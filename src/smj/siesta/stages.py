"""SIESTA 5.4.2 molecular, Au principal-layer, device and open-boundary stages."""
from __future__ import annotations

import math
from pathlib import Path
import re
import shlex
import shutil
from uuid import uuid4

import numpy as np

from ..checks import ATOMIC_NUMBERS
from ..errors import ConfigError, StageError
from ..junction.geometry import read_xyz, write_xyz
from ..stages._helpers import resolve_anchor_indices_atoms
from ..stages.base import Stage
from ..stages.molecule_opt import MoleculeOptStage
from .geometry import align_device_positions, build_device, build_electrode, inspect_device
from .io import (NUM, VERSION, convergence, dump, inspect_hsx, inspect_restart, mulliken_populations,
                 open_charge, read_json, relaxed_geometry, sha, verify_files)

DEFAULTS = {"mesh_cutoff_ry": 300., "smearing_ev": 0.025, "scf_tolerance": 1e-4,
            "force_ev_ang": 0.02, "electrode_kgrid": [1, 1, 12], "device_kgrid": [1, 1, 1],
            "spin": "unpolarized", "molecule_vacuum_ang": 10., "relax_steps": 200,
            "scf_steps": 300}
MODULE = "siesta/5.4.2-gpu-netcdf"
# Pre-registered before seeing the 5x5 result: a device Au layer counts as
# bulk-like when its Mulliken mean is within 0.05 e of the periodic electrode
# mean, and the periodic electrode's own layer spread is below the same value.
SCREENING_TOLERANCE_E = 0.05
SCREENING_CRITERION = ("outermost two device Au layers within %.2f e of the periodic "
                       "electrode Mulliken mean; electrode layer spread <= %.2f e"
                       % (SCREENING_TOLERANCE_E, SCREENING_TOLERANCE_E))


def validate_settings(settings):
    opts = DEFAULTS | settings
    if "device_mixing_weight" in opts:
        w = opts["device_mixing_weight"]
        if type(w) not in (int, float) or not 0 < w <= 1:
            raise ConfigError("siesta.device_mixing_weight must be in (0, 1]")
    for key in ("mesh_cutoff_ry", "smearing_ev", "scf_tolerance", "force_ev_ang", "molecule_vacuum_ang"):
        v = opts[key]
        if isinstance(v, bool) or not isinstance(v, (float, int)) or not math.isfinite(v) or v <= 0:
            raise ConfigError(f"siesta.{key} must be positive and finite")
    for key in ("scf_steps", "relax_steps"):
        if type(opts[key]) is not int or opts[key] < 1:
            raise ConfigError(f"siesta.{key} must be a positive integer")
    for key in ("electrode_kgrid", "device_kgrid"):
        grid = opts[key]
        if (not isinstance(grid, (tuple, list)) or len(grid) != 3 or
                any(type(x) is not int or x < 1 for x in grid)):
            raise ConfigError(f"siesta.{key} must contain three positive integers")
    if opts["device_kgrid"][2] != 1:
        raise ConfigError("siesta.device_kgrid needs a single transport-axis k point")
    if list(opts["device_kgrid"][:2]) != list(opts["electrode_kgrid"][:2]):
        raise ConfigError("electrode_kgrid and device_kgrid require matching transverse sampling")
    if opts["spin"] not in ("unpolarized", "collinear"):
        raise ConfigError("siesta.spin must be unpolarized or collinear")
    for key in ("basis", "pseudopotentials"):
        values = opts.get(key, {})
        if not isinstance(values, dict):
            raise ConfigError(f"siesta.{key} must be an element mapping")
        for symbol, val in values.items():
            if symbol not in ATOMIC_NUMBERS or not isinstance(val, str) or not val.strip():
                raise ConfigError(f"invalid siesta.{key} entry {symbol!r}")
            if key == "basis" and val not in ("SZ", "DZ", "SZP", "DZP", "TZ", "TZP"):
                raise ConfigError(f"unsupported PAO size for {symbol}: {val}")
            if key == "pseudopotentials" and Path(val).suffix.lower() not in (".psf", ".psml"):
                raise ConfigError(f"expected .psf/.psml pseudopotential for {symbol}")
    if opts.get("electrode_reuse") is not None and not isinstance(opts["electrode_reuse"], str):
        raise ConfigError("siesta.electrode_reuse must be a collected electrode directory")


def pseudo_paths(ctx, elements):
    entries = ctx.config.siesta.get("pseudopotentials", {})
    paths = {}
    for element in sorted(set(elements)):
        path = ctx.resolve(entries.get(element))
        if path is None or not path.is_file() or not path.stat().st_size:
            raise ConfigError(f"missing readable siesta.pseudopotentials.{element}: {path}")
        if element not in ctx.config.siesta.get("basis", {}):
            raise ConfigError(f"missing siesta.basis.{element}")
        paths[element] = path
    return paths


def launch(ctx, executable, fdf, *options):
    args = shlex.join([executable, *options, fdf])
    return 'mpirun -np "${SLURM_NTASKS:-%d}" %s' % (ctx.config.cluster.ntasks or 1, args)


def slurm_directives(ctx, stage_name):
    c = ctx.config.cluster
    args = {"job-name": stage_name, "ntasks": c.ntasks or 1, "partition": c.partition,
            "time": c.time_limit, "mem": c.mem, "gres": c.gres, "account": c.account, "qos": c.qos}
    lines = []
    for key, val in args.items():
        if val is None:
            continue
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_:,.+-]*", str(val)):
            raise ConfigError(f"invalid Slurm --{key}: {val!r}")
        lines.append(f"#SBATCH --{key}={val}")
    return lines


def write_script(ctx, stage_dir, stage_name, command_lines, clean):
    lines = ["#!/bin/bash -l", *slurm_directives(ctx, stage_name)]
    lines += ["set -eo pipefail", """trap 'rc=$?; printf "%s\\n" "$rc" > .dsh_done' EXIT""",
              "rm -f -- " + shlex.join([".dsh_done", *clean]),
              f"module load {MODULE}", "ulimit -s unlimited", *command_lines]
    path = Path(stage_dir) / "sub_siesta.sh"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def _block(name, rows):
    return "\n".join([f"%block {name}", *rows, f"%endblock {name}"])


def _xyz(atoms):
    return {"elements": [a.element for a in atoms], "positions": [list(a.coord()) for a in atoms]}


def _format_vec(row):
    return " ".join(f"{v:.10f}" for v in row)


def _fdf(opts, geometry, pseudo_names, label, grid, spin):
    elements = geometry["elements"]
    ids = {name: i for i, name in enumerate(sorted(set(elements)), 1)}
    return [f"SystemLabel {label}", f"NumberOfAtoms {len(elements)}", f"NumberOfSpecies {len(ids)}",
            _block("ChemicalSpeciesLabel",
                   [f"{idx} {ATOMIC_NUMBERS[e]} {e} {pseudo_names[e]}" for e, idx in ids.items()]),
            _block("PAO.BasisSizes", [f"{e} {opts['basis'][e]}" for e in ids]),
            "XC.functional GGA", "XC.authors PBE", f"MeshCutoff {opts['mesh_cutoff_ry']} Ry",
            f"ElectronicTemperature {opts['smearing_ev']} eV",
            f"Spin {'polarized' if spin == 'collinear' else 'non-polarized'}",
            "Charge.Spin end",
            "SCF.MustConverge true", f"SCF.DM.Tolerance {opts['scf_tolerance']}",
            f"MaxSCFIterations {opts['scf_steps']}", "DM.UseSaveDM false",
            "LatticeConstant 1.0 Ang", _block("LatticeVectors", map(_format_vec, geometry["cell"])),
            "AtomicCoordinatesFormat Ang",
            _block("AtomicCoordinatesAndAtomicSpecies",
                   [f"{_format_vec(pos)} {ids[e]}" for pos, e in zip(geometry["positions"], elements)]),
            _block("kgrid.MonkhorstPack", [
                " ".join(str(n if i == j else 0) for j in range(3)) + " 0.0"
                for i, n in enumerate(grid)])]


def _relax(opts):
    return ["MD.TypeOfRun CG", f"MD.Steps {opts['relax_steps']}", "MD.VariableCell false",
            "MD.UseSaveXV false", "MD.UseSaveCG false", "GeometryMustConverge true",
            f"MD.MaxForceTol {opts['force_ev_ang']} eV/Ang"]


def _file_digests(directory, names):
    return {name: sha(directory / name) for name in names}


def _sha_manifest(path, expected):
    """Parse the restricted basename-only format generated by sha256sum."""
    try:
        data = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            match = re.fullmatch(r"([0-9a-f]{64})  ([\w.-]+)", line, flags=re.ASCII)
            if not match or match[2] in data:
                raise ValueError("invalid sha256sum record")
            data[match[2]] = match[1]
        if set(data) != set(expected):
            raise ValueError("unexpected or missing file records")
        verify_files(path.parent, {"files": data})
        return data
    except (OSError, ValueError) as exc:
        raise StageError(f"{path}: {exc}") from exc


def _upstream(ctx, name, manifest, required):
    directory = ctx.stage_dir(name)
    data = read_json(directory / manifest)
    if not set(required).issubset(data.get("files", {})):
        raise StageError(f"{name}/{manifest}: missing required file digests")
    verify_files(directory, data)
    return directory, data


class SiestaStage(Stage):
    script = "sub_siesta.sh"
    options = ()
    raw = ()

    @property
    def remote_files(self):
        return (*self.raw, "output.sha256")

    def write_inputs(self, ctx, d, geom, grid, spin, extra, protocol, parents=()):
        d.mkdir(parents=True, exist_ok=True)
        opts = DEFAULTS | ctx.config.siesta
        pseudos = pseudo_paths(ctx, geom["elements"])
        names = {e: e + p.suffix.lower() for e, p in pseudos.items()}
        for e, p in pseudos.items():
            if p.resolve() != (d / names[e]).resolve():
                shutil.copyfile(p, d / names[e])
        dump(d / "geometry.json", geom)
        fdf = self.fdf
        lines = _fdf(opts, geom, names, self.label, grid, spin) + extra
        if self.label == "device" and "device_mixing_weight" in opts:
            lines.extend([f"SCF.Mixer.Weight {opts['device_mixing_weight']}",
                          "SCF.Mixer.History 6", "SCF.Mix.First true",
                          "SCF.Mix.First.Force true"])
        (d / fdf).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        files = [fdf, "geometry.json", *names.values(), *parents]
        protocol.update(version=VERSION, spin=spin, attempt=uuid4().hex,
                        pseudo_sources={e: str(p) for e, p in pseudos.items()},
                        files=_file_digests(d, files))
        dump(d / "protocol.json", protocol)
        files.append("protocol.json")
        (d / "input.sha256").write_text(
            "".join(f"{sha(d / n)}  {n}\n" for n in files), encoding="ascii", newline="\n")
        files.append("input.sha256")
        command = launch(ctx, "siesta", fdf, *self.options) + " > siesta.out 2>&1"
        lines = ["sha256sum -c input.sha256", command,
                 "sha256sum " + shlex.join(["input.sha256", *self.raw]) + " > output.sha256"]
        clean = [*self.remote_files, self.label + ".DM", self.label + ".XV",
                 self.label + ".CG", self.label + ".TSDE", "0_NORMAL_EXIT"]
        script = write_script(ctx, d, self.name, lines, clean)
        return self.prepared(ctx, d, inputs=[d / n for n in files] + [script],
                             commands=[f"sbatch {script.name}"])

    def read_outputs(self, d):
        protocol = read_json(d / "protocol.json")
        inputs = _sha_manifest(d / "input.sha256", [*protocol["files"], "protocol.json"])
        verify_files(d, protocol)
        outputs = _sha_manifest(d / "output.sha256", ["input.sha256", *self.raw])
        log = (d / "siesta.out").read_text(encoding="utf-8", errors="replace")
        result = convergence(log, open_system=self.name == "transiesta_scf")
        result.update(spin=protocol["spin"], files=inputs | outputs,
                      attempt=protocol["attempt"])
        if result["spin"] == "collinear" and result["final_moment_muB"] is None:
            raise StageError(f"{self.name}: missing final spin moment")
        if result["spin"] == "unpolarized":
            result["final_moment_muB"] = 0.0
        return protocol, result, log

    def finish(self, ctx, d, result, manifest, extras=()):
        result["files"].update(_file_digests(d, extras))
        dump(d / manifest, result)
        return self.collected(ctx, d, artifacts=[manifest, *result["files"]], results=result)

    def read_relax(self, d):
        protocol, result, log = self.read_outputs(d)
        matches = re.findall(rf"(?im)^\s*Max\s+({NUM})\s+constrained\s*$", log)
        if not matches or not (0 <= float(matches[-1]) <= protocol["force_ev_ang"] + 5e-7):
            raise StageError(f"{self.name}: final constrained force missing or above threshold")
        result.update(max_force_ev_ang=float(matches[-1]), geometry_converged=True,
                      stop_reason="GeometryMustConverge with final force within MD.MaxForceTol")
        geom = read_json(d / "geometry.json")
        atoms, cell = relaxed_geometry(d / (self.label + ".XV"), geom)
        return protocol, result, geom, atoms, cell


class Molecule(SiestaStage):
    name = "molecule_opt"
    label = "molecule"
    fdf = "molecule.fdf"
    raw = ("siesta.out", "molecule.XV")
    outputs = ("opt.xyz", "molecule.json")
    input_hint = "molecule.xyz and SIESTA species"

    def prepare(self, ctx, stage_dir):
        opts = DEFAULTS | ctx.config.siesta
        xyz = ctx.molecule_xyz()
        mol = ctx.config.molecule
        atoms = read_xyz(xyz)
        if int(xyz.read_text(encoding="utf-8").splitlines()[0]) != len(atoms):
            raise StageError("XYZ atom count differs from parsed rows")
        protocol = MoleculeOptStage._electron_state(xyz, mol)
        anchors = resolve_anchor_indices_atoms(atoms, ctx.config.electrode.left_anchor,
                                               ctx.config.electrode.right_anchor, stage=self.name)
        geom = _xyz(atoms)
        pos = np.asarray(geom["positions"])
        length = float(np.ptp(pos, axis=0).max()) + 2 * opts["molecule_vacuum_ang"]
        geom.update(cell=(np.eye(3)*length).tolist(),
                    positions=(pos - (pos.min(0)+pos.max(0))/2 + length/2).tolist())
        open_spin = mol.multiplicity != 1 or mol.broken_symmetry or bool(mol.initial_moments)
        spin = "collinear" if open_spin else "unpolarized"
        extra = [*_relax(opts), f"NetCharge {mol.charge}"]
        if open_spin and mol.fixed_spin_moment:
            extra += ["Spin.Fix true", f"Spin.Total {mol.multiplicity-1}"]
        moments = mol.initial_moments or {}
        if moments:
            if max(moments) > len(atoms):
                raise ConfigError("initial magnetic moment atom outside XYZ")
            extra.append(_block("DM.InitSpin", [f"{i} {v}" for i, v in moments.items()]))
        protocol.update(xyz_sha256=sha(xyz), anchor_ids=list(anchors),
                        initial_moments_muB=moments, fixed_spin=bool(open_spin and mol.fixed_spin_moment),
                        force_ev_ang=opts["force_ev_ang"], vacuum_ang=opts["molecule_vacuum_ang"],
                        boundary="periodic cubic supercell at Gamma",
                        charge_treatment="periodic compensating background; cubic molecular Madelung correction requested"
                        if mol.charge else "neutral cell")
        return self.write_inputs(ctx, stage_dir, geom, [1,1,1], spin, extra, protocol)

    def collect(self, ctx, stage_dir):
        protocol, result, _, atoms, _ = self.read_relax(stage_dir)
        if protocol["fixed_spin"] and abs(result["final_moment_muB"] - (protocol["multiplicity"]-1)) > 0.02:
            raise StageError("molecule_opt: final spin incompatible with requested multiplicity")
        if protocol["charge"]:
            log = (stage_dir / "siesta.out").read_text(encoding="utf-8", errors="replace")
            result["cubic_madelung_correction_confirmed"] = "madelung: Charged system with SC cell" in log
        result.update(anchor_ids=protocol["anchor_ids"], original_ids=list(range(1, len(atoms)+1)))
        write_xyz(stage_dir / "opt.xyz", atoms, "SIESTA converged molecule, original atom order")
        return self.finish(ctx, stage_dir, result, "molecule.json", ["opt.xyz"])


class Electrode(SiestaStage):
    name = "siesta_electrode"
    label = "electrode"
    options = ("--electrode",)
    fdf = "electrode.fdf"
    raw = ("siesta.out", "electrode.HSX", "electrode.TSDE")
    outputs = ("electrode.HSX", "electrode.TSDE", "electrode.json")
    input_hint = "periodic Au(111) principal layer or validated electrode_reuse directory"

    def prepare(self, ctx, stage_dir):
        opts = DEFAULTS | ctx.config.siesta
        source = ctx.resolve(opts.get("electrode_reuse"))
        if source and source.resolve() == stage_dir.resolve():
            raise StageError("siesta.electrode_reuse must reference a separate collected stage")
        geom = build_electrode(opts)
        prepared = self.write_inputs(ctx, stage_dir, geom, opts["electrode_kgrid"], opts["spin"],
                                     ["MD.Steps 0", "SaveHS true", "TS.DE.Save true",
                                      "WriteMullikenPop 1"], {})
        if source is None:
            return prepared
        saved = read_json(source / "electrode.json")
        geom_files = read_json(stage_dir / "protocol.json")["files"]
        if (saved.get("version") != VERSION or saved.get("scf_converged") is not True
                or any(saved.get("files", {}).get(n) != v for n, v in geom_files.items())):
            raise StageError("incompatible or unconverged electrode_reuse: resolved electrode inputs differ")
        verify_files(source, saved)
        # Re-read the source output and its paired input/attempt receipt. A checksum
        # of an unrelated HSX or a manually edited convergence flag is insufficient.
        _, result, _ = self.read_outputs(source)
        hsx = inspect_hsx(source / "electrode.HSX", geom, result["spin"], electrode=True)
        inspect_restart(source / "electrode.TSDE", hsx, result["spin"])
        for name in [*self.remote_files, "input.sha256", "protocol.json"]:
            shutil.copyfile(source / name, stage_dir / name)
        prepared.local_only = True
        prepared.script = None
        prepared.commands = []
        prepared.notes = [f"Reused verified electrode from {source}"]
        return prepared

    def collect(self, ctx, stage_dir):
        _, result, log = self.read_outputs(stage_dir)
        geom = read_json(stage_dir / "geometry.json")
        result.update(inspect_hsx(stage_dir / "electrode.HSX", geom, result["spin"], electrode=True))
        inspect_restart(stage_dir / "electrode.TSDE", result, result["spin"])
        source = ctx.resolve(ctx.config.siesta.get("electrode_reuse"))
        result.update(principal_layer_checked=True, source=str(source or stage_dir))
        try:
            pops = mulliken_populations(log)
        except StageError:
            result.update(mulliken_bulk_charge_e=None, mulliken_layer_spread_e=None)
        else:
            _mulliken_charges(geom, pops)
            positions = np.asarray(geom["positions"], dtype=float)
            layers = {}
            for row in pops:
                layers.setdefault(round(float(positions[row[0] - 1][2]), 3), []).append(row[1])
            means = [float(np.mean(values)) for _, values in sorted(layers.items())]
            core = means[1:-1] if len(means) > 3 else means
            result.update(mulliken_bulk_charge_e=float(np.mean(core)),
                          mulliken_layer_spread_e=max(core) - min(core))
        return self.finish(ctx, stage_dir, result, "electrode.json")


def _require_electrode(ctx):
    return _upstream(ctx, "siesta_electrode", "electrode.json",
                     ["geometry.json", "electrode.HSX", "electrode.TSDE", "electrode.fdf"])


def _moments(ctx, device, spin):
    moments = ctx.config.molecule.initial_moments or {}
    if spin != "collinear" or not moments:
        return []
    mapping = device["mapping"]["original_to_device"]
    missing = [i for i in moments if str(i) not in mapping]
    if missing:
        raise ConfigError(f"initial magnetic moments refer to removed/absent atoms: {missing}")
    return [_block("DM.InitSpin", [f"{mapping[str(i)]} {v}" for i, v in moments.items()])]


def _mulliken_charges(geom, pops):
    elements = geom["elements"]
    if len(pops) != len(elements) or [row[0] for row in pops] != list(range(1, len(elements) + 1)):
        raise StageError("Mulliken table does not follow the device atom order")
    if [row[3] for row in pops] != elements:
        raise StageError("Mulliken table species do not match the device geometry")
    return {row[0]: row[1] for row in pops}


def _layer_screening(geom, charges, bulk, spread):
    report = {"status": "unverified", "criterion": SCREENING_CRITERION,
              "tolerance_e": SCREENING_TOLERANCE_E,
              "electrode_bulk_charge_e": bulk,
              "electrode_layer_spread_e": spread, "sides": {}}
    if not isinstance(bulk, (int, float)) or not math.isfinite(float(bulk)):
        report["reason"] = "periodic electrode has no Mulliken reference"
        return report
    if not isinstance(spread, (int, float)) or not math.isfinite(float(spread)):
        report["reason"] = "periodic electrode has no layer-spread reference"
        return report
    if spread > SCREENING_TOLERANCE_E:
        report["reason"] = "periodic electrode layer charges are not uniform"
        return report
    positions = np.asarray(geom["positions"], dtype=float)
    contacts = geom.get("contacts") or {}
    for side in ("left", "right"):
        layer_ids = (contacts.get(side) or {}).get("layer_ids") or []
        if not layer_ids:
            continue
        profile = []
        for ids in layer_ids:
            profile.append({"atoms": len(ids), "z_A": float(np.mean([positions[i - 1][2] for i in ids])),
                            "mean_charge_e": float(np.mean([charges[i] for i in ids]))})
        profile.sort(key=lambda row: row["z_A"])
        if side == "right":
            profile.reverse()  # both sides: index 0 is the outermost (lead-facing) layer
        deltas = [row["mean_charge_e"] - float(bulk) for row in profile]
        for row, delta in zip(profile, deltas):
            row["delta_e"] = delta
        outer = deltas[0]
        second = deltas[1] if len(deltas) > 1 else None
        passed = abs(outer) <= SCREENING_TOLERANCE_E and (
            second is None or abs(second) <= SCREENING_TOLERANCE_E)
        report["sides"][side] = {"passed": passed, "outermost_delta_e": outer,
                                 "second_outermost_delta_e": second,
                                 "innermost_delta_e": deltas[-1], "layers": profile}
    if not report["sides"]:
        report["reason"] = "device contact layer map is missing"
        return report
    report["status"] = "passed" if all(row["passed"] for row in report["sides"].values()) else "failed"
    report["reason"] = ("outer two Au layers are bulk-like within the pre-registered tolerance"
                        if report["status"] == "passed" else
                        "outer Au layer deviates from the periodic electrode beyond tolerance")
    return report


def _relax_reuse_dir(ctx):
    raw = ctx.config.siesta.get("reuse_relax_from")
    if not raw:
        return None
    root = ctx.resolve(raw)
    if root is None:
        return None
    candidates = [root] if root.name == "siesta_junction_relax" else [root, root / "siesta_junction_relax"]
    for path in candidates:
        if (path / "relax.json").is_file():
            return path
    raise StageError("siesta.reuse_relax_from does not contain a collected siesta_junction_relax/relax.json")


def _scf_reuse_dir(ctx):
    raw = ctx.config.siesta.get("reuse_scf_from")
    if not raw:
        return None
    root = ctx.resolve(raw)
    if root is None:
        return None
    candidates = [root] if root.name == "transiesta_scf" else [root, root / "transiesta_scf"]
    for path in candidates:
        if (path / "scf.json").is_file():
            return path
    raise StageError("siesta.reuse_scf_from does not contain a collected transiesta_scf/scf.json")


class JunctionRelax(SiestaStage):
    name = "siesta_junction_relax"
    label = "device"
    fdf = "device.fdf"
    raw = ("siesta.out", "device.XV")
    depends = ("molecule_opt", "siesta_electrode")
    outputs = ("final_geometry.json", "contacts.json", "relax.json")
    input_hint = "optimized molecule and verified periodic electrode"

    @classmethod
    def dependencies(cls, ctx):
        config = getattr(ctx, "config", ctx)
        if not config.siesta.get("reuse_relax_from"):
            return cls.depends
        state = getattr(ctx, "state", None)
        if state is None:
            return cls.depends
        molecule = state.stage("molecule_opt")
        if getattr(molecule, "status", None) == "collected":
            return cls.depends
        return ("siesta_electrode",)

    def prepare(self, ctx, stage_dir):
        source = _relax_reuse_dir(ctx)
        if source is not None:
            return self._prepare_reuse(ctx, stage_dir, source)
        opts = DEFAULTS | ctx.config.siesta
        mol_dir, mol = _upstream(ctx, "molecule_opt", "molecule.json", ["opt.xyz"])
        lead_dir, lead = _require_electrode(ctx)
        atoms = read_xyz(mol_dir / "opt.xyz")
        anchors = tuple(mol["anchor_ids"])
        strip = ctx.config.electrode.strip_thiol_h
        if strip is None:
            strip = all(atoms[i-1].element == "S" for i in anchors)
        geom = build_device(atoms, anchors, read_json(lead_dir / "geometry.json"), opts, strip_h=strip)
        extra = _relax(opts) + [_block("Geometry.Constraints", [f"atom {i}" for i in geom["fixed_ids"]])]
        extra += _moments(ctx, geom, opts["spin"])
        protocol = {"force_ev_ang": opts["force_ev_ang"], "opt_xyz_sha256": mol["files"]["opt.xyz"],
                    "electrode_hsx_sha256": lead["files"]["electrode.HSX"],
                    "electrode_geometry_sha256": lead["files"]["geometry.json"]}
        return self.write_inputs(ctx, stage_dir, geom, opts["device_kgrid"], opts["spin"], extra, protocol)

    def _prepare_reuse(self, ctx, stage_dir, source):
        if source.resolve() == stage_dir.resolve():
            raise StageError("siesta.reuse_relax_from must reference a separate collected stage")
        relax = read_json(source / "relax.json")
        verify_files(source, relax)
        contacts = read_json(source / "contacts.json")
        if not contacts.get("passed"):
            raise StageError("reused siesta_junction_relax does not have passed contacts")
        protocol_file = source / "protocol.json"
        protocol = read_json(protocol_file) if protocol_file.is_file() else {}
        electrode_geometry = (protocol.get("electrode_geometry_sha256")
                              or relax.get("electrode_geometry_sha256"))
        if not electrode_geometry:
            raise StageError("reused siesta_junction_relax has no electrode geometry digest")
        local_mol = ctx.stage_dir("molecule_opt") / "molecule.json"
        if local_mol.is_file():
            _, mol = _upstream(ctx, "molecule_opt", "molecule.json", ["opt.xyz"])
        else:
            # A TBtrans-only run may not carry its own molecule stage; the source
            # run's molecule receipt still proves which optimized geometry was used.
            source_mol = source.parent / "molecule_opt"
            mol = read_json(source_mol / "molecule.json")
            verify_files(source_mol, mol)
        if relax.get("opt_xyz_sha256") != mol["files"]["opt.xyz"]:
            raise StageError("reused relaxation was prepared for a different optimized molecule")
        lead_dir, lead = _require_electrode(ctx)
        if electrode_geometry != lead["files"]["geometry.json"]:
            raise StageError("reused relaxation was prepared for a different electrode geometry")
        stage_dir.mkdir(parents=True, exist_ok=True)
        for name in ("final_geometry.json", "contacts.json"):
            shutil.copyfile(source / name, stage_dir / name)
        marker = {"source": str(source), "electrode_geometry_sha256": electrode_geometry,
                  "electrode_hsx_sha256": relax.get("electrode_hsx_sha256"),
                  "opt_xyz_sha256": relax.get("opt_xyz_sha256"),
                  "contacts": contacts.get("contacts") or relax.get("contacts") or {},
                  "files": {name: sha(stage_dir / name)
                            for name in ("final_geometry.json", "contacts.json")}}
        dump(stage_dir / "reuse.json", marker)
        return self.prepared(ctx, stage_dir, local_only=True,
                             inputs=[stage_dir / name for name in marker["files"]] + [stage_dir / "reuse.json"],
                             notes=[f"Reused relaxed junction geometry from {source}"])

    def collect(self, ctx, stage_dir):
        marker = stage_dir / "reuse.json"
        if marker.is_file():
            info = read_json(marker)
            for name, expected in info["files"].items():
                if sha(stage_dir / name) != expected:
                    raise StageError(f"reused siesta_junction_relax: {name} changed after prepare")
            result = {"reused_from": info["source"],
                      "electrode_geometry_sha256": info["electrode_geometry_sha256"],
                      "electrode_hsx_sha256": info["electrode_hsx_sha256"],
                      "opt_xyz_sha256": info["opt_xyz_sha256"],
                      "contacts": info["contacts"], "files": dict(info["files"])}
            return self.finish(ctx, stage_dir, result, "relax.json")
        protocol, result, geom, atoms, cell = self.read_relax(stage_dir)
        lead_dir, lead = _require_electrode(ctx)
        if lead["files"]["geometry.json"] != protocol["electrode_geometry_sha256"]:
            raise StageError("electrode geometry changed since junction preparation")
        report = inspect_device(geom, read_json(lead_dir / "geometry.json"), atoms, cell)
        dump(stage_dir / "contacts.json", report)
        if not report["passed"]:
            raise StageError(f"{self.name}: " + ", ".join(report["issues"]))
        final = geom | {"positions": align_device_positions(geom, atoms, cell), "cell": cell}
        dump(stage_dir / "final_geometry.json", final)
        result.update(electrode_hsx_sha256=protocol["electrode_hsx_sha256"],
                      opt_xyz_sha256=protocol["opt_xyz_sha256"], contacts=report["contacts"])
        return self.finish(ctx, stage_dir, result, "relax.json", ["contacts.json", "final_geometry.json"])


def _open_boundary(geom):
    lines = ["SolutionMethod transiesta", "TS.Voltage 0.0 eV", "SaveHS true", "TS.DE.Save true",
             "TS.Elecs.Neglect.Principal false", "TS.Elecs.Gf.Reuse false",
             _block("TS.Elecs", ["Left", "Right"]), _block("TS.ChemPots", ["Equilibrium"]),
             _block("TS.ChemPot.Equilibrium", ["mu 0.0 eV", "contour.eq.pole.N 100"])]
    for side, direction in [("Left", "-a3"), ("Right", "+a3")]:
        first = geom["electrode_ids"][side.lower()][0]
        lines.append(_block(f"TS.Elec.{side}", ["HS electrode.HSX", "DE electrode.TSDE",
                            f"semi-inf-dir {direction}", f"elec-pos {first}", "chem-pot Equilibrium"]))
    return lines


class OpenSCF(SiestaStage):
    name = "transiesta_scf"
    label = "device"
    fdf = "transiesta.fdf"
    raw = ("siesta.out", "device.TS.HSX", "device.TSDE")
    depends = ("siesta_electrode", "siesta_junction_relax")
    outputs = ("device.TS.HSX", "device.TSDE", "scf.json")
    input_hint = "relaxed matched device and both electrode references"

    def prepare(self, ctx, stage_dir):
        source = _scf_reuse_dir(ctx)
        if source is not None:
            return self._prepare_reuse(ctx, stage_dir, source)
        opts = DEFAULTS | ctx.config.siesta
        lead_dir, lead = _require_electrode(ctx)
        device_dir, relax = _upstream(ctx, "siesta_junction_relax", "relax.json",
                                     ["final_geometry.json", "contacts.json"])
        contacts_passed = read_json(device_dir / "contacts.json").get("passed")
        if relax.get("reused_from"):
            # A reused relaxation keeps the electrode geometry but deliberately
            # changes the transverse k mesh, so its HSX digest must differ.
            geometry_matches = (relax.get("electrode_geometry_sha256")
                                == lead["files"]["geometry.json"])
        else:
            geometry_matches = relax.get("electrode_hsx_sha256") == lead["files"]["electrode.HSX"]
        if not contacts_passed or not geometry_matches:
            raise StageError("open device does not match the collected electrode or contacts")
        geom = read_json(device_dir / "final_geometry.json")
        stage_dir.mkdir(parents=True, exist_ok=True)
        parents = ["electrode.HSX", "electrode.TSDE"]
        for name in parents:
            shutil.copyfile(lead_dir / name, stage_dir / name)
        extra = ["MD.Steps 0", "WriteMullikenPop 1", *_open_boundary(geom),
                 *_moments(ctx, geom, opts["spin"])]
        protocol = {"geometry_sha256": relax["files"]["final_geometry.json"],
                    "contacts_sha256": relax["files"]["contacts.json"]}
        return self.write_inputs(ctx, stage_dir, geom, opts["device_kgrid"], opts["spin"], extra,
                                 protocol, parents)

    def _prepare_reuse(self, ctx, stage_dir, source):
        if source.resolve() == stage_dir.resolve():
            raise StageError("siesta.reuse_scf_from must reference a separate collected stage")
        scf = read_json(source / "scf.json")
        verify_files(source, scf)
        if scf.get("scf_converged") is not True or scf.get("open_boundary_converged") is not True:
            raise StageError("reused transiesta_scf is not converged or lacks open-boundary evidence")
        opts = DEFAULTS | ctx.config.siesta
        if scf.get("spin") != opts["spin"]:
            raise StageError("reused transiesta_scf was computed for a different spin setting")
        device_dir, relax = _upstream(ctx, "siesta_junction_relax", "relax.json",
                                      ["final_geometry.json", "contacts.json"])
        if read_json(device_dir / "contacts.json").get("passed") is not True:
            raise StageError("reused transiesta_scf requires passed contacts")
        if scf.get("geometry_sha256") != relax["files"]["final_geometry.json"]:
            raise StageError("reused transiesta_scf was prepared for a different junction geometry")
        if scf.get("contacts_sha256") != relax["files"]["contacts.json"]:
            raise StageError("reused transiesta_scf was prepared for different contacts")
        lead_dir, lead = _require_electrode(ctx)
        required = {"transiesta.fdf", "device.TS.HSX", "device.TSDE", "electrode.HSX",
                    "electrode.TSDE", "siesta.out"}
        missing = required - set(scf.get("files", {}))
        if missing:
            raise StageError("reused transiesta_scf is missing: " + ", ".join(sorted(missing)))
        if scf["files"]["electrode.HSX"] != lead["files"]["electrode.HSX"]:
            raise StageError("reused transiesta_scf was prepared with a different electrode")
        names = sorted(scf["files"])
        for name in names:
            if Path(name).name != name:
                raise StageError("reused transiesta_scf contains a non-local artifact name")
        stage_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            shutil.copyfile(source / name, stage_dir / name)
        marker = {"source": str(source), "scf_sha256": sha(source / "scf.json"), "scf": scf,
                  "files": {name: sha(stage_dir / name) for name in names}}
        dump(stage_dir / "reuse.json", marker)
        return self.prepared(ctx, stage_dir, local_only=True,
                             inputs=[stage_dir / name for name in names] + [stage_dir / "reuse.json"],
                             notes=[f"Reused converged TranSIESTA SCF from {source}"])

    def collect(self, ctx, stage_dir):
        marker = stage_dir / "reuse.json"
        if marker.is_file():
            info = read_json(marker)
            source = _scf_reuse_dir(ctx)
            if source is None or source.resolve() != Path(info["source"]).resolve():
                raise StageError("siesta.reuse_scf_from changed after prepare")
            for name, expected in info["files"].items():
                if sha(stage_dir / name) != expected:
                    raise StageError(f"reused transiesta_scf: {name} changed after prepare")
            result = dict(info["scf"])
            result["reused_from"] = info["source"]
            result["files"] = dict(info["files"])
            return self.finish(ctx, stage_dir, result, "scf.json")
        protocol, result, log = self.read_outputs(stage_dir)
        geom = read_json(stage_dir / "geometry.json")
        result.update(inspect_hsx(stage_dir / "device.TS.HSX", geom, result["spin"]))
        inspect_restart(stage_dir / "device.TSDE", result, result["spin"])
        target, charge = open_charge(log, result["spin"])
        lead_dir, lead = _require_electrode(ctx)
        try:
            pops = mulliken_populations(log)
        except StageError as exc:
            screening = {"status": "unverified", "reason": str(exc),
                         "criterion": SCREENING_CRITERION, "tolerance_e": SCREENING_TOLERANCE_E,
                         "sides": {}}
        else:
            screening = _layer_screening(geom, _mulliken_charges(geom, pops),
                                         lead.get("mulliken_bulk_charge_e"),
                                         lead.get("mulliken_layer_spread_e"))
        result.update(open_boundary_converged=True, target_charge_e=target, charge_e=charge,
                      charge_deviation_e=charge-target, screening=screening,
                      geometry_sha256=protocol["geometry_sha256"], contacts_sha256=protocol["contacts_sha256"],
                      chemical_potentials_ev={"Left": result["fermi_ev"], "Right": result["fermi_ev"]})
        return self.finish(ctx, stage_dir, result, "scf.json")


STAGES = {stage.name: stage for stage in (Molecule(), Electrode(), JunctionRelax(), OpenSCF())}


def remote_files(ctx, name):
    return STAGES[name].remote_files
