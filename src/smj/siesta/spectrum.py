"""TBtrans NetCDF extraction and zero-bias, two-terminal Landauer conductance."""
from __future__ import annotations

from pathlib import Path
import math
import re
import shutil
import uuid

import numpy as np

from ..config import HARTREE_EV
from ..errors import DependencyMissing, StageError
from ..junction.criteria import MAX_FERMI_MOMENT_KBT, MIN_QUADRATURE_RATIO, MIN_THERMAL_COVERAGE
from ..junction.landauer import fermi_quadrature, minus_fermi_derivative, thermal_window, trapz
from ..stages.base import Stage
from .io import dump, read_json, sha, verify_files


# SI defining constants, G0 = 2 e^2/h.
G0_S = 2 * 1.602176634e-19 ** 2 / 6.62607015e-34


def _data(variable):
    data = variable[:]
    if np.ma.is_masked(data) or not np.all(np.isfinite(data)):
        raise ValueError(f"incomplete or nonfinite values in {variable.name}")
    return np.asarray(data, dtype=float)


def _energy(variable):
    unit = getattr(variable, "unit", getattr(variable, "units", None))
    factors = {"ry": HARTREE_EV / 2, "ev": 1.0, "ha": HARTREE_EV}
    if str(unit).lower() not in factors:
        raise ValueError(f"unknown energy unit {unit!r} for {variable.name}")
    return _data(variable) * factors[unit.lower()], unit


def _read_channel(path):
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise DependencyMissing("TBtrans NetCDF analysis requires pip install 'smj-toolkit[siesta]'") from exc
    try:
        with Dataset(path) as ds:
            energy, unit = _energy(ds.variables["E"])
            if energy.ndim != 1 or len(energy) < 2 or np.any(np.diff(energy) <= 0):
                raise ValueError("expected a strictly increasing TBtrans energy grid")
            potentials = []
            for name in ("Left", "Right"):
                mu, _ = _energy(ds.groups[name].variables["mu"])
                if mu.size != 1:
                    raise ValueError(f"invalid {name} chemical potential")
                potentials.append(mu.item())
            if not np.allclose(potentials, 0., rtol=0, atol=1e-8):
                raise ValueError(f"expected equilibrium electrodes at zero reference, got {potentials} eV")
            variable = ds.groups["Left"].variables["Right.T"]
            data = _data(variable)
            if np.any(data < -1e-12):
                raise ValueError("negative left-to-right transmission")
            data = np.maximum(data, 0.)
            e_dim = ds.variables["E"].dimensions[0]
            weights = points = None
            if "wkpt" in ds.variables or "kpt" in ds.variables:
                weights = _data(ds.variables["wkpt"])
                points = _data(ds.variables["kpt"])
                if (weights.ndim != 1 or points.shape != (len(weights), 3)
                        or np.any(weights < 0) or not np.isclose(weights.sum(), 1., rtol=0, atol=1e-7)):
                    raise ValueError("invalid k-points or unnormalized k weights")
            if variable.dimensions == (e_dim,):
                if not path.name.endswith(".AV.nc"):
                    raise ValueError("k-independent transmission requires a TBtrans .AV.nc source")
                aggregation = "already averaged"
            else:
                if weights is None:
                    raise ValueError("missing k weights for k-resolved transmission")
                k_dim = ds.variables["wkpt"].dimensions[0]
                if variable.dimensions != (k_dim, e_dim):
                    raise ValueError(f"unexpected transmission dimensions {variable.dimensions}")
                data = weights @ data
                aggregation = "weighted once"
            return energy, data, {
                "filename": path.name, "raw_energy_unit": unit,
                "raw_energy_reference": "TranSIESTA equilibrium Fermi energy",
                "variable": "Left/Right.T", "chemical_potentials_relative_eV": potentials,
                "k_aggregation": aggregation,
                "k_weights": None if weights is None else weights.tolist(),
                "k_points": None if points is None else points.tolist(),
            }
    except (OSError, RuntimeError, ValueError, KeyError, IndexError) as exc:
        raise StageError(f"tbtrans_te: {path.name}: {exc}") from exc


def _raw_files(spin):
    if spin == "unpolarized":
        return ["device.TBT.nc"]
    if spin == "collinear":
        return ["device.TBT_UP.nc", "device.TBT_DN.nc"]
    raise StageError(f"tbtrans_te: unsupported spin mode {spin!r}")


def normalize_spectra(paths, spin):
    """Read raw Left -> Right physical channels and apply transverse weights once."""
    expected = _raw_files(spin)
    if len(paths) != len(expected):
        raise StageError(f"tbtrans_te: expected spin files {expected}")
    channels, sources = [], []
    for path, name in zip(paths, expected):
        path = Path(path)
        if path.name not in (name, name[:-3] + ".AV.nc"):
            raise StageError(f"tbtrans_te: expected {name}, got {path.name}")
        energy, channel, meta = _read_channel(path)
        if channels:
            if energy.shape != grid.shape or not np.allclose(energy, grid, rtol=0, atol=1e-9):
                raise StageError("tbtrans_te: spin energy grids differ")
            for key in ("k_points", "k_weights", "k_aggregation"):
                if sources[0][key] != meta[key]:
                    raise StageError(f"tbtrans_te: spin channels differ in {key}")
        else:
            grid = energy
        channels.append(channel)
        sources.append(meta)
    return grid, np.column_stack(channels), sources


def conductance(energy, channels, temperature):
    """Use a sampled/interpolated Fermi value at 0 K or the thermal integral."""
    energy, channels = np.asarray(energy, dtype=float), np.asarray(channels, dtype=float)
    if (not math.isfinite(temperature) or temperature < 0 or energy.ndim != 1 or len(energy) < 2
            or channels.shape not in ((len(energy), 1), (len(energy), 2))
            or not np.all(np.isfinite(energy)) or not np.all(np.isfinite(channels))
            or np.any(channels < 0) or np.any(np.diff(energy) <= 0)):
        raise StageError("junction_analyze: invalid temperature, energy grid or spin transmission")
    if energy[0] > 0 or energy[-1] < 0:
        raise StageError("junction_analyze: energy grid must bracket the common Fermi energy")
    diagnostics = {}
    warnings = []
    if temperature == 0.:
        l0 = [float(np.interp(0., energy, channel)) for channel in channels.T]
        diagnostics["ef_evaluation"] = "sampled" if np.any(energy == 0.) else "linear interpolation"
    else:
        rows = [(None, value) for value in energy.tolist()]
        window = thermal_window(rows, temperature)
        quad = fermi_quadrature(rows, temperature)
        if window["coverage"] < MIN_THERMAL_COVERAGE:
            warnings.append("energy window truncates the Fermi derivative")
        ratio, moment = quad["ratio"], quad["first_moment_kbt"]
        if (ratio is None or abs(ratio - 1.) > 1 - MIN_QUADRATURE_RATIO or
                moment is None or abs(moment) > MAX_FERMI_MOMENT_KBT):
            warnings.append("energy grid does not resolve the Fermi derivative")
        weights = np.array([minus_fermi_derivative(e, temperature) for e in energy])
        l0 = [float(trapz(channel * weights, energy)) for channel in channels.T]
        diagnostics.update(thermal_window=window, fermi_quadrature=quad)
    g = sum(l0) / len(l0)
    return {"temperature_K": temperature, "points": len(energy), "channel_L0": l0,
            "G_over_G0": g, "G_S": g * G0_S, "G0_S": G0_S,
            "spin_convention": ("per-spin, twofold degenerate" if len(l0) == 1
                                else "up and down each counted once"),
            "warnings": warnings, **diagnostics}


class TbtransStage(Stage):
    name = "tbtrans_te"
    depends = ("siesta_electrode", "transiesta_scf")
    script = "sub_siesta.sh"
    outputs = ("transmission.tsv", "spectrum.json")
    remote_files = ("tbtrans.out", "output.sha256")
    input_hint = "TranSIESTA device.TS.HSX and periodic electrode.HSX"

    def remote_files_for(self, ctx):
        manifest = ctx.artifact(self.name, "tbtrans-input.json")
        spin = (read_json(manifest)["spin"] if manifest.is_file()
                else ctx.config.siesta.get("spin", "unpolarized"))
        return [*self.remote_files, *_raw_files(spin)]

    def prepare(self, ctx, stage_dir):
        from .stages import launch, write_script

        scf_dir, electrode_dir, geometry_dir = (
            ctx.stage_dir(stage) for stage in ("transiesta_scf", "siesta_electrode", "siesta_junction_relax"))
        scf = read_json(scf_dir / "scf.json")
        electrode = read_json(electrode_dir / "electrode.json")
        if scf.get("open_boundary_converged") is not True:
            raise StageError("tbtrans_te: missing converged open-boundary SCF")
        verify_files(scf_dir, scf)
        verify_files(electrode_dir, electrode)
        required = {"transiesta.fdf", "device.TS.HSX", "electrode.HSX"}
        if (not required.issubset(scf["files"]) or scf["files"]["electrode.HSX"]
                != electrode["files"].get("electrode.HSX")):
            raise StageError("tbtrans_te: electrode/device inputs differ from the collected SCF record")
        raw = _raw_files(scf.get("spin"))
        geometry = geometry_dir / "final_geometry.json"
        contacts = geometry_dir / "contacts.json"
        if (not geometry.is_file() or scf.get("geometry_sha256") != sha(geometry)
                or not contacts.is_file() or scf.get("contacts_sha256") != sha(contacts)):
            raise StageError("tbtrans_te: missing or mismatched final junction geometry/contacts")
        fermi_ev = scf.get("fermi_ev")
        if not isinstance(fermi_ev, (float, int)) or not math.isfinite(fermi_ev):
            raise StageError("tbtrans_te: missing finite equilibrium Fermi energy")
        smearing = re.search(r"(?im)^ElectronicTemperature\s+(\S+)\s+eV\s*$",
                             (scf_dir / "transiesta.fdf").read_text(encoding="utf-8"))
        if not smearing:
            raise StageError("tbtrans_te: missing DFT occupation smearing in TranSIESTA FDF")
        stage_dir.mkdir(parents=True, exist_ok=True)
        for name in required:
            shutil.copyfile(scf_dir / name, stage_dir / name)
        start, step, end = ctx.config.transport.relative_window()
        if not all(map(math.isfinite, (start, step, end))) or not (end > start and step > 0):
            raise StageError("tbtrans_te: invalid transmission energy window")
        grid = start + np.arange(math.floor((end - start) / step + 1e-9) + 1) * step
        if len(grid) < 2:
            raise StageError("tbtrans_te: energy window requires at least two points")
        # Use explicit points: the default mid-rule samples half a step inside each endpoint.
        grid_path = stage_dir / "energies.dat"
        np.savetxt(grid_path, np.column_stack((grid, np.full(len(grid), step))),
                   fmt="%.14g", header="E-E_F [eV], weight [eV]")
        opts = ctx.config.siesta
        kgrid = opts.get("tbtrans_kgrid", opts.get("device_kgrid", [1, 1, 1]))
        if (not isinstance(kgrid, (list, tuple)) or len(kgrid) != 3 or
                any(type(k) is not int or k <= 0 for k in kgrid) or kgrid[2] != 1):
            raise StageError("tbtrans_te: tbtrans_kgrid needs three positive integers with kz = 1")
        fdf = stage_dir / "tbtrans.fdf"
        fdf.write_text(
            f"TBT.HS device.TS.HSX\nTBT.Voltage 0 eV\nTBT.k [{' '.join(map(str, kgrid))}]\n"
            "TBT.Elecs.Gf.Reuse false\n"
            "%block TBT.Contours\n spectrum\n%endblock\n%block TBT.Contour.spectrum\n"
            f" from {grid[0]:.14g} eV to {grid[-1]:.14g} eV\n file energies.dat\n%endblock\n"
            "%include transiesta.fdf\n", encoding="ascii", newline="\n")
        files = {name: sha(stage_dir / name) for name in sorted(required | {fdf.name, grid_path.name})}
        attempt = uuid.uuid4().hex
        manifest = dump(stage_dir / "tbtrans-input.json", {
            "files": files, "attempt": attempt, "fermi_eV": fermi_ev, "spin": scf["spin"],
            "kgrid": list(kgrid), "smearing_eV": float(smearing[1]),
            "upstream_sha256": {"siesta_electrode/electrode.json": sha(electrode_dir / "electrode.json"),
                                "transiesta_scf/scf.json": sha(scf_dir / "scf.json"),
                                "siesta_junction_relax/final_geometry.json": sha(geometry),
                                "siesta_junction_relax/contacts.json": sha(contacts)},
            "contacts": read_json(contacts)})
        check = stage_dir / "tbtrans-input.sha256"
        check.write_text("".join(f"{digest}  {name}\n" for name, digest in
                                 (files | {manifest.name: sha(manifest)}).items()),
                         encoding="ascii", newline="\n")
        commands = ["sha256sum -c tbtrans-input.sha256",
                    launch(ctx, "tbtrans", "tbtrans.fdf") + " > tbtrans.out 2>&1",
                    "sha256sum tbtrans-input.json tbtrans.out " + " ".join(raw) + " > output.sha256"]
        script = write_script(ctx, stage_dir, self.name, commands,
                              [*self.remote_files, *_raw_files("unpolarized"), *_raw_files("collinear")])
        return self.prepared(ctx, stage_dir, inputs=[*[stage_dir / n for n in files], check, manifest, script],
                             artifacts=[*raw, "tbtrans.out", *self.outputs])

    def collect(self, ctx, stage_dir):
        from .stages import _sha_manifest

        input_path = stage_dir / "tbtrans-input.json"
        inputs = read_json(input_path)
        verify_files(stage_dir, inputs)
        raw = _raw_files(inputs["spin"])
        files = _sha_manifest(stage_dir / "output.sha256", [*raw, "tbtrans.out", input_path.name])
        energy, channels, sources = normalize_spectra([stage_dir / f for f in raw], inputs["spin"])
        expected = np.loadtxt(stage_dir / "energies.dat", ndmin=2)[:, 0]
        if energy.shape != expected.shape or not np.allclose(energy, expected, atol=1e-8, rtol=1e-6):
            raise StageError("tbtrans_te: NetCDF energies differ from requested grid")
        table = stage_dir / "transmission.tsv"
        labels = ["T_per_spin"] if channels.shape[1] == 1 else ["T_up", "T_down"]
        np.savetxt(table, np.column_stack((energy, channels)), delimiter="\t", fmt="%.17g",
                   header="\t".join(["E_minus_EF_eV", *labels]))
        files[table.name] = sha(table)
        files["output.sha256"] = sha(stage_dir / "output.sha256")
        meta = {key: inputs[key] for key in ("fermi_eV", "spin", "contacts", "upstream_sha256",
                                             "smearing_eV", "kgrid")}
        meta.update(files=files, sources=sources, columns=["E_minus_EF_eV", *labels],
                    transmission_unit="dimensionless", electrodes=["Left", "Right"])
        dump(stage_dir / "spectrum.json", meta)
        return self.collected(ctx, stage_dir, artifacts=[*files, "spectrum.json"],
                              results={"points": len(energy), "spin": inputs["spin"]})


class JunctionAnalyzeStage(Stage):
    name = "junction_analyze"
    depends = ("tbtrans_te",)
    needs_submit = False
    outputs = ("results.json", "T_E.png")
    input_hint = "tbtrans_te/transmission.tsv and spectrum.json"

    def prepare(self, ctx, stage_dir):
        source = ctx.override_path("te") or ctx.stage_dir("tbtrans_te")
        directory = source if source.is_dir() else source.parent
        table = directory / "transmission.tsv"
        if not source.is_dir() and source != table:
            raise StageError("junction_analyze: expected TBtrans transmission.tsv")
        metadata = directory / "spectrum.json"
        meta = read_json(metadata)
        verify_files(directory, meta)
        if table.name not in meta["files"]:
            raise StageError("junction_analyze: missing digest for normalized transmission")
        data = np.loadtxt(table, ndmin=2)
        channels = data[:, 1:]
        if channels.shape[1] != len(_raw_files(meta["spin"])):
            raise StageError("junction_analyze: normalized transmission has missing spin channel")
        result = conductance(data[:, 0], channels, ctx.config.transport.temperature)
        result.update({k: meta[k] for k in ("fermi_eV", "smearing_eV", "contacts", "spin", "upstream_sha256")})
        result.update(backend="transiesta", physical_verdict="inconclusive",
                      spectrum=str(metadata), spectrum_sha256=sha(metadata), source_files=meta["files"])
        try:
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_agg import FigureCanvasAgg
        except ImportError as exc:
            raise DependencyMissing("TBtrans plotting requires pip install 'smj-toolkit[siesta]'") from exc
        out = ctx.override_path("out") or stage_dir
        out.mkdir(parents=True, exist_ok=True)
        fig = Figure(figsize=(5.2, 3.5), layout="constrained")
        FigureCanvasAgg(fig)
        ax = fig.subplots()
        for channel, label in zip(channels.T, meta["columns"][1:]):
            ax.plot(data[:, 0], channel, label=label)
        ax.axvline(0., color="0.55", linestyle="--", linewidth=.8)
        ax.set(xlabel=r"$E-E_F$ (eV)", ylabel="Transmission")
        ax.legend(frameon=False)
        fig.savefig(out / "T_E.png", dpi=210)
        dump(out / "results.json", result)
        return self.prepared(ctx, stage_dir, inputs=[table, metadata], local_only=True)

    def collect(self, ctx, stage_dir):
        out = ctx.override_path("out") or stage_dir
        data = read_json(out / "results.json")
        if not (out / "T_E.png").is_file():
            raise StageError("junction_analyze: missing T_E.png")
        return self.collected(ctx, stage_dir, artifacts=list(self.outputs), output_dir=out,
                              results={"junction": data}, notes=data["warnings"])


STAGES = {stage.name: stage for stage in (TbtransStage(), JunctionAnalyzeStage())}
