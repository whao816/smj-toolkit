"""Producer-side tests with SIESTA 5.4.2 text and HSX v2 record layouts."""
import struct
from types import SimpleNamespace

import numpy as np
import pytest

from smj.errors import ConfigError, StageError


def geom():
    return {"cell": np.eye(3).tolist(), "elements": ["Au"], "positions": [[0, 0, 0]]}


def tiny_hsx(path, *, element="Au", z=0., spin=1):
    from smj.siesta.io import BOHR
    def rec(fmt, *args):
        data = struct.pack("<" + fmt, *args)
        length = struct.pack("<i", len(data))
        return length + data + length
    data = [rec("i", 2), rec("i", 1), rec("7i", 1, 1, spin, 1, 1, 1, 1),
            rec("12d", *(np.eye(3) / BOHR).ravel(), -0.25, 1, 0.01),
            rec("3i3d2i", 0, 0, 0, 0, 0, z / BOHR, 1, 1),
            rec("20sdi", element.encode().ljust(20), 1., 1), rec("3i", 1, 0, 1),
            rec("9i3d", *np.eye(3, dtype=int).ravel(), 0, 0, 0), rec("i", 1),
            rec("i", 1)] + [rec("d", 1.)] * (spin + 1)
    path.write_bytes(b"".join(data))


def test_hsx_rejects_cut_overlap_and_wrong_species(tmp_path):
    from smj.siesta.io import inspect_hsx
    p = tmp_path / "electrode.HSX"
    tiny_hsx(p)
    assert inspect_hsx(p, geom(), "unpolarized")["fermi_ev"] == pytest.approx(-3.4014, abs=1e-4)
    with pytest.raises(StageError, match="uncoupled"):
        inspect_hsx(p, geom(), "unpolarized", electrode=True)
    p.write_bytes(p.read_bytes()[:-16])
    with pytest.raises(StageError, match="HSX"):
        inspect_hsx(p, geom(), "unpolarized")
    tiny_hsx(p, element="Ag")
    with pytest.raises(StageError, match="HSX"):
        inspect_hsx(p, geom(), "unpolarized")


def test_transiesta_convergence_requires_final_open_cycle():
    from smj.siesta.io import convergence
    log = ("Siesta Version  : 5.4.2\nscf:   1  0  0  0  0 -3.1234 0\n"
           "SCF cycle converged after 1 iterations\n")
    end = ">> End of run:   10-JAN-2025\n"
    with pytest.raises(StageError, match="ts-scf"):
        convergence(log + end, open_system=True)
    ts = "ts-scf:   1  0 0 0 0 -3.1234 0\n"
    with pytest.raises(StageError, match="converg"):
        convergence(log + ts + end, open_system=True)
    assert convergence(log + ts + "SCF cycle converged after 1 iterations\n" + end,
                       open_system=True)["scf_converged"]


def test_open_charge_uses_sum_of_both_spin_channels():
    from smj.siesta.io import open_charge
    log = "transiesta: Charge distribution, target = 10.0\nTotal charge [Q] : 5.51 4.5 10.01\n"
    assert open_charge(log, "collinear") == (10., 10.01)
    with pytest.raises(StageError):
        open_charge(log.replace("4.5 10.01", ""), "collinear")


def test_rejects_mismatched_periodic_transverse_meshes():
    from smj.siesta.stages import validate_settings
    with pytest.raises(ConfigError, match="transverse"):
        validate_settings({"electrode_kgrid": [2,3,12], "device_kgrid": [2,2,1]})


def test_stage_script_records_failure_and_validates_scheduler_values(tmp_path):
    from smj.siesta.stages import write_script
    cluster = SimpleNamespace(ntasks=2, partition="gpu1", mem=None, gres="gpu:1",
                              time_limit=None, account=None, qos=None)
    ctx = SimpleNamespace(config=SimpleNamespace(cluster=cluster))
    script = write_script(ctx, tmp_path, "molecule_opt", ["siesta molecule.fdf"], ["molecule.XV"])
    assert ".dsh_done" in script.read_text() and "set -e" in script.read_text()
    cluster.partition = "gpu1\nexit 0"
    with pytest.raises(ConfigError):
        write_script(ctx, tmp_path, "molecule_opt", [], [])


def _context(root):
    from smj.config import MoleculeConfig, ElectrodeConfig, ClusterConfig
    cfg = SimpleNamespace(molecule=MoleculeConfig(), electrode=ElectrodeConfig(strip_thiol_h=False),
                          cluster=ClusterConfig(), siesta={
                              "basis": {"S": "DZP", "Au": "SZ"},
                              "pseudopotentials": {"Au": "Au.psf", "S": "S.psf"},
                              "contact": {"left": {"site": "top"}, "right": {"site": "top"}}})
    for e in ("Au", "S"):
        (root / f"{e}.psf").write_text(e)
    xyz = root / "source.xyz"
    xyz.write_text("2\ninput\nS 0 0 0\nS 0 0 4\n")
    return SimpleNamespace(config=cfg, resolve=lambda p: root / p if p else None,
                           molecule_xyz=lambda: xyz, stage_dir=lambda n: root / n)


def _fixture_output(stage, d, *, mulliken=False):
    """Synthetic completed job for the actual geometry prepared in this directory."""
    from smj.checks import ATOMIC_NUMBERS
    from smj.siesta.io import BOHR, read_json, sha
    g = read_json(d / "geometry.json")
    elements = g["elements"]
    labels = sorted(set(elements))
    ids = [1 + labels.index(e) for e in elements]
    n = len(elements)
    xyz = np.array(g["positions"]) / BOHR
    cell = np.array(g["cell"]) / BOHR
    log = ("Version         : 5.4.2\nscf: 1 0 0 0 0 -3.40 0\n"
           "SCF cycle converged after 1 iterations\n")
    if mulliken:
        rows = "\n".join(f"{i:6d}     0.000000    1.000000  {e}"
                         for i, e in enumerate(elements, 1))
        log += ("Mulliken Atomic Populations:\n"
                "Atom #   charge [q] valence [e]  Species\n" + rows +
                "\n--------------------\n Total     0.000000\n")
    if stage.name == "transiesta_scf":
        log += ("ts-scf: 1 0 0 0 0 -3.40 0\nSCF cycle converged after 1 iterations\n"
                "transiesta: Charge distribution, target = 55.0\nTotal charge [Q] : 55.01\n")
    (d / "siesta.out").write_text(log + " Max 0.001 constrained\n>> End of run: 10-JAN-2025\n")
    for filename in stage.raw:
        p = d / filename
        if filename.endswith(".XV"):
            rows = [" ".join(map(str, row)) + "  0 0 0" for row in cell]
            rows += [str(n)] + [f"{i} {ATOMIC_NUMBERS[e]} " + " ".join(map(str, row)) + " 0 0 0"
                                for i, e, row in zip(ids, elements, xyz)]
            p.write_text("\n".join(rows)+"\n")
        elif filename.endswith(".TSDE"):
            rows = [struct.pack("<5i", n, 1, 1, 1, 1), struct.pack("<i", 1)*n]
            rows += [struct.pack("<i", i) for i in range(1, n+1)] + [struct.pack("<d", 1.)]*(2*n)
            rows += [struct.pack("<d", -.25)]
            p.write_bytes(b"".join(struct.pack("<i", len(r))+r+struct.pack("<i", len(r)) for r in rows))
        elif filename.endswith(".HSX"):
            def pack(fmt, *vals):
                return struct.pack("<" + fmt, *vals)
            translations = [0, 1, -1] if stage.name == "siesta_electrode" else [0]
            nz = len(translations)
            rows = [pack("i", 2), pack("i", 1), pack("7i", n,n,1,len(labels),1,1,nz),
                    pack("12d", *cell.ravel(), -.25, 1., .01),
                    b"".join(pack("3i", 0,0,z) for z in translations)
                    + xyz.astype("<f8").tobytes() + pack(f"{n}i", *ids)
                    + pack(f"{n}i", *range(1, n+1)),
                    b"".join(pack("20sdi", e.encode().ljust(20), 1.,1) for e in labels)]
            rows += [pack("3i", 1,0,1)]*len(labels) + [pack("9i3d", *np.eye(3,dtype=int).ravel(), 0,0,0)]
            rows += [pack("i", nz)*n]
            rows += [pack(f"{nz}i", *(i + j*n for j in range(nz))) for i in range(1,n+1)]
            rows += [pack("d", 1.)*nz]*(2*n)
            p.write_bytes(b"".join(pack("i", len(r))+r+pack("i", len(r)) for r in rows))
    (d / "output.sha256").write_text("".join(f"{sha(d / p)}  {p}\n" for p in ["input.sha256", *stage.raw]))


def test_mulliken_populations_uses_last_complete_table():
    from smj.siesta.io import mulliken_populations
    table = ("Mulliken Atomic Populations:\n"
             "Atom #   charge [q] valence [e]  Species\n"
             "     1     0.100000    4.900000  N\n"
             "     2     0.200000    3.800000  C\n"
             "--------------------\n"
             " Total    -0.000000\n")
    later = table.replace("0.100000", "0.300000").replace("0.200000", "0.400000")
    rows = mulliken_populations("noise\n" + table + "trailer\n" + later)
    assert [row[1] for row in rows] == [0.3, 0.4]
    assert rows[0][2] == pytest.approx(4.9)
    with pytest.raises(StageError, match="Mulliken"):
        mulliken_populations("no table here")


def test_mulliken_charges_require_matching_geometry_order():
    from smj.siesta.stages import _mulliken_charges
    g = {"elements": ["Au", "Au"]}
    pops = [(1, 0.1, 4.9, "Au"), (2, 0.2, 4.8, "Au")]
    assert _mulliken_charges(g, pops) == {1: 0.1, 2: 0.2}
    with pytest.raises(StageError, match="order"):
        _mulliken_charges(g, pops[:1])
    with pytest.raises(StageError, match="species"):
        _mulliken_charges(g, [(1, 0.1, 4.9, "Ag"), pops[1]])


def test_layer_screening_orders_both_sides_from_the_lead_surface():
    from smj.siesta.stages import _layer_screening
    geom = {"positions": [[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 0, 3]],
            "contacts": {"left": {"layer_ids": [[1], [2]]},
                         "right": {"layer_ids": [[3], [4]]}}}
    charges = {1: 0.02, 2: 0.03, 3: 0.04, 4: 0.01}
    report = _layer_screening(geom, charges, 0.0, 0.0)
    assert report["status"] == "passed"
    assert report["sides"]["left"]["outermost_delta_e"] == pytest.approx(0.02)
    assert report["sides"]["right"]["outermost_delta_e"] == pytest.approx(0.01)
    failed = _layer_screening(geom, {1: 0.2, 2: 0.03, 3: 0.04, 4: 0.01}, 0.0, 0.0)
    assert failed["status"] == "failed"
    missing = _layer_screening({"positions": [[0, 0, 0]]}, {1: 0.0}, 0.0, 0.0)
    assert missing["status"] == "unverified" and "layer map" in missing["reason"]


def test_electrode_collect_extracts_mulliken_bulk_reference(tmp_path):
    from smj.siesta.stages import STAGES
    ctx = _context(tmp_path)
    stage = STAGES["siesta_electrode"]
    d = ctx.stage_dir(stage.name)
    stage.prepare(ctx, d)
    _fixture_output(stage, d, mulliken=True)
    result = stage.collect(ctx, d)
    assert result.results["mulliken_bulk_charge_e"] == pytest.approx(0.0)
    assert result.results["mulliken_layer_spread_e"] == pytest.approx(0.0)


def test_junction_relax_reuses_a_collected_geometry(tmp_path):
    import shutil
    from smj.siesta.stages import STAGES
    ctx = _context(tmp_path)
    for name in ("molecule_opt", "siesta_electrode", "siesta_junction_relax"):
        stage = STAGES[name]
        d = ctx.stage_dir(name)
        stage.prepare(ctx, d)
        _fixture_output(stage, d)
        stage.collect(ctx, d)
    source = tmp_path / "siesta_junction_relax"
    root = tmp_path / "reuse"
    for name in ("molecule_opt", "siesta_electrode"):
        shutil.copytree(tmp_path / name, root / name)
    reuse = _context(root)
    reuse.config.siesta["reuse_relax_from"] = str(source)
    stage = STAGES["siesta_junction_relax"]
    d = reuse.stage_dir("siesta_junction_relax")
    prepared = stage.prepare(reuse, d)
    assert prepared.local_only
    result = stage.collect(reuse, d)
    assert result.results["reused_from"] == str(source)
    assert (d / "final_geometry.json").is_file()


def test_transiesta_scf_reuses_a_collected_device(tmp_path):
    import shutil
    from smj.siesta.stages import STAGES
    from smj.siesta.io import read_json
    ctx = _context(tmp_path)
    for name in ("molecule_opt", "siesta_electrode", "siesta_junction_relax", "transiesta_scf"):
        stage = STAGES[name]
        d = ctx.stage_dir(name)
        stage.prepare(ctx, d)
        _fixture_output(stage, d)
        stage.collect(ctx, d)
    source = tmp_path / "transiesta_scf"
    root = tmp_path / "reuse_scf"
    for name in ("molecule_opt", "siesta_electrode", "siesta_junction_relax"):
        shutil.copytree(tmp_path / name, root / name)
    reuse = _context(root)
    reuse.config.siesta["reuse_scf_from"] = str(source)
    stage = STAGES["transiesta_scf"]
    d = reuse.stage_dir("transiesta_scf")
    prepared = stage.prepare(reuse, d)
    assert prepared.local_only
    result = stage.collect(reuse, d)
    assert result.results["reused_from"] == str(source)
    assert read_json(d / "scf.json")["open_boundary_converged"] is True
    assert (d / "device.TS.HSX").is_file()
    assert (d / "reuse.json").is_file()

    reuse.config.siesta["reuse_scf_from"] = str(tmp_path / "missing")
    with pytest.raises(StageError, match="does not contain"):
        stage.collect(reuse, d)


def test_complete_producers_reuse_and_receipts(tmp_path):
    from smj.siesta.stages import STAGES
    from smj.siesta.io import read_json
    ctx = _context(tmp_path)
    for name, stage in STAGES.items():
        d = ctx.stage_dir(name)
        prepared = stage.prepare(ctx, d)
        assert prepared.script and all(p.is_file() for p in prepared.inputs)
        _fixture_output(stage, d)
        stage.collect(ctx, d)
    scf = ctx.stage_dir("transiesta_scf")
    assert read_json(scf / "scf.json")["charge_deviation_e"] == pytest.approx(.01)
    assert "elec-pos 1" in (scf / "transiesta.fdf").read_text()
    assert "DE electrode.TSDE" in (scf / "transiesta.fdf").read_text()
    assert "contour.eq.pole " not in (scf / "transiesta.fdf").read_text().lower()
    assert (scf / "electrode.TSDE").is_file()
    # A new attempt may not claim the local output of the previous one.
    stage.prepare(ctx, scf)
    with pytest.raises(StageError, match="sha256"):
        stage.collect(ctx, scf)
    lead = STAGES["siesta_electrode"]
    ctx.config.siesta["electrode_reuse"] = "siesta_electrode"
    reused = tmp_path / "reused"
    assert lead.prepare(ctx, reused).local_only
    assert lead.collect(ctx, reused).results["source"] == str(ctx.stage_dir("siesta_electrode"))
    ctx.config.siesta["mesh_cutoff_ry"] = 501
    with pytest.raises(StageError, match="incompatible"):
        lead.prepare(ctx, reused)


def test_missing_restart_data_rejected_even_with_valid_receipt(tmp_path):
    from smj.siesta.stages import STAGES
    from smj.siesta.io import sha
    ctx = _context(tmp_path)
    stage = STAGES["siesta_electrode"]
    d = ctx.stage_dir(stage.name)
    stage.prepare(ctx, d)
    _fixture_output(stage, d)
    restart = d / "electrode.TSDE"
    restart.write_bytes(restart.read_bytes()[:-16])
    (d / "output.sha256").write_text("".join(f"{sha(d / p)}  {p}\n"
                                             for p in ["input.sha256", *stage.raw]))
    with pytest.raises(StageError, match="TSDE"):
        stage.collect(ctx, d)
