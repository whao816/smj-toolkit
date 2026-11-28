"""Offline TBtrans NetCDF and zero-bias conductance regressions."""
import numpy as np
import pytest

from smj.errors import StageError
from smj.siesta.spectrum import conductance, normalize_spectra


def netcdf(path, energy, transmission, *, weights=(0.25, 0.75), unit="Ry",
           mu=(0., 0.), k_shift=0.):
    nc = pytest.importorskip("netCDF4")
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("ne", len(energy))
        ds.createDimension("nkpt", len(weights))
        ds.createDimension("xyz", 3)
        e = ds.createVariable("E", "f8", ("ne",))
        e.unit = unit
        e[:] = energy
        ds.createVariable("wkpt", "f8", ("nkpt",))[:] = weights
        points = np.zeros((len(weights), 3))
        points[:, 0] = np.arange(len(weights)) / len(weights) + k_shift
        ds.createVariable("kpt", "f8", ("nkpt", "xyz"))[:] = points
        left = ds.createGroup("Left")
        right = ds.createGroup("Right")
        for grp, chemical in zip((left, right), mu):
            pot = grp.createVariable("mu", "f8", ())
            pot.unit = unit
            pot[...] = chemical
        axes = ("nkpt", "ne") if np.ndim(transmission) == 2 else ("ne",)
        left.createVariable("Right.T", "f8", axes)[:] = transmission


def test_asymmetric_peak_ry_relative_and_weighted_once(tmp_path):
    path = tmp_path / "device.TBT.nc"
    netcdf(path, [-0.1, 0., 0.1], [[0., 0., 1.], [0., 0.8, 0.]])
    energy, channels, meta = normalize_spectra([path], "unpolarized")
    np.testing.assert_allclose(energy, [-1.36057, 0., 1.36057], rtol=2e-5)
    np.testing.assert_allclose(channels[:, 0], [0., 0.6, 0.25])
    assert meta[0]["k_weights"] == [0.25, 0.75]


def test_preaveraged_spectrum_is_not_weighted_again(tmp_path):
    path = tmp_path / "device.TBT.AV.nc"
    netcdf(path, [-1., 0., 1.], [0.2, 0.8, 0.4], unit="eV")
    energy, channels, meta = normalize_spectra([path], "unpolarized")
    assert energy.tolist() == [-1., 0., 1.]
    assert channels[:, 0].tolist() == [0.2, 0.8, 0.4]
    assert meta[0]["k_aggregation"] == "already averaged"


def test_invalid_negative_transmission_at_one_k_is_detected_before_averaging(tmp_path):
    path = tmp_path / "device.TBT.nc"
    netcdf(path, [-0.1, 0., 0.1], [[-1., 0., 1.], [1., 0.8, 0.]])
    with pytest.raises(StageError, match="negative"):
        normalize_spectra([path], "unpolarized")


@pytest.mark.parametrize("error", ["weights", "kgrid", "energy", "missing", "bias"])
def test_rejects_unmatched_spin_sources(tmp_path, error):
    up = tmp_path / "device.TBT_UP.nc"
    down = tmp_path / "device.TBT_DN.nc"
    netcdf(up, [-.1, .1], [[1., 1.], [1., 1.]])
    if error != "missing":
        netcdf(down, [-.1, .2 if error == "energy" else .1], np.zeros((2, 2)),
               weights=(.5, .5) if error == "weights" else (.25, .75),
               k_shift=.1 if error == "kgrid" else 0.,
               mu=(.1, 0.) if error == "bias" else (0., 0.))
    with pytest.raises(StageError):
        normalize_spectra([up, down], "collinear")


@pytest.mark.parametrize("temperature", [0., 300.])
@pytest.mark.parametrize(("trans", "expected"), [([1.], 1.), ([1., 0.], .5)])
def test_conductance_physical_units_and_spin(temperature, trans, expected):
    e = np.linspace(-1., 1., 801)
    result = conductance(e, np.tile(trans, (len(e), 1)), temperature)
    assert result["G_over_G0"] == pytest.approx(expected)
    assert result["G_S"] == pytest.approx(7.74809172986e-5 * expected)


def test_zero_temperature_interpolation_and_extrapolation():
    result = conductance(np.array([-1., .5]), np.array([[0.], [.9]]), 0.)
    assert result["G_over_G0"] == pytest.approx(.6)
    assert "interpol" in result["ef_evaluation"]
    with pytest.raises(StageError):
        conductance(np.array([.1, .5]), np.ones((2, 1)), 0.)


def test_underresolved_temperature_window_is_flagged():
    result = conductance(np.array([-.1, 0., .1]), np.ones((3, 1)), 10.)
    assert result["fermi_quadrature"]["ratio"] > 1.05
    assert result["warnings"]


def test_stages_preserve_raw_and_contact_provenance_and_refuse_stale_outputs(tmp_path):
    pytest.importorskip("matplotlib")
    from smj.config import RunConfig, RunInfo
    from smj.context import StageContext
    from smj.siesta.io import dump, read_json, sha
    from smj.siesta.spectrum import STAGES
    from smj.state import RunLayout

    config = RunConfig(run=RunInfo("offline", mode="local"))
    config.transport.backend = "transiesta"
    config.transport.temperature = 0
    config.transport.energy_window = (-1., 1., 1.)
    ctx = StageContext(config, RunLayout(tmp_path), None)
    geom, elec, scf, tbt, out = [
        ctx.ensure_stage_dir(s) for s in (
            "siesta_junction_relax", "siesta_electrode", "transiesta_scf",
            "tbtrans_te", "junction_analyze")]
    geometry = dump(geom / "final_geometry.json", {"positions": [[0., 0., 0.]]})
    contacts = {"left": {"motif": "flat"}, "right": {"motif": "adatom"}}
    dump(geom / "contacts.json", contacts)
    (elec / "electrode.HSX").write_bytes(b"electrode hamiltonian")
    dump(elec / "electrode.json", {"files": {"electrode.HSX": sha(elec / "electrode.HSX")}})
    (scf / "device.TS.HSX").write_bytes(b"open-system hamiltonian")
    (scf / "electrode.HSX").write_bytes((elec / "electrode.HSX").read_bytes())
    (scf / "transiesta.fdf").write_text("SystemLabel device\nElectronicTemperature 0.013 eV\n",
                                        encoding="ascii")
    dump(scf / "scf.json", {
        "files": {f.name: sha(f) for f in scf.iterdir()}, "fermi_ev": -4.7,
        "spin": "unpolarized", "open_boundary_converged": True, "geometry_sha256": sha(geometry),
        "contacts_sha256": sha(geom / "contacts.json")})
    stage = STAGES["tbtrans_te"]
    stage.prepare(ctx, tbt)
    assert "TBT.HS device.TS.HSX" in (tbt / "tbtrans.fdf").read_text()
    assert "device.TBT.nc" in stage.remote_files_for(ctx)
    netcdf(tbt / "device.TBT.nc", [-1., 0., 1.], [[.1, .2, .6], [.2, .8, .3]], unit="eV")
    (tbt / "tbtrans.out").write_text("TBtrans\n")
    raw = ["device.TBT.nc", "tbtrans.out", "tbtrans-input.json"]
    (tbt / "output.sha256").write_text("".join(f"{sha(tbt / f)}  {f}\n" for f in raw))
    stage.collect(ctx, tbt)
    analyze = STAGES["junction_analyze"]
    analyze.prepare(ctx, out)
    result = analyze.collect(ctx, out).results["junction"]
    assert result["G_over_G0"] == pytest.approx(.65)
    assert result["contacts"] == contacts
    assert result["fermi_eV"] == -4.7
    assert result["smearing_eV"] == .013
    assert result["physical_verdict"] == "inconclusive"
    assert (out / "T_E.png").read_bytes().startswith(b"\x89PNG")
    (tbt / "transmission.tsv").write_text("0 1\n1 1\n")
    with pytest.raises(StageError):
        analyze.prepare(ctx, out)
    stage.prepare(ctx, tbt)
    with pytest.raises(StageError):
        stage.collect(ctx, tbt)
    dump(geom / "contacts.json", {})
    with pytest.raises(StageError):
        stage.prepare(ctx, tbt)
