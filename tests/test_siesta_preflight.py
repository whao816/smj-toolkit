import pytest

from smj.config import load_config
from smj.preflight import run_preflight
from smj.transport import CommandResult


def make_config(workdir):
    cfg = workdir / "smj.yaml"
    cfg.write_text("run: {name: au-test, mode: dry-run}\n"
                   "transport: {backend: transiesta}\n"
                   "molecule: {xyz: molecule.xyz}\n"
                   "siesta: {contact: {left: {site: top}, right: {motif: adatom}}}\n",
                   encoding="utf-8")
    (workdir / "molecule.xyz").write_text("2\nfixture\nS 0 0 0\nS 0 0 6\n")
    return load_config(cfg)


def test_route_checks_actual_geometry_and_pseudopotentials(workdir):
    cfg = make_config(workdir)
    result = run_preflight(cfg, check_env=False)
    assert "siesta_geometry_ok" in result.codes()
    assert "siesta_pseudos_invalid" in result.codes()
    assert "species_dir_missing" not in result.codes()
    cfg.siesta["contact"]["right"]["site"] = "hcp"  # adatom has only an apex
    assert "siesta_geometry_invalid" in run_preflight(cfg, check_env=False).codes()
    del cfg.siesta["contact"]
    assert "siesta_geometry_ok" in run_preflight(cfg, check_env=False).codes()


def test_electrode_preflight_does_not_require_molecule(workdir):
    cfg = make_config(workdir)
    (workdir / "molecule.xyz").unlink()
    result = run_preflight(cfg, stage="siesta_electrode", check_env=False)
    assert "molecule_not_found" not in result.codes()
    assert "siesta_geometry_ok" in result.codes()


@pytest.mark.parametrize("version,cdf,ok", [
    ("5.4.2", "NetCDF-4 support", True),
    ("5.4.2", "NetCDF support", False),
    ("5.4.22", "NetCDF-4 support", False),
])
def test_read_only_cluster_version_and_build_flags(workdir, version, cdf, ok):
    cfg = make_config(workdir)
    cfg.run.mode = "slurm-ssh"
    cfg.cluster.host = "gpu1"
    cfg.cluster.partition = "gpu1"
    class Fake:
        def run(self, command, **kw):
            if "realpath -m" in command:
                return CommandResult(0, "/home/test/smj-runs/au-test\n")
            if "siesta -v" in command or "tbtrans -v" in command:
                return CommandResult(0, f"Version         : {version}\n{cdf}\n")
            return CommandResult(0)
    result = run_preflight(cfg, transport=Fake())
    assert result.remote_checked
    assert ("siesta_binaries_ok" in result.codes()) is ok
    assert ("siesta_binary_invalid" in result.codes()) is not ok
