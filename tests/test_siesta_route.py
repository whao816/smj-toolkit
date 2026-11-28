import json

import pytest

from smj.cli import main
from smj.config import load_config
from smj.errors import SmjError
from smj.runner import StageRunner


def config_path(tmp_path):
    path = tmp_path / "smj.yaml"
    path.write_text(
        "run: {name: siesta-test, dir: runs, mode: dry-run}\n"
        "transport: {backend: transiesta, temperature: 0}\n"
        "siesta: {contact: {left: {site: fcc}, right: {site: fcc}}}\n",
        encoding="utf-8",
    )
    return path


def test_plan_and_status_follow_siesta_route(tmp_path, capsys):
    config = config_path(tmp_path)
    expected = ["molecule_opt", "siesta_electrode", "siesta_junction_relax",
                "transiesta_scf", "tbtrans_te", "junction_analyze"]
    for command in ("plan", "status"):
        assert main([command, "-c", str(config), "--json"]) == 0
        result = json.loads(capsys.readouterr().out)
        assert [s["stage"] for s in result["stages"]] == expected


def test_route_bound_to_run_state(tmp_path):
    config = load_config(config_path(tmp_path))
    runner = StageRunner(config)
    runner.state.set("siesta_electrode", status="prepared")
    runner.state.save()
    assert StageRunner(config).state.status("siesta_electrode") == "prepared"
    with pytest.raises(SmjError):
        runner.get_stage("aims_scf")
    config.transport.backend = "aitranss"
    with pytest.raises(SmjError):
        StageRunner(config)

def test_six_stage_runner_and_changed_contact(tmp_path, capsys):
    from test_siesta_stages import _context, _fixture_output
    from test_siesta_spectrum import netcdf
    from smj.siesta.io import sha
    import yaml

    fixture = _context(tmp_path)
    path = config_path(tmp_path)
    data = yaml.safe_load(path.read_text())
    data.update(siesta=fixture.config.siesta, molecule={"xyz": "source.xyz"},
                electrode={"strip_thiol_h": False})
    data["transport"]["energy_window"] = [-1, 1, 1]
    path.write_text(yaml.safe_dump(data))
    runner = StageRunner(load_config(path))
    for name in runner.state.order:
        runner.prepare(name)
        d = runner.layout.stage_dir(name)
        if name == "tbtrans_te":
            netcdf(d / "device.TBT.nc", [-1., 0., 1.], [[.2]*3, [.6]*3], unit="eV")
            (d / "tbtrans.out").write_text("TBtrans")
            files = ["device.TBT.nc", "tbtrans.out", "tbtrans-input.json"]
            (d / "output.sha256").write_text("".join(f"{sha(d / f)}  {f}\n" for f in files))
        elif name != "junction_analyze":
            _fixture_output(runner.get_stage(name), d)
        runner.collect(name)
    assert all(r.ok for r in runner.run())
    assert main(["analyze", "-c", str(path), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["G_over_G0"] == pytest.approx(.5)
    from smj.report import build_manifest
    manifest = build_manifest(runner.config, include_checks=True)
    assert manifest["verdict"]["engineering_passed"], manifest["stage_issues"]
    assert manifest["transport"]["G_over_G0"] == pytest.approx(.5)
    runner.config.siesta["contact"]["right"] = {"motif": "adatom"}
    assert build_manifest(runner.config)["transport"] is None
    with pytest.raises(SmjError, match="siesta_junction_relax"):
        runner.run("junction_analyze", "junction_analyze")
    runner.prepare("siesta_junction_relax", force=True)
    assert runner.state.status("siesta_electrode") == "collected"
    assert runner.state.status("junction_analyze") == "pending"
