"""随包电极力常数资产（``smj.phonon.assets``）与 ``thermal.lead_fc: auto``。

这套测试守的是"第二个使用者能不能只装包就跑"：以前 `thermal_analyze` 要求
先自己算一次体相 lead FC，否则直接 `StageBlocked`。现在资产随包发布，
``auto`` 是默认值，这里钉住三件事：资产在、sha 对、坏资产会被拒绝。
"""

from __future__ import annotations

import json

import pytest

from smj.phonon import assets


def test_bundled_asset_is_present_and_self_consistent():
    info = assets.lead_fc_info()
    assert info["a_lattice_A"] == pytest.approx(4.07823)
    assert tuple(info["in_plane"]) == (3, 3)
    assert info["layers"] == 2
    assert info["n_atoms_per_layer"] == 9
    # 2 层 × 9 原子 = 18 原子 → 54 自由度
    assert info["sha256"] == info["sha256_actual"]
    import numpy as np

    fc = np.load(info["path"])
    assert fc.shape == (54, 54)


def test_describe_lines_mention_the_key_numbers():
    text = "\n".join(assets.describe_lines())
    assert "Au" in text
    assert "3×3" in text
    assert "sha256" in text


def test_export_round_trips(workdir):
    written = assets.export_lead_fc(workdir / "lead")
    assert len(written) == 2
    info = assets.lead_fc_info(workdir / "lead")
    assert info["sha256_actual"] == assets.lead_fc_info()["sha256"]


def test_corrupted_asset_is_rejected(workdir):
    """sha 校验不是装饰:力常数被截断/被文本转换动过时必须报错。"""
    written = assets.export_lead_fc(workdir / "lead")
    npy = [p for p in written if p.suffix == ".npy"][0]
    payload = bytearray(npy.read_bytes())
    payload[0] ^= 0xFF
    npy.write_bytes(bytes(payload))
    with pytest.raises(ValueError) as excinfo:
        assets.lead_fc_info(workdir / "lead")
    assert "sha256 不符" in str(excinfo.value)


def test_meta_file_without_npy_is_an_error(workdir):
    directory = workdir / "empty"
    directory.mkdir()
    (directory / "lead_fc.json").write_text(json.dumps({"sha256": ""}),
                                            encoding="utf-8")
    with pytest.raises(OSError):
        assets.lead_fc_info(directory)


def test_auto_is_not_treated_as_a_relative_path(workdir):
    """``auto`` 是关键字:不能被 resolve 拼成 yaml 目录下的一个文件。"""
    from smj.config import load_config

    yaml_path = workdir / "smj.yaml"
    yaml_path.write_text(
        "run:\n  name: t\n  mode: dry-run\nthermal:\n  enabled: true\n",
        encoding="utf-8")
    config = load_config(str(yaml_path))
    assert config.thermal.lead_fc == "auto"


def test_explicit_lead_fc_path_is_resolved(workdir):
    from smj.config import load_config

    yaml_path = workdir / "smj.yaml"
    yaml_path.write_text(
        "run:\n  name: t\n  mode: dry-run\n"
        "thermal:\n  enabled: true\n  lead_fc: my_lead\n",
        encoding="utf-8")
    config = load_config(str(yaml_path))
    assert config.thermal.lead_fc.endswith("my_lead")
    # 相对路径按 yaml 所在目录解析:必须是绝对路径，且落在 workdir 下
    from pathlib import Path

    resolved = Path(config.thermal.lead_fc)
    assert resolved.is_absolute()
    assert resolved.parent == Path(workdir).resolve()


def test_stage_resolves_auto_to_the_bundled_asset():
    from smj.stages.thermal_analyze import _lead_fc_path

    path = _lead_fc_path("auto")
    assert path.name == "lead_fc.json"
    assert path.parent == assets.lead_fc_dir()
    assert _lead_fc_path(None) == path
