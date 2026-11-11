"""F1 回归：mos.aims 只取表头（head），且回取后 aitranss_te 能真的算出 nsaos。

背景（2026-11-11 本地补丁 F1）：``artifacts.retrieval_plan`` 早就把 ``mos.aims`` 归进
``head``（只取前 2 行），但 ``RemoteSession.download`` 只消费 ``must``，``head`` 从未
被读取。结果是默认 slurm-ssh 链路在 aims_scf 之后本地没有 mos.aims，
``aitranss_te.prepare`` 在 ``parse_mos_nsaos`` 处 FileNotFoundError -> StageBlocked，
第 5 阶段必断（P0）。
"""
from __future__ import annotations

import json
import shutil
import uuid
from pathlib import Path

import pytest

from smj._vendor.junction import artifacts
from smj._vendor.junction.inputs import make_aitranss_input
from smj._vendor.junction.tcontrol import parse_mos_nsaos
from smj.remote import RemoteSession
from smj.testing import FakeTransport

REMOTE = "/remote/runs/r1"


@pytest.fixture()
def workdir() -> Path:
    base = Path(".tmp") / ("f1-" + uuid.uuid4().hex[:8])
    base.mkdir(parents=True)
    try:
        yield base
    finally:
        shutil.rmtree(base, ignore_errors=True)


class _Ctx:
    def __init__(self, transport: FakeTransport) -> None:
        self.transport = transport

    def remote_stage_dir(self, stage: str) -> str:
        # aitranss_te 复用 aims_scf 的远端目录（与真实链路一致）
        return "%s/%s" % (REMOTE, "aims_scf" if stage == "aitranss_te" else stage)


def _seed_aims_scf(transport: FakeTransport) -> None:
    d = "%s/aims_scf" % REMOTE
    for name in ("geometry.in", "control.in", "aims.out", "basis-indices.out",
                 "sub_aims.sh", ".dsh_done"):
        transport.add_text("%s/%s" % (d, name), "placeholder\n")
    # 0.9 GB 级的真实 mos.aims：表头含 nsaos，后面是大量轨道系数
    transport.add_text("%s/mos.aims" % d,
                       "nsaos = 4210\n" + ("0.123456 0.234567 0.345678\n" * 20000))
    transport.add_text("%s/omat.aims" % d, "omat header\n" + ("1.0\n" * 20000))


def test_plan_head_only_for_aims_stages():
    scf = artifacts.retrieval_plan("aims_scf")
    assert "mos.aims" in scf["head"]
    assert "mos.aims" not in scf["must"]
    # gaussian_opt 不产 mos.aims，不该凭空多出必然 missing 的 head 条目
    assert artifacts.retrieval_plan("gaussian_opt")["head"] == {}


def test_download_uses_head_for_mos_aims(workdir: Path):
    transport = FakeTransport()
    _seed_aims_scf(transport)
    session = RemoteSession(_Ctx(transport), transport=transport)
    out = session.download("aims_scf", workdir / "scf")

    assert out["missing"] == []
    assert "mos.aims" in out["head_fetched"]
    local = (workdir / "scf" / "mos.aims").read_text(encoding="utf-8")
    assert local.startswith("nsaos = 4210")
    assert len(local) < 200, "head 模式必须只取表头，不得整份拉回"
    assert (workdir / "scf" / "aims.out").is_file()
    # 关键：head 下来的表头足以解析 nsaos
    assert parse_mos_nsaos(workdir / "scf" / "mos.aims") == 4210


def test_missing_must_is_reported_not_silently_collected(workdir: Path):
    transport = FakeTransport()
    d = "%s/aims_scf" % REMOTE
    transport.add_text("%s/geometry.in" % d, "placeholder\n")   # 其余 must 故意缺席
    session = RemoteSession(_Ctx(transport), transport=transport)
    out = session.download("aims_scf", workdir / "scf")
    assert "aims.out" in out["missing"]
    assert "mos.aims" not in out["head_fetched"]


def test_aitranss_input_builds_after_head_retrieval(workdir: Path):
    """P0 闭环：head 回取 -> aitranss_te.prepare 的 make_tcontrol 不再 FileNotFoundError。"""
    transport = FakeTransport()
    _seed_aims_scf(transport)
    session = RemoteSession(_Ctx(transport), transport=transport)
    scf = workdir / "scf"
    session.download("aims_scf", scf)

    # 把下载到的占位文件换成 make_tcontrol 真正会读的内容
    atoms = "\n".join("atom %.3f %.3f %.3f Au" % (1.0 * i, 0.0, 0.0) for i in range(8))
    (scf / "geometry.in").write_text(atoms + "\n", encoding="utf-8")
    (scf / "aims.out").write_text("Fermi level: -0.2000 Ha\n", encoding="utf-8")
    (scf / "mos.aims").write_text("nsaos = 4210\n", encoding="utf-8")
    (scf / "electrode_metadata.json").write_text(json.dumps(
        {"left_surface_indices": [1, 2, 3], "right_surface_indices": [4, 5, 6]}),
        encoding="utf-8")

    result = make_aitranss_input(
        scf, workdir / "te", metadata=scf / "electrode_metadata.json")
    tcontrol = Path(result["tcontrol"]).read_text(encoding="utf-8")
    assert "nsaos" in tcontrol.lower()
    assert "4210" in tcontrol
    assert result["values"]["nsaos"] == 4210
