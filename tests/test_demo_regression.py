"""demo 谱的数值回归（2026-11-20 复核 R9 / 后续①）。

``demo/`` 下的 TE.dat 是仓库里唯一随发布走的真实透射谱。这里对
``demo/demo-reference.json`` 逐个比对：文件哈希（谱有没有被悄悄换掉）、
点数与温度口径、以及 ``analyze_te`` 的积分量（log10 G/G0、S、k_e、ZT）。
参考值由 ``scripts/gen_demo_reference.py`` 生成，CI 另跑一次 ``--check``
守住"参考值文件本身与代码一致"。
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from smj.junction.landauer import analyze_te
from smj.analysis import compare_reference

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "demo"
REFERENCE = DEMO / "demo-reference.json"

#: 明确进 reference 的逐量比对键；其余（点数/温度/哈希）单独断言。
SCALAR_KEYS = ("log10_G_over_G0_like", "S_uV_per_K", "k_e_pW_per_K", "ZT_electronic")


def _reference() -> dict:
    return json.loads(REFERENCE.read_text(encoding="utf-8"))


def _entries() -> list[tuple[str, dict]]:
    return sorted(_reference()["demos"].items())


def _analyze(rel: str, entry: dict) -> dict:
    return analyze_te(DEMO / rel, temperature=float(entry["temperature_K"]))


def test_reference_exists_and_has_tolerances():
    meta = _reference()
    assert meta["schema"] == "smj.demo-reference.v1"
    assert meta["note"]
    for key in SCALAR_KEYS:
        assert meta["tolerance"][key], "缺少 %s 的容差" % key
    assert meta["demos"], "参考值文件里没有任何 demo"


def test_reference_covers_every_demo_spectrum():
    """新增 demo 谱必须同时补参考值，否则这份谱就没有回归保护。"""
    on_disk = {
        p.relative_to(DEMO).as_posix() for p in sorted(DEMO.glob("**/aitranss_te/TE.dat"))
    }
    assert on_disk, "demo/ 下一个 TE.dat 都没有，测试环境不对"
    assert on_disk == set(_reference()["demos"]), "参考值与磁盘上的 demo 谱不一致"


@pytest.mark.parametrize("rel,entry", _entries(), ids=[r.rsplit("/", 3)[0] for r, _ in _entries()])
def test_demo_spectrum_unchanged(rel: str, entry: dict):
    path = DEMO / rel
    assert path.is_file(), "参考值指向的谱不存在: %s" % rel
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == entry["te_sha256"], "谱文件被改动但没有更新参考值: %s" % rel
    result = _analyze(rel, entry)
    assert result["points"] == entry["points"]
    assert result["temperature_K"] == pytest.approx(entry["temperature_K"], abs=1e-9)
    # 谱必须是升序：R1 之后倒序/乱序会被排序并告警，这里顺带守住 demo 自身干净。
    assert result["input_energy_order"] == entry["input_energy_order"] == "ascending"
    assert result["reordered"] is False
    assert result["warnings"] == []


@pytest.mark.parametrize("rel,entry", _entries(), ids=[r.rsplit("/", 3)[0] for r, _ in _entries()])
def test_demo_integrals_within_tolerance(rel: str, entry: dict):
    result = _analyze(rel, entry)
    expected = {key: entry[key] for key in SCALAR_KEYS}
    report = compare_reference(result, expected, tolerance=_reference()["tolerance"])
    assert report["ok"], "%s 的积分量超出容差: %s" % (rel, report["items"])


@pytest.mark.parametrize("rel,entry", _entries(), ids=[r.rsplit("/", 3)[0] for r, _ in _entries()])
def test_demo_max_transmission_within_tolerance(rel: str, entry: dict):
    result = _analyze(rel, entry)
    peak = result["max_transmission"]
    assert peak["T"] == pytest.approx(entry["max_transmission_T"], rel=1e-6, abs=1e-9)
    assert peak["E_minus_Ef_eV"] == pytest.approx(entry["max_transmission_at_eV"], abs=1e-6)


def test_tampered_demo_spectrum_fails_the_comparison(workdir):
    """验收要求：故意改 TE.dat 一个数，回归必须红。

    两点都要成立：
    ① 文件哈希变 → ``test_demo_spectrum_unchanged`` 直接失败（逐字节守住谱）；
    ② 改动落在峰值上时，``max_transmission`` 的 1e-6 级比对也能抓住——
       积分量（S/ZT）对单点改动不敏感，这是物理上应有的，不能指望它们报警。
    """
    rel, entry = _entries()[-1]
    lines = (DEMO / rel).read_text(encoding="utf-8").splitlines(keepends=True)
    data_rows = [i for i, line in enumerate(lines) if line.strip() and not line.lstrip().startswith("#")]
    assert data_rows, "demo 谱没有数据行，测试前提不成立"
    peak = max(data_rows, key=lambda i: float(lines[i].split()[2]))
    columns = lines[peak].split()
    columns[2] = "%.10f" % (float(columns[2]) * 1.10)
    lines[peak] = " ".join(columns) + "\n"

    tampered = workdir / "TE.dat"
    tampered.write_text("".join(lines), encoding="utf-8")

    digest = hashlib.sha256(tampered.read_bytes()).hexdigest()
    assert digest != entry["te_sha256"], "改动后哈希应当变化"

    result = analyze_te(tampered, temperature=float(entry["temperature_K"]))
    peak_after = result["max_transmission"]["T"]
    assert abs(peak_after - entry["max_transmission_T"]) > 1e-6 * abs(entry["max_transmission_T"]), (
        "单个透射值被改 10%%，峰值比对应随之变化，否则 max_transmission 比对形同虚设"
    )

    # 整体口径漂移（例如单位或温度写错）由积分量比对负责抓住。
    scaled = workdir / "TE_scaled.dat"
    scaled_lines = list(lines)
    for i in data_rows:
        cols = scaled_lines[i].split()
        cols[2] = "%.10f" % (float(cols[2]) * 1.5)
        scaled_lines[i] = " ".join(cols) + "\n"
    scaled.write_text("".join(scaled_lines), encoding="utf-8")
    report = compare_reference(
        analyze_te(scaled, temperature=float(entry["temperature_K"])),
        {key: entry[key] for key in SCALAR_KEYS},
        tolerance=_reference()["tolerance"],
    )
    assert report["ok"] is False, "整体透射放大 50% 后积分量仍在容差内，回归样本没有意义"


def test_reference_values_are_actually_reproducible():
    """compare_reference 必须能对"还没写进 reference 的量"报警。

    故意把 S 改大 50%：模型层面这是错的，比对必须红；否则容差比对本身是空转。
    """
    rel, entry = _entries()[0]
    result = _analyze(rel, entry)
    drifted = dict(entry)
    drifted["S_uV_per_K"] = entry["S_uV_per_K"] * 1.5
    report = compare_reference(
        result,
        {"S_uV_per_K": drifted["S_uV_per_K"]},
        tolerance=_reference()["tolerance"],
    )
    assert report["ok"] is False
    item = report["items"][0]
    assert item["ok"] is False
    assert "limit" in item
