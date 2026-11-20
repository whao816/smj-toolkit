"""重新生成 ``src/smj/_vendor/junction/VENDOR.json``（本地补丁台账）。

为什么需要它：0.1.0 的 wheel **没有**带 VENDOR.json（README 却说它记录 sha256 与
主仓 commit），而 0.1.1/0.1.2 在 vendored 目录里做了 4 类本地补丁（F2 判据、F1 头回取、
F4 提示文案、F-LIB 库电极与 au_relax），必须留下可追溯的记录，否则下次从主仓 sync 会静默覆盖。

本脚本只用标准库，输出确定性内容（文件按路径排序；不写时间戳到 files 之外时
可选择 --no-timestamp 便于 diff）。用法::

    python scripts/gen_vendor_json.py           # 写入 VENDOR.json
    python scripts/gen_vendor_json.py --check   # 只比对，不一致返回 1
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VENDOR_DIR = ROOT / "src" / "smj" / "_vendor" / "junction"
VENDOR_JSON = VENDOR_DIR / "VENDOR.json"

# 0.1.1 的本地补丁台账：每条写明改了什么、为什么、被哪条测试守住。
LOCAL_PATCHES = [
    {
        "id": "F2-te-present",
        "commit": "bc7febf",
        "files": ["criteria.py", "validate.py"],
        "what": "无 TE.dat 时登记 te_present(engineering/error, missing=True)，"
                "physics_verdict([]) 返回 inconclusive 而非通过；validate_junction_run "
                "接受 temperature 参数。",
        "why": "原实现里空 aitranss_te/ 目录会得到 passed=True / 支持 H0，"
               "缺证据被当成通过（P0 fail-open）。",
        "tests": ["tests/test_criteria_failopen.py"],
    },
    {
        "id": "F1-mos-aims-head",
        "commit": "f7d3970",
        "files": ["artifacts.py"],
        "what": "新增 AIMS_HEAD_STAGES=('aims_scf',)，只有 aims_scf 的 retrieval_plan "
                "才带 mos.aims/omat.aims 的 head 条目。",
        "why": "原先所有阶段都收到 head 条目，gaussian_opt 等阶段凭空多出必然 missing "
               "的项；真正消费 head 的下载逻辑在 smj/remote.py（非 vendored）。",
        "tests": ["tests/test_remote_mos_head.py"],
    },
    {
        "id": "F4-pending-message",
        "commit": "fe48c29",
        "files": ["review.py"],
        "what": "pending_message 不再提示主仓的 `agenticdft smj-review`，改为说明独立包"
                "的真实入口 `smj review --approve --who <名字>`。",
        "why": "该提示指向独立包里不存在的命令（且该检查点模块本身未接线）。",
        "tests": ["tests/test_connection_check_approve.py"],
    },
    {
        "id": "F-LIB-library-electrode",
        "commit": "00294d9",
        "files": ["electrode.py", "inputs.py"],
        "what": "新增 LibraryElectrode（直接吃 FHI-aims electrodes.library 的预定义团簇）、"
                "backplane_indices（按输运轴取背层整层）；append_au60 接受 cluster/provenance "
                "并写 electrode_source / library_cluster / left|right_backplane_indices；"
                "新增 make_au_relax_input（背层整层 constrain_relaxation + bfgs）。",
        "why": "手册 §5.3.1 要求用 electrodes.library 的预定义团簇，并在贴装后固定两个背层"
               "做一次弛豫；原独立包只能自建三角岛、且没有弛豫阶段。",
        "tests": ["tests/test_library_au_relax.py"],
    },
    {
        "id": "F-RELAX-au-relax-stage",
        "commit": "00294d9",
        "files": ["artifacts.py", "review.py"],
        "what": "artifacts 登记 au_relax 的回取清单；review.find_build_dir 优先 "
                "run_dir/au_relax（连接核对（原 A2/表面核对）审计对象改为真正送 SCF 的弛豫后几何）。",
        "why": "0.1.2 新增 au_relax 阶段；不登记会被 download 视为 unknown，"
               "连接核对也会继续审计未弛豫的 au60_build 几何。",
        "tests": ["tests/test_library_au_relax.py", "tests/test_connection_check_approve.py"],
    },
    {
        "id": "F-ELEC-element-lattice",
        "commit": "2026-11-15",
        "files": ["electrode.py", "geometry.py", "inputs.py", "self_energy.py",
                  "validate.py", "a2_viz.py"],
        "what": "电极金属可参数化：add_tip_au(element=…)、make_tip_relax_input(tip_element=…)、"
                "append_au60(element=, a_lattice=) 在非 Au 时写 electrode_element/"
                "electrode_a_lattice、make_self_energy(element=…)、validate/a2_viz 从 metadata "
                "读电极元素。",
        "why": "阶段 4「构建电极」要求 source/library 与电极类型（Au/Ag/Cu/Ni/Pt/Pd/Al 或显式 "
               "a_lattice）都是参数；原实现只认 Au，tip/自能/校验全硬编码。",
        "tests": ["tests/test_electrode_metals.py", "tests/test_no_sulfur_anchor.py"],
    },
    {
        "id": "F-CCHECK-rename",
        "commit": "2026-11-15",
        "files": ["validate.py", "a2_viz.py", "fingerprint.py", "criteria.py", "review.py"],
        "what": "「A2 人工确认」→「表面核对」→「连接核对」更名：check 名 "
                "junction_connection_check；a2_viz 文案与默认文件名改 connection_check_*"
                "（保留数据键 area_A2）；fingerprint 写 connection_check.json、仍读旧 "
                "surface_check.json / a2_confirmation.json；criteria 登记判据族 BPY。",
        "why": "用户可理解的命名 + 默认自动登记（auto_connection_check）+ BPY 不再被 "
               "Au–BDT–Au 标定的判据带 fail-safe 误判。",
        "tests": ["tests/test_connection_check_approve.py", "tests/test_cli_connection_check.py",
                  "tests/test_runner_connection_check.py", "tests/test_review_family_drift.py",
                  "tests/test_family.py"],
    },
    {
        "id": "F-FAMILY-sme",
        "commit": "2026-11-15",
        "files": ["criteria.py"],
        "what": "SYSTEM_FAMILIES 登记判据族 SME（Au-(1,4-bis(methylthio)benzene)-Au），"
                "让 t_ef_magnitude / g_in_theory_band / seebeck_sign 记「不适用（仅记录）」"
                "而不是闸门。",
        "why": "这三条判据带都按 Au–BDT–Au 硫醇盐标定（0.03–0.40 G0、Seebeck 正号）；"
               "甲硫基（硫醚）是配位型接触，S 上无 H、不脱质子，带内数值对它没有参考意义，"
               "不登记族名会把一次物理上正常的算例判成「否证 H0」。",
        "tests": ["tests/test_family.py"],
    },
    {
        "id": "F-FAMILY-alkane",
        "commit": "2026-11-15",
        "files": ["criteria.py"],
        "what": "SYSTEM_FAMILIES 登记判据族 Au-alkane-dithiol-Au（键 HDT/ODT/DDT/C6/C8/C10），"
                "让 t_ef_magnitude / g_in_theory_band / seebeck_sign 记「不适用（仅记录）」"
                "而不是闸门。",
        "why": "烷烃二硫醇走饱和 σ 链的隧穿通道（Xu & Tao 2003 的链长衰减序列），"
               "与 Au–BDT–Au 的 π 通道不同族；三条判据带都按 BDT 硫醇盐标定，"
               "不登记族名会把 C6/C8/C10 正常算例（G 比 BDT 低 2–3 个量级）判成「否证 H0」。",
        "tests": ["tests/test_family.py"],
    },
    {
        "id": "F-CRIT-2026-11-15",
        "commit": "2026-11-15",
        "files": ["criteria.py", "validate.py"],
        "what": "判据口径修订：① t_ef_no_peak_near_ef → t_ef_resonance_attribution"
                "（无峰 pass；峰按 EXPECTED_RESONANCE 族表归因到 HOMO/LUMO 侧，"
                "|ΔE|∈[0.15,1.0] eV 记「已归因」；未登记族按实测 Seebeck 符号自洽回退）；"
                "② 删除 t_ef_magnitude（0.03–0.40 G0 非 Thygesen 2005 原文区间，"
                "派生包络并入 g_in_theory_band，后者 gate error→warn 并标 evidence=derived）；"
                "③ au_s_bond 带 2.40–2.49 → 2.40–2.44 Å（Jackson 2009 S–Au=2.42±0.02 Å；"
                "2.49 无一手源），validate.py 的 AU_S_BAND 与 4 处 detail 同步；"
                "④ te_present 在 TE.dat 不可解析时也记 missing。",
        "why": "原共振判据把「E_F 附近有峰」一律当错误信号，会把 2CN（LUMO 传导）与 AZ"
               "（HOMO 侧）的正常算例判成 inconclusive，BDT 自身 −0.34 eV 的 HOMO 峰也中招；"
               "有源 ≠ 有阈值，派生包络不足以否证 H0。证据见 docs/criteria-source-research.md。",
        "tests": ["tests/test_resonance_attribution.py", "tests/test_family.py",
                  "tests/test_criteria_failopen.py"],
    },
    {
        "id": "F-RECHECK-2026-11-20-window",
        "commit": "2026-11-20",
        "files": ["criteria.py"],
        "what": "新增 physical/error 判据 te_window_covers_efermi（E−E_F 必须跨 0）；"
                "physics_verdict 增加闸门覆盖判定：所有带 applies_to 的族闸门都 "
                "applicable=False 时判 inconclusive（inconclusive_criteria 记 "
                "no_applicable_gate），不再因为「没有闸门可判」而 passed。",
        "why": "复核 P1-2：整窗在 E_F 一侧的合成谱、以及族闸门全部不适用的体系（如薁），"
               "原实现给出 physics_verdict=passed / 支持 H0，属于把「没有证据」包装成"
               "「证据支持」。",
        "tests": ["tests/test_verdict_gates.py"],
    },
    {
        "id": "F-UPSTREAM-PENDING-2026-11-20",
        "commit": "2026-11-20",
        "files": ["fixtures.py", "molecules.py", "surface_index.py", "tcontrol.py",
                  "templates/control_relax.in", "templates/gaussian_opt_freq.gjf",
                  "templates/sub_g16.sh"],
        "what": "副本与主仓工作区一致、但与主仓 HEAD（26c197b）不同：tcontrol 的 "
                "_spin_moment_line、surface_index.pick_surface_triple(away=…)、三份模板的 "
                "SCF 精度/scratch 目录/spin 行、fixtures 的 2026-11-15 判据注记，都是主仓"
                "尚未提交的改动；molecules.py 是等价实现（按 cid 逐个取 conformer）。",
        "why": "副本当年是从主仓**工作区**同步的，主仓这些改动一直没提交，导致 "
               "sync_core --check 每次都报漂移。先登记避免 CI 永久红；主仓提交后应用 "
               "sync_core.py 全量同步并删除本条。",
        "tests": ["scripts/sync_core.py --check"],
    },
    {
        "id": "F-RECHECK-2026-11-20-aims-out-stage",
        "commit": "2026-11-20",
        "files": ["validate.py"],
        "what": "_ARTIFACT_STAGE_HOME['aims.out'] 从 ('aims_scf','aims_tip_relax') 收紧为 "
                "('aims_scf',)：静态 SCF 判据只认 aims_scf 阶段（扁平布局仍由 "
                "run_dir/aims.out 命中），不再回落到弛豫日志。",
        "why": "复核 R13b：弛豫日志每一步都写 Self-consistency cycle converged，"
               "回落到 tip_relax/aims.out 时，一个「静态 SCF 阶段根本没跑」的 run 也会拿到 "
               "aims_scf=OK 判据，报告把弛豫里的某一步 SCF 当成了送 aitranss 的静态 SCF。",
        "tests": ["tests/test_validate_artifacts.py"],
    },
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):  # pragma: no cover
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def build_manifest() -> dict:
    files = {}
    for path in sorted(VENDOR_DIR.rglob("*")):
        if not path.is_file() or path.name == "VENDOR.json":
            continue
        if "__pycache__" in path.parts:
            continue
        rel = path.relative_to(VENDOR_DIR).as_posix()
        files[rel] = {"sha256": _sha256(path), "bytes": path.stat().st_size}
    return {
        "schema": "smj.vendor.v1",
        "upstream": "agenticdft/agenticdft/tools/junction",
        "upstream_commit": None,
        "rebuilt_from": "smj_toolkit-0.1.0-py3-none-any.whl（原 wheel 未带 VENDOR.json）",
        "repo_commit": _git_head(),
        "files": files,
        "local_patches": LOCAL_PATCHES,
        "note": "0.1.1 本地补丁的权威记录；回迁主仓后应删除对应条目并重新 sync。"
                "本文件由 scripts/gen_vendor_json.py 生成，勿手改。"
                "已知偏离（2026-11-20 复核 R12）：vendored criteria.py 是 0.2.1 版注册表，"
                "k_e_zt_order 的元数据仍写 implemented=False（实际评估器会 emit 该判据），"
                "与主仓新注册表不一致；该判据不参与闸门，本轮按要求不改 vendored 判定语义，"
                "仅在此登记，待回迁主仓时一并对齐。",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="只比对现有 VENDOR.json，不一致返回 1（不写文件）")
    args = parser.parse_args(argv)

    manifest = build_manifest()
    text = json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    if args.check:
        if not VENDOR_JSON.is_file():
            print("缺少 VENDOR.json：请运行 python scripts/gen_vendor_json.py")
            return 1
        try:
            current = json.loads(VENDOR_JSON.read_text(encoding="utf-8"))
        except ValueError as exc:
            print("VENDOR.json 不是合法 JSON：%s" % exc)
            return 1
        # repo_commit 记的是"生成时所在提交"，每次 commit 都会变；它不该让
        # 校验变红（2026-11-20 复核 R10：旧实现整份文本比较，CI 必然红）。
        # 其余字段（files 的 sha256/size、local_patches、schema、upstream 等）
        # 逐字段比较，任何真实漂移仍然报错。
        generated = dict(manifest)
        recorded = dict(current)
        generated.pop("repo_commit", None)
        recorded.pop("repo_commit", None)
        if recorded != generated:
            print("VENDOR.json 与磁盘内容不一致：请运行 python scripts/gen_vendor_json.py")
            return 1
        print("VENDOR.json 一致：%d 个文件（repo_commit 不参与比对）" % len(manifest["files"]))
        return 0

    VENDOR_JSON.write_text(text, encoding="utf-8", newline="\n")
    print("写入 %s（%d 个文件，%d 条本地补丁）"
          % (VENDOR_JSON, len(manifest["files"]), len(manifest["local_patches"])))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
