"""单分子结 aitranss 自能文件（self.energy.in）生成器。

用途: 为 aitranss 透射计算生成 `self.energy.in` —— 逐原子列出虚部自能
      (imaginary self-energy), 电极原子 (left/right) 按离 Au tip 距离分层:
      近 3 个 (z=1 层) = 0.05, 远 6 个 (z=2 层) = 0.1; 分子原子与 Au tip
      为 empty (0.0)。

开发背景 (AGD-2026-09-22-SMJ, P4d/P5): 该生成器是本会话为 2CN-656 大分子
流程开发的新工具, 来源工具包 (tmp/自动化 qmtransport) 无此生成器。

**left/right 纪律 (重要)**:
- 左右电极划分**必须以 `electrode_metadata.json` 的
  `left_tip_index / right_tip_index` 为准** (append_au60 的输出, 与 tcontrol
  表面索引一致), 不得用"距 tip 距离"自动猜测 —— 2CN 流程实测按距离判定
  会把 left/right 反置, 导致电极虚部加到错误一侧。
- 生成后应校验: `left_surface_indices` 全部落在 left 电极、`right_surface_indices`
  全部落在 right 电极 (见 tests/test_self_energy.py)。
- 输出数值用 Fortran (D20.14) 风格 (首数 0.ddd × 10^e), 与 aitranss 参考
  self.energy.in 一致。

纯本地工具: 只做输入文件生成, 不做真实计算。
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from .geometry import read_aims_geometry

__all__ = ["make_self_energy", "fortran_d"]


def fortran_d(v: float) -> str:
    """Fortran (D20.14) 风格: 首数 0.ddd × 10^e, 指数 2 位带符号。

    例: 0.05 -> "0.50000000000000D-01", 0.0 -> "0.00000000000000D+00"。
    """
    if v == 0.0:
        return "0.00000000000000D+00"
    e = math.floor(math.log10(abs(v)))
    m = v / 10.0 ** e
    if abs(m) >= 10.0:
        m /= 10.0
        e += 1
    m /= 10.0
    e += 1
    sign = "+" if e >= 0 else "-"
    return "%.14fD%s%02d" % (m, sign, abs(e))


DEFAULT_SIGMA_NEAR = 0.05
DEFAULT_SIGMA_FAR = 0.10


def make_self_energy(
    geometry_in: str | Path,
    output: str | Path,
    au_per_electrode: int = 10,
    metadata: str | Path | None = None,
    sigma_near: float = DEFAULT_SIGMA_NEAR,
    sigma_far: float = DEFAULT_SIGMA_FAR,
) -> dict:
    """由结 geometry.in 生成 self.energy.in。

    Args:
        geometry_in: append_au60 之后的完整结 geometry.in (分子 + Au tip + 电极)。
        output: 写出路径。
        au_per_electrode: 每侧电极原子数 (含 tip), 默认 10 → 每侧 9 个电极 Au。
        metadata: 可选 electrode_metadata.json。提供时以其中的
            left_tip_index/right_tip_index 划分左右 (推荐, 保证与 tcontrol
            表面索引一致); 缺省时按"距最近非 Au 原子距离"猜测 (仅诊断用)。
        sigma_near: 距 tip 最近 3 个电极原子的虚部自能 (Ha), 默认 0.05。
        sigma_far: 其余 6 个电极原子的虚部自能 (Ha), 默认 0.10。

    `sigma_near/sigma_far` 是 **模型常数**, 不是收敛参数: 生产默认值必须原样使用,
    否则结果不再与本仓冻结的 BDT/C4/2CN 基准可比 (molecular-junction-plan §5.2.6)。
    它们只供 dev plan §5.4.4 的 σ 判决实验 (W0-6, `scripts/smj_sigma_scan.py`) 扫描
    ——该实验的全部意义就是改动这对常数看结论怎么动。缺省值下输出**逐字节不变**。

    Returns:
        摘要 dict: {output, atom_count, left_tip, right_tip,
                    left_electrode, right_electrode}。

    Raises:
        ValueError: Au 原子不足 2 个, 或电极数 ≠ 2*(au_per_electrode-1),
            或左右划分失败。
    """
    geom = Path(geometry_in)
    out = Path(output)
    atoms = read_aims_geometry(geom)
    au_indices = [i for i, a in enumerate(atoms) if a.element == "Au"]
    if len(au_indices) < 2:
        raise ValueError(f"Expected at least two Au atoms, found {len(au_indices)}")

    def _d2(i: int, j: int) -> float:
        a, b = atoms[i], atoms[j]
        return (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2

    data: dict | None = None
    if metadata:
        data = json.loads(Path(metadata).read_text(encoding="utf-8"))
        left_tip = int(data["left_tip_index"]) - 1
        right_tip = int(data["right_tip_index"]) - 1
    else:
        # 诊断回退: tip = 与最近非 Au 原子距离最小的两个 Au (顺序不可靠)。
        def _min_nonau(i: int) -> float:
            return min(_d2(i, j) for j, b in enumerate(atoms) if b.element != "Au")

        ranked = sorted(au_indices, key=_min_nonau)
        left_tip, right_tip = ranked[0], ranked[1]

    electrode = [i for i in au_indices if i not in (left_tip, right_tip)]
    if len(electrode) != 2 * (au_per_electrode - 1):
        raise ValueError(
            f"Expected {2 * (au_per_electrode - 1)} electrode Au, found {len(electrode)}"
        )

    # W0-2: 厚电极的贴装轴随分子倾斜, 远层原子可能离对侧 tip 更近 —— 按距离
    # 划分会把它们分错 (260 原子结实测左 122 / 右 124)。metadata 若给出精确
    # 归属 (append_au60 写入的 left/right_electrode_indices, 1-based, 含 tip),
    # 直接采用; 否则回退到按距离划分 (旧行为, 逐字节不变)。
    left: list[int] = []
    right: list[int] = []
    if data is not None:
        lk = data.get("left_electrode_indices")
        rk = data.get("right_electrode_indices")
        if lk and rk:
            left = [
                int(i) - 1 for i in lk if int(i) - 1 not in (left_tip, right_tip)
            ]
            right = [
                int(i) - 1 for i in rk if int(i) - 1 not in (left_tip, right_tip)
            ]
            if sorted(left + right) != sorted(electrode):
                raise ValueError(
                    "metadata left/right_electrode_indices 与几何不符: "
                    f"metadata {len(left) + len(right)} 个, 几何 {len(electrode)} 个"
                )
    if not left or not right:
        left = [i for i in electrode if _d2(i, left_tip) < _d2(i, right_tip)]
        right = [i for i in electrode if _d2(i, right_tip) < _d2(i, left_tip)]
    if len(left) != au_per_electrode - 1 or len(right) != au_per_electrode - 1:
        raise ValueError(f"Left/right split failed: left={len(left)} right={len(right)}")

    def _sigma(electrode_ids: list[int], tip: int) -> dict[int, float]:
        ids = sorted(electrode_ids, key=lambda i: _d2(i, tip))
        return {i: (sigma_near if k < 3 else sigma_far) for k, i in enumerate(ids)}

    left_sig = _sigma(left, left_tip)
    right_sig = _sigma(right, right_tip)

    lines = ["$self.energy: imaginary piece per atom"]
    for idx, a in enumerate(atoms, 1):
        if idx - 1 in (left_tip, right_tip) or a.element != "Au":
            region, sigma = "empty", 0.0
        elif idx - 1 in left_sig:
            region, sigma = "left", left_sig[idx - 1]
        elif idx - 1 in right_sig:
            region, sigma = "right", right_sig[idx - 1]
        else:
            raise ValueError(f"Unclassified atom {idx}")
        lines.append(
            "%5d %15.7f %15.7f %15.7f  %-5s %-9s %s"
            % (idx, a.x, a.y, a.z, a.element.lower(), region, fortran_d(sigma))
        )
    lines.append("$end")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {
        "output": str(out),
        "atom_count": len(atoms),
        "left_tip": left_tip + 1,
        "right_tip": right_tip + 1,
        "left_electrode": [i + 1 for i in left],
        "right_electrode": [i + 1 for i in right],
        "sigma_near": sigma_near,
        "sigma_far": sigma_far,
    }


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: python -m agenticdft.tools.junction.self_energy "
              "<geometry.in> <self.energy.in> [electrode_metadata.json]")
        sys.exit(2)
    meta = sys.argv[3] if len(sys.argv) > 3 else None
    result = make_self_energy(sys.argv[1], sys.argv[2], metadata=meta)
    print(json.dumps(result, indent=2))