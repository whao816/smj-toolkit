"""按原子分组投影的 PDOS（2026-11-20，规划 §6 P2-8）。

输入是 FHI-aims 静态 SCF 的三件产物：

* ``basis-indices.out`` —— 把"基函数编号"映射到"原子编号"（``# fn. type at. n l m``）；
* ``mos.aims``（开壳层是 ``alpha.aims`` / ``beta.aims``）—— 每个 KS 态的本征值与本征向量；
* ``TE.dat`` —— 只为取 ``efermi`` 与那条已经用于 T(E) 的能量轴。

投影口径（写进结果与报告，不隐藏）：

    w(group, state) = occupation * Σ_{ν ∈ group} c_{ν,state}²

即把该态在本组基函数上的系数平方和当作本组权重（Löwdin/正交基下的 Mulliken 型投影）。
**不做重叠矩阵修正**（``omat.aims`` 不参与），所以它不是严格的 Mulliken 布居；
同族的能量积分恒等于总权重，跨组求和恒等于 1（乘占据数后为总电子数），可以用来校验。

再按高斯展宽（默认 σ = 0.1 eV）落到 TE.dat 的 ``E − E_F`` 网格上，因此 PDOS 与 T(E)
可以直接画在同一能量轴。

已知边界（会写进 ``notes``，不静默降级）：

* ``mos.aims`` 在回取清单里属于 ``head`` 文件，很多 run 本地只有头部（``nsaos`` 行 +
  第一个态的本征值）。这时 :func:`parse_mos` 抛 :class:`PdosError` 并说清缺什么，
  不会拿"只有本征值"的数据凑一份 PDOS。
* 权重是"系数平方和"，不是分轨道（s/p/d）或分 m 的投影；需要分 l 时按
  ``basis-indices.out`` 的 ``l`` 列再切一刀即可（本轮不做）。
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .errors import SmjError

HARTREE_TO_EV = 27.211386245988
PDOS_SCHEMA = "smj.pdos.v1"
DEFAULT_SIGMA_EV = 0.1

_MOS_HEADER_RE = re.compile(
    r"^\s*\$(?P<kind>scfmo|uhfmo_alpha|uhfmo_beta)\.aims(?:.*?iteration:\s*(?P<iteration>\d+))?"
    r"(?:.*?format\(\s*(?P<per_line>\d+)[dD](?P<width>\d+)\.(?P<decimals>\d+)\s*\))?")
_STATE_RE = re.compile(
    r"^\s*(\d+)\s+([ab])\s+eigenvalue=\s*(\S+)\s+"
    r"(?:occupation=\s*(\S+)\s+)?nsaos=\s*(\d+)\s*$"
)
_BASIS_COUNT_RE = re.compile(r"amount of basis functions used:\s*(\d+)")
_BASIS_ROW_RE = re.compile(r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+([-+]?\d+)\s+([-+]?\d+)\s*$")
_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[dDeE][-+]?\d+)?")


class PdosError(SmjError, ValueError):
    """PDOS 无法从现有产物计算（缺文件、文件被截断、编号对不上）。"""


def _float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def _read(path: str | Path) -> str:
    p = Path(path)
    if not p.is_file():
        raise PdosError("找不到文件：%s" % p)
    return p.read_text(encoding="utf-8", errors="ignore")


# --------------------------------------------------------------------------- #
# 解析
# --------------------------------------------------------------------------- #
def parse_basis_indices(path: str | Path) -> dict[str, Any]:
    """解析 ``basis-indices.out``，返回基函数 → 原子的映射。

    Returns
    -------
    dict
        ``{"path", "nsaos", "atom_of_fn"}``；``atom_of_fn[i]`` 对应基函数 ``i + 1``
        （1-based，与 tcontrol 的 ``$nsaos`` 一致）。
    """
    text = _read(path)
    nsaos: int | None = None
    atom_of_fn: list[int] = []
    for line in text.splitlines():
        if nsaos is None:
            m = _BASIS_COUNT_RE.search(line)
            if m:
                nsaos = int(m.group(1))
        m = _BASIS_ROW_RE.match(line)
        if not m:
            continue
        fn, _kind, atom = int(m.group(1)), m.group(2), int(m.group(3))
        if fn != len(atom_of_fn) + 1:
            raise PdosError(
                "basis-indices.out 的基函数编号不连续：期望 %d，实际 %d（%s）"
                % (len(atom_of_fn) + 1, fn, path)
            )
        atom_of_fn.append(atom)
    if not atom_of_fn:
        raise PdosError("basis-indices.out 里没有基函数行（%s）" % path)
    declared = nsaos if nsaos is not None else len(atom_of_fn)
    if declared != len(atom_of_fn):
        raise PdosError(
            "basis-indices.out 声明 %d 个基函数，实际解析到 %d 行（%s）"
            % (declared, len(atom_of_fn), path)
        )
    return {"path": str(path), "nsaos": declared, "atom_of_fn": atom_of_fn}


def _channel_of(kind: str) -> str:
    return {"scfmo": "closed", "uhfmo_alpha": "alpha", "uhfmo_beta": "beta"}.get(kind, "closed")


def _split_coefficients(line: str, width: int | None) -> list[float]:
    """切一行系数。

    FHI-aims 用 ``format(4d20.14)`` 写定宽字段（无分隔符），相邻负数会粘成
    ``…D-08-.752…``；按 20 字符定宽切才安全。测试夹具用空白分隔，``width``
    为 ``None`` 时退回 token 切分。
    """
    if width:
        values: list[float] = []
        for start in range(0, len(line), width):
            chunk = line[start:start + width].strip()
            if chunk:
                values.append(_float(chunk))
        return values
    return [_float(token) for token in _FLOAT_RE.findall(line)]


def _scan_states(path: str | Path, *, group_of_fn: Sequence[int | None] | None = None):
    """流式扫描波函数文件，逐态产出 meta + 按分组的投影权重。

    FHI-aims 的本征向量满足 ``cᵀ S c = 1``（NAO 基不是正交基），``Σc²`` 本身
    不是归一化权重，也没有 Mulliken 意义。这里按每个态的 ``Σ_all c²`` 归一：
    ``w_g = Σ_{ν∈g} c² / Σ_all c²``（Mulliken 需要 ``omat.aims`` 的重叠矩阵，
    本地只有 60 字节的 head）。组间比例可读，绝对量不能当 Mulliken 布居。

    不保留系数，383 MB 的 ``alpha.aims`` 也不会把内存撑爆。
    """
    p = Path(path)
    if not p.is_file():
        raise PdosError("找不到文件：%s" % p)
    groups = len({g for g in (group_of_fn or []) if g is not None})
    current: dict[str, Any] | None = None
    weights: list[float] = [0.0] * groups
    count = 0
    total_sq = 0.0
    width: int | None = None
    channel = "closed"
    iteration: int | None = None
    nsaos_declared: int | None = None

    def _close():
        nonlocal current, weights, count, total_sq
        if current is None:
            return
        want = int(current["nsaos"])
        if count != want:
            raise PdosError(
                "%s 的第 %d 个态只有 %d/%d 个系数：这是被截断的波函数文件"
                "（本地通常只回取了 head）。PDOS 需要完整的本征向量；请用 "
                "`smj pdos --fetch` 按显式清单整取，或重跑 aims_scf 后保留完整文件。"
                % (path, current["index"], count, want)
            )
        if group_of_fn is not None:
            scale = 1.0 / total_sq if total_sq > 0 else 0.0
            current["weights"] = [value * scale for value in weights]
            current["sum_c_squared"] = total_sq
        else:
            current["weights"] = None
        yield_state = dict(current)
        current, weights, count, total_sq = None, [0.0] * groups, 0, 0.0
        return yield_state

    with p.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            header = _MOS_HEADER_RE.match(line)
            if header:
                channel = _channel_of(header.group("kind"))
                if header.group("iteration"):
                    iteration = int(header.group("iteration"))
                if header.group("width"):
                    width = int(header.group("width"))
                continue
            state = _STATE_RE.match(line)
            if state:
                ready = _close()
                if ready is not None:
                    yield ready
                idx, _spin, eig, occ, nsaos = state.groups()
                nsaos_declared = int(nsaos)
                current = {
                    "index": int(idx),
                    "spin": channel,
                    "eigenvalue_Ha": _float(eig),
                    "occupation": _float(occ) if occ is not None else None,
                    "nsaos": int(nsaos),
                    "iteration": iteration,
                }
                continue
            if current is None:
                continue
            if line.lstrip().startswith("$"):
                break
            values = _split_coefficients(line, width)
            if not values:
                continue
            for offset, value in enumerate(values):
                fn = count + offset + 1              # 1-based
                if nsaos_declared is not None and fn > nsaos_declared:
                    break
                total_sq += value * value
                if group_of_fn is not None:
                    gi = group_of_fn[fn - 1]
                    if gi is not None:
                        weights[gi] += value * value
            count += len(values)
    ready = _close()
    if ready is not None:
        yield ready


def parse_mos(path: str | Path) -> dict[str, Any]:
    """解析 ``mos.aims`` / ``alpha.aims`` / ``beta.aims``（测试与调试用）。

    生产路径走 :func:`_scan_states` 的流式投影，不会把整个文件读进内存。
    """
    states: list[dict[str, Any]] = []
    channel = "closed"
    nsaos: int | None = None
    # 需要系数本身时走一遍定宽/分词解析（文件通常很小；大文件请用 _scan_states）
    text = _read(path)
    width: int | None = None
    current: dict[str, Any] | None = None
    coefficients: list[float] = []

    def _flush() -> None:
        nonlocal current, coefficients
        if current is None:
            return
        want = int(current["nsaos"])
        if len(coefficients) != want:
            raise PdosError(
                "%s 的第 %d 个态只有 %d/%d 个系数：这是被截断的波函数文件"
                "（本地通常只回取了 head）。PDOS 需要完整的本征向量；请用 "
                "`smj pdos --fetch` 按显式清单整取，或重跑 aims_scf 后保留完整文件。"
                % (path, current["index"], len(coefficients), want)
            )
        current["coefficients"] = list(coefficients)
        states.append(current)
        current, coefficients = None, []

    for line in text.splitlines():
        header = _MOS_HEADER_RE.match(line)
        if header:
            channel = _channel_of(header.group("kind"))
            if header.group("width"):
                width = int(header.group("width"))
            continue
        state = _STATE_RE.match(line)
        if state:
            _flush()
            idx, _spin, eig, occ, nsaos_text = state.groups()
            nsaos = int(nsaos_text)
            current = {"index": int(idx), "spin": channel,
                       "eigenvalue_Ha": _float(eig),
                       "occupation": _float(occ) if occ is not None else None,
                       "nsaos": nsaos}
            continue
        if current is None:
            continue
        if line.lstrip().startswith("$"):
            break
        coefficients.extend(_split_coefficients(line, width))
    _flush()
    if not states:
        raise PdosError("mos.aims 里没有可识别的 KS 态（%s）" % path)
    return {"path": str(path), "iteration": None, "nsaos": nsaos,
            "states": states, "channel": channel}


# --------------------------------------------------------------------------- #
# 分组与展宽
# --------------------------------------------------------------------------- #
def parse_groups(spec: str | Sequence[Any], *, natoms: int) -> list[dict[str, Any]]:
    """把 ``"molecule=1-20,left=21-52"`` 解析成 ``[{"name","atoms"}]``。

    分隔符 ``=`` 或 ``:`` 都接受；原子编号 1-based、闭区间，可写 ``3,7-9``。
    """
    if isinstance(spec, str):
        items: list[tuple[str, str]] = []
        for chunk in spec.split(","):
            chunk = chunk.strip()
            if not chunk:
                continue
            if "=" in chunk:
                name, body = chunk.split("=", 1)
            elif ":" in chunk:
                name, body = chunk.split(":", 1)
            else:
                raise PdosError("分组写法应为 name=1-20 或 name:1-20，实际是 %r" % chunk)
            items.append((name.strip(), body.strip()))
    else:
        items = [(str(entry[0]), str(entry[1])) for entry in spec]  # type: ignore[index]
    if not items:
        raise PdosError("至少要给一个原子分组（例：molecule=1-20,left=21-52）")

    groups: list[dict[str, Any]] = []
    seen: dict[int, str] = {}
    for name, body in items:
        atoms: list[int] = []
        for part in body.split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part.lstrip("-"):
                lo_s, hi_s = part.split("-", 1)
                lo, hi = int(lo_s), int(hi_s)
            else:
                lo = hi = int(part)
            if lo > hi:
                lo, hi = hi, lo
            if lo < 1 or hi > natoms:
                raise PdosError(
                    "分组 %s 的原子区间 %s 超出 1..%d" % (name, part, natoms)
                )
            atoms.extend(range(lo, hi + 1))
        if not atoms:
            raise PdosError("分组 %s 没有原子" % name)
        for atom in atoms:
            if atom in seen:
                raise PdosError("原子 %d 同时出现在分组 %s 与 %s" % (atom, seen[atom], name))
            seen[atom] = name
        groups.append({"name": name, "atoms": atoms})
    return groups


def broaden(values_ev: Iterable[float], weights: Iterable[float], grid_ev: Sequence[float],
            sigma_ev: float) -> list[float]:
    """高斯展宽：``Σ_s w_s · N(E − ε_s; σ)``，``N`` 归一化到积分 1。"""
    if sigma_ev <= 0:
        raise PdosError("σ 必须为正（当前 %r eV）" % sigma_ev)
    norm = 1.0 / (sigma_ev * math.sqrt(2.0 * math.pi))
    out = [0.0] * len(grid_ev)
    for value, weight in zip(values_ev, weights):
        if weight == 0.0:
            continue
        for i, energy in enumerate(grid_ev):
            dx = (energy - value) / sigma_ev
            out[i] += weight * norm * math.exp(-0.5 * dx * dx)
    return out


_OCC_HEADER_RE = re.compile(r"State\s+Occupation\s+Eigenvalue \[Ha\]\s+Eigenvalue \[eV\]")
_OCC_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?)\s+"
    r"[-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?\s+[-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?\s*$")


def parse_occupations(aims_out: str | Path) -> dict[str, Any]:
    """从 ``aims.out`` 的 ``State Occupation Eigenvalue`` 表读占据数。

    FHI-aims 的 ``.aims`` 波函数文件只写本征值、不写占据数；自旋极化时
    ``aims.out`` 在末尾按 spin 各打一张表（先 α 后 β）。取最后两张等长的表
    作为 α/β，只有一张则视为闭壳层。
    """
    p = Path(aims_out)
    if not p.is_file():
        raise PdosError("找不到 aims.out：%s（PDOS 需要它给出态占据数）" % p)
    tables: list[list[float]] = []
    current: list[float] | None = None
    with p.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if _OCC_HEADER_RE.search(line):
                if current:
                    tables.append(current)
                current = []
                continue
            if current is None:
                continue
            row = _OCC_ROW_RE.match(line)
            if row:
                current.append(_float(row.group(2)))
                continue
            if current:
                tables.append(current)
            current = None
    if current:
        tables.append(current)
    if not tables:
        raise PdosError("aims.out 里没有 State/Occupation 表：%s" % p)
    same = [t for t in tables if len(t) == len(tables[-1])]
    if len(same) >= 2:
        return {"source": str(p), "tables": len(tables),
                "alpha": same[-2], "beta": same[-1]}
    return {"source": str(p), "tables": len(tables),
            "alpha": tables[-1], "beta": None}


def energy_grid(te_path: str | Path) -> dict[str, Any]:
    """从 TE.dat 取能量轴（``E − E_F`` eV，已按升序）与 ``efermi``。"""
    from smj.analysis import te_header  # 表头口径与 report/analysis 共用
    from smj.junction.landauer import analyze_te as _analyze  # 复用已排序、已校验的读取
    from smj.junction.landauer import read_te_ordered

    rows = _analyze(te_path)
    data, _ = read_te_ordered(te_path)      # 2026-11-20 复审 P1-3: 单一能量序
    grid = [float(row[1]) for row in data]
    header = te_header(te_path)
    return {
        "te_path": str(te_path),
        "grid_eV": grid,
        "points": len(grid),
        "efermi_Ha": header.get("efermi_Ha"),
        "energy_order": rows.get("input_energy_order"),
    }


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #
def _find(run_dir: Path, *names: str) -> Path | None:
    for name in names:
        candidate = run_dir / name
        if candidate.is_file():
            return candidate
    return None


def _channel_files(run_dir: Path) -> dict[str, Path | None]:
    scf = run_dir / "aims_scf" if (run_dir / "aims_scf").is_dir() else run_dir
    alpha = _find(scf, "alpha.aims")
    beta = _find(scf, "beta.aims")
    if alpha is not None or beta is not None:
        # 自旋极化 SCF 不写 mos.aims；即使目录里残留一份，也以 alpha/beta 为准
        return {"alpha": alpha, "beta": beta}
    return {"mos": _find(scf, "mos.aims"), "alpha": None, "beta": None}


def analyze_pdos(
    run_dir: str | Path,
    *,
    groups: str | Sequence[Any],
    sigma_ev: float = DEFAULT_SIGMA_EV,
    te: str | Path | None = None,
    output: str | Path | None = None,
    plot: str | Path | bool | None = None,
) -> dict[str, Any]:
    """算一个 run 的分组 PDOS，返回 ``smj.pdos.v1`` 字典。

    ``run_dir`` 里应能找到 ``aims_scf/{mos.aims|alpha.aims|beta.aims}``、
    ``aims_scf/basis-indices.out``（扁平布局也接受）与 ``aitranss_te/TE.dat``。
    """
    run = Path(run_dir)
    te_path = Path(te) if te else (run / "aitranss_te" / "TE.dat")
    if not te_path.is_file():
        alt = run / "TE.dat"
        if alt.is_file():
            te_path = alt
        else:
            raise PdosError("找不到 TE.dat：%s（PDOS 需要它的能量轴与 efermi）" % te_path)

    scf = run / "aims_scf" if (run / "aims_scf").is_dir() else run
    basis_path = scf / "basis-indices.out"
    basis = parse_basis_indices(basis_path)
    natoms = max(basis["atom_of_fn"])
    groups_parsed = parse_groups(groups, natoms=natoms)

    # 基函数 → 分组
    group_functions: list[list[int]] = [[] for _ in groups_parsed]
    group_of_atom: dict[int, int] = {}
    for gi, group in enumerate(groups_parsed):
        for atom in group["atoms"]:
            group_of_atom[atom] = gi
    for fn, atom in enumerate(basis["atom_of_fn"], start=1):
        gi = group_of_atom.get(atom)
        if gi is not None:
            group_functions[gi].append(fn)

    files = _channel_files(run)
    channels = {name: path for name, path in files.items() if path is not None}
    if not channels:
        raise PdosError(
            "aims_scf 下没有 mos.aims / alpha.aims / beta.aims：PDOS 需要本征向量（%s）" % scf
        )

    grid = energy_grid(te_path)
    efermi = grid["efermi_Ha"]
    if efermi is None:
        raise PdosError(
            "TE.dat 表头没有 efermi，无法把本征值换到 E − E_F 轴（%s）" % te_path
        )

    group_of_fn: list[int | None] = [group_of_atom.get(atom) for atom in basis["atom_of_fn"]]
    occupation_tables = None
    aims_out = _find(scf, "aims.out")
    if aims_out is not None:
        occupation_tables = parse_occupations(aims_out)

    notes: list[str] = [
        "权重按每个态的 Σ_all c² 归一（w_group = Σ_{ν∈group} c² / Σ_all c²）：NAO 基非正交，"
        "真正的 Mulliken 布居需要 omat.aims 的重叠矩阵（本地只有 head）；组间比例可读，"
        "绝对量不能当 Mulliken 电荷",
        "PDOS 与 T(E) 共用 TE.dat 的 E − E_F 网格（%d 点）；σ = %.3g eV" % (grid["points"], sigma_ev),
    ]
    if occupation_tables is not None:
        notes.append(
            "占据数取自 %s 的 State Occupation 表（共 %d 张%s）"
            % (occupation_tables["source"], occupation_tables["tables"],
               "，末两张为 α/β 终态" if occupation_tables["beta"] else "，单张=闭壳层"))
    else:
        notes.append("没有 aims.out：按每个态占据 1 计权（只能看态密度形状，不能当占据数解释）")
    result: dict[str, Any] = {
        "schema": PDOS_SCHEMA,
        "run_dir": str(run),
        "te_path": str(te_path),
        "efermi_Ha": efermi,
        "sigma_eV": float(sigma_ev),
        "natoms": natoms,
        "nsaos": basis["nsaos"],
        "groups": [
            {"name": group["name"], "atoms": group["atoms"],
             "functions": len(group_functions[gi])}
            for gi, group in enumerate(groups_parsed)
        ],
        "grid_eV": grid["grid_eV"],
        "channels": {},
        "notes": notes,
    }
    for name, path in channels.items():
        table = None
        if occupation_tables is not None:
            table = occupation_tables["beta"] if name == "beta" else occupation_tables["alpha"]
        values_ev: list[float] = []
        weights: dict[str, list[float]] = {group["name"]: [] for group in groups_parsed}
        states = 0
        for index, state in enumerate(_scan_states(path, group_of_fn=group_of_fn)):
            occupation = state["occupation"]
            if occupation is None and table is not None and index < len(table):
                occupation = table[index]
            if occupation is None:
                occupation = 1.0
            values_ev.append((state["eigenvalue_Ha"] - float(efermi)) * HARTREE_TO_EV)
            per_group = state["weights"] or [0.0] * len(groups_parsed)
            for gi, group in enumerate(groups_parsed):
                weights[group["name"]].append(occupation * per_group[gi])
            states += 1
        if not states:
            raise PdosError("%s 里没有可识别的 KS 态" % path)
        pdos = {key: broaden(values_ev, weight, result["grid_eV"], sigma_ev)
                for key, weight in weights.items()}
        result["channels"][name] = {
            "file": str(path),
            "states": states,
            "spin": "closed" if name == "mos" else name,
            "integrated_weight": {key: sum(weight) for key, weight in weights.items()},
            "state_energies_eV": values_ev,
            "pdos": pdos,
        }

    if output:
        target = Path(output)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                          encoding="utf-8", newline="\n")
        result["json_path"] = str(target)

    if plot:
        target = Path(plot) if plot is not True else (run / "PDOS.png")
        try:
            result["plot"] = plot_pdos(result, te_path, target)
        except Exception as exc:  # noqa: BLE001 —— matplotlib 是可选依赖
            result["plot"] = {"skipped": True,
                              "reason": "matplotlib 不可用（extra \"plot\" 未安装）：%s" % exc}
    return result


def plot_pdos(result: Mapping[str, Any], te_path: str | Path, target: str | Path) -> dict:
    """把 PDOS 与 T(E) 画在同一能量轴上（上图 T(E)，下图分组 PDOS）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from smj.junction.landauer import read_te_ordered

    grid = list(result["grid_eV"])
    rows, _ = read_te_ordered(te_path)
    by_energy = {float(row[1]): float(row[2]) for row in rows}
    transmission = [by_energy.get(energy) for energy in grid]

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True,
                                     gridspec_kw={"height_ratios": [1.0, 1.4]})
    top.semilogy(grid, [max(value, 1e-12) if value else 1e-12 for value in transmission],
                 color="#1f77b4", lw=1.2)
    top.set_ylabel("T(E) [G/G$_0$]")
    top.set_title("PDOS %s" % result.get("run_dir", ""))
    colours = ["#d62728", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]
    for index, channel in enumerate(result["channels"]):
        style = "-" if channel in ("mos", "alpha") else "--"
        label_suffix = "" if channel == "mos" else " (%s)" % channel
        for gi, group in enumerate(result["groups"]):
            series = result["channels"][channel]["pdos"][group["name"]]
            bottom.plot(grid, series, style, color=colours[gi % len(colours)], lw=1.2,
                        label="%s%s" % (group["name"], label_suffix))
    bottom.axvline(0.0, color="0.4", ls=":", lw=1.0)
    bottom.set_xlabel("E − E$_F$ [eV]")
    bottom.set_ylabel("PDOS [states/eV]")
    bottom.legend(frameon=False, fontsize=8)
    for axis in (top, bottom):
        axis.tick_params(direction="in", which="both", top=True, right=True)
    fig.tight_layout()
    path = Path(target)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return {"path": str(path), "skipped": False}
