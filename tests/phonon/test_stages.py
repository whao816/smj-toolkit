"""声子热导链的两个 control.in 模板(thermal_relax / thermal_force)。

为什么值得单独钉死(方案附录 A-3):FHI-aims 里打印力的关键字是
``compute_forces``,写成 ``calculate_forces`` **不会报错**,只会静默不给力 ——
整个 324 次位移扫描算完才发现 aims.out 里没有力块。所以这里数的是关键字出现的
次数,不是「包含」。
"""
from __future__ import annotations

import pathlib

import pytest

from smj.config import ConfigError, load_config
from smj.junction.geometry import Atom, write_aims_geometry
from smj.junction.tcontrol import make_control

_ATOMS = [
    Atom("Au", 0.0, 0.0, 0.0),
    Atom("S", 0.0, 0.0, 2.4),
    Atom("C", 1.4, 0.0, 3.9),
]


@pytest.fixture()
def setup(workdir):
    """几何 + 假的 species 目录(每个元素一个 ``<El>.in``)。"""
    geometry = workdir / "geometry.in"
    write_aims_geometry(geometry, _ATOMS,
                        lattice=((4.07823, 0.0, 0.0), (0.0, 4.07823, 0.0),
                                 (0.0, 0.0, 25.0)))
    species = workdir / "species"
    species.mkdir()
    for element in ("Au", "S", "C"):
        (species / f"{element}.in").write_text(
            "# %s block\n  species %s\n" % (element, element), encoding="utf-8")
    return geometry, species


def _render(setup, mode: str, **kwargs) -> str:
    geometry, species = setup
    out = geometry.parent / ("control_%s.in" % mode)
    make_control(geometry, out, mode, species,
                 k_grid=(3, 3, 1) if mode.startswith("thermal") else None,
                 **kwargs)
    return out.read_text(encoding="utf-8")


def _active(text: str) -> str:
    """只留 FHI-aims 会读的行:去掉 ``#`` 注释。

    模板里的注释**故意**提到 ``calculate_forces`` / ``output aitranss`` 这些
    名字(解释为什么不写它们)。直接在全文里搜关键字会被自己的注释绊倒。
    """
    return "\n".join(
        line for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#"))


def test_thermal_force_has_exactly_one_compute_forces(setup):
    text = _active(_render(setup, "thermal_force"))
    assert text.count("compute_forces .true.") == 1
    # 拼错的关键字会被 FHI-aims 静默忽略 —— 一个都不能出现
    assert "calculate_forces" not in text


def test_thermal_force_has_k_grid_and_no_aitranss_outputs(setup):
    text = _active(_render(setup, "thermal_force"))
    assert "k_grid 3 3 1" in text
    assert "output aitranss" not in text
    assert "output mulliken" not in text
    assert "restart aims.restart" not in text
    assert "final_forces_cleaned .true." in text
    assert "sc_accuracy_forces 5E-3" in text


def test_thermal_relax_adds_bfgs(setup):
    text = _active(_render(setup, "thermal_relax"))
    assert "relax_geometry bfgs 1.e-2" in text
    assert text.count("compute_forces .true.") == 1
    assert "k_grid 3 3 1" in text
    assert "output aitranss" not in text


def test_thermal_modes_require_k_grid(setup):
    geometry, species = setup
    with pytest.raises(ValueError):
        make_control(geometry, geometry.parent / "c.in", "thermal_force", species)


def test_kz_must_be_one(setup):
    geometry, species = setup
    with pytest.raises(ValueError):
        make_control(geometry, geometry.parent / "c.in", "thermal_force", species,
                     k_grid=(3, 3, 4))


def test_k_grid_must_be_three_positive_ints(setup):
    geometry, species = setup
    for bad in ((3, 3), (3, 3, 0), None):
        with pytest.raises(ValueError):
            make_control(geometry, geometry.parent / "c.in", "thermal_force",
                         species, k_grid=bad)


def test_existing_scf_mode_is_unchanged(setup):
    """既有 scf 口径不能因为新增 mode 而漂:aitranss/mulliken/restart 仍在。"""
    text = _render(setup, "scf")
    assert "output aitranss" in text
    assert "output mulliken" in text
    assert "restart aims.restart" in text
    assert "compute_forces" not in text
    assert "k_grid" not in text


def test_unknown_mode_still_rejected(setup):
    geometry, species = setup
    with pytest.raises(ValueError):
        make_control(geometry, geometry.parent / "c.in", "thermal", species)


# --- thermal 配置段 ---------------------------------------------------------

_MINIMAL = "run:\n  name: t\n  mode: dry-run\n"


def _load(workdir, extra: str = ""):
    path = workdir / "smj.yaml"
    path.write_text(_MINIMAL + extra, encoding="utf-8")
    return load_config(path)


def test_thermal_disabled_by_default(workdir):
    """不写 thermal 段时电子链行为零变化:全部默认值,enabled=False。"""
    config = _load(workdir)
    thermal = config.thermal
    assert thermal.enabled is False
    assert thermal.temperature_K is None
    assert thermal.lead_fc is None
    assert thermal.electrode.in_plane == (3, 3)
    assert thermal.electrode.layers == 4
    assert thermal.electrode.device_layers == 2
    assert thermal.electrode.strip_thiol_h is True
    assert thermal.fc.displacement == pytest.approx(0.01)
    assert thermal.fc.cutoff == pytest.approx(10.0)
    assert thermal.fc.array_parallel == 8
    assert thermal.fc.k_grid == (3, 3, 1)
    # D15 冻结协议：η 与步长都取 0.05（见 ThermalNEGFConfig 的 docstring）
    assert thermal.negf.eta_meV == pytest.approx(0.05)
    assert thermal.negf.energy_step_meV == pytest.approx(0.05)
    assert thermal.negf.band_top_meV == pytest.approx(19.5)
    # 温度缺省继承 transport.temperature(默认 300 K)
    assert thermal.resolved_temperature(config.transport.temperature) == 300.0


def test_thermal_section_parsed_and_paths_resolved(workdir):
    config = _load(workdir, """
thermal:
  enabled: true
  temperature: 250.0
  lead_fc: leads/au111.npz
  electrode:
    in_plane: [4, 4]
    layers: 6
    device_layers: 3
    vacuum: 18.0
  fc:
    displacement: 0.005
    k_grid: [2, 2, 1]
  negf:
    eta_mev: 0.25
""")
    thermal = config.thermal
    assert thermal.enabled is True
    assert thermal.resolved_temperature(300.0) == pytest.approx(250.0)
    assert thermal.electrode.in_plane == (4, 4)
    assert thermal.electrode.layers == 6
    assert thermal.electrode.device_layers == 3
    assert thermal.electrode.vacuum == pytest.approx(18.0)
    assert thermal.fc.displacement == pytest.approx(0.005)
    assert thermal.fc.k_grid == (2, 2, 1)
    assert thermal.negf.eta_meV == pytest.approx(0.25)
    # lead_fc 相对 yaml 所在目录解析成绝对路径
    assert pathlib.Path(thermal.lead_fc).is_absolute()
    assert (pathlib.Path(thermal.lead_fc)
            == (workdir / "leads" / "au111.npz").resolve())


@pytest.mark.parametrize("extra", [
    "thermal:\n  fc:\n    k_grid: [3, 3, 4]\n",       # kz 必须是 1
    "thermal:\n  electrode:\n    layers: 4\n    device_layers: 4\n",
    "thermal:\n  fc:\n    displacement: 0.0\n",
    "thermal:\n  fc:\n    displacement: 0.2\n",
    "thermal:\n  fc:\n    cutoff: 6.5\n",              # 4 层 lead 跨度就是 7.06 Å
    "thermal:\n  negf:\n    eta_mev: 0.0\n",
    "thermal:\n  temperature: -5\n",
])
def test_thermal_invalid_values_rejected(workdir, extra):
    with pytest.raises(ConfigError):
        _load(workdir, extra)


def test_thermal_unknown_keys_warn_not_fail(workdir):
    """未知字段只 warning(向前兼容),与其它段一致。"""
    config = _load(workdir, "thermal:\n  enabled: true\n  bogus: 1\n"
                            "  fc:\n    bogus2: 2\n")
    assert config.thermal.enabled is True
    joined = " ".join(config.warnings)
    assert "thermal" in joined
    assert "bogus" in joined
    assert "bogus2" in joined


# --------------------------------------------------------------------------
# 模板的字节卫生:E/FHI-aims 与 sbatch 都会把控制字符当成语法的一部分
# --------------------------------------------------------------------------

def test_thermal_templates_have_no_control_characters():
    """模板里不得出现 CR / FF 之类的控制字符。

    为什么钉这条(2026-11-23 真实事故):写模板时文本里混进了 ``\\r`` 与 ``\\f``,
    它们落成真实控制字符之后,``control.in`` 里一条注释被 CR 拆成两行,后半行
    以 ``restart aims.restart`` 开头,于是 FHI-aims 把它当真关键字,报
    ``Unknown keyword estart`` 并停止 —— 作业白跑,而且报错信息指向的关键字在
    模板里根本不存在。sbatch 也会因为 CRLF 直接拒收脚本。
    """
    from smj.junction.tcontrol import TEMPLATE_DIR

    for name in ("control_thermal_force.in", "control_thermal_relax.in",
                 "sub_thermal_relax.sh", "sub_thermal_fc.sh"):
        path = pathlib.Path(TEMPLATE_DIR) / name
        raw = path.read_bytes()
        assert b"\r" not in raw, "%s 含 CR" % name
        bad = [c for c in raw if c < 0x20 and c not in (0x0A, 0x09)]
        assert not bad, "%s 含控制字符 %r" % (name, [hex(c) for c in bad[:4]])
        for line in raw.decode("utf-8").split("\n"):
            assert not line.lstrip().startswith("restart"), (
                "%s 有一行以 restart 开头: %r" % (name, line))

def test_thermal_scripts_request_no_gpu(workdir):
    """gpu1 节点带 Gres=gpu:4:不显式要 gpu:0 时整个节点只容得下一个作业。

    实测（2026-11-23）``sbatch --test-only -n 1 --mem=100M`` 会把新作业排到
    当前作业时限结束,加 ``--gres=gpu:0`` 后立刻可跑;热导链的 FHI-aims 是纯 CPU
    (input_summary.json 的 aims_gpus 为 "cpu"),不要卡。
    """
    from smj.stages._thermal_common import render_thermal_sub

    env = {"partition": "main", "aims_mem": "12G", "aims_walltime": "06:00:00",
           "aims_bin": "aims.x"}
    relax = render_thermal_sub("sub_thermal_relax.sh", env, ntasks=8)
    fc = render_thermal_sub("sub_thermal_fc.sh", env, ntasks=8,
                            extra={"last": 3, "parallel": 2, "count": 4})
    for text in (relax, fc):
        assert "#SBATCH --gres=gpu:0" in text
        assert "--gres=gpu:4" not in text
    # 显式配了 gres 就照用
    other = render_thermal_sub("sub_thermal_relax.sh",
                               dict(env, gres="gpu:1"), ntasks=8)
    assert "#SBATCH --gres=gpu:1" in other
    assert "gpu:0" not in other
def test_fc_script_calls_extract_forces_with_two_arguments():
    """extract_forces 的 CLI 是 ``extract_forces.py <aims.out> <forces.json>``。

    模板里曾写成 ``python ../../extract_forces.py aims.out > forces.json``：
    重定向把用法提示写进了 forces.json、脚本返回 2，数组任务因此全部失败
    (job 26443,54 个 aims.out 都算完了但一个 done 都没写)。
    """
    from smj.stages._thermal_common import render_thermal_sub

    env = {"partition": "main", "aims_mem": "12G", "aims_walltime": "06:00:00",
           "aims_bin": "aims.x"}
    text = render_thermal_sub("sub_thermal_fc.sh", env, ntasks=8,
                              extra={"last": 3, "parallel": 2, "count": 4})
    assert "extract_forces.py aims.out forces.json" in text
    assert "extract_forces.py aims.out >" not in text

def test_negf_defaults_are_the_frozen_protocol():
    """D15：η 与步长是冻结的生产协议，改动等于改方法。

    实测依据(同一份力常数)：η=0.5 meV 时 G_ph(300 K)=10.9 pW/K、0.3→12.37、
    0.1→13.61、0.05→13.72，η→0 外推 13.8。跨分子比较不同 η 会凭空造出 25%
    的系统差别，所以默认值必须是协议值。
    """
    from smj.config import ThermalNEGFConfig

    cfg = ThermalNEGFConfig()
    assert cfg.eta_meV == 0.05
    assert cfg.energy_step_meV == 0.05
