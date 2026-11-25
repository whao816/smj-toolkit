"""诊断:thermal_analyze 在带顶之上的透射到底有多大、lead 模型的真实带顶在哪。"""
import json
from pathlib import Path

import numpy as np

from smj.phonon.dynamical import dynamical_matrix, imaginary_modes, masses_from_elements, mode_frequencies_cm1
from smj.phonon.fc import FullFCMeta, asr_diagnostics, assemble_full_fc, load_lead_fc
from smj.phonon.checks import MAX_ASR_RELATIVE
from smj.phonon.negf import transmission
from smj.phonon.units import EV_TO_CM1

run = Path("runs/nap-sh-thermal/nap-sh-thermal")
meta = json.loads((run / "thermal_device/device_metadata.json").read_text())
k_dev = np.load(run / "thermal_fc/force_constants.npy")
lead = load_lead_fc("runs/lead-fc/lead")
masses = masses_from_elements(
    [line.split()[-1] for line in
     (run / "thermal_fc/geometry.in").read_text(encoding="utf-8").splitlines()
     if line.strip().startswith("atom")])

full = assemble_full_fc(
    k_dev, lead,
    FullFCMeta(len(meta["layer_index"]), meta["device_indices"],
               left_indices=meta["left_plot"] if "left_plot" in meta else
               [i for i, (s, l) in enumerate(zip(meta["side"], meta["layer_index"]))
                if s == "left" and l in (0, 1)],
               right_indices=[i for i, (s, l) in enumerate(zip(meta["side"], meta["layer_index"]))
                              if s == "right" and l in (0, 1)]))
D_raw = dynamical_matrix(full, masses)

layer, side = meta["layer_index"], meta["side"]
def dof(atoms):
    return np.array([3 * a + k for a in atoms for k in range(3)])

centre = dof(meta["device_indices"])
asr = asr_diagnostics(D_raw, masses, threshold=MAX_ASR_RELATIVE, rows=centre)
D = asr["dynamical"]
print("ASR: 绝对 %.4e 相对 %.4e 是否投影 %s 投影范数 %.3e"
      % (asr["residual_ev_per_A2"], asr["relative_residual"],
         asr["projected"], asr["projection_norm"]))
left_pl = dof([i for i in range(len(layer)) if side[i] == "left" and layer[i] == 1])
left_out = dof([i for i in range(len(layer)) if side[i] == "left" and layer[i] == 0])
right_pl = dof([i for i in range(len(layer)) if side[i] == "right" and layer[i] == 1])
right_out = dof([i for i in range(len(layer)) if side[i] == "right" and layer[i] == 0])

d_c = D[np.ix_(centre, centre)]
d_l, d_r = D[np.ix_(left_pl, left_pl)], D[np.ix_(right_pl, right_pl)]
v_l, v_r = D[np.ix_(centre, left_pl)], D[np.ix_(centre, right_pl)]
v_ll, v_lr = D[np.ix_(left_pl, left_out)], D[np.ix_(right_pl, right_out)]

# lead 模型(层1 为主层、层0 为下一层)的带顶:扫 k
kp = np.linspace(0, np.pi, 121)
top = 0.0
for k in kp:
    blk = d_l + v_ll * np.exp(1j * k) + v_ll.T * np.exp(-1j * k)
    f = mode_frequencies_cm1(np.real(blk))
    top = max(top, f.max())
print("lead 模型带顶: %.1f cm-1 = %.3f meV" % (top, top / EV_TO_CM1 * 1000))

# 器件中心区本征模的最高频
print("器件区最高模: %.1f cm-1" % mode_frequencies_cm1(d_c).max())

eta = 0.5e-3
w = np.arange(0.0, 24.4e-3, 1e-4)  # 0.1 meV 步长
T = np.array([transmission(x, d_c, d_l, d_r, v_l, v_r, eta ** 2,
                           v_lead_l=v_ll, v_lead_r=v_lr) for x in w])
above = w > 19.5e-3
print("带顶之上(>19.5 meV): max T = %.3e  条数=%d" % (T[above].max(), above.sum()))
print("  T>1e-6 的最高能量: %.3f meV" % (w[above][T[above] > 1e-6].max() * 1000)
      if (T[above] > 1e-6).any() else "  none")
print("  T>1e-6 的最高能量: n/a" if not (T[above] > 1e-6).any() else "")
for e in (19.0, 19.4, 19.5, 19.55, 19.6, 19.7, 20.0):
    i = int(np.argmin(np.abs(w - e * 1e-3)))
    print("   T(%.2f meV) = %.3e" % (w[i] * 1000, T[i]))
print("整条谱 max T = %.4f @ %.2f meV" % (T.max(), w[np.argmax(T)] * 1000))
print("虚频:", imaginary_modes(D)["count"])
