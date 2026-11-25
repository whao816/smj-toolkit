"""一次性诊断：检查已完成的位移单点是否都正常收敛（只读，不改任何东西）。"""
import json
import re
import statistics
import sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")

cyc = re.compile(r"Number of self-consistency cycles\s*:\s*(\d+)")
row = re.compile(r"^\s*\|\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")

rows = []
bad = []
warn = []
for d in sorted(root.glob("fc/[0-9][0-9][0-9][0-9]")):
    out = d / "aims.out"
    if not out.exists():
        bad.append((d.name, "no aims.out"))
        continue
    text = out.read_text(errors="replace")
    ok = "Have a nice day." in text
    cycles = [int(x) for x in cyc.findall(text)]
    n_conv = text.count("Self-consistency cycle converged")
    n_notconv = text.count("Self-consistency cycle NOT converged")

    # 最后一个力块的 max|F|
    fmax = None
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if "Total atomic forces" in line:
            start = i
    if start is not None:
        vals = []
        for line in lines[start + 1:]:
            m = row.match(line)
            if not m:
                break
            vals.append(max(abs(float(m.group(k).replace("D", "E")))
                            for k in (2, 3, 4)))
        if vals:
            fmax = max(vals)

    for line in text.splitlines():
        low = line.lower()
        if "not converged" in low or "error" in low and "error in hartree" not in low:
            warn.append((d.name, line.strip()[:140]))

    if not ok:
        bad.append((d.name, "no 'Have a nice day'"))
    if n_notconv:
        bad.append((d.name, "SCF NOT converged"))
    rows.append((d.name, ok, cycles[-1] if cycles else None, n_conv, fmax,
                 (d / "forces.json").exists()))

done = [r for r in rows if r[1]]
print("== 位移目录数：%d，正常结束：%d" % (len(rows), len(done)))
c = [r[2] for r in rows if r[2] is not None]
if c:
    print("== SCF 圈数：min %d median %.0f max %d（全部单次收敛：%s）"
          % (min(c), statistics.median(c), max(c),
             all(r[3] >= 1 for r in rows if r[2] is not None)))
f = [r[4] for r in done if r[4] is not None]
if f:
    print("== 末个力块 max|F|：min %.4f median %.4f max %.4f eV/A"
          % (min(f), statistics.median(f), max(f)))
print("== forces.json 齐备：%d/%d" % (sum(1 for r in rows if r[5]), len(rows)))
if bad:
    print("== 异常：")
    for name, why in bad:
        print("   %s: %s" % (name, why))
seen = set()
odd = []
for name, line in warn:
    k = line[:60]
    if k in seen:
        continue
    seen.add(k)
    odd.append((name, line))
print("== 可疑行（去重后 %d 类）：" % len(odd))
for name, line in odd[:10]:
    print("   %s: %s" % (name, line))
