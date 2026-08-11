#!/usr/bin/env python3
"""Final audit: recompute every derived number in verification.md from the raw
JSON by an independent route, and re-verify the theory claims the report makes.
Exits non-zero if anything disagrees.
"""
import json, math, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "..", "data", "verification.json")))
txt = open(os.path.join(HERE, "..", "verification.md"), encoding="utf-8").read()
bad = []
def check(name, ok, detail=""):
    print(("  OK   " if ok else "  FAIL ") + name + ("  " + detail if detail else ""))
    if not ok: bad.append(name)

print("[1] 报告中的数字是否都能在 JSON 中找到")
s = d["setup"]
check("lambda_c", f"{s['lambda_c']:.5f}" in txt, f"{s['lambda_c']:.5f}")
check("lambda", f"{s['lambda']:.5f}" in txt, f"{s['lambda']:.5f}")
for r in d["scaling"]:
    check(f"scaling N={r['N']}", f"{r['closure']:.5f}" in txt and f"{r['meanfield']:.5f}" in txt)
for r in d["threshold_crosscheck"]:
    check(f"threshold {r['case'][:18]}", f"{r['jacobian']:.8f}" in txt)

print("\n[2] 独立重算报告中的导出量")
inv = d["invariants"]
ratios = [inv[i]["identity"]/inv[i+1]["identity"] for i in range(len(inv)-1)]
check("RK4 四阶：三个比值均在 15.5-16.5", all(15.5 < r < 16.5 for r in ratios),
      " ".join(f"{r:.2f}" for r in ratios))
check("求和规则残差 < 1e-18", all(r["sumrule"] < 1e-18 for r in inv))
zs = [(r["recursion"]-r["mc"])/r["mc_se"] for r in d["cascade"]]
mz = sum(zs)/len(zs)
check("级联 z 均值与零相容 (<2 sigma)", abs(mz)*math.sqrt(len(zs)) < 2,
      f"mean z={mz:+.3f} -> {abs(mz)*math.sqrt(len(zs)):.2f} sigma")
check("级联 |z|>3 的格数为 0", sum(1 for z in zs if abs(z) > 3) == 0)
# A whole-dict write by a concurrently running script once reverted this block
# to a superseded version. The two superseded versions are identifiable: the
# standard error was a Poisson approximation, and the seeds were sequential.
poisson = [abs(r["mc_se"] - 2*math.sqrt(r["mc"]/400000)) < 1e-9 for r in d["cascade"]]
check("级联 SE 用样本方差而非 Poisson 近似", not any(poisson),
      f"{sum(poisson)}/{len(poisson)} 格匹配 Poisson 形式")
check("m=2 时 C 恰为 T", all(abs(r["recursion"]-r["naive"]) < 1e-12
                            for r in d["cascade"] if r["m"] == 2))
check("m>=3 时 C > (m-1)T", all(r["recursion"] > r["naive"]
                                for r in d["cascade"] if r["m"] >= 3))
sc = d["scaling"]
check("闭合偏差单调下降（在 CI 内）",
      all(sc[i]["closure_lo"] <= sc[i-1]["closure_hi"] for i in range(1, len(sc))))
check("最大 N 处闭合落入零假设水平", sc[-1]["closure"] < sc[-1]["noise"],
      f"{sc[-1]['closure']:.5f} < {sc[-1]['noise']:.5f}")
check("均场在最大 N 处远高于零假设水平", sc[-1]["meanfield"] > 50*sc[-1]["noise"],
      f"{sc[-1]['meanfield']/sc[-1]['noise']:.0f}x")
check("每个 N 的 CI 都包含点估计",
      all(r["closure_lo"] <= r["closure"] <= r["closure_hi"] and
          r["meanfield_lo"] <= r["meanfield"] <= r["meanfield_hi"] for r in sc))

print("\n[3] 重新验证理论断言（不读 JSON，直接算）")
from audit_theory import Nmat, rho, C
from verify_closure import Closure
import numpy as np
P = {(2, 2): .5, (3, 3): .5}; MS = tuple(s["ms"]); W = tuple(s["w"])
lo, hi = 1e-4, 6.0
for _ in range(70):
    m = .5*(lo+hi)
    if rho(Nmat(P, MS, [m*x for x in W], (1, 1))) > 1: hi = m
    else: lo = m
check("lambda_c 可复现", abs(.5*(lo+hi)-s["lambda_c"])/s["lambda_c"] < 1e-6)
A = Nmat(P, MS, [s["lambda"]*x for x in W], (1, 1))
check("N 逐元非负（命题 4 的前提）", all(A[i][j] >= 0 for i in range(2) for j in range(2)))
check("rho(N) >= max 对角元（命题 4）", rho(A) >= max(A[0][0], A[1][1]) - 1e-12,
      f"rho={rho(A):.4f} vs max diag {max(A[0][0], A[1][1]):.4f}")
check("超临界（lambda=1.6 lambda_c）", rho(A) > 1, f"rho={rho(A):.4f}")
c2 = Closure(P, MS, [s["lambda"]*x for x in W])
idn, sm = c2.invariants(s["eps"], 12.0, 0.005)
check("恒等式在本文参数下成立", idn < 1e-11, f"{idn:.2e}")
check("求和规则在本文参数下成立", sm < 1e-18, f"{sm:.2e}")

print("\n[4] 图与数据是否同源")
for f in ("fig1_trajectory", "fig2_verification"):
    for ext in ("pdf", "png"):
        p = os.path.join(HERE, "..", "figures", f+"."+ext)
        check(f"{f}.{ext} 存在", os.path.exists(p))
jt = os.path.getmtime(os.path.join(HERE, "..", "data", "verification.json"))
for f in ("fig1_trajectory.png", "fig2_verification.png"):
    check(f"{f} 不早于数据",
          os.path.getmtime(os.path.join(HERE, "..", "figures", f)) >= jt)
check("verification.md 不早于数据",
      os.path.getmtime(os.path.join(HERE, "..", "verification.md")) >= jt)

print()
if bad:
    print(f"FAILED: {len(bad)} 项 -> {bad}"); sys.exit(1)
print("全部通过")
