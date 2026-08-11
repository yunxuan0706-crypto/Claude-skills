#!/usr/bin/env python3
"""Is the 2.5% gap between the closure's own final-size bisection and
rho(N)=1 an integration-window artefact?

Near threshold the growth rate vanishes, so a finite t_max fails to see growth
and biases the measured threshold upward -- which is the observed sign. If that
is the whole story the gap must fall monotonically with t_max.
"""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from verify_closure import Closure
from audit_theory import Nmat, rho

P = {(2, 4): .5, (4, 2): .5}; ms = (3, 4); w = (1.0, 1.0); th = (1, 1)
lo, hi = 1e-3, 4.0
for _ in range(60):
    m = .5*(lo+hi)
    if rho(Nmat(P, ms, [m*x for x in w], th)) > 1: hi = m
    else: lo = m
pred = .5*(lo+hi)
print(f"rho(N)=1 prediction: {pred:.6f}")
print(f"{'t_max':>7} {'ODE threshold':>15} {'rel. deviation':>15}")
for tmax in (100.0, 200.0, 400.0, 800.0, 1600.0):
    lo, hi = 1e-3, 1.0
    for _ in range(22):
        mid = .5*(lo+hi)
        c = Closure(P, ms, [mid*x for x in w], list(th))
        _, y = c.run(1e-9, tmax, 0.05)
        if y[-1] > 1e-3: hi = mid
        else: lo = mid
    m = .5*(lo+hi)
    print(f"{tmax:7.0f} {m:15.6f} {abs(m-pred)/pred:15.2e}")
