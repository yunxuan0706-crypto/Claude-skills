"""Is the simulation really alpha-independent at fixed bi-degree sequence?

Section III.C says the tree-like closure cannot see the overlap alpha, and
offers that as a falsifiable prediction. The closure side is trivial -- it only
reads the degree sequence. The simulation side has never been tested.
"""
import math, random, statistics as st, sys
LOG=open('/tmp/claude-0/-home-user-Claude-skills/4b148b58-a40b-5105-80bb-fdabe381c458/scratchpad/alpha_result.txt','w',buffering=1)
sys.path.insert(0,'tools')
from ebcm_directed import PGF, bidegree, config_model, degrees, ebcm, gillespie_traj
from feasibility_check import index

def alphas(E, tau, eta):
    TT=HH=0; mTT=mHH=0
    for i in range(len(E)):
        for j in range(len(E)):
            if i==j: continue
            c=len(E[i][0]&E[j][0])
            if c: TT+=c/tau; mTT+=1
            c=len(E[i][1]&E[j][1])
            if c: HH+=c/eta; mHH+=1
    return (TT/mTT if mTT else 0), (HH/mHH if mHH else 0)

def drive(E, rng, n, sign):
    """degree-preserving double swaps that raise (sign=+1) or lower alpha_par"""
    cur = sum(alphas(E,2,2))/2
    for _ in range(n):
        i,j = rng.sample(range(len(E)),2)
        A,B = (E[i][0],E[j][0]) if rng.random()<0.5 else (E[i][1],E[j][1])
        v1,v2 = rng.choice(sorted(A)), rng.choice(sorted(B))
        if v1==v2 or v1 in E[j][0]|E[j][1] or v2 in E[i][0]|E[i][1]: continue
        A.remove(v1); A.add(v2); B.remove(v2); B.add(v1)
        new = sum(alphas(E,2,2))/2
        if sign*(new-cur) > 0: cur = new
        else:
            A.remove(v2); A.add(v1); B.remove(v1); B.add(v2)
    return cur

N, tau, eta, runs, eps = 500, 2, 2, 400, 0.006
rng = random.Random(555)
kin, kout = bidegree(N, tau, eta, 2, 0.0, rng)
base = config_model(kin, kout, tau, eta, rng)
ki, ko = degrees(N, base); p = PGF(ki, ko)
lc = 1.0/(tau*p.kappa-1)
print(f"N={N} kappa={p.kappa:.4f} lambda_c={lc:.4f}  (bi-degree sequence identical throughout)")

out=[]
for tag, sign, nsw in (("low  alpha", -1, 3000), ("base     ", 0, 0), ("high alpha", +1, 3000)):
    E=[(set(T),set(H)) for T,H in base]
    if nsw: drive(E, random.Random(9+nsw+sign), nsw, sign)
    ki2,ko2 = degrees(N,E)
    assert ki2==ki and ko2==ko, "degree sequence changed!"
    aTT,aHH = alphas(E,tau,eta); apar=(aTT+aHH)/2
    tails,heads,tails_of = index(N,E)
    for r in (1.0, 1.6):
        lam=r*lc
        vals=[]
        for s in range(runs):
            rr=random.Random(31000+77*s+int(100*r)+nsw*sign)
            sd=[v for v in range(N) if rr.random()<eps]
            vals.append(gillespie_traj(N,tails,heads,tails_of,lam,1,rr,sd,[400.0])[-1][2])
        m=st.mean(vals); se=math.sqrt(st.variance(vals)/runs)
        e=ebcm(p,tau,lam,eps=eps,tmax=250.0,dt=0.01)[-1][3]
        out.append((tag,apar,r,m,se,e))
        print(f"  {tag}  alpha_par={apar:.4f}  lam={r:.1f}lc:  sim R={m:.5f}+/-{se:.5f}   EBCM={e:.5f}")
print()
for r in (1.0,1.6):
    rows=[o for o in out if o[2]==r]
    lo,hi = rows[0], rows[-1]
    d=hi[3]-lo[3]; sd=math.sqrt(lo[4]**2+hi[4]**2)
    print(f"lam={r:.1f}lc: alpha_par {lo[1]:.4f} -> {hi[1]:.4f}   sim 变化 {d:+.5f} +/- {sd:.5f}  ({abs(d)/sd:.2f} sigma, 相对 {100*d/lo[3]:+.1f}%)")
