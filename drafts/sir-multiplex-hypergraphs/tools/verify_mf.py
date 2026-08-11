"""Check Prop 3.4 (the mean-field matrix) and Prop 3.3 (synergy)."""
import numpy as np
from functools import lru_cache
def C(m, lam):
    @lru_cache(None)
    def u(i,s):
        if i==0 or s==0: return 0.0
        t=lam*s+i
        return (lam*s/t)*(1+u(i+1,s-1))+(i/t)*u(i-1,s)
    return u(1,m-1)

def mf_growth(P, ms, ws, lam):
    """numerically linearise the degree-based mean field, mu=1"""
    M=len(ms); betas=[lam*w for w in ws]
    ks=list(P); kbar=[sum(p*k[a] for k,p in P.items()) for a in range(M)]
    n=len(ks); eps=1e-7
    def rhs(i_vec):
        phi=[sum(P[k]*k[a]*i_vec[j] for j,k in enumerate(ks))/kbar[a] for a in range(M)]
        Th=[1-(1-phi[a])**(ms[a]-1) for a in range(M)]
        return np.array([sum(k[a]*betas[a]*Th[a] for a in range(M))*(1-i_vec[j]) - i_vec[j]
                         for j,k in enumerate(ks)])
    x0=np.zeros(n); f0=rhs(x0); J=np.zeros((n,n))
    for j in range(n):
        xp=x0.copy(); xp[j]+=eps; J[:,j]=(rhs(xp)-f0)/eps
    return max(np.linalg.eigvals(J).real)

def N_mf(P, ms, ws, lam):
    M=len(ms); kbar=[sum(p*k[a] for k,p in P.items()) for a in range(M)]
    A=np.zeros((M,M))
    for a in range(M):
        for b in range(M):
            kab=sum(p*k[a]*k[b] for k,p in P.items())
            A[a,b]=(ms[b]-1)*lam*ws[b]*kab/kbar[a]
    return A

def bisect(f, lo=1e-4, hi=5.0):
    for _ in range(80):
        m=.5*(lo+hi)
        if f(m)>0: hi=m
        else: lo=m
    return .5*(lo+hi)

P={(2,2):.5,(4,4):.5}; ms=(3,4); ws=(1.0,1.0)
a=bisect(lambda L: mf_growth(P,ms,ws,L))
b=bisect(lambda L: max(abs(np.linalg.eigvals(N_mf(P,ms,ws,L))))-1)
print(f"MF threshold: ODE Jacobian {a:.6f}   matrix (eq 3.4) {b:.6f}   rel {abs(a-b)/b:.1e}")

# Prop 3.4: sign of the net mean-field error across (m, lambda)
print("\nnet MF error at the true threshold  (lc_MF/lc_true - 1; <0 means MF underestimates lc)")
print(f"{'m1=m2':>6} {'lc_true':>9} {'lc_MF':>9} {'ratio':>8}")
for m in (2,3,4,6,8):
    Pu={(3,3):.5,(5,5):.5}; mm=(m,m)
    def Ntrue(L):
        kb=[sum(p*k[a] for k,p in Pu.items()) for a in range(2)]
        A=np.zeros((2,2))
        for x in range(2):
            for y in range(2):
                kxy=sum(p*k[x]*k[y] for k,p in Pu.items())
                A[x,y]=C(mm[y],L)*(kxy/kb[x]-(1 if x==y else 0))
        return A
    lt=bisect(lambda L: max(abs(np.linalg.eigvals(Ntrue(L))))-1)
    lm=bisect(lambda L: max(abs(np.linalg.eigvals(N_mf(Pu,mm,(1,1),L))))-1)
    print(f"{m:>6} {lt:9.5f} {lm:9.5f} {lm/lt:8.4f}")

# Prop 3.3: two subcritical layers, jointly supercritical?
print("\nsynergy check (Prop 3.3)")
Ps={(2,2):.5,(3,3):.5}; ms2=(3,3)
def N2(L, cross):
    kb=[sum(p*k[a] for k,p in Ps.items()) for a in range(2)]
    A=np.zeros((2,2))
    for x in range(2):
        for y in range(2):
            kxy=sum(p*k[x]*k[y] for k,p in Ps.items())
            A[x,y]=C(ms2[y],L)*(kxy/kb[x]-(1 if x==y else 0))
            if x!=y and not cross: A[x,y]=0.0
    return A
for L in (0.10, 0.13, 0.16):
    Aiso=N2(L, False); Aful=N2(L, True)
    print(f"  lam={L}: single-layer R = {Aiso[0,0]:.4f}, {Aiso[1,1]:.4f}"
          f"   |  multiplex rho = {max(abs(np.linalg.eigvals(Aful))):.4f}")
