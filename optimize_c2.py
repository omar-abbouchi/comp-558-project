import numpy as np
from scipy.optimize import minimize
import time, cirq

N = 256
T2 = np.exp(-0.5*(np.arange(N)-N//2)**2/(N//4)**2); T2 /= T2.sum()

def cnot_mat(ctrl, tgt, n=8):
    dim = 2**n; mat = np.eye(dim, dtype=complex)
    for i in range(dim):
        bits = [(i >> (n-1-q)) & 1 for q in range(n)]
        if bits[ctrl] == 1:
            bits[tgt] ^= 1
            j = sum(b << (n-1-q) for q, b in enumerate(bits))
            mat[i,i] = 0; mat[j,i] = 1
    return mat

def apply_ry(sv, params, offset):
    sv = sv.reshape([2]*8)
    for qi in range(8):
        t = params[offset+qi]; c, s = np.cos(t/2), np.sin(t/2)
        i0 = tuple(0 if k==qi else slice(None) for k in range(8))
        i1 = tuple(1 if k==qi else slice(None) for k in range(8))
        a, b = sv[i0].copy(), sv[i1].copy()
        sv[i0] = c*a - s*b; sv[i1] = s*a + c*b
    return sv.reshape(256)

def make_objective(cnot_pairs):
    chain = np.eye(256, dtype=complex)
    for c, t in cnot_pairs: chain = cnot_mat(c,t) @ chain
    def run(p):
        sv = np.zeros(256, dtype=complex); sv[0] = 1.0
        sv = apply_ry(sv,p,0); sv = chain@sv; sv = apply_ry(sv,p,8)
        return 0.5*np.sum(np.abs(np.abs(sv)**2 - T2))
    def f(params):
        v = run(params); g = np.zeros(16)
        for i in range(16):
            pp=params.copy(); pp[i]+=np.pi/2
            pm=params.copy(); pm[i]-=np.pi/2
            g[i] = (run(pp) - run(pm)) / 2
        return v, g
    return f

configs = [
    ([(0,1),(0,2),(0,3)],    500),
    ([(0,1),(0,2)],           400),
    ([(0,1),(0,2),(1,3)],    400),
    ([(0,1),(1,2),(2,3)],    400),
    ([(0,1),(0,2),(0,3),(0,4)], 300),
]

best_g = {'tvd': 999, 'params': None, 'cnots': None}

for cnots, n_restarts in configs:
    n2=len(cnots); n1=16; d=n2+3
    f = make_objective(cnots)
    best=999; bp=None; t0=time.time()
    print(f"\npattern={cnots}, {n_restarts} restarts", flush=True)
    for i in range(n_restarts):
        np.random.seed(i)
        x0 = np.random.uniform(0, 2*np.pi, 16)
        r = minimize(f, x0, method='L-BFGS-B', jac=True,
                     options={'maxiter':500,'ftol':1e-13,'gtol':1e-9})
        if r.fun < best: best=r.fun; bp=list(r.x)
        if (i+1) % 50 == 0 or (i+1) == n_restarts:
            elapsed=time.time()-t0; eta=(n_restarts-i-1)/max((i+1)/elapsed,1e-9)
            s=1/(0.75*best+0.05*n1/100+0.1*n2/100+0.1*d/100)
            print(f"  [{i+1:3d}/{n_restarts}] TVD={best:.5f} score={s:.2f} | {elapsed:.0f}s ~{eta:.0f}s left", flush=True)
    if best < best_g['tvd']:
        best_g = {'tvd':best,'params':bp,'cnots':cnots}
        print("  *** NEW BEST ***", flush=True)

# Polish
f = make_objective(best_g['cnots'])
r = minimize(f, np.array(best_g['params']), method='L-BFGS-B', jac=True,
             options={'maxiter':50000,'ftol':1e-15,'gtol':1e-12})
n2=len(best_g['cnots'])
print(f"\nFINAL: TVD={r.fun:.8f} score_est={1/(0.75*r.fun+0.05*16/100+0.1*n2/100+0.1*(n2+3)/100):.2f}")
print(f"CNOTS: {best_g['cnots']}")
print(f"PARAMS: {list(r.x)}")

# Verify
q=cirq.LineQubit.range(8); sim=cirq.Simulator(); scores=[]
for run in range(3):
    ops=[cirq.ry(float(r.x[i]))(q[i]) for i in range(8)]
    for c,t in best_g['cnots']: ops.append(cirq.CX(q[c],q[t]))
    ops+=[cirq.ry(float(r.x[8+i]))(q[i]) for i in range(8)]
    ops.append(cirq.measure(*q,key='m'))
    c2=cirq.Circuit(ops)
    res=sim.run(c2,repetitions=100000)
    counts=res.histogram(key='m')
    p_out=np.zeros(256)
    for k,v in counts.items(): p_out[k]=v/100000
    tvd_s=0.5*np.sum(np.abs(p_out-T2))
    all_ops=[op for op in c2.all_operations() if not isinstance(op.gate,cirq.MeasurementGate)]
    n1s=sum(1 for op in all_ops if len(op.qubits)==1)
    n2s=sum(1 for op in all_ops if len(op.qubits)==2)
    depth=len(cirq.Circuit(c2.all_operations()))
    sc=1/(0.75*tvd_s+0.05*n1s/100+0.1*n2s/100+0.1*depth/100)
    scores.append(sc)
    print(f"  run {run+1}: TVD={tvd_s:.4f} 1q={n1s} 2q={n2s} depth={depth} score={sc:.2f}")
print(f"Mean sampled score: {np.mean(scores):.2f}")