import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hh1d import *
n = 31; kh = float(sys.argv[1])
A, S = helmholtz(n, kh)
cs, s7 = supports(n, 7)
P5 = np.load(f"Pnorm_n{n}_kh{kh}_w=5.npy")
v, P = true_opt_norm(A, S, s7, [pack(P5, s7)])
old = np.linalg.norm(twogrid(A, S, np.load(f"Pnorm_n{n}_kh{kh}_w=7.npy")), 2)
print(f"kh={kh} w=7 restarted from w=5 optimum: ||E||2 = {v:.4f} (rho {rho(twogrid(A,S,P)):.4f}, cond P {np.linalg.cond(P):.1f}); previous {old:.4f}")
if v < old: np.save(f"Pnorm_n{n}_kh{kh}_w=7.npy", P)
