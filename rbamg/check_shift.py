"""E_j e0 = P~(I - B_c A_c) c_j + (I - P~ B_c R A) Q[(S - S~) c_k + K^k q0]  -- check on 2 levels."""
import numpy as np
from ml_twins import ML
from robust_n import recirc
rng = np.random.default_rng(1)
A = recirc(63, 2.0); n = A.shape[0]; w = 0.8
for tau in (0.0, 0.5):
    for sh in (0.0, 0.5, 1.0):
        ml = ML(A, omega=w, tauP=tau, mode="lump", max_levels=2, tau_filter=1e-2, shift=sh)
        L, Lc = ml.levels
        x0 = rng.standard_normal(n); b = np.zeros(n)          # x* = 0, e = -x
        xs = [x0]
        for _ in range(L.k): xs.append(xs[-1] + w * (b - A @ xs[-1]) / L.D)
        j = int(round(sh * L.k))
        cj, ck = -xs[j][L.C], -xs[-1][L.C]
        ek = -xs[-1]; qk_part = np.zeros(n)
        # Q[(S - S~) c_k + K^k q0]: F-part of eps_N minus S~ c_k  (eps_N = P' c_k + Q K^k q0)
        q0 = -x0[L.F] - 0                                      # Q0^sharp e0 = e0_F - 0 (P0 = [I;0]: Q^sharp = F rows - S? no: Q0^sharp = [0 I] - ?)
        Kq = q0.copy()
        for _ in range(L.k): Kq = L.K @ Kq
        # for P0/Q0: Psh = [I 0], Qsh = [0 I]; e0 = P0 c0 + Q0 q0 with q0 = e0_F
        v = np.zeros(n); v[L.F] = (L.Sex - L.S) @ ck + Kq
        Av = A @ v; ev = Lc.lu.solve(Av[L.C] + L.ST @ Av[L.F])
        v[L.C] -= ev; v[L.F] -= L.S @ ev
        u = cj - Lc.lu.solve(ml.RAP(L, cj))                     # (I - B_c A_c) c_j
        pred = v.copy(); pred[L.C] += u; pred[L.F] += L.S @ u
        meas = -ml.cycle(0, b, x0, "N")
        print(f"tau={tau}, j={j}/{L.k}: ‖E_j e0 - prediction‖/‖E_j e0‖ = {np.linalg.norm(meas-pred)/np.linalg.norm(meas):.1e}")
