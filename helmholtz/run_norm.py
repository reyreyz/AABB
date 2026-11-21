import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hh1d import *
n = int(sys.argv[1]); kh = float(sys.argv[2]); rng = np.random.default_rng(0)
A, S = helmholtz(n, kh)
for w in [int(x) for x in sys.argv[3].split(',')]:
    if w == 0:
        cs, _ = supports(n, 3); sups = [np.arange(n)] * len(cs); label = "dense"
        _, s0 = supports(n, 3); B, mask = waves_bc(n, kh, 1)
        P0, _ = constrained_min_masked(n, s0, cs, A @ A, B, mask)
    else:
        cs, sups = supports(n, w); label = f"w={w}"
        _, s3 = supports(n, 3); B, mask = waves_bc(n, kh, 1)
        P0, _ = constrained_min_masked(n, s3, cs, A @ A, B, mask)   # width-3 rho-optimum, embedded
    starts = [pack(P0, sups)]
    B, mask = waves_bc(n, kh, (w - 1) // 2 if w else 1)
    if w > 3:
        starts.append(pack(constrained_min_masked(n, sups, cs, A @ A, B, mask)[0], sups))
    starts += [starts[0] + 0.2 * rng.standard_normal(len(starts[0])) for _ in range(2)]
    t = time.time()
    v, P = true_opt_norm(A, S, sups, starts)
    E = twogrid(A, S, P)
    E0 = twogrid(A, S, P0)
    print(f"n={n} kh={kh} {label}: min ||E||2 = {v:.4f} (rho there {rho(E):.4f}, cond P {np.linalg.cond(P):.1f}); "
          f"rho-optimal w=3: ||E||2 = {np.linalg.norm(E0, 2):.4f}  [{time.time()-t:.0f}s]", flush=True)
    np.save(f"Pnorm_n{n}_kh{kh}_{label}.npy", P)
