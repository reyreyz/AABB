# Two-level scan, constant-wind centered CD, k = 4, q = 4

| Pe | delta | m | nc/n | rho(K) Q0 | rho(K_Z) | M | N (Q0) | N_Z |
|---|---|---|---|---|---|---|---|---|
| 4 | 0.6 | 31 | 0.48 | 0.486 | 0.447 | **0.67** | **0.00126** | **0.0159** |
| 4 | 0.6 | 63 | 0.49 | 0.495 | 0.450 | 1.42 | **0.0968** | **0.0611** |
| 4 | 0.6 | 127 | 0.50 | 0.498 | 0.448 | 1.68 | **0.199** | **0.0613** |
| 4 | 1.0 | 31 | 0.47 | 0.759 | 0.758 | **0.893** | **0.308** | **0.308** |
| 4 | 1.0 | 63 | 0.48 | 0.750 | 0.749 | 3.32 | **0.482** | **0.483** |
| 4 | 1.0 | 127 | 0.49 | 0.740 | 0.738 | 3.3 | **0.538** | **0.539** |
TypeError: float() argument must be a string or a real number, not 'csr_matrix'

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/home/user/AABB/rbamg/scan_o_n.py", line 46, in <module>
    main()
    ~~~~^^
  File "/home/user/AABB/rbamg/scan_o_n.py", line 35, in main
    c = mzz.Alg73Z(A, k=k, delta=dl, max_levels=2, **kw)
  File "/home/user/AABB/rbamg/ml_z73.py", line 115, in __init__
    Afc, Aff = A[F][:, C].tocsr(), A[F][:, F].tocsr()
            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/user/AABB/rbamg/ml_nz.py", line 140, in z_kk
    res = minimize(fg, z0, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter))
  File "/root/.local/lib/python3.13/site-packages/scipy/optimize/_minimize.py", line 784, in minimize
    res = _minimize_lbfgsb(fun, x0, args, jac, bounds,
                           callback=callback, **options)
  File "/root/.local/lib/python3.13/site-packages/scipy/optimize/_lbfgsb_py.py", line 362, in _minimize_lbfgsb
    sf = _prepare_scalar_function(fun, x0, jac=jac, args=args, epsilon=eps,
                                  bounds=bounds,
                                  finite_diff_rel_step=finite_diff_rel_step,
                                  workers=workers)
  File "/root/.local/lib/python3.13/site-packages/scipy/optimize/_optimize.py", line 310, in _prepare_scalar_function
    sf = ScalarFunction(fun, x0, args, grad, hess,
                        finite_diff_rel_step, bounds, epsilon=epsilon,
                        workers=workers)
  File "/root/.local/lib/python3.13/site-packages/scipy/optimize/_differentiable_functions.py", line 252, in __init__
    self.x = xp.astype(_x, _dtype)
             ~~~~~~~~~^^^^^^^^^^^^
  File "/root/.local/lib/python3.13/site-packages/scipy/_external/array_api_compat/numpy/_aliases.py", line 106, in astype
    return x.astype(dtype=dtype, copy=copy)
           ~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^
ValueError: setting an array element with a sequence.
