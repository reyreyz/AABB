"""Build a first-difference (bidiagonal) operator on a uniform grid.

The matrix has `a` on the main diagonal and `b` on the first superdiagonal:

    [ a  b  0  ...  0 ]
    [ 0  a  b  ...  0 ]
    [ ...             ]
    [ 0  0  ...  0  a ]
"""

import numpy as np
import matplotlib.pyplot as plt

# --- grid -------------------------------------------------------------------
h = 0.001                 # grid spacing
L = 3 * np.pi             # domain length
n = int(round(L / h))     # number of grid points (9425)

x = np.arange(n) * h      # spacing is exactly h; linspace(0, L, n) would not be

a = 1.0
b = -1.0


# --- operator ---------------------------------------------------------------
def build_dense(n, a, b):
    """Dense n x n version. Costs n**2 * 8 bytes (~710 MB at n=9425), so only
    call this for small n."""
    m = np.zeros((n, n))
    np.fill_diagonal(m, a)            # in place: returns None, don't assign it
    np.fill_diagonal(m[:-1, 1:], b)   # view onto the first superdiagonal
    return m


# Equivalent allocation-free representation: just the two nonzero bands.
main = np.full(n, a)
upper = np.full(n - 1, b)

print(f"n = {n}, dense size would be {n * n * 8 / 1e9:.2f} GB")
print("main diagonal :", main[:5], "...")
print("superdiagonal :", upper[:5], "...")

# Small dense version, safe to print and display.
D_small = build_dense(8, a, b)
print("\nD (n=8):\n", D_small)

# Applying the operator without ever forming the matrix:
#   (D @ u)[i] = a*u[i] + b*u[i+1],  and (D @ u)[n-1] = a*u[n-1]
u = np.sin(x)
Du = a * u.copy()
Du[:-1] += b * u[1:]

# Divide by h to get an actual finite difference. With a=1, b=-1 this is
# -(u[i+1] - u[i])/h, so negate it to get the forward difference du/dx.
# The last entry has no u[i+1] to pair with, so it is not a derivative: drop it.
deriv = -Du[:-1] / h

# --- plots ------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

ax1.imshow(D_small, cmap="coolwarm")
ax1.set_title(f"D structure (n=8, a={a}, b={b})")

ax2.plot(x, u, label="u = sin(x)")
ax2.plot(x[:-1], deriv, "--", label="forward difference ≈ cos(x)")
ax2.set_xlabel("x")
ax2.legend()
ax2.set_title("operator applied to sin(x)")

fig.tight_layout()
plt.show()
