"""SVG small multiples of the eps-pseudospectra of E (|A|-norm) from psE.json. No matplotlib needed."""
import json, os, sys
import numpy as np

D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "psE.json")))
NORM = sys.argv[1] if len(sys.argv) > 1 else "|A|"
KH = ["0.0", "0.3", "0.6", "0.9"]
PS = ["linear", "rho-opt (w=3)", "norm-opt (w=3)", "eigvec (dense)"]
TITLES = ["linear (hat)", "ρ-optimal, w=3", "‖E‖₂-optimal, w=3", "eigenvector space (dense)"]
LEVELS = [-0.5, -1.0, -1.5, -2.0]
COLORS = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]       # ordinal blue ramp, outer (large eps) light
INK, MUTED, GRID, ACC = "#1f1f1e", "#6b6a66", "#d9d8d4", "#e0721f"

def segments(Z, X, Y, lev):
    """Marching squares on Z[y, x]; returns list of ((x0,y0),(x1,y1)) in data coords."""
    segs = []
    for i in range(len(Y) - 1):
        for j in range(len(X) - 1):
            c = [Z[i, j], Z[i, j + 1], Z[i + 1, j + 1], Z[i + 1, j]]
            p = [(X[j], Y[i]), (X[j + 1], Y[i]), (X[j + 1], Y[i + 1]), (X[j], Y[i + 1])]
            pts = []
            for k in range(4):
                a, b = c[k], c[(k + 1) % 4]
                if (a < lev) != (b < lev):
                    t = (lev - a) / (b - a)
                    pts.append((p[k][0] + t * (p[(k + 1) % 4][0] - p[k][0]), p[k][1] + t * (p[(k + 1) % 4][1] - p[k][1])))
            if len(pts) == 2: segs.append(tuple(pts))
            elif len(pts) == 4: segs += [(pts[0], pts[1]), (pts[2], pts[3])]
    return segs

W, H, M = 200, 200, 14                     # panel size, inner margin
LEFT, TOP, GAPX, GAPY = 70, 106, 22, 46
fw = LEFT + 4 * W + 3 * GAPX + 20
fh = TOP + 4 * H + 3 * GAPY + 30
o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{fw}" height="{fh}" viewBox="0 0 {fw} {fh}" font-family="Helvetica, Arial, sans-serif">',
     f'<rect width="100%" height="100%" fill="#ffffff"/>',
     f'<text x="{LEFT}" y="24" font-size="16" font-weight="600" fill="{INK}">ε-pseudospectra of the two-grid operator E ({NORM}-norm), 1D Helmholtz, n = 31, damped Jacobi ω = 2/3</text>']
# legend
lx = LEFT
for lev, col in zip(LEVELS, COLORS):
    o.append(f'<line x1="{lx}" y1="48" x2="{lx+22}" y2="48" stroke="{col}" stroke-width="2"/>')
    o.append(f'<text x="{lx+28}" y="52" font-size="12" fill="{INK}">ε = 10<tspan dy="-5" font-size="9">{lev:g}</tspan></text>'); lx += 100
o.append(f'<circle cx="{lx+6}" cy="48" r="3" fill="{INK}"/><text x="{lx+16}" y="52" font-size="12" fill="{INK}">eigenvalues</text>'); lx += 100
o.append(f'<line x1="{lx}" y1="48" x2="{lx+22}" y2="48" stroke="{ACC}" stroke-width="1.5" stroke-dasharray="4 3"/><text x="{lx+28}" y="52" font-size="12" fill="{INK}">|z| = ρ(E)</text>'); lx += 100
o.append(f'<line x1="{lx}" y1="48" x2="{lx+22}" y2="48" stroke="{MUTED}" stroke-width="1" stroke-dasharray="2 3"/><text x="{lx+28}" y="52" font-size="12" fill="{INK}">|z| = 1</text>')
o.append(f'<text x="{LEFT}" y="72" font-size="11" fill="{MUTED}">Each panel has its own scale (half-width R shown). E is real, so only the upper half-plane was computed and it is mirrored.</text>')
for c, t in enumerate(TITLES):
    o.append(f'<text x="{LEFT + c*(W+GAPX) + W/2}" y="{TOP-8}" font-size="13" font-weight="600" text-anchor="middle" fill="{INK}">{t}</text>')
for r, kh in enumerate(KH):
    y0 = TOP + r * (H + GAPY)
    o.append(f'<text x="{LEFT-14}" y="{y0+H/2}" font-size="13" font-weight="600" text-anchor="end" fill="{INK}">kh = {kh}</text>')
    for c, pn in enumerate(PS):
        d = D[kh][pn][NORM]; x0 = LEFT + c * (W + GAPX)
        R = d["R"]; X = np.array(d["X"]); Y = np.array(d["Y"]); Z = np.array(d["logs"])
        Zf = np.vstack([Z[::-1][:-1], Z]); Yf = np.concatenate([-Y[::-1][:-1], Y])
        sx = lambda x: x0 + M + (x + R) / (2 * R) * (W - 2 * M)
        sy = lambda y: y0 + M + (R - y) / (2 * R) * (H - 2 * M)
        o.append(f'<rect x="{x0}" y="{y0}" width="{W}" height="{H}" fill="#fbfbfa" stroke="{GRID}"/>')
        o.append(f'<line x1="{sx(-R)}" y1="{sy(0)}" x2="{sx(R)}" y2="{sy(0)}" stroke="{GRID}" stroke-width="1"/>')
        o.append(f'<line x1="{sx(0)}" y1="{sy(-R)}" x2="{sx(0)}" y2="{sy(R)}" stroke="{GRID}" stroke-width="1"/>')
        o.append(f'<defs><clipPath id="c{r}{c}"><rect x="{x0}" y="{y0}" width="{W}" height="{H}"/></clipPath></defs><g clip-path="url(#c{r}{c})">')
        k = (W - 2 * M) / (2 * R)
        o.append(f'<circle cx="{sx(0)}" cy="{sy(0)}" r="{k:.2f}" fill="none" stroke="{MUTED}" stroke-width="1" stroke-dasharray="2 3"/>')
        for lev, col in zip(LEVELS, COLORS):
            path = "".join(f"M{sx(a[0]):.1f},{sy(a[1]):.1f}L{sx(b[0]):.1f},{sy(b[1]):.1f}" for a, b in segments(Zf, X, Yf, lev))
            if path: o.append(f'<path d="{path}" fill="none" stroke="{col}" stroke-width="1.6" stroke-linecap="round"/>')
        o.append(f'<circle cx="{sx(0)}" cy="{sy(0)}" r="{d["rho"]*k:.2f}" fill="none" stroke="{ACC}" stroke-width="1.5" stroke-dasharray="4 3"/>')
        for ex, ey in d["ev"]:
            o.append(f'<circle cx="{sx(ex):.1f}" cy="{sy(ey):.1f}" r="2.2" fill="{INK}" stroke="#ffffff" stroke-width="0.8"/>')
        o.append('</g>')
        o.append(f'<text x="{x0+6}" y="{y0+H+14}" font-size="11" fill="{MUTED}">R={R:.2g}</text>')
        o.append(f'<text x="{x0+W-4}" y="{y0+H+14}" font-size="11" text-anchor="end" fill="{INK}">ρ={d["rho"]:.3f}  ‖E‖={d["normE"]:.2f}</text>')
        o.append(f'<text x="{x0+W-4}" y="{y0+H+28}" font-size="11" text-anchor="end" fill="{MUTED}">α(0.1)={d["alpha"][0]:.2f}  α(0.01)={d["alpha"][1]:.3f}</text>')
o.append('</svg>')
name = "psE_absA.svg" if NORM == "|A|" else "psE_2norm.svg"
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), name), "w").write("\n".join(o))
print("wrote", name)
