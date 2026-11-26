# Section 6.3 reproduction (PyAMG components, own work count)
Recirculating flow: stand-in discretization (not identified).

## Table SM15 (fixed configurations): total WU to 1e-10, ours / paper (cycles, oc)

| problem | method | config | m = 127 | m = 255 | m = 511 |
|---|---|---|---|---|---|
| Advection | lAIR | 1pt lAIR2 -|FF theta=0.25 tauA=0.0001 | 43.5 / 29.7 (6, 3.53) | 46.2 / 30.5 (6, 3.75) | 47.1 / 31.2 (6, 3.85) |
| Advection | N (P') | dir.+2F inj. FFF|- theta=0.25 tauA=0.0001 | 17.9 / 17.4 (4, 3.58) | 18.6 / 17.9 (4, 3.78) | 19.2 / 17.9 (4, 3.94) |
| Recirc. 1e-4 | lAIR | 1pt lAIR2 -|FC theta=0.25 tauA=0.001 | 103.6 / 77.5 (14, 3.18) | 138.4 / 79.2 (17, 3.45) | 1189.5 / 97.2 (129, 3.86) |
| Recirc. 1e-4 | N (P') | dir.+2F inj. FF|- theta=0.5 tauA=0.0001 | 57.7 / 57.4 (9, 4.46) | 75.2 / 72.9 (11, 4.72) | 98.4 / 86.2 (13, 5.09) |
| Recirc. 1e-3 | lAIR | 1pt lAIR2 -|FFC theta=0.5 tauA=0.0001 | 116.2 / 118.4 (9, 5.22) | 152.8 / 104.4 (10, 6.13) | 181.9 / 124.6 (11, 6.56) |
| Recirc. 1e-3 | N (P') [dir.+2F instead of dir. Z2] | dir.+2F inj. -|FF theta=0.5 tauA=0.0001 | 346.2 / 82.5 (34, 5.56) | 574.3 / 105.1 (47, 6.40) | 969.3 / 118.3 (73, 6.74) |

## Table 3 (reduced search): best total WU to 1e-10 per family, ours / paper

| problem | m | lAIR | lAIR+P' | RS | N (P') |
|---|---|---|---|---|---|
