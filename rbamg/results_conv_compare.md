m=15, 432 cases

| outcome | cases |
|---|---|
| both | 324 |
| only N | 48 |
| only M | 0 |
| neither | 60 |

When both converge: rho(E_N) <= rho(E_M) in 242 of 324; median rho(E_N)/rho(E_M) = 0.834

## Cases where M converges and N does not, or M is faster by more than 1.5x

| problem | smoother | split | k | rho(T) | rho(K) | rho(E_N) | rho(E_M) |
|---|---|---|---|---|---|---|---|

## Full table

| problem | smoother | split | k | rho(T) | rho(K) | rho(E_N) | rho(E_M) | outcome |
|---|---|---|---|---|---|---|---|---|
| Poisson | Jacobi w=0.67 | greedy d=0.3 | 1 | 0.987 | 0.5 | 0.895 | 0.868 | both |
| Poisson | Jacobi w=0.67 | greedy d=0.3 | 2 | 0.987 | 0.5 | 0.595 | 0.586 | both |
| Poisson | Jacobi w=0.67 | greedy d=0.3 | 4 | 0.987 | 0.5 | 0.106 | 0.105 | both |
| Poisson | Jacobi w=0.67 | greedy d=0.6 | 1 | 0.987 | 0.667 | 0.948 | 0.942 | both |
| Poisson | Jacobi w=0.67 | greedy d=0.6 | 2 | 0.987 | 0.667 | 0.835 | 0.822 | both |
| Poisson | Jacobi w=0.67 | greedy d=0.6 | 4 | 0.987 | 0.667 | 0.449 | 0.437 | both |
| Poisson | Jacobi w=0.67 | greedy d=0.9 | 1 | 0.987 | 0.829 | 0.974 | 0.972 | both |
| Poisson | Jacobi w=0.67 | greedy d=0.9 | 2 | 0.987 | 0.829 | 0.934 | 0.931 | both |
| Poisson | Jacobi w=0.67 | greedy d=0.9 | 4 | 0.987 | 0.829 | 0.81 | 0.805 | both |
| Poisson | Jacobi w=0.67 | red-black (C = black) | 1 | 0.987 | 0.333 | 0.825 | 0.75 | both |
| Poisson | Jacobi w=0.67 | red-black (C = black) | 2 | 0.987 | 0.333 | 0.323 | 0.387 | both |
| Poisson | Jacobi w=0.67 | red-black (C = black) | 4 | 0.987 | 0.333 | 0.0161 | 0.0656 | both |
| Poisson | Jacobi w=0.67 | full coarsening (C = even,even) | 1 | 0.987 | 0.796 | 0.974 | 0.971 | both |
| Poisson | Jacobi w=0.67 | full coarsening (C = even,even) | 2 | 0.987 | 0.796 | 0.931 | 0.928 | both |
| Poisson | Jacobi w=0.67 | full coarsening (C = even,even) | 4 | 0.987 | 0.796 | 0.78 | 0.782 | both |
| Poisson | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 1 | 0.987 | 0.916 | 0.983 | 0.982 | both |
| Poisson | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 2 | 0.987 | 0.916 | 0.963 | 0.962 | both |
| Poisson | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 4 | 0.987 | 0.916 | 0.912 | 0.913 | both |
| Poisson | Jacobi w=1 | greedy d=0.3 | 1 | 0.981 | 0.25 | 0.579 | 0.898 | both |
| Poisson | Jacobi w=1 | greedy d=0.3 | 2 | 0.981 | 0.25 | 0.106 | 0.875 | both |
| Poisson | Jacobi w=1 | greedy d=0.3 | 4 | 0.981 | 0.25 | 0.0041 | 0.827 | both |
| Poisson | Jacobi w=1 | greedy d=0.6 | 1 | 0.981 | 0.5 | 0.875 | 0.845 | both |
| Poisson | Jacobi w=1 | greedy d=0.6 | 2 | 0.981 | 0.5 | 0.567 | 0.57 | both |
| Poisson | Jacobi w=1 | greedy d=0.6 | 4 | 0.981 | 0.5 | 0.101 | 0.599 | both |
| Poisson | Jacobi w=1 | greedy d=0.9 | 1 | 0.981 | 0.743 | 0.953 | 0.947 | both |
| Poisson | Jacobi w=1 | greedy d=0.9 | 2 | 0.981 | 0.743 | 0.87 | 0.861 | both |
| Poisson | Jacobi w=1 | greedy d=0.9 | 4 | 0.981 | 0.743 | 0.601 | 0.658 | both |
| Poisson | Jacobi w=1 | red-black (C = black) | 1 | 0.981 | 0 | 0 | 0.962 | both |
| Poisson | Jacobi w=1 | red-black (C = black) | 2 | 0.981 | 0 | 0 | 0.962 | both |
| Poisson | Jacobi w=1 | red-black (C = black) | 4 | 0.981 | 0 | 0 | 0.925 | both |
| Poisson | Jacobi w=1 | full coarsening (C = even,even) | 1 | 0.981 | 0.694 | 0.951 | 0.98 | both |
| Poisson | Jacobi w=1 | full coarsening (C = even,even) | 2 | 0.981 | 0.694 | 0.868 | 0.962 | both |
| Poisson | Jacobi w=1 | full coarsening (C = even,even) | 4 | 0.981 | 0.694 | 0.541 | 0.925 | both |
| Poisson | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.981 | 0.874 | 0.973 | 0.981 | both |
| Poisson | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.981 | 0.874 | 0.937 | 0.939 | both |
| Poisson | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.981 | 0.874 | 0.841 | 0.851 | both |
| Poisson | GS w=1 | greedy d=0.3 | 1 | 0.962 | 0.297 | 0.814 | 0.798 | both |
| Poisson | GS w=1 | greedy d=0.3 | 2 | 0.962 | 0.297 | 0.39 | 0.454 | both |
| Poisson | GS w=1 | greedy d=0.3 | 4 | 0.962 | 0.297 | 0.0334 | 0.115 | both |
| Poisson | GS w=1 | greedy d=0.6 | 1 | 0.962 | 0.445 | 0.848 | 0.836 | both |
| Poisson | GS w=1 | greedy d=0.6 | 2 | 0.962 | 0.445 | 0.541 | 0.558 | both |
| Poisson | GS w=1 | greedy d=0.6 | 4 | 0.962 | 0.445 | 0.103 | 0.162 | both |
| Poisson | GS w=1 | greedy d=0.9 | 1 | 0.962 | 0.668 | 0.925 | 0.924 | both |
| Poisson | GS w=1 | greedy d=0.9 | 2 | 0.962 | 0.668 | 0.807 | 0.815 | both |
| Poisson | GS w=1 | greedy d=0.9 | 4 | 0.962 | 0.668 | 0.475 | 0.532 | both |
| Poisson | GS w=1 | red-black (C = black) | 1 | 0.962 | 0.247 | 0.775 | 0.742 | both |
| Poisson | GS w=1 | red-black (C = black) | 2 | 0.962 | 0.247 | 0.307 | 0.432 | both |
| Poisson | GS w=1 | red-black (C = black) | 4 | 0.962 | 0.247 | 0.0197 | 0.116 | both |
| Poisson | GS w=1 | full coarsening (C = even,even) | 1 | 0.962 | 0.701 | 0.932 | 0.93 | both |
| Poisson | GS w=1 | full coarsening (C = even,even) | 2 | 0.962 | 0.701 | 0.826 | 0.849 | both |
| Poisson | GS w=1 | full coarsening (C = even,even) | 4 | 0.962 | 0.701 | 0.528 | 0.605 | both |
| Poisson | GS w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.962 | 0.855 | 0.953 | 0.953 | both |
| Poisson | GS w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.962 | 0.855 | 0.899 | 0.905 | both |
| Poisson | GS w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.962 | 0.855 | 0.77 | 0.8 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.3 | 1 | 0.987 | 0.34 | 0.83 | 0.759 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.3 | 2 | 0.987 | 0.34 | 0.337 | 0.394 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.3 | 4 | 0.987 | 0.34 | 0.0177 | 0.0663 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.6 | 1 | 0.987 | 0.67 | 0.96 | 0.952 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.6 | 2 | 0.987 | 0.67 | 0.87 | 0.864 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.6 | 4 | 0.987 | 0.67 | 0.518 | 0.528 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.9 | 1 | 0.987 | 0.67 | 0.96 | 0.952 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.9 | 2 | 0.987 | 0.67 | 0.87 | 0.864 | both |
| aniso 0.01 | Jacobi w=0.67 | greedy d=0.9 | 4 | 0.987 | 0.67 | 0.518 | 0.528 | both |
| aniso 0.01 | Jacobi w=0.67 | red-black (C = black) | 1 | 0.987 | 0.333 | 0.825 | 0.75 | both |
| aniso 0.01 | Jacobi w=0.67 | red-black (C = black) | 2 | 0.987 | 0.333 | 0.323 | 0.387 | both |
| aniso 0.01 | Jacobi w=0.67 | red-black (C = black) | 4 | 0.987 | 0.333 | 0.0161 | 0.0655 | both |
| aniso 0.01 | Jacobi w=0.67 | full coarsening (C = even,even) | 1 | 0.987 | 0.981 | 0.983 | 0.983 | both |
| aniso 0.01 | Jacobi w=0.67 | full coarsening (C = even,even) | 2 | 0.987 | 0.981 | 0.966 | 0.966 | both |
| aniso 0.01 | Jacobi w=0.67 | full coarsening (C = even,even) | 4 | 0.987 | 0.981 | 0.933 | 0.933 | both |
| aniso 0.01 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 1 | 0.987 | 0.984 | 0.985 | 0.985 | both |
| aniso 0.01 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 2 | 0.987 | 0.984 | 0.971 | 0.971 | both |
| aniso 0.01 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 4 | 0.987 | 0.984 | 0.942 | 0.942 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.3 | 1 | 0.981 | 0.00971 | 0.012 | 0.96 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.3 | 2 | 0.981 | 0.00971 | 9.45e-05 | 0.961 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.3 | 4 | 0.981 | 0.00971 | 8.89e-09 | 0.924 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.6 | 1 | 0.981 | 0.505 | 0.904 | 0.98 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.6 | 2 | 0.981 | 0.505 | 0.64 | 0.805 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.6 | 4 | 0.981 | 0.505 | 0.126 | 0.813 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.9 | 1 | 0.981 | 0.505 | 0.904 | 0.98 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.9 | 2 | 0.981 | 0.505 | 0.64 | 0.805 | both |
| aniso 0.01 | Jacobi w=1 | greedy d=0.9 | 4 | 0.981 | 0.505 | 0.126 | 0.813 | both |
| aniso 0.01 | Jacobi w=1 | red-black (C = black) | 1 | 0.981 | 0 | 0 | 0.962 | both |
| aniso 0.01 | Jacobi w=1 | red-black (C = black) | 2 | 0.981 | 0 | 0 | 0.962 | both |
| aniso 0.01 | Jacobi w=1 | red-black (C = black) | 4 | 0.981 | 0 | 0 | 0.925 | both |
| aniso 0.01 | Jacobi w=1 | full coarsening (C = even,even) | 1 | 0.981 | 0.971 | 0.974 | 0.978 | both |
| aniso 0.01 | Jacobi w=1 | full coarsening (C = even,even) | 2 | 0.981 | 0.971 | 0.949 | 0.962 | both |
| aniso 0.01 | Jacobi w=1 | full coarsening (C = even,even) | 4 | 0.981 | 0.971 | 0.9 | 0.925 | both |
| aniso 0.01 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.981 | 0.976 | 0.978 | 0.981 | both |
| aniso 0.01 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.981 | 0.976 | 0.956 | 0.956 | both |
| aniso 0.01 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.981 | 0.976 | 0.913 | 0.913 | both |
| aniso 0.01 | GS w=1 | greedy d=0.3 | 1 | 0.962 | 0.252 | 0.779 | 0.746 | both |
| aniso 0.01 | GS w=1 | greedy d=0.3 | 2 | 0.962 | 0.252 | 0.316 | 0.432 | both |
| aniso 0.01 | GS w=1 | greedy d=0.3 | 4 | 0.962 | 0.252 | 0.0209 | 0.11 | both |
| aniso 0.01 | GS w=1 | greedy d=0.6 | 1 | 0.962 | 0.502 | 0.923 | 0.927 | both |
| aniso 0.01 | GS w=1 | greedy d=0.6 | 2 | 0.962 | 0.502 | 0.76 | 0.789 | both |
| aniso 0.01 | GS w=1 | greedy d=0.6 | 4 | 0.962 | 0.502 | 0.297 | 0.416 | both |
| aniso 0.01 | GS w=1 | greedy d=0.9 | 1 | 0.962 | 0.502 | 0.923 | 0.927 | both |
| aniso 0.01 | GS w=1 | greedy d=0.9 | 2 | 0.962 | 0.502 | 0.76 | 0.789 | both |
| aniso 0.01 | GS w=1 | greedy d=0.9 | 4 | 0.962 | 0.502 | 0.297 | 0.416 | both |
| aniso 0.01 | GS w=1 | red-black (C = black) | 1 | 0.962 | 0.252 | 0.778 | 0.746 | both |
| aniso 0.01 | GS w=1 | red-black (C = black) | 2 | 0.962 | 0.252 | 0.312 | 0.431 | both |
| aniso 0.01 | GS w=1 | red-black (C = black) | 4 | 0.962 | 0.252 | 0.0203 | 0.109 | both |
| aniso 0.01 | GS w=1 | full coarsening (C = even,even) | 1 | 0.962 | 0.943 | 0.95 | 0.95 | both |
| aniso 0.01 | GS w=1 | full coarsening (C = even,even) | 2 | 0.962 | 0.943 | 0.901 | 0.901 | both |
| aniso 0.01 | GS w=1 | full coarsening (C = even,even) | 4 | 0.962 | 0.943 | 0.809 | 0.809 | both |
| aniso 0.01 | GS w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.962 | 0.953 | 0.957 | 0.958 | both |
| aniso 0.01 | GS w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.962 | 0.953 | 0.915 | 0.915 | both |
| aniso 0.01 | GS w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.962 | 0.953 | 0.834 | 0.834 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.3 | 1 | 0.987 | 0.333 | 0.825 | 0.75 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.3 | 2 | 0.987 | 0.333 | 0.322 | 0.387 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.3 | 4 | 0.987 | 0.333 | 0.0161 | 0.0972 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.6 | 1 | 0.987 | 0.658 | 0.949 | 0.942 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.6 | 2 | 0.987 | 0.658 | 0.836 | 0.823 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.6 | 4 | 0.987 | 0.658 | 0.445 | 0.434 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.9 | 1 | 0.987 | 0.817 | 0.973 | 0.972 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.9 | 2 | 0.987 | 0.817 | 0.933 | 0.93 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | greedy d=0.9 | 4 | 0.987 | 0.817 | 0.802 | 0.798 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | red-black (C = black) | 1 | 0.987 | 0.333 | 0.825 | 0.75 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | red-black (C = black) | 2 | 0.987 | 0.333 | 0.322 | 0.387 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | red-black (C = black) | 4 | 0.987 | 0.333 | 0.0161 | 0.0972 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | full coarsening (C = even,even) | 1 | 0.987 | 0.795 | 0.974 | 0.971 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | full coarsening (C = even,even) | 2 | 0.987 | 0.795 | 0.931 | 0.928 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | full coarsening (C = even,even) | 4 | 0.987 | 0.795 | 0.779 | 0.782 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 1 | 0.987 | 0.913 | 0.983 | 0.983 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 2 | 0.987 | 0.913 | 0.963 | 0.962 | both |
| recirc centered Pe 2 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 4 | 0.987 | 0.913 | 0.911 | 0.912 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.3 | 1 | 0.981 | 0 | 0 | 0.962 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.3 | 2 | 0.981 | 0 | 0 | 0.962 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.3 | 4 | 0.981 | 0 | 0 | 0.925 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.6 | 1 | 0.981 | 0.487 | 0.876 | 0.845 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.6 | 2 | 0.981 | 0.487 | 0.563 | 0.612 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.6 | 4 | 0.981 | 0.487 | 0.0962 | 0.628 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.9 | 1 | 0.981 | 0.725 | 0.952 | 0.946 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.9 | 2 | 0.981 | 0.725 | 0.865 | 0.858 | both |
| recirc centered Pe 2 | Jacobi w=1 | greedy d=0.9 | 4 | 0.981 | 0.725 | 0.577 | 0.67 | both |
| recirc centered Pe 2 | Jacobi w=1 | red-black (C = black) | 1 | 0.981 | 0 | 0 | 0.962 | both |
| recirc centered Pe 2 | Jacobi w=1 | red-black (C = black) | 2 | 0.981 | 0 | 0 | 0.962 | both |
| recirc centered Pe 2 | Jacobi w=1 | red-black (C = black) | 4 | 0.981 | 0 | 0 | 0.925 | both |
| recirc centered Pe 2 | Jacobi w=1 | full coarsening (C = even,even) | 1 | 0.981 | 0.692 | 0.951 | 0.98 | both |
| recirc centered Pe 2 | Jacobi w=1 | full coarsening (C = even,even) | 2 | 0.981 | 0.692 | 0.868 | 0.962 | both |
| recirc centered Pe 2 | Jacobi w=1 | full coarsening (C = even,even) | 4 | 0.981 | 0.692 | 0.54 | 0.925 | both |
| recirc centered Pe 2 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.981 | 0.869 | 0.973 | 0.981 | both |
| recirc centered Pe 2 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.981 | 0.869 | 0.937 | 0.938 | both |
| recirc centered Pe 2 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.981 | 0.869 | 0.838 | 0.847 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.3 | 1 | 0.962 | 0.234 | 0.764 | 0.743 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.3 | 2 | 0.962 | 0.234 | 0.302 | 0.449 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.3 | 4 | 0.962 | 0.234 | 0.0191 | 0.192 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.6 | 1 | 0.962 | 0.423 | 0.855 | 0.848 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.6 | 2 | 0.962 | 0.423 | 0.556 | 0.581 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.6 | 4 | 0.962 | 0.423 | 0.105 | 0.178 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.9 | 1 | 0.962 | 0.655 | 0.924 | 0.924 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.9 | 2 | 0.962 | 0.655 | 0.801 | 0.811 | both |
| recirc centered Pe 2 | GS w=1 | greedy d=0.9 | 4 | 0.962 | 0.655 | 0.455 | 0.519 | both |
| recirc centered Pe 2 | GS w=1 | red-black (C = black) | 1 | 0.962 | 0.234 | 0.764 | 0.743 | both |
| recirc centered Pe 2 | GS w=1 | red-black (C = black) | 2 | 0.962 | 0.234 | 0.302 | 0.449 | both |
| recirc centered Pe 2 | GS w=1 | red-black (C = black) | 4 | 0.962 | 0.234 | 0.0191 | 0.192 | both |
| recirc centered Pe 2 | GS w=1 | full coarsening (C = even,even) | 1 | 0.962 | 0.701 | 0.931 | 0.93 | both |
| recirc centered Pe 2 | GS w=1 | full coarsening (C = even,even) | 2 | 0.962 | 0.701 | 0.824 | 0.848 | both |
| recirc centered Pe 2 | GS w=1 | full coarsening (C = even,even) | 4 | 0.962 | 0.701 | 0.522 | 0.599 | both |
| recirc centered Pe 2 | GS w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.962 | 0.852 | 0.953 | 0.953 | both |
| recirc centered Pe 2 | GS w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.962 | 0.852 | 0.897 | 0.904 | both |
| recirc centered Pe 2 | GS w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.962 | 0.852 | 0.765 | 0.795 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.3 | 1 | 0.93 | 0.468 | 0.524 | 0.469 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.3 | 2 | 0.93 | 0.468 | 0.22 | 0.219 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.3 | 4 | 0.93 | 0.468 | 0.0479 | 0.107 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.6 | 1 | 0.93 | 0.678 | 0.801 | 0.783 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.6 | 2 | 0.93 | 0.678 | 0.538 | 0.513 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.6 | 4 | 0.93 | 0.678 | 0.214 | 0.228 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.9 | 1 | 0.93 | 0.802 | 0.878 | 0.873 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.9 | 2 | 0.93 | 0.802 | 0.731 | 0.713 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | greedy d=0.9 | 4 | 0.93 | 0.802 | 0.421 | 0.415 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | red-black (C = black) | 1 | 0.93 | 0.333 | 0.401 | 0.333 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | red-black (C = black) | 2 | 0.93 | 0.333 | 0.111 | 0.2 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | red-black (C = black) | 4 | 0.93 | 0.333 | 0.0123 | 0.0914 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | full coarsening (C = even,even) | 1 | 0.93 | 0.758 | 0.865 | 0.853 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | full coarsening (C = even,even) | 2 | 0.93 | 0.758 | 0.675 | 0.642 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | full coarsening (C = even,even) | 4 | 0.93 | 0.758 | 0.292 | 0.308 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 1 | 0.93 | 0.866 | 0.91 | 0.906 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 2 | 0.93 | 0.866 | 0.808 | 0.799 | both |
| const-flow upwind Pe 2 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 4 | 0.93 | 0.866 | 0.562 | 0.506 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.3 | 1 | 0.895 | 0.201 | 0.203 | 0.886 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.3 | 2 | 0.895 | 0.201 | 0.0412 | 0.772 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.3 | 4 | 0.895 | 0.201 | 0.00167 | 0.621 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.6 | 1 | 0.895 | 0.516 | 0.627 | 0.875 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.6 | 2 | 0.895 | 0.516 | 0.295 | 0.617 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.6 | 4 | 0.895 | 0.516 | 0.0728 | 0.482 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.9 | 1 | 0.895 | 0.702 | 0.797 | 0.814 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.9 | 2 | 0.895 | 0.702 | 0.528 | 0.522 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | greedy d=0.9 | 4 | 0.895 | 0.702 | 0.25 | 0.393 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | red-black (C = black) | 1 | 0.895 | 0 | 0 | 0.965 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | red-black (C = black) | 2 | 0.895 | 0 | 0 | 0.8 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | red-black (C = black) | 4 | 0.895 | 0 | 0 | 0.641 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | full coarsening (C = even,even) | 1 | 0.895 | 0.637 | 0.762 | 0.913 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | full coarsening (C = even,even) | 2 | 0.895 | 0.637 | 0.467 | 0.8 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | full coarsening (C = even,even) | 4 | 0.895 | 0.637 | 0.165 | 0.641 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.895 | 0.798 | 0.856 | 0.901 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.895 | 0.798 | 0.685 | 0.672 | both |
| const-flow upwind Pe 2 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.895 | 0.798 | 0.429 | 0.499 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.3 | 1 | 0.8 | 0.241 | 0.465 | 0.561 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.3 | 2 | 0.8 | 0.241 | 0.169 | 0.37 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.3 | 4 | 0.8 | 0.241 | 0.0187 | 0.206 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.6 | 1 | 0.8 | 0.454 | 0.538 | 0.552 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.6 | 2 | 0.8 | 0.454 | 0.242 | 0.358 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.6 | 4 | 0.8 | 0.454 | 0.0446 | 0.194 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.9 | 1 | 0.8 | 0.57 | 0.649 | 0.646 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.9 | 2 | 0.8 | 0.57 | 0.385 | 0.436 | both |
| const-flow upwind Pe 2 | GS w=1 | greedy d=0.9 | 4 | 0.8 | 0.57 | 0.126 | 0.237 | both |
| const-flow upwind Pe 2 | GS w=1 | red-black (C = black) | 1 | 0.8 | 0.207 | 0.44 | 0.559 | both |
| const-flow upwind Pe 2 | GS w=1 | red-black (C = black) | 2 | 0.8 | 0.207 | 0.148 | 0.364 | both |
| const-flow upwind Pe 2 | GS w=1 | red-black (C = black) | 4 | 0.8 | 0.207 | 0.0139 | 0.21 | both |
| const-flow upwind Pe 2 | GS w=1 | full coarsening (C = even,even) | 1 | 0.8 | 0.585 | 0.669 | 0.679 | both |
| const-flow upwind Pe 2 | GS w=1 | full coarsening (C = even,even) | 2 | 0.8 | 0.585 | 0.414 | 0.5 | both |
| const-flow upwind Pe 2 | GS w=1 | full coarsening (C = even,even) | 4 | 0.8 | 0.585 | 0.15 | 0.294 | both |
| const-flow upwind Pe 2 | GS w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.8 | 0.712 | 0.757 | 0.739 | both |
| const-flow upwind Pe 2 | GS w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.8 | 0.712 | 0.529 | 0.556 | both |
| const-flow upwind Pe 2 | GS w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.8 | 0.712 | 0.266 | 0.33 | both |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.3 | 1 | 3.02 | 0.333 | 0.825 | 13.6 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.3 | 2 | 3.02 | 0.333 | 0.322 | 18.5 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.3 | 4 | 3.02 | 0.333 | 0.0161 | 245 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.6 | 1 | 3.02 | 0.333 | 0.825 | 13.6 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.6 | 2 | 3.02 | 0.333 | 0.322 | 18.5 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.6 | 4 | 3.02 | 0.333 | 0.0161 | 245 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.9 | 1 | 3.02 | 0.652 | 0.908 | 15.3 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.9 | 2 | 3.02 | 0.652 | 0.683 | 24.4 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | greedy d=0.9 | 4 | 3.02 | 0.652 | 0.178 | 327 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | red-black (C = black) | 1 | 3.02 | 0.333 | 0.825 | 13.6 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | red-black (C = black) | 2 | 3.02 | 0.333 | 0.322 | 18.5 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | red-black (C = black) | 4 | 3.02 | 0.333 | 0.0161 | 245 | only N |
| recirc centered Pe 20 | Jacobi w=0.67 | full coarsening (C = even,even) | 1 | 3.02 | 3.03 | 3.32 | 10.5 | neither |
| recirc centered Pe 20 | Jacobi w=0.67 | full coarsening (C = even,even) | 2 | 3.02 | 3.03 | 12.4 | 19.7 | neither |
| recirc centered Pe 20 | Jacobi w=0.67 | full coarsening (C = even,even) | 4 | 3.02 | 3.03 | 202 | 193 | neither |
| recirc centered Pe 20 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 1 | 3.02 | 3.03 | 3.96 | 8.11 | neither |
| recirc centered Pe 20 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 2 | 3.02 | 3.03 | 10 | 14.9 | neither |
| recirc centered Pe 20 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 4 | 3.02 | 3.03 | 104 | 128 | neither |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.3 | 1 | 4.5 | 0 | 0 | 21.1 | only N |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.3 | 2 | 4.5 | 0 | 0 | 20.2 | only N |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.3 | 4 | 4.5 | 0 | 0 | 409 | only N |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.6 | 1 | 4.5 | 0 | 0 | 21.1 | only N |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.6 | 2 | 4.5 | 0 | 0 | 20.2 | only N |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.6 | 4 | 4.5 | 0 | 0 | 409 | only N |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.9 | 1 | 4.5 | 0.84 | 0.773 | 23.7 | only N |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.9 | 2 | 4.5 | 0.84 | 0.638 | 30.3 | only N |
| recirc centered Pe 20 | Jacobi w=1 | greedy d=0.9 | 4 | 4.5 | 0.84 | 0.46 | 766 | only N |
| recirc centered Pe 20 | Jacobi w=1 | red-black (C = black) | 1 | 4.5 | 0 | 0 | 21.1 | only N |
| recirc centered Pe 20 | Jacobi w=1 | red-black (C = black) | 2 | 4.5 | 0 | 0 | 20.2 | only N |
| recirc centered Pe 20 | Jacobi w=1 | red-black (C = black) | 4 | 4.5 | 0 | 0 | 409 | only N |
| recirc centered Pe 20 | Jacobi w=1 | full coarsening (C = even,even) | 1 | 4.5 | 4.52 | 5.39 | 16.3 | neither |
| recirc centered Pe 20 | Jacobi w=1 | full coarsening (C = even,even) | 2 | 4.5 | 4.52 | 20.5 | 31.9 | neither |
| recirc centered Pe 20 | Jacobi w=1 | full coarsening (C = even,even) | 4 | 4.5 | 4.52 | 419 | 409 | neither |
| recirc centered Pe 20 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 1 | 4.5 | 4.52 | 6.25 | 13.4 | neither |
| recirc centered Pe 20 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 2 | 4.5 | 4.52 | 20.3 | 22.3 | neither |
| recirc centered Pe 20 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 4 | 4.5 | 4.52 | 299 | 397 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.3 | 1 | 20.2 | 6.18 | 5.62 | 3.31e+06 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.3 | 2 | 20.2 | 6.18 | 72.8 | 1.5e+09 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.3 | 4 | 20.2 | 6.18 | 2.83e+04 | 3.48e+11 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.6 | 1 | 20.2 | 6.18 | 5.62 | 3.31e+06 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.6 | 2 | 20.2 | 6.18 | 72.8 | 1.5e+09 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.6 | 4 | 20.2 | 6.18 | 2.83e+04 | 3.48e+11 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.9 | 1 | 20.2 | 6.19 | 5.42 | 5.12e+06 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.9 | 2 | 20.2 | 6.19 | 114 | 5.49e+08 | neither |
| recirc centered Pe 20 | GS w=1 | greedy d=0.9 | 4 | 20.2 | 6.19 | 3.49e+04 | 3.04e+13 | neither |
| recirc centered Pe 20 | GS w=1 | red-black (C = black) | 1 | 20.2 | 6.18 | 5.62 | 3.31e+06 | neither |
| recirc centered Pe 20 | GS w=1 | red-black (C = black) | 2 | 20.2 | 6.18 | 72.8 | 1.5e+09 | neither |
| recirc centered Pe 20 | GS w=1 | red-black (C = black) | 4 | 20.2 | 6.18 | 2.83e+04 | 3.48e+11 | neither |
| recirc centered Pe 20 | GS w=1 | full coarsening (C = even,even) | 1 | 20.2 | 20.5 | 5.9 | 6.29e+04 | neither |
| recirc centered Pe 20 | GS w=1 | full coarsening (C = even,even) | 2 | 20.2 | 20.5 | 167 | 1.15e+07 | neither |
| recirc centered Pe 20 | GS w=1 | full coarsening (C = even,even) | 4 | 20.2 | 20.5 | 6.91e+04 | 1.34e+09 | neither |
| recirc centered Pe 20 | GS w=1 | aggressive (C = every 3rd, 3rd) | 1 | 20.2 | 20.4 | 15.1 | 8.19e+03 | neither |
| recirc centered Pe 20 | GS w=1 | aggressive (C = every 3rd, 3rd) | 2 | 20.2 | 20.4 | 179 | 1.71e+06 | neither |
| recirc centered Pe 20 | GS w=1 | aggressive (C = every 3rd, 3rd) | 4 | 20.2 | 20.4 | 3.7e+04 | 3.16e+09 | neither |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.3 | 1 | 0.652 | 0.399 | 0.399 | 0.399 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.3 | 2 | 0.652 | 0.399 | 0.159 | 0.262 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.3 | 4 | 0.652 | 0.399 | 0.0253 | 0.135 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.6 | 1 | 0.652 | 0.462 | 0.462 | 0.462 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.6 | 2 | 0.652 | 0.462 | 0.213 | 0.213 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.6 | 4 | 0.652 | 0.462 | 0.0454 | 0.172 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.9 | 1 | 0.652 | 0.535 | 0.528 | 0.525 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.9 | 2 | 0.652 | 0.535 | 0.276 | 0.282 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | greedy d=0.9 | 4 | 0.652 | 0.535 | 0.0754 | 0.117 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | red-black (C = black) | 1 | 0.652 | 0.333 | 0.333 | 0.333 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | red-black (C = black) | 2 | 0.652 | 0.333 | 0.111 | 0.284 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | red-black (C = black) | 4 | 0.652 | 0.333 | 0.0123 | 0.118 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | full coarsening (C = even,even) | 1 | 0.652 | 0.563 | 0.524 | 0.524 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | full coarsening (C = even,even) | 2 | 0.652 | 0.563 | 0.275 | 0.275 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | full coarsening (C = even,even) | 4 | 0.652 | 0.563 | 0.0911 | 0.202 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 1 | 0.652 | 0.619 | 0.602 | 0.591 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 2 | 0.652 | 0.619 | 0.338 | 0.338 | both |
| const-flow upwind Pe 20 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 4 | 0.652 | 0.619 | 0.115 | 0.171 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.3 | 1 | 0.479 | 0.0981 | 0.0981 | 0.867 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.3 | 2 | 0.479 | 0.0981 | 0.0097 | 0.279 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.3 | 4 | 0.479 | 0.0981 | 9.35e-05 | 0.0716 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.6 | 1 | 0.479 | 0.192 | 0.192 | 0.658 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.6 | 2 | 0.479 | 0.192 | 0.0371 | 0.243 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.6 | 4 | 0.479 | 0.192 | 0.00137 | 0.0979 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.9 | 1 | 0.479 | 0.303 | 0.287 | 0.976 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.9 | 2 | 0.479 | 0.303 | 0.0915 | 0.229 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | greedy d=0.9 | 4 | 0.479 | 0.303 | 0.00838 | 0.0525 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | red-black (C = black) | 1 | 0.479 | 0 | 0 | 0.979 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | red-black (C = black) | 2 | 0.479 | 0 | 0 | 0.229 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | red-black (C = black) | 4 | 0.479 | 0 | 0 | 0.0525 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | full coarsening (C = even,even) | 1 | 0.479 | 0.345 | 0.286 | 0.743 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | full coarsening (C = even,even) | 2 | 0.479 | 0.345 | 0.119 | 0.229 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | full coarsening (C = even,even) | 4 | 0.479 | 0.345 | 0.0198 | 0.175 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.479 | 0.428 | 0.392 | 0.645 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.479 | 0.428 | 0.231 | 0.214 | both |
| const-flow upwind Pe 20 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.479 | 0.428 | 0.0738 | 0.156 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.3 | 1 | 0.229 | 0.0682 | 0.126 | 0.25 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.3 | 2 | 0.229 | 0.0682 | 0.0144 | 0.0692 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.3 | 4 | 0.229 | 0.0682 | 0.000161 | 0.00481 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.6 | 1 | 0.229 | 0.0871 | 0.14 | 0.246 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.6 | 2 | 0.229 | 0.0871 | 0.0182 | 0.0699 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.6 | 4 | 0.229 | 0.0871 | 0.00027 | 0.00503 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.9 | 1 | 0.229 | 0.107 | 0.119 | 0.249 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.9 | 2 | 0.229 | 0.107 | 0.0128 | 0.0673 | both |
| const-flow upwind Pe 20 | GS w=1 | greedy d=0.9 | 4 | 0.229 | 0.107 | 0.000125 | 0.00457 | both |
| const-flow upwind Pe 20 | GS w=1 | red-black (C = black) | 1 | 0.229 | 0.06 | 0.119 | 0.251 | both |
| const-flow upwind Pe 20 | GS w=1 | red-black (C = black) | 2 | 0.229 | 0.06 | 0.0127 | 0.0683 | both |
| const-flow upwind Pe 20 | GS w=1 | red-black (C = black) | 4 | 0.229 | 0.06 | 0.000123 | 0.0047 | both |
| const-flow upwind Pe 20 | GS w=1 | full coarsening (C = even,even) | 1 | 0.229 | 0.168 | 0.19 | 0.25 | both |
| const-flow upwind Pe 20 | GS w=1 | full coarsening (C = even,even) | 2 | 0.229 | 0.168 | 0.0353 | 0.0668 | both |
| const-flow upwind Pe 20 | GS w=1 | full coarsening (C = even,even) | 4 | 0.229 | 0.168 | 0.00112 | 0.00439 | both |
| const-flow upwind Pe 20 | GS w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.229 | 0.204 | 0.216 | 0.248 | both |
| const-flow upwind Pe 20 | GS w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.229 | 0.204 | 0.0467 | 0.0641 | both |
| const-flow upwind Pe 20 | GS w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.229 | 0.204 | 0.00206 | 0.00399 | both |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.3 | 1 | 15.3 | 0.333 | 0.825 | 350 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.3 | 2 | 15.3 | 0.333 | 0.322 | 468 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.3 | 4 | 15.3 | 0.333 | 0.0161 | 1.64e+05 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.6 | 1 | 15.3 | 0.333 | 0.825 | 350 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.6 | 2 | 15.3 | 0.333 | 0.322 | 468 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.6 | 4 | 15.3 | 0.333 | 0.0161 | 1.64e+05 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.9 | 1 | 15.3 | 0.333 | 0.825 | 350 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.9 | 2 | 15.3 | 0.333 | 0.322 | 468 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | greedy d=0.9 | 4 | 15.3 | 0.333 | 0.0161 | 1.64e+05 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | red-black (C = black) | 1 | 15.3 | 0.333 | 0.825 | 350 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | red-black (C = black) | 2 | 15.3 | 0.333 | 0.322 | 468 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | red-black (C = black) | 4 | 15.3 | 0.333 | 0.0161 | 1.64e+05 | only N |
| recirc centered Pe 100 | Jacobi w=0.67 | full coarsening (C = even,even) | 1 | 15.3 | 15.2 | 95.3 | 273 | neither |
| recirc centered Pe 100 | Jacobi w=0.67 | full coarsening (C = even,even) | 2 | 15.3 | 15.2 | 345 | 457 | neither |
| recirc centered Pe 100 | Jacobi w=0.67 | full coarsening (C = even,even) | 4 | 15.3 | 15.2 | 1.58e+05 | 1.6e+05 | neither |
| recirc centered Pe 100 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 1 | 15.3 | 15.2 | 109 | 244 | neither |
| recirc centered Pe 100 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 2 | 15.3 | 15.2 | 269 | 332 | neither |
| recirc centered Pe 100 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 4 | 15.3 | 15.2 | 5.93e+04 | 6.21e+04 | neither |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.3 | 1 | 22.9 | 0 | 0 | 527 | only N |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.3 | 2 | 22.9 | 0 | 0 | 525 | only N |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.3 | 4 | 22.9 | 0 | 0 | 2.75e+05 | only N |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.6 | 1 | 22.9 | 0 | 0 | 527 | only N |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.6 | 2 | 22.9 | 0 | 0 | 525 | only N |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.6 | 4 | 22.9 | 0 | 0 | 2.75e+05 | only N |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.9 | 1 | 22.9 | 0 | 0 | 527 | only N |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.9 | 2 | 22.9 | 0 | 0 | 525 | only N |
| recirc centered Pe 100 | Jacobi w=1 | greedy d=0.9 | 4 | 22.9 | 0 | 0 | 2.75e+05 | only N |
| recirc centered Pe 100 | Jacobi w=1 | red-black (C = black) | 1 | 22.9 | 0 | 0 | 527 | only N |
| recirc centered Pe 100 | Jacobi w=1 | red-black (C = black) | 2 | 22.9 | 0 | 0 | 525 | only N |
| recirc centered Pe 100 | Jacobi w=1 | red-black (C = black) | 4 | 22.9 | 0 | 0 | 2.75e+05 | only N |
| recirc centered Pe 100 | Jacobi w=1 | full coarsening (C = even,even) | 1 | 22.9 | 22.8 | 144 | 410 | neither |
| recirc centered Pe 100 | Jacobi w=1 | full coarsening (C = even,even) | 2 | 22.9 | 22.8 | 519 | 812 | neither |
| recirc centered Pe 100 | Jacobi w=1 | full coarsening (C = even,even) | 4 | 22.9 | 22.8 | 2.69e+05 | 2.74e+05 | neither |
| recirc centered Pe 100 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 1 | 22.9 | 22.8 | 157 | 357 | neither |
| recirc centered Pe 100 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 2 | 22.9 | 22.8 | 513 | 531 | neither |
| recirc centered Pe 100 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 4 | 22.9 | 22.8 | 2.36e+05 | 2.41e+05 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.3 | 1 | 524 | 166 | 3.18e+03 | 2.23e+10 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.3 | 2 | 524 | 166 | 8.26e+05 | 4.9e+11 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.3 | 4 | 524 | 166 | 1.87e+10 | 8.96e+18 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.6 | 1 | 524 | 166 | 3.18e+03 | 2.23e+10 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.6 | 2 | 524 | 166 | 8.26e+05 | 4.9e+11 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.6 | 4 | 524 | 166 | 1.87e+10 | 8.96e+18 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.9 | 1 | 524 | 166 | 3.18e+03 | 2.23e+10 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.9 | 2 | 524 | 166 | 8.26e+05 | 4.9e+11 | neither |
| recirc centered Pe 100 | GS w=1 | greedy d=0.9 | 4 | 524 | 166 | 1.87e+10 | 8.96e+18 | neither |
| recirc centered Pe 100 | GS w=1 | red-black (C = black) | 1 | 524 | 166 | 3.18e+03 | 2.23e+10 | neither |
| recirc centered Pe 100 | GS w=1 | red-black (C = black) | 2 | 524 | 166 | 8.26e+05 | 4.9e+11 | neither |
| recirc centered Pe 100 | GS w=1 | red-black (C = black) | 4 | 524 | 166 | 1.87e+10 | 8.96e+18 | neither |
| recirc centered Pe 100 | GS w=1 | full coarsening (C = even,even) | 1 | 524 | 522 | 9.65e+03 | 4.99e+08 | neither |
| recirc centered Pe 100 | GS w=1 | full coarsening (C = even,even) | 2 | 524 | 522 | 4.87e+06 | 8.79e+11 | neither |
| recirc centered Pe 100 | GS w=1 | full coarsening (C = even,even) | 4 | 524 | 522 | 1.33e+12 | 8.25e+15 | neither |
| recirc centered Pe 100 | GS w=1 | aggressive (C = every 3rd, 3rd) | 1 | 524 | 522 | 4.56e+03 | 3.2e+07 | neither |
| recirc centered Pe 100 | GS w=1 | aggressive (C = every 3rd, 3rd) | 2 | 524 | 522 | 1.47e+06 | 1.54e+10 | neither |
| recirc centered Pe 100 | GS w=1 | aggressive (C = every 3rd, 3rd) | 4 | 524 | 522 | 1.48e+11 | 4.84e+16 | neither |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.3 | 1 | 0.528 | 0.365 | 0.365 | 0.365 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.3 | 2 | 0.528 | 0.365 | 0.133 | 0.295 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.3 | 4 | 0.528 | 0.365 | 0.0177 | 0.123 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.6 | 1 | 0.528 | 0.395 | 0.395 | 0.395 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.6 | 2 | 0.528 | 0.395 | 0.156 | 0.164 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.6 | 4 | 0.528 | 0.395 | 0.0242 | 0.163 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.9 | 1 | 0.528 | 0.432 | 0.425 | 0.424 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.9 | 2 | 0.528 | 0.432 | 0.18 | 0.309 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | greedy d=0.9 | 4 | 0.528 | 0.432 | 0.0325 | 0.104 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | red-black (C = black) | 1 | 0.528 | 0.333 | 0.333 | 0.333 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | red-black (C = black) | 2 | 0.528 | 0.333 | 0.111 | 0.311 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | red-black (C = black) | 4 | 0.528 | 0.333 | 0.0123 | 0.105 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | full coarsening (C = even,even) | 1 | 0.528 | 0.475 | 0.428 | 0.435 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | full coarsening (C = even,even) | 2 | 0.528 | 0.475 | 0.182 | 0.182 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | full coarsening (C = even,even) | 4 | 0.528 | 0.475 | 0.0646 | 0.187 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 1 | 0.528 | 0.494 | 0.478 | 0.475 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 2 | 0.528 | 0.494 | 0.213 | 0.216 | both |
| const-flow upwind Pe 100 | Jacobi w=0.67 | aggressive (C = every 3rd, 3rd) | 4 | 0.528 | 0.494 | 0.0705 | 0.166 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.3 | 1 | 0.236 | 0.0468 | 0.0468 | 0.864 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.3 | 2 | 0.236 | 0.0468 | 0.00219 | 0.105 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.3 | 4 | 0.236 | 0.0468 | 4.81e-06 | 0.0101 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.6 | 1 | 0.236 | 0.0918 | 0.0918 | 0.575 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.6 | 2 | 0.236 | 0.0918 | 0.0113 | 0.2 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.6 | 4 | 0.236 | 0.0918 | 0.000214 | 0.0356 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.9 | 1 | 0.236 | 0.148 | 0.135 | 0.984 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.9 | 2 | 0.236 | 0.148 | 0.0218 | 0.0538 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | greedy d=0.9 | 4 | 0.236 | 0.148 | 0.000477 | 0.00289 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | red-black (C = black) | 1 | 0.236 | 0 | 0 | 0.986 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | red-black (C = black) | 2 | 0.236 | 0 | 0 | 0.0538 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | red-black (C = black) | 4 | 0.236 | 0 | 0 | 0.00289 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | full coarsening (C = even,even) | 1 | 0.236 | 0.177 | 0.14 | 0.665 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | full coarsening (C = even,even) | 2 | 0.236 | 0.177 | 0.0281 | 0.175 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | full coarsening (C = even,even) | 4 | 0.236 | 0.177 | 0.0071 | 0.157 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.236 | 0.212 | 0.187 | 0.592 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.236 | 0.212 | 0.072 | 0.374 | both |
| const-flow upwind Pe 100 | Jacobi w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.236 | 0.212 | 0.0192 | 0.217 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.3 | 1 | 0.0538 | 0.0159 | 0.0248 | 0.0578 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.3 | 2 | 0.0538 | 0.0159 | 0.000568 | 0.00353 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.3 | 4 | 0.0538 | 0.0159 | 2.75e-07 | 1.23e-05 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.6 | 1 | 0.0538 | 0.0203 | 0.0284 | 0.0572 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.6 | 2 | 0.0538 | 0.0203 | 0.000752 | 0.00356 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.6 | 4 | 0.0538 | 0.0203 | 4.93e-07 | 1.28e-05 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.9 | 1 | 0.0538 | 0.0253 | 0.0248 | 0.0573 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.9 | 2 | 0.0538 | 0.0253 | 0.000624 | 0.00345 | both |
| const-flow upwind Pe 100 | GS w=1 | greedy d=0.9 | 4 | 0.0538 | 0.0253 | 4e-07 | 1.18e-05 | both |
| const-flow upwind Pe 100 | GS w=1 | red-black (C = black) | 1 | 0.0538 | 0.014 | 0.0231 | 0.0578 | both |
| const-flow upwind Pe 100 | GS w=1 | red-black (C = black) | 2 | 0.0538 | 0.014 | 0.000491 | 0.00349 | both |
| const-flow upwind Pe 100 | GS w=1 | red-black (C = black) | 4 | 0.0538 | 0.014 | 2.03e-07 | 1.2e-05 | both |
| const-flow upwind Pe 100 | GS w=1 | full coarsening (C = even,even) | 1 | 0.0538 | 0.0395 | 0.0414 | 0.059 | both |
| const-flow upwind Pe 100 | GS w=1 | full coarsening (C = even,even) | 2 | 0.0538 | 0.0395 | 0.00166 | 0.00351 | both |
| const-flow upwind Pe 100 | GS w=1 | full coarsening (C = even,even) | 4 | 0.0538 | 0.0395 | 2.56e-06 | 1.16e-05 | both |
| const-flow upwind Pe 100 | GS w=1 | aggressive (C = every 3rd, 3rd) | 1 | 0.0538 | 0.0479 | 0.0487 | 0.059 | both |
| const-flow upwind Pe 100 | GS w=1 | aggressive (C = every 3rd, 3rd) | 2 | 0.0538 | 0.0479 | 0.00232 | 0.00342 | both |
| const-flow upwind Pe 100 | GS w=1 | aggressive (C = every 3rd, 3rd) | 4 | 0.0538 | 0.0479 | 5.15e-06 | 1.08e-05 | both |
