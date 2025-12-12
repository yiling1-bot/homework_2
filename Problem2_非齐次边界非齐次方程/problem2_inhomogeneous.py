"""
Problem 2: inhomogeneous PDEs in 1D Cartesian coordinates
Compare analytic/benchmark solutions with finite difference results for:
  1) Poisson (steady, inhomogeneous boundary)
  2) Heat with inhomogeneous Dirichlet boundary (Crank-Nicolson)
  3) Wave with boundary forcing (centered differences; benchmark is simplified modal approximation)
"""

import os
import numpy as np
import matplotlib.pyplot as plt

output_dir = os.path.dirname(os.path.abspath(__file__))

# ---------- 1. Poisson ----------
# PDE: u_xx = -2, u(0)=0, u(1)=1
# Analytic: u(x) = -x^2 + 2x
print("\n" + "=" * 60)
print("Problem 2: Poisson (inhomogeneous boundary)")
print("=" * 60)

L = 1.0
nx = 101
dx = L / (nx - 1)
x = np.linspace(0.0, L, nx)

u_analytic_poisson = -x**2 + 2 * x

A = np.zeros((nx, nx))
b = np.zeros(nx)
A[0, 0] = 1.0
A[-1, -1] = 1.0
b[-1] = 1.0

for i in range(1, nx - 1):
    A[i, i - 1] = 1.0
    A[i, i] = -2.0
    A[i, i + 1] = 1.0
    b[i] = -2.0 * dx**2

u_num_poisson = np.linalg.solve(A, b)
error_poisson = np.max(np.abs(u_num_poisson - u_analytic_poisson))
print(f"Max error (Poisson): {error_poisson:.2e}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(x, u_analytic_poisson, "b-", lw=2, label="Analytic")
axes[0].plot(x, u_num_poisson, "ro", markevery=5, ms=4, label="FD")
axes[0].set_xlabel("x")
axes[0].set_ylabel("u")
axes[0].set_title("Poisson: analytic vs FD")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].semilogy(x, np.abs(u_num_poisson - u_analytic_poisson) + 1e-16, "g-")
axes[1].set_xlabel("x")
axes[1].set_ylabel("|error|")
axes[1].set_title("Poisson: absolute error")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "01_poisson_compare.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved: 01_poisson_compare.png")

# ---------- 2. Heat with inhomogeneous Dirichlet ----------
# PDE: u_t = D u_xx, u(0,t)=0, u(1,t)=1, u(x,0)=sin(pi x)
# Steady solution is u_s = x; transient decays. Crank-Nicolson.
print("\n" + "=" * 60)
print("Problem 2: Heat (inhomogeneous Dirichlet, Crank-Nicolson)")
print("=" * 60)

D = 0.1
nt = 2000
dt = 0.001
t_final = nt * dt

x_heat = np.linspace(0.0, L, nx)
dx_heat = x_heat[1] - x_heat[0]
u_heat = np.sin(np.pi * x_heat)

alpha = D * dt / (2 * dx_heat**2)
n_int = nx - 2

diag_main = (1 + 2 * alpha) * np.ones(n_int)
diag_off = -alpha * np.ones(n_int - 1)
A_cn = np.diag(diag_main) + np.diag(diag_off, 1) + np.diag(diag_off, -1)

diag_main_rhs = (1 - 2 * alpha) * np.ones(n_int)
diag_off_rhs = alpha * np.ones(n_int - 1)
B_cn = np.diag(diag_main_rhs) + np.diag(diag_off_rhs, 1) + np.diag(diag_off_rhs, -1)

for _ in range(nt):
    rhs = B_cn @ u_heat[1:-1]
    rhs[0] += alpha * (u_heat[0] + 0.0)   # left boundary = 0
    rhs[-1] += alpha * (u_heat[-1] + 1.0) # right boundary = 1
    u_heat[1:-1] = np.linalg.solve(A_cn, rhs)
    u_heat[0] = 0.0
    u_heat[-1] = 1.0

u_analytic_heat = x_heat  # steady-state reference
error_heat = np.max(np.abs(u_heat - u_analytic_heat))
print(f"Max error (heat, vs steady x): {error_heat:.2e}, t_final={t_final:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(x_heat, u_analytic_heat, "b-", lw=2, label="Steady analytic u=x")
axes[0].plot(x_heat, u_heat, "ro", markevery=5, ms=4, label="CN numeric")
axes[0].set_xlabel("x")
axes[0].set_ylabel("u")
axes[0].set_title(f"Heat (inhomogeneous) t={t_final:.4f}")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].semilogy(x_heat, np.abs(u_heat - u_analytic_heat) + 1e-16, "g-")
axes[1].set_xlabel("x")
axes[1].set_ylabel("|error|")
axes[1].set_title("Heat (inhomogeneous): deviation from steady")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "02_heat_inhomogeneous_compare.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved: 02_heat_inhomogeneous_compare.png")

# ---------- 3. Wave with boundary forcing (manufactured analytic solution) ----------
# Choose analytic: u(x,t) = x + sin(pi x) cos(omega t)
# This satisfies u_tt = c^2 u_xx when omega=pi (source term zero), and boundaries u(0,t)=0, u(1,t)=1 (non-hom).
# Initial conditions: u(x,0)=x+sin(pi x), u_t(x,0)=0.
print("\n" + "=" * 60)
print("Problem 2: Wave (boundary forcing)")
print("=" * 60)

c = 1.0
omega = np.pi
nt_wave = 300
dt_wave = 0.001
t_final_wave = nt_wave * dt_wave

x_wave = np.linspace(0.0, L, nx)
dx_wave = x_wave[1] - x_wave[0]
r = c * dt_wave / dx_wave
print(f"CFL = {r:.4f} (<=1 recommended)")
if r > 1.0:
    raise RuntimeError("CFL too large; reduce dt_wave or increase nx.")

def u_exact_wave(x, t):
    return x + np.sin(np.pi * x) * np.cos(omega * t)

u_prev = u_exact_wave(x_wave, 0.0)
u_curr = u_prev.copy()  # u_t(0)=0
u_next = np.zeros_like(u_curr)

for step in range(nt_wave):
    t_next = (step + 1) * dt_wave
    if step == 0:
        u_next[1:-1] = u_curr[1:-1] + 0.5 * r**2 * (
            u_curr[2:] - 2 * u_curr[1:-1] + u_curr[:-2]
        )
    else:
        u_next[1:-1] = (
            2 * u_curr[1:-1]
            - u_prev[1:-1]
            + r**2 * (u_curr[2:] - 2 * u_curr[1:-1] + u_curr[:-2])
        )
    # Non-hom boundary from analytic solution
    u_next[0] = u_exact_wave(0.0, t_next)
    u_next[-1] = u_exact_wave(1.0, t_next)
    u_prev, u_curr = u_curr, u_next
    u_next = u_prev.copy()

u_num_wave = u_curr
u_ref_wave = u_exact_wave(x_wave, t_final_wave)
error_wave = np.max(np.abs(u_num_wave - u_ref_wave))
print(f"Max error (wave, manufactured exact): {error_wave:.2e}, t_final={t_final_wave:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(x_wave, u_ref_wave, "b-", lw=2, label="Analytic (manufactured)")
axes[0].plot(x_wave, u_num_wave, "ro", markevery=5, ms=4, label="FD")
axes[0].set_xlabel("x")
axes[0].set_ylabel("u")
axes[0].set_title(f"Wave (boundary forcing) t={t_final_wave:.4f}")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].semilogy(x_wave, np.abs(u_num_wave - u_ref_wave) + 1e-16, "g-")
axes[1].set_xlabel("x")
axes[1].set_ylabel("|error|")
axes[1].set_title("Wave (boundary forcing): absolute error")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "03_wave_forcing_compare.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved: 03_wave_forcing_compare.png")

# ---------- Error summary ----------
print("\n" + "=" * 60)
print("Problem 2: error summary")
print("=" * 60)
print(f"Poisson max error                  : {error_poisson:.2e}")
print(f"Heat (inhomogeneous) max error     : {error_heat:.2e}")
print(f"Wave (boundary forcing, vs approx) : {error_wave:.2e}")
print("=" * 60)
