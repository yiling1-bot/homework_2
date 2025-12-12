"""
Problem 1: homogeneous PDEs in 1D Cartesian coordinates
Compare analytic solutions with finite difference results for:
  1) Steady Laplace equation
  2) Heat equation (parabolic)
  3) Wave equation (hyperbolic)
Boundary: u(0,t)=0, u(1,t)=0 (homogeneous Dirichlet)
Initial: u(x,0)=sin(pi x) (for time-dependent problems)
"""

import os
import numpy as np
import matplotlib.pyplot as plt

output_dir = os.path.dirname(os.path.abspath(__file__))

# ---------- 1. Steady Laplace ----------
# PDE: d^2 u / dx^2 = 0, u(0)=0, u(1)=0
# Analytic: u(x)=0
print("\n" + "=" * 60)
print("Problem 1: Steady Laplace (homogeneous)")
print("=" * 60)

L = 1.0
nx = 101
dx = L / (nx - 1)
x = np.linspace(0.0, L, nx)

u_analytic_steady = np.zeros_like(x)

# Finite difference discretization: A u = b
A = np.zeros((nx, nx))
b = np.zeros(nx)

A[0, 0] = 1.0
A[-1, -1] = 1.0

for i in range(1, nx - 1):
    A[i, i - 1] = 1.0
    A[i, i] = -2.0
    A[i, i + 1] = 1.0

u_num_steady = np.linalg.solve(A, b)
error_steady = np.max(np.abs(u_num_steady - u_analytic_steady))
print(f"Max error (steady): {error_steady:.2e}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(x, u_analytic_steady, "b-", lw=2, label="Analytic")
axes[0].plot(x, u_num_steady, "ro", markevery=5, ms=4, label="FD")
axes[0].set_xlabel("x")
axes[0].set_ylabel("u")
axes[0].set_title("Steady Laplace: analytic vs FD")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].semilogy(x, np.abs(u_num_steady - u_analytic_steady) + 1e-16, "g-")
axes[1].set_xlabel("x")
axes[1].set_ylabel("|error|")
axes[1].set_title("Steady Laplace: absolute error")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "01_steady_laplace_compare.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved: 01_steady_laplace_compare.png")

# ---------- 2. Heat equation ----------
# PDE: u_t = D u_xx, u(0,t)=u(1,t)=0, u(x,0)=sin(pi x)
# Analytic: u(x,t) = exp(-pi^2 D t) sin(pi x)
print("\n" + "=" * 60)
print("Problem 1: Heat equation (homogeneous)")
print("=" * 60)

D = 0.1
nt = 200
dt = 0.0005  # reduce to satisfy Courant<=0.5 at nx=101
t_final = nt * dt

x_heat = np.linspace(0.0, L, nx)
dx_heat = x_heat[1] - x_heat[0]
u_heat = np.sin(np.pi * x_heat)

courant = D * dt / dx_heat**2
print(f"Courant number = {courant:.4f} (<=0.5 for stability)")

for _ in range(nt):
    u_new = u_heat.copy()
    u_new[1:-1] = (
        u_heat[1:-1]
        + D * dt / dx_heat**2 * (u_heat[0:-2] - 2 * u_heat[1:-1] + u_heat[2:])
    )
    u_new[0] = 0.0
    u_new[-1] = 0.0
    u_heat = u_new

u_analytic_heat = np.exp(-np.pi**2 * D * t_final) * np.sin(np.pi * x_heat)
error_heat = np.max(np.abs(u_heat - u_analytic_heat))
print(f"Max error (heat): {error_heat:.2e}, t_final={t_final:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(x_heat, u_analytic_heat, "b-", lw=2, label="Analytic")
axes[0].plot(x_heat, u_heat, "ro", markevery=5, ms=4, label="FD")
axes[0].set_xlabel("x")
axes[0].set_ylabel("u")
axes[0].set_title(f"Heat: analytic vs FD (t={t_final:.4f})")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].semilogy(x_heat, np.abs(u_heat - u_analytic_heat) + 1e-16, "g-")
axes[1].set_xlabel("x")
axes[1].set_ylabel("|error|")
axes[1].set_title("Heat: absolute error")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "02_heat_equation_compare.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved: 02_heat_equation_compare.png")

# ---------- 3. Wave equation ----------
# PDE: u_tt = c^2 u_xx, u(0,t)=u(1,t)=0, u(x,0)=sin(pi x), u_t(x,0)=0
# Analytic: u(x,t)=sin(pi x) cos(pi c t)
print("\n" + "=" * 60)
print("Problem 1: Wave equation (homogeneous)")
print("=" * 60)

c = 1.0
nt_wave = 200
dt_wave = 0.0005
t_final_wave = nt_wave * dt_wave

x_wave = np.linspace(0.0, L, nx)
dx_wave = x_wave[1] - x_wave[0]
r = c * dt_wave / dx_wave
print(f"CFL = {r:.4f} (<=1 for stability)")
if r > 1.0:
    raise RuntimeError("CFL too large; reduce dt_wave or increase nx.")

u_prev = np.sin(np.pi * x_wave)
u_curr = u_prev.copy()
u_next = np.zeros_like(u_curr)

for step in range(nt_wave):
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
    u_next[0] = 0.0
    u_next[-1] = 0.0
    u_prev, u_curr = u_curr, u_next
    u_next = u_prev.copy()

u_num_wave = u_curr
u_analytic_wave = np.sin(np.pi * x_wave) * np.cos(c * np.pi * t_final_wave)
error_wave = np.max(np.abs(u_num_wave - u_analytic_wave))
print(f"Max error (wave): {error_wave:.2e}, t_final={t_final_wave:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(x_wave, u_analytic_wave, "b-", lw=2, label="Analytic")
axes[0].plot(x_wave, u_num_wave, "ro", markevery=5, ms=4, label="FD")
axes[0].set_xlabel("x")
axes[0].set_ylabel("u")
axes[0].set_title(f"Wave: analytic vs FD (t={t_final_wave:.4f})")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].semilogy(x_wave, np.abs(u_num_wave - u_analytic_wave) + 1e-16, "g-")
axes[1].set_xlabel("x")
axes[1].set_ylabel("|error|")
axes[1].set_title("Wave: absolute error")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "03_wave_equation_compare.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved: 03_wave_equation_compare.png")

# ---------- Error summary ----------
print("\n" + "=" * 60)
print("Problem 1: error summary")
print("=" * 60)
print(f"Steady (Laplace) max error   : {error_steady:.2e}")
print(f"Heat equation max error      : {error_heat:.2e}")
print(f"Wave equation max error      : {error_wave:.2e}")
print("=" * 60)
