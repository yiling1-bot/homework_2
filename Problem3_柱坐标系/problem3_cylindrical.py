"""
Problem 3: Axisymmetric PDEs in cylindrical coordinates (r, z)
Compare finite-difference solutions with reference fields for:
  1) Steady Laplace (axisymmetric)
  2) Heat equation (axisymmetric)
  3) Wave equation (axisymmetric)
Note: Reference steady field here is a manufactured solution for trend comparison;
for strict validation, replace with a Bessel-based separated solution.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

output_dir = os.path.dirname(os.path.abspath(__file__))

# ---------- 1. Steady Poisson with manufactured solution (axisymmetric) ----------
# Exact (manufactured): u(r,z) = (1 - r^2) * sin(pi z)
# Domain: r in [0,1], z in [0,1]; Dirichlet u=0 at r=1, z=0, z=1; symmetry at r=0.
# Source f = -Laplace(u) so that Laplace(u) + f = 0 is satisfied.
print("\n" + "=" * 60)
print("Problem 3: Axisymmetric Poisson (steady, manufactured solution)")
print("=" * 60)

nr = 51
nz = 51
dr = 1.0 / (nr - 1)
dz = 1.0 / (nz - 1)
r = np.linspace(0.0, 1.0, nr)
z = np.linspace(0.0, 1.0, nz)
R, Z = np.meshgrid(r, z)

def manufactured_u(r, z):
    return (1.0 - r**2) * np.sin(np.pi * z)

def compute_rhs(R, Z):
    # Laplace(u) for u=(1-r^2) sin(pi z) is -4 sin(pi z) - pi^2 (1-r^2) sin(pi z)
    lap = -4.0 * np.sin(np.pi * Z) - (np.pi**2) * (1 - R**2) * np.sin(np.pi * Z)
    return lap  # f = Laplace(u_exact)

def solve_poisson_gs(nr, nz, dr, dz, f, max_iter=12000, tol=1e-6):
    """Gauss-Seidel for axisymmetric Poisson: Lap u = f."""
    u = np.zeros((nz, nr))
    # boundary: u=0 at r=1, z=0, z=1; symmetry at r=0
    for it in range(max_iter):
        u_old = u.copy()
        for i in range(1, nz - 1):
            for j in range(1, nr - 1):
                r_j = j * dr
                denom = 2 / dr**2 + 2 / dz**2
                if r_j > 1e-8:
                    num = (
                        (u[i, j + 1] + u[i, j - 1]) / dr**2
                        + (u[i + 1, j] + u[i - 1, j]) / dz**2
                        + (u[i, j + 1] - u[i, j - 1]) / (2 * r_j * dr)
                        - f[i, j]
                    )
                else:
                    num = (
                        2 * u[i, 1] / dr**2
                        + (u[i + 1, j] + u[i - 1, j]) / dz**2
                        - f[i, j]
                    )
                u[i, j] = num / denom
        # boundaries
        u[:, -1] = 0.0
        u[-1, :] = 0.0
        u[0, :] = 0.0
        u[1:-1, 0] = u[1:-1, 1]  # symmetry at r=0
        res = np.max(np.abs(u - u_old))
        if res < tol:
            print(f"Converged in {it} iterations, residual={res:.2e}")
            break
    return u

f_rhs = compute_rhs(R, Z)
u_cyl_steady = solve_poisson_gs(nr, nz, dr, dz, f_rhs)
u_ref_steady = manufactured_u(R, Z)
error_cyl_steady = np.max(np.abs(u_cyl_steady - u_ref_steady))
print(f"Max error (steady, vs manufactured exact): {error_cyl_steady:.2e}")

fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(121, projection="3d")
ax1.plot_surface(R, Z, u_cyl_steady, cmap="viridis", alpha=0.9)
ax1.set_xlabel("r")
ax1.set_ylabel("z")
ax1.set_zlabel("u")
ax1.set_title("Numeric steady solution")
ax2 = fig.add_subplot(122, projection="3d")
ax2.plot_surface(R, Z, u_ref_steady, cmap="plasma", alpha=0.9)
ax2.set_xlabel("r")
ax2.set_ylabel("z")
ax2.set_zlabel("u")
ax2.set_title("Reference field (manufactured)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "01_cyl_steady_3d_compare.png"), dpi=150, bbox_inches="tight")
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
levels = np.linspace(np.min(u_cyl_steady), np.max(u_cyl_steady), 20)
cs1 = axes[0].contourf(R, Z, u_cyl_steady, levels=levels, cmap="viridis")
axes[0].set_xlabel("r")
axes[0].set_ylabel("z")
axes[0].set_title("Numeric steady (contour)")
plt.colorbar(cs1, ax=axes[0])
cs2 = axes[1].contourf(R, Z, np.abs(u_cyl_steady - u_ref_steady), levels=20, cmap="Reds")
axes[1].set_xlabel("r")
axes[1].set_ylabel("z")
axes[1].set_title("Absolute error vs reference")
plt.colorbar(cs2, ax=axes[1])
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "01_cyl_steady_contour.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved steady plots.")

# ---------- 2. Heat equation (axisymmetric) ----------
print("\n" + "=" * 60)
print("Problem 3: Axisymmetric heat equation")
print("=" * 60)

D = 0.1
dt = 1e-4
nt = 60

nr_h = 41
nz_h = 41
dr_h = 1.0 / (nr_h - 1)
dz_h = 1.0 / (nz_h - 1)
r_h = np.linspace(0.0, 1.0, nr_h)
z_h = np.linspace(0.0, 1.0, nz_h)
R_h, Z_h = np.meshgrid(r_h, z_h)

u_heat = np.sin(np.pi * R_h) * np.sin(np.pi * Z_h)
courant_h = D * dt * (1 / dr_h**2 + 1 / dz_h**2)
print(f"Generalized Courant = {courant_h:.4f} (reduce dt if >>1)")

for _ in range(nt):
    u_new = u_heat.copy()
    for i in range(1, nz_h - 1):
        for j in range(1, nr_h - 1):
            r_j = j * dr_h
            if r_j > 1e-8:
                lap_r = (
                    (u_heat[i, j + 1] - 2 * u_heat[i, j] + u_heat[i, j - 1]) / dr_h**2
                    + (u_heat[i, j + 1] - u_heat[i, j - 1]) / (2 * r_j * dr_h)
                )
            else:
                lap_r = 2 * (u_heat[i, 1] - u_heat[i, 0]) / dr_h**2
            lap_z = (u_heat[i + 1, j] - 2 * u_heat[i, j] + u_heat[i - 1, j]) / dz_h**2
            u_new[i, j] = u_heat[i, j] + D * dt * (lap_r + lap_z)
    u_new[:, -1] = 0.0
    u_new[-1, :] = 0.0
    u_new[:, 0] = u_new[:, 1]
    u_heat = u_new

t_final_heat = dt * nt
u_analytic_heat = np.exp(-2 * np.pi**2 * D * t_final_heat) * np.sin(np.pi * R_h) * np.sin(np.pi * Z_h)
error_heat = np.max(np.abs(u_heat - u_analytic_heat))
print(f"Max error (heat): {error_heat:.2e}, t_final={t_final_heat:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
cs1 = axes[0].contourf(R_h, Z_h, u_heat, levels=20, cmap="viridis")
axes[0].set_xlabel("r")
axes[0].set_ylabel("z")
axes[0].set_title(f"Heat numeric (t={t_final_heat:.4f})")
plt.colorbar(cs1, ax=axes[0])
cs2 = axes[1].contourf(R_h, Z_h, np.abs(u_heat - u_analytic_heat), levels=20, cmap="Reds")
axes[1].set_xlabel("r")
axes[1].set_ylabel("z")
axes[1].set_title("Heat absolute error")
plt.colorbar(cs2, ax=axes[1])
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "02_cyl_heat_contour.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved heat plots.")

# ---------- 3. Wave equation (axisymmetric) ----------
print("\n" + "=" * 60)
print("Problem 3: Axisymmetric wave equation")
print("=" * 60)

c = 1.0
dt_w = 1e-4
nt_w = 120

nr_w = 41
nz_w = 41
dr_w = 1.0 / (nr_w - 1)
dz_w = 1.0 / (nz_w - 1)
r_w = np.linspace(0.0, 1.0, nr_w)
z_w = np.linspace(0.0, 1.0, nz_w)
R_w, Z_w = np.meshgrid(r_w, z_w)

u_prev = np.sin(np.pi * R_w) * np.sin(np.pi * Z_w)
u_curr = u_prev.copy()
u_next = np.zeros_like(u_curr)

cfl = c * dt_w / np.sqrt(dr_w**2 * dz_w**2 / (dr_w**2 + dz_w**2))
print(f"CFL (approx) = {cfl:.4f} (<=1 recommended)")
for step in range(nt_w):
    if step == 0:
        for i in range(1, nz_w - 1):
            for j in range(1, nr_w - 1):
                r_j = j * dr_w
                if r_j > 1e-8:
                    lap = (
                        (u_curr[i, j + 1] - 2 * u_curr[i, j] + u_curr[i, j - 1]) / dr_w**2
                        + (u_curr[i, j + 1] - u_curr[i, j - 1]) / (2 * r_j * dr_w)
                        + (u_curr[i + 1, j] - 2 * u_curr[i, j] + u_curr[i - 1, j]) / dz_w**2
                    )
                else:
                    lap = (
                        2 * (u_curr[i, 1] - u_curr[i, 0]) / dr_w**2
                        + (u_curr[i + 1, j] - 2 * u_curr[i, j] + u_curr[i - 1, j]) / dz_w**2
                    )
                u_next[i, j] = u_curr[i, j] + 0.5 * c**2 * dt_w**2 * lap
    else:
        for i in range(1, nz_w - 1):
            for j in range(1, nr_w - 1):
                r_j = j * dr_w
                if r_j > 1e-8:
                    lap = (
                        (u_curr[i, j + 1] - 2 * u_curr[i, j] + u_curr[i, j - 1]) / dr_w**2
                        + (u_curr[i, j + 1] - u_curr[i, j - 1]) / (2 * r_j * dr_w)
                        + (u_curr[i + 1, j] - 2 * u_curr[i, j] + u_curr[i - 1, j]) / dz_w**2
                    )
                else:
                    lap = (
                        2 * (u_curr[i, 1] - u_curr[i, 0]) / dr_w**2
                        + (u_curr[i + 1, j] - 2 * u_curr[i, j] + u_curr[i - 1, j]) / dz_w**2
                    )
                u_next[i, j] = 2 * u_curr[i, j] - u_prev[i, j] + c**2 * dt_w**2 * lap
    u_next[:, -1] = 0.0
    u_next[-1, :] = 0.0
    u_next[:, 0] = u_next[:, 1]
    u_prev, u_curr = u_curr, u_next
    u_next = u_prev.copy()

t_final_wave = dt_w * nt_w
lambda_cyl = np.pi * np.sqrt(2)
u_analytic_wave = (
    np.sin(np.pi * R_w) * np.sin(np.pi * Z_w) * np.cos(lambda_cyl * c * t_final_wave)
)
error_wave = np.max(np.abs(u_curr - u_analytic_wave))
print(f"Max error (wave): {error_wave:.2e}, t_final={t_final_wave:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
cs1 = axes[0].contourf(R_w, Z_w, u_curr, levels=20, cmap="viridis")
axes[0].set_xlabel("r")
axes[0].set_ylabel("z")
axes[0].set_title(f"Wave numeric (t={t_final_wave:.4f})")
plt.colorbar(cs1, ax=axes[0])
cs2 = axes[1].contourf(R_w, Z_w, np.abs(u_curr - u_analytic_wave), levels=20, cmap="Reds")
axes[1].set_xlabel("r")
axes[1].set_ylabel("z")
axes[1].set_title("Wave absolute error")
plt.colorbar(cs2, ax=axes[1])
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "03_cyl_wave_contour.png"), dpi=150, bbox_inches="tight")
plt.close()
print("Saved wave plots.")

# ---------- Error summary ----------
print("\n" + "=" * 60)
print("Problem 3: error summary")
print("=" * 60)
print(f"Steady Laplace (vs manufactured ref) : {error_cyl_steady:.2e}")
print(f"Heat equation max error              : {error_heat:.2e}")
print(f"Wave equation max error              : {error_wave:.2e}")
print("=" * 60)
