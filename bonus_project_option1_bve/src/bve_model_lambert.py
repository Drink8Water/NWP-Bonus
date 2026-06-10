"""
Lambert-grid finite-area barotropic vorticity equation model.

The prognostic equation is written in conformal projection coordinates as

    ∂ζ/∂t = -m² J(ψ, ζ + f) + ν m² ∇²ζ - α(ζ - ζ_ref),
    ζ = m² ∇²ψ,

with zero streamfunction on the lateral boundaries.  This is the finite-area
counterpart to the original beta-plane/periodic BVE model.
"""

import numpy as np

from .operators import arakawa_jacobian
from .poisson_dirichlet import solve_poisson_dirichlet, laplacian_dirichlet


class LambertBVEModel:
    """BVE model on a Lambert conformal grid with fixed psi boundaries."""

    def __init__(self, grid, dt=600.0, nu=0.0, sponge_width=0,
                 sponge_tau_hours=6.0, zeta_ref=None):
        self.grid = grid
        self.dt = float(dt)
        self.dx = grid.dx
        self.dy = grid.dy
        self.m2 = grid.m ** 2
        self.f = grid.f
        self.ny = grid.ny
        self.nx = grid.nx
        self.nu = float(nu)
        self.sponge_width = int(sponge_width)
        self.sponge_tau = float(sponge_tau_hours) * 3600.0
        self.zeta_ref = zeta_ref
        self._sponge_mask = None

        if self.nu > 0:
            print(f"  Diffusion: nu = {self.nu:.1e} m²/s")
        if self.sponge_width > 0:
            print(f"  Sponge: width = {self.sponge_width} pts, tau = {sponge_tau_hours:.1f} h")

    def _build_sponge_mask(self):
        if self._sponge_mask is not None:
            return self._sponge_mask

        ny, nx = self.ny, self.nx
        width = self.sponge_width
        alpha_max = 1.0 / self.sponge_tau
        j = np.arange(ny).reshape(ny, 1)
        i = np.arange(nx).reshape(1, nx)
        dist_y = np.minimum(j, ny - 1 - j)
        dist_x = np.minimum(i, nx - 1 - i)
        dist = np.minimum(dist_y, dist_x)
        alpha = np.zeros((ny, nx), dtype=float)
        in_layer = dist < width
        alpha[in_layer] = alpha_max * 0.5 * (1.0 + np.cos(np.pi * dist[in_layer] / width))
        self._sponge_mask = alpha
        return alpha

    def psi_from_zeta(self, zeta):
        """Invert ζ = m²∇²ψ with ψ = 0 on the boundary."""
        return solve_poisson_dirichlet(zeta / self.m2, self.dx, self.dy)

    def tendency(self, psi, zeta):
        q = zeta + self.f
        jac = arakawa_jacobian(psi, q, self.dx, self.dy)
        dzdt = -self.m2 * jac

        if self.nu > 0:
            dzdt += self.nu * self.m2 * laplacian_dirichlet(zeta, self.dx, self.dy)

        if self.sponge_width > 0:
            alpha = self._build_sponge_mask()
            dzdt -= alpha * (zeta - self.zeta_ref)

        dzdt[0, :] = 0.0
        dzdt[-1, :] = 0.0
        dzdt[:, 0] = 0.0
        dzdt[:, -1] = 0.0
        return dzdt

    def step_rk4(self, psi, zeta):
        k1 = self.tendency(psi, zeta)

        zeta2 = zeta + 0.5 * self.dt * k1
        psi2 = self.psi_from_zeta(zeta2)
        k2 = self.tendency(psi2, zeta2)

        zeta3 = zeta + 0.5 * self.dt * k2
        psi3 = self.psi_from_zeta(zeta3)
        k3 = self.tendency(psi3, zeta3)

        zeta4 = zeta + self.dt * k3
        psi4 = self.psi_from_zeta(zeta4)
        k4 = self.tendency(psi4, zeta4)

        zeta_new = zeta + (self.dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        zeta_new[0, :] = self.zeta_ref[0, :] if self.zeta_ref is not None else zeta[0, :]
        zeta_new[-1, :] = self.zeta_ref[-1, :] if self.zeta_ref is not None else zeta[-1, :]
        zeta_new[:, 0] = self.zeta_ref[:, 0] if self.zeta_ref is not None else zeta[:, 0]
        zeta_new[:, -1] = self.zeta_ref[:, -1] if self.zeta_ref is not None else zeta[:, -1]
        psi_new = self.psi_from_zeta(zeta_new)
        return psi_new, zeta_new

    def forecast(self, psi0, zeta0, hours, save_hours=(12.0, 24.0)):
        if self.zeta_ref is None:
            self.zeta_ref = zeta0.copy()

        n_steps = int(round(hours * 3600.0 / self.dt))
        save_steps = {int(round(h * 3600.0 / self.dt)): float(h) for h in save_hours}
        history = {
            'time_hours': [0.0],
            'psi': [psi0.copy()],
            'zeta': [zeta0.copy()],
        }

        psi = psi0.copy()
        zeta = zeta0.copy()
        print(f"Integrating Lambert BVE {hours:g} h with dt = {self.dt:g} s ({n_steps} steps)")

        for step in range(1, n_steps + 1):
            psi, zeta = self.step_rk4(psi, zeta)

            if step in save_steps:
                history['time_hours'].append(save_steps[step])
                history['psi'].append(psi.copy())
                history['zeta'].append(zeta.copy())

            if step % max(1, n_steps // 6) == 0:
                print(f"  Step {step}/{n_steps}  t={step*self.dt/3600:.1f} h"
                      f"  |ζ|_max={np.max(np.abs(zeta)):.2e} s⁻¹")

        print("Lambert integration complete.")
        return history


def initial_condition_from_height(z0, grid):
    """
    Build ψ and ζ from geopotential height using variable-f geostrophic balance.

    g Z' = f ψ, with Z' defined relative to the domain mean.  Boundary ψ is set
    to zero to match the finite-area Poisson inversion.
    """
    z_anom = z0 - np.mean(z0)
    psi = grid.G * z_anom / grid.f
    psi = psi - np.mean(psi[1:-1, 1:-1])
    psi[0, :] = 0.0
    psi[-1, :] = 0.0
    psi[:, 0] = 0.0
    psi[:, -1] = 0.0
    zeta = (grid.m ** 2) * laplacian_dirichlet(psi, grid.dx, grid.dy)

    print(f"Lambert initial condition: mean(Z) = {np.mean(z0):.1f} m")
    print(f"  |ψ0|_max = {np.max(np.abs(psi)):.2e} m²/s")
    print(f"  |ζ0|_max = {np.max(np.abs(zeta)):.2e} s⁻¹")
    return psi, zeta


def height_anomaly_from_psi(psi, grid):
    """Recover balanced height anomaly with variable f: Z' = f ψ / g."""
    z_anom = grid.f * psi / grid.G
    return z_anom - np.mean(z_anom)


def vorticity_from_height(z, grid):
    """Diagnostic vorticity from height using the same variable-f balance."""
    psi, _ = initial_condition_from_height(z, grid)
    return (grid.m ** 2) * laplacian_dirichlet(psi, grid.dx, grid.dy)
