"""
Barotropic Vorticity Equation (BVE) model on a beta-plane.

Governing equation:
    ∂ζ/∂t + J(ψ, ζ) + β ∂ψ/∂x = 0                    (CTRL)
    ∂ζ/∂t + J(ψ, ζ) + β ∂ψ/∂x = ν ∇²ζ                (DIFF)
    ∂ζ/∂t + J(ψ, ζ) + β ∂ψ/∂x = -α(x,y)·(ζ − ζ₀)     (SPONGE)
    ∂ζ/∂t + J(ψ, ζ) + β ∂ψ/∂x = ν ∇²ζ − α(x,y)·(ζ − ζ₀)  (DIFF_SPONGE)

    ∇²ψ = ζ

Time integration: RK4 (4th-order Runge–Kutta)
"""

import numpy as np
from .poisson_fft import solve_poisson_fft, laplacian_fft
from .operators import arakawa_jacobian, ddx


class BVEModel:
    """Barotropic vorticity equation model on a beta-plane.

    Supports optional Laplacian diffusion (nu > 0) and
    boundary sponge relaxation (sponge_width > 0).
    """

    def __init__(self, grid, dt=600.0, use_arakawa=True,
                 nu=0.0, sponge_width=0, sponge_tau_hours=6.0,
                 zeta_ref=None):
        """
        Parameters
        ----------
        grid : BetaPlaneGrid
            The computational grid.
        dt : float
            Time step in seconds.
        use_arakawa : bool
            If True, use the Arakawa Jacobian (conserves energy + enstrophy).
        nu : float
            Laplacian diffusion coefficient (m² s⁻¹). 0 = no diffusion.
        sponge_width : int
            Number of grid points for the sponge layer near boundaries.
            0 = no sponge.
        sponge_tau_hours : float
            Sponge relaxation timescale in hours.
        zeta_ref : ndarray (ny, nx) or None
            Reference vorticity field for sponge relaxation.
            If None, uses the initial zeta (must be set before forecast).
        """
        self.grid = grid
        self.dt = dt
        self.use_arakawa = use_arakawa
        self.dx = grid.dx
        self.dy = grid.dy
        self.beta = grid.beta
        self.f0 = grid.f0
        self.ny = grid.ny
        self.nx = grid.nx

        # Diffusion
        self.nu = nu

        # Sponge
        self.sponge_width = int(sponge_width)
        self.sponge_tau = sponge_tau_hours * 3600.0  # convert to seconds
        self.zeta_ref = zeta_ref
        self._sponge_mask = None  # lazy init

        if self.nu > 0:
            print(f"  Diffusion: nu = {self.nu:.1e} m²/s")
        if self.sponge_width > 0:
            print(f"  Sponge: width = {self.sponge_width} pts, tau = {sponge_tau_hours:.1f} h")

    def _build_sponge_mask(self):
        """Build the smooth sponge damping coefficient α(x, y).

        Uses a cosine ramp: α(d) = α_max * (1 − cos(π·d / W)) / 2
        where d is the distance from the nearest boundary (in grid points),
        and W = sponge_width.
        """
        if self._sponge_mask is not None:
            return self._sponge_mask

        ny, nx = self.ny, self.nx
        W = self.sponge_width
        alpha_max = 1.0 / self.sponge_tau  # s⁻¹

        # Distance from each boundary
        dist_n = np.arange(ny).reshape(ny, 1)            # 0 at top (north)
        dist_s = (ny - 1 - np.arange(ny)).reshape(ny, 1) # 0 at bottom (south)
        dist_w = np.arange(nx).reshape(1, nx)             # 0 at left (west)
        dist_e = (nx - 1 - np.arange(nx)).reshape(1, nx)  # 0 at right (east)

        # Minimum distance to any boundary
        dist = np.minimum(np.minimum(dist_n, dist_s),
                          np.minimum(dist_w, dist_e))

        # Cosine ramp
        alpha = np.zeros((ny, nx))
        interior = dist < W
        alpha[interior] = alpha_max * 0.5 * (1.0 - np.cos(np.pi * dist[interior] / W))

        self._sponge_mask = alpha
        return self._sponge_mask

    def tendency(self, psi, zeta):
        """
        Compute ∂ζ/∂t = −J(ψ, ζ) − β ∂ψ/∂x + ν ∇²ζ − α(x,y)·(ζ − ζ₀).

        Parameters
        ----------
        psi : ndarray (ny, nx)
            Streamfunction.
        zeta : ndarray (ny, nx)
            Relative vorticity.

        Returns
        -------
        dzeta_dt : ndarray (ny, nx)
        """
        if self.use_arakawa:
            jac = arakawa_jacobian(psi, zeta, self.dx, self.dy)
        else:
            from .operators import simple_jacobian
            jac = simple_jacobian(psi, zeta, self.dx, self.dy)

        beta_term = self.beta * ddx(psi, self.dx)

        dzetadt = -jac - beta_term

        # ── Optional: Laplacian diffusion ──
        if self.nu > 0:
            dzetadt += self.nu * laplacian_fft(zeta, self.dx, self.dy)

        # ── Optional: Boundary sponge relaxation ──
        if self.sponge_width > 0:
            alpha = self._build_sponge_mask()
            dzetadt -= alpha * (zeta - self.zeta_ref)

        return dzetadt

    def step_rk4(self, psi, zeta):
        """
        Advance one time step using classical RK4.

        Returns (psi_new, zeta_new).
        """
        # Stage 1
        k1 = self.tendency(psi, zeta)

        # Stage 2
        zeta2 = zeta + 0.5 * self.dt * k1
        psi2 = solve_poisson_fft(zeta2, self.dx, self.dy)
        k2 = self.tendency(psi2, zeta2)

        # Stage 3
        zeta3 = zeta + 0.5 * self.dt * k2
        psi3 = solve_poisson_fft(zeta3, self.dx, self.dy)
        k3 = self.tendency(psi3, zeta3)

        # Stage 4
        zeta4 = zeta + self.dt * k3
        psi4 = solve_poisson_fft(zeta4, self.dx, self.dy)
        k4 = self.tendency(psi4, zeta4)

        zeta_new = zeta + (self.dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        psi_new = solve_poisson_fft(zeta_new, self.dx, self.dy)

        return psi_new, zeta_new

    def forecast(self, psi0, zeta0, hours, save_interval_hours=1):
        """
        Integrate the BVE forward in time.

        Parameters
        ----------
        psi0 : ndarray (ny, nx)
            Initial streamfunction.
        zeta0 : ndarray (ny, nx)
            Initial relative vorticity.
        hours : float
            Total forecast length in hours.
        save_interval_hours : float
            How often to save output (default 1 h).

        Returns
        -------
        history : dict
            Keys: 'time_hours', 'psi', 'zeta' — lists of saved states.
        """
        # If sponge is active and no reference zeta was provided, use zeta0
        if self.sponge_width > 0 and self.zeta_ref is None:
            self.zeta_ref = zeta0.copy()
            print(f"  Sponge reference: using initial zeta")

        n_steps = int(hours * 3600.0 / self.dt)
        save_every = max(1, int(save_interval_hours * 3600.0 / self.dt))

        psi = psi0.copy()
        zeta = zeta0.copy()

        history = {
            'time_hours': [0.0],
            'psi': [psi.copy()],
            'zeta': [zeta.copy()],
        }

        print(f"Integrating {hours} h with dt = {self.dt} s ({n_steps} steps)")
        print(f"  Saving every {save_every} steps ({save_interval_hours} h)")

        for step in range(1, n_steps + 1):
            psi, zeta = self.step_rk4(psi, zeta)

            if step % save_every == 0:
                t_hours = step * self.dt / 3600.0
                history['time_hours'].append(t_hours)
                history['psi'].append(psi.copy())
                history['zeta'].append(zeta.copy())

            if step % max(1, n_steps // 10) == 0:
                zeta_max = np.max(np.abs(zeta))
                print(f"  Step {step}/{n_steps}  (t = {step*self.dt/3600:.1f} h)"
                      f"  |ζ|_max = {zeta_max:.2e} s⁻¹")

        print("Integration complete.")
        return history


def initial_condition(z0, grid):
    """
    Compute initial ψ and ζ from 500 hPa geopotential height.

    Φ0 = g * Z0
    ψ0 = (Φ0 - mean(Φ0)) / f0
    ζ0 = ∇²ψ0  (computed via FFT for consistency with Poisson solver)

    Parameters
    ----------
    z0 : ndarray (ny, nx)
        500 hPa geopotential height (m).
    grid : BetaPlaneGrid

    Returns
    -------
    psi0, zeta0 : ndarray (ny, nx)
    """
    phi0 = grid.G * z0  # geopotential
    psi0 = (phi0 - np.mean(phi0)) / grid.f0
    zeta0 = laplacian_fft(psi0, grid.dx, grid.dy)

    print(f"Initial condition: mean(Z) = {np.mean(z0):.1f} m")
    print(f"  |ψ0|_max = {np.max(np.abs(psi0)):.2e} m²/s")
    print(f"  |ζ0|_max = {np.max(np.abs(zeta0)):.2e} s⁻¹")

    return psi0, zeta0
