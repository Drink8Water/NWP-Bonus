"""
Finite-difference operators for a regular 2-D Cartesian grid.

Uses second-order centred differences throughout.
"""

import numpy as np


def ddx(f, dx):
    """Second-order centred difference in x (periodic in x)."""
    result = np.empty_like(f)
    result[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dx)
    # periodic in x
    result[:, 0] = (f[:, 1] - f[:, -1]) / (2.0 * dx)
    result[:, -1] = (f[:, 0] - f[:, -2]) / (2.0 * dx)
    return result


def ddy(f, dy):
    """Second-order centred difference in y (zero-gradient at y-boundaries)."""
    result = np.empty_like(f)
    result[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2.0 * dy)
    # one-sided at y-boundaries
    result[0, :] = (f[1, :] - f[0, :]) / dy
    result[-1, :] = (f[-1, :] - f[-2, :]) / dy
    return result


def laplacian(f, dx, dy):
    """
    Second-order 5-point Laplacian with periodic-x, one-sided y-boundaries.

    Interior uses standard centred differences.
    x-boundaries use periodic wraparound.
    y-boundaries (first and last row) use one-sided second differences.
    """
    ny, nx = f.shape
    result = np.zeros_like(f)

    # Interior: 2nd-order centred in both directions
    # d²f/dx² for all rows, interior columns
    d2x = (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / (dx ** 2)
    # d²f/dy² for interior rows, all columns
    d2y = (f[2:, :] - 2.0 * f[1:-1, :] + f[:-2, :]) / (dy ** 2)

    # Interior points (1:-1 in both directions)
    result[1:-1, 1:-1] = d2x[1:-1, :] + d2y[:, 1:-1]

    # x-periodic boundaries for interior rows
    for j in range(1, ny - 1):
        # Left boundary (periodic: use rightmost as left neighbour)
        result[j, 0] = ((f[j, 1] - 2.0 * f[j, 0] + f[j, -1]) / (dx ** 2)
                        + (f[j + 1, 0] - 2.0 * f[j, 0] + f[j - 1, 0]) / (dy ** 2))
        # Right boundary (periodic: use leftmost as right neighbour)
        result[j, -1] = ((f[j, 0] - 2.0 * f[j, -1] + f[j, -2]) / (dx ** 2)
                         + (f[j + 1, -1] - 2.0 * f[j, -1] + f[j - 1, -1]) / (dy ** 2))

    # y-boundary rows (one-sided d²f/dy² + periodic d²f/dx²)
    for i in range(nx):
        # Top row (j=0)
        d2x_top = (f[0, (i + 1) % nx] - 2.0 * f[0, i] + f[0, (i - 1) % nx]) / (dx ** 2)
        d2y_top = (f[0, i] - 2.0 * f[1, i] + f[2, i]) / (dy ** 2) if ny > 2 else 0.0
        result[0, i] = d2x_top + d2y_top

        # Bottom row (j=ny-1)
        d2x_bot = (f[-1, (i + 1) % nx] - 2.0 * f[-1, i] + f[-1, (i - 1) % nx]) / (dx ** 2)
        d2y_bot = (f[-1, i] - 2.0 * f[-2, i] + f[-3, i]) / (dy ** 2) if ny > 2 else 0.0
        result[-1, i] = d2x_bot + d2y_bot

    return result


def arakawa_jacobian(psi, zeta, dx, dy):
    """
    Arakawa Jacobian J(psi, zeta) = dpsi/dx * dzeta/dy - dpsi/dy * dzeta/dx.

    Uses the Arakawa (1966) discretisation that conserves both kinetic energy
    and enstrophy when integrated over a closed or periodic domain.

    J = (J++ + J+x + Jx+) / 3
    """
    ny, nx = psi.shape

    J_pp = np.zeros_like(psi)
    J_xp = np.zeros_like(psi)
    J_px = np.zeros_like(psi)

    # J++ : standard centred Jacobian
    J_pp[1:-1, 1:-1] = (
        (psi[1:-1, 2:] - psi[1:-1, :-2]) * (zeta[2:, 1:-1] - zeta[:-2, 1:-1])
        - (psi[2:, 1:-1] - psi[:-2, 1:-1]) * (zeta[1:-1, 2:] - zeta[1:-1, :-2])
    ) / (4.0 * dx * dy)

    # J+x : flux form with d(psi, dzeta/dx)/dy - d(psi, dzeta/dy)/dx
    # terms rearranged
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            # J+x contribution
            term1 = (psi[j, i + 1] * (zeta[j + 1, i + 1] - zeta[j - 1, i + 1])
                     - psi[j, i - 1] * (zeta[j + 1, i - 1] - zeta[j - 1, i - 1])) / (4.0 * dx * dy)
            term2 = (psi[j + 1, i] * (zeta[j + 1, i + 1] - zeta[j + 1, i - 1])
                     - psi[j - 1, i] * (zeta[j - 1, i + 1] - zeta[j - 1, i - 1])) / (4.0 * dx * dy)
            J_xp[j, i] = term1 - term2

            # Jx+ contribution
            term3 = (zeta[j, i + 1] * (psi[j + 1, i + 1] - psi[j - 1, i + 1])
                     - zeta[j, i - 1] * (psi[j + 1, i - 1] - psi[j - 1, i - 1])) / (4.0 * dx * dy)
            term4 = (zeta[j + 1, i] * (psi[j + 1, i + 1] - psi[j + 1, i - 1])
                     - zeta[j - 1, i] * (psi[j - 1, i + 1] - psi[j - 1, i - 1])) / (4.0 * dx * dy)
            J_px[j, i] = term3 - term4

    return (J_pp + J_xp + J_px) / 3.0


def simple_jacobian(psi, zeta, dx, dy):
    """Simple second-order centred Jacobian J(psi, zeta)."""
    psi_x = ddx(psi, dx)
    psi_y = ddy(psi, dy)
    zeta_x = ddx(zeta, dx)
    zeta_y = ddy(zeta, dy)
    return psi_x * zeta_y - psi_y * zeta_x
