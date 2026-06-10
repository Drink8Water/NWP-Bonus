"""
Fast Poisson solver for zero-Dirichlet boundaries.

Solves

    ∇² psi = rhs,    psi = 0 on all four boundaries

using a discrete sine transform implemented with numpy FFTs.  This removes the
doubly periodic assumption in the original FFT solver without requiring scipy.
"""

import numpy as np


def _dst1(x, axis=-1):
    """Unnormalised DST-I along one axis."""
    x = np.asarray(x)
    n = x.shape[axis]
    shape = list(x.shape)
    shape[axis] = 2 * (n + 1)
    work = np.zeros(shape, dtype=float)

    sl = [slice(None)] * x.ndim
    sl[axis] = slice(1, n + 1)
    work[tuple(sl)] = x

    sl[axis] = slice(n + 2, None)
    work[tuple(sl)] = -np.flip(x, axis=axis)

    transformed = np.fft.fft(work, axis=axis)
    sl[axis] = slice(1, n + 1)
    return -transformed[tuple(sl)].imag


def _idst1(x, axis=-1):
    """Inverse of the unnormalised DST-I along one axis."""
    n = x.shape[axis]
    return _dst1(x, axis=axis) / (2.0 * (n + 1))


def solve_poisson_dirichlet(rhs, dx, dy):
    """
    Solve ∇² psi = rhs with psi = 0 on the boundary.

    Parameters
    ----------
    rhs : ndarray (ny, nx)
        Right-hand side on the full grid. Boundary values are ignored.
    dx, dy : float
        Grid spacing in projected metres.

    Returns
    -------
    psi : ndarray (ny, nx)
        Solution with zero boundary values.
    """
    rhs = np.asarray(rhs, dtype=float)
    ny, nx = rhs.shape
    psi = np.zeros_like(rhs)

    if ny <= 2 or nx <= 2:
        return psi

    interior = rhs[1:-1, 1:-1]
    rhs_hat = _dst1(_dst1(interior, axis=0), axis=1)

    jj = np.arange(1, ny - 1).reshape(-1, 1)
    ii = np.arange(1, nx - 1).reshape(1, -1)
    lam_y = 2.0 * (np.cos(np.pi * jj / (ny - 1)) - 1.0) / (dy * dy)
    lam_x = 2.0 * (np.cos(np.pi * ii / (nx - 1)) - 1.0) / (dx * dx)
    psi_hat = rhs_hat / (lam_x + lam_y)

    psi[1:-1, 1:-1] = _idst1(_idst1(psi_hat, axis=0), axis=1)
    return psi


def laplacian_dirichlet(field, dx, dy):
    """Five-point Laplacian on the interior; boundary values are set to zero."""
    field = np.asarray(field)
    out = np.zeros_like(field)
    out[1:-1, 1:-1] = (
        (field[1:-1, 2:] - 2.0 * field[1:-1, 1:-1] + field[1:-1, :-2]) / (dx * dx)
        + (field[2:, 1:-1] - 2.0 * field[1:-1, 1:-1] + field[:-2, 1:-1]) / (dy * dy)
    )
    return out
