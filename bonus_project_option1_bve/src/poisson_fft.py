"""
FFT-based Poisson solver for ∇²ψ = ζ on a doubly periodic domain.

The continuous problem:
    ∂²ψ/∂x² + ∂²ψ/∂y² = ζ

With periodic boundary conditions in both x and y, the solution in
Fourier space is:
    ψ̂(kx, ky) = -ζ̂(kx, ky) / (kx² + ky²)   for (kx, ky) ≠ (0, 0)
    ψ̂(0, 0) = 0                             (zero-mean condition)
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq


def solve_poisson_fft(zeta, dx, dy):
    """
    Solve ∇²ψ = ζ on a doubly periodic domain using FFT.

    Parameters
    ----------
    zeta : ndarray (ny, nx)
        Vorticity field.
    dx, dy : float
        Grid spacing in x and y (m).

    Returns
    -------
    psi : ndarray (ny, nx)
        Streamfunction with zero mean.
    """
    ny, nx = zeta.shape

    # Wavenumbers (radians per metre)
    kx = 2.0 * np.pi * fftfreq(nx, d=dx)
    ky = 2.0 * np.pi * fftfreq(ny, d=dy)

    KX, KY = np.meshgrid(kx, ky)
    k_sq = KX ** 2 + KY ** 2

    # Forward FFT
    zeta_hat = fft2(zeta)

    # Solve in spectral space
    psi_hat = np.zeros_like(zeta_hat, dtype=complex)
    mask = k_sq > 0
    psi_hat[mask] = -zeta_hat[mask] / k_sq[mask]
    # psi_hat[0, 0] remains 0 → zero mean

    # Inverse FFT
    psi = np.real(ifft2(psi_hat))

    return psi


def laplacian_fft(f, dx, dy):
    """
    Compute ∇²f using FFT spectral differentiation (consistent with solve_poisson_fft).

    This is the exact inverse of solve_poisson_fft: solve_poisson_fft(laplacian_fft(f)) == f
    (up to a constant, with zero-mean enforced).
    """
    ny, nx = f.shape
    kx = 2.0 * np.pi * fftfreq(nx, d=dx)
    ky = 2.0 * np.pi * fftfreq(ny, d=dy)
    KX, KY = np.meshgrid(kx, ky)
    k_sq = KX ** 2 + KY ** 2

    f_hat = fft2(f)
    lap_hat = -k_sq * f_hat
    return np.real(ifft2(lap_hat))
    """
    Test the FFT Poisson solver with a known analytic solution.

    Uses ψ_exact = sin(2πx/Lx) * sin(2πy/Ly) and checks that
    solving ∇²ψ = ζ recovers ψ.
    """
    Lx = nx * dx
    Ly = ny * dy

    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    X, Y = np.meshgrid(x, y)

    kx0 = 2.0 * np.pi / Lx
    ky0 = 2.0 * np.pi / Ly

    psi_exact = np.sin(kx0 * X) * np.sin(ky0 * Y)
    zeta_exact = -(kx0 ** 2 + ky0 ** 2) * psi_exact

    psi_solved = solve_poisson_fft(zeta_exact, dx, dy)

    error = np.max(np.abs(psi_solved - psi_exact))
    rel_error = error / np.max(np.abs(psi_exact))

    print(f"Poisson solver test: max error = {error:.2e}, relative = {rel_error:.2e}")
    if rel_error < 1e-10:
        print("  PASSED")
        return True
    else:
        print("  WARNING: error larger than expected for double precision")
        return False


if __name__ == "__main__":
    test_poisson_solver()
