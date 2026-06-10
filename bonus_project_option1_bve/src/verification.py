"""
Verification module for the BVE forecast.

Computes RMSE, spatial pattern correlation (anomaly correlation),
and vorticity correlation between forecast and analysis fields.
"""

import numpy as np
from .poisson_fft import solve_poisson_fft, laplacian_fft


def anomaly_correlation(fcst, anal, mask=None):
    """
    Spatial anomaly pattern correlation.

    ACC = Σ[(F_f - μ_f)(F_a - μ_a)] / sqrt[Σ(F_f - μ_f)² Σ(F_a - μ_a)²]

    Parameters
    ----------
    fcst, anal : ndarray
        Forecast and analysis fields.
    mask : ndarray or None
        Boolean mask for the verification region.

    Returns
    -------
    corr : float
    """
    f = fcst.copy()
    a = anal.copy()

    if mask is not None:
        f = f[mask]
        a = a[mask]

    f_anom = f - np.mean(f)
    a_anom = a - np.mean(a)

    num = np.sum(f_anom * a_anom)
    den = np.sqrt(np.sum(f_anom ** 2) * np.sum(a_anom ** 2))

    if den < 1e-30:
        return 0.0

    return num / den


def rmse(fcst, anal, mask=None):
    """
    Root-mean-square error.

    Parameters
    ----------
    fcst, anal : ndarray
    mask : ndarray or None

    Returns
    -------
    rmse_val : float
    """
    f = fcst.copy()
    a = anal.copy()

    if mask is not None:
        f = f[mask]
        a = a[mask]

    return np.sqrt(np.mean((f - a) ** 2))


def bias(fcst, anal, mask=None):
    """Mean bias (forecast - analysis)."""
    f = fcst.copy()
    a = anal.copy()
    if mask is not None:
        f = f[mask]
        a = a[mask]
    return np.mean(f - a)


def debiased_rmse(fcst, anal, mask=None):
    """
    Root-mean-square error after removing the domain-mean bias.

    Debiased RMSE = sqrt(mean(((fcst - bias) - anal)^2))

    This separates the domain-mean drift from the spatial-pattern error.
    """
    f = fcst.copy()
    a = anal.copy()
    if mask is not None:
        f = f[mask]
        a = a[mask]
    b = np.mean(f - a)
    return np.sqrt(np.mean(((f - b) - a) ** 2))


def compute_psi_from_z(z, grid):
    """
    Compute streamfunction from geopotential height.

    Φ = g * Z
    ψ = (Φ - mean(Φ)) / f0

    Parameters
    ----------
    z : ndarray (ny, nx)
        Geopotential height (m).
    grid : BetaPlaneGrid

    Returns
    -------
    psi : ndarray (ny, nx)
    """
    phi = grid.G * z
    return (phi - np.mean(phi)) / grid.f0


def compute_zeta_from_z(z, grid):
    """
    Compute relative vorticity from geopotential height via streamfunction.

    Uses FFT-based Laplacian for consistency with the Poisson solver.

    Parameters
    ----------
    z : ndarray (ny, nx)
        Geopotential height (m).
    grid : BetaPlaneGrid

    Returns
    -------
    zeta : ndarray (ny, nx)
    """
    psi = compute_psi_from_z(z, grid)
    return laplacian_fft(psi, grid.dx, grid.dy)


def forecast_height_anomaly(psi_fcst, grid):
    """
    Convert forecast streamfunction anomaly to geopotential height anomaly.

    Z'_fcst = f0 * ψ_fcst / g

    Parameters
    ----------
    psi_fcst : ndarray (ny, nx)
    grid : BetaPlaneGrid

    Returns
    -------
    Z_anom : ndarray (ny, nx)
    """
    return grid.f0 * psi_fcst / grid.G


def analysis_height_anomaly(z_anal):
    """
    Compute analysis height anomaly by removing the domain mean.

    Parameters
    ----------
    z_anal : ndarray (ny, nx)

    Returns
    -------
    Z_anom : ndarray (ny, nx)
    """
    return z_anal - np.mean(z_anal)


def verify(psi_fcst, z_anal, grid, verification_mask=None):
    """
    Comprehensive verification of a BVE forecast.

    Parameters
    ----------
    psi_fcst : ndarray (ny, nx)
        Forecast streamfunction.
    z_anal : ndarray (ny, nx)
        Verifying analysis geopotential height (m).
    grid : BetaPlaneGrid
    verification_mask : ndarray or None

    Returns
    -------
    scores : dict
    """
    # Forecast height anomaly
    zf_anom = forecast_height_anomaly(psi_fcst, grid)

    # Analysis height anomaly
    za_anom = analysis_height_anomaly(z_anal)

    # Analysis streamfunction and vorticity
    psi_anal = compute_psi_from_z(z_anal, grid)
    zeta_fcst = laplacian_fft(psi_fcst, grid.dx, grid.dy)
    zeta_anal = compute_zeta_from_z(z_anal, grid)

    scores = {}

    # Height anomaly metrics
    scores['height_rmse_m'] = rmse(zf_anom, za_anom, verification_mask)
    scores['height_bias_m'] = bias(zf_anom, za_anom, verification_mask)
    scores['height_debiased_rmse_m'] = debiased_rmse(zf_anom, za_anom, verification_mask)
    scores['height_acc'] = anomaly_correlation(zf_anom, za_anom, verification_mask)

    # Vorticity metrics
    scores['vorticity_rmse'] = rmse(zeta_fcst, zeta_anal, verification_mask)
    scores['vorticity_acc'] = anomaly_correlation(zeta_fcst, zeta_anal, verification_mask)

    return {
        'scores': scores,
        'zf_anom': zf_anom,
        'za_anom': za_anom,
        'zeta_fcst': zeta_fcst,
        'zeta_anal': zeta_anal,
        'psi_anal': psi_anal,
    }
