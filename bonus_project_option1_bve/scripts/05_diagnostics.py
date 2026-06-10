#!/usr/bin/env python3
"""
05 — Numerical diagnostics for the Lambert BVE forecast experiment.

Runs four diagnostic checks that complement the main forecast verification:

  1. Poisson solver residual test     — verify inversion accuracy
  2. Energy & enstrophy time series   — confirm bounded integration
  3. Scale-decomposed RMSE / ACC      — large-scale vs small-scale skill
  4. CFL / time-step sensitivity      — dt = 300, 600, 900 s comparison

Also computes the skill score (SS) relative to persistence.
"""

import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.lambert_grid import load_lambert_grid
from src.bve_model_lambert import (
    LambertBVEModel, initial_condition_from_height, height_anomaly_from_psi
)
from src.poisson_dirichlet import solve_poisson_dirichlet, laplacian_dirichlet
from src.verification import rmse, bias, debiased_rmse, anomaly_correlation
from scipy.ndimage import gaussian_filter

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed_lambert')
OUT_DIR  = os.path.join(os.path.dirname(__file__), '..', 'outputs')

GRID_PATH   = os.path.join(DATA_DIR, 'grid_info_lambert.npz')
INIT_PATH   = os.path.join(DATA_DIR, 'initial_20251230_00.npz')
ANAL12_PATH = os.path.join(DATA_DIR, 'analysis_20251230_12.npz')
ANAL24_PATH = os.path.join(DATA_DIR, 'analysis_20251231_00.npz')

VERIF_LAT = (25.0, 55.0)
VERIF_LON = (90.0, 145.0)


# ═══════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════

def make_mask(grid):
    return ((grid.lat2d >= VERIF_LAT[0]) & (grid.lat2d <= VERIF_LAT[1])
            & (grid.lon2d >= VERIF_LON[0]) & (grid.lon2d <= VERIF_LON[1]))


def zeta_from_height(z, grid):
    z_anom = z - np.mean(z)
    psi = grid.G * z_anom / grid.f
    psi = psi - np.mean(psi[1:-1, 1:-1])
    psi[0, :] = psi[-1, :] = psi[:, 0] = psi[:, -1] = 0.0
    return (grid.m ** 2) * laplacian_dirichlet(psi, grid.dx, grid.dy)


def scale_decompose(field, sigma=2.0):
    """Split a 2-D field into large-scale (smoothed) and small-scale (residual)."""
    large = gaussian_filter(field, sigma=sigma)
    small = field - large
    return large, small


# ═══════════════════════════════════════════════════════════════════════════
#  1. Poisson solver residual test
# ═══════════════════════════════════════════════════════════════════════════

def test_poisson_residual(grid, psi0, zeta0):
    """Verify that ∇²ψ = ζ/m² is solved to numerical precision."""
    rhs = zeta0 / (grid.m ** 2)
    psi_solved = solve_poisson_dirichlet(rhs, grid.dx, grid.dy)
    laplacian_psi = laplacian_dirichlet(psi_solved, grid.dx, grid.dy)
    residual = (grid.m ** 2) * laplacian_psi - zeta0
    rel_res = np.linalg.norm(residual) / np.linalg.norm(zeta0)
    print(f"\n  Relative residual ‖m²∇²ψ − ζ‖₂ / ‖ζ‖₂ = {rel_res:.2e}")
    print(f"  Max |residual| = {np.max(np.abs(residual)):.2e} s⁻¹")
    return rel_res


# ═══════════════════════════════════════════════════════════════════════════
#  2. Energy & enstrophy time series
# ═══════════════════════════════════════════════════════════════════════════

def compute_ke_enstrophy(psi, zeta, grid):
    """Approximate kinetic energy and enstrophy on the Lambert grid."""
    dA = grid.dx * grid.dy
    mi = grid.m[1:-1, 1:-1]  # interior map factor
    u = -mi * (psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2.0 * grid.dy)
    v =  mi * (psi[1:-1, 2:] - psi[1:-1, :-2]) / (2.0 * grid.dx)
    ke = 0.5 * np.sum(u**2 + v**2) * dA
    ens = 0.5 * np.sum(zeta[1:-1, 1:-1]**2) * dA
    return float(ke), float(ens)


def run_energy_diagnostics(grid, psi0, zeta0):
    """Integrate CTRL for 24 h, saving KE and enstrophy every hour."""
    print("\n  Integrating CTRL with hourly diagnostics ...")
    model = LambertBVEModel(grid, dt=600.0, nu=0.0, sponge_width=0)
    n_steps = int(round(24.0 * 3600.0 / model.dt))
    save_every = int(round(3600.0 / model.dt))  # every hour

    psi = psi0.copy(); zeta = zeta0.copy()
    times_h = [0.0]
    ke_vals = []; ens_vals = []
    ke0, ens0 = compute_ke_enstrophy(psi, zeta, grid)
    ke_vals.append(ke0); ens_vals.append(ens0)

    for step in range(1, n_steps + 1):
        psi, zeta = model.step_rk4(psi, zeta)
        if step % save_every == 0:
            times_h.append(step * model.dt / 3600.0)
            ke, ens = compute_ke_enstrophy(psi, zeta, grid)
            ke_vals.append(ke); ens_vals.append(ens)

    ke_arr = np.array(ke_vals); ens_arr = np.array(ens_vals)
    ke_drift = (ke_arr[-1] - ke_arr[0]) / abs(ke_arr[0]) * 100
    ens_drift = (ens_arr[-1] - ens_arr[0]) / abs(ens_arr[0]) * 100
    print(f"  KE drift over 24 h: {ke_drift:+.2f} %")
    print(f"  Enstrophy drift over 24 h: {ens_drift:+.2f} %")
    print(f"  KE  min/max ratio: {ke_arr.min()/ke_arr[0]:.4f} – {ke_arr.max()/ke_arr[0]:.4f}")
    print(f"  Ens min/max ratio: {ens_arr.min()/ens_arr[0]:.4f} – {ens_arr.max()/ens_arr[0]:.4f}")
    return {'time_h': times_h, 'KE': ke_arr.tolist(), 'enstrophy': ens_arr.tolist(),
            'KE_drift_pct': float(ke_drift), 'Ens_drift_pct': float(ens_drift)}


# ═══════════════════════════════════════════════════════════════════════════
#  3. Scale-decomposed error
# ═══════════════════════════════════════════════════════════════════════════

def run_scale_diagnostics(grid, psi0, zeta0, Z12, Z24, mask):
    """Score CTRL with large-scale and small-scale decomposition."""
    print("\n  Integrating CTRL for scale-decomposition scoring ...")
    model = LambertBVEModel(grid, dt=600.0, nu=0.0, sponge_width=0)
    hist = model.forecast(psi0.copy(), zeta0.copy(), 24.0, (12.0, 24.0))

    results = {}
    for i, lead in enumerate([12.0, 24.0]):
        psi_f = hist['psi'][i + 1]
        z_anal = Z12 if lead == 12.0 else Z24

        zf_anom = height_anomaly_from_psi(psi_f, grid)
        za_anom = z_anal - np.mean(z_anal)

        for scale_name, sigma in [('full', None), ('large', 2.0), ('small', 2.0)]:
            if sigma is None:
                zf_s, za_s = zf_anom, za_anom
            else:
                zf_l, zf_s = scale_decompose(zf_anom, sigma)
                za_l, za_s = scale_decompose(za_anom, sigma)
                zf_s, za_s = (zf_l, za_l) if scale_name == 'large' else (zf_s, za_s)

            key = f'{scale_name}_{int(lead)}h'
            results[key] = {
                'rmse_m': float(rmse(zf_s, za_s, mask)),
                'acc': float(anomaly_correlation(zf_s, za_s, mask)),
            }

        # Print a clean line
        full = results[f'full_{int(lead)}h']
        large = results[f'large_{int(lead)}h']
        small = results[f'small_{int(lead)}h']
        print(f"  +{int(lead)}h  full: RMSE={full['rmse_m']:.1f} ACC={full['acc']:.3f}"
              f"  |  large-scale: RMSE={large['rmse_m']:.1f} ACC={large['acc']:.3f}"
              f"  |  small-scale: RMSE={small['rmse_m']:.1f} ACC={small['acc']:.3f}")
    return results


# ═══════════════════════════════════════════════════════════════════════════
#  4. Time-step sensitivity
# ═══════════════════════════════════════════════════════════════════════════

def run_dt_sensitivity(grid, psi0, zeta0, Z12, Z24, mask):
    """Run CTRL at dt = 300, 600, 900 s and compare +24 h scores."""
    print("\n  Time-step sensitivity (CTRL):")
    results = {}
    for dt in [300, 600, 900]:
        model = LambertBVEModel(grid, dt=float(dt), nu=0.0, sponge_width=0)
        hist = model.forecast(psi0.copy(), zeta0.copy(), 24.0, (12.0, 24.0))
        for i, lead in enumerate([12.0, 24.0]):
            psi_f = hist['psi'][i + 1]
            z_anal = Z12 if lead == 12.0 else Z24
            zf = height_anomaly_from_psi(psi_f, grid)
            za = z_anal - np.mean(z_anal)
            key = f'dt{int(dt)}_{int(lead)}h'
            results[key] = {
                'rmse_m': float(rmse(zf, za, mask)),
                'acc': float(anomaly_correlation(zf, za, mask)),
            }
        print(f"  dt={int(dt):4d}s  +12h RMSE={results[f'dt{int(dt)}_12h']['rmse_m']:.1f}"
              f"  +24h RMSE={results[f'dt{int(dt)}_24h']['rmse_m']:.1f}"
              f"  +24h ACC={results[f'dt{int(dt)}_24h']['acc']:.3f}")
    return results


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 68)
    print("  Numerical Diagnostics — Lambert BVE")
    print("=" * 68)

    os.makedirs(OUT_DIR, exist_ok=True)
    grid = load_lambert_grid(GRID_PATH)
    init  = np.load(INIT_PATH)
    anal12 = np.load(ANAL12_PATH)
    anal24 = np.load(ANAL24_PATH)
    Z0, Z12, Z24 = init['Z'], anal12['Z'], anal24['Z']
    psi0, zeta0 = initial_condition_from_height(Z0, grid)
    mask = make_mask(grid)

    all_diag = {}

    # ---- 1. Poisson residual ----
    print("\n── 1. Poisson solver residual ──")
    all_diag['poisson_rel_residual'] = float(test_poisson_residual(grid, psi0, zeta0))

    # ---- 2. Energy & enstrophy ----
    print("\n── 2. Energy & enstrophy time series ──")
    all_diag['energy_enstrophy'] = run_energy_diagnostics(grid, psi0, zeta0)

    # ---- 3. Scale decomposition ----
    print("\n── 3. Scale-decomposed verification ──")
    all_diag['scale_decomp'] = run_scale_diagnostics(grid, psi0, zeta0, Z12, Z24, mask)

    # ---- 4. Time-step sensitivity ----
    print("\n── 4. Time-step sensitivity ──")
    all_diag['dt_sensitivity'] = run_dt_sensitivity(grid, psi0, zeta0, Z12, Z24, mask)

    # ---- 5. Skill Score ----
    print("\n── 5. Skill Score vs Persistence ──")
    # Compute PERSIST scores directly (initial anomaly unchanged)
    z0_anom_fcst = height_anomaly_from_psi(psi0, grid)
    z12_anom = Z12 - np.mean(Z12)
    z24_anom = Z24 - np.mean(Z24)
    rmse_p12 = float(rmse(z0_anom_fcst, z12_anom, mask))
    rmse_p24 = float(rmse(z0_anom_fcst, z24_anom, mask))
    # Use CTRL scores from the fresh scale-decomp run (section 3)
    rmse_c12 = all_diag['scale_decomp']['full_12h']['rmse_m']
    rmse_c24 = all_diag['scale_decomp']['full_24h']['rmse_m']
    for lead, rp, rc in [('12h', rmse_p12, rmse_c12), ('24h', rmse_p24, rmse_c24)]:
        ss = 1.0 - rc / rp
        all_diag[f'skill_score_{lead}'] = float(ss)
        print(f"  +{lead}  RMSE_persist={rp:.1f}  RMSE_ctrl={rc:.1f}  SS={ss:+.3f}")

    # ---- Save ----
    out_path = os.path.join(OUT_DIR, 'diagnostics.json')
    with open(out_path, 'w') as f:
        json.dump(all_diag, f, indent=2)
    print(f"\nSaved: {out_path}")

    print("\n" + "=" * 68)
    print("  All diagnostics complete.")
    print("=" * 68)


if __name__ == '__main__':
    main()
