#!/usr/bin/env python3
"""
05 — Sensitivity experiments for the BVE forecast.

Runs a matrix of experiments to diagnose the +24 h forecast degradation:

  PERSIST   : Persistence baseline (anomaly unchanged from t=0)
  CTRL      : BVE only (no diffusion, no sponge)
  DIFF      : BVE + Laplacian vorticity diffusion (ν = 2.5e4 m²/s)
  SPONGE    : BVE + lateral sponge relaxation (width=8 pts, tau=6 h)
  DIFF_SPONGE : BVE + diffusion + sponge combined

All forecasts start from the same initial condition.
Scores are computed on the inner verification domain 25°N–55°N, 90°E–145°E.
"""

import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.grid import BetaPlaneGrid
from src.bve_model import BVEModel, initial_condition
from src.verification import verify, forecast_height_anomaly, analysis_height_anomaly

# ── Paths ──────────────────────────────────────────────────────────────────

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'outputs')

# ── Configuration ──────────────────────────────────────────────────────────

DT = 600.0
FORECAST_HOURS = 24.0
SAVE_HOURS = [12.0, 24.0]

# Verification region
VERIF_LAT_RANGE = (25.0, 55.0)
VERIF_LON_RANGE = (90.0, 145.0)

# Experiment definitions
EXPERIMENTS = [
    {
        'name': 'PERSIST',
        'desc': 'Persistence baseline',
        'nu': None,           # not integrated
        'sponge_width': None,
    },
    {
        'name': 'CTRL',
        'desc': 'BVE only',
        'nu': 0.0,
        'sponge_width': 0,
    },
    {
        'name': 'DIFF',
        'desc': 'BVE + diffusion (ν = 2.5e4 m²/s)',
        'nu': 2.5e4,
        'sponge_width': 0,
    },
    {
        'name': 'SPONGE',
        'desc': 'BVE + sponge (width=8, tau=6 h)',
        'nu': 0.0,
        'sponge_width': 8,
    },
    {
        'name': 'DIFF_SPONGE',
        'desc': 'BVE + diffusion + sponge',
        'nu': 2.5e4,
        'sponge_width': 8,
    },
]


# ── Helpers ────────────────────────────────────────────────────────────────

def make_verification_mask(lat_1d, lon_1d, lat_range, lon_range):
    """Create boolean mask for the verification region."""
    lat_2d, lon_2d = np.meshgrid(lat_1d, lon_1d, indexing='ij')
    mask = ((lat_2d >= lat_range[0]) & (lat_2d <= lat_range[1]) &
            (lon_2d >= lon_range[0]) & (lon_2d <= lon_range[1]))
    return mask


def compute_metrics(psi_fcst, z_anal, grid, verif_mask):
    """Return a flat dict of verification metrics."""
    result = verify(psi_fcst, z_anal, grid, verif_mask)
    s = result['scores']
    return {
        'rmse_m': float(s['height_rmse_m']),
        'bias_m': float(s['height_bias_m']),
        'debiased_rmse_m': float(s['height_debiased_rmse_m']),
        'acc': float(s['height_acc']),
        'vorticity_acc': float(s['vorticity_acc']),
        'vorticity_rmse_s-1': float(s['vorticity_rmse']),
    }, result


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    print("=" * 64)
    print("  Sensitivity Experiments — BVE Forecast")
    print("=" * 64)

    # ── 1. Load data ───────────────────────────────────────────────────
    init_path = os.path.join(DATA_DIR, 'initial_20251230_00.npz')
    grid_path = os.path.join(DATA_DIR, 'grid_info.npz')
    anal_12_path = os.path.join(DATA_DIR, 'analysis_20251230_12.npz')
    anal_24_path = os.path.join(DATA_DIR, 'analysis_20251231_00.npz')

    init_data = np.load(init_path)
    grid_data = np.load(grid_path)
    anal_12 = np.load(anal_12_path)
    anal_24 = np.load(anal_24_path)

    Z0 = init_data['Z']
    Z_anal_12 = anal_12['Z']
    Z_anal_24 = anal_24['Z']

    grid = BetaPlaneGrid(
        lat0=float(grid_data['lat0']),
        lon0=float(grid_data['lon0']),
        lat_range=tuple(grid_data['lat_range']),
        lon_range=tuple(grid_data['lon_range']),
        dx_deg=float(grid_data['dx_deg']),
        dy_deg=float(grid_data['dy_deg']),
    )

    # Initial condition
    print("\n── Computing initial condition ──")
    psi0, zeta0 = initial_condition(Z0, grid)

    # Verification mask
    verif_mask = make_verification_mask(
        grid.lat_1d, grid.lon_1d,
        VERIF_LAT_RANGE, VERIF_LON_RANGE
    )
    print(f"Verification region: {VERIF_LAT_RANGE[0]}–{VERIF_LAT_RANGE[1]}°N, "
          f"{VERIF_LON_RANGE[0]}–{VERIF_LON_RANGE[1]}°E  "
          f"({np.sum(verif_mask)} pts)")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── 2. Run experiments ──────────────────────────────────────────────
    all_scores = {}

    for exp in EXPERIMENTS:
        name = exp['name']
        print(f"\n{'─' * 50}")
        print(f"  Experiment: {name} — {exp['desc']}")
        print(f"{'─' * 50}")

        # ── PERSIST: no integration ──
        if name == 'PERSIST':
            # Forecast height anomaly = initial height anomaly at all times
            Z_init_anom = Z0 - np.mean(Z0)
            # Convert to psi via psi = g * Z_anom / f0 (inverse of forecast_height_anomaly)
            psi_persist = grid.G * Z_init_anom / grid.f0

            for lead_h in SAVE_HOURS:
                z_anal = Z_anal_12 if lead_h == 12.0 else Z_anal_24
                metrics, result = compute_metrics(psi_persist, z_anal, grid, verif_mask)
                all_scores[f'{name}_{int(lead_h)}h'] = metrics

                # Save
                out_path = os.path.join(OUTPUT_DIR, f'{name}_forecast_{int(lead_h)}h.npz')
                np.savez(out_path,
                         psi=psi_persist,
                         zeta=zeta0,  # initial zeta, unchanged
                         time_hours=lead_h, experiment=name)
                print(f"  Saved: {out_path}")

                # Save verification arrays
                verif_path = os.path.join(OUTPUT_DIR, f'{name}_verification_{int(lead_h)}h.npz')
                np.savez(verif_path,
                         zf_anom=result['zf_anom'],
                         za_anom=result['za_anom'],
                         zeta_fcst=result['zeta_fcst'],
                         zeta_anal=result['zeta_anal'])
                print(f"  Saved: {verif_path}")

            continue

        # ── Integrated experiments (CTRL, DIFF, SPONGE, DIFF_SPONGE) ──
        model = BVEModel(
            grid,
            dt=DT,
            use_arakawa=True,
            nu=exp['nu'],
            sponge_width=exp['sponge_width'],
            sponge_tau_hours=6.0,
            zeta_ref=zeta0.copy(),
        )

        history = model.forecast(
            psi0.copy(), zeta0.copy(),
            FORECAST_HOURS,
            save_interval_hours=min(SAVE_HOURS),
        )

        # Find the saved states at +12 h and +24 h
        for lead_h in SAVE_HOURS:
            # Find closest saved time
            times = np.array(history['time_hours'])
            idx = np.argmin(np.abs(times - lead_h))
            t_actual = times[idx]

            psi_fcst = history['psi'][idx]
            zeta_fcst = history['zeta'][idx]

            z_anal = Z_anal_12 if lead_h == 12.0 else Z_anal_24
            metrics, result = compute_metrics(psi_fcst, z_anal, grid, verif_mask)
            all_scores[f'{name}_{int(lead_h)}h'] = metrics

            # Save forecast
            out_path = os.path.join(OUTPUT_DIR, f'{name}_forecast_{int(lead_h)}h.npz')
            np.savez(out_path,
                     psi=psi_fcst, zeta=zeta_fcst,
                     time_hours=float(t_actual), experiment=name)
            print(f"  Saved: {out_path}")

            # Save verification arrays
            verif_path = os.path.join(OUTPUT_DIR, f'{name}_verification_{int(lead_h)}h.npz')
            np.savez(verif_path,
                     zf_anom=result['zf_anom'],
                     za_anom=result['za_anom'],
                     zeta_fcst=result['zeta_fcst'],
                     zeta_anal=result['zeta_anal'])
            print(f"  Saved: {verif_path}")

    # ── 3. Print score table and save ──────────────────────────────────
    print(f"\n\n{'═' * 80}")
    print("  SCORE SUMMARY — All Experiments")
    print(f"{'═' * 80}")
    print(f"{'Experiment':<14} {'Lead':>6} {'RMSE':>8} {'Bias':>8} "
          f"{'DebRMSE':>8} {'ACC':>7} {'ζ Corr':>7}")
    print(f"{'─' * 14} {'─' * 6} {'─' * 8} {'─' * 8} {'─' * 8} {'─' * 7} {'─' * 7}")

    for exp_name in [e['name'] for e in EXPERIMENTS]:
        for lead in [12, 24]:
            key = f'{exp_name}_{lead}h'
            if key in all_scores:
                s = all_scores[key]
                print(f"{exp_name:<14} {f'+{lead}h':>6} "
                      f"{s['rmse_m']:8.1f} {s['bias_m']:8.1f} "
                      f"{s['debiased_rmse_m']:8.1f} {s['acc']:7.4f} "
                      f"{s['vorticity_acc']:7.4f}")

    print(f"{'═' * 80}")

    # ── 4. Save scores ─────────────────────────────────────────────────
    # JSON
    json_path = os.path.join(OUTPUT_DIR, 'scores_experiment_matrix.json')
    with open(json_path, 'w') as f:
        json.dump(all_scores, f, indent=2)
    print(f"\nScores saved: {json_path}")

    # CSV
    csv_path = os.path.join(OUTPUT_DIR, 'scores_experiment_matrix.csv')
    with open(csv_path, 'w') as f:
        cols = ['rmse_m', 'bias_m', 'debiased_rmse_m', 'acc', 'vorticity_acc']
        # Parse key into experiment and lead
        rows = []
        for key, s in all_scores.items():
            parts = key.rsplit('_', 1)
            exp_name = parts[0]
            lead = parts[1]
            rows.append([exp_name, lead] + [s[c] for c in cols])

        f.write('experiment,lead,' + ','.join(cols) + '\n')
        for row in sorted(rows):
            f.write(','.join(str(x) for x in row) + '\n')
    print(f"CSV saved: {csv_path}")

    print("\n" + "=" * 64)
    print("  Sensitivity experiments complete!")
    print("=" * 64)


if __name__ == "__main__":
    main()
