#!/usr/bin/env python3
"""
03 — Compute verification scores for the BVE forecast.

Compares forecast at +12 h and +24 h against the verifying ERA5 analysis.
Saves scores as scores.json.
"""

import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.grid import BetaPlaneGrid
from src.verification import verify


# ── Configuration ──────────────────────────────────────────────────────────

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'outputs')

# Verification region (inner domain)
VERIF_LAT_RANGE = (25.0, 55.0)
VERIF_LON_RANGE = (90.0, 145.0)

# Forecast and analysis pairs
PAIRS = [
    {
        'label': '+12h',
        'fcst_file': 'forecast_12_0h.npz',
        'anal_file': 'analysis_20251230_12.npz',
    },
    {
        'label': '+24h',
        'fcst_file': 'forecast_24_0h.npz',
        'anal_file': 'analysis_20251231_00.npz',
    },
]


def make_verification_mask(lat, lon, lat_range, lon_range):
    """Create boolean mask for verification region."""
    lat_2d, lon_2d = np.meshgrid(lat, lon, indexing='ij')
    mask = ((lat_2d >= lat_range[0]) & (lat_2d <= lat_range[1]) &
            (lon_2d >= lon_range[0]) & (lon_2d <= lon_range[1]))
    return mask


def main():
    print("=" * 60)
    print("  Step 03: Compute Verification Scores")
    print("=" * 60)

    # ── 1. Load grid info ──────────────────────────────────────────────
    grid_path = os.path.join(DATA_DIR, 'grid_info.npz')
    grid_data = np.load(grid_path)
    grid = BetaPlaneGrid(
        lat0=float(grid_data['lat0']),
        lon0=float(grid_data['lon0']),
        lat_range=tuple(grid_data['lat_range']),
        lon_range=tuple(grid_data['lon_range']),
        dx_deg=float(grid_data['dx_deg']),
        dy_deg=float(grid_data['dy_deg']),
    )

    # Create verification mask on the model grid
    verif_mask = make_verification_mask(
        grid.lat_1d, grid.lon_1d,
        VERIF_LAT_RANGE, VERIF_LON_RANGE
    )
    n_verif = np.sum(verif_mask)
    print(f"Verification region: {VERIF_LAT_RANGE[0]}–{VERIF_LAT_RANGE[1]}°N, "
          f"{VERIF_LON_RANGE[0]}–{VERIF_LON_RANGE[1]}°E")
    print(f"  {n_verif} grid points in verification region")

    all_scores = {}

    for pair in PAIRS:
        label = pair['label']
        print(f"\n{'─'*40}")
        print(f"  {label} verification")

        # Load forecast
        fcst_path = os.path.join(OUTPUT_DIR, pair['fcst_file'])
        anal_path = os.path.join(DATA_DIR, pair['anal_file'])

        if not os.path.exists(fcst_path):
            print(f"  ERROR: Forecast file not found: {fcst_path}")
            continue
        if not os.path.exists(anal_path):
            print(f"  ERROR: Analysis file not found: {anal_path}")
            continue

        fcst_data = np.load(fcst_path)
        anal_data = np.load(anal_path)

        psi_fcst = fcst_data['psi']
        z_anal = anal_data['Z']

        # Verify
        result = verify(psi_fcst, z_anal, grid, verif_mask)
        scores = result['scores']

        print(f"  Height anomaly RMSE:    {scores['height_rmse_m']:.2f} m")
        print(f"  Height full RMSE:       {scores['height_full_rmse_m']:.2f} m")
        print(f"  Height anomaly corr:    {scores['height_acc']:.4f}")
        print(f"  Height full corr:       {scores['height_full_acc']:.4f}")
        print(f"  Vorticity corr:         {scores['vorticity_acc']:.4f}")
        print(f"  Height bias:            {scores['height_bias_m']:.2f} m")

        all_scores[label] = {
            'height_rmse_m': float(scores['height_rmse_m']),
            'height_full_rmse_m': float(scores['height_full_rmse_m']),
            'height_acc': float(scores['height_acc']),
            'height_full_acc': float(scores['height_full_acc']),
            'vorticity_rmse_s-1': float(scores['vorticity_rmse']),
            'vorticity_acc': float(scores['vorticity_acc']),
            'height_bias_m': float(scores['height_bias_m']),
            'verification_region': {
                'lat_range': VERIF_LAT_RANGE,
                'lon_range': VERIF_LON_RANGE,
                'n_points': int(n_verif),
            },
        }

        # Also save the result arrays for plotting
        result_path = os.path.join(OUTPUT_DIR, f'verification_{label}.npz')
        np.savez(result_path,
                 zf_anom=result['zf_anom'],
                 za_anom=result['za_anom'],
                 zf_full=result['zf_full'],
                 zeta_fcst=result['zeta_fcst'],
                 zeta_anal=result['zeta_anal'],
                 psi_anal=result['psi_anal'])
        print(f"  Saved verification arrays: {result_path}")

    # ── 3. Save scores as JSON ────────────────────────────────────────
    scores_path = os.path.join(OUTPUT_DIR, 'scores.json')
    with open(scores_path, 'w') as f:
        json.dump(all_scores, f, indent=2)
    print(f"\n  Scores saved: {scores_path}")

    print("\n" + "=" * 60)
    print("  Score computation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
