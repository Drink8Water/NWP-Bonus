#!/usr/bin/env python3
"""
04 — Generate figures for the BVE forecast experiment.

Produces:
  1. Initial ERA5 500 hPa geopotential height
  2. +12 h forecast/analysis/error 3-panel figure
  3. +24 h forecast/analysis/error 3-panel figure
  4. Initial relative vorticity
  5. Bar chart / table of verification scores
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.grid import BetaPlaneGrid
from src.plotting import (
    plot_geopotential_height, plot_forecast_verification,
    plot_vorticity, plot_scores
)
from src.operators import laplacian

# ── Configuration ──────────────────────────────────────────────────────────

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'outputs')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')


def main():
    print("=" * 60)
    print("  Step 04: Generate Figures")
    print("=" * 60)

    os.makedirs(FIG_DIR, exist_ok=True)

    # ── Load grid info ─────────────────────────────────────────────────
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

    # ── 1. Initial ERA5 500 hPa geopotential height ────────────────────
    print("\n[1/5] Initial geopotential height")
    init_data = np.load(os.path.join(DATA_DIR, 'initial_20251230_00.npz'))
    Z0 = init_data['Z']
    plot_geopotential_height(
        Z0, grid.lon2d, grid.lat2d,
        title='ERA5 500 hPa Geopotential Height — 2025-12-30 00 UTC',
        filepath=os.path.join(FIG_DIR, 'fig01_initial_z500.png'),
        contour_interval=60,
    )

    # ── 2. +12 h verification figure ──────────────────────────────────
    print("\n[2/5] +12 h forecast verification")
    verif_12 = np.load(os.path.join(OUTPUT_DIR, 'verification_+12h.npz'))
    plot_forecast_verification(
        verif_12['zf_anom'], verif_12['za_anom'],
        grid.lon2d, grid.lat2d,
        title_prefix='+12 h',
        filepath=os.path.join(FIG_DIR, 'fig02_forecast_verification_12h.png'),
        contour_interval=30,
        error_interval=15,
    )

    # ── 3. +24 h verification figure ──────────────────────────────────
    print("\n[3/5] +24 h forecast verification")
    verif_24 = np.load(os.path.join(OUTPUT_DIR, 'verification_+24h.npz'))
    plot_forecast_verification(
        verif_24['zf_anom'], verif_24['za_anom'],
        grid.lon2d, grid.lat2d,
        title_prefix='+24 h',
        filepath=os.path.join(FIG_DIR, 'fig03_forecast_verification_24h.png'),
        contour_interval=30,
        error_interval=15,
    )

    # ── 4. Initial relative vorticity ─────────────────────────────────
    print("\n[4/5] Initial relative vorticity")
    from src.bve_model import initial_condition
    psi0, zeta0 = initial_condition(Z0, grid)
    plot_vorticity(
        zeta0, grid.lon2d, grid.lat2d,
        title='Initial Relative Vorticity — 2025-12-30 00 UTC',
        filepath=os.path.join(FIG_DIR, 'fig04_initial_vorticity.png'),
    )

    # ── 5. Score summary figure ───────────────────────────────────────
    print("\n[5/5] Verification score summary")
    import json
    scores_path = os.path.join(OUTPUT_DIR, 'scores.json')
    with open(scores_path) as f:
        all_scores = json.load(f)

    plot_scores(
        all_scores.get('+12h', {}),
        all_scores.get('+24h', {}),
        filepath=os.path.join(FIG_DIR, 'fig05_scores_summary.png'),
    )

    print("\n" + "=" * 60)
    print("  All figures saved to figures/")
    print("=" * 60)


if __name__ == "__main__":
    main()
