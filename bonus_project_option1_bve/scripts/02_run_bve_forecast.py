#!/usr/bin/env python3
"""
02 — Run the barotropic vorticity equation forecast.

Reads the processed initial condition and integrates the BVE forward
for 24 hours.  Saves forecast output at +12 h and +24 h.
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.grid import BetaPlaneGrid
from src.bve_model import BVEModel, initial_condition
from src.operators import laplacian

# ── Configuration ──────────────────────────────────────────────────────────

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'outputs')

DT = 600.0        # time step (s)
FORECAST_HOURS = 24.0
SAVE_INTERVAL = 12.0  # save every 12 h


def main():
    print("=" * 60)
    print("  Step 02: Run BVE Forecast")
    print("=" * 60)

    # ── 1. Load initial condition ──────────────────────────────────────
    init_path = os.path.join(DATA_DIR, 'initial_20251230_00.npz')
    grid_path = os.path.join(DATA_DIR, 'grid_info.npz')

    if not os.path.exists(init_path):
        print(f"ERROR: Initial condition file not found: {init_path}")
        print("Run 01_preprocess_local_data.py first.")
        sys.exit(1)

    init_data = np.load(init_path)
    grid_data = np.load(grid_path)

    Z0 = init_data['Z']
    print(f"Loaded initial field: shape={Z0.shape}")
    print(f"  Z range: {Z0.min():.1f} – {Z0.max():.1f} m")

    # ── 2. Set up grid ─────────────────────────────────────────────────
    grid = BetaPlaneGrid(
        lat0=float(grid_data['lat0']),
        lon0=float(grid_data['lon0']),
        lat_range=tuple(grid_data['lat_range']),
        lon_range=tuple(grid_data['lon_range']),
        dx_deg=float(grid_data['dx_deg']),
        dy_deg=float(grid_data['dy_deg']),
    )

    # ── 3. Compute initial ψ and ζ ─────────────────────────────────────
    psi0, zeta0 = initial_condition(Z0, grid)

    # Check CFL
    # From geostrophic balance, |u| ~ (g/f) * |∇Z|
    # Compute derivatives on the interior (2:end-2 in both dims)
    dZdx = (Z0[2:-2, 3:-1] - Z0[2:-2, 1:-3]) / (2.0 * grid.dx)
    dZdy = (Z0[3:-1, 2:-2] - Z0[1:-3, 2:-2]) / (2.0 * grid.dy)
    u_geo = -grid.G * dZdy / grid.f0
    v_geo = grid.G * dZdx / grid.f0
    u_max = np.max(np.sqrt(u_geo**2 + v_geo**2))
    cfl = u_max * DT / grid.dx
    print(f"  |u|_max (geostrophic) ≈ {u_max:.1f} m/s")
    print(f"  CFL (advection) = {cfl:.3f}")
    if cfl > 0.8:
        print(f"  WARNING: CFL > 0.8. Consider reducing dt to {int(0.8 * grid.dx / u_max)} s.")

    # ── 4. Run model ───────────────────────────────────────────────────
    model = BVEModel(grid, dt=DT, use_arakawa=True)

    history = model.forecast(psi0, zeta0, FORECAST_HOURS,
                              save_interval_hours=SAVE_INTERVAL)

    # ── 5. Save outputs ────────────────────────────────────────────────
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # The history has entries at t=0, t=12, t=24 h
    for idx, t_hours in enumerate(history['time_hours']):
        tag = f"{t_hours:04.1f}h".replace('.', '_')
        out_path = os.path.join(OUTPUT_DIR, f'forecast_{tag}.npz')
        psi = history['psi'][idx]
        zeta = history['zeta'][idx]
        np.savez(out_path,
                 psi=psi, zeta=zeta,
                 time_hours=t_hours,
                 dt=DT,
                 description=f"BVE forecast at t = {t_hours} h")
        print(f"  Saved forecast: {out_path}")

    print("\n" + "=" * 60)
    print("  Forecast complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
