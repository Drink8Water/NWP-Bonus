#!/usr/bin/env python3
"""
07 — Regrid processed ERA5 fields to a Lambert conformal model grid.

This uses the already processed lat-lon files from data/processed/ as source
data, then writes data/processed_lambert/ for the finite-area Lambert BVE.
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.lambert_grid import LambertGrid
from src.interpolation import interp_regular_latlon


SRC_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed_lambert')

FILES = [
    'initial_20251230_00.npz',
    'analysis_20251230_12.npz',
    'analysis_20251231_00.npz',
]


def main():
    print("=" * 64)
    print("  Step 07: Prepare Lambert Model Data")
    print("=" * 64)

    os.makedirs(OUT_DIR, exist_ok=True)

    grid = LambertGrid(
        center_lon=115.0,
        center_lat=40.0,
        central_lon=115.0,
        central_lat=35.0,
        standard_parallels=(25.0, 45.0),
        width_m=8.0e6,
        height_m=5.0e6,
        dx=150_000.0,
    )

    grid_path = os.path.join(OUT_DIR, 'grid_info_lambert.npz')
    np.savez(grid_path, **grid.to_npz_kwargs())
    print(f"  Saved grid: {grid_path}")

    for name in FILES:
        src_path = os.path.join(SRC_DIR, name)
        src = np.load(src_path)
        lon2d, lat2d = np.meshgrid(src['lon'], src['lat'])

        Z = interp_regular_latlon(src['Z'], lon2d, lat2d, grid.lon2d, grid.lat2d)
        U = interp_regular_latlon(src['U'], lon2d, lat2d, grid.lon2d, grid.lat2d)
        V = interp_regular_latlon(src['V'], lon2d, lat2d, grid.lon2d, grid.lat2d)

        out_path = os.path.join(OUT_DIR, name)
        np.savez(out_path,
                 Z=Z, U=U, V=V,
                 lon2d=grid.lon2d,
                 lat2d=grid.lat2d,
                 X=grid.X,
                 Y=grid.Y,
                 m=grid.m,
                 f=grid.f,
                 time=src['time'],
                 pressure_level=src['pressure_level'],
                 description=f"Lambert-grid ERA5 500 hPa — {src['time']}")
        print(f"  {name}: Z mean={np.nanmean(Z):.1f} m,"
              f" range={np.nanmin(Z):.1f}–{np.nanmax(Z):.1f} m")
        print(f"    Saved: {out_path}")

    print("\nLambert preprocessing complete.")


if __name__ == '__main__':
    main()
