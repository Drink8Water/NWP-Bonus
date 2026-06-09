#!/usr/bin/env python3
"""
01 — Preprocess local ERA5 data for the BVE forecast experiment.

Reads NetCDF files from data/raw/, selects 500 hPa, extracts three
time steps (initial + two verification times), subsets the East Asia
domain, coarsens to ~1° resolution, converts geopotential to height,
and saves processed .npz files.
"""

import os
import sys
import numpy as np

# Add parent to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.preprocess import (
    find_data_files, open_dataset, detect_variable_names,
    select_pressure_level, select_times, subset_region,
    coarsen_to_resolution, convert_geopotential, extract_field
)
from src.grid import BetaPlaneGrid

# ── Configuration ──────────────────────────────────────────────────────────

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'raw')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')

PRESSURE_LEVEL = 500.0  # hPa

TIMES = [
    '2025-12-30T00:00:00',
    '2025-12-30T12:00:00',
    '2025-12-31T00:00:00',
]

OUTPUT_NAMES = [
    'initial_20251230_00.npz',
    'analysis_20251230_12.npz',
    'analysis_20251231_00.npz',
]

# Domain and grid settings
LAT_RANGE = (15.0, 65.0)
LON_RANGE = (60.0, 170.0)
TARGET_DX_DEG = 1.0
TARGET_DY_DEG = 1.0

# Beta-plane reference
LAT0 = 40.0
LON0 = 115.0  # approximate centre of East Asia domain


def main():
    print("=" * 60)
    print("  Step 01: Preprocess Local ERA5 Data")
    print("=" * 60)

    # ── 1. Find data files ─────────────────────────────────────────────
    files = find_data_files(DATA_DIR, search_parent=True)
    if not files:
        print(f"\nERROR: No .nc files found in {DATA_DIR} or its parent.")
        print("Please place your ERA5 NetCDF file in data/raw/")
        sys.exit(1)

    print(f"\nFound {len(files)} NetCDF file(s):")
    for f in files:
        print(f"  {f}")

    # ── 2. Open dataset ────────────────────────────────────────────────
    ds = open_dataset(files[0])

    # ── 3. Detect variable names ───────────────────────────────────────
    z_name, u_name, v_name = detect_variable_names(ds)
    print(f"  Geopotential variable: '{z_name}'")
    print(f"  U-wind variable:       '{u_name}'")
    print(f"  V-wind variable:       '{v_name}'")

    # ── 4. Select pressure level ───────────────────────────────────────
    ds = select_pressure_level(ds, PRESSURE_LEVEL)

    # ── 5. Select times ────────────────────────────────────────────────
    ds_times = select_times(ds, TIMES)

    # ── 6. Subset region ───────────────────────────────────────────────
    ds_times = [subset_region(d, LAT_RANGE, LON_RANGE) for d in ds_times]

    # ── 7. Coarsen to ~1° resolution ───────────────────────────────────
    ds_times = [coarsen_to_resolution(d, TARGET_DX_DEG, TARGET_DY_DEG)
                for d in ds_times]

    # ── 8. Convert geopotential to height ──────────────────────────────
    ds_times = [convert_geopotential(d, z_name) for d in ds_times]
    # Each element is now (ds, z_name_used)
    ds_times, z_names = zip(*ds_times)
    z_name = z_names[0]  # Should be same for all

    # ── 9. Extract fields from first dataset to determine grid shape ──
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Get the actual lat/lon from the coarsened data
    ds0 = ds_times[0]
    lat_1d = ds0['latitude'].values if 'latitude' in ds0.coords \
             else ds0['lat'].values
    lon_1d = ds0['longitude'].values if 'longitude' in ds0.coords \
             else ds0['lon'].values

    # Detect if latitudes are descending (N-to-S, typical ERA5)
    lat_descending = lat_1d[0] > lat_1d[-1]
    if lat_descending:
        print("  Latitudes are N→S; will flip to S→N for model grid consistency.")
        # Sort to increasing order
        lat_1d = np.sort(lat_1d)
        # We will flip the extracted arrays below

    actual_lat_range = (float(lat_1d[0]), float(lat_1d[-1]))
    actual_lon_range = (float(lon_1d[0]), float(lon_1d[-1]))
    actual_dx_deg = float(abs(lon_1d[1] - lon_1d[0])) if len(lon_1d) > 1 else TARGET_DX_DEG
    actual_dy_deg = float(abs(lat_1d[1] - lat_1d[0])) if len(lat_1d) > 1 else TARGET_DY_DEG

    print(f"\n  Actual data grid: lat={actual_lat_range}, lon={actual_lon_range}")
    print(f"  Actual resolution: dy={actual_dy_deg:.3f}°, dx={actual_dx_deg:.3f}°")
    print(f"  Data shape: {len(lat_1d)} x {len(lon_1d)}")

    # Create model grid matching the actual data
    grid = BetaPlaneGrid(lat0=LAT0, lon0=LON0,
                          lat_range=actual_lat_range, lon_range=actual_lon_range,
                          dx_deg=actual_dx_deg, dy_deg=actual_dy_deg)

    for idx, (ds_t, out_name) in enumerate(zip(ds_times, OUTPUT_NAMES)):
        # Extract 2-D fields (flip if data was N→S)
        Z = extract_field(ds_t, z_name)
        U = extract_field(ds_t, u_name)
        V = extract_field(ds_t, v_name)
        if lat_descending:
            Z = np.flipud(Z)
            U = np.flipud(U)
            V = np.flipud(V)

        print(f"\n  Time {idx} ({TIMES[idx]}):")
        print(f"    Z: shape={Z.shape}, min={Z.min():.1f}, max={Z.max():.1f}, "
              f"mean={Z.mean():.1f} m")
        print(f"    U: shape={U.shape}, min={U.min():.1f}, max={U.max():.1f}, "
              f"mean={U.mean():.1f} m/s")
        print(f"    V: shape={V.shape}, min={V.min():.1f}, max={V.max():.1f}, "
              f"mean={V.mean():.1f} m/s")

        out_path = os.path.join(OUTPUT_DIR, out_name)
        np.savez(out_path,
                 Z=Z, U=U, V=V,
                 lat=lat_1d, lon=lon_1d,
                 time=TIMES[idx],
                 pressure_level=PRESSURE_LEVEL,
                 description=f"ERA5 500 hPa — {TIMES[idx]}")
        print(f"    Saved: {out_path}")

    # Save grid info too
    grid_path = os.path.join(OUTPUT_DIR, 'grid_info.npz')
    np.savez(grid_path,
             lat0=LAT0, lon0=LON0,
             lat_range=np.array(actual_lat_range), lon_range=np.array(actual_lon_range),
             dx_deg=actual_dx_deg, dy_deg=actual_dy_deg,
             f0=grid.f0, beta=grid.beta,
             dx_m=grid.dx, dy_m=grid.dy,
             ny=grid.ny, nx=grid.nx)
    print(f"\n  Saved grid info: {grid_path}")

    print("\n" + "=" * 60)
    print("  Preprocessing complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
