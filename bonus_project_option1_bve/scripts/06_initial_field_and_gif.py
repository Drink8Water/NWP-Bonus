#!/usr/bin/env python3
"""
06 — Generate a static initial-field figure (ERA5, t = 0 h) and a 5-frame
GIF of the ERA5 500 hPa height-anomaly evolution at 6-hourly intervals
(0 h → 6 h → 12 h → 18 h → 24 h), both on the Lambert conformal grid
with the same curved-boundary Cartopy style as fig13.

Requires: xarray, numpy, scipy, matplotlib, cartopy, pillow
"""

import os
import sys
import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.preprocess import (
    find_data_files, open_dataset, detect_variable_names,
    select_pressure_level, select_times, subset_region,
    coarsen_to_resolution, convert_geopotential, extract_field
)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Reuse the same plotting machinery as fig13
HAS_CARTOPY = False
try:
    import cartopy.crs as ccrs
    HAS_CARTOPY = True
except ImportError:
    pass

if HAS_CARTOPY:
    from src.plotting import (
        make_lambert_rect_axes, plot_lambert_latlon_field,
        LAMBERT_PROJ,
    )

# ── Paths ──────────────────────────────────────────────────────────────────

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'raw')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')

PRESSURE_LEVEL = 500.0
LAT_RANGE = (15.0, 65.0)
LON_RANGE = (60.0, 170.0)
TARGET_DX_DEG = 1.0
TARGET_DY_DEG = 1.0

TIMES = [
    '2025-12-30T00:00:00',
    '2025-12-30T06:00:00',
    '2025-12-30T12:00:00',
    '2025-12-30T18:00:00',
    '2025-12-31T00:00:00',
]

TIME_LABELS = ['0 h', '6 h', '12 h', '18 h', '24 h']


# ── Helpers ─────────────────────────────────────────────────────────────────

def _lonlat_extent(lon2d, lat2d):
    """Return (lon_min, lon_max, lat_min, lat_max)."""
    return (float(np.nanmin(lon2d)), float(np.nanmax(lon2d)),
            float(np.nanmin(lat2d)), float(np.nanmax(lat2d)))


def plot_one_frame(za_anom, lon2d, lat2d, title, filepath,
                   anomaly_interval=30):
    """Plot a single ERA5 height-anomaly field with curved Lambert boundary."""
    max_abs = np.nanmax(np.abs(za_anom))
    max_anom = np.ceil(max_abs / anomaly_interval) * anomaly_interval
    if max_anom < anomaly_interval:
        max_anom = anomaly_interval
    anomaly_levels = np.arange(-max_anom, max_anom + anomaly_interval,
                                anomaly_interval)

    if not HAS_CARTOPY:
        fig, ax = plt.subplots(figsize=(10, 6.5))
        cf = ax.contourf(lon2d, lat2d, za_anom, levels=anomaly_levels,
                         cmap='RdBu_r', extend='both')
        ax.contour(lon2d, lat2d, za_anom, levels=anomaly_levels,
                   colors='k', linewidths=0.3)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.set_xlim(LON_RANGE[0], LON_RANGE[1])
        ax.set_ylim(LAT_RANGE[0], LAT_RANGE[1])
        cbar = fig.colorbar(cf, ax=ax, orientation='horizontal', pad=0.08)
        cbar.set_label('Height anomaly (m)')
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return

    # Cartopy path — use make_lambert_rect_axes exactly like fig13
    lon_min, lon_max, lat_min, lat_max = _lonlat_extent(lon2d, lat2d)

    fig = plt.figure(figsize=(10, 6.5))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[12, 1],
                  hspace=0.08, left=0.06, right=0.94, bottom=0.10, top=0.92)

    # make_lambert_rect_axes handles: projection, boundary, xlim/ylim,
    # coastlines, borders, gridlines, verification box, and clipping
    ax = make_lambert_rect_axes(
        fig, gs[0, 0],
        lon_min=lon_min, lon_max=lon_max,
        lat_min=lat_min, lat_max=lat_max)

    # Plot field — same function as fig13, auto-clips to boundary
    cf = plot_lambert_latlon_field(
        ax, lon2d, lat2d, za_anom, anomaly_levels,
        cmap='RdBu_r', title=title, contour=True, clabel=False)

    # Colorbar
    cax = fig.add_subplot(gs[1, 0])
    cbar = fig.colorbar(cf, cax=cax, orientation='horizontal')
    cbar.set_label('500 hPa height anomaly (m)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    cax.set_facecolor('none')
    for s in cax.spines.values():
        s.set_visible(False)
    cax.set_xticks([])
    cax.set_yticks([])

    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  06 — Initial Field Figure & ERA5 Evolution GIF")
    print("=" * 60)

    os.makedirs(FIG_DIR, exist_ok=True)

    # ── 1. Load ERA5 ────────────────────────────────────────────────────
    files = find_data_files(DATA_DIR, search_parent=True)
    if not files:
        print("ERROR: No .nc file found.")
        sys.exit(1)
    ds = open_dataset(files[0])
    z_name, u_name, v_name = detect_variable_names(ds)
    ds = select_pressure_level(ds, PRESSURE_LEVEL)
    ds_times = select_times(ds, TIMES)
    ds_times = [subset_region(d, LAT_RANGE, LON_RANGE) for d in ds_times]
    ds_times = [coarsen_to_resolution(d, TARGET_DX_DEG, TARGET_DY_DEG)
                for d in ds_times]
    ds_times = [convert_geopotential(d, z_name) for d in ds_times]
    ds_times, z_names = zip(*ds_times)
    z_name = z_names[0]

    # ── 2. Process each time step ───────────────────────────────────────
    #     Plot ERA5 lat-lon data DIRECTLY on Lambert projection — Cartopy
    #     warps the regular geographic grid to fill the curved domain.
    frames = []
    frame_paths = []

    for idx, (ds_t, tlabel) in enumerate(zip(ds_times, TIME_LABELS)):
        print(f"\n  Processing {TIMES[idx]} ({tlabel}) …")

        # Extract and handle N→S ordering
        lat_1d = (ds_t['latitude'].values if 'latitude' in ds_t.coords
                  else ds_t['lat'].values)
        lat_descending = lat_1d[0] > lat_1d[-1]
        if lat_descending:
            lat_1d = np.sort(lat_1d)

        Z_raw = extract_field(ds_t, z_name)
        if lat_descending:
            Z_raw = np.flipud(Z_raw)

        lon_1d = (ds_t['longitude'].values if 'longitude' in ds_t.coords
                  else ds_t['lon'].values)
        lon2d_src, lat2d_src = np.meshgrid(lon_1d, lat_1d)

        # Height anomaly on the native lat-lon grid
        za_anom = Z_raw - np.nanmean(Z_raw)
        frames.append((za_anom, lon2d_src, lat2d_src))

        print(f"    Z anomaly: min={za_anom.min():.1f}, "
              f"max={za_anom.max():.1f}, mean={za_anom.mean():.4f}")

        # ── 4. Plot each frame ────────────────────────────────────────
        za_anom, lon2d_plt, lat2d_plt = frames[idx]
        timestamp = TIMES[idx].replace('T00:00:00', ' 00Z') \
                              .replace('T06:00:00', ' 06Z') \
                              .replace('T12:00:00', ' 12Z') \
                              .replace('T18:00:00', ' 18Z')
        title = f'ERA5 500 hPa Height Anomaly  —  {timestamp}  ({tlabel})'
        frame_path = os.path.join(FIG_DIR,
                                  f'era5_analysis_{tlabel.replace(" ", "h")}.png')
        plot_one_frame(za_anom, lon2d_plt, lat2d_plt, title, frame_path)
        frame_paths.append(frame_path)
        print(f"    Saved: {frame_path}")

    # ── 5. Save static initial-field figure (t = 0 h) ──────────────────
    initial_path = os.path.join(FIG_DIR, 'fig_initial_field_lambert.png')
    # Copy the t=0 frame with a cleaner filename
    import shutil
    shutil.copy(frame_paths[0], initial_path)
    print(f"\n  Static initial-field figure: {initial_path}")

    # ── 6. Build GIF ────────────────────────────────────────────────────
    gif_path = os.path.join(FIG_DIR, 'fig_era5_evolution_0h_24h.gif')
    images = [Image.open(p) for p in frame_paths]
    images[0].save(
        gif_path, save_all=True, append_images=images[1:],
        duration=800, loop=0,
        optimize=False,
    )
    print(f"  GIF saved: {gif_path}")
    print(f"  Frames: {len(images)}, duration: 800 ms each")

    # Clean up individual frames (keep static initial field)
    for p in frame_paths:
        if p != initial_path:
            os.remove(p)

    print("\n" + "=" * 60)
    print("  Done — initial field figure + GIF generated.")
    print("=" * 60)


if __name__ == '__main__':
    main()
