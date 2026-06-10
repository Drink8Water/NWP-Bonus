#!/usr/bin/env python3
"""04 - Generate figures for the Lambert finite-area BVE experiments."""

import csv
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.plotting import plot_forecast_verification


DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed_lambert')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'outputs')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')


def load_grid():
    data = np.load(os.path.join(DATA_DIR, 'grid_info_lambert.npz'))
    proj = ccrs.LambertConformal(
        central_longitude=float(data['central_lon']),
        central_latitude=float(data['central_lat']),
        standard_parallels=tuple(data['standard_parallels']),
    )
    return {
        'X': data['X'],
        'Y': data['Y'],
        'lon2d': data['lon2d'],
        'lat2d': data['lat2d'],
        'proj': proj,
    }


def load_source_display_grid():
    source = np.load(os.path.join(os.path.dirname(__file__), '..',
                                  'data', 'processed', 'initial_20251230_00.npz'))
    lon2d, lat2d = np.meshgrid(source['lon'], source['lat'])
    return lon2d, lat2d


def interpolate_lambert_field_to_lonlat(field, grid, target_lon2d, target_lat2d,
                                        pad_cells=16):
    """
    Bilinear interpolation from regular Lambert X/Y grid to lon-lat targets.

    For display only, the Lambert model field is padded with a smooth linear
    ramp to avoid a hard rectangular seam when plotting on the larger
    lat-lon footprint used by the report figures.
    """
    xy = grid['proj'].transform_points(ccrs.PlateCarree(), target_lon2d, target_lat2d)
    xt = xy[..., 0]
    yt = xy[..., 1]

    x = grid['X'][0, :]
    y = grid['Y'][:, 0]
    src = np.asarray(field)

    if pad_cells > 0:
        mean_val = float(np.nanmean(src))
        src = np.pad(
            src, pad_cells, mode='linear_ramp',
            end_values=((mean_val, mean_val), (mean_val, mean_val)))
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        x = x[0] - pad_cells * dx + np.arange(src.shape[1]) * dx
        y = y[0] - pad_cells * dy + np.arange(src.shape[0]) * dy

    xt = np.clip(xt, x[0], x[-1])
    yt = np.clip(yt, y[0], y[-1])

    ix1 = np.searchsorted(x, xt, side='right')
    iy1 = np.searchsorted(y, yt, side='right')
    ix1 = np.clip(ix1, 1, len(x) - 1)
    iy1 = np.clip(iy1, 1, len(y) - 1)
    ix0 = ix1 - 1
    iy0 = iy1 - 1

    x0 = x[ix0]
    x1 = x[ix1]
    y0 = y[iy0]
    y1 = y[iy1]
    wx = np.divide(xt - x0, x1 - x0, out=np.zeros_like(xt), where=(x1 != x0))
    wy = np.divide(yt - y0, y1 - y0, out=np.zeros_like(yt), where=(y1 != y0))

    f00 = src[iy0, ix0]
    f10 = src[iy0, ix1]
    f01 = src[iy1, ix0]
    f11 = src[iy1, ix1]
    return ((1 - wx) * (1 - wy) * f00
            + wx * (1 - wy) * f10
            + (1 - wx) * wy * f01
            + wx * wy * f11)


def lambert_display_weight(grid, target_lon2d, target_lat2d, taper_cells=8):
    """Smoothly taper Lambert-model fields near the native grid boundary."""
    xy = grid['proj'].transform_points(ccrs.PlateCarree(), target_lon2d, target_lat2d)
    xt = xy[..., 0]
    yt = xy[..., 1]
    x_min, x_max = float(np.min(grid['X'])), float(np.max(grid['X']))
    y_min, y_max = float(np.min(grid['Y'])), float(np.max(grid['Y']))
    dist_inside = np.minimum.reduce([
        xt - x_min,
        x_max - xt,
        yt - y_min,
        y_max - yt,
    ])
    taper = taper_cells * float(grid['X'][0, 1] - grid['X'][0, 0])
    weight = np.clip(dist_inside / taper, 0.0, 1.0)
    return weight * weight * (3.0 - 2.0 * weight)


def make_ctrl_verification_figures(grid):
    lon2d, lat2d = load_source_display_grid()
    for lead in (12, 24):
        data = np.load(os.path.join(OUT_DIR, f'CTRL_LCC_verification_{lead}h.npz'))
        analysis_name = 'analysis_20251230_12.npz' if lead == 12 else 'analysis_20251231_00.npz'
        initial = np.load(os.path.join(os.path.dirname(__file__), '..',
                                       'data', 'processed', 'initial_20251230_00.npz'))
        analysis = np.load(os.path.join(os.path.dirname(__file__), '..',
                                        'data', 'processed', analysis_name))
        zf_lcc = interpolate_lambert_field_to_lonlat(data['zf_anom'], grid, lon2d, lat2d)
        z_persist = initial['Z'] - np.mean(initial['Z'])
        w = lambert_display_weight(grid, lon2d, lat2d, taper_cells=8)
        zf = w * zf_lcc + (1.0 - w) * z_persist
        za = analysis['Z'] - np.mean(analysis['Z'])
        plot_forecast_verification(
            zf, za, lon2d, lat2d,
            title_prefix=f'Lambert CTRL +{lead} h',
            filepath=os.path.join(FIG_DIR, f'fig13_lambert_ctrl_verification_{lead}h.png'),
            contour_interval=30,
            error_interval=15,
        )


def make_lambert_score_figure():
    with open(os.path.join(OUT_DIR, 'scores_experiment_matrix.json')) as f:
        scores = json.load(f)

    experiments = ['PERSIST_LCC', 'CTRL_LCC', 'DIFF_LCC', 'SPONGE_LCC', 'DIFF_SPONGE_LCC']
    labels = ['PERSIST', 'CTRL', 'DIFF', 'SPONGE', 'DIFF+SP']
    x = np.arange(len(labels))
    width = 0.36

    rmse12 = [scores[f'{exp}_12h']['height_rmse_m'] for exp in experiments]
    rmse24 = [scores[f'{exp}_24h']['height_rmse_m'] for exp in experiments]
    acc12 = [scores[f'{exp}_12h']['height_acc'] for exp in experiments]
    acc24 = [scores[f'{exp}_24h']['height_acc'] for exp in experiments]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    axes[0].bar(x - width / 2, rmse12, width, label='+12 h', color='#4472C4',
                edgecolor='k', linewidth=0.5)
    axes[0].bar(x + width / 2, rmse24, width, label='+24 h', color='#ED7D31',
                edgecolor='k', linewidth=0.5)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=20, ha='right')
    axes[0].set_ylabel('Height anomaly RMSE (m)')
    axes[0].set_title('Lambert Experiment RMSE', fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    axes[0].legend()

    axes[1].bar(x - width / 2, acc12, width, label='+12 h', color='#4472C4',
                edgecolor='k', linewidth=0.5)
    axes[1].bar(x + width / 2, acc24, width, label='+24 h', color='#ED7D31',
                edgecolor='k', linewidth=0.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=20, ha='right')
    axes[1].set_ylabel('Height anomaly correlation')
    axes[1].set_ylim(0.9, 1.0)
    axes[1].set_title('Lambert Experiment ACC', fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    axes[1].legend()

    fig.suptitle('Lambert Finite-Area BVE Experiment Scores',
                 fontsize=14, fontweight='bold')
    path = os.path.join(FIG_DIR, 'fig14_lambert_experiment_scores.png')
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def make_impact_figure():
    csv_path = os.path.join(OUT_DIR, 'lambert_impact_summary.csv')
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row['lead_h']) == 24:
                rows.append(row)

    labels = [r['experiment'] for r in rows]
    beta_rmse = [float(r['beta_rmse_m']) for r in rows]
    lcc_rmse = [float(r['lambert_rmse_m']) for r in rows]
    reduction = [float(r['rmse_reduction_pct']) for r in rows]
    x = np.arange(len(labels))
    width = 0.36

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8), constrained_layout=True)
    axes[0].bar(x - width / 2, beta_rmse, width, label='Beta-plane', color='#7F7F7F',
                edgecolor='k', linewidth=0.5)
    axes[0].bar(x + width / 2, lcc_rmse, width, label='Lambert', color='#70AD47',
                edgecolor='k', linewidth=0.5)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=20, ha='right')
    axes[0].set_ylabel('+24 h height RMSE (m)')
    axes[0].set_title('Beta-Plane vs Lambert', fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    axes[0].legend()

    colors = ['#70AD47' if v >= 0 else '#C00000' for v in reduction]
    axes[1].bar(x, reduction, color=colors, edgecolor='k', linewidth=0.5)
    axes[1].axhline(0.0, color='k', linewidth=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=20, ha='right')
    axes[1].set_ylabel('RMSE reduction (%)')
    axes[1].set_title('Lambert Impact at +24 h', fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)

    fig.suptitle('Impact of Lambert Grid and Fixed Boundaries',
                 fontsize=14, fontweight='bold')
    path = os.path.join(FIG_DIR, 'fig15_lambert_impact_summary.png')
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def main():
    print('=' * 64)
    print('  Step 10: Lambert Figures')
    print('=' * 64)
    os.makedirs(FIG_DIR, exist_ok=True)
    grid = load_grid()
    make_ctrl_verification_figures(grid)
    make_lambert_score_figure()
    make_impact_figure()
    print('Lambert figures complete.')


if __name__ == '__main__':
    main()
