#!/usr/bin/env python3
"""
06 — Sensitivity experiment comparison figures.

Uses the same rectangular Lambert domain as the main figures.
"""

import os, sys, json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from src.grid import BetaPlaneGrid
from src.plotting import (
    make_lambert_rect_axes, _get_shared_domain,
    plot_lambert_gridded_field,
    get_lambert_projection,
    interpolate_fields_to_lambert_rectangle,
    _save_and_close, _symmetric_levels,
)
LAMBERT_PROJ = get_lambert_projection()

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'outputs')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')

EXPERIMENT_NAMES = ['PERSIST', 'CTRL', 'DIFF', 'SPONGE', 'DIFF_SPONGE']
EXPERIMENT_LABELS = ['PERSIST', 'CTRL', 'DIFF', 'SPONGE', 'DIFF+SPONGE']


def load_grid():
    gd = np.load(os.path.join(DATA_DIR, 'grid_info.npz'))
    return BetaPlaneGrid(
        lat0=float(gd['lat0']), lon0=float(gd['lon0']),
        lat_range=tuple(gd['lat_range']), lon_range=tuple(gd['lon_range']),
        dx_deg=float(gd['dx_deg']), dy_deg=float(gd['dy_deg']))


def load_verification(exp_name, lead_h):
    p = os.path.join(OUTPUT_DIR, f'{exp_name}_verification_{lead_h}h.npz')
    return np.load(p) if os.path.exists(p) else None


# ══════════════════════════════════════════════════════════════════════
#  Figure 1: Score comparison
# ══════════════════════════════════════════════════════════════════════

def plot_score_comparison(scores, filepath):
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    metrics_config = [
        ('rmse_m', 'Height Anomaly RMSE (m)', axes[0, 0], False),
        ('debiased_rmse_m', 'Debiased RMSE (m)', axes[0, 1], False),
        ('acc', 'Anomaly Correlation', axes[1, 0], True),
        ('bias_m', 'Height Bias (m)', axes[1, 1], False),
    ]
    for metric, ylabel, ax, higher_better in metrics_config:
        x = np.arange(len(EXPERIMENT_NAMES)); w = 0.35
        v12 = [scores.get(f'{e}_12h', {}).get(metric, np.nan) for e in EXPERIMENT_NAMES]
        v24 = [scores.get(f'{e}_24h', {}).get(metric, np.nan) for e in EXPERIMENT_NAMES]
        b12 = ax.bar(x - w/2, v12, w, label='+12 h', color='#4472C4', edgecolor='k', linewidth=0.5)
        b24 = ax.bar(x + w/2, v24, w, label='+24 h', color='#ED7D31', edgecolor='k', linewidth=0.5)
        for b in [b12, b24]:
            for bar in b:
                h = bar.get_height()
                if not np.isnan(h):
                    ax.text(bar.get_x() + bar.get_width()/2, h + 0.02*max(1, abs(h)),
                            f'{h:.1f}' if not higher_better else f'{h:.3f}',
                            ha='center', va='bottom' if h>=0 else 'top', fontsize=6, rotation=90)
        ax.set_xticks(x); ax.set_xticklabels(EXPERIMENT_LABELS, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=11); ax.legend(fontsize=9)
        ax.grid(axis='y', alpha=0.3)
        if higher_better: ax.set_ylim(0, 1.05)
        ax.set_title(ylabel, fontsize=12, fontweight='bold')
    fig.suptitle('Sensitivity Experiment: Score Comparison', fontsize=14, fontweight='bold')
    _save_and_close(fig, filepath)


# ══════════════════════════════════════════════════════════════════════
#  Figure 2: +24 h error comparison
# ══════════════════════════════════════════════════════════════════════

def plot_error_comparison(grid, filepath):
    exps = ['CTRL', 'DIFF', 'SPONGE', 'DIFF_SPONGE']
    descs = ['CTRL: BVE only', 'DIFF: BVE + diffusion',
             'SPONGE: BVE + sponge', 'DIFF+SPONGE: combined']
    X_lcc, Y_lcc, lon_t, lat_t, xlim, ylim = _get_shared_domain()

    # Collect errors and interpolate
    err_dict = {}
    for exp in exps:
        v = load_verification(exp, 24)
        if v is not None:
            err_dict[exp] = v['zf_anom'] - v['za_anom']
    lcc_errs = interpolate_fields_to_lambert_rectangle(
        err_dict, grid.lon2d, grid.lat2d, lon_t, lat_t)

    all_errs = [v for v in lcc_errs.values() if np.any(np.isfinite(v))]
    max_abs = max(np.nanmax(np.abs(v)) for v in all_errs)
    interval = 20
    emax = np.ceil(max_abs / interval) * interval
    error_levels = np.arange(-emax, emax + interval, interval)

    fig = plt.figure(figsize=(19, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.22, wspace=0.10,
                  left=0.03, right=0.97, bottom=0.08, top=0.92)
    for idx, (exp, desc) in enumerate(zip(exps, descs)):
        ax = make_lambert_rect_axes(fig, gs[idx//2, idx%2], xlim, ylim)
        if exp not in lcc_errs:
            ax.set_title(f'{exp}: (not found)', fontsize=10); continue
        cf = plot_lambert_gridded_field(ax, X_lcc, Y_lcc, lcc_errs[exp],
                                        levels=error_levels, cmap='RdBu_r',
                                        title=desc, contour=True, clabel=False)
        cbar = fig.colorbar(cf, ax=ax, orientation='horizontal', pad=0.06, shrink=0.82)
        cbar.set_label('Error (m)', fontsize=9); cbar.ax.tick_params(labelsize=8)
    fig.suptitle('+24 h Error Comparison: CTRL vs Sensitivity Experiments',
                 fontsize=14, fontweight='bold')
    _save_and_close(fig, filepath)


# ══════════════════════════════════════════════════════════════════════
#  Figure 3: +24 h improvement maps
# ══════════════════════════════════════════════════════════════════════

def plot_improvement_maps(grid, filepath):
    ctrl_v = load_verification('CTRL', 24)
    if ctrl_v is None: return
    ctrl_err = np.abs(ctrl_v['zf_anom'] - ctrl_v['za_anom'])

    exps = ['DIFF', 'SPONGE', 'DIFF_SPONGE']
    descs = ['DIFF − CTRL', 'SPONGE − CTRL', 'DIFF+SPONGE − CTRL']
    X_lcc, Y_lcc, lon_t, lat_t, xlim, ylim = _get_shared_domain()

    imp_dict = {}
    for exp in exps:
        v = load_verification(exp, 24)
        if v is not None:
            imp_dict[exp] = ctrl_err - np.abs(v['zf_anom'] - v['za_anom'])
    lcc_imps = interpolate_fields_to_lambert_rectangle(
        imp_dict, grid.lon2d, grid.lat2d, lon_t, lat_t)

    all_vals = [v for v in lcc_imps.values() if np.any(np.isfinite(v))]
    max_imp = max(np.nanmax(np.abs(v)) for v in all_vals) if all_vals else 100
    interval = max(5, int(np.ceil(max_imp / 10)))
    imp_levels = np.arange(-max_imp, max_imp + interval, interval)

    fig = plt.figure(figsize=(17, 5.5))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[10, 1],
                  hspace=0.15, wspace=0.12,
                  left=0.03, right=0.97, bottom=0.12, top=0.90)

    for idx, (exp, desc) in enumerate(zip(exps, descs)):
        ax = make_lambert_rect_axes(fig, gs[0, idx], xlim, ylim)
        if exp not in lcc_imps:
            ax.set_title(f'{exp}: (not found)', fontsize=10); continue
        cf = plot_lambert_gridded_field(ax, X_lcc, Y_lcc, lcc_imps[exp],
                                        levels=imp_levels, cmap='BrBG',
                                        title=desc, contour=True, clabel=False)
        ax.contour(X_lcc, Y_lcc, lcc_imps[exp], levels=[0],
                   colors='k', linewidths=0.8, transform=LAMBERT_PROJ)
        cax = fig.add_subplot(gs[1, idx])
        cbar = fig.colorbar(cf, cax=cax, orientation='horizontal')
        cbar.set_label('|Error CTRL| − |Error EXP| (m)', fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        for sp in cax.spines.values(): sp.set_visible(False)
        cax.set_xticks([]); cax.set_yticks([])

    fig.suptitle('+24 h Improvement over CTRL (positive = better)',
                 fontsize=14, fontweight='bold', y=0.97)
    _save_and_close(fig, filepath)


# ══════════════════════════════════════════════════════════════════════
#  Figure 4: +24 h vorticity CTRL vs DIFF
# ══════════════════════════════════════════════════════════════════════

def plot_vorticity_comparison(grid, filepath):
    ctrl_v = load_verification('CTRL', 24)
    diff_v = load_verification('DIFF', 24)
    if ctrl_v is None or diff_v is None: return

    X_lcc, Y_lcc, lon_t, lat_t, xlim, ylim = _get_shared_domain()
    zeta_dict = {
        'ctrl': ctrl_v['zeta_fcst'] * 1e5,
        'diff': diff_v['zeta_fcst'] * 1e5,
    }
    lcc_zeta = interpolate_fields_to_lambert_rectangle(
        zeta_dict, grid.lon2d, grid.lat2d, lon_t, lat_t)

    zeta_max = max(np.nanmax(np.abs(v)) for v in lcc_zeta.values())
    zeta_max = max(np.ceil(zeta_max * 2) / 2, 1)
    zeta_levels = np.linspace(-zeta_max, zeta_max, 21)

    titles = ['CTRL: BVE only', 'DIFF: BVE + diffusion (ν = 2.5×10⁴ m²/s)']
    keys = ['ctrl', 'diff']

    fig = plt.figure(figsize=(14, 5.2))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[10, 1],
                  hspace=0.15, wspace=0.15,
                  left=0.03, right=0.97, bottom=0.15, top=0.88)

    for idx, (key, title) in enumerate(zip(keys, titles)):
        ax = make_lambert_rect_axes(fig, gs[0, idx], xlim, ylim)
        cf = plot_lambert_gridded_field(ax, X_lcc, Y_lcc, lcc_zeta[key],
                                        levels=zeta_levels, cmap='RdBu_r',
                                        title=title, contour=True, clabel=False)
        cax = fig.add_subplot(gs[1, idx])
        cbar = fig.colorbar(cf, cax=cax, orientation='horizontal')
        cbar.set_label('Relative vorticity (10⁻⁵ s⁻¹)', fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        for sp in cax.spines.values(): sp.set_visible(False)
        cax.set_xticks([]); cax.set_yticks([])

    fig.suptitle('+24 h Relative Vorticity: CTRL vs DIFF',
                 fontsize=14, fontweight='bold', y=0.97)
    _save_and_close(fig, filepath)


# ══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  Sensitivity Figures")
    print("=" * 60)
    grid = load_grid()
    scores_path = os.path.join(OUTPUT_DIR, 'scores_experiment_matrix.json')
    if not os.path.exists(scores_path):
        print(f"ERROR: {scores_path} not found. Run 05 first."); sys.exit(1)
    with open(scores_path) as f:
        scores = json.load(f)
    os.makedirs(FIG_DIR, exist_ok=True)

    print("\n[1/4] Score comparison bar chart")
    plot_score_comparison(scores, os.path.join(FIG_DIR, 'fig06_sensitivity_scores.png'))
    print("\n[2/4] +24 h error comparison maps")
    plot_error_comparison(grid, os.path.join(FIG_DIR, 'fig07_error_comparison_24h.png'))
    print("\n[3/4] +24 h improvement over CTRL")
    plot_improvement_maps(grid, os.path.join(FIG_DIR, 'fig08_improvement_24h.png'))
    print("\n[4/4] +24 h vorticity CTRL vs DIFF")
    plot_vorticity_comparison(grid, os.path.join(FIG_DIR, 'fig09_vorticity_ctrl_vs_diff_24h.png'))

    print("\n" + "=" * 60)
    print("  Sensitivity figures complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
