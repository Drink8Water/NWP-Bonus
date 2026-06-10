#!/usr/bin/env python3
"""
08 — Run Lambert-grid BVE experiments and score their optimisation impact.
"""

import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.lambert_grid import load_lambert_grid
from src.bve_model_lambert import (
    LambertBVEModel, initial_condition_from_height, height_anomaly_from_psi
)
from src.poisson_dirichlet import laplacian_dirichlet
from src.verification import rmse, bias, debiased_rmse, anomaly_correlation


DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed_lambert')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'outputs')

DT = 600.0
FORECAST_HOURS = 24.0
SAVE_HOURS = (12.0, 24.0)
VERIF_LAT_RANGE = (25.0, 55.0)
VERIF_LON_RANGE = (90.0, 145.0)

EXPERIMENTS = [
    {'name': 'PERSIST_LCC', 'nu': None, 'sponge_width': None},
    {'name': 'CTRL_LCC', 'nu': 0.0, 'sponge_width': 0},
    {'name': 'DIFF_LCC', 'nu': 2.5e4, 'sponge_width': 0},
    {'name': 'SPONGE_LCC', 'nu': 0.0, 'sponge_width': 5},
    {'name': 'DIFF_SPONGE_LCC', 'nu': 2.5e4, 'sponge_width': 5},
]


def make_verification_mask(grid):
    return ((grid.lat2d >= VERIF_LAT_RANGE[0]) & (grid.lat2d <= VERIF_LAT_RANGE[1])
            & (grid.lon2d >= VERIF_LON_RANGE[0]) & (grid.lon2d <= VERIF_LON_RANGE[1]))


def zeta_from_height_quiet(z, grid):
    z_anom = z - np.mean(z)
    psi = grid.G * z_anom / grid.f
    psi = psi - np.mean(psi[1:-1, 1:-1])
    psi[0, :] = 0.0
    psi[-1, :] = 0.0
    psi[:, 0] = 0.0
    psi[:, -1] = 0.0
    return (grid.m ** 2) * laplacian_dirichlet(psi, grid.dx, grid.dy)


def score_forecast(psi_fcst, zeta_fcst, z_anal, grid, mask):
    zf_anom = height_anomaly_from_psi(psi_fcst, grid)
    za_anom = z_anal - np.mean(z_anal)
    zeta_anal = zeta_from_height_quiet(z_anal, grid)
    scores = {
        'height_rmse_m': float(rmse(zf_anom, za_anom, mask)),
        'height_debiased_rmse_m': float(debiased_rmse(zf_anom, za_anom, mask)),
        'height_acc': float(anomaly_correlation(zf_anom, za_anom, mask)),
        'vorticity_rmse_s-1': float(rmse(zeta_fcst, zeta_anal, mask)),
        'vorticity_acc': float(anomaly_correlation(zeta_fcst, zeta_anal, mask)),
        'height_bias_m': float(bias(zf_anom, za_anom, mask)),
    }
    return scores, {
        'zf_anom': zf_anom,
        'za_anom': za_anom,
        'zeta_fcst': zeta_fcst,
        'zeta_anal': zeta_anal,
    }


def main():
    print("=" * 72)
    print("  Step 08: Lambert BVE Experiment Matrix")
    print("=" * 72)

    os.makedirs(OUT_DIR, exist_ok=True)
    grid = load_lambert_grid(os.path.join(DATA_DIR, 'grid_info_lambert.npz'))
    print(f"Loaded Lambert grid: {grid.nx} x {grid.ny}, dx={grid.dx/1000:.1f} km")

    init = np.load(os.path.join(DATA_DIR, 'initial_20251230_00.npz'))
    anal12 = np.load(os.path.join(DATA_DIR, 'analysis_20251230_12.npz'))
    anal24 = np.load(os.path.join(DATA_DIR, 'analysis_20251231_00.npz'))
    Z0 = init['Z']
    Z_ANAL = {12.0: anal12['Z'], 24.0: anal24['Z']}

    psi0, zeta0 = initial_condition_from_height(Z0, grid)
    mask = make_verification_mask(grid)
    print(f"Verification points: {np.sum(mask)}")

    all_scores = {}

    for exp in EXPERIMENTS:
        name = exp['name']
        print(f"\n{'-' * 58}")
        print(f"Experiment: {name}")
        print(f"{'-' * 58}")

        if name == 'PERSIST_LCC':
            states = {12.0: (psi0, zeta0), 24.0: (psi0, zeta0)}
        else:
            model = LambertBVEModel(
                grid, dt=DT,
                nu=exp['nu'],
                sponge_width=exp['sponge_width'],
                sponge_tau_hours=6.0,
                zeta_ref=zeta0.copy(),
            )
            history = model.forecast(psi0.copy(), zeta0.copy(), FORECAST_HOURS, SAVE_HOURS)
            states = {
                float(t): (history['psi'][i], history['zeta'][i])
                for i, t in enumerate(history['time_hours'])
            }

        for lead in SAVE_HOURS:
            psi_fcst, zeta_fcst = states[lead]
            scores, arrays = score_forecast(psi_fcst, zeta_fcst, Z_ANAL[lead], grid, mask)
            key = f"{name}_{int(lead)}h"
            all_scores[key] = scores
            print(f"  +{int(lead)}h RMSE={scores['height_rmse_m']:.1f} m,"
                  f" debRMSE={scores['height_debiased_rmse_m']:.1f} m,"
                  f" bias={scores['height_bias_m']:.1f} m,"
                  f" ACC={scores['height_acc']:.3f}")

            np.savez(os.path.join(OUT_DIR, f'{name}_verification_{int(lead)}h.npz'),
                     **arrays)
            np.savez(os.path.join(OUT_DIR, f'{name}_forecast_{int(lead)}h.npz'),
                     psi=psi_fcst, zeta=zeta_fcst, time_hours=lead, experiment=name)

    json_path = os.path.join(OUT_DIR, 'scores_lambert.json')
    with open(json_path, 'w') as f:
        json.dump(all_scores, f, indent=2)
    print(f"\nSaved: {json_path}")

    print("\nSummary:")
    print(f"{'Experiment':<18} {'Lead':>5} {'RMSE':>8} {'DebRMSE':>8} {'Bias':>8} {'ACC':>7}")
    for exp in EXPERIMENTS:
        for lead in SAVE_HOURS:
            key = f"{exp['name']}_{int(lead)}h"
            s = all_scores[key]
            print(f"{exp['name']:<18} +{int(lead):<4d}"
                  f" {s['height_rmse_m']:8.1f}"
                  f" {s['height_debiased_rmse_m']:8.1f}"
                  f" {s['height_bias_m']:8.1f}"
                  f" {s['height_acc']:7.3f}")


if __name__ == '__main__':
    main()
