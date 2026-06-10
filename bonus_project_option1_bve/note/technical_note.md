# Technical Note: Lambert Finite-Area 500 hPa BVE Forecast

**Bonus Project Option 1 — Numerical Weather Prediction Course**

---

## 1. Case And Data

The case is a winter East Asian 500 hPa flow pattern initialized at
2025-12-30 00 UTC. ERA5 pressure-level reanalysis is used for the initial
condition and verifying analyses.

Three times are used:

- 2025-12-30 00 UTC: initial condition
- 2025-12-30 12 UTC: +12 h verifying analysis
- 2025-12-31 00 UTC: +24 h verifying analysis

The raw ERA5 fields are subset to 15°N–65°N, 60°E–170°E and coarsened to about
1° resolution. Geopotential is converted to geopotential height by
`Z = Φ/g`.

---

## 2. Lambert Finite-Area Model

### 2.1 Governing Equation

The current dynamical core is a finite-area barotropic vorticity equation
model on a Lambert conformal projection grid:

$$
\frac{\partial \zeta}{\partial t}
= -m^2 J(\psi,\zeta+f)
  + \nu m^2\nabla^2\zeta
  - \alpha(\zeta-\zeta_0),
\qquad
\zeta = m^2\nabla^2\psi.
$$

Here `ψ` is streamfunction, `ζ` is relative vorticity, `m(i,j)` is the Lambert
map factor, and `f(i,j)=2Ωsinφ` is the Coriolis parameter evaluated on the
model grid. The CTRL experiment sets `ν=0` and `α=0`; DIFF and SPONGE switch on
the optional diffusion and boundary relaxation terms.

### 2.2 Grid

The Lambert grid is generated in `src/lambert_grid.py`.

| Quantity | Value |
| --- | --- |
| Projection | Lambert conformal |
| Standard parallels | 25°N / 45°N |
| Grid spacing | 150 km |
| Grid size | 54 × 34 |
| Approximate coverage | 56°E–174°E, 12°N–62°N |
| Verification region | 25°N–55°N, 90°E–145°E |

ERA5 fields are first processed on a regular latitude-longitude grid and then
bilinearly interpolated to the Lambert model grid by
`scripts/02_prepare_lambert_grid.py`.

### 2.3 Initial Condition

The initial streamfunction is diagnosed from height anomaly with a variable-f
geostrophic relation:

$$
\psi_0 = \frac{g(Z_0-\overline{Z_0})}{f}.
$$

The boundary value of `ψ` is set to zero to match the fixed-boundary Poisson
solver. Initial vorticity is then computed as

$$
\zeta_0 = m^2\nabla^2\psi_0.
$$

### 2.4 Poisson Solver

The inversion

$$
\nabla^2\psi = \frac{\zeta}{m^2}
$$

is solved with zero streamfunction on all four boundaries. The implementation
uses a discrete sine transform written with NumPy FFTs in
`src/poisson_dirichlet.py`. This avoids the doubly periodic assumption of a
standard FFT Poisson solver.

---

## 3. Numerical Method

- Spatial derivatives use second-order centered differences.
- The nonlinear Jacobian uses the Arakawa scheme.
- Time integration uses fourth-order Runge-Kutta with `Δt = 600 s`.
- Optional diffusion uses Laplacian vorticity diffusion with
  `ν = 2.5 × 10^4 m² s⁻¹`.
- Optional sponge relaxation damps vorticity toward the initial vorticity near
  the lateral boundaries.

---

## 4. Height Recovery And Verification

The model predicts vorticity and streamfunction. Forecast height anomaly is
diagnosed by

$$
Z'_\text{fcst} = \frac{f\psi_\text{fcst}}{g}.
$$

The model does not predict absolute domain-mean height, so full-field RMSE is
not reported. Verification focuses on:

| Metric | Meaning |
| --- | --- |
| Height anomaly RMSE | RMSE of `Z'_fcst - Z'_anal` |
| Height debiased RMSE | RMSE after removing mean forecast-analysis bias |
| Height ACC | Spatial anomaly correlation |
| Height bias | Mean forecast-analysis height anomaly difference |
| Vorticity correlation | Spatial correlation of forecast and analysis vorticity |

All scores are computed on the inner verification region
25°N–55°N, 90°E–145°E.

---

## 5. Experiments

The Lambert experiment matrix is run by `scripts/03_run_experiments.py`.

| Experiment | Description |
| --- | --- |
| PERSIST_LCC | Persistence baseline: initial height anomaly remains fixed |
| CTRL_LCC | Lambert BVE without diffusion or sponge |
| DIFF_LCC | Lambert BVE + Laplacian vorticity diffusion |
| SPONGE_LCC | Lambert BVE + lateral sponge relaxation |
| DIFF_SPONGE_LCC | Lambert BVE + diffusion + sponge |

### 5.1 Results

| Experiment | Lead | RMSE (m) | Debiased RMSE (m) | Bias (m) | ACC |
| --- | --- | ---: | ---: | ---: | ---: |
| PERSIST_LCC | +12 h | 42.9 | 41.4 | -11.1 | 0.985 |
| PERSIST_LCC | +24 h | 70.3 | 69.2 | -12.5 | 0.953 |
| CTRL_LCC | +12 h | 55.4 | 54.1 | +12.1 | 0.975 |
| CTRL_LCC | +24 h | 56.6 | 55.7 | +10.0 | 0.963 |
| DIFF_LCC | +12 h | 55.9 | 55.1 | +9.7 | 0.973 |
| DIFF_LCC | +24 h | 57.9 | 57.9 | -1.0 | 0.960 |
| SPONGE_LCC | +12 h | 50.8 | 50.7 | +2.6 | 0.980 |
| SPONGE_LCC | +24 h | 69.0 | 69.0 | +1.6 | 0.960 |
| DIFF_SPONGE_LCC | +12 h | 52.2 | 52.1 | +1.4 | 0.978 |
| DIFF_SPONGE_LCC | +24 h | 71.9 | 71.8 | -1.7 | 0.956 |

### 5.2 Interpretation

The case is highly persistent: PERSIST_LCC already gives +24 h ACC of 0.953.
CTRL_LCC nevertheless improves +24 h RMSE to 56.6 m and ACC to 0.963. In this
Lambert finite-area configuration, diffusion and sponge relaxation no longer
provide the large RMSE reduction seen in earlier periodic-boundary experiments;
they mainly alter bias and boundary behaviour.

The best +24 h RMSE among integrated Lambert experiments is CTRL_LCC. The best
+12 h RMSE is SPONGE_LCC, but the sponge does not improve +24 h RMSE.

---

## 6. Figures

`src/plotting.py` defines a curved Lambert map boundary with
`set_lambert_curved_boundary()`. The boundary is constructed from the projected
latitude-longitude source-domain outline:

- bottom edge: fixed `lat = lat_min`
- right edge: fixed `lon = lon_max`
- top edge: fixed `lat = lat_max`
- left edge: fixed `lon = lon_min`

The resulting boundary is applied with

```python
ax.set_boundary(boundary_path, transform=ax.transData)
```

and the default rectangular `geo` spine is hidden. Coastlines, borders,
gridlines, contours, and filled contours are clipped to the curved boundary.

`fig13_lambert_ctrl_verification_12h.png` and
`fig13_lambert_ctrl_verification_24h.png` show Lambert CTRL forecast,
analysis, and error maps. For visual consistency with the report footprint, the
Lambert forecast field is interpolated to the original lat-lon display grid and
smoothly tapered near the native model boundary. This display processing does
not affect the quantitative scores reported in the experiment matrix.

---

## 7. Limitations

1. The BVE is non-divergent and single-level; it cannot represent baroclinic
   development, vertical motion, latent heating, or boundary-layer processes.
2. The fixed boundary condition `ψ=0` is still idealized and does not provide
   time-dependent lateral forcing.
3. Height recovery uses `Z'=fψ/g`; a stricter balanced-height diagnostic would
   solve the linear balance equation.
4. Only one case is tested, so the results are not statistically representative.
5. Figure display interpolation is separate from verification scoring.

---

## References

- Arakawa, A. (1966). Computational design for long-term numerical integration
  of the equations of fluid motion. *J. Comput. Phys.*, 1, 119–143.
- Haltiner, G. J., and R. T. Williams (1980). *Numerical Prediction and
  Dynamic Meteorology*, 2nd ed., Wiley.
- Hersbach, H., et al. (2020). The ERA5 global reanalysis. *Quart. J. Roy.
  Meteorol. Soc.*, 146, 1999–2049.
- Holton, J. R., and G. J. Hakim (2013). *An Introduction to Dynamic
  Meteorology*, 5th ed., Academic Press.
