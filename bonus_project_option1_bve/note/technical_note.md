# Technical Note: 12–24 h 500 hPa Barotropic Vorticity Equation Forecast

**Bonus Project Option 1 — Numerical Weather Prediction Course**

---

## 1. Case and Data

### 1.1 Case Description

The experiment targets a winter circulation pattern over East Asia.  The initial
time is 2025-12-30 00 UTC.  December 2025 featured a well-organised mid-latitude
wave pattern across Eurasia, with a distinct trough over East Asia and an
upstream ridge near the Urals — a configuration that is dominantly barotropic
and well-suited to a single-level non-divergent model.

### 1.2 Data

ERA5 reanalysis data on pressure levels (Hersbach et al., 2020) is used.  The
dataset was obtained from the Copernicus Climate Data Store and stored locally
as a NetCDF file.  Variables used are 500 hPa geopotential (z), u-component of
wind (u), and v-component of wind (v).  The native resolution is 0.25° × 0.25°
on a regular latitude–longitude grid with 6-hourly temporal resolution.

Three time steps are extracted:
- 2025-12-30 00 UTC (initial condition)
- 2025-12-30 12 UTC (+12 h verifying analysis)
- 2025-12-31 00 UTC (+24 h verifying analysis)

### 1.3 Domain and Preprocessing

The data are subset to the East Asian region 15°N–65°N, 60°E–170°E and coarsened
to approximately 1° × 1° resolution via block averaging.  This reduces the grid
size from 721 × 1440 (global) to ~51 × 111 (regional), making the model tractable
on a laptop while retaining synoptic-scale features.  Geopotential (in m² s⁻²) is
converted to geopotential height Z = Φ/g (in metres).

---

## 2. Model Equation and Assumptions

### 2.1 Governing Equation

The barotropic vorticity equation (BVE) on a local beta-plane is used:

$$
\frac{\partial \zeta}{\partial t} + J(\psi, \zeta) + \beta \frac{\partial \psi}{\partial x} = 0,
\quad \nabla^2 \psi = \zeta,
$$

where ζ is relative vorticity (s⁻¹), ψ is streamfunction (m² s⁻¹), β = df/dy,
and J(ψ, ζ) = ∂ψ/∂x · ∂ζ/∂y − ∂ψ/∂y · ∂ζ/∂x is the Jacobian operator.

### 2.2 Beta-Plane Geometry

A local Cartesian coordinate system is defined with reference latitude
φ₀ = 40°N:

$$
x = a \cos\phi_0 (\lambda - \lambda_0), \quad
y = a (\phi - \phi_0), \quad
f_0 = 2\Omega \sin\phi_0, \quad
\beta = \frac{2\Omega \cos\phi_0}{a},
$$

with a = 6.371 × 10⁶ m and Ω = 7.292 × 10⁻⁵ s⁻¹.

### 2.3 Initial Condition

The initial streamfunction is derived from the ERA5 500 hPa geopotential Φ₀
assuming geostrophic balance:

$$
\psi_0 = \frac{\Phi_0 - \overline{\Phi_0}}{f_0}, \quad
\zeta_0 = \nabla^2 \psi_0.
$$

### 2.4 Key Assumptions

1. **Non-divergent flow.** Horizontal divergence is neglected; vertical motion
   and associated stretching of vortex tubes are not represented.
2. **Barotropic structure.** The flow is assumed to have no vertical shear.
3. **Beta-plane approximation.** The spherical geometry is approximated locally;
   errors grow with distance from the reference latitude.
4. **No diabatic or frictional processes.** The model is a dry, inviscid
   dynamical core.
5. **Periodic lateral boundaries.** The Poisson solver assumes doubly periodic
   boundary conditions.

---

## 3. Numerical Method

### 3.1 Spatial Discretisation

Second-order centred finite differences are used on a regular Cartesian grid
with spacing dx ≈ dy ≈ 80–110 km (≈ 1° at mid-latitudes).  The Jacobian is
evaluated using the Arakawa (1966) scheme, which is the arithmetic mean of
three algebraically equivalent forms of J(ψ, ζ).  This discretisation conserves
both domain-integrated kinetic energy and enstrophy when integrated over a
closed or periodic domain, suppressing nonlinear computational instability.

### 3.2 Time Integration

The classical fourth-order Runge–Kutta (RK4) scheme is used with a time step
of Δt = 600 s.  RK4 provides fourth-order accuracy in time with good stability
properties.  The advective CFL number is monitored at the initial time to
ensure numerical stability.

### 3.3 Poisson Solver

At each RK4 substage, the Poisson equation ∇²ψ = ζ is solved using a
two-dimensional Fast Fourier Transform (FFT).  In Fourier space the solution is

$$
\hat{\psi}(k_x, k_y) = -\frac{\hat{\zeta}(k_x, k_y)}{k_x^2 + k_y^2},
\quad (k_x, k_y) \neq (0, 0),
$$

with $\hat{\psi}(0, 0) = 0$, enforcing a zero-mean streamfunction.  The FFT
solver implicitly imposes periodic boundary conditions in both x and y.

---

## 4. Forecast Verification

### 4.1 Verification Metrics

Forecasts at +12 h and +24 h are verified against the corresponding ERA5
analyses on the same grid over the inner verification domain
25°N–55°N, 90°E–145°E:

| Metric                          | Formula |
| ------------------------------- | ------- |
| RMSE (height anomaly)           | √[mean ((Z′_fcst − Z′_anal)²)] |
| Pattern correlation (anomaly)   | Anomaly correlation coefficient |
| Vorticity correlation           | Spatial correlation of relative vorticity |

### 4.2 Height Recovery

Forecast height anomaly is recovered diagnostically from the streamfunction
using the linear balance relation:

$$
Z'_{\text{fcst}} = \frac{f_0}{g} \psi_{\text{fcst}}.
$$

Because the FFT Poisson solver enforces a zero-mean streamfunction, the BVE
forecast does not prognose the absolute domain-mean height.  Verification
therefore focuses on height anomaly, bias, and debiased RMSE rather than a
separate full-field RMSE.

---

## 5. Results

The BVE model was integrated for 24 hours from 2025-12-30 00 UTC.
Verification was performed against ERA5 analyses over the inner domain
25°N–55°N, 90°E–145°E.

### 5.1 Verification Scores

| Metric                          | +12 h    | +24 h    |
| ------------------------------- | -------- | -------- |
| Height anomaly RMSE (m)         | 80.8     | 547.0    |
| Height debiased RMSE (m)        | 53.8     | 227.3    |
| Height anomaly correlation      | 0.973    | 0.515    |
| Vorticity correlation           | 0.332    | 0.225    |
| Height bias (m)                 | +60.3    | +497.5   |

### 5.2 Discussion

At +12 h, the anomaly correlation coefficient (ACC) of 0.973 indicates
excellent skill in capturing the large-scale height pattern.  The RMSE
of 81 m reflects amplitude errors, particularly near the domain
boundaries where the periodic FFT Poisson solver creates artificial
interactions between the northern and southern edges.

By +24 h, the ACC drops to 0.515 with a large positive bias (+498 m)
and RMSE of 547 m.  The centred (bias-corrected) RMSE is approximately
227 m, which is more indicative of the pattern error.  The systematic
drift toward higher heights reflects the accumulation of boundary
errors and the missing divergent damping that would constrain
geopotential tendencies in a more complete model.

The model captures the first-order Rossby-wave propagation with high
fidelity at 12 h but loses predictive skill by 24 h, consistent with
expectations for an unforced single-level barotropic model with
periodic lateral boundaries on a beta-plane.

The key results are:

1. Initial 500 hPa geopotential height field showing the synoptic pattern
   at 2025-12-30 00 UTC.
2. +12 h and +24 h forecast height anomaly, verifying analysis anomaly,
   and their difference (error map).
3. Initial relative vorticity field with physically reasonable interior
   values (mean |ζ| ≈ 3.9 × 10⁻⁵ s⁻¹).
4. RMSE and pattern correlation at +12 h and +24 h.

The BVE model captures the large-scale trough and ridge movement reasonably
well over the first 12 hours, with forecast skill degrading by 24 hours as
the lack of divergent and baroclinic processes becomes significant.
Boundary contamination from the periodic Poisson solver is most visible
near the northern and southern edges of the domain.

### 5.3 Lambert Finite-Area Upgrade

To test whether the large +24 h degradation is primarily a boundary-geometry
problem, a Lambert conformal finite-area version was added. ERA5 fields are
interpolated to a regular Lambert grid with spacing 150 km. The model uses
spatially varying map factor $m(i,j)$ and Coriolis parameter $f(i,j)$:

$$
\frac{\partial \zeta}{\partial t}
= -m^2 J(\psi,\zeta+f) + \nu m^2\nabla^2\zeta - \alpha(\zeta-\zeta_0),
\qquad
\zeta = m^2\nabla^2\psi.
$$

The Poisson inversion is solved with zero streamfunction on all boundaries,
using a discrete sine-transform solver rather than a doubly periodic FFT.

| Experiment       | +12 h RMSE | +24 h RMSE | +24 h ACC |
| ---------------- | ---------- | ---------- | --------- |
| PERSIST_LCC      | 42.9 m     | 70.3 m     | 0.953     |
| CTRL_LCC         | 55.4 m     | 56.6 m     | 0.963     |
| DIFF_LCC         | 55.9 m     | 57.9 m     | 0.960     |
| SPONGE_LCC       | 50.8 m     | 69.0 m     | 0.960     |
| DIFF_SPONGE_LCC  | 52.2 m     | 71.9 m     | 0.956     |

The upgrade reduces CTRL +24 h RMSE from 547.0 m to 56.6 m and removes the
large positive bias. In this new configuration, diffusion and sponge layers
are no longer the dominant source of improvement; their role is limited to
small-scale noise control and boundary relaxation.

---

## 6. Sensitivity Experiments

### 6.1 Motivation

The CTRL (BVE-only) forecast exhibits a strong positive bias (+498 m) and
reduced pattern correlation (ACC = 0.52) at +24 h.  To diagnose the cause
of this degradation, a matrix of sensitivity experiments was conducted:

| Experiment    | Description                                              |
| ------------- | -------------------------------------------------------- |
| PERSIST       | Persistence: forecast anomaly = initial anomaly at all times |
| CTRL          | BVE: ∂ζ/∂t + J(ψ, ζ) + β ψₓ = 0                        |
| DIFF          | BVE + Laplacian vorticity diffusion (ν = 2.5 × 10⁴ m² s⁻¹) |
| SPONGE        | BVE + boundary sponge relaxation (width = 8 pts, τ = 6 h) |
| DIFF\_SPONGE  | BVE + diffusion + sponge combined                        |

### 6.2 Results

All scores computed over the inner verification domain 25°N–55°N, 90°E–145°E.

| Experiment   | Lead | RMSE (m) | Bias (m) | Deb. RMSE (m) | ACC   |
| ------------ | ---- | -------- | -------- | ------------- | ----- |
| PERSIST      | +12h | 35.7     | −1.2     | 35.6          | 0.986 |
| PERSIST      | +24h | 65.3     | −0.4     | 65.3          | 0.953 |
| CTRL         | +12h | 80.8     | +60.3    | 53.8          | 0.974 |
| CTRL         | +24h | 547.0    | +497.5   | 227.3         | 0.515 |
| DIFF         | +12h | 54.2     | +21.0    | 49.9          | 0.973 |
| DIFF         | +24h | 193.7    | +173.1   | 86.9          | 0.918 |
| SPONGE       | +12h | 67.4     | +44.7    | 50.5          | 0.977 |
| SPONGE       | +24h | 193.2    | +160.7   | 107.3         | 0.857 |
| DIFF\_SPONGE | +12h | 53.1     | +21.2    | 48.7          | 0.978 |
| DIFF\_SPONGE | +24h | 103.5    | +54.9    | 87.8          | 0.921 |

### 6.3 Interpretation

**Persistence as a strong baseline.**  The December 2025 East Asian 500 hPa
pattern was highly persistent: the initial height anomaly alone achieves an
ACC of 0.95 and RMSE of only 65 m at +24 h.  No integrated model beats this
baseline — a striking but physically meaningful result for a winter blocking-like
pattern with slow evolution.  This demonstrates that persistence is a non-trivial
benchmark that barotropic models must be evaluated against.

**Diffusion (DIFF) dramatically improves CTRL.**  At +24 h, adding weak
Laplacian vorticity diffusion (ν = 2.5 × 10⁴ m² s⁻¹) reduces the RMSE from
547 m to 194 m (65 % reduction), the bias from +498 m to +173 m, and raises
the ACC from 0.52 to 0.92.  This indicates that a significant fraction of
the CTRL error growth is driven by grid-scale vorticity noise — likely
generated at the boundaries by the periodic FFT Poisson solver and then
amplified by the nonlinear Jacobian.  The diffusion term selectively damps
small scales while preserving the synoptic-scale Rossby wave signal.

**Sponge relaxation (SPONGE) provides moderate improvement.**  The boundary
sponge reduces +24 h RMSE to 193 m and ACC to 0.86.  It is less effective
than diffusion at preserving the spatial pattern, likely because the sponge
only constrains the boundary region and does not damp internally-generated
noise.  However, its bias reduction (+161 m vs CTRL's +498 m) confirms that
boundary effects are a real contributor to the domain-mean drift.

**Combined (DIFF\_SPONGE) gives the best integrated-model forecast.**
The combined experiment achieves +24 h RMSE = 104 m, ACC = 0.92, and bias
= +55 m — the best among all integrated models.  The diffusion handles
interior noise and the sponge handles boundary contamination, with
complementary effects.

**Debiased RMSE isolates pattern error.**  For CTRL, 78 % of the +24 h RMSE
(497.5² / 547.0² ≈ 0.83) comes from the domain-mean bias.  The debiased
RMSE of 227 m is the true spatial-pattern error.  DIFF reduces this to
87 m and DIFF\_SPONGE to 88 m, confirming that the spatial pattern is
substantially improved by the added damping.

**Vorticity diagnostics.**  PERSIST achieves the highest vorticity
correlation (0.37 at +24 h), while all integrated models have lower values
(0.10–0.23).  This suggests that the BVE model — even with diffusion —
introduces small-scale vorticity features not present in the initial
condition, reducing point-wise vorticity correlation even when the
large-scale height pattern is well-captured (ACC > 0.9 for DIFF).

### 6.4 Scientific Conclusions

1. **The +24 h CTRL degradation is caused by both grid-scale noise and
   boundary effects**, with boundary-generated noise propagating into the
   interior and being amplified by nonlinear advection.
2. **Laplacian diffusion is an effective and simple remedy** for this
   class of periodic-boundary BVE model, reducing RMSE by ~65 % at +24 h.
3. **Persistence is a strong baseline for a slowly-evolving winter pattern.**
   The BVE model should be expected to beat persistence only for cases with
   significant Rossby-wave propagation; the 2025-12-30 case was inherently
   persistent.
4. **These experiments are diagnostic, not tuning.**  The goal is to
   understand error sources, not to optimise scores.  The same ν and sponge
   parameters would need re-evaluation for a different case.

---

## 7. Limitations

1. **Model physics.** The BVE is a single prognostic equation; it cannot
   represent divergence, vertical motion, baroclinic instability, latent
   heat release, or boundary-layer processes.
2. **Boundary conditions.** The periodic FFT Poisson solver is inappropriate
   for a limited-area domain.  Boundary-induced errors increase with
   integration time.
3. **Resolution.** The 1° grid resolves only synoptic scales.  Forecast
   error at smaller scales grows rapidly.
4. **Initialisation.** The geostrophic initial condition neglects ageostrophic
   components of the flow.
5. **Forecast horizon.** Skill beyond 24–48 h is not expected from a
   single-level, unforced model.
6. **Single case.** Results from one case study are not statistically
   representative of model performance.

Despite these limitations, the experiment provides a clear demonstration
of the behaviour of a simplified dynamical core with real initial data —
a useful pedagogical exercise in numerical weather prediction.

---

## 8. Individual Contributions

*(To be completed by the student.)*

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
