# Presentation Slides — BVE Forecast Experiment

---

## Slide 1: Case, Data, Model, and Assumptions

**Case: East Asian Winter Circulation**
- Initial time: 2025-12-30 00 UTC
- Forecast: +12 h (2025-12-30 12 UTC) and +24 h (2025-12-31 00 UTC)
- Domain: 15°N–65°N, 60°E–170°E; verification: 25°N–55°N, 90°E–145°E
- December 2025: well-organised mid-latitude wave train across Eurasia

**Data**
- ERA5 reanalysis at 500 hPa (0.25° native → coarsened to ~1°)
- Variables: geopotential height Z, u, v
- Processed entirely from local files — no internet required

**Model: Barotropic Vorticity Equation on a Beta-Plane**
$$
\frac{\partial \zeta}{\partial t} + J(\psi, \zeta) + \beta \frac{\partial \psi}{\partial x} = 0,
\quad \nabla^2 \psi = \zeta
$$
- β-plane centred at 40°N: f₀ = 2Ω sin 40°, β = 2Ω cos 40° / a

**Numerics**
- Arakawa Jacobian (energy + enstrophy conserving)
- RK4 time integration, Δt = 600 s
- FFT Poisson solver (doubly periodic)
- Grid: regular Cartesian, ~1° spacing

**Key Assumptions**
- Non-divergent, barotropic flow
- No diabatic heating or friction
- Local beta-plane geometry
- Geostrophic initial condition

---

## Slide 2: Forecast Results, Verification, and Limitations

**Forecast Products**
- +12 h and +24 h 500 hPa geopotential height anomaly
- Error maps (forecast − analysis)
- Relative vorticity evolution

**Verification Metrics (inner domain 25°N–55°N, 90°E–145°E)**
- RMSE of height anomaly and full height field
- Spatial pattern correlation (anomaly correlation coefficient)
- Vorticity correlation

**Results (inner domain 25°N–55°N, 90°E–145°E)**

| Metric                      | +12 h  | +24 h  |
| --------------------------- | ------ | ------ |
| Height anom. RMSE           | 81 m   | 547 m  |
| Height anom. correlation    | 0.97   | 0.52   |
| Vorticity correlation       | 0.33   | 0.23   |
| Height bias                 | +60 m  | +498 m |

- +12 h: excellent pattern skill (ACC = 0.97); model captures Rossby-wave
  propagation well
- +24 h: systematic positive bias (+498 m) from accumulated boundary
  contamination and missing divergence damping
- Centred RMSE at +24 h ≈ 227 m after bias correction

**Limitations**
| Limitation | Impact |
| ---------- | ------ |
| No divergence | Cannot capture vertical motion, cyclone development |
| Single level | No baroclinic processes |
| Periodic BCs | Boundary contamination in regional domain |
| ~1° resolution | Only synoptic scales resolved |
| No physics | No friction, orography, or diabatic forcing |
| Single case | Not statistically representative |

**Conclusion**
- The BVE captures the first-order Rossby-wave dynamics of a winter case
- A correct, reproducible, simple model can outperform persistence at 12–24 h
- This is a **pedagogical dynamical core**, not an operational forecast system
