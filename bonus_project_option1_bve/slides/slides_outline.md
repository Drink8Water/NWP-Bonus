# Presentation Slides — Lambert-Grid BVE Forecast

---

## Slide 1: Workflow and Model Setup

**Case:** East Asian winter 500 hPa — 2025-12-30 00 UTC
(+12 h / +24 h verification)

**Workflow:**

```
ERA5 Z500, u500, v500 (0.25°)
  → subset 15°N–65°N, 60°E–170°E, coarsen to ~1°
  → interpolate to Lambert conformal grid (54×34, d = 150 km)
  → initial ζ, ψ from geostrophic balance (ψ = 0 on boundaries)
  → integrate finite-area BVE 12–24 h (RK4, Δt = 600 s)
  → verify height anomalies vs ERA5 analysis
```

**Model equation:**

$$\frac{\partial\zeta}{\partial t} = -m^2 J(\psi,\zeta+f) + \nu m^2\nabla^2\zeta - \alpha(\zeta-\zeta_0), \qquad \zeta = m^2\nabla^2\psi$$

**Key configuration:**

| Item | Setting |
| --- | --- |
| Projection | Lambert conformal (25°N / 45°N) |
| Grid | 54 × 34, d = 150 km |
| Poisson solver | DST, `ψ = 0` on all boundaries |
| Jacobian | Arakawa (energy + enstrophy conserving) |
| Time stepping | RK4, Δt = 600 s |
| Verification region | 25°N–55°N, 90°E–145°E |

---

## Slide 2: Results and Interpretation

**Experiment matrix:**

| Experiment | +12 h RMSE | +24 h RMSE | +24 h ACC |
| --- | ---: | ---: | ---: |
| PERSIST | 42.9 m | 70.3 m | 0.953 |
| **CTRL** | 55.4 m | **56.6 m** | **0.963** |
| DIFF | 55.9 m | 57.9 m | 0.960 |
| SPONGE | 50.8 m | 69.0 m | 0.960 |
| DIFF+SPONGE | 52.2 m | 71.9 m | 0.956 |

**Key findings:**

- Case is highly persistent (PERSIST ACC = 0.953 at +24 h).
- At +24 h, CTRL improves RMSE (70.3 → 56.6 m) and ACC (0.953 → 0.963) over
  persistence, demonstrating useful large-scale phase-evolution skill.
- Diffusion and sponge are diagnostic sensitivity runs; neither further
  reduces +24 h RMSE. The Dirichlet boundary configuration is already stable
  for 24 h forecasts.

**Limitations:**

| Limitation | Impact |
| --- | --- |
| Non-divergent, single-level | No cyclone development or baroclinic processes |
| `Z' = fψ/g` diagnostic | First-order geostrophic approximation |
| Fixed `ψ = 0` boundaries | No time-varying lateral forcing |
| One case only | Not statistically representative |
| ~150 km resolution | Synoptic scales only |

**Conclusion:** A clean Lambert finite-area BVE demonstrates modest but
physically meaningful +24 h skill over persistence for a slowly evolving
winter case. Diffusion and sponge are diagnostic confirmations of numerical
stability, not forecast improvements.
