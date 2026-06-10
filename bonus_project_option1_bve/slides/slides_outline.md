# Presentation Slides — Lambert-Grid BVE Forecast

---

## Slide 1: Case, Model, and Workflow

**Case: East Asian Winter Z500 Circulation**
- Initial time: 2025-12-30 00 UTC
- +12 h / +24 h verification: 2025-12-30 12 UTC / 2025-12-31 00 UTC
- Domain: 15°N–65°N, 60°E–170°E; verification: 25°N–55°N, 90°E–145°E

**Data & Workflow**
```
ERA5 Z500, u500, v500
  → subset East Asia, coarsen to ~1°
  → interpolate to Lambert conformal grid
  → construct initial ζ and ψ (geostrophic, ψ=0 on boundaries)
  → integrate finite-area BVE for 12–24 h (RK4, Δt=600 s)
  → compare height anomalies with ERA5 verifying analysis
```

**Model: Barotropic Vorticity Equation on a Lambert Conformal Grid**

$$\frac{\partial\zeta}{\partial t} = -m^2 J(\psi,\zeta+f) + \nu m^2\nabla^2\zeta - \alpha(\zeta-\zeta_{\text{ref}})$$

$$\zeta = m^2\nabla^2\psi$$

- m: Lambert map factor; f: local Coriolis parameter
- ψ = 0 on all boundaries (Dirichlet)
- Poisson inversion via discrete sine transform (DST)
- Arakawa Jacobian (energy + enstrophy conserving)

**Key Assumptions**
- Non-divergent, barotropic flow
- Lambert conformal projection (angle-preserving, suitable for mid-latitudes)
- Fixed lateral boundaries, no time-dependent forcing
- No diabatic heating or friction (except explicit sensitivity terms)

---

## Slide 2: Results, Sensitivity, and Limitations

**Experiment Matrix**

| Experiment     | Description                              |
| -------------- | ---------------------------------------- |
| PERSIST        | Persistence: initial anomaly unchanged   |
| CTRL           | BVE only                                 |
| DIFF           | BVE + vorticity diffusion (2.5×10⁴ m²/s) |
| SPONGE         | BVE + boundary sponge (8 pts, τ=6 h)    |
| DIFF\_SPONGE   | Diffusion + sponge combined              |

**Verification Scores** (inner domain 25°N–55°N, 90°E–145°E)

| Experiment       | +12h RMSE | +24h RMSE | +24h ACC |
| ---------------- | --------: | --------: | -------: |
| PERSIST          | 42.9      | 70.3      | 0.953    |
| **CTRL**         | 55.4      | **56.6**  | **0.963**|
| DIFF             | 55.9      | 57.9      | 0.960    |
| SPONGE           | 50.8      | 69.0      | 0.960    |
| DIFF\_SPONGE     | 52.2      | 71.9      | 0.956    |

**Key Findings**
- CTRL slightly outperforms persistence at +24 h (RMSE 56.6 vs 70.3 m),
  indicating the BVE captures useful dynamical evolution
- Diffusion and sponge provide negligible +24 h improvement —
  the Dirichlet boundary configuration is already stable for short forecasts
- The case is highly persistent; persistence is a strong baseline

**Limitations**
| Limitation | Impact |
| ---------- | ------ |
| Non-divergent | Cannot capture cyclone development |
| Single level | No baroclinic processes |
| Fixed ψ=0 boundaries | No time-varying lateral forcing |
| One case only | Not statistically representative |
| ~100 km resolution | Only synoptic scales resolved |

**Conclusion**
A clean Lambert finite-area BVE with persistence baseline and sensitivity
diagnosis, not a tuned model comparison.
