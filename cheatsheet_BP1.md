# Bonus Project Option 1 — Cheatsheet

## A Simple Forecast Model with Real Atmospheric Data

This cheatsheet collects the equations and numerical recipes you may use for
**Bonus Project Option 1**. The main body describes a **barotropic primitive
equation (BPE) model** integrated for 12–24 h on a 500 hPa isobaric surface.
Brief alternatives — the **barotropic vorticity equation (BVE)** and the
**shallow-water (SW)** model — are included at the end. You may follow this
recipe verbatim, or use it as a reference while implementing your own
simplified model. A simple but well-tested model can receive full credit.

---

## 0. Notation

| Symbol                                 | Meaning                                                       |
| -------------------------------------- | ------------------------------------------------------------- |
| $u, v$                                 | horizontal wind on a map-projection grid (m s⁻¹)              |
| $z$                                    | 500 hPa geopotential height (geopotential metres)             |
| $\phi = gz$                            | geopotential (m² s⁻²)                                         |
| $f$                                    | Coriolis parameter, $f = 2\Omega\sin\varphi$                  |
| $m$                                    | map-projection scale factor (image-distance / earth-distance) |
| $\varphi, \lambda$                     | latitude, longitude                                           |
| $(X, Y)$                               | image-plane coordinates (m); $(i, j)$ are grid indices        |
| $d = \Delta X = \Delta Y$              | grid spacing on the image plane (m)                           |
| $\Delta t$                             | time step (s); $NB$ counts time steps from $t = 0$            |
| $g = 9.80665$ m s⁻²                    | gravity                                                       |
| $\Omega = 7.292\times 10^{-5}$ rad s⁻¹ | earth rotation rate                                           |

Recommended configuration: $d = 100\text{–}300$ km, a Lambert conformal
grid over East Asia (reference point 35°N 115°E, standard parallels 25°N/45°N).
A $60 \times 50$ grid at $d = 100$ km or a $20 \times 17$ grid at $d = 300$ km
are both manageable on a laptop.
Choose $\Delta t$ to satisfy CFL (see §3.6); typical values are
$\Delta t = 300$ s (BPE/SW, $d = 100$ km) or $\Delta t = 600$ s (BVE or $d = 300$ km).
$NB = 1 \ldots 288$ for a 24 h run at $\Delta t = 300$ s;
$NB = 1 \ldots 144$ at $\Delta t = 600$ s.

---

## 1. Workflow

```
   real reanalysis  ──►  initial condition  ──►  numerical model
   (ERA5 / FNL)         (z, u, v on grid)        (12–24 h integration)
                                                         │
                                                         ▼
                                              verifying analysis
                                              + error map / RMSE / ACC
```

You need to produce, at minimum, three figures: **initial field**, **forecast
field**, and **verifying analysis** at the forecast time. An **error map** or
a quantitative score (RMSE or spatial correlation) is strongly encouraged.

---

## 2. Data acquisition (post-2025 East Asia case)

Pick a mid-latitude weather case **not earlier than 2025**. Suggested datasets:

| Dataset          | Native resolution | Access                              |
| ---------------- | ----------------- | ----------------------------------- |
| ERA5 reanalysis  | 0.25° hourly      | Copernicus CDS (`cdsapi` in Python) |
| GDAS / FNL 0.25° | 0.25° 6-hourly    | NCAR RDA [ds083.3](https://rda.ucar.edu/datasets/d083003/) (archive); NCEP NOMADS (last ~10 days only) |
| GFS analysis     | 0.25° 6-hourly    | NCEP NOMADS (last ~10 days); NCAR RDA [ds084.1](https://rda.ucar.edu/datasets/d084001/) (archive) |

Recommended variables on the 500 hPa surface:

- `z` — geopotential height (gpm)
- `u`, `v` — wind components (m s⁻¹)

**Regridding.** ERA5 / FNL are on a regular lat–lon grid. You must interpolate
them onto your **model grid** in the map-projection plane (Lambert conformal
conic or polar-stereographic). Bilinear interpolation is sufficient. Save the
map-projection scale factor $m(i, j)$ and Coriolis parameter $f(i, j)$ as
**PARAM** on the same grid — they are needed every time step.

### Recommended case types

The barotropic models in this cheatsheet work best when the synoptic-scale
flow is **dominated by the rotational component** and **changes little with
height** through the troposphere. The following situations satisfy these
conditions and are recommended:

| Case type | Why it works | Typical season |
| --------- | ------------ | -------------- |
| **Strong East Asian trough + upstream ridge** (e.g. Ural or Siberian blocking high with a deep trough over East Asia) | Highly barotropic, large-scale vorticity advection dominates; BVE captures the trough movement well | Dec–Feb |
| **Rossby wave train propagation** across Eurasia | Long-wave phase speed is well represented by quasi-geostrophic theory; BVE/BPE both work | Nov–Mar |
| **Amplifying long-wave pattern** prior to a major cold surge | The pre-surge build-up is quasi-barotropic; 12–24 h forecast window is within the useful range of a barotropic model | Dec–Feb |

**Situations to avoid:**

- Rapidly deepening extratropical cyclones (baroclinic instability dominates — the model has no temperature field)
- Western Pacific typhoons or tropical systems (strong divergence, latent heat release)
- Active East Asian summer monsoon or Meiyu situations (strong divergent flow, mesoscale convective systems)
- Cases where the 500 hPa flow shows significant tilt with height (strong baroclinic structure)

**Practical tip for case selection.** Browse 500 hPa height analyses from ERA5
or GDAS and look for a date when the mid-latitude belt shows a clear,
smooth long-wave pattern (2–4 large troughs and ridges across the hemisphere).
Avoid dates when the height contours are tightly packed or highly irregular
at small scales — these indicate active baroclinic development that
barotropic models cannot capture.

---

## 3. Barotropic Primitive Equation (BPE) model — main

The BPE atmosphere supports both **slowly-moving long waves** and
**fast-moving inertia–gravity waves**, so the model can represent
quasi-geostrophic evolution **and** geostrophic adjustment. The trade-off is
that the high-frequency waves force a short time step, and the model is
sensitive to **initial imbalance** and **boundary conditions**.

### 3.1 Governing equations in $p$-coordinates

Starting point — momentum and mass continuity on an isobaric surface, with
$z$ standing for the free-surface geopotential height:

$$
\begin{cases}
\dfrac{\partial u}{\partial t} + u\dfrac{\partial u}{\partial x} + v\dfrac{\partial u}{\partial y} - f v + g\dfrac{\partial z}{\partial x} = 0 \\[6pt]
\dfrac{\partial v}{\partial t} + u\dfrac{\partial v}{\partial x} + v\dfrac{\partial v}{\partial y} + f u + g\dfrac{\partial z}{\partial y} = 0 \\[6pt]
\dfrac{\partial z}{\partial t} + u\dfrac{\partial z}{\partial x} + v\dfrac{\partial z}{\partial y} + z\left(\dfrac{\partial u}{\partial x} + \dfrac{\partial v}{\partial y}\right) = 0
\end{cases}
$$

### 3.2 Map-projection coordinates

Rewrite the system in conformal image-plane coordinates $(X, Y, P, t)$ with
map factor $m$:

$$
\begin{cases}
\dfrac{\partial u}{\partial t} = -m\left(u\dfrac{\partial u}{\partial x} + v\dfrac{\partial u}{\partial y}\right) + f^{*} v - m g\dfrac{\partial z}{\partial x} \\[6pt]
\dfrac{\partial v}{\partial t} = -m\left(u\dfrac{\partial v}{\partial x} + v\dfrac{\partial v}{\partial y}\right) - f^{*} u - m g\dfrac{\partial z}{\partial y} \\[6pt]
\dfrac{\partial z}{\partial t} = -m^{2}\left(u\dfrac{\partial}{\partial x}\!\left(\dfrac{z}{m}\right) + v\dfrac{\partial}{\partial y}\!\left(\dfrac{z}{m}\right) + \dfrac{z}{m}\left(\dfrac{\partial u}{\partial x} + \dfrac{\partial v}{\partial y}\right)\right)
\end{cases}
$$

with the **modified Coriolis parameter**

$$
f^{*} = f + u\,\dfrac{\partial m}{\partial y} - v\,\dfrac{\partial m}{\partial x}.
$$

For a Lambert conic projection, $m$ and $f$ are functions of latitude and can
be tabulated once into a `PARAM` array on the model grid.

### 3.3 Finite-difference operators (quadratic-conservative advection)

For any grid field $A$, define the half-step average and centred derivative:

$$
\bar{A}^{x} = \tfrac{1}{2}\left(A_{i+\frac12,j} + A_{i-\frac12,j}\right), \qquad
A_{x} = \tfrac{1}{d}\left(A_{i+\frac12,j} - A_{i-\frac12,j}\right).
$$

On an unstaggered grid the composite operators evaluated at integer indices are

$$
\overline{\bar{A}^{x}}^{x}_{i,j} = \tfrac{1}{2}(A_{i+1,j} + A_{i-1,j}), \qquad
\bar{A}_{x}^{x}\big|_{i,j} = \tfrac{1}{2d}(A_{i+1,j} - A_{i-1,j}).
$$

The **quadratic-conservative** form of $\overline{A B_{x}}$ is

$$
\overline{\bar{A}^{x} B_{x}}^{x}
= \tfrac{1}{2}\!\left(
\dfrac{(A_{i+1,j}+A_{i,j})(B_{i+1,j}-B_{i,j})}{2d}
+ \dfrac{(A_{i,j}+A_{i-1,j})(B_{i,j}-B_{i-1,j})}{2d}
\right).
$$

Analogous formulas hold in the $y$-direction. Applying these operators to
the three governing equations gives the discrete system

$$
\begin{cases}
\dfrac{\partial u_{i,j}}{\partial t} = -m_{i,j}\!\left(\overline{\bar{u}^{x} u_{x}}^{x} + \overline{\bar{v}^{y} u_{y}}^{y} + g\,\bar{z}_{x}^{x}\right) + \widetilde{f^{*}}_{i,j}\, v_{i,j} \equiv E_{i,j} \\[6pt]
\dfrac{\partial v_{i,j}}{\partial t} = -m_{i,j}\!\left(\overline{\bar{u}^{x} v_{x}}^{x} + \overline{\bar{v}^{y} v_{y}}^{y} + g\,\bar{z}_{y}^{y}\right) - \widetilde{f^{*}}_{i,j}\, u_{i,j} \equiv G_{i,j} \\[6pt]
\dfrac{\partial z_{i,j}}{\partial t} = -m_{i,j}^{\,2}\!\left(\overline{\bar{u}^{x}\!\left(\tfrac{z}{m}\right)_{\!x}}^{x} + \overline{\bar{v}^{y}\!\left(\tfrac{z}{m}\right)_{\!y}}^{y} + \dfrac{z_{i,j}}{m_{i,j}}(\bar{u}_{x}^{x} + \bar{v}_{y}^{y})\right) \equiv H_{i,j}
\end{cases}
$$

where the discrete modified Coriolis term is

$$
\widetilde{f^{*}}_{i,j} = f_{i,j} + u_{i,j}\,\bar{m}_{y}^{y} - v_{i,j}\,\bar{m}_{x}^{x}.
$$

> **Implementation tip.** Note $m^{2}\overline{\bar{u}^{x}(z/m)_{x}}^{x} \neq m\,\bar{z}_{x}^{x}$.
> Keep $z/m$ inside the advection operator; do not factor $m$ out.

In matrix form, write the prognostic vector $F = (u, v, z/m)^{T}$ and

$$
\dfrac{\partial F_{i,j}}{\partial t} = A_{i,j}\, F_{i,j}, \qquad
A_{i,j} =
\begin{bmatrix}
L & \widetilde{f^{*}}_{i,j} & -m_{i,j}^{2}\, g\,\overline{(\,\cdot\,)_{x}}^{x} \\
-\widetilde{f^{*}}_{i,j} & L & -m_{i,j}^{2}\, g\,\overline{(\,\cdot\,)_{y}}^{y} \\
-z_{i,j}\,\overline{(\,\cdot\,)_{x}}^{x} & -z_{i,j}\,\overline{(\,\cdot\,)_{y}}^{y} & L
\end{bmatrix}
$$

with $L = -m_{i,j}\!\left(\overline{\bar{u}^{x}(\,\cdot\,)_{x}}^{x} + \overline{\bar{v}^{y}(\,\cdot\,)_{y}}^{y}\right)$.

### 3.4 Geostrophic initial condition (required)

Compute the initial wind from **geostrophic balance** with the interpolated
500 hPa height field. From the reanalysis $z^{0}_{i,j}$ on the model grid:

$$
u^{0}_{i,j} = -\dfrac{m_{i,j}\, g}{f_{i,j}}\,\dfrac{\partial z^{0}_{i,j}}{\partial y}, \qquad
v^{0}_{i,j} = +\dfrac{m_{i,j}\, g}{f_{i,j}}\,\dfrac{\partial z^{0}_{i,j}}{\partial x}.
$$

This is the **recommended default**. Because $u^{0}, v^{0}$ are derived from
$z^{0}$ by construction, the initial state is in geostrophic balance and
high-frequency inertia–gravity waves are not spuriously excited at $t = 0$.
Using raw reanalysis winds instead introduces imbalance that typically
manifests as $2\Delta x$ ringing in $z$ within the first few time steps
(Lynch and Huang 1992).

If you prefer to use the raw reanalysis wind components as the initial
condition, the imbalance can be suppressed by **digital filter
initialization (DFI)**: run the model forward and backward for a short
spin-up window (e.g. 3–6 h each way), then apply a low-pass time filter
(e.g. Lanczos or Dolph-Chebyshev) to the sequence of model states at
$t = 0$. The filtered state retains Rossby-wave scales while removing
inertia–gravity wave oscillations. See Lynch and Huang (1992) for the
method and Lynch (1997) for a computationally efficient implementation.

### 3.5 Lateral boundary condition

Use **fixed lateral boundary conditions** on the perimeter $\beta$ of the
domain:

$$
\left.\dfrac{\partial u_{i,j}}{\partial t}\right|_{\beta} =
\left.\dfrac{\partial v_{i,j}}{\partial t}\right|_{\beta} =
\left.\dfrac{\partial z_{i,j}}{\partial t}\right|_{\beta} = 0.
$$

That is, the four boundary rows/columns retain their initial values for the
entire run. This is the simplest choice and the source of most BPE errors —
the BPE model is **highly sensitive to boundary conditions**. Three
approaches can mitigate boundary-induced errors, in order of implementation
cost: (1) **enlarge the domain** so that boundary reflections do not reach
the region of interest within the integration period; (2) **Davies
relaxation** — nudge the fields in a sponge layer near the boundary back
toward a reference state, so that waves are absorbed rather than reflected;
(3) **time-varying boundary conditions** — update boundary values at each
output interval from interpolated reanalysis fields. The boundary smoothing
applied in §3.8 also partially alleviates the sharp gradient between the
frozen boundary and the evolving interior.

### 3.6 Time integration scheme

The model is started with **Euler-backward** to damp spurious high-frequency
oscillations from initial imbalance, then switched to a **three-step leapfrog
start + centred time difference**.

**Stage 1 — Euler-backward** ($NB = 1 \ldots 6$, i.e. the first hour with
$\Delta t = 600$ s):

$$
\begin{cases}
F^{*\,n+1}_{i,j} = F^{n}_{i,j} + \Delta t\, A^{n}_{i,j}\, F^{n}_{i,j} \\
F^{n+1}_{i,j} = F^{n}_{i,j} + \Delta t\, A^{*\,n+1}_{i,j}\, F^{*\,n+1}_{i,j}
\end{cases}
$$

**Stage 2 — three-step centred-difference start, then centred-difference
leapfrog** ($NB = 7 \ldots 72$, hours 1 through 12):

$$
\begin{cases}
F^{n+\frac12}_{i,j} = F^{n}_{i,j} + \tfrac{1}{2}\Delta t\, A^{n}_{i,j}\, F^{n}_{i,j} \\[4pt]
F^{n+1}_{i,j} = F^{n}_{i,j} + \Delta t\, A^{n+\frac12}_{i,j}\, F^{n+\frac12}_{i,j} \\[4pt]
F^{n+2}_{i,j} = F^{n}_{i,j} + 2\Delta t\, A^{n+1}_{i,j}\, F^{n+1}_{i,j}
\end{cases}
$$

The first two lines run once at the leapfrog initialisation; the third is the
standard centred (leapfrog) step that you iterate.

**Time line.**

```
NB:   0    1                 6   7                            72
      |----| Euler-backward  |---| 3-step start + leapfrog    |
      0 h          1 h                                       12 h
```

After 12 h, repeat the loop (with leapfrog only) to reach 24 h. CFL must hold
on the **fastest gravity wave**:
$\Delta t \le d / (\,|\mathbf{u}|_{\max} + \sqrt{g\bar z}\,)$.
With $\bar z \approx 5500$ gpm, $\sqrt{g\bar z}\approx 230$ m s⁻¹ and
typical $|\mathbf{u}|_{\max} \approx 60$ m s⁻¹:

- $d = 300$ km $\Rightarrow$ $\Delta t \lesssim 1\,030$ s; a **600 s** step is conservative.
- $d = 100$ km $\Rightarrow$ $\Delta t \lesssim 340$ s; a **300 s** step is conservative.

### 3.7 Time smoothing — Robert–Asselin filter

The leapfrog scheme supports a **computational mode** (a spurious $2\Delta t$
oscillation) that grows slowly but never dissipates without explicit damping.
The standard remedy is the **Robert–Asselin filter** applied every step:

$$
\widetilde{F}^{\,n}_{i,j} = F^{n}_{i,j} + \dfrac{\alpha}{2}\,(F^{n+1}_{i,j} - 2F^{n}_{i,j} + F^{n-1}_{i,j}),
$$

with $\alpha \in [0.05, 0.2]$ (typical operational value: $\alpha = 0.1$).
The filtered $\widetilde{F}^n$ replaces $F^n$ before the next leapfrog step.
This is the **recommended approach** for any production-quality leapfrog
integration and should be applied starting from the first leapfrog step
($NB = 7$ onward).

**Simplified alternative.** If you prefer a minimal implementation, apply
the three-point filter only at two consecutive steps after 6 h ($NB = 36, 37$)
with $S = \tfrac{1}{2}$:

$$
\widetilde{F}_{i,j}^{\,n} = (1 - S)\, F^{n}_{i,j} + \dfrac{S}{2}\,(F^{n+1}_{i,j} + F^{n-1}_{i,j}).
$$

With $S = \tfrac{1}{2}$ this completely removes the $2\Delta t$ computational
mode in a single pass. It is less rigorous than the every-step filter but
sufficient for a short (12–24 h) run with geostrophic initial conditions.

### 3.8 Spatial smoothing

Two filters are applied:

**(a) 9-point boundary smoothing** — once per hour during the first 11 h,
applied to the **first interior ring** of grid points (one row/column in
from $\beta$). This reduces the artificial gradient between the frozen
boundary and the evolving interior:

$$
\widetilde{F}_{i,j}^{\,x,y} = F_{i,j}
+ \dfrac{S}{2}(1 - S)\,(F_{i+1,j} + F_{i,j+1} + F_{i-1,j} + F_{i,j-1} - 4F_{i,j})
+ \dfrac{S^{2}}{4}\,(F_{i+1,j+1} + F_{i+1,j-1} + F_{i-1,j+1} + F_{i-1,j-1} - 4F_{i,j}).
$$

**(b) 5-point interior smoothing** — once at $NB = 72$ (12 h), applied to
**all interior points**, to suppress short-wave noise and nonlinear
computational instability:

$$
\widetilde{F}_{i,j}^{\,x,y} = F_{i,j}
+ \dfrac{S}{4}(F_{i+1,j} + F_{i,j+1} + F_{i-1,j} + F_{i,j-1} - 4F_{i,j}).
$$

Use $S = \tfrac{1}{2}$ unless you have a reason to do otherwise.

### 3.9 Algorithm flowchart

```
   Read initial Z500 field  z_ij^0
              │
   Load PARAM:  m_ij,  f_ij  (lat/lon orientation must match!)
              │
   Compute u_ij^0, v_ij^0  from geostrophic balance
              │
   ┌──────────► Compute E_ij, G_ij, H_ij  (tendencies)
   │          │
   │          Euler-backward for 1 h   (NB = 1..6)
   │          │
   │          Leapfrog (centred-diff) for 11 h   (NB = 7..72)
   │          │
   │          NB == 36, 37 ?  → time smoothing (3-point)
   │          │
   │          NB == 72       ?  → spatial smoothing (5-point interior)
   │          │
   │          every hour     ?  → spatial smoothing (9-point boundary ring)
   │          │
   └──────────  NB < 144  ?  → loop, else print & stop
              │
   Done at 24 h.
```

---

## 4. Alternative model A: Barotropic Vorticity Equation (BVE)

A simpler model. Only one prognostic variable, no gravity waves, very stable.

### 4.1 Equations

$$
\dfrac{\partial \zeta}{\partial t} + J(\psi, \zeta + f) = 0, \qquad
\nabla^{2}\psi = \zeta.
$$

- $\zeta$: relative vorticity
- $\psi$: streamfunction, recovered each step by solving the Poisson equation
  $\nabla^{2}\psi = \zeta$ (FFT, multigrid, or SOR work fine).
- $J(\psi, \zeta) = \psi_{x}\zeta_{y} - \psi_{y}\zeta_{x}$ — use the
  **Arakawa Jacobian** for energy + enstrophy conservation.

### 4.2 Numerical scheme

- Time: leapfrog with Robert–Asselin filter ($\alpha \approx 0.1$, every step).
- Space: 2nd-order centred differences on a regular $(x, y)$ grid.
- CFL: $\Delta t \le d / |\mathbf{u}|_{\max}$ — no gravity-wave constraint,
  so $\Delta t$ can be much larger than in the BPE model.
- BC: $\psi = 0$ on $\beta$ (or fixed $\zeta$ on $\beta$).

**Poisson solver.** At every time step you must solve $\nabla^{2}\psi = \zeta$
for $\psi$ with $\psi = 0$ on $\beta$. The recommended approach for Python
users is to assemble the five-point finite-difference Laplacian as a sparse
matrix once before the time loop, then call `scipy.sparse.linalg.spsolve`
(direct solver, exact, fast for grids up to ~$200\times200$) each step:

```python
from scipy.sparse import diags, kron, eye
from scipy.sparse.linalg import spsolve

# Build 1-D second-difference matrix for nx interior points
e = np.ones(nx); D = diags([e, -2*e, e], [-1,0,1], shape=(nx,nx)) / dx**2
L2d = kron(eye(ny), D) + kron(diags([e,-2*e,e],[-1,0,1],shape=(ny,ny))/dy**2, eye(nx))
# Each step: psi_flat = spsolve(L2d, zeta_flat)
```

This direct solver is exact and fast for grids up to $\sim200\times200$.

### 4.3 Initial condition from reanalysis

$\zeta^{0} = \partial_{x} v - \partial_{y} u$ from the reanalysis 500 hPa wind.
Solve $\nabla^{2}\psi^{0} = \zeta^{0}$ on your grid for the initial
streamfunction. To recover the height field afterwards, the
**linear balance equation** is a good diagnostic:

$$
\nabla^{2}(g z) = \nabla \cdot (f \nabla \psi).
$$

---

## 5. Alternative model B: Cartesian Shallow-Water (CartSW) model

A non-trivial dynamical core, simpler than the BPE because there is one
prognostic mass variable instead of three.

> **Note on naming.** The version of the shallow-water equations discussed in
> lecture retains the map-projection scale factor $m$ in the advection and
> pressure-gradient terms. The model here uses plain Cartesian coordinates 
> ($m \equiv 1$), which further simplifies the implementation at the cost 
> of some map-projection accuracy.

### 5.1 Equations

$$
\begin{cases}
\dfrac{\partial u}{\partial t} + u\dfrac{\partial u}{\partial x} + v\dfrac{\partial u}{\partial y} - f v + g\dfrac{\partial h}{\partial x} = 0 \\[4pt]
\dfrac{\partial v}{\partial t} + u\dfrac{\partial v}{\partial x} + v\dfrac{\partial v}{\partial y} + f u + g\dfrac{\partial h}{\partial y} = 0 \\[4pt]
\dfrac{\partial h}{\partial t} + \dfrac{\partial (uh)}{\partial x} + \dfrac{\partial (vh)}{\partial y} = 0
\end{cases}
$$

Interpret $h$ as the local 500 hPa height. Setting $m \equiv 1$ in the BPE
system (§3.2) recovers these equations exactly.

### 5.2 Numerical scheme

- Time: leapfrog with **Robert–Asselin** filter
  $\bar{F}^{n} = F^{n} + \alpha(F^{n+1} - 2F^{n} + F^{n-1})$, $\alpha \approx 0.05$.
- Space: 2nd-order centred differences, optionally on an Arakawa C-grid.
- CFL: $\Delta t \le d / (|\mathbf{u}|_{\max} + \sqrt{g\bar{h}})$.
- BC: **sponge** or **fixed** lateral BC, same options as the BPE model.

### 5.3 Initial condition

Use $h^{0} = z^{0}$ (500 hPa height from reanalysis) and either geostrophic
$u^{0}, v^{0}$ as in §3.4, or the raw reanalysis winds.

---

## 6. Verification

At the verifying time $T \in \{12, 24\}$ h, compare your forecast field
$F^{\text{fcst}}$ with the corresponding reanalysis $F^{\text{anal}}$ on the
**same grid**:

| Metric | Formula                                                                                                                                                                   |
| ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Bias   | $\frac{1}{N}\sum (F^{\text{fcst}} - F^{\text{anal}})$                                                                                                                     |
| RMSE   | $\sqrt{\frac{1}{N}\sum (F^{\text{fcst}} - F^{\text{anal}})^{2}}$                                                                                                          |
| ACC    | $\dfrac{\sum (F^{\text{fcst}} - \bar F^{c})(F^{\text{anal}} - \bar F^{c})}{\sqrt{\sum (F^{\text{fcst}} - \bar F^{c})^{2}}\sqrt{\sum (F^{\text{anal}} - \bar F^{c})^{2}}}$ |

where $\bar F^{c}$ is a climatological reference (e.g. monthly mean Z500 from
ERA5 1991–2020). For a single case study, anomaly correlation against
climatology is more informative than raw RMSE.

A simple visual **error map** $F^{\text{fcst}} - F^{\text{anal}}$ on the
model domain is usually the most revealing diagnostic.

---

## 7. Practical tips

1. **Verify the IC before integrating.** Plot $z^{0}$, $u^{0}$, $v^{0}$ and
   compare against the reanalysis 500 hPa charts on the same date. A
   90° rotation or hemisphere flip from a `lat`/`lon` axis-order bug is the
   single most common error.
2. **Check CFL** before you change $d$ or $\Delta t$. Print $|\mathbf{u}|_{\max}$
   each hour to make sure the wind is not exploding.
3. **Initial imbalance shows up as ringing in $z$**. If the height field
   develops spurious $2\Delta x$ ripples within the first hour, your time step
   is too long, your initial winds are not balanced enough, or your boundary
   gradient is too sharp.
4. **Fixed lateral BCs distort the interior** by hour 12–24. If the
   forecast looks pinned to the boundary, switch to a sponge or relaxation BC.
5. **Reproducibility.** Save: the case date/time, the data source, the
   regridding method, the domain, $\Delta t$, the projection parameters
   (reference point, standard parallels), and the smoothing schedule. Put
   them in your technical note.
6. **Keep it simple.** A clean BVE or shallow-water run that beats persistence
   is more impressive than a buggy BPE run. A correct, reproducible, and
   well-interpreted simple model receives full credit.

---

## 8. References

- Haltiner, G. J., and R. T. Williams (1980), _Numerical Prediction and
  Dynamic Meteorology_, 2nd ed., Wiley — standard reference for BPE/BVE
  discretisation, quadratic-conservative advection, Arakawa Jacobian.
- Robert, A. J. (1966), The integration of a low order spectral form of the
  primitive meteorological equations, _J. Meteorol. Soc. Japan_, **44**, 237–245
  — original Robert time filter.
- Asselin, R. (1972), Frequency filter for time integrations, _Mon. Wea.
  Rev._, **100**, 487–490 — Robert–Asselin parameter $\alpha$.
- Arakawa, A. (1966), Computational design for long-term numerical
  integration of the equations of fluid motion, _J. Comput. Phys._, **1**, 119–143
  — Arakawa Jacobian for energy + enstrophy conservation.
- Lynch, P., and X.-Y. Huang (1992), Initialization of the HIRLAM model
  using a digital filter, _Mon. Wea. Rev._, **120**, 1019–1034 — foundational
  paper showing that uninitialized analyses excite spurious high-frequency
  gravity waves that dominate early forecast error; motivates §3.4.
- Lynch, P. (1997), The Dolph-Chebyshev window: A simple filter for
  efficient initialization, _Mon. Wea. Rev._, **125**, 655–666 — practical
  digital-filter initialization (DFI) implementation.
- ECMWF, _IFS Documentation_, Cy49r1 (2024) — modern operational reference
  for time-integration and smoothing choices.
