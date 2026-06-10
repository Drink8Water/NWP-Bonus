# BVE Forecast Experiment — Bonus Project Option 1

## 项目目标

用正压涡度方程（BVE）在局地 β 平面上做一个 12–24 h 的 500 hPa 预报实验。
个例是 2025 年 12 月底的一次东亚冬季环流过程，初始场和验证场都来自 ERA5 再分析资料。

## 个例

| 项目                  | 值                           |
| --------------------- | ---------------------------- |
| 初始时间              | 2025-12-30 00 UTC            |
| +12 h 验证            | 2025-12-30 12 UTC            |
| +24 h 验证            | 2025-12-31 00 UTC            |
| 等压面                | 500 hPa                      |
| 全区域                | 15°N–65°N, 60°E–170°E        |
| 验证子区域            | 25°N–55°N, 90°E–145°E        |
| β 平面参考纬度        | 40°N                         |

2025 年 12 月底东亚—西伯利亚中纬度有明显的长波槽脊结构，500 hPa 环流以旋转风为主，我们采用正压模型进行预测。

## 数据

需要 December 2025 的 ERA5 再分析资料，包含 500 hPa 的位势高度（或 geopotential）、
u 风和 v 风，放在规则经纬网格上。把 NetCDF 文件放到：

```
data/raw/
```

预处理脚本会自动搜索 `data/raw/` 及其父目录下的 `.nc` 文件。
数据下载后不需要联网。

## 模型方程

局地 β 平面上的正压涡度方程：

$$
\frac{\partial \zeta}{\partial t} + J(\psi, \zeta) + \beta \frac{\partial \psi}{\partial x} = 0,
\qquad
\zeta = \nabla^2 \psi
$$

| 符号   | 含义                      | 单位     |
| ------ | ------------------------- | -------- |
| ζ      | 相对涡度                  | s⁻¹      |
| ψ      | 流函数                    | m² s⁻¹   |
| β      | f 的经向梯度 (df/dy)      | m⁻¹ s⁻¹  |
| J(·,·) | Jacobian 算子             | —        |

β 平面中心取 φ₀ = 40°N：

- `f₀ = 2Ω sin(φ₀)`
- `β = 2Ω cos(φ₀) / a`
- `x = a cos(φ₀) (λ − λ₀)`, `y = a (φ − φ₀)`

常数：`a = 6.371×10⁶ m`, `Ω = 7.292×10⁻⁵ s⁻¹`, `g = 9.80665 m s⁻²`。

初始场由 ERA5 的 500 hPa geopotential Φ₀ = g Z₀ 导出：

$$
\psi_0 = \frac{\Phi_0 - \overline{\Phi_0}}{f_0}, \qquad
\zeta_0 = \nabla^2 \psi_0
$$

## 数值方法

- **空间离散**：规则 Cartesian 网格，二阶中心差分。Jacobian 用 Arakawa (1966) 格式，
  同时守恒动能和拟能。
- **时间积分**：经典四阶 Runge–Kutta (RK4)。
- **Poisson 求解**：FFT 谱方法，双周期边界条件。零波数分量设为零（零平均流函数）。
- **网格**：从原始 0.25° 粗化到约 1°（中纬度约 80–110 km），50 × 110 个格点。
- **时间步长**：Δt = 600 s，满足平流 CFL 条件。

## 安装依赖

```bash
conda env create -f environment.yml
conda activate bve-nwp
```

或者用 pip：

```bash
pip install numpy scipy xarray netCDF4 matplotlib cartopy
```

## 怎么跑

数据放好之后，按顺序执行：

```bash
# 1. 预处理本地数据
python scripts/01_preprocess_local_data.py

# 2. 运行 BVE 预报
python scripts/02_run_bve_forecast.py

# 3. 计算评分
python scripts/03_compute_scores.py

# 4. 画图
python scripts/04_make_figures.py

# 5. 灵敏度实验（可选）
python scripts/05_run_sensitivity_experiments.py

# 6. 灵敏度实验画图（可选）
python scripts/06_sensitivity_figures.py

# 7. 插值到 Lambert 模式网格
python scripts/07_prepare_lambert_data.py

# 8. 运行 Lambert / 固定边界实验矩阵
python scripts/08_run_lambert_experiments.py

# 9. 对比旧 beta-plane 与新 Lambert 结果
python scripts/09_compare_lambert_impact.py
```

## 输出文件

### 预处理和预报

| 文件                                          | 内容                             |
| --------------------------------------------- | -------------------------------- |
| `data/processed/initial_20251230_00.npz`      | 初始场 (Z, U, V)                 |
| `data/processed/analysis_20251230_12.npz`     | +12 h 分析场                     |
| `data/processed/analysis_20251231_00.npz`     | +24 h 分析场                     |
| `data/processed/grid_info.npz`                | 网格参数                         |
| `outputs/forecast_12_0h.npz`                  | CTRL 预报 ψ, ζ（+12 h）          |
| `outputs/forecast_24_0h.npz`                  | CTRL 预报 ψ, ζ（+24 h）          |
| `outputs/verification_+12h.npz`               | CTRL 验证场（+12 h）             |
| `outputs/verification_+24h.npz`               | CTRL 验证场（+24 h）             |
| `outputs/scores.json`                         | CTRL 的 RMSE、ACC、bias          |
| `outputs/scores_lambert.json`                 | Lambert 实验矩阵评分             |
| `outputs/lambert_impact_summary.csv`          | beta-plane vs Lambert 对比表      |

### 图片

| 文件                                            | 内容                                  |
| ----------------------------------------------- | ------------------------------------- |
| `figures/fig01_initial_z500.png`                | 初始 500 hPa 位势高度                 |
| `figures/fig02_forecast_verification_12h.png`   | +12 h 预报-分析-误差 三面板           |
| `figures/fig03_forecast_verification_24h.png`   | +24 h 预报-分析-误差 三面板           |
| `figures/fig04_initial_vorticity.png`           | 初始相对涡度                          |
| `figures/fig05_scores_summary.png`              | CTRL 评分汇总                         |
| `figures/fig06_sensitivity_scores.png`          | 灵敏度实验评分对比                    |
| `figures/fig07_error_comparison_24h.png`        | +24 h 误差四面板对比                  |
| `figures/fig08_improvement_24h.png`             | +24 h 相对 CTRL 改进量                |
| `figures/fig09_vorticity_ctrl_vs_diff_24h.png`  | +24 h 涡度 CTRL vs DIFF               |
| `figures/fig10_diff_sponge_verification_12h.png`| DIFF_SPONGE +12 h 验证                |
| `figures/fig11_diff_sponge_verification_24h.png`| DIFF_SPONGE +24 h 验证                |
| `figures/fig12_diff_vorticity_24h.png`          | DIFF +24 h 涡度                       |

数值积分使用局地 beta 平面 Cartesian 网格；Lambert Conformal 投影只用于地图绘制
（标准纬线 25°N/47°N，中心 105°E/35°N）。填色区域是投影空间的规则矩形，
不是经纬度矩形投影后的弯曲扇形。

## CTRL 预报结果

验证区域 25°N–55°N, 90°E–145°E：

| 指标               | +12 h | +24 h  |
| ------------------ | ----- | ------ |
| 高度异常 RMSE      | 81 m  | 547 m  |
| 高度异常 ACC       | 0.97  | 0.52   |
| 涡度相关           | 0.33  | 0.23   |
| Bias               | +60 m | +498 m |

+12 h 的空间形态几乎完全正确（ACC 0.97）。+24 h 出现了很强的正偏差（接近 +500 m），
主要来自 FFT 周期性边界产生的虚假涡度被非线性平流放大。去掉系统性偏差后，
+24 h 的去偏 RMSE 约 227 m——这才是真正的空间形态误差。

## 灵敏度实验

为了诊断 +24 h 退化的原因，跑了五个实验：

| 实验         | 说明                                              |
| ------------ | ------------------------------------------------- |
| PERSIST      | 持续性：假设初始异常完全不变                       |
| CTRL         | 原始 BVE                                          |
| DIFF         | BVE + 涡度扩散（ν = 2.5×10⁴ m²/s）                |
| SPONGE       | BVE + 边界海绵松弛（宽度 8 格点, τ = 6 h）         |
| DIFF_SPONGE  | 扩散 + 海绵联合                                   |

### 灵敏度实验结果

| 实验         | 时效  | RMSE (m) | Bias (m) | 去偏 RMSE (m) | ACC   |
| ------------ | ----- | -------- | -------- | ------------- | ----- |
| PERSIST      | +12h  | 36       | −1       | 36            | 0.99  |
| PERSIST      | +24h  | 65       | 0        | 65            | 0.95  |
| CTRL         | +12h  | 81       | +60      | 54            | 0.97  |
| CTRL         | +24h  | 547      | +498     | 227           | 0.52  |
| DIFF         | +12h  | 54       | +21      | 50            | 0.97  |
| DIFF         | +24h  | 194      | +173     | 87            | 0.92  |
| SPONGE       | +12h  | 67       | +45      | 51            | 0.98  |
| SPONGE       | +24h  | 193      | +161     | 107           | 0.86  |
| DIFF_SPONGE  | +12h  | 53       | +21      | 49            | 0.98  |
| DIFF_SPONGE  | +24h  | 104      | +55      | 88            | 0.92  |

### 几点发现

1. **持续性模型是个很强的 baseline。** 2025-12-30 这次过程的 500 hPa 高度场非常稳定，PERSIST 的 +24 h ACC 高达 0.95、RMSE 仅 65 m。所有积分的 BVE 模型在 RMSE 上都低于此模型。对于缓慢演变的冬季阻塞型环流，正压模型可能引入比它消除还多的误差。

2. **扩散大幅改善了 CTRL。** +24 h RMSE 从 547 m 降到 194 m（降了 65%），ACC 从 0.52 恢复到 0.92。原因是 FFT 周期性边界产生的格点尺度涡度噪声被非线性 Jacobian 放大，弱扩散选择性阻尼了小尺度噪声，保留了天气尺度的 Rossby 波信号。

3. **海绵单独用效果中等。** +24 h RMSE 193 m，ACC 0.86。海绵只在边界附近把涡度往回拽，没法处理内部的格点噪声。

4. **联合（DIFF_SPONGE）是积分模型里最好的。** +24 h RMSE 104 m，ACC 0.92，bias 仅 +55 m。扩散处理内部噪声，海绵处理边界污染，互补。

5. **CTRL 的 RMSE 大约 83% 来自系统性偏差。** 去偏后的 RMSE 只有 227 m（对比原始 547 m）。DIFF 和 DIFF_SPONGE 把去偏 RMSE 降到了 87–88 m，说明扩散确实改善了空间形态，不只是压制漂移。

灵敏度实验的输出文件：

| 文件                                        | 内容                     |
| ------------------------------------------- | ------------------------ |
| `outputs/scores_experiment_matrix.json`     | 全部实验的评分（JSON）    |
| `outputs/scores_experiment_matrix.csv`      | 全部实验的评分（CSV）     |
| `outputs/{EXP}_forecast_{LEAD}h.npz`        | 各实验的预报场            |
| `outputs/{EXP}_verification_{LEAD}h.npz`    | 各实验的验证场            |

## Lambert 有限区域动力核心

新增 Lambert 版本使用投影平面规则网格、空间变化的 `m(i,j)` 与 `f(i,j)`，
并把 Poisson 求解从双周期 FFT 改为固定边界 `ψ = 0` 的正弦变换解法。预报方程写成：

$$
\frac{\partial \zeta}{\partial t}
= -m^2 J(\psi,\zeta+f) + \nu m^2\nabla^2\zeta - \alpha(\zeta-\zeta_0),
\qquad
\zeta = m^2\nabla^2\psi.
$$

Lambert 网格为 54 × 34，`d = 150 km`，覆盖约 56°E–174°E, 12°N–62°N。
验证区仍为 25°N–55°N, 90°E–145°E。

| 实验             | +12 h RMSE | +24 h RMSE | +24 h ACC | 说明 |
| ---------------- | ---------- | ---------- | --------- | ---- |
| PERSIST_LCC      | 42.9 m     | 70.3 m     | 0.953     | Lambert 网格上的持续性基准 |
| CTRL_LCC         | 55.4 m     | 56.6 m     | 0.963     | Lambert + 固定边界主实验 |
| DIFF_LCC         | 55.9 m     | 57.9 m     | 0.960     | 扩散在新核心中不再明显改善 RMSE |
| SPONGE_LCC       | 50.8 m     | 69.0 m     | 0.960     | 降低 bias，但 +24 h RMSE 略差 |
| DIFF_SPONGE_LCC  | 52.2 m     | 71.9 m     | 0.956     | 联合方案不优于 CTRL_LCC |

最主要的改进来自边界和几何处理：旧 beta-plane 周期边界 CTRL 的 +24 h RMSE 为
547.0 m，新 Lambert/固定边界 CTRL 降到 56.6 m，下降约 89.7%。这说明旧 CTRL 的
大误差主要不是缺少 diffusion/sponge，而是周期边界和简化几何导致的系统性漂移。

## 已知局限

1. **无辐散。** BVE 假设大气无辐散，槽前辐散和脊前辐合都不在方程里。
2. **没有斜压过程。** 单层模型没法表示温度平流、垂直耦合和斜压不稳定。
3. **旧 beta-plane 路径仍是周期性边界。** FFT Poisson 求解器强制双周期边界，对区域模式不现实。
   南北边界从 65°N"连"到 16°N，产生虚假梯度并向内传播。
4. **Lambert 高度恢复仍是诊断近似。** 当前用 $Z'=f\psi/g$ 恢复高度距平；更严格版本应解线性平衡方程。
5. **粗分辨率。** ~1° 网格只能分辨天气尺度（>500 km），中尺度细节全丢了。
6. **没有物理参数化。** 无摩擦、无地形、无非绝热加热——纯粹干动力学核心。
7. **预报时效短。** 没有外强迫、边界误差累积，24–48 h 之后基本没法看。
8. **只用了一个个例。** 结果不能代表模型在别的环流背景下的表现。
9. **持续性碾压了所有积分模型。** 对于这次高度稳定的冬季环流个例，初始场本身就包含了
   24 h 内大部分有用信息。对于 Rossby 波传播更快的过程，BVE 模型可能会相对更有价值。

## 参考文献

- Arakawa, A. (1966). Computational design for long-term numerical integration of the equations of fluid motion. *J. Comput. Phys.*, 1, 119–143.
- Haltiner, G. J., and R. T. Williams (1980). *Numerical Prediction and Dynamic Meteorology*, 2nd ed., Wiley.
- Holton, J. R., and G. J. Hakim (2013). *An Introduction to Dynamic Meteorology*, 5th ed., Academic Press.
- Hersbach, H., et al. (2020). The ERA5 global reanalysis. *Quart. J. Roy. Meteorol. Soc.*, 146, 1999–2049.
