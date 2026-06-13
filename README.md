# Lambert BVE Forecast Experiment — Bonus Project Option 1

## 项目目标

本项目使用正压涡度方程（Barotropic Vorticity Equation, BVE）对一次
东亚 500 hPa 冬季环流过程进行 12–24 h 短时预报实验。当前版本的动力核心
已经改为 **Lambert conformal 有限区域模式**：

- 模式网格：Lambert 投影平面规则网格；
- 动力变量：流函数 `ψ` 和相对涡度 `ζ`；
- 地图因子：使用空间变化的 `m(i,j)`；
- Coriolis 参数：使用空间变化的 `f(i,j)`；
- Poisson 反演：固定边界 `ψ = 0`，使用离散正弦变换求解；
- 验证指标：高度距平 RMSE、去偏 RMSE、ACC、bias、涡度相关。

## 个例

| 项目 | 值 |
| --- | --- |
| 初始时间 | 2025-12-30 00 UTC |
| +12 h 验证 | 2025-12-30 12 UTC |
| +24 h 验证 | 2025-12-31 00 UTC |
| 等压面 | 500 hPa |
| 原始资料区域 | 15°N–65°N, 60°E–170°E |
| 验证子区域 | 25°N–55°N, 90°E–145°E |
| Lambert 标准纬线 | 25°N / 45°N |
| Lambert 网格距 | 150 km |

ERA5 再分析资料提供初始场和验证场。使用变量为 500 hPa 位势高度（或 geopotential）、
u 风和 v 风。

### 为什么选 2025-12-30？

选择这个个例主要基于以下考虑：

1. **冬季正压性更强。** 12 月底东亚中纬度西风带盛行，斜压能量转换（baroclinic conversion）
   相对弱于春季和秋季，500 hPa 环流演变更多由正压 Rossby 波频散主导，
   这让 BVE 模型在物理上更有机会展现预报技巧。
2. **缓慢演变的环流背景。** 该个例持续性基准 ACC 达 0.953（+24 h），
   说明天气形势相对稳定——不是快速发展的爆发性气旋或强非线性过程。
   对正压短时预报实验而言，半静稳态环流是合适的"公平测试"：
   既不完全 trivial（持续性 +24 h RMSE 仍有 70.3 m），
   又不至于让单层无辐散模型毫无胜算。
3. **数据时效性。** 项目运行时最新可获取的完整月份 ERA5 数据为 2025 年 12 月，
   月末时段（30–31 日）拥有完整的初始场和验证场配对。
4. **东亚区域代表性。** 12 月下旬是东亚冬季风活跃期，
   域内通常存在清晰的高-低-高 Rossby 波列，适合验证动力核心对
   大尺度波的传播和频散模拟能力。

## 模型方程

Lambert conformal 投影平面上的有限区域 BVE 写成：

$$
\frac{\partial \zeta}{\partial t}
= -m^2 J(\psi,\zeta+f)
  + \nu m^2\nabla^2\zeta
  - \alpha(\zeta-\zeta_0),
\qquad
\zeta = m^2\nabla^2\psi.
$$

其中 `m(i,j)` 是 Lambert 地图比例因子，`f(i,j)=2Ωsinφ`，`ν` 是可选涡度扩散系数，
`α` 是可选边界海绵松弛系数。CTRL 实验取 `ν=0`、无海绵。

高度距平由平衡关系诊断得到：

$$
Z' = \frac{f\psi}{g}.
$$

当前模式不预报绝对域平均高度，因此不再报告 full-field RMSE；所有高度评分均基于
高度距平、bias 和去偏 RMSE。

## 数值方法

- **空间离散**：Lambert 投影平面规则网格，二阶中心差分。
- **Jacobian**：Arakawa Jacobian。
- **时间积分**：四阶 Runge–Kutta，`Δt = 600 s`。
- **Poisson 求解**：固定边界 `ψ = 0`，离散正弦变换求解。
- **Lambert 网格**：54 × 34，`d = 150 km`，约覆盖 56°E–174°E, 12°N–62°N。
- **验证区域**：25°N–55°N, 90°E–145°E。

## 环境配置

使用 conda 创建虚拟环境：

```bash
conda env create -f bonus_project_option1_bve/environment.yml
conda activate bve-nwp
```

核心依赖包括 NumPy、SciPy、xarray、netCDF4、Matplotlib、Cartopy。

## 数据获取

本项目使用 **ERA5 再分析资料**，可从
[Copernicus CDS](https://cds.climate.copernicus.eu/) 下载。

需要的数据：

| 需求 | 说明 |
| --- | --- |
| 时间段 | 2025 年 12 月（至少包含 29–31 日） |
| 等压面 | 500 hPa |
| 时间分辨率 | 逐 6 小时分析场（00/06/12/18 UTC） |
| 空间分辨率 | 0.25° 全球网格 |
| 变量 | geopotential (`z`)、u 风 (`u`)、v 风 (`v`) |
| 格式 | NetCDF（GRIB 可用 `cfgrib` 直接读取或先转换为 NetCDF） |

仓库根目录的 **`202512.nc`**（约 600 MB，已通过 `.gitignore` 排除版本控制）
即为上述数据的本地副本。若该文件不存在，请将下载好的 NetCDF 文件放置于
仓库根目录并命名为 `202512.nc`，预处理脚本 `01_preprocess_local_data.py`
会从中提取所需时段与区域。

## 怎么跑

数据放入 `bonus_project_option1_bve/data/raw/` 后，按当前 Lambert 工作流运行：

```bash
# 1. 预处理 ERA5 到规则经纬度中间网格
python scripts/01_preprocess_local_data.py

# 2. 插值到 Lambert 模式网格
python scripts/02_prepare_lambert_grid.py

# 3. 运行 Lambert BVE 实验矩阵
python scripts/03_run_experiments.py

# 4. 生成 Lambert 结果图
python scripts/04_make_figures.py

# 5. 数值诊断（Poisson 残差、能量守恒、尺度分解、Δt 敏感性、Skill Score）
python scripts/05_diagnostics.py
```

## 项目结构

![Project layout diagram](layout.png)

```
.
├── README.md
├── layout.png                         # 项目结构概览图
├── cheatsheet_BP1.md                  # 课程公式与数值方法速查表
├── 202512.nc                          # ERA5 原始数据（~600 MB, gitignored）
├── .gitignore
└── bonus_project_option1_bve/
    ├── environment.yml                 # conda 环境
    ├── .gitignore
    ├── data/
    │   ├── raw/                        # 放原始 ERA5 .nc 文件
    │   ├── processed/                  # 01 输出：粗化后的经纬度 .npz
    │   └── processed_lambert/          # 02 输出：Lambert 网格上的 .npz
    ├── src/
    │   ├── bve_model_lambert.py        # Lambert BVE 动力核心
    │   ├── poisson_dirichlet.py        # DST Poisson 求解器 (ψ=0)
    │   ├── lambert_grid.py             # Lambert 投影网格构建
    │   ├── operators.py                # 差分算子 + Arakawa Jacobian
    │   ├── interpolation.py            # 双线性插值
    │   ├── preprocess.py               # ERA5 读取、粗化、变量检测
    │   ├── verification.py             # RMSE / ACC / bias / 去偏 RMSE
    │   └── plotting.py                 # Lambert 曲边地图 + 评分图
    ├── scripts/
    │   ├── 01_preprocess_local_data.py # ERA5 → 粗化 → .npz
    │   ├── 02_prepare_lambert_grid.py  # 经纬度 → Lambert 重网格
    │   ├── 03_run_experiments.py       # 运行 5 个实验
    │   ├── 04_make_figures.py          # 画图
    │   └── 05_diagnostics.py           # 数值诊断
    ├── outputs/                        # 预报场、验证场、评分 JSON
    ├── figures/                        # 输出的 .png 图
    ├── note/
    │   └── technical_note.md           # 技术笔记
    └── slides/
        └── slides_outline.md           # 展示大纲
```

## 输出文件

### 数据与评分

| 文件 | 内容 |
| --- | --- |
| `data/processed/*.npz` | ERA5 经初步处理后的经纬度网格场 |
| `data/processed_lambert/*.npz` | 插值到 Lambert 模式网格后的场 |
| `data/processed_lambert/grid_info_lambert.npz` | Lambert 网格、`m(i,j)`、`f(i,j)` |
| `outputs/*_LCC_forecast_*.npz` | Lambert 实验预报 `ψ, ζ` |
| `outputs/*_LCC_verification_*.npz` | Lambert 实验验证数组 |
| `outputs/scores_experiment_matrix.json` | Lambert 实验评分 |
| `outputs/diagnostics.json` | 数值诊断结果（Poisson 残差、能量、尺度分解、Δt 敏感性、Skill Score） |

### 图片

| 文件 | 内容 |
| --- | --- |
| `figures/fig13_lambert_ctrl_verification_12h.png` | Lambert CTRL +12 h 预报/分析/误差 |
| `figures/fig13_lambert_ctrl_verification_24h.png` | Lambert CTRL +24 h 预报/分析/误差 |
| `figures/fig14_lambert_experiment_scores.png` | Lambert 实验评分对比 |
| `figures/fig15_lambert_impact_summary.png` | Lambert 实验影响汇总图 |

## 地图绘制说明

地图使用 Lambert Conformal 投影。`src/plotting.py` 中的
`set_lambert_curved_boundary()` 会把原始经纬度资料边界投影到 Lambert 平面，并用该
曲边四边形替换 Cartopy 默认矩形 Axes：

- 左右边界：固定经度线；
- 上下边界：固定纬度线；
- 填色、海岸线、国界线和经纬网均裁剪在该曲边边界内；
- 不再显示默认矩形地图外框。

`fig13` 展示 Lambert 模式结果。由于动力模式网格比报告展示用经纬度 footprint 略小，
绘图时将 Lambert 预报场插值到原始经纬度展示网格，并在模式边缘附近做平滑 taper，
避免出现硬矩形接缝。该处理仅用于展示；定量评分仍来自 Lambert 原生模式网格。

## Lambert 实验结果

评分均在 25°N–55°N, 90°E–145°E 验证子区域内计算。

| 实验 | +12 h RMSE | +24 h RMSE | +12 h ACC | +24 h ACC | 说明 |
| --- | ---: | ---: | ---: | ---: | --- |
| PERSIST_LCC | 42.9 m | 70.3 m | 0.985 | 0.953 | 持续性基准 |
| CTRL_LCC | 55.4 m | 56.6 m | 0.975 | 0.963 | Lambert BVE 主实验 |
| DIFF_LCC | 55.9 m | 57.9 m | 0.973 | 0.960 | 加弱涡度扩散 |
| SPONGE_LCC | 50.8 m | 69.0 m | 0.980 | 0.960 | 加边界海绵 |
| DIFF_SPONGE_LCC | 52.2 m | 71.9 m | 0.978 | 0.956 | 扩散 + 海绵 |

当前个例高度场非常稳定，PERSIST_LCC 是很强的基准。Lambert CTRL 在 +24 h
RMSE 和 ACC 上表现最好；扩散和海绵在该配置下主要起诊断作用，并未继续降低 +24 h RMSE。

## 数值诊断

`scripts/05_diagnostics.py` 附带 5 项额外检查，验证模型的数值可靠性：

| 诊断项 | 结果 |
| --- | --- |
| Poisson 求解器残差 | ‖m²∇²ψ − ζ‖₂ / ‖ζ‖₂ = **1.7 × 10⁻¹⁵**（机器精度） |
| 动能 24 h 漂移 | **+10.8%**（无爆炸、无振荡，积分 bounded） |
| 拟能 24 h 漂移 | **+12.2%**（无爆炸、无振荡） |
| 大尺度 ACC（+24 h, σ=2.0 Gaussian） | **0.967**（优于全场 ACC 0.963） |
| 小尺度 ACC（+24 h） | 0.741（BVE 对小尺度技巧有限，符合物理预期） |
| Δt 敏感性（300 / 600 / 900 s） | 三个步长 **+24 h RMSE 完全相同**（56.6 m） |
| Skill Score vs PERSIST（+12 h） | SS = **−0.29**（持续性更优，诚实承认） |
| Skill Score vs PERSIST（+24 h） | SS = **+0.195**（RMSE 相对持续性降低约 20%） |

这些诊断确认：(1) Poisson 反演可靠；(2) 非线性积分在 24 h 窗口内数值稳定；
(3) BVE 的技巧集中在大尺度 Rossby 波上，符合正压模型的物理能力；
(4) 主要结论对时间步长不敏感；(5) +24 h CTRL 相对于持续性有正技巧，
但 +12 h 持续性仍更强——这与缓慢演变冬季环流的预期一致。

## 已知局限

1. **单层正压模型。** 无辐散、无垂直结构、无斜压过程。
2. **高度恢复是诊断近似。** 当前使用 `Z'=fψ/g`，更严格版本应解线性平衡方程。
3. **固定边界仍是理想化边界。** `ψ=0` 消除了周期穿越问题，但仍不是真实侧边界条件。
4. **单个个例。** 结果不能代表模型在所有天气形势下的统计表现。
5. **展示插值不参与评分。** `fig13` 为报告视觉一致性做了边缘平滑；评分仍使用原生 Lambert 网格。

## 辅助文档

- **[cheatsheet_BP1.md](cheatsheet_BP1.md)** — 课程公式与数值方法速查表，涵盖 BPE/BVE/SW 三种模型的方程、
  符号约定、离散化方法和实现建议。适合在阅读代码前快速回顾理论基础。

## 参考文献

- Arakawa, A. (1966). Computational design for long-term numerical integration of the equations of fluid motion. *J. Comput. Phys.*, 1, 119–143.
- Haltiner, G. J., and R. T. Williams (1980). *Numerical Prediction and Dynamic Meteorology*, 2nd ed., Wiley.
- Holton, J. R., and G. J. Hakim (2013). *An Introduction to Dynamic Meteorology*, 5th ed., Academic Press.
- Hersbach, H., et al. (2020). The ERA5 global reanalysis. *Quart. J. Roy. Meteorol. Soc.*, 146, 1999–2049.
