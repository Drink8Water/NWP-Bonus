# A 12–24 h Lambert-Grid Barotropic Vorticity Forecast Experiment for an East Asian Winter Z500 Case

用 Lambert 保角投影有限区域网格上的正压涡度方程，对 2025 年 12 月底一次东亚冬季
500 hPa 环流过程做 12–24 h 简化动力预报实验。初始场和验证场均来自本地 ERA5
再分析资料。

## 个例

| 项目               | 值                           |
| ------------------ | ---------------------------- |
| 初始时间           | 2025-12-30 00 UTC            |
| +12 h / +24 h 验证 | 2025-12-30 12 UTC / 2025-12-31 00 UTC |
| 等压面             | 500 hPa                      |
| 全区域             | 15°N–65°N, 60°E–170°E        |
| 验证子区域         | 25°N–55°N, 90°E–145°E        |

2025 年 12 月底东亚—西伯利亚中纬度有明显长波槽脊，500 hPa 环流以旋转风为主，
适合用正压模型。但该过程演变缓慢，是高度持续型环流，持续性本身就是一个很强的
baseline。

## 数据

December 2025 的 ERA5 再分析资料，500 hPa geopotential、u 风和 v 风。
NetCDF 文件放到 `data/raw/`，不需要联网。

## 模型

Lambert 保角投影有限区域网格上的正压涡度方程（BVE）：

$$
\frac{\partial\zeta}{\partial t} = -m^2 J(\psi, \zeta+f) + \nu m^2\nabla^2\zeta - \alpha(\zeta-\zeta_{\text{ref}})
$$

$$
\zeta = m^2\nabla^2\psi
$$

其中 m 是 Lambert 投影的 map factor，f 是随纬度变化的 Coriolis 参数。
ψ = 0 在侧边界上（Dirichlet 条件），用离散正弦变换（DST）求解 Poisson 方程。

- 空间离散：二阶中心差分，Arakawa (1966) Jacobian（能量 + 拟能守恒）
- 时间积分：RK4，Δt = 600 s
- Poisson 求解：DST，天然满足 ψ = 0 边界条件
- 可选阻尼：涡度扩散（ν）和边界海绵松弛（α）

## 安装

```bash
conda env create -f environment.yml
conda activate bve-nwp
```

或 pip：

```bash
pip install numpy scipy xarray netCDF4 matplotlib cartopy
```

## 怎么跑

```bash
python scripts/01_preprocess_local_data.py    # ERA5 → 粗化 → .npz
python scripts/02_prepare_lambert_grid.py     # 重网格到 Lambert 投影
python scripts/03_run_experiments.py          # 运行实验
python scripts/04_make_figures.py             # 画图
```

## 实验设计

| 实验         | 说明                                        |
| ------------ | ------------------------------------------- |
| PERSIST      | 持续性 baseline：初始高度异常保持不变        |
| CTRL         | 原始 BVE，无扩散、无海绵                    |
| DIFF         | BVE + 涡度扩散（ν = 2.5×10⁴ m²/s）          |
| SPONGE       | BVE + 边界海绵松弛（8 格点, τ = 6 h）       |
| DIFF_SPONGE  | 扩散 + 海绵联合                             |

评分在验证子区域 25°N–55°N, 90°E–145°E 上计算：RMSE、bias、去偏 RMSE、
高度异常相关系数（ACC）。

## 结果

| 实验            | +12h RMSE (m) | +24h RMSE (m) | +24h ACC |
| --------------- | ------------: | ------------: | -------: |
| PERSIST_LCC     | 42.9          | 70.3          | 0.953    |
| CTRL_LCC        | 55.4          | 56.6          | 0.963    |
| DIFF_LCC        | 55.9          | 57.9          | 0.960    |
| SPONGE_LCC      | 50.8          | 69.0          | 0.960    |
| DIFF_SPONGE_LCC | 52.2          | 71.9          | 0.956    |

几个发现：

1. **CTRL 在 +24 h 略优于持续性。** PERSIST 的 +24 h RMSE 为 70.3 m，CTRL 降到
   56.6 m，ACC 也从 0.953 提升到 0.963。对这个缓慢演变的冬季个例，持续性本身已经很
   强，但 Lambert BVE 仍然捕捉到了一部分有用的 Rossby 波演变信息。

2. **扩散和海绵对 Lambert 配置改善有限。** 在 Dirichlet 边界条件下，边界噪声
   不是主导误差来源。DIFF 的 +24 h RMSE（57.9 m）与 CTRL（56.6 m）接近。
   SPONGE 在 +12 h 表现最好（50.8 m），但 +24 h 退回到 69.0 m。

3. **CTRL 是积分模型中 +24 h 综合表现最好的。** 不需要额外调参，干净的有限区域
   BVE 在这个短时积分窗口中已经足够稳定。

## 输出

| 路径                                | 内容                      |
| ----------------------------------- | ------------------------- |
| `data/processed_lambert/`           | Lambert 网格上的 .npz 场  |
| `outputs/{EXP}_LCC_*.npz`           | 各实验预报场和验证场       |
| `outputs/scores_experiment_matrix.json` | 全部实验评分           |
| `figures/`                          | 预报验证图和评分对比图     |

## 已知局限

1. 无辐散，无斜压过程，无物理参数化
2. 预报时效短（24–48 h 后退化）
3. 仅一个个例，结果不能推广
4. 持续性在这个高度稳定的冬季个例中是最强 baseline

## 参考文献

- Arakawa, A. (1966). Computational design for long-term numerical integration
  of the equations of fluid motion. *J. Comput. Phys.*, 1, 119–143.
- Haltiner, G. J. & R. T. Williams (1980). *Numerical Prediction and Dynamic
  Meteorology*, 2nd ed., Wiley.
- Holton, J. R. & G. J. Hakim (2013). *An Introduction to Dynamic Meteorology*,
  5th ed., Academic Press.
- Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quart. J. Roy.
  Meteorol. Soc.*, 146, 1999–2049.
