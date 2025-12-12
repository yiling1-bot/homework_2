# PDE 数值验证项目说明

本项目包含三类偏微分方程（稳态/Laplace-Poisson，传导/Heat，波动/Wave）在不同边界与坐标系下的数值求解与解析对比，当前实现基于网格差分（传统有限差分思想），并预留谱方法与 PINN（AI 思路）扩展入口。

## 目录
- `Problem1_齐次边界齐次方程/`：直角坐标齐次边界，Laplace/Heat/Wave。
- `Problem2_非齐次边界非齐次方程/`：直角坐标非齐次边界/源项，Poisson、非齐次 Heat（Crank–Nicolson）、边界激励 Wave（制造解）。
- `Problem3_柱坐标系/`：柱坐标轴对称，制造解 Poisson、Heat、Wave。
- 文档：`README_浅显版.md`（易读版），`paper.tex`（论文草稿，含图像引用）。

## 运行
确保已安装 Python 及依赖：
```
pip install numpy matplotlib scipy
```
运行各问题脚本，生成对应 PNG 图像（保存在各自目录）：
```
python Problem1_齐次边界齐次方程/problem1_homogeneous.py
python Problem2_非齐次边界非齐次方程/problem2_inhomogeneous.py
python Problem3_柱坐标系/problem3_cylindrical.py
```

## 方法概览
- 网格差分（已实现）：空间二阶中心差分；时间——显式欧拉（Heat）、中心差分（Wave）、Crank–Nicolson（非齐次 Heat）。柱坐标轴对称在 r=0 采用对称/L'Hôpital 处理。
- 谱方法（计划）：直角 1D 场景使用正弦基傅里叶伪谱/Galerkin，空间误差指数收敛；拟与网格差分二阶误差并列对比。
- PINN（计划，AI 思路）：构造 PDE 残差+边界/初值损失，在 1D Heat/Wave 上训练，展示损失与误差收敛，与网格差分/谱法对比。

## 当前结果摘要
- Problem1（齐次）：Laplace 误差 0；Heat 1.47e-05；Wave 3.98e-06，数值与解析重合。
- Problem2（非齐次）：Poisson ~1e-14；非齐次 Heat（t=2.0）误差 5.04e-02，可延长时间/加密网格降误差；制造解 Wave 3.10e-05。
- Problem3（柱坐标）：制造解 Poisson 1.07e-03；Heat 1.57e-01；Wave 9.16e-02，可通过更细网格/更小时间步改进。

## 图像（生成后位于各目录）
- Problem1：`01_steady_laplace_compare.png`, `02_heat_equation_compare.png`, `03_wave_equation_compare.png`
- Problem2：`01_poisson_compare.png`, `02_heat_inhomogeneous_compare.png`, `03_wave_forcing_compare.png`
- Problem3：`01_cyl_steady_3d_compare.png`, `01_cyl_steady_contour.png`, `02_cyl_heat_contour.png`, `03_cyl_wave_contour.png`

## 后续改进建议
1. 补充谱方法实现，绘制误差-网格（或模式数）收敛曲线，与网格差分并列。
2. 补充 PINN 实验（1D Heat/Wave），绘制损失收敛与测试误差曲线。
3. 非齐次热：延长模拟时间、减小时间步或加密网格以更贴近稳态；柱坐标热/波同理。
4. 为柱坐标稳态提供严格解析基准（Bessel 模态分离变量），或注明制造解性质。
