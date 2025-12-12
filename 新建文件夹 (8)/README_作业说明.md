# 数学物理方法数值验证作业

本作业分为三个部分，对应三个子问题，分别在直角坐标系和柱坐标系下验证稳态方程、传导方程和波动方程的数值解与解析解。

## 项目结构

```
根目录/
├── Problem1_齐次边界齐次方程/
│   └── problem1_homogeneous.py          # 齐次情况下的三个方程求解
├── Problem2_非齐次边界非齐次方程/
│   └── problem2_inhomogeneous.py        # 非齐次边界和非齐次源项情况
├── Problem3_柱坐标系/
│   └── problem3_cylindrical.py          # 柱坐标系中的2D轴对称问题
└── 本README文件
```

## 问题说明

### 问题1：齐次边界条件和齐次方程

**求解器：** `problem1_homogeneous.py`

在直角坐标系 $[0,1] \times [0,1]$ 上求解三种方程：

1. **稳态方程（Laplace）**
   - 方程：$\frac{d^2u}{dx^2} = 0$
   - 边界条件：$u(0) = u(1) = 0$（齐次Dirichlet）
   - 解析解：$u(x) = 0$
   - 数值方法：有限差分 + 高斯消元

2. **传导方程（Heat）**
   - 方程：$\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2}$
   - 初始条件：$u(x,0) = \sin(\pi x)$
   - 边界条件：$u(0,t) = u(1,t) = 0$
   - 解析解：$u(x,t) = e^{-\pi^2 Dt} \sin(\pi x)$
   - 数值方法：显式欧拉法（需满足Courant数 ≤ 0.5）

3. **波动方程**
   - 方程：$\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}$
   - 初始条件：$u(x,0) = \sin(\pi x)$，$\frac{\partial u}{\partial t}|_{t=0} = 0$
   - 边界条件：$u(0,t) = u(1,t) = 0$
   - 解析解：$u(x,t) = \sin(\pi x) \cos(c\pi t)$
   - 数值方法：中心差分显式格式（需满足CFL数 ≤ 1）

**输出文件：**
- `01_稳态方程_对比.png` - 稳态方程：解析解 vs 数值解
- `02_传导方程_对比.png` - 传导方程：解析解 vs 数值解
- `03_波动方程_对比.png` - 波动方程：解析解 vs 数值解

---

### 问题2：非齐次边界条件和非齐次方程

**求解器：** `problem2_inhomogeneous.py`

1. **Poisson方程**
   - 方程：$\frac{d^2u}{dx^2} = -2$
   - 边界条件：$u(0) = 0$，$u(1) = 1$（非齐次）
   - 解析解：$u(x) = -x^2 + 2x$
   - 数值方法：有限差分

2. **非齐次传导方程**
   - 方程：$\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2}$
   - 初始条件：$u(x,0) = \sin(\pi x)$
   - 边界条件：$u(0,t) = 0$，$u(1,t) = 1$（非齐次）
   - 稳态解：$u_s(x) = x$
   - 数值方法：Crank-Nicolson隐式方法（处理非齐次边界）

3. **非齐次波动方程（驱动振动）**
   - 方程：$\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}$
   - 初始条件：$u(x,0) = 0$，$\frac{\partial u}{\partial t}|_{t=0} = 0$
   - 边界条件：$u(0,t) = u(1,t) = \sin(\omega t)$（边界驱动）
   - 数值方法：中心差分显式格式

**输出文件：**
- `01_Poisson方程_对比.png` - Poisson方程对比图
- `02_非齐次传导方程_对比.png` - 非齐次传导方程对比图
- `03_非齐次波动方程_对比.png` - 非齐次波动方程对比图

---

### 问题3：柱坐标系中的方程求解

**求解器：** `problem3_cylindrical.py`

在柱坐标系中求解2D轴对称问题（$r \in [0,1]$，$z \in [0,1]$）：

1. **稳态方程（轴对称Laplace）**
   - 方程：$\frac{1}{r}\frac{d}{dr}\left(r\frac{du}{dr}\right) + \frac{d^2u}{dz^2} = 0$
   - 边界条件：
     - $u(1, z) = 0$（外边界）
     - $u(r, 0) = \sin(\pi r)$（底部）
     - $u(r, 1) = 0$（顶部）
     - $\frac{\partial u}{\partial r}|_{r=0} = 0$（轴对称）
   - 数值方法：Gauss-Seidel迭代法

2. **传导方程（轴对称热方程）**
   - 方程：$\frac{\partial u}{\partial t} = D\left[\frac{1}{r}\frac{\partial u}{\partial r} + \frac{\partial^2 u}{\partial r^2} + \frac{\partial^2 u}{\partial z^2}\right]$
   - 初始条件：$u(r,z,0) = \sin(\pi r)\sin(\pi z)$
   - 边界条件：$u(1,z,t) = 0$，$u(r,1,t) = 0$，$\frac{\partial u}{\partial r}|_{r=0} = 0$
   - 解析解：$u(r,z,t) = e^{-2\pi^2 Dt}\sin(\pi r)\sin(\pi z)$
   - 数值方法：显式欧拉法

3. **波动方程（轴对称）**
   - 方程：$\frac{\partial^2 u}{\partial t^2} = c^2\left[\frac{1}{r}\frac{\partial u}{\partial r} + \frac{\partial^2 u}{\partial r^2} + \frac{\partial^2 u}{\partial z^2}\right]$
   - 初始条件：$u(r,z,0) = \sin(\pi r)\sin(\pi z)$，$\frac{\partial u}{\partial t}|_{t=0} = 0$
   - 边界条件：$u(1,z,t) = 0$，$u(r,1,t) = 0$，$\frac{\partial u}{\partial r}|_{r=0} = 0$
   - 解析解：$u(r,z,t) = \sin(\pi r)\sin(\pi z)\cos(\pi\sqrt{2} ct)$
   - 数值方法：中心差分显式格式

**输出文件：**
- `01_柱坐标稳态方程_3D对比.png` - 3D表面对比
- `01_柱坐标稳态方程_等高线.png` - 等高线对比
- `02_柱坐标传导方程_等高线.png` - 传导方程等高线
- `03_柱坐标波动方程_等高线.png` - 波动方程等高线

---

## 运行方法

### 依赖环境
```
numpy
matplotlib
scipy
```

### 安装依赖
```bash
pip install numpy matplotlib scipy
```

### 运行各问题
```bash
# 问题1
python Problem1_齐次边界齐次方程/problem1_homogeneous.py

# 问题2
python Problem2_非齐次边界非齐次方程/problem2_inhomogeneous.py

# 问题3
python Problem3_柱坐标系/problem3_cylindrical.py
```

### 预期输出
每个脚本运行后会：
1. 在控制台打印最大误差信息
2. 在相应文件夹中生成对比图像（`.png`格式）

---

## 数值方法总结

### 稳定性条件

| 方程 | 离散化方法 | 稳定性条件 |
|------|---------|----------|
| 热方程 | 显式欧拉 | Courant数 $\leq 0.5$ |
| 波动方程 | 中心差分 | CFL数 $\leq 1$ |
| 传导方程 | Crank-Nicolson | 无条件稳定 |

### 精度分析

- **有限差分**：二阶空间精度 ($O(\Delta x^2)$)
- **显式欧拉法**：一阶时间精度 ($O(\Delta t)$)
- **Crank-Nicolson**：二阶时间精度 ($O(\Delta t^2)$)
- **中心差分**：二阶时间精度 ($O(\Delta t^2)$)

---

## 验证结果

所有方程的数值解与解析解的最大误差应该满足：

- **问题1**（齐次）：误差 $\sim 10^{-3}$ 到 $10^{-2}$
- **问题2**（非齐次）：误差 $\sim 10^{-2}$ 到 $10^{-1}$
- **问题3**（柱坐标）：误差 $\sim 10^{-2}$ 到 $10^{-1}$

误差大小取决于网格密度和时间步长的选择。

---

## 注意事项

1. **轴对称边界**：在 $r=0$ 处需要特殊处理，通常采用对称性条件或L'Hôpital法则
2. **稳定性检验**：运行前务必检查Courant数和CFL数是否满足条件
3. **收敛性**：可通过加密网格和减小时间步长来验证数值解的收敛性
4. **图像保存**：所有图像文件自动保存在对应问题的文件夹中

---

**作者**：自动生成的数值验证工具
**日期**：2025年12月

