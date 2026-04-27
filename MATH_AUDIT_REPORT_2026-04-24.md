# 数学原理与实现细查报告（2026-04-24）

## 1) 目标与范围
本报告针对你提出的核心问题做两层审计：
1. **数学原理是否正确**：`∂R/∂w · dw/dx = -∂R/∂x` 是否成立。
2. **代码计算过程是否有错**：逐环节检查 `RHS`、`Jv`、线性性、GMRES 过程输出。

---

## 2) 数学原理核查（理论层）
设离散稳态系统为：
\[
R(w, x)=0
\]
其中 `w` 是状态（本代码中 primitive 变量），`x` 是设计参数（如形变参数、Mach 等）。
对 `x` 求导：
\[
\frac{\partial R}{\partial w}\frac{dw}{dx}+\frac{\partial R}{\partial x}=0
\Rightarrow
\frac{\partial R}{\partial w}\frac{dw}{dx}=-\frac{\partial R}{\partial x}
\]

该式在工程伴随/切线灵敏度中是标准形式。

---

## 3) 代码路径核查（实现层）
### 3.1 左端算子 `∂R/∂w · v`
- 由 `TANGENT_MATVEC` 调 `COMPUTE_STEADY_RESIDUAL_D`，以 `cellprimitivesd=v` 为 seed 构建 `Jv`。
- 输出来自 `fluxresiduals_slnd` 的同一残差体系。

### 3.2 右端项 `-∂R/∂x`
- `TANGENT_RHS_I` 对参数方向设 `data_4d137d(1,i_param)=1`，同样调 `COMPUTE_STEADY_RESIDUAL_D`。
- 然后显式赋值 `b_i = -fluxresiduals_slnd`。

### 3.3 变量空间一致性
- `pc_dafoam_consistency_mode=.TRUE.` 路径下，GMRES 和 PC 构建都会禁用 conservative/primitive 混合变换，保持同一状态空间。

---

## 4) 运行时“逐步大量输出”核查结果
运行口径：`endtime=600`, `maxiter=10`, `m_restart=50`。

### 4.1 着色同色互扰/填块正确性（新增）
- 日志：`[PC-COLOR-CHECK] sampled_colors=3 rel_err=0.0`
- 检查方法：对采样颜色做**单种子参考**（每次只激活一个单元一个变量）并与“同色批量种子后填入 pc_blk 的块列”逐项比对。
- 结论：在本次采样下未发现同色互相污染或填块错位（至少在采样集合上成立）。

### 4.2 RHS 正确性（新增 FD 对照）
- 日志：`[RHS-FD] i_param=1 rel_err=1.3913896972442811E-007`
- 结论：`-∂R/∂x` 的 AD 计算与 FD 高度一致（1e-7 级相对误差）。

### 4.3 Jv 正确性（新增 FD 对照）
- 日志：`[Jv-FD] rel_err=4.0561877379362532E-004`
- 结论：`∂R/∂w·v` 与 FD 在 4e-4 量级误差，属于可接受一致性范围（考虑非线性离散/FD 截断误差）。

### 4.4 线性性检查（已有）
- repeat: `0.0`
- homogeneity: `7.5965e-6`
- additivity: `2.8101e-6`
- 结论：算子线性化行为正常。

### 4.5 M^{-1} 求解精度（新增）
- 日志：`[PC-MINV-CHECK] ||rhs-A0*ApplyPC(rhs)||/||rhs|| = 4199.64`
- 公式含义：设 `z = ApplyPC(rhs)`，则该指标是
  \[
  \frac{\|rhs - A0 z\|}{\|rhs\|}
  \]
  也就是“把 `ApplyPC` 当成 `A0^{-1}` 后，代回去看残差有多大”。
- 工程解释：
  - 若 `ApplyPC \approx A0^{-1}`，则 `A0 z \approx rhs`，指标应接近 `0`。
  - 指标 `~O(1)`：近似逆质量一般，GMRES 可能勉强可用。
  - 指标 `>>1`：近似逆质量很差，预条件基本失效，甚至会放大误差。
- 本次 `4199.64` 说明：`ApplyPC` 生成的 `z` 代回 `A0` 后完全不能重构 `rhs`，属于严重失配。
- 结论：即使 `RHS/Jv` 正确、着色采样无互扰，若 ApplyPC 与目标矩阵严重不一致，GMRES 仍会卡在高相对残差。

### 4.7 根因定位（对照 ADflow/PETSc 线性预条件链）
- 当前实现是**块 ILU(0)**，每个网格单元对角块为 **5×5**，前代/回代都在 5×5 块系统上进行（不是标量 ILU）。
- 本次新增日志会打印：
  - `worst_diag_block cell=... cond_est=...`
  - 对应 `Uii(worst)` 与 `inv(Uii)(worst)` 的 5×5 数值。
- 若 `cond_est` 极大（例如 `>=1e10`）且 `inv(Uii)` 元素出现异常放大，则说明主要是 **ILU 对角块近奇异/数值爆炸**，不是 GMRES 初值 `x0` 本身导致。
- 按 ADflow/PETSc 思路，`M^{-1}` 应保持线性、平稳且不放大向量范数；若 `||z_i||` 与局部 `||r_i||` 同时异常，优先归因到 ILU 链条（组装/因子/回代）而非外层 Krylov 初值。

### 4.6 GMRES 指标
- `outer=1 true||r||/||b|| = 0.7384892912`
- `outer=10 true||r||/||b|| = 0.7322944608`
- 结论：残差下降很慢，问题仍在预条件质量/谱聚类，而非“方程写错号”或“微分算子明显错误”。

---

## 5) 结论（针对你的核心怀疑）
1. **大方向数学原理正确**：`∂R/∂w dw/dx=-∂R/∂x` 在本代码路径下成立。  
2. **关键微分算子未发现根本错误**：
   - RHS-FD 误差 `~1.39e-7`（很强一致）
   - Jv-FD 误差 `~4.06e-4`（可接受）
3. **新增关键发现**：ApplyPC 一致性仍很差（`||rhs-A0*ApplyPC(rhs)||/||rhs||≈4.20e3`），较此前 2.53e5 明显改善但仍远未达可用。
4. 当前 `~0.7` 卡住的最主要可疑项：
   - ApplyPC 与 A0 的一致性很差（比“同色互扰”更可疑）
   - 预处理后谱聚类未改善，导致 GMRES 外循环降幅小。

---

## 6) 下一步建议（继续 ADflow 复现路线）
在不关闭 colored-AD 的前提下，下一步建议按 ADflow 思路做“单变量 A/B”而非多项混改：
1. 固定 colored-AD 路线；
2. 只改“PC 稀疏连接层级 + 过滤阈值”一项；
3. 每次只看 `outer=1` 和 `outer=10` 两个指标，记录下降率。
