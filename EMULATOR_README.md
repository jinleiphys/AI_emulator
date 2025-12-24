# DBMM Emulator 使用说明

## 概述

DBMM Emulator 是一个基于降维的快速预测器，用于加速CDCC（Continuum-Discretized Coupled-Channels）计算。通过离线训练学习解空间的低维结构，在线预测时可以大幅减少计算时间。

## 文件结构

```
pmm/
├── emulator_train.F90      # 训练程序
├── emulator_predict.F90    # 预测程序
├── dbmm_emulator.F90       # 核心emulator模块 (在 cookie/ 目录)
├── emulator_io.F90         # 输入输出模块
├── param_mapping.F90       # 参数映射模块
└── EMULATOR_README.md      # 本文档
```

---

## 1. 降维方法

### 1.1 SVD/PCA方法 (`method = 'svd'`)

**函数**: `emulator_finalize_training()`

**原理**:
- 对训练快照矩阵 S = [c₁, c₂, ..., cₙ] 进行奇异值分解 (SVD)
- S = U Σ V^H
- 选取前 n_basis 个左奇异向量作为reduced basis: X = U(:, 1:n_basis)

**特点**:
- ✅ 最优降维（Eckart-Young定理）
- ✅ 给出明确的截断误差估计 (σₖ₊₁/σ₁)
- ❌ 需要预先指定 n_basis
- ❌ 无法增量更新

**配置示例**:
```fortran
&EMULATOR_TRAIN
  n_basis = 20        ! 必须指定基函数数目
  method = 'svd'
/
```

### 1.2 Greedy Gram-Schmidt方法 (`method = 'gs'`)

**函数**: `emulator_finalize_training_gs()`

**原理**:
1. 找到范数最大的快照 → φ₁
2. 将所有快照投影到φ₁，找残差最大的 → 正交化得到φ₂
3. 重复直到残差小于容差tol

**特点**:
- ✅ 自动确定基函数数目
- ✅ 可增量更新 (`emulator_add_snapshot_gs`)
- ✅ 每个基向量对应一个具体的训练样本
- ❌ 不是全局最优

**配置示例**:
```fortran
&EMULATOR_TRAIN
  n_basis = 0         ! 0 = 自动确定
  method = 'gs'
  tol = 1.0d-10       ! 残差容差
/
```

### 1.3 方法对比

| 特性 | SVD | Gram-Schmidt |
|------|-----|--------------|
| 最优性 | 全局最优 | 贪婪近似 |
| 基数目 | 需预设 | 自动确定 |
| 增量更新 | 不支持 | 支持 |
| 计算复杂度 | O(mn²) | O(mn·k) |

---

## 2. 预测方法

### 2.1 Matrix Reconstruction (Galerkin投影)

**原理**:
```
M_red = X^H · M(θ) · X    # 投影矩阵
b_red = X^H · b           # 投影源项
c_red = M_red⁻¹ · b_red   # 解reduced系统
c ≈ X · c_red             # 重构解
```

**特点**:
- ✅ 物理上准确（在子空间内）
- ❌ 假设解在span{X}中，远离训练样本时误差大

### 2.2 RBF Coefficient Interpolation (系数插值)

**全称**: Radial Basis Function (径向基函数) 插值

**原理**:
```
wᵢ(θ) = exp(-|θ - θᵢ|² / 2σ²)   # 高斯RBF核权重
α(θ) = Σᵢ wᵢ(θ) · αᵢ / Σᵢ wᵢ    # 归一化加权平均
c ≈ X · α(θ)                     # 重构解
```

**特点**:
- ✅ 对参数变化更robust
- ✅ 误差更稳定
- ❌ 不使用新参数的物理信息
- ❌ 不提供不确定性估计（与高斯过程回归不同）

**注**: 使用Nadaraya-Watson核回归，也称RBF插值。

### 2.3 实测对比

| 参数 | Matrix | RBF插值 | 说明 |
|------|--------|---------|------|
| nominal | 4.0E-06 | 5.3E-03 | Matrix在训练点附近极准 |
| 偏离10% | 2.7E-02 | 6.3E-03 | RBF更稳定 |
| 偏离20% | 7.1E-02 | 7.6E-03 | RBF明显更好 |

**建议**: 参数扫描用RBF方法；精确计算用Matrix方法+密集训练。

---

## 3. 训练程序 (emulator_train)

### 3.1 使用方法

```bash
cd cookie/pmm
echo "path/to/train_config.in" | ./emulator_train
```

### 3.2 输入文件格式

```fortran
NAMELIST
! 注释行

&EMULATOR_TRAIN
  base_input = '../test/test_d58Ni_cdcc.in'  ! 基础CDCC输入文件
  n_samples = 30                              ! LHS采样数
  n_basis = 0                                 ! 基函数数 (0=自动)
  method = 'gs'                               ! 'svd' 或 'gs'
  tol = 1.0d-10                               ! GS方法的残差容差
  output_file = 'emulator.dat'                ! 输出文件名
/

! 参数范围定义 (可多个)
! 格式: name='kp1_kp2_param'
!   kp1: 't'=proton-Target, 'x'=neutron-Target, 'b'=breakup, 'p'=p-n
!   kp2: 势能索引 (通常为1)
!   param: uv, rv, av, uw, rw, aw, wd, rwd, awd, vd, rvd, avd, rc

&PARAM_RANGE  name='t_1_uv', min_val=45.0, max_val=60.0 /
&PARAM_RANGE  name='t_1_wd', min_val=5.0, max_val=12.0 /
&PARAM_RANGE /   ! 空的PARAM_RANGE表示结束
```

### 3.3 参数说明

| 参数 | 类型 | 说明 |
|------|------|------|
| `base_input` | string | 基础CDCC计算的输入文件路径 |
| `n_samples` | int | Latin Hypercube采样点数 |
| `n_basis` | int | 基函数数目，0=自动确定(仅GS) |
| `method` | string | 'svd' 或 'gs' |
| `tol` | real | GS方法残差容差 |
| `output_file` | string | 输出emulator数据文件 |

### 3.4 输出

- **屏幕输出**: 训练进度、奇异值谱、基函数数目
- **文件输出**: 二进制emulator数据文件 (`.dat`)

---

## 4. 预测程序 (emulator_predict)

### 4.1 使用方法

```bash
cd cookie/pmm
echo "path/to/predict_config.in" | ./emulator_predict
```

### 4.2 输入文件格式

```fortran
NAMELIST
! 注释行

&EMULATOR_PREDICT
  emulator_file = 'emulator.dat'              ! 训练好的emulator文件
  base_input = '../test/test_d58Ni_cdcc.in'   ! 基础输入(用于初始化物理系统)
  predict_method = 'matrix'                    ! 'matrix' 或 'rbf'
  output_format = 'cross_section'              ! 输出格式
/

! 预测参数集 (可多个)
&PARAM_SET  t_1_uv=53.3, t_1_wd=7.8 /
&PARAM_SET  t_1_uv=52.0, t_1_wd=8.0 /
&PARAM_SET  t_1_uv=54.0, t_1_wd=7.5 /
&PARAM_SET /   ! 空的PARAM_SET表示结束
```

### 4.3 参数说明

| 参数 | 类型 | 说明 |
|------|------|------|
| `emulator_file` | string | 训练好的emulator文件 |
| `base_input` | string | 基础CDCC输入(初始化物理系统) |
| `predict_method` | string | 'matrix'=Galerkin, 'rbf'=系数插值 |
| `output_format` | string | 输出格式 |

### 4.4 输出

- **屏幕输出**: 每个参数集的预测结果和误差
- 当前版本同时计算Matrix和RBF两种方法的误差进行对比

---

## 5. 完整示例

### 5.1 训练

```bash
# 创建训练配置
cat > train.in << 'EOF'
NAMELIST
&EMULATOR_TRAIN
  base_input = '../test/test_d58Ni_cdcc.in'
  n_samples = 50
  n_basis = 0
  method = 'gs'
  tol = 1.0d-10
  output_file = 'my_emulator.dat'
/
&PARAM_RANGE  name='t_1_uv', min_val=45.0, max_val=60.0 /
&PARAM_RANGE  name='t_1_wd', min_val=5.0, max_val=12.0 /
&PARAM_RANGE /
EOF

# 运行训练
echo "train.in" | ./emulator_train
```

### 5.2 预测

```bash
# 创建预测配置
cat > predict.in << 'EOF'
NAMELIST
&EMULATOR_PREDICT
  emulator_file = 'my_emulator.dat'
  base_input = '../test/test_d58Ni_cdcc.in'
  predict_method = 'matrix'
  output_format = 'cross_section'
/
&PARAM_SET  t_1_uv=53.3, t_1_wd=7.8 /
&PARAM_SET  t_1_uv=50.0, t_1_wd=10.0 /
&PARAM_SET /
EOF

# 运行预测
echo "predict.in" | ./emulator_predict
```

---

## 6. 性能参考

基于 d+58Ni CDCC测试 (nlag=80, nch=19, ntot=1520):

| 指标 | 值 |
|------|-----|
| 训练时间 (30样本) | ~8秒 |
| 单次DBMM求解 | ~280ms |
| Emulator预测 | ~80ms |
| 加速比 | ~3.5x |
| 降维比 | 1520 → 20 (76x) |

---

## 7. 注意事项

1. **缓存问题**: 修改势能参数后必须调用 `cleanup_potential_cache()` 清除缓存
2. **参数范围**: 预测参数应在训练范围内，外推精度会显著下降
3. **基函数数目**: SVD方法需要预设，建议先用GS方法自动确定合适数目
4. **内存**: 大规模计算时注意快照存储 (~16 × ntot × n_samples 字节)

---

## 8. 参考文献

Liu, Lei, Ren, "Reduced-basis emulator for CDCC calculations", Phys. Lett. B 858 (2024) 139070
