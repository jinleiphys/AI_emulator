# Unified Basis Multi-J Emulator

## Overview

Attempted to implement a unified basis approach for the multi-J DBMM emulator, allowing a single set of basis functions to work across ALL J values.

**STATUS: ❌ 方案不可行 - 已放弃，改用Per-J Emulator**

## 失败原因分析

### 根本问题：边界条件不兼容

不同J值的解住在**不同的向量空间**：
```
c(θ, J=0)  ∈ ℂ^{80×19}   (1520维)
c(θ, J=20) ∈ ℂ^{80×37}   (2960维)
```

更关键的是，不同channel有本质不同的边界条件：
- **Open channels**: oscillatory (F + iG Coulomb functions)
- **Closed channels**: exponentially decaying (Whittaker W)

这两种radial behavior在任何L²基下都无法同时高效表示。

### Zero-Padding方法的问题

1. 引入人工结构，破坏物理意义
2. 不同J的channel数差异大（18-37）
3. SVD无法区分真实物理结构和padding artifacts

### 测试结果

| J值 | 平均相对误差 |
|-----|-------------|
| J=0  | 90.4% |
| J=1  | 83.3% |
| J=5  | 101% |
| J=20 | 232% |

**误差~100%，完全不可接受**

## 其他方案分析

| 方案 | 可行性 | 说明 |
|------|--------|------|
| Zero-padding unified basis | ❌ | 边界条件不兼容 |
| Channel-wise低秩分解 | ❌ | 低秩假设不成立 |
| Interior Matching | ⚠️ | 减轻但不解决维度不匹配 |
| S-matrix Emulation | ⚠️ | 变成ML surrogate，非physics emulator |
| Complex Scaling | ✓ | 需要大改DBMM核心代码 |
| **Per-J Emulator** | ✓ | **推荐方案** |

### 为什么Complex Scaling能work

参考PLB论文的Complex Scaling方法：
```
r → r·e^(iθ)  for r > R₀

Before: e^(ikr)     → oscillatory, non-L²
After:  e^(ikr·e^(iθ)) → exponentially decaying, L²
```

Complex scaling把所有partial wave的边界条件统一成exponentially decaying，因此可以用同一套L²基展开。DBMM没有这个变换，所以unified basis不可行。

## 最终决定

**采用Per-J Emulator方案**

优点：
- 每个J独立的reduced basis，物理正确
- 精度已验证良好
- 代码已实现并测试
- 实际使用足够高效

## 已清理的代码

以下unified basis相关代码已从代码库中移除：
- `test_unified_accuracy.F90` - 测试程序
- `emulator_train_unified.in` - 配置文件
- `emulator_train_unified_small.in` - 配置文件
- `emulator_unified_small.dat` - 训练数据
- `dbmm_emulator.F90`中的unified basis函数（保留per-J功能）

## 参考

- PLB论文：Complex Scaling Emulator for scattering
- Eigenvector Continuation / Reduced Basis Method文献
