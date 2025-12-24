# Emulator Input File Design

## 1. 训练输入文件格式 (`emulator_train.in`)

```fortran
NAMELIST
! Emulator Training Configuration

&EMULATOR_TRAIN
  base_input = 'test_d58Ni_cdcc.in'   ! 基础CDCC输入文件
  n_samples = 50                       ! LHS采样数
  n_basis = 0                          ! 基函数数 (0=GS自动确定)
  method = 'gs'                        ! 'svd' 或 'gs'
  tol = 1.0d-8                         ! GS收敛精度
  output_file = 'emulator_dNi.dat'     ! 保存文件名
/

! 参数范围定义
! 格式: kp1_kp2_param (势类型_编号_参数名)
! kp1: 't'=target, 'x'=core-target, 'b'=breakup, 'p'=p-n
! param: uv, rv, av, uw, rw, aw, wd, rwd, awd

&PARAM_RANGE  name='t_1_uv', min=45.0, max=60.0 /    ! proton-Target Vreal
&PARAM_RANGE  name='t_1_wd', min=5.0, max=12.0 /     ! proton-Target Wsurf
&PARAM_RANGE  name='x_1_uv', min=40.0, max=55.0 /    ! neutron-Target Vreal
&PARAM_RANGE  name='x_1_wd', min=6.0, max=12.0 /     ! neutron-Target Wsurf
&PARAM_RANGE /
```

## 2. 预测输入文件格式 (`emulator_predict.in`)

```fortran
NAMELIST
! Emulator Prediction Configuration

&EMULATOR_PREDICT
  emulator_file = 'emulator_dNi.dat'   ! 训练好的emulator
  base_input = 'test_d58Ni_cdcc.in'    ! 基础输入
  predict_method = 'matrix'            ! 'rbf' 或 'matrix'
  output_format = 'cross_section'      ! 输出类型
/

! 要预测的参数组合 (与PARAM_RANGE对应)
&PARAM_SET  t_1_uv=52.0, t_1_wd=8.0, x_1_uv=48.0, x_1_wd=9.0 /
&PARAM_SET  t_1_uv=55.0, t_1_wd=7.5, x_1_uv=50.0, x_1_wd=8.5 /
&PARAM_SET  t_1_uv=50.0, t_1_wd=9.0, x_1_uv=46.0, x_1_wd=10.0 /
&PARAM_SET /
```

## 3. 参数命名规则

| 势类型 (kp1) | 说明 | 示例参数名 |
|-------------|------|-----------|
| `t` | Target (proton-Target) | `t_1_uv`, `t_1_wd` |
| `x` | Core-Target (neutron-Target) | `x_1_uv`, `x_1_wd` |
| `b` | Breakup channel | `b_1_uv`, `b_2_uv` |
| `p` | p-n potential | `p_1_uv` |
| `a` | Adiabatic | `a_1_uv` |

| 参数名 | 说明 |
|-------|------|
| `uv` | 实部势深 (Woods-Saxon) |
| `rv` | 实部半径参数 |
| `av` | 实部扩散参数 |
| `uw` | 虚部体积项势深 |
| `rw` | 虚部体积项半径 |
| `aw` | 虚部体积项扩散 |
| `wd` | 虚部表面项势深 |
| `rwd` | 虚部表面项半径 |
| `awd` | 虚部表面项扩散 |

## 4. 数据结构

### 4.1 param_range_t
```fortran
type :: param_range_t
    character(len=32) :: name      ! 参数名 (如 't_1_uv')
    real(dpreal) :: min_val        ! 最小值
    real(dpreal) :: max_val        ! 最大值
    ! 解析后的信息
    character(len=4) :: kp1        ! 势类型
    integer :: kp2                 ! 势编号
    character(len=8) :: param      ! 参数名
end type
```

### 4.2 param_set_t
```fortran
type :: param_set_t
    integer :: n_params
    character(len=32), allocatable :: names(:)
    real(dpreal), allocatable :: values(:)
end type
```

## 5. 工作流程

### 训练阶段
```
1. 读取 emulator_train.in
2. 解析 EMULATOR_TRAIN namelist
3. 读取所有 PARAM_RANGE
4. 读取 base_input (完整CDCC设置)
5. LHS采样 n_samples 个参数组合
6. 对每个样本:
   a. 修改势参数 (根据param_range映射)
   b. 调用完整CDCC计算 (dbmm_coupled_channel_cc)
   c. 存储系数向量 c
7. 训练emulator (SVD或GS)
8. 保存到 output_file
```

### 预测阶段
```
1. 读取 emulator_predict.in
2. 加载 emulator_file
3. 读取所有 PARAM_SET
4. 对每个参数组合:
   a. 构建势 (根据参数值)
   b. 调用 emulator_predict
   c. 计算并输出截面
```

## 6. 文件结构

```
cookie/
├── pmm/
│   ├── emulator_io.F90        ! 输入解析模块
│   ├── param_mapping.F90      ! 参数映射模块
│   ├── emulator_train.F90     ! 训练主程序
│   ├── emulator_predict.F90   ! 预测主程序
│   └── Makefile
├── dbmm_emulator.F90          ! 核心emulator模块 (已有)
└── test/
    ├── emulator_train_dNi.in  ! 示例训练输入
    └── emulator_predict_dNi.in ! 示例预测输入
```
