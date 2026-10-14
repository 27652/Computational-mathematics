# RPN-CG 仓库完成度检查报告

检查依据：`ROPTLIB/RPNCG - latex/req.txt`。检查范围包括 ROPTLIB 中的 C++ 实现、测试入口、已有 LaTeX 文档，以及仓库根目录下的 MATLAB RPN-CG 实现。

## 结论

当前仓库已经完成 **RPN-CG 的 C++ 初步移植、基本测试入口和实现说明文档**，但没有完成需求中最关键的 **同口径效率/效果对比实验及其结果写入文档**。因此整体完成度建议评定为：

> **部分完成。算法代码和说明文档基本具备，实验对比和统计闭环尚未完成。**

需求没有规定各项权重，因此不以未经定义的百分比作为验收结论。本报告区分源码存在、历史记录与本次实际验证；代码已存在不等于移植数值等价性已得到验证。

不能据当前仓库材料得出“C++ RPN-CG 优于现有 C++ 一阶算法”或“优于 MATLAB 二阶方法”的结论。

## 按需求逐项核对

| 需求 | 状态 | 证据与说明 |
|---|---|---|
| 完成 RPN-CG 的 C++ 版本 | 基本完成 | `ROPTLIB/Solvers/RPNCG.cpp` 和 `ROPTLIB/Solvers/RPNCG.h` 已实现外层迭代、半光滑 Newton-CG (`FindDir`)、约束投影 tCG (`tCG`)、二阶算子 (`ApplyBB`) 和步长/回溯逻辑。 |
| 将 MATLAB 思路翻译为 C++ | 基本完成 | `RPN-CG/SPCA/solvers/driver_rpncg.m` 与 `ROPTLIB/Solvers/RPNCG.cpp` 存在对应的近端、乘子、tCG、回缩和线搜索流程；C++ 中也有 MATLAB 列主序索引和活动集处理的移植代码。 |
| 代码格式与现有风格统一 | 部分完成 | 代码可接入 ROPTLIB 的 `Problem`、`Manifold`、`SolversNSM` 接口，且使用项目的 `Variable/Vector` 类型；但 `RPNCG.cpp` 仍有较多 `MODIFIED` 注释、重复空行、旧式缩进和局部裸指针接口，建议后续按项目风格清理。 |
| 新建 Overleaf 文档 | 本地文档完成，线上状态未验证 | `ROPTLIB/RPNCG - latex/main.tex` 已建立，且已有 `main.pdf`；`README.md` 说明了 `compile.sh` 的编译方式。仓库材料不能证明已创建或上传至 Overleaf 线上项目。 |
| 说明类的关系 | 已完成 | `main.tex` 第 23--41 行描述 `RPNCG`、`SolversNSM`、`Problem` 和 `Stiefel` 的继承/使用关系。 |
| 说明成员函数和成员变量 | 已完成 | `main.tex` 第 43--75 行列出了主要成员函数、算法参数、继承字段和 tCG 状态码，并指出 `RPNCGParams` 当前没有被 `Run()` 使用。 |
| 说明运行逻辑 | 已完成 | `main.tex` 第 77--106 行说明外层迭代、`FindDir`、活动集、tCG、回缩和回溯过程。 |
| 说明内存管理 | 已完成 | `main.tex` 第 108--111 行说明非拥有指针、值对象、局部临时对象、`BuildBbar` 的内存开销和计数边界。 |
| 运行 C++ 测试 | 已完成但证据不完整 | `test/mainRPNCG.cpp` 固定随机种子 2，调用 `testRPNCGSPCA()`；已有文档记录一次结果：147 次外层迭代、目标值约 `-38.50512`、耗时 `0.558956 s`。当前环境尝试运行仓库已有二进制时缺少 `libmkl_intel_lp64.so.2`，因此该结果不能在本环境重新复现。 |
| 与现有 C++ 一阶算法比较效果和效率 | 未完成 | 仓库中有 `IARPG`（`Solvers/IARPG.cpp/.h`）及 `test/TestStieSPCA.cpp` 的测试，但其默认问题规模、`lambda`、初值和停止设置与 RPN-CG 测试不同，也没有在同一入口输出可比较的完整统计表。 |
| 与 MATLAB 二阶方法比较效果和效率 | 未完成 | 仓库中确实存在 MATLAB RPN-CG 代码：`RPN-CG/SPCA/solvers/driver_rpncg.m`；同时存在 `driver_rpn_cgh.m` 等二阶相关驱动。但没有找到与 C++ 测试共享数据、共享初值、统一停止条件并输出对比结果的实验记录。 |
| 输出时间开销和主要计算次数 | 未完成 | C++ 代码累计了 `nf/ng/nR` 的部分基类计数，但 `FindDir` 中的近端/Jacobian 调用、`tCG` 中的目标函数和 `ApplyBB` 中的 Hessian-vector 调用没有完整计数，且 `nH` 没有递增。MATLAB 驱动主要输出迭代数、目标值、CPU 时间、稀疏度和残差，没有形成与 C++ 对齐的调用次数表。 |
| 将对比结果写入 Overleaf 文档 | 未完成 | `main.tex` 第 126--140 行明确把 C++ 一阶和 MATLAB 对比列为“未取得”，只有 RPN-CG C++ 单次样例，没有最终对比表或图。 |

## 已有 C++ 实现和文档的具体内容

### 1. 算法实现

目前实现的主要结构如下：

- `RPNCG` 继承 `SolversNSM`。
- `StieSPCA` 提供稀疏 PCA 的目标函数、欧氏梯度、欧氏 Hessian 作用、近端映射和广义 Jacobian。
- `Stiefel` 提供 Stiefel 流形初值和 polar 型回缩。
- `Run()` 负责外层迭代、活动集划分、方向组合、步长调整和回溯。
- `FindDir()` 负责乘子方程的半光滑 Newton-CG 求解。
- `tCG()` 负责约束子空间中的二阶修正。
- `ApplyBB()` 调用 `Problem::EucHessianEta()` 实现 Hessian-vector 作用。

现有 `test/testRPNCG.cpp` 的测试配置为：`n=800`、`m=50`、`p=5`、`lambda=0.8`、随机种子为 2、最大迭代数为 1000、容差为 `1e-6`。输入矩阵列经过去均值和归一化，流形回缩选择 `ChooseParamsSet4()`。

### 2. 已有文档

`main.tex` 已覆盖需求中的四类说明内容，内容不仅是接口摘要，还记录了以下实现事实：

- `Prob` 和 `Mani` 是非拥有型指针，调用者负责保证其生命周期。
- `Variable`/`Vector` 和 STL 容器负责多数临时数据的自动生命周期管理。
- `BuildBbar` 显式构造活动集矩阵及其 Gram 矩阵/伪逆，存在额外内存开销。
- `Run()` 使用方向范数停止条件，而不是父类的全部停止机制。
- 当前计数器不能代表完整的 Hessian、近端映射和 Jacobian 调用次数。

这些说明足以作为当前实现的技术文档基础，但仍应补充编译版本、硬件、BLAS 后端、重复实验方式和正式对比结果。

## 实验和复现检查

### C++ 编译

已检查 `ROPTLIB/CMakeLists.txt`：目标为 `testRPNCG`，依赖 BLAS/LAPACK，C++ 标准为 C++11。当前环境可以完成 CMake 配置并发现 BLAS/LAPACK，但完整编译在限定时间内未结束，未形成新的可执行文件。

本次构建使用独立目录 `/tmp/opencode/rpncg-completeness-build`，Release 构建在 120 秒后超时，随后按 `run.sh` 使用的 Debug 配置重试，在 600 秒后超时。两次均未观察到编译诊断错误，但不能据此判定完整构建成功。XeLaTeX 输出同样写入 `/tmp/opencode`。

仓库已有的 `build/bin/testRPNCG` 在当前环境启动失败，原因是运行时找不到：

```text
libmkl_intel_lp64.so.2
```

这属于已有二进制的运行时依赖问题，不能作为算法测试失败的证据。应使用当前机器可用的 BLAS/LAPACK 重新构建，或补齐 MKL 运行库路径后再运行。

### LaTeX 编译

使用 XeLaTeX 编译 `main.tex` 成功，生成 3 页 PDF。编译过程中有字体替代和少量 overfull/underfull box 警告，但没有 LaTeX 错误。建议后续调整表格和长句排版，并确认 Overleaf 环境安装了 `Noto Sans CJK SC` 字体。

### 已记录的单次 C++ 样例

| 方法 | 外层迭代 | 最终目标值 | 程序内耗时 | 备注 |
|---|---:|---:|---:|---|
| C++ RPN-CG | 147 | 约 `-3.850512e+01` | `0.558956 s` | 来自现有文档记录的单次运行 |

该样例只能证明测试曾经产生过合理的收敛输出，不能用于稳定的速度排名。文档没有记录完整的 `nf`、`ng`、`nR`、Hessian-vector、近端、Newton 和 CG 总调用次数。

### 文档与移植的待确认点

- `main.tex` 第 134、138 行称仓库未找到 MATLAB RPN-CG 驱动，这对整个仓库而言不准确：驱动位于 `RPN-CG/SPCA/solvers/driver_rpncg.m`。应将缺口描述改为“未完成匹配的跨语言实验”。
- C++ 将 tCG 迭代上限设为 `max(1,min(200,round(1.2*|id1|)))`，MATLAB 驱动直接使用 `round(1.2*|id1|)`；这会影响运行路径和计算量，应解释或统一。
- C++ 的 `Tolerance` 同时影响初始内层精度和外层停止，MATLAB 分别使用 `option.tol` 和 `option.stop`，比较时不能只设置一个相同数字。
- C++ `tCG()` 在循环尾更新 `inner_it`，提前退出可能少计已经执行的 CG 步。`Run()` 结束时只保存 `ComTime`，没有调用已定义的 `PrintFinalInfo()` 输出完整汇总。
- `RPNCGParams` 未接入实际参数读取，`max_inner_iter` 未初始化也未使用；头文件和实现的缩进、注释仍需整理。
- 当前测试是单个随机 SPCA 示例，没有结果断言、最终正交约束检查或与 MATLAB 逐项核对的记录；CMake 也没有注册 CTest 测试。现有文件不足以证明所有边界情况及内存安全均已验证。

## 尚需完成的工作

### A. 统一测试数据和初值

至少固定以下内容，并让 C++ 与 MATLAB 读取同一份文件：

1. 数据矩阵 `B`；
2. 初值 `X0`；
3. `n,m,p`、`lambda`；
4. 回缩方式；
5. 容差、最大迭代数和停止判据；
6. BLAS/MATLAB 线程数和硬件环境。

MATLAB 原始实验中常见的设置与 C++ 当前用例并不完全一致，例如 `table_rand.m` 使用多个规模、20 个随机种子和不同算法；不能直接把这些结果与当前 C++ 单次 `seed=2` 结果相减比较。

### B. 增加可比的 C++ 一阶基线

建议新建一个统一测试入口，同时运行：

- RPN-CG；
- IARPG 或 IRPG；
- 如需与 MATLAB 论文实验对应，也可加入 ManPG/ManPG-Ada 等对应方法。

每种方法使用完全相同的 `B`、`X0` 和配置，并重复至少 5 次，报告中位数和最好/最差值。

### C. 增加 MATLAB 二阶基线

仓库已有 `RPN-CG/SPCA/solvers/driver_rpncg.m` 以及 `driver_rpn_cgh.m` 等驱动，应先确认需求所称“MATLAB 二阶方法”的准确方法名，再选定一个基线。必须导出与 C++ 相同的数据和初值，而不是在 MATLAB 中重新随机生成。

### D. 补齐计算次数统计

建议在 `RPNCG` 和 `StieSPCA` 中增加明确的实验计数器，至少统计：

- 目标函数 `f/g` 调用次数；
- 梯度调用次数；
- Hessian-vector 调用次数；
- 近端映射调用次数；
- 广义 Jacobian 调用次数；
- Newton 外层迭代数和 Newton 内部 CG 次数；
- tCG 迭代数；
- 回缩次数和回溯次数；
- 总外层迭代数和总耗时。

计数器应在算法内部调用点递增，不能用外层迭代数估算。尤其要避免把 `nH=0` 解释为没有二阶计算。

### E. 将实验结果补入 `main.tex`

最终文档至少应增加一张表，包含：

| 方法 | 最终目标值 | 最终残差/方向范数 | 外层迭代 | 总时间 | 梯度次数 | Hessian 次数 | 近端次数 | Newton/tCG 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|

同时注明机器、编译器、BLAS/MKL、线程数、随机种子、重复次数和停止条件。只有这些信息齐全后，才可以对效果和效率作出可信结论。

## 最终判定

1. **代码部分：基本完成。** 已有可读的 RPN-CG C++ 实现和 ROPTLIB 测试入口。
2. **说明部分：基本完成。** `main.tex` 已覆盖类关系、接口/成员、运行逻辑和内存管理。
3. **C++ 测试：部分完成。** 有单次测试样例，但当前环境无法运行已有 MKL 二进制，且缺少稳定重复实验。
4. **C++ 一阶对比：未完成。** 有候选算法和旧测试代码，但没有同数据、同初值、同指标的结果。
5. **MATLAB 二阶对比：未完成。** MATLAB RPN-CG 及相关驱动存在，但没有完成可复现的跨语言对比。
6. **时间和主要计算次数：未完成。** 现有统计不足以支撑效率比较。
7. **Overleaf 最终结果：未完成。** 文档框架完成，实验对比表和结论仍需补写。

因此，按 `req.txt` 的完整要求，当前仓库应标记为：

> **RPN-CG C++ 实现与技术说明基本完成；实验对比、完整计算统计和最终 Overleaf 结果记录未完成。**
