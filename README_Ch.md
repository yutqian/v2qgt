# v2qgt_v5 使用说明

[English README](README.md)

v2qgt_v5 从 VASP 的 `WAVECAR` 读取波函数、能量及 k 点，使用 Kubo 求和计算 Berry 曲率、量子度规或双带量子几何张量（QGT）不等式。k 点数量与顺序取自 `WAVECAR`，因此改变 `KPOINTS` 后需要用对应的新 `WAVECAR`。

## `-kubo` 编号

| 参数 | 当前代码计算的量 |
| --- | --- |
| `-kubo 0`（默认） | 跳过 Kubo 计算。当前 v5 主程序没有启用旧的重叠法 Berry 曲率分支。 |
| `-kubo 1` | 在 Cartesian `x,y` 方向计算逐带量。用 `-qm 0` 输出 Berry 曲率，用 `-qm 1` 输出量子度规。选定多条带时还会输出总量。 |
| `-kubo 2`、`3`、`4` | 与 `-kubo 1` 进入同一个分支，没有独立公式。 |
| `-kubo 5` | 双带 QGT 不等式。导数方向是 Cartesian `kx,ky`。每个 k 点同时输出度规、无迹 Berry 曲率矩阵及不等式两侧。 |
| `-kubo 6` | 同样的双带 QGT 不等式。导数沿倒易基矢，由 `QGT_KDIR1`、`QGT_KDIR2` 选取，默认为 `1,2`。 |
| `-kubo 7` | 同样的双带 QGT 不等式。导数沿两个任意的倒易空间方向，由 `QGT_DVEC1`、`QGT_DVEC2` 指定。 |

当前源码还会把其他正整数（例如 `8`）送进与 `1` 相同的分支。使用 `1` 表示该计算即可。

## 命令行参数

能带编号从 1 开始，对应 `WAVECAR` 里的顺序。参数按 `选项 值` 成对传入。

| 参数 | `-kubo 1–4` | `-kubo 5–7` |
| --- | --- | --- |
| `-nn N` | 逐带循环由 `-is` 或 `-ii/-if` 选择第一条带，`-nn` 在此分支不参与选带。 | 目标双带子空间的第一条带。 |
| `-mm M` | `-abelian 2` 时指定第二条带。 | 目标双带子空间的第二条带。请给出与 `N` 不同的有效编号。 |
| `-qm 0/1` | `0`（默认）输出曲率，主要看中间的 `xy` 分量。 `1` 输出 `g_xx, g_xy, g_yy`。 | 同时计算度规和曲率，此开关不影响结果。 |
| `-abelian 1/2` | `1`（默认）为单带 Abelian。 `2` 为双带子空间中的 non-Abelian 逐带计算，需提供 `-mm`。 | 始终调用双带例程。设为 `2` 可避免提示该模式预期使用 non-Abelian 双带。 |
| `-Qnm 0/1` | 设计上只在 `-abelian 2` 的逐带输出中生效。 `0` 选能带空间对角元 `Q_nn`，`1` 选非对角元 `Q_nm`。见下方源码问题。 | 双带例程内部计算完整的 2×2 能带矩阵，此开关不影响结果。 |
| `-is N` | 只输出第 `N` 条带的逐带结果。 | 双带由 `-nn/-mm` 指定。 |
| `-ii N -if M` | 逐带输出 `N` 到 `M`（含端点）。区间包含多条带时另写总量。 | 双带由 `-nn/-mm` 指定。 |

`-kubo 1–4` 的总量来自单独的例程。它对选定区间内的目标带求和，求和中的中间态取自 `-if` 以上的能带。因此总量通常不等于把各逐带文件直接相加，也不采用逐带分支的 `-abelian/-mm/-Qnm` 设置。

当前源码中，`main.f90` 将 `Q_nm` 声明为整数，但 `parse.f90` 和 `kubo_berry_curvature_nonabelian.f90` 的对应形参被隐式声明为双精度实数。这会影响 `-Qnm` 的读入及传递。使用 `-kubo 1 -abelian 2` 的逐带结果前，应先统一三处的类型并验证结果。

| 通用参数 | 用途 |
| --- | --- |
| `-f FILE` | 输入波函数文件，默认 `WAVECAR`。 |
| `-o TAG` | 设置输出标识。常规 QGT 调用中，`-o TAG` 对应文件前缀 `BERRYCURV.TAG`，省略时为 `BERRYCURV_KUBO`。 |
| `-s 2` / `-s 1` | 按 spinor（例如 SOC，默认）/ 非 spinor 格式读取波函数，应与生成 `WAVECAR` 的设置一致。 |
| `-kp N` | 控制普通 Kubo 输出的布里渊区扩展，默认 `2`。模式 5–7 直接输出原始 k 点。 |

## 模式 5–7 的方向和外带选择

- 模式 6 中，`QGT_KDIR1`、`QGT_KDIR2` 分别取 `1`、`2` 或 `3`，对应倒易基矢 `b1`、`b2`、`b3`。两个编号应不同。
- 模式 7 中，每个方向写作带引号的三个系数，例如 `QGT_DVEC1="1 1 0"` 表示沿 `b1+b2` 求导。默认两个方向分别为 `"1 0 0"` 和 `"0 1 0"`。方向必须非零且不平行。代码保留向量长度，缩放向量会改变对应的导数及张量分量。导数方向和 `WAVECAR` 中 k 点的取样路径可独立设置。
- `QGT_EXT_MODE=all`（默认）使用双带之外的所有中间带。 `near` 只取 `QGT_EXT_BMIN` 到 `QGT_EXT_BMAX`（含端点）内的中间带，`remote` 只取区间外的中间带。三种设置始终排除 `-nn/-mm`。
- 代码用两条目标带的平均能量作为求和分母的参考能量，并跳过能量差小于或等于 `1e-6 eV` 的项。用 `min_gap_eV` 和 `nskip_gap` 检查双带子空间与所选外带之间的小能隙。

## 输出

- 模式 1–4 的逐带文件名包含 `.EIG-N`，普通 SOC 输出为 `.dat`，共线双自旋可带 `.UP.dat`、`.DN.dat` 后缀。数据行依次为三个 Cartesian k 坐标、三个结果分量（`xx,xy,yy`）、三个分数 k 坐标。文件头沿用 Berry 曲率的旧标签，读结果时按 `-qm` 判断这三个数的含义。文件头的 Chern number 在当前分支未实际验证。
- 模式 5–7 的后缀分别为 `.qgtineq2d_kxky_MODE.dat`、`.qgtineq2d_k1k2_MODE.dat`（数字按所选轴变化）、`.qgtineq2d_dvec_MODE.dat`。`MODE` 为 `all`、`near` 或 `remote`。同目录中重复扫描不同能带对或模式 7 的不同方向时，用不同的 `-o` 值区分文件。
- 模式 5–7 的每行包含 k 点编号、三个分数坐标、半迹度规 `g11,g22,g12`、`detg`、`lhs`、`omega_norm2`、`rhs`、`margin`、无迹曲率矩阵 `Omega_SU2` 四个元素的实部和虚部，以及 `min_gap_eV,nskip_gap`。
- 不等式采用 `Omega = i(Q12-Q21)`，`Omega_SU2 = Omega - Tr(Omega) I/2`。输出的 `omega_norm2 = Tr(Omega_SU2²)/2`，`lhs = sqrt(max(detg,0))`，`rhs = sqrt(max(omega_norm2,0))/2`，`margin = lhs - rhs`。按代码的判据，`margin ≥ 0` 表示该点满足不等式。

## 示例

```bash
# 第 10 条带的 Abelian Berry 曲率
./v2qgt.x -kubo 1 -is 10 -qm 0 -o bc10

# 第 10 条带的 Abelian 量子度规
./v2qgt.x -kubo 1 -is 10 -qm 1 -o metric10

# 第 10、11 条带在 Cartesian kx,ky 方向的 QGT 不等式
./v2qgt.x -kubo 5 -nn 10 -mm 11 -abelian 2 -o pair10_11_xy

# 沿 b2,b3，且只让 9–12 中的外带进入求和（即第 9、12 条带）
QGT_KDIR1=2 QGT_KDIR2=3 \
QGT_EXT_MODE=near QGT_EXT_BMIN=9 QGT_EXT_BMAX=12 \
./v2qgt.x -kubo 6 -nn 10 -mm 11 -abelian 2 -o pair10_11_k2k3

# 沿 b1+b2 和 b3
QGT_DVEC1="1 1 0" QGT_DVEC2="0 0 1" \
./v2qgt.x -kubo 7 -nn 10 -mm 11 -abelian 2 -o pair10_11_dvec
```

在 Linux 上用 Intel MPI 和 MKL 编译 MPI 版时，在 Makefile 的 `F90` 参数中加 `-DMPI_USE`，再执行 `make clean && make`。以下命令让 4 个 MPI 进程分摊 k 点，由 rank 0 写出模式 7 的结果。

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 4 ./v2qgt.x -kubo 7 -nn 10 -mm 11 -abelian 2
```
