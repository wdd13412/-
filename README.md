效率

## 云端编译说明（Fortran）
仓库已包含主程序依赖的 `.f90` 源码与 `mesh/OFairfoilMesh` 网格文件。

如果你希望尽量与 Ubuntu 22.04 本地环境（GNU Fortran 11.4.0）一致，建议安装 `gfortran-11`：

```bash
sudo apt update
sudo apt install -y gfortran-11 liblapack-dev libblas-dev liblapacke-dev
```

> 说明：`sudo apt upgrade -y` 不是必须步骤，通常可省略以避免在云端环境引入额外系统变更。

编译：

```bash
./run_cloud_compile.sh
```

编译并立即运行：

```bash
./run_cloud_compile.sh --run
```

指定编译器（例如强制用 11 版）：

```bash
FC=gfortran-11 ./run_cloud_compile.sh --run
```

该脚本会：
1. 检查 5 个 Fortran 源文件与网格文件是否齐全；
2. 自动优先使用 `gfortran-11`（若存在），否则回退到 `gfortran`；
3. 输出当前编译器版本；
4. 执行以下编译命令：

```bash
gfortran -o buflow_test_d TypesModule.f90 meshdeformationn_d.f90 BuFlow_test_d.f90 main_d.f90 run_parameter_d_AAA.f90 -llapack
```
