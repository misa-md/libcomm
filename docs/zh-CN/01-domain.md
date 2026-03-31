## domain API

## 三维空间区域分解的概念
三维空间下，假设模拟的区域大小为 $X \times Y \times Z$ (其中 X,Y,Z为正整数，表示各个维度的格点的数量)，
并假设启动了P个MPI进程，进程可以从一维到三维笛卡尔拓扑映射，使得再x，y，z维度分别有 $P_x$, $P_y$， $P_z$ 个进程
（其中 $P_x$, $P_y$， $P_z$ 均为正整数，且$P_x \times P_y \times P_z = P$）。

## 1 domain 的构建
### 1.1 一般的 domain 的构建

模拟 domain 的构建过程，会依据传入的各个维度的格点数和MPI进程信息，进行：  
（1）将MPI进程进行三维笛卡尔拓扑映射，并记录各个维度的进程数、三维下的进程坐标、三维的邻居进程号等。  
（2）将格点进行三维划分（让每个进程上格点数尽量可能均等），并记录：每个进程上的负责的格点数、进程负责的格点再全局格点坐标下的起始范围、ghost区域大小、加上ghost区域后的格点起始范围等信息。

**需要注意：划分和邻居的索引，均是周期性边界条件下的**。

对应的头文件如下：
```cpp
#include <comm/domain/domain.h>
```

依据格点构建模拟的 domain 的调用代码：
```cpp
comm::Domain* domain =
        comm::Domain::Builder()
            .setComm(comm::mpi_process process, MPI_Comm *new_comm)
            .setPhaseSpace(const int64_t lattice_num[3])
            .setCutoffRadius(const double cutoff_radius_factor)
            .setLatticeConst(const double lattice_const)
            .setGhostSize(const unsigned int ghost_size)
            .setMPIMap3dSubDim(const int mpi_map_3d_sub_dim[3])
            .build()
```
可以看出，构建 domain 采用链式的调用进行构建。
各个函数接收的参数解释如下表：

|API | 参数类型 | 函数说明 | 参数说明| 是否必须 |
|--|--|--|--|--|
|`setComm`           |`(comm::mpi_process process, MPI_Comm *new_comm)` |设置MPI通信参数| （进程属性，进程划分后的新的通信域） | Y |
|`setPhaseSpace`     |`const int64_t lattice_num[3]` | 设置三个维度的格点数量，构建阶段依据此进行格点划分 |全局的、xyz三个维度的格点数| Y |
|`setCutoffRadius`   |`const double cutoff_radius`| 设置相互作业的截断半径。如果没有显示设置ghost size，则用截断半径向上取整作为ghost size| 截断半径大小 | Y |
|`setLatticeConst`   |`const double latticeConst`| 设置格点间常数（即：格点阵列之间的最小距离） | 格点阵列之间的最小距离（如材料中的晶格常数） | Y |
|`setGhostSize`      |`const unsigned int ghost_size` | 设置xyz三个维度的 ghost 区域的大小 | ghost 区域大小（单位：格点间常数。）| N |
|`setMPIMap3dSubDim` |`const int mpi_map_3d_sub_dim[3]` | 设置三维MPI进程映射的子维度（如计算节点内的若干进程如何进行三维映射） | 三维MPI进程映射的子维度 | N |

### 1.2 BccDomain 的构建
除了通用的 Domain 外，为了方便，还支持 BccDomain。此domain是为了支持 BCC 晶格而设计的，与一般性的 domain的区域是，格点的坐标、数量在x维度均乘以了2
（使用的时候，需要用 BccDomain 类中的 `dbx_*` 的变量）。

BccDomain 扩展自 Domain类，其构建方式和 Domain 类似。

### 1.3 ColoredDomain 的构建
此外，libcomm 还支持 ColoredDomain。ColoredDomain 主要是用于 sub-lattice 算法（并行 KMC 中会用到）中的8个 sector 的区域范围的标记。

ColoredDomain 扩展自 Domain类，其构建方式和 Domain 类似。

## 2. comm::Domain 类的成员
`comm::Domain` 类中包含若干成员，主要包括进程信息和格点信息，如下：

### 体系属性类
|成员| 类型 |说明 |
|--|--|--|
|`lattice_const`| `double` | 格点间常数（即：格点阵列之间的最小距离）|
|`cutoff_radius_factor`| `const double` | 相互作业的截断半径|
|`cut_lattice`| `comm::_type_lattice_size` | 截断半径转换为格点单位（向上取整） |
|`phase_space`| `std::array<uint64_t, 3>` | 全局的xyz三个维度的格点数 |
|`ghost_size`| `int [3]` | xyz三个维度的 ghost 区域的大小（单位：格点间常数。）|

### MPI 进程相关属性
|成员| 类型 |说明 |
|--|--| -- |
|`grid_coord`| `int [3]` | 当前 MPI 进程在三维拓扑划分后的进程坐标 |
|`grid_size`| `int [3]` | 当前 MPI 进程在三维拓扑划分后的每个方向的MPI进程数 |
|`rank_id_neighbours`| `int [3][2]` | 当前MPI进程的6个面相邻的邻居进程的进程 rank。二维数组的低维度表示方向（前/后,右/左,上/下），高维度表示xyz的维度。 |

注意：邻居进程符合周期性边界条件。

### 模拟box：实际长度和范围
|成员 | 类型 | 说明 |
|--|--| --|
|`meas_global_length`| `const double [3]` | 全局模拟 box在三个维度上的长度 |
|`meas_global_region`| `const comm::Region<double>` | 全局模拟 box的模拟区域的范围 |
|`meas_sub_box_region`| `const comm::Region<double>` | 本MPI进程对应的实际模拟区域的范围 |
|`meas_ghost_length`| `const double [3]` | 本进程的ghost区域在xyz三个维度上的宽度 |
|`meas_ghost_ext_region`| `const comm::Region<double>` | 本进程模拟区域叠加外层的ghost区域后形成扩展区域，该成员表示该扩展区域的范围 |

### 格点范围和数量
|成员| 类型 |说明 |
|--|--|--|
|`sub_box_lattice_size`| `comm::_type_lattice_size [3]` | 本MPI进程负责的区域所对应的格点数（xyz三个方向） |
|`ghost_extended_lattice_size`|  `comm::_type_lattice_size [3]` | 本MPI进程负责的区域叠加外层ghost区域后形成扩展区域，该成员表示该扩展区域在xyz三个维度上所对应的格点数 |
|`lattice_size_ghost`| `comm::_type_lattice_size [3]` | ghost 区域的范围 |
|`sub_box_lattice_region`| `comm::Region<_type_lattice_coord>` | 本MPI进程负责的区域所对应的格点范围 |
|`ghost_ext_lattice_region`| `comm::Region<_type_lattice_coord>` | 本MPI进程负责的区域叠加外层ghost区域后形成扩展区域，该成员表示该扩展区域的格点范围 |
|`local_sub_box_lattice_region`| `comm::Region<_type_lattice_coord>` | 本MPI进程本地负责的区域所对应的格点范围 |
|`local_ghost_ext_lattice_region`| `comm::Region<_type_lattice_coord>` | 本MPI进程本地负责的区域叠加外层ghost区域后形成扩展区域，该成员表示该扩展区域的格点范围 |

注意，`sub_box_lattice_region` 和 `ghost_ext_lattice_region` 是全局坐标系下的格点范围，而 `local_sub_box_lattice_region` 和 `local_ghost_ext_lattice_region` 是本地坐标系下的格点范围，如下图是本地坐标系的示意图：
```
| 0:ghost_lower     | sub_box_lower       | sub_box_upper  | ghost_upper
|-------------------|---------...---------|----------------|
```

其中，`sub_box_lower` 和 `sub_box_upper` 分别是 `local_sub_box_lattice_region` 的对应区域起始坐标和结束坐标，而 `ghost_lower` 和 `ghost_upper` 分别是 `local_ghost_ext_lattice_region` 对应区域的起始坐标和结束坐标。  
对于全局坐标系，也是类似的，只是是放到整个模拟Box中的格点而言的。
