## Forward Communication API

### 1. 转发通信
转发通信可以实现27个邻居进程的通信，只用和6个面相邻的邻居进程通信即可完成（下称为 “全三维通信”），大大地降低了通信的开销。
相关论文可以参考：
[1] Thompson A P, Aktulga H M, Berger R, et al. LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales [J]. Computer Physics Communications, 2022, 271 : 108171.
[2] Plimpton S. Fast Parallel Algorithms for Short–Range Molecular Dynamics [J]. 1995 : 42.

此外，由于 KMC 中，需要和7个邻居进程进行 sector 的通信，同样可以采用转发通信算法（下称为“sector 通信”），通过和三个面相邻的邻居进程通信以达成目的。

### 2. 转发通信 API
头文件：
```cpp
#include <comm/comm.hpp>
```

#### 2.1 全三维通信 API:
```cpp
template <typename T, bool F = false>
void neiSendReceive(Packer<T> *packer, const mpi_process processes, const MPI_Datatype data_type,
                    const _MPI_Rank (&neighbours_rank)[DIMENSION_SIZE][2]);
```

其中，模板参数：
- `T` 是发送或接收的数据类型
- `F` 表示是否进行反向通信（即正向通信是按X，Y，Z维度的顺序进行通信，反向通信是按 Z，Y，X 维度的顺序进行通信）。 

函数参数：
- `packer` 是一个 `comm::Packer` 对象，用于数据打包和解包；
- `processes` 是一个 `comm::mpi_process` 对象，表示当前 MPI 进程的通信域；
- `data_type` 是 `MPI_Datatype` 类型的变量；
- `neighbours_rank` 是中的邻居进程的 rank，6个面相邻的邻居进程按维度、方向进行组织。

#### 2.2 sector 通信 API:
```cpp
template <typename T, typename RT, bool F = false>
void singleSideForwardComm(RegionPacker<T, RT> *packer, const mpi_process processes,
                            const MPI_Datatype data_type,
                            const std::array<std::vector<comm::Region<RT>>, DIMENSION_SIZE> send_regions,
                            const std::array<std::vector<comm::Region<RT>>, DIMENSION_SIZE> recv_regions,
                            const std::array<unsigned int, DIMENSION_SIZE> send_ranks,
                            const std::array<unsigned int, DIMENSION_SIZE> recv_ranks);
```
singleSideForwardComm 中的 singleSide 的含义是指，在某个维度的通信时，只考虑一个方向的通信（例如，在做 X 维度的通信时，只考虑和左侧的邻居进程进行通信，不和右侧的邻居进程进行通信）。
这个适合于 KMC 中的 sector 通信，因为在 KMC 中，某个维度的通信时，只需要和一个方向的邻居进程进行通信即可完成 sector 的通信。

其中，模板参数：
- `T` 是发送或接收的数据类型；
- `RT` 是 Region 中的坐标类型；
- `F` 表示是否进行反向通信（即正向通信是按X，Y，Z维度的顺序进行通信，反向通信是按 Z，Y，X 维度的顺序进行通信）。

函数参数：
- `packer` 是一个 `comm::RegionPacker` 对象（和 `comm::Packer` 类似），用于数据打包和解包；
- `processes` 是一个 `comm::mpi_process` 对象，表示当前 MPI 进程的通信域；
- `data_type` 是 `MPI_Datatype` 类型的变量，用于MPI通信函数；
- `send_regions` 表示每个维度的发送区域，数组长度对应于三个维度；
- `recv_regions` 表示每个维度的接收区域，数组长度对应于三个维度；
- `send_ranks` 表示每个维度的发送邻居 MPI 进程的 rank；
- `recv_ranks` 表示每个维度的接收邻居 MPI 进程的 rank。

### 3. Packer 框架
libcomm 中有三种 Packer，分别是 `comm::Packer`，`comm::RegionPacker` 和 `comm::NativePacker`。
其中，`comm::NativePacker` 继承自 `comm::Packer`，这两个均可用于全三维通信，而 `comm::RegionPacker` 只能用于 sector 通信。

`comm::NativePacker` 是旧的Packer API，需有在打包之前自己手动计算数据的长度；而`comm::Packer`由于采用 std::vector 进行数据存储，无需提前计算数据长度。

为了实现自定义的数据打包和解包，开发者需要继承以上三个类之一，然后实现类中的如下 API。

对于 `comm::Packer<T>`，需要实现以下 API：
| API | 参数说明 | API 说明|
|--|--|--|
| `void initialize()` | / | 通信前的初始化 |
| `const std::size_t pack(std::vector<T> &data, const int dimension, const int direction)` | data：打包的数据存放的buffer；dimension和direction：分别是通信的维度和方向（下同）| 数据打包函数。数据追加放入data中，返回打包后数据的长度。 |
| `void unpack(const std::vector<T> &data, const int dimension, const int direction)` | data：解包的数据buffer，可以通过 data.size() 获取数据长度；dimension和direction：分别是通信的维度和方向 | 数据解包函数；|
| `done()` | / | 通信完成后的函数（例如可以进行数据清理和内存释放）|

对于 `comm::NativePacker<T>`，需要实现以下 API：
| API | 参数说明 | API 说明|
|--|--|--|
| `const std::size_t sendLength(const int dimension, const int direction)` | dimension和direction：分别是通信的维度和方向 | 数据打包前，计算数据长度并返回。单位是数据类型 T（不是字节）。 |
| `void onSend(T buffer[], const std::size_t send_len, const int dimension, const int direction)` | buffer: 已经申请好内存的buffer，打包的数据存储在此，其长度为 `send_len`； send_len: 要打包的数据的长度（一般和 sendLength的返回值一样），单位T； dimension和direction：分别是通信的维度和方向 | 进行数据打包，buffer 内存已经提前申请好。 |
| `void onReceive(const T buffer[], const std::size_t receive_len, const int dimension, const int direction)` | buffer: 已经接收完毕的数据的buffer； receive_len: 已经接收到的戴解包的数据的长度，单位T； dimension和direction：分别是通信的维度和方向 | 进行数据解包。 |
| `void onFinish()` | / | 通信完成后的函数（例如可以进行数据清理和内存释放）|

对于 `comm::RegionPacker<T>`，需要实现以下 API 和 `comm::NativePacker<T>` 基本类似，
只是其 sendLength、onSend、onReceive这三个 API会额外多一个 region 的参数，表示发送或者接收的区域。

### 4. 通信区域的预设

为了便于通信的区域的计算，libcomm 结合 Domain API，内置了用于上述通信的区域，包括：
- `fwCommLocalSendRegion`、`fwCommLocalRecvRegion`：可以用于全三维通信，获取对应的发送区域或者接收区域的范围。坐标单位可以用格点单位或者实际的坐标单位。
- `fwCommSectorSendRegion`、`fwCommSectorRecvRegion`：可以用于sector 通信，获取对应的发送区域或者接收区域的范围。坐标单位只能是格点。

弃用的API：
- `fwCommLocalRegion`, `fwCommLocalMeaRegion`。