# 数据获取与放置

本仓库不包含GPS定位数据和逐段ERA5匹配表（约209 MB，含逐点坐标）。

## 需要的本地文件

审计脚本读取阶段数据包 `GWFG_Spring_DayNight_Wind_StageReport_Package/` 中的以下文件：

| 文件 | 行数 | SHA-256（前16位） |
|---|---:|---|
| `02_Data/01_Frozen_Inputs/Phase5a_Frozen_Movement_Data.csv` | 17,784 | `f374784b32630375` |
| `02_Data/01_Frozen_Inputs/Phase5c_GEE_Full_Matched_Movement_Data.csv` | 17,784 | `a30f748c60dfae3f` |
| `02_Data/02_Wind_Reconstruction/Phase5d2A_Relative_Wind_Data.csv` | 17,784 | `436323eca5cd7107` |

完整哈希见 `manifest/manifest_sha256.csv`。

## 尚未取得、但下一阶段必需的文件

| 文件 | 用途 |
|---|---|
| `2015 Tsinghua waterfowl (Yangtze).csv`（原始GPS，SHA-256 `9aca38d3…a7d21`） | 重建全日连续轨迹、利用设备瞬时地速与GPS高度 |
| 4X辅助晨昏/昼夜边界轨迹 | 恢复被晨昏剔除截断的出发过程 |
| Phase 5a、5c（GEE提取）、5d-1、5d-2A、5d-2C处理脚本 | 完整复现流程 |

## 放置方式

可以把数据包放在任意位置，用 `--package` 指定路径；如果放在本目录下（`data/GWFG_Spring_DayNight_Wind_StageReport_Package/`），`.gitignore` 会自动将其排除在版本控制之外。
