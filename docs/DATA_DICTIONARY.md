# 数据字典（文件级与关键字段）

| 文件 | 行数/对象 | 内容与用途 |
|---|---:|---|
| `Phase5a_Frozen_Movement_Data.csv` | 17,784个GPS移动段 | 冻结的GPS段级基础数据；包含所有阶段、昼夜标签、起止时间/坐标、净位移和时长。 |
| `Phase5c_GEE_Full_Matched_Movement_Data.csv` | 17,784个GPS移动段 | 与冻结GPS段按Segment_ID匹配的GEE/ERA5环境数据；17,784段均为matched。 |
| `Phase5d2A_Relative_Wind_Data.csv` | 17,784个GPS移动段 | 方向候选、春迁路线轴LOO、起中终点U/V和相对风变量。 |
| `Phase5d2A_Heading_Method_Audit.csv` | 方向审计记录 | H1/H2/H3及春迁路线轴的覆盖率、阈值敏感性和阶段×昼夜覆盖。 |
| `Phase5d2A_Coverage_QC.csv` | QC汇总记录 | 样本覆盖、距离分布、个体与春迁年份覆盖等质量审计。 |
| `Phase5d2A_Spring_Route_Axis_Audit.csv` | 26个鸟×春迁年份 | 各鸟年路线的首末时间、净位移和方位。 |
| `Phase5d2B_Long_Displacement_Sequences.csv` | 261条序列 | ≥30 km长位移段的连续序列摘要；不是261次独立起飞。 |
| `Phase5d2C_Sustained_Long_Movement_Events.csv` | 224条候选事件 | 严格冻结样本中至少由两个连续≥10 km段组成的候选事件。 |
| `Phase5d2C_Event_Onset_Segments.csv` | 218条有效起始段 | 具有有效事件起始风场的候选事件首段。 |
| `Phase5d2C_Risk_Set_Control_Segments.csv` | 6,976条对照段 | 与事件起始段按鸟、春迁年份、月份和昼夜匹配的非事件风险集对照。 |
| `Phase5d2C_Matched_Stratum_Wind_Contrasts.csv` | 92个匹配层 | 匹配层内事件起始风与对照风的差异。 |
| `Phase5d2_Model_Comparison.csv` | 670行模型记录 | 越冬/春迁昼夜移动响应候选模型比较与诊断记录。 |
| `Phase5d2_Robustness_Results.csv` | 5,020行 | 阈值、Bootstrap、留一鸟、预测和稳健性结果。 |
| `Phase5d2_Response_Variable_Audit.csv` | 18,069行审计记录 | 原始GPS响应变量的公式复核、分布统计和诊断记录。 |

## 关键字段

| 字段 | 含义 | 解释注意事项 |
|---|---|---|
| `Segment_ID` | GPS移动段唯一编号 | 跨Phase 5a/5c/5d-2A的主连接键。 |
| `Bird_ID` / `Tag_ID` | 个体与设备标识 | 重复观测结构；模型必须考虑个体层级。 |
| `Start_Time_UTC` / `End_Time_UTC` | 移动段起止UTC | 保留原始时间，不用固定钟点替代太阳时。 |
| `Start_Longitude` / `Start_Latitude`、`End_Longitude` / `End_Latitude` | GPS起止坐标 | 两端点之间的地表净位移不等于飞行路径。 |
| `Distance_km` | 两次GPS定位间净位移 | 30 km/10 km/50 km均为操作性位移阈值，不是飞行真值。 |
| `Duration_h` | 两次定位间实际时长 | 不等于连续飞行时间。 |
| `Speed_kmh` | 净位移/观测时长 | 不等于鸟的空速。 |
| `Migration_Stage` / `Day_Night` | 冻结阶段与昼夜标签 | 沿用既有太阳高度、太阳日期和阶段规则。 |
| `Spring_Year` | 春迁年份 | 用于鸟×春迁年份分层与路线轴重构。 |
| `GEE_Start_U10_mps` / `GEE_Start_V10_mps` | GPS段起点10 m风U/V分量 | 近地风，不是飞行高度风。 |
| `Spring_Seasonal_Route_Axis_LOO_Deg` | 春迁季节性路线轴LOO方位 | 回顾性路线参照，不是实时局地方位。 |
| `SpringRoute_AlongTrackWind_mps` | 沿春迁路线轴风分量 | 正值为助推分量，负值为阻碍分量。 |
| `SpringRoute_Crosswind_mps` | 侧风强度 | 非负值；用于表征横向风分量。 |
| `Bout_ID` | 连续长位移序列编号 | 一个序列可能含多个GPS段，不是多个独立事件。 |

## 事件尺度字段提示

`Phase5d2C_Sustained_Long_Movement_Events.csv`每行代表一条候选连续事件；`N_segments`、事件首末UTC、累计段距和首末Segment_ID用于追溯。事件累计距离是连续GPS段净位移之和，不是真实飞行路径。`Phase5d2C_Event_Onset_Segments.csv`保留事件起始段的GPS与风场字段；事件风场暴露优先采用起始段起点10 m U/V。
