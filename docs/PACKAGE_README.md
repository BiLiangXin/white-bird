# 白额雁春迁昼夜移动与风场响应：阶段数据包

本数据包归档当前可访问的Phase 5a/5c冻结分析输入、Phase 5d-2A风向轴与相对风数据、Phase 5d-2B移动段尺度探索性模型、Phase 5d-2C持续长距离移动候选事件分析，以及阶段审计报告、图件和现有分析脚本。

## 重要解释边界

- GPS `Distance_km`是相邻定位点间净位移，不是真实飞行路径长度；`Speed_kmh`不是空速。
- 10 m ERA5风不是鸟实际飞行高度的风。
- 春迁路线轴是排除当前段端点后的回顾性季节参照，不是严格实时航向。
- 224条事件是严格冻结样本中的持续长距离移动候选事件，不是已确认的224次完整迁飞。
- 当前逐段模型存在时间自相关；风向结果属于探索性关联，不是因果证明。
- 本脚本仅进行文件整理、哈希记录和结构/行数核验，不重新下载ERA5、不重跑GEE、不重新拟合生态模型。

## 目录结构

- `01_Report/`：阶段研究进展报告（Markdown/Word，如提供）。
- `02_Data/`：冻结输入、风向重构、移动序列、事件尺度和阶段比较数据。
- `03_Analysis_Artifacts/`：Phase 5d-2B/2C的模型输出、诊断JSON和源文件哈希记录。
- `04_Figures/`：Phase 5d-2、2A、2B、2C配套图件。
- `05_Code/`：数据打包脚本和已有Phase 5d-2B分析脚本。
- `06_Audit_Reports/`：Phase 5d-1至5d-2C审计报告。
- `08_Source_Archives/`：本阶段原始分析子包ZIP，供核对原始交付结构。
- `07_Manifest/`：文件SHA-256清单和自动核验结果。

## 重建数据包

需要Python 3.10或更高版本；打包脚本仅使用Python标准库。

```bash
python build_gwfg_wind_package.py --input-dir ./source --output-dir ./GWFG_Spring_DayNight_Wind_StageReport_Package --report-md ./GWFG_Spring_DayNight_Wind_Stage_Progress_Report.md --report-docx ./GWFG_Spring_DayNight_Wind_Stage_Progress_Report.docx
```

若本地有原始GPS源文件，可增加：

```bash
--raw-gps "路径/2015 Tsinghua waterfowl (Yangtze).csv"
```

该原始文件的既有SHA-256应为：`9aca38d37f9bed3e02cbb67a7d0df963c2c38df9371aaf469ec48bee851a7d21`。如有4X辅助晨昏/昼夜边界轨迹原文件，可通过`--aux-4x`加入归档；当前阶段库中未取得该文件。

## 数据完整性

自动核验包括：Phase 5a、5c、5d-2A行数与Segment_ID唯一性；Phase 5a/5c/5d-2A Segment_ID集合一致性；GEE匹配状态；春迁昼夜段数；路线轴有效覆盖；≥30 km段数；Phase 5d-2C事件、case、control和匹配层数量；已知冻结文件SHA-256。

原始GPS源文件和4X辅助轨迹原件如未通过可选参数提供，将在核验报告中明确标记为未包含，不会伪装成已打包文件。
