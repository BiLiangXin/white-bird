#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build and verify the Greater White-fronted Goose spring day/night wind package.

This script packages existing frozen inputs, derived tables, reports, figures and
analysis scripts. It does NOT rerun GEE/ERA5 extraction or refit ecological models.
Only Python's standard library is required.

Example:
  python build_gwfg_wind_package.py \
      --input-dir ./source \
      --output-dir ./GWFG_Spring_DayNight_Wind_StageReport_Package \
      --report-md ./GWFG_Spring_DayNight_Wind_Stage_Progress_Report.md \
      --report-docx ./GWFG_Spring_DayNight_Wind_Stage_Progress_Report.docx

Optional original source files:
  --raw-gps /path/to/2015_Tsinghua_waterfowl_Yangtze.csv
  --aux-4x /path/to/QA_checked_4X_auxiliary_tracks.csv
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path, PurePosixPath
from datetime import datetime, timezone

EXPECTED_SHA256 = {
    "Phase5a_Frozen_Movement_Data.csv": "f374784b32630375ee8e79444577a2ffe2a9319e8dd6016d50bef4a049455812",
    "Phase5c_GEE_Full_Matched_Movement_Data.csv": "a30f748c60dfae3f206bead43053b9c333db4007ade520df127d224979de3490",
    "2015 Tsinghua waterfowl (Yangtze).csv": "9aca38d37f9bed3e02cbb67a7d0df963c2c38df9371aaf469ec48bee851a7d21",
}

FILE_MAP = {
    "Phase5a_Frozen_Movement_Data.csv": "02_Data/01_Frozen_Inputs/Phase5a_Frozen_Movement_Data.csv",
    "Phase5c_GEE_Full_Matched_Movement_Data.csv": "02_Data/01_Frozen_Inputs/Phase5c_GEE_Full_Matched_Movement_Data.csv",
    "Phase5d2A_Relative_Wind_Data.csv": "02_Data/02_Wind_Reconstruction/Phase5d2A_Relative_Wind_Data.csv",
    "Phase5d2A_Heading_Method_Audit.csv": "02_Data/02_Wind_Reconstruction/Phase5d2A_Heading_Method_Audit.csv",
    "Phase5d2A_Coverage_QC.csv": "02_Data/02_Wind_Reconstruction/Phase5d2A_Coverage_QC.csv",
    "Phase5d2A_Spring_Route_Axis_Audit.csv": "02_Data/02_Wind_Reconstruction/Phase5d2A_Spring_Route_Axis_Audit.csv",
    "Phase5d2B_Long_Displacement_Sequences.csv": "02_Data/03_Movement_Sequences/Phase5d2B_Long_Displacement_Sequences.csv",
    "Phase5d2_Model_Comparison.csv": "02_Data/05_Stage_Comparison/Phase5d2_Model_Comparison.csv",
    "Phase5d2_Robustness_Results.csv": "02_Data/05_Stage_Comparison/Phase5d2_Robustness_Results.csv",
    "Phase5d2_Response_Variable_Audit.csv": "02_Data/05_Stage_Comparison/Phase5d2_Response_Variable_Audit.csv",
}
REPORT_FILES = {
    "Phase5d1_Movement_State_Audit.md": "06_Audit_Reports/Phase5d1_Movement_State_Audit.md",
    "Phase5d2_Manual_Review_Report.md": "06_Audit_Reports/Phase5d2_Manual_Review_Report.md",
    "Phase5d2A_Manual_Audit_Report.md": "06_Audit_Reports/Phase5d2A_Manual_Audit_Report.md",
    "Phase5d2B_Manual_Model_Report.md": "06_Audit_Reports/Phase5d2B_Manual_Model_Report.md",
    "Phase5d2C_Manual_Audit_Report.md": "06_Audit_Reports/Phase5d2C_Manual_Audit_Report.md",
}
ARCHIVES = [
    "Phase5d2_Manual_Review_Package.zip",
    "Phase5d2A_Manual_Review_Package.zip",
    "Phase5d2B_Manual_Analysis.zip",
    "Phase5d2C_Manual_Review_Package.zip",
]

README = r"""# 白额雁春迁昼夜移动与风场响应：阶段数据包

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
"""

DATA_DICTIONARY = r"""# 数据字典（文件级与关键字段）

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
"""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def find_exact(root: Path, name: str) -> Path | None:
    matches = [p for p in root.rglob(name) if p.is_file()]
    if not matches:
        return None
    # Prefer a direct child or the shallowest matching path.
    matches.sort(key=lambda p: (len(p.relative_to(root).parts), str(p)))
    if len(matches) > 1:
        hashes = {sha256(p) for p in matches}
        if len(hashes) > 1:
            raise RuntimeError(f"Multiple non-identical files named {name}: {matches}")
    return matches[0]


def copy_named(input_root: Path, name: str, output_root: Path, relpath: str, required: bool = True) -> bool:
    src = find_exact(input_root, name)
    if src is None:
        if required:
            raise FileNotFoundError(f"Required file not found: {name}")
        return False
    dst = output_root / relpath
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def safe_extract_selected(zpath: Path, output_root: Path, archive_kind: str) -> None:
    with zipfile.ZipFile(zpath) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            name = PurePosixPath(info.filename).name
            suffix = Path(name).suffix.lower()
            if archive_kind == "5d2":
                if suffix == ".png":
                    rel = Path("04_Figures/Phase5d2") / name
                else:
                    continue
            elif archive_kind == "5d2a":
                if suffix == ".png":
                    rel = Path("04_Figures/Phase5d2A") / name
                else:
                    continue
            elif archive_kind == "5d2b":
                if suffix == ".png":
                    rel = Path("04_Figures/Phase5d2B") / name
                elif suffix == ".py":
                    rel = Path("05_Code/Phase5d2B") / name
                elif suffix in (".csv", ".json") and name not in ("Phase5d2B_Long_Displacement_Sequences.csv",):
                    rel = Path("03_Analysis_Artifacts/Phase5d2B") / name
                else:
                    continue
            elif archive_kind == "5d2c":
                if suffix == ".png":
                    rel = Path("04_Figures/Phase5d2C") / name
                elif suffix == ".csv":
                    rel = Path("02_Data/04_Event_Scale") / name
                elif suffix == ".json":
                    rel = Path("03_Analysis_Artifacts/Phase5d2C") / name
                else:
                    continue
            else:
                continue
            dst = output_root / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            # Avoid unsafe archive paths: only the basename is used above.
            with zf.open(info) as src, dst.open("wb") as out:
                shutil.copyfileobj(src, out)


def stream_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yield row


def count_and_ids(path: Path, id_col: str = "Segment_ID"):
    n = 0
    ids = set()
    duplicates = 0
    for row in stream_csv(path):
        n += 1
        if id_col in row:
            value = row.get(id_col, "")
            if value in ids:
                duplicates += 1
            ids.add(value)
    return n, ids, duplicates


def verify_package(root: Path, raw_gps_included: bool, aux_4x_included: bool) -> dict:
    checks = {}
    p5a = root / FILE_MAP["Phase5a_Frozen_Movement_Data.csv"]
    p5c = root / FILE_MAP["Phase5c_GEE_Full_Matched_Movement_Data.csv"]
    p5w = root / FILE_MAP["Phase5d2A_Relative_Wind_Data.csv"]
    for label, p in (("Phase5a", p5a), ("Phase5c", p5c), ("Phase5d2A", p5w)):
        n, ids, dup = count_and_ids(p)
        checks[label] = {"rows": n, "unique_Segment_ID": len(ids), "duplicate_Segment_ID_rows": dup, "sha256": sha256(p)}
        if n != 17784 or len(ids) != 17784 or dup != 0:
            raise RuntimeError(f"{label} row/key check failed: {checks[label]}")
    ids_a = count_and_ids(p5a)[1]
    ids_c = count_and_ids(p5c)[1]
    ids_w = count_and_ids(p5w)[1]
    checks["Segment_ID_set_equality"] = {"Phase5a_equals_Phase5c": ids_a == ids_c, "Phase5a_equals_Phase5d2A": ids_a == ids_w}
    if ids_a != ids_c or ids_a != ids_w:
        raise RuntimeError("Segment_ID sets differ across Phase 5a/5c/5d-2A")

    # Phase 5c match status check.
    matched = 0
    total = 0
    for row in stream_csv(p5c):
        total += 1
        if row.get("GEE_Match_Status") == "matched":
            matched += 1
    checks["GEE_match_status"] = {"total": total, "matched": matched, "missing_or_other": total - matched}
    if matched != 17784:
        raise RuntimeError("GEE matched-row check failed")

    spring_n = day_n = night_n = route_valid = long30 = long30_day = long30_night = 0
    spring_birds = set()
    route_birds = set()
    for row in stream_csv(p5w):
        if row.get("Migration_Stage") != "Spring":
            continue
        spring_n += 1
        bird = row.get("Bird_ID", "")
        spring_birds.add(bird)
        dn = row.get("Day_Night")
        if dn == "Day": day_n += 1
        elif dn == "Night": night_n += 1
        valid = str(row.get("Spring_Route_Axis_Valid", "")).strip().lower() in ("true", "1", "yes")
        if valid:
            route_valid += 1
            route_birds.add(bird)
            try: dist = float(row.get("Distance_km", ""))
            except (TypeError, ValueError): dist = float("nan")
            if dist >= 30:
                long30 += 1
                if dn == "Day": long30_day += 1
                elif dn == "Night": long30_night += 1
    checks["Spring_route_axis"] = {
        "spring_segments": spring_n, "spring_birds": len(spring_birds),
        "day_segments": day_n, "night_segments": night_n,
        "route_axis_valid_segments": route_valid, "route_axis_valid_birds": len(route_birds),
        "valid_axis_and_ge30km": long30, "ge30km_day": long30_day, "ge30km_night": long30_night,
    }
    expected = (10092, 15, 6617, 3475, 9779, 14, 540, 400, 140)
    observed = (spring_n, len(spring_birds), day_n, night_n, route_valid, len(route_birds), long30, long30_day, long30_night)
    if observed != expected:
        raise RuntimeError(f"Spring coverage check failed: expected {expected}, observed {observed}")

    expected_event_files = {
        "Phase5d2C_Sustained_Long_Movement_Events.csv": 224,
        "Phase5d2C_Event_Onset_Segments.csv": 218,
        "Phase5d2C_Risk_Set_Control_Segments.csv": 6976,
        "Phase5d2C_Matched_Stratum_Wind_Contrasts.csv": 92,
        "Phase5d2C_Spring_Segment_Continuity_Audit.csv": 10092,
    }
    checks["Phase5d2C_event_tables"] = {}
    for name, expected_n in expected_event_files.items():
        p = root / "02_Data/04_Event_Scale" / name
        if not p.exists():
            raise FileNotFoundError(f"Event table missing from archive: {name}")
        n = sum(1 for _ in stream_csv(p))
        checks["Phase5d2C_event_tables"][name] = n
        if n != expected_n:
            raise RuntimeError(f"{name} expected {expected_n} rows, found {n}")

    seq_path = root / FILE_MAP["Phase5d2B_Long_Displacement_Sequences.csv"]
    seq_n = sum(1 for _ in stream_csv(seq_path))
    checks["Phase5d2B_sequences"] = {"rows": seq_n}
    if seq_n != 261:
        raise RuntimeError(f"Expected 261 Phase5d2B sequence rows, found {seq_n}")

    # Check immutable source hashes.
    for name, expected_hash in EXPECTED_SHA256.items():
        p = root / FILE_MAP.get(name, "02_Data/00_Original_Source/" + name)
        if p.exists():
            actual = sha256(p)
            if actual != expected_hash:
                raise RuntimeError(f"SHA-256 mismatch for {name}: {actual}")
            checks.setdefault("known_SHA256", {})[name] = {"expected": expected_hash, "actual": actual, "passed": True}
        elif name == "2015 Tsinghua waterfowl (Yangtze).csv":
            checks.setdefault("known_SHA256", {})[name] = {"expected": expected_hash, "actual": None, "passed": False, "status": "not_included_in_current_accessible_files"}

    checks["optional_sources"] = {
        "original_raw_GPS_included": raw_gps_included,
        "4X_auxiliary_tracks_included": aux_4x_included,
    }
    checks["interpretation"] = {
        "wind_effect_status": "exploratory_association_only",
        "continuous_event_status": "candidate_events_not_confirmed_flight_truth",
        "script_scope": "package_and_verify_existing_files_only; no GEE extraction or ecological model refitting",
    }
    return checks


def write_manifest(root: Path) -> None:
    manifest_dir = root / "07_Manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest = manifest_dir / "manifest_sha256.csv"
    rows = []
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p == manifest:
            continue
        rows.append((p.relative_to(root).as_posix(), p.stat().st_size, sha256(p)))
    with manifest.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["relative_path", "size_bytes", "sha256"])
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description="Build and verify the GWFG spring day/night wind package")
    ap.add_argument("--input-dir", required=True, help="Folder containing the Phase 5 source files and ZIP archives")
    ap.add_argument("--output-dir", required=True, help="Package folder to create")
    ap.add_argument("--report-md", help="Stage report Markdown to include")
    ap.add_argument("--report-docx", help="Stage report DOCX to include")
    ap.add_argument("--raw-gps", help="Optional original GPS source CSV")
    ap.add_argument("--aux-4x", help="Optional QA-checked 4X auxiliary trajectory file")
    ap.add_argument("--zip", dest="zip_path", help="Optional output ZIP path")
    args = ap.parse_args()

    input_root = Path(args.input_dir).expanduser().resolve()
    output_root = Path(args.output_dir).expanduser().resolve()
    if not input_root.exists():
        raise FileNotFoundError(input_root)
    if output_root == input_root or input_root in output_root.parents:
        raise RuntimeError("Output directory must not be inside the input directory")
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True)

    # Core frozen and derived tables.
    for name, rel in FILE_MAP.items():
        copy_named(input_root, name, output_root, rel, required=True)

    # Audit reports.
    for name, rel in REPORT_FILES.items():
        copy_named(input_root, name, output_root, rel, required=True)

    # Stage report files.
    report_dir = output_root / "01_Report"
    report_dir.mkdir(parents=True, exist_ok=True)
    for arg_value, default_name in ((args.report_md, "GWFG_Spring_DayNight_Wind_Stage_Progress_Report.md"),
                                    (args.report_docx, "GWFG_Spring_DayNight_Wind_Stage_Progress_Report.docx")):
        if arg_value:
            src = Path(arg_value).expanduser().resolve()
            if not src.exists():
                raise FileNotFoundError(src)
            shutil.copy2(src, report_dir / default_name)

    # Optional original source files.
    raw_gps_included = False
    if args.raw_gps:
        raw = Path(args.raw_gps).expanduser().resolve()
        if not raw.exists(): raise FileNotFoundError(raw)
        actual = sha256(raw)
        expected = EXPECTED_SHA256["2015 Tsinghua waterfowl (Yangtze).csv"]
        if actual != expected:
            raise RuntimeError(f"Original GPS source hash mismatch: expected {expected}, got {actual}")
        dst = output_root / "02_Data/00_Original_Source/2015 Tsinghua waterfowl (Yangtze).csv"
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(raw, dst)
        raw_gps_included = True
    aux_4x_included = False
    if args.aux_4x:
        aux = Path(args.aux_4x).expanduser().resolve()
        if not aux.exists(): raise FileNotFoundError(aux)
        dst = output_root / "02_Data/00_Original_Source" / aux.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(aux, dst)
        aux_4x_included = True

    # Preserve original sub-package ZIPs and unpack relevant outputs.
    archive_specs = [
        ("Phase5d2_Manual_Review_Package.zip", "5d2"),
        ("Phase5d2A_Manual_Review_Package.zip", "5d2a"),
        ("Phase5d2B_Manual_Analysis.zip", "5d2b"),
        ("Phase5d2C_Manual_Review_Package.zip", "5d2c"),
    ]
    archive_dir = output_root / "08_Source_Archives"
    archive_dir.mkdir(parents=True, exist_ok=True)
    for name, kind in archive_specs:
        src = find_exact(input_root, name)
        if src is None:
            raise FileNotFoundError(f"Required archive missing: {name}")
        shutil.copy2(src, archive_dir / name)
        safe_extract_selected(src, output_root, kind)

    # Add README and data dictionary.
    (output_root / "00_README.md").write_text(README, encoding="utf-8")
    (output_root / "00_DATA_DICTIONARY.md").write_text(DATA_DICTIONARY, encoding="utf-8")

    # Copy this builder script into the package.
    this_script = Path(__file__).resolve()
    code_dst = output_root / "05_Code/build_gwfg_wind_package.py"
    code_dst.parent.mkdir(parents=True, exist_ok=True)
    if this_script != code_dst:
        shutil.copy2(this_script, code_dst)

    # Automated structural verification.
    (output_root / "07_Manifest").mkdir(parents=True, exist_ok=True)
    checks = verify_package(output_root, raw_gps_included, aux_4x_included)
    checks["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    (output_root / "07_Manifest/verification_report.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8")
    write_manifest(output_root)

    # Build a ZIP outside the package directory.
    zip_path = Path(args.zip_path).expanduser().resolve() if args.zip_path else output_root.with_suffix(".zip")
    if zip_path.exists(): zip_path.unlink()
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for p in sorted(output_root.rglob("*")):
            if p.is_file():
                zf.write(p, p.relative_to(output_root.parent).as_posix())
    print(json.dumps({"package_dir": str(output_root), "zip": str(zip_path), "checks": checks}, ensure_ascii=False, indent=2))
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
