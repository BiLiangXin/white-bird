"""探索性检查：月光、迁徙段长度与风选择、飞行定位点高度。

仅用于筛选下一阶段研究方向，不作为确认性结论。复用独立审计脚本的数据读取与事件构建。

用法：
  python -I src/audit/explore_new_factors.py --package "<阶段数据包路径>" --out docs
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, get_body, get_sun
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt  # noqa: E402

from gwfg_independent_audit import (ACCENT, DAY, INK, INK2, MUTED, NIGHT, bird_boot,  # noqa: E402
                                    build_onsets, contrast, design_month_strata, load_spring)

# 研究方向优先级（主观评估：创新性、现有数据可行性，1–5）
DIRECTIONS = [
    ("D1 准备度×触发的出发机制", 4.0, 3.0),
    ("D2 航段长度/屏障与风选择", 3.5, 4.0),
    ("D3 夜间出发与晴空（天体线索）", 4.0, 3.5),
    ("D4 月光×阶段的夜间活动", 3.5, 4.5),
    ("D5 人类压力与夜行化", 4.5, 3.5),
    ("D6 昼夜飞行速度与高度", 3.0, 2.5),
    ("D7 越冬携带效应与个体重复性", 3.0, 2.0),
    ("D8 农田景观与晨昏通勤", 3.5, 2.5),
    ("D9 高纬白夜下的节律", 3.0, 1.5),
]

warnings.filterwarnings("ignore")


def moon_light(df: pd.DataFrame) -> pd.DataFrame:
    """段起点处月球高度与月面照亮比例。"""
    t = Time(df.t0.dt.tz_convert(None).to_numpy())
    loc = EarthLocation(lon=df.Start_Longitude.to_numpy() * u.deg, lat=df.Start_Latitude.to_numpy() * u.deg)
    moon, sun = get_body("moon", t), get_sun(t)
    alt = moon.transform_to(AltAz(obstime=t, location=loc)).alt.deg
    elong = sun.separation(moon).rad
    illum = (1 - np.cos(elong)) / 2
    return df.assign(moon_alt=alt, moon_illum=illum,
                     moonlight=np.where(alt > 0, illum, 0.0))


def probe_moon(sp: pd.DataFrame) -> dict:
    night = moon_light(sp[sp.Day_Night == "Night"].copy())
    night["move1"] = (night.Distance_km >= 1).astype(int)
    night["long10"] = (night.Distance_km >= 10).astype(int)
    out = {"night_segments": int(len(night)), "moon_up_share": float((night.moon_alt > 0).mean())}
    bands = pd.cut(night.moonlight, [-0.01, 0.1, 0.6, 1.0], labels=["暗夜（<0.1）", "中等", "亮夜（>0.6）"])
    out["move1_by_band"] = {str(k): {"rate": float(v["mean"]), "n": int(v["size"])}
                            for k, v in night.groupby(bands, observed=True).move1.agg(["mean", "size"]).iterrows()}
    for y in ["move1", "long10"]:
        m = smf.glm(f"{y} ~ moonlight + C(bird_year) + C(month)", data=night,
                    family=sm.families.Binomial()).fit(cov_type="cluster",
                                                       cov_kwds={"groups": pd.factorize(night.Bird_ID)[0]})
        ci = m.conf_int().loc["moonlight"]
        out[y] = {"beta_full_vs_dark": float(m.params["moonlight"]), "ci": [float(ci[0]), float(ci[1])],
                  "rate_dark": float(night.loc[night.moonlight < 0.1, y].mean()),
                  "rate_bright": float(night.loc[night.moonlight > 0.6, y].mean()),
                  "n_dark": int((night.moonlight < 0.1).sum()), "n_bright": int((night.moonlight > 0.6).sum())}
    # 鸟等权：亮夜与暗夜的≥1 km移动率差
    g = night.assign(bright=night.moonlight > 0.5).groupby(["Bird_ID", "bright"]).move1.agg(["mean", "size"]).unstack()
    g = g[(g[("size", True)] >= 20) & (g[("size", False)] >= 20)]
    est, lo, hi = bird_boot(g[("mean", True)] - g[("mean", False)])
    out["bird_equal_move1_bright_minus_dark"] = {"estimate": est, "ci": [lo, hi], "birds": int(len(g))}
    return out


def probe_leg_length(sp: pd.DataFrame) -> dict:
    """离开东北停歇区（起点44–50°N、事件终点>50°N）的北上事件 vs 其他事件。"""
    ev = build_onsets(sp, 10)
    ev = ev[ev.Spring_Route_Axis_Valid == True]
    # 事件终点纬度：同一连续串最后一段的终点
    long_ = sp.Distance_km >= 10
    starts = long_ & ~(sp.cont_prev & long_.shift(1, fill_value=False))
    run = starts.cumsum().where(long_)
    end_lat = sp.groupby(run).End_Latitude.last()
    run_of = run.loc[ev.index]
    ev = ev.assign(end_lat=run_of.map(end_lat).to_numpy())
    ev["leave_ne"] = ev.Start_Latitude.between(44, 50) & (ev.end_lat > 50)
    valid = sp[sp.Spring_Route_Axis_Valid == True]
    pool = valid[valid.Distance_km < 10]
    out = {}
    for name, sub in [("离开东北停歇区北上", ev[ev.leave_ne]), ("其他持续事件", ev[~ev.leave_ne])]:
        a, c = design_month_strata(sub, pool)
        r = contrast(a, c, "along")
        out[name] = {"events": int(len(a)), "birds": r["birds"], "along_diff": r["estimate"], "ci": r["ci"],
                     "night_share": float((a.Day_Night == "Night").mean())}
    return out


def probe_flight_fixes(pkg: Path) -> dict:
    cols = ["Migration_Stage", "Day_Night", "Bird_ID", "End_Device_Ground_Speed_mps",
            "End_Height_Above_Ellipsoid_m", "End_Latitude"]
    a = pd.read_csv(pkg / "02_Data/01_Frozen_Inputs/Phase5a_Frozen_Movement_Data.csv", usecols=cols, low_memory=False)
    sp = a[a.Migration_Stage == "Spring"]
    fly = sp[sp.End_Device_Ground_Speed_mps >= 5]
    return {
        "speed_available_share": float(sp.End_Device_Ground_Speed_mps.notna().mean()),
        "birds_with_speed": int(sp.loc[sp.End_Device_Ground_Speed_mps.notna(), "Bird_ID"].nunique()),
        "flight_fixes": int(len(fly)), "flight_birds": int(fly.Bird_ID.nunique()),
        "flight_fixes_by_daynight": fly.Day_Night.value_counts().to_dict(),
        "ellipsoid_height_median_flight": fly.groupby("Day_Night").End_Height_Above_Ellipsoid_m.median().to_dict(),
        "flight_ground_speed_median_mps": fly.groupby("Day_Night").End_Device_Ground_Speed_mps.median().to_dict(),
    }


def plot_signals(res: dict, figdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax = axes[0]
    bands = res["moon"]["move1_by_band"]
    x = np.arange(len(bands))
    ax.bar(x, [v["rate"] * 100 for v in bands.values()], width=0.6, color=NIGHT)
    for xi, v in zip(x, bands.values()):
        ax.text(xi, v["rate"] * 100 + 0.4, f"{v['rate']:.1%}\nn={v['n']:,}", ha="center", color=INK2, fontsize=9)
    ax.set_xticks(x, list(bands))
    ax.set_ylim(0, 24)
    ax.set_ylabel("夜间≥1 km移动段比例（%）")
    ax.set_title("(a) 亮夜的夜间局地移动更多（未校正）", loc="left", color=INK)
    ax.grid(axis="x", visible=False)

    ax = axes[1]
    rows = list(res["leg_length"].items())
    y = np.arange(len(rows))[::-1]
    for yi, (name, r), col in zip(y, rows, [ACCENT, DAY]):
        ax.hlines(yi, r["ci"][0], r["ci"][1], color=col, linewidth=2)
        ax.plot(r["along_diff"], yi, "o", color=col, markersize=8, markeredgecolor="white", markeredgewidth=1.5)
        ax.text(r["ci"][1] + 0.08, yi, f"{r['events']}个事件/{r['birds']}只鸟", va="center", color=INK2, fontsize=9)
    ax.axvline(0, color=INK2, linewidth=1, linestyle="--")
    ax.set_yticks(y, [n for n, _ in rows])
    ax.set_xlim(-0.3, 3.2)
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("事件起始 − 对照的沿轴风（m/s，鸟等权95%区间）")
    ax.set_title("(b) 跨越泰加林的大跳跃更挑风", loc="left", color=INK)
    ax.grid(axis="y", visible=False)
    fig.suptitle("图B2  气象以外因素的初步信号（探索性，未经确认性检验）", x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "B2_preliminary_signals.png")
    plt.close(fig)


def plot_priority(figdir: Path) -> None:
    fig, (ax, side) = plt.subplots(1, 2, figsize=(10.5, 5.2), gridspec_kw={"width_ratios": [1.5, 1]})
    ax.axvspan(2.9, 5, ymin=(3.4 - 2) / 3, ymax=1, color=DAY, alpha=0.06, linewidth=0)
    ax.text(4.95, 4.92, "优先开展区", ha="right", va="top", color=DAY, fontsize=10)
    for name, innov, feas in DIRECTIONS:
        hi = innov >= 3.5 and feas >= 3.0
        ax.scatter(feas, innov, s=110, color=DAY if hi else MUTED, edgecolors="white", linewidths=1.5, zorder=3)
        ax.annotate(name[:2], (feas, innov), xytext=(0, 9), textcoords="offset points", ha="center",
                    fontsize=9, color=INK if hi else INK2)
    ax.set_xlim(1, 5)
    ax.set_ylim(2, 5)
    side.axis("off")
    for i, (name, innov, feas) in enumerate(DIRECTIONS):
        hi = innov >= 3.5 and feas >= 3.0
        side.text(0, 1 - i * 0.105, name, transform=side.transAxes, va="top", fontsize=10,
                  color=INK if hi else INK2, fontweight="bold" if hi else "normal")
    ax.set_xlabel("现有数据可行性（1–5）")
    ax.set_ylabel("创新性（1–5）")
    fig.suptitle("图B1  候选研究方向的优先级（主观评估）", x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "B1_direction_priority.png")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--package", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    sp = load_spring(args.package)
    res = {"moon": probe_moon(sp), "leg_length": probe_leg_length(sp), "flight_fixes": probe_flight_fixes(args.package)}
    (args.out / "explore_new_factors.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    figdir = args.out / "figures"
    plot_signals(res, figdir)
    plot_priority(figdir)
    print(json.dumps(res["moon"]["move1_by_band"], ensure_ascii=False))


if __name__ == "__main__":
    main()
