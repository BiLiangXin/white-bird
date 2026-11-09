"""白额雁春迁昼夜移动—风场项目：独立数据审计与事件起始风场再分析。

只读取阶段数据包中的冻结数据，不修改任何输入文件。输出：
  - <out>/audit_results.json  全部审计数值
  - <out>/figures/*.png       中文图件

用法：
  python -I src/audit/gwfg_independent_audit.py \
      --package "<阶段数据包路径>" \
      --out docs
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.discrete.conditional_models import ConditionalLogit

warnings.filterwarnings("ignore")

# 颜色：固定顺序的分类色（白天=槽位1，夜间=槽位2），文字与网格使用中性色
DAY, NIGHT = "#2a78d6", "#eb6834"
ACCENT, MUTED = "#6250d6", "#a3a29c"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "font.size": 10,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.axisbelow": True,
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})

N_BOOT = 5000
SEED = 20261109


# ---------------------------------------------------------------- 数据读取
def load_spring(pkg: Path) -> pd.DataFrame:
    w = pd.read_csv(pkg / "02_Data/02_Wind_Reconstruction/Phase5d2A_Relative_Wind_Data.csv", low_memory=False)
    c = pd.read_csv(pkg / "02_Data/01_Frozen_Inputs/Phase5c_GEE_Full_Matched_Movement_Data.csv",
                    low_memory=False, usecols=["Segment_ID", "GEE_Start_Temperature2m_C"])
    assert len(w) == 17784 and w.Segment_ID.is_unique
    sp = w[w.Migration_Stage == "Spring"].merge(c, on="Segment_ID", how="left").copy()
    assert len(sp) == 10092
    sp["t0"] = pd.to_datetime(sp.Start_Time_UTC, utc=True)
    sp["t1"] = pd.to_datetime(sp.End_Time_UTC, utc=True)
    sp = sp.sort_values(["Bird_ID", "Tag_ID", "t0"]).reset_index(drop=True)

    # 起点10 m风投影到季节性路线轴：有符号沿轴风（正=助推）与侧风强度
    th = np.deg2rad(sp.Spring_Seasonal_Route_Axis_LOO_Deg)
    sp["along"] = sp.GEE_Start_U10_mps * np.sin(th) + sp.GEE_Start_V10_mps * np.cos(th)
    sp["cross"] = np.abs(sp.GEE_Start_U10_mps * np.cos(th) - sp.GEE_Start_V10_mps * np.sin(th))
    u_mean = (sp.GEE_Start_U10_mps + 4 * sp.GEE_Mid_U10_mps + sp.GEE_End_U10_mps) / 6
    v_mean = (sp.GEE_Start_V10_mps + 4 * sp.GEE_Mid_V10_mps + sp.GEE_End_V10_mps) / 6
    sp["along_simpson"] = u_mean * np.sin(th) + v_mean * np.cos(th)
    sp["ws"] = sp.GEE_Start_Wind_Speed10_mps
    sp["T"] = sp.GEE_Start_Temperature2m_C
    sp["bird_year"] = sp.Bird_ID + "_" + sp.Spring_Year.astype(int).astype(str)
    sp["night"] = (sp.Day_Night == "Night").astype(int)
    sp["month"] = sp.t0.dt.month
    sp["solar_hour"] = np.floor(sp.Start_Apparent_Solar_Hour).astype(int) % 24

    # 与Phase 5d-2C相同的严格连续性规则
    g = sp.groupby(["Bird_ID", "Tag_ID"])
    cont = ((sp.Start_Event_ID == g.End_Event_ID.shift(1))
            & (sp.t0 == g.t1.shift(1))
            & np.isclose(sp.Start_Longitude, g.End_Longitude.shift(1), atol=1e-8, rtol=0)
            & np.isclose(sp.Start_Latitude, g.End_Latitude.shift(1), atol=1e-8, rtol=0)
            & (sp.Spring_Year == g.Spring_Year.shift(1)))
    sp["cont_prev"] = cont.fillna(False).astype(bool)
    sp["cont_next"] = sp.groupby(["Bird_ID", "Tag_ID"]).cont_prev.shift(-1).fillna(False).astype(bool)
    return sp


def bird_boot(values: pd.Series, seed: int = SEED) -> tuple[float, float, float]:
    """鸟等权均值及整鸟重抽样95%区间。"""
    arr = values.dropna().to_numpy()
    rng = np.random.default_rng(seed)
    bs = rng.choice(arr, size=(N_BOOT, len(arr)), replace=True).mean(axis=1)
    return float(arr.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


# ---------------------------------------------------------------- 1 样本结构
def audit_composition(sp: pd.DataFrame, figdir: Path) -> dict:
    bins = [30, 40, 45, 50, 55, 60, 65]
    labels = ["30–40°N", "40–45°N", "45–50°N", "50–55°N", "55–60°N", "60–65°N"]
    sp = sp.assign(lat_band=pd.cut(sp.Start_Latitude, bins, labels=labels))
    tab = sp.groupby(["lat_band", "Day_Night"], observed=True).size().unstack(fill_value=0)
    tab["night_share"] = tab.Night / (tab.Day + tab.Night)
    tab["birds"] = sp.groupby("lat_band", observed=True).Bird_ID.nunique()
    tab["median_doy"] = sp.groupby("lat_band", observed=True).t0.apply(lambda x: float(x.dt.dayofyear.median()))

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9), gridspec_kw={"width_ratios": [1.25, 1]})
    y = np.arange(len(tab))
    ax = axes[0]
    ax.barh(y, tab.Day, color=DAY, height=0.62, label="白天")
    ax.barh(y, tab.Night, left=tab.Day, color=NIGHT, height=0.62, label="夜间",
            edgecolor="white", linewidth=1.5)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("GPS移动段数")
    ax.set_title("(a) 春迁移动段的纬度分布", loc="left", color=INK)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    ax.hlines(y, 0, tab.night_share * 100, color=MUTED, linewidth=1.2)
    ax.plot(tab.night_share * 100, y, "o", color=NIGHT, markersize=8)
    for yi, (share, n) in enumerate(zip(tab.night_share, tab.Day + tab.Night)):
        ax.text(share * 100 + 2, yi, f"{share:.0%}（n={n:,}）", va="center", color=INK2, fontsize=9)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 62)
    ax.set_xlabel("夜间段占比（%）")
    ax.set_title("(b) 夜间样本随纬度迅速消失", loc="left", color=INK)
    ax.grid(axis="y", visible=False)
    fig.suptitle("图A1  春迁“昼夜”比较与纬度/季节进程混杂", x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "A1_latitude_daynight_confound.png")
    plt.close(fig)

    return {
        "spring_segments": int(len(sp)),
        "spring_birds": int(sp.Bird_ID.nunique()),
        "bird_years": int(sp.bird_year.nunique()),
        "by_latitude_band": {k: {kk: (float(vv) if isinstance(vv, (float, np.floating)) else int(vv))
                                 for kk, vv in row.items()} for k, row in tab.to_dict("index").items()},
        "segments_per_bird_year": sp.groupby("bird_year").size().describe().round(1).to_dict(),
    }


# ---------------------------------------------------------------- 2 路线轴与风变量
def audit_axis(sp: pd.DataFrame) -> dict:
    v = sp[sp.Spring_Route_Axis_Valid == True]
    uniq = v.groupby("bird_year").Spring_Seasonal_Route_Axis_LOO_Deg.nunique()
    axis_signed = v.Spring_Seasonal_Route_Axis_LOO_Deg.where(v.Spring_Seasonal_Route_Axis_LOO_Deg < 180,
                                                             v.Spring_Seasonal_Route_Axis_LOO_Deg - 360)
    stored = v.SpringRoute_AlongTrackWind_mps
    L = v[v.Distance_km >= 30]
    within = v.assign(a=v.along - v.groupby("bird_year").along.transform("mean"),
                      t=v["T"] - v.groupby("bird_year")["T"].transform("mean"))
    return {
        "loo_axis_unique_values_per_bird_year_max": int(uniq.max()),
        "axis_deg_range": [float(axis_signed.min()), float(axis_signed.max())],
        "corr_along_vs_V10": float(np.corrcoef(v.along, v.GEE_Start_V10_mps)[0, 1]),
        "corr_within_along_vs_T2m": float(np.corrcoef(within.a, within.t)[0, 1]),
        "stored_SpringRoute_AlongTrackWind_equals": (
            "simpson_mean" if np.allclose(stored, v.along_simpson, atol=1e-9) else "start"),
        "ge30_sign_disagreement_start_vs_simpson": int((np.sign(L.along) != np.sign(L.along_simpson)).sum()),
        "ge30_n": int(len(L)),
    }


# ---------------------------------------------------------------- 3 晨昏剔除导致的连续性断点
def audit_continuity(sp: pd.DataFrame, figdir: Path) -> dict:
    sp = sp.assign(L10=sp.Distance_km >= 10)
    t = sp.groupby("solar_hour").agg(n=("L10", "size"), L10=("L10", "sum"), next_ok=("cont_next", "mean"))
    t["L10_next_ok"] = sp[sp.L10].groupby("solar_hour").cont_next.mean()
    even = t[t.index % 2 == 0]

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 5.6), sharex=True, gridspec_kw={"height_ratios": [1.1, 1]})
    hot = [16, 2]
    for ax in axes:
        ax.axvspan(17, 19, color=MUTED, alpha=0.15, linewidth=0)
    axes[0].text(18, even.L10.max() * 0.92, f"晨昏窗口\n多数段被剔除\n（18时仅{int(even.n[18])}段）",
                 ha="center", va="top", color=INK2, fontsize=8.5)
    colors = [ACCENT if h in hot else DAY for h in even.index]
    axes[0].bar(even.index, even.L10, width=1.4, color=colors)
    for h in hot:
        axes[0].text(h, even.L10[h] + 6, f"{int(even.L10[h])}", ha="center", color=INK, fontsize=9)
    axes[0].set_ylabel("≥10 km段数")
    axes[0].set_title("(a) 各太阳时起点的≥10 km净位移段", loc="left", color=INK)
    axes[0].grid(axis="x", visible=False)

    axes[1].plot(even.index, even.L10_next_ok * 100, "-o", color=DAY, linewidth=2, markersize=8)
    for h in hot:
        axes[1].plot(h, even.L10_next_ok[h] * 100, "o", color=ACCENT, markersize=9)
        axes[1].annotate(f"{even.L10_next_ok[h]:.0%}", (h, even.L10_next_ok[h] * 100), xytext=(8, 4),
                         textcoords="offset points", color=INK, fontsize=9)
    axes[1].set_ylabel("后继段连续（%）")
    axes[1].set_ylim(-5, 105)
    axes[1].set_xlabel("段起点当地视太阳时（h）")
    axes[1].set_xticks(range(0, 24, 2))
    axes[1].set_title("(b) 该≥10 km段之后存在严格连续后继段的比例", loc="left", color=INK)
    fig.suptitle("图A2  晨昏段被剔除后，16时（黄昏前）与2时出发的位移链被系统性截断",
                 x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "A2_twilight_truncation.png")
    plt.close(fig)

    return {
        "by_solar_hour": {int(h): {"n": int(r.n), "L10": int(r.L10), "next_ok": round(float(r.next_ok), 3),
                                   "L10_next_ok": (None if pd.isna(r.L10_next_ok) else round(float(r.L10_next_ok), 3))}
                          for h, r in t.iterrows()},
        "share_of_L10_starting_16h": float(t.L10.get(16, 0) / t.L10.sum()),
    }


# ---------------------------------------------------------------- 4 昼夜效应的纬度敏感性
def audit_latitude_sensitivity(sp: pd.DataFrame) -> dict:
    v = sp[sp.Spring_Route_Axis_Valid == True].copy()
    v["y"] = (v.Distance_km >= 30).astype(int)
    for x in ["along", "cross"]:
        v[x + "_w"] = v[x] - v.groupby("bird_year")[x].transform("mean")
    form = "y ~ night + along_w + cross_w + night:along_w + C(bird_year) + C(month)"
    subsets = {
        "全部有效路线轴段": v,
        "仅<50°N": v[v.Start_Latitude < 50],
        "仅<50°N且≤5月10日": v[(v.Start_Latitude < 50) & (v.t0.dt.dayofyear <= 130)],
    }
    out = {}
    for name, d in subsets.items():
        m = smf.glm(form, data=d, family=sm.families.Binomial()).fit(
            cov_type="cluster", cov_kwds={"groups": pd.factorize(d.Bird_ID)[0]})
        ci = m.conf_int()
        out[name] = {"n": int(len(d)), "events": int(d.y.sum()), "birds": int(d.Bird_ID.nunique()),
                     **{term: [float(m.params[term]), float(ci.loc[term, 0]), float(ci.loc[term, 1])]
                        for term in ["night", "along_w", "night:along_w"]}}
    return out


# ---------------------------------------------------------------- 5 事件起始风场再分析
def build_onsets(sp: pd.DataFrame, thr: float) -> pd.DataFrame:
    long_ = sp.Distance_km >= thr
    starts = long_ & ~(sp.cont_prev & long_.shift(1, fill_value=False))
    run_id = starts.cumsum().where(long_)
    run_len = run_id.map(run_id.value_counts())
    return sp[starts & (run_len >= 2)].copy()


def design_month_strata(ev, pool):
    key = lambda d: d.Bird_ID + "|" + d.Spring_Year.astype(str) + "|" + d.month.astype(str) + "|" + d.Day_Night
    ev, pool = ev.assign(stratum=key(ev)), pool.assign(stratum=key(pool))
    pool = pool[pool.stratum.isin(ev.stratum)]
    return ev[ev.stratum.isin(pool.stratum)], pool


def design_time_local(ev, pool, days=3, hours=1.0):
    """同鸟同年、±days天、当地太阳时±hours小时的时间局地风险集。"""
    rows = []
    for idx, e in ev.iterrows():
        c = pool[(pool.bird_year == e.bird_year) & ((pool.t0 - e.t0).abs() <= pd.Timedelta(days=days))]
        dh = np.abs(((c.Start_Apparent_Solar_Hour - e.Start_Apparent_Solar_Hour) + 12) % 24 - 12)
        c = c[dh <= hours]
        if len(c):
            rows.append(c.assign(stratum=f"e{idx}"))
    ctrl = pd.concat(rows)
    ev = ev.assign(stratum=[f"e{i}" for i in ev.index])
    return ev[ev.stratum.isin(ctrl.stratum)], ctrl


def contrast(ev, ctrl, var):
    s = (pd.concat([ev.assign(y=1), ctrl.assign(y=0)])
         .groupby(["stratum", "Bird_ID", "y"])[var].mean().unstack("y").dropna())
    per_bird = (s[1] - s[0]).groupby("Bird_ID").mean()
    est, lo, hi = bird_boot(per_bird)
    return {"estimate": est, "ci": [lo, hi], "birds": int(len(per_bird)),
            "birds_positive": int((per_bird > 0).sum()), "per_bird": per_bird.round(3).to_dict()}


def clogit(ev, ctrl, cols):
    d = pd.concat([ev.assign(y=1), ctrl.assign(y=0)]).dropna(subset=cols)
    m = ConditionalLogit(d.y, d[cols], groups=d.stratum).fit(disp=0)
    return {c: [float(m.params[c]), float(m.bse[c])] for c in cols}


def audit_events(sp: pd.DataFrame, figdir: Path) -> dict:
    valid = sp[sp.Spring_Route_Axis_Valid == True]
    pool = valid[valid.Distance_km < 10]
    out = {}
    for thr in [10, 30, 50]:
        ev_all = build_onsets(sp, thr)
        ev = ev_all[ev_all.Spring_Route_Axis_Valid == True]
        evA, cA = design_month_strata(ev, pool)
        evB, cB = design_time_local(ev, pool)
        res = {"events_ge2_segments": int(len(ev_all)), "events_valid_axis": int(len(ev)),
               "A": {"events": int(len(evA)), "controls": int(len(cA)), "strata": int(evA.stratum.nunique())},
               "B": {"events": int(len(evB)), "controls": int(len(cB))}}
        for var in ["along", "cross", "ws", "T"]:
            res["A"][var] = contrast(evA, cA, var)
            res["B"][var] = contrast(evB, cB, var)
        if thr == 10:
            res["A"]["original_truncated_simpson_tailwind"] = contrast(evA, cA, "SpringRoute_Tailwind_mps")
            res["A"]["simpson_signed"] = contrast(evA, cA, "along_simpson")
            evB = evB.assign(along_night=evB.along * evB.night)
            cB = cB.assign(along_night=cB.along * cB.night)
            res["B"]["clogit_with_night_interaction"] = clogit(evB, cB, ["along", "along_night", "cross", "ws"])
        res["B"]["clogit"] = clogit(evB, cB, ["along", "cross", "ws"])
        out[thr] = res

    _plot_event_forest(out, figdir)
    _plot_per_bird(out[10], figdir)
    return out


def _plot_event_forest(out: dict, figdir: Path) -> None:
    rows = [("原5d-2C：截断助推风（段内平均）", out[10]["A"]["original_truncated_simpson_tailwind"], MUTED)]
    for thr in [10, 30, 50]:
        rows.append((f"有符号起点沿轴风 · 月份层 · ≥{thr} km", out[thr]["A"]["along"], DAY))
    for thr in [10, 30, 50]:
        rows.append((f"有符号起点沿轴风 · 时间局地层 · ≥{thr} km", out[thr]["B"]["along"], ACCENT))
    ws_rows = [(f"10 m风速 · 月份层 · ≥{thr} km", out[thr]["A"]["ws"], DAY) for thr in [10, 30, 50]]
    ws_rows += [(f"10 m风速 · 时间局地层 · ≥{thr} km", out[thr]["B"]["ws"], ACCENT) for thr in [10, 30, 50]]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=False)
    for ax, data, title in [(axes[0], rows, "(a) 沿迁徙轴风（事件起始 − 对照，m/s）"),
                            (axes[1], ws_rows, "(b) 10 m风速（事件起始 − 对照，m/s）")]:
        y = np.arange(len(data))[::-1]
        for yi, (lab, r, col) in zip(y, data):
            ax.hlines(yi, r["ci"][0], r["ci"][1], color=col, linewidth=2)
            ax.plot(r["estimate"], yi, "o", color=col, markersize=8, markeredgecolor="white", markeredgewidth=1.5)
        ax.axvline(0, color=INK2, linewidth=1, linestyle="--")
        ax.set_yticks(y, [d[0] for d in data], fontsize=9)
        ax.set_title(title, loc="left", color=INK)
        ax.grid(axis="y", visible=False)
    axes[0].set_xlabel("鸟等权差值及整鸟Bootstrap 95%区间")
    axes[1].set_xlabel("鸟等权差值及整鸟Bootstrap 95%区间")
    fig.suptitle("图A3  事件起始风场：原“无效应”结论来自截断变量；有符号沿轴风在各阈值和设计下均为正",
                 x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "A3_event_onset_wind_reanalysis.png")
    plt.close(fig)


def _plot_per_bird(res10: dict, figdir: Path) -> None:
    pb = pd.Series(res10["A"]["along"]["per_bird"]).sort_values()
    names = [b.replace("(GWFG)", "") for b in pb.index]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = np.arange(len(pb))
    cols = [DAY if v > 0 else NIGHT for v in pb.values]
    ax.hlines(y, 0, pb.values, color=cols, linewidth=2)
    ax.plot(pb.values, y, "o", color=INK2, markersize=0)
    ax.scatter(pb.values, y, s=64, c=cols, edgecolors="white", linewidths=1.5, zorder=3)
    ax.axvline(0, color=INK2, linewidth=1)
    est, (lo, hi) = res10["A"]["along"]["estimate"], res10["A"]["along"]["ci"]
    ax.axvspan(lo, hi, color=DAY, alpha=0.08, linewidth=0)
    ax.axvline(est, color=DAY, linewidth=1, linestyle="--")
    ax.text(hi + 0.08, 1.0, f"鸟等权均值 {est:+.2f} m/s\n95%区间 [{lo:+.2f}, {hi:+.2f}]\n（阴影与虚线）",
            color=INK, fontsize=9, va="center")
    ax.set_yticks(y, names)
    ax.set_xlabel("事件起始 − 同层对照的起点沿轴风（m/s）")
    ax.set_title(f"图A4  14只鸟中{res10['A']['along']['birds_positive']}只的事件起始沿轴风高于对照",
                 loc="left", color=INK, fontsize=12)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(figdir / "A4_event_onset_per_bird.png")
    plt.close(fig)


def _plot_latitude(lat: dict, figdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.2))
    names = list(lat)
    y = np.arange(len(names))[::-1]
    for ax, term, title in [(axes[0], "night", "(a) 夜间主效应（log-odds）"),
                            (axes[1], "along_w", "(b) 白天沿轴风斜率（log-odds，每+1 m/s）")]:
        for yi, n in zip(y, names):
            est, lo, hi = lat[n][term]
            col = NIGHT if term == "night" else DAY
            ax.hlines(yi, lo, hi, color=col, linewidth=2)
            ax.plot(est, yi, "o", color=col, markersize=8, markeredgecolor="white", markeredgewidth=1.5)
        ax.axvline(0, color=INK2, linewidth=1, linestyle="--")
        ax.set_yticks(y, [f"{n}（{lat[n]['events']}/{lat[n]['n']}）" for n in names], fontsize=9)
        ax.set_title(title, loc="left", color=INK)
        ax.grid(axis="y", visible=False)
    fig.suptitle("图A5  ≥30 km发生模型：剔除高纬晚季样本后夜间效应区间跨零，沿轴风效应保持",
                 x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "A5_latitude_sensitivity.png")
    plt.close(fig)


# ---------------------------------------------------------------- 主程序
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--package", required=True, type=Path, help="阶段数据包根目录")
    ap.add_argument("--out", required=True, type=Path, help="输出目录（写入audit_results.json与figures/）")
    args = ap.parse_args()
    figdir = args.out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    sp = load_spring(args.package)
    results = {
        "composition": audit_composition(sp, figdir),
        "route_axis": audit_axis(sp),
        "continuity": audit_continuity(sp, figdir),
        "latitude_sensitivity": audit_latitude_sensitivity(sp),
        "event_onset": audit_events(sp, figdir),
    }
    _plot_latitude(results["latitude_sensitivity"], figdir)
    (args.out / "audit_results.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: results["event_onset"][10]["A"][k]["estimate"]
                      for k in ["along", "ws", "original_truncated_simpson_tailwind"]}, indent=1))


if __name__ == "__main__":
    main()
