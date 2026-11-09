"""越冬与春季停歇期夜间活动的预检：月光效应、局地昼夜差、个体一致性。

为“夜间活动机制”研究方案提供可行性依据；属于探索性筛查，不是确认性分析。

用法：
  python -I src/audit/precheck_nocturnal_activity.py --package "<阶段数据包路径>" --out docs
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
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt  # noqa: E402

from gwfg_independent_audit import DAY, INK, INK2, NIGHT, bird_boot  # noqa: E402

warnings.filterwarnings("ignore")

COLS = ["Bird_ID", "Migration_Stage", "Day_Night", "Distance_km", "Start_Time_UTC", "Start_Longitude",
        "Start_Latitude", "Winter_Season", "Spring_Year", "GEE_Start_Temperature2m_C"]


def load(pkg: Path) -> pd.DataFrame:
    d = pd.read_csv(pkg / "02_Data/01_Frozen_Inputs/Phase5c_GEE_Full_Matched_Movement_Data.csv",
                    usecols=COLS, low_memory=False)
    assert len(d) == 17784
    d["t0"] = pd.to_datetime(d.Start_Time_UTC, utc=True)
    d["season"] = np.where(d.Migration_Stage == "Winter", "W" + d.Winter_Season.astype(str),
                           "S" + d.Spring_Year.astype(str))
    d["bird_season"] = d.Bird_ID + "|" + d.season
    d["month"] = d.t0.dt.month
    d["move1"] = (d.Distance_km >= 1).astype(int)
    # 两个驻留期：鄱阳湖越冬、春季东北停歇区（45–50°N）
    d["site"] = np.select([d.Migration_Stage == "Winter",
                           (d.Migration_Stage == "Spring") & d.Start_Latitude.between(45, 50)],
                          ["越冬期（鄱阳湖）", "春季停歇期（45–50°N）"], default="春季其他")
    return d


def add_moonlight(nt: pd.DataFrame) -> pd.DataFrame:
    t = Time(nt.t0.dt.tz_convert(None).to_numpy())
    loc = EarthLocation(lon=nt.Start_Longitude.to_numpy() * u.deg, lat=nt.Start_Latitude.to_numpy() * u.deg)
    moon, sun = get_body("moon", t), get_sun(t)
    alt = moon.transform_to(AltAz(obstime=t, location=loc)).alt.deg
    illum = (1 - np.cos(sun.separation(moon).rad)) / 2
    return nt.assign(moonlight=np.where(alt > 0, illum, 0.0))


def moon_glm(df: pd.DataFrame, form: str = "move1 ~ moonlight + C(bird_season) + C(month)"):
    m = smf.glm(form, data=df, family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(df.Bird_ID)[0]})
    ci = m.conf_int().loc["moonlight"]
    return float(m.params["moonlight"]), [float(ci[0]), float(ci[1])]


def moon_effects(d: pd.DataFrame) -> dict:
    nt = add_moonlight(d[d.Day_Night == "Night"].copy())
    nt["band"] = pd.cut(nt.moonlight, [-0.01, 0.1, 0.6, 1.0], labels=["暗夜", "中等", "亮夜"])
    out = {}
    for site in ["越冬期（鄱阳湖）", "春季停歇期（45–50°N）"]:
        s = nt[nt.site == site]
        est, ci = moon_glm(s)
        rates = s.groupby("band", observed=True).move1.agg(["mean", "size"])
        g = s.assign(bright=s.moonlight > 0.5).groupby(["Bird_ID", "bright"]).move1.agg(["mean", "size"]).unstack()
        g = g[(g[("size", True)] >= 15) & (g[("size", False)] >= 15)]
        diff = g[("mean", True)] - g[("mean", False)]
        be, lo, hi = bird_boot(diff)
        out[site] = {"n": int(len(s)), "birds": int(s.Bird_ID.nunique()), "beta": est, "ci": ci,
                     "rates": {str(k): {"rate": float(r["mean"]), "n": int(r["size"])} for k, r in rates.iterrows()},
                     "bird_equal_bright_minus_dark": {"estimate": be, "ci": [lo, hi], "birds": int(len(diff)),
                                                      "birds_positive": int((diff > 0).sum())}}

    w = nt[nt.site == "越冬期（鄱阳湖）"]
    loo = []
    for b in w.Bird_ID.unique():
        s = w[w.Bird_ID != b]
        loo.append(float(smf.glm("move1 ~ moonlight + C(bird_season) + C(month)", data=s,
                                 family=sm.families.Binomial()).fit().params["moonlight"]))
    out["winter_leave_one_bird_out_range"] = [min(loo), max(loo)]
    out["winter_by_season"] = {}
    for ws, s in w.groupby("season"):
        if s.Bird_ID.nunique() >= 3:
            est, ci = moon_glm(s)
            out["winter_by_season"][ws] = {"n": int(len(s)), "birds": int(s.Bird_ID.nunique()), "beta": est, "ci": ci}

    both = nt[nt.site != "春季其他"].assign(spring=lambda x: (x.Migration_Stage == "Spring").astype(int))
    # 月份与阶段完全共线，交互模型只保留鸟×季节固定效应
    m = smf.glm("move1 ~ moonlight + moonlight:spring + C(bird_season)", data=both,
                family=sm.families.Binomial()).fit(cov_type="cluster",
                                                   cov_kwds={"groups": pd.factorize(both.Bird_ID)[0]})
    ci = m.conf_int().loc["moonlight:spring"]
    out["interaction_staging_vs_winter"] = {"beta": float(m.params["moonlight:spring"]),
                                            "ci": [float(ci[0]), float(ci[1])]}
    return out


def local_daynight(d: pd.DataFrame) -> dict:
    """排除≥10 km迁飞段后，两个驻留期的局地昼夜差及同鸟配对差中差。"""
    gaps = {}
    out = {}
    for site in ["越冬期（鄱阳湖）", "春季停歇期（45–50°N）"]:
        x = d[(d.site == site) & (d.Distance_km < 10)]
        r = x.groupby(["Bird_ID", "Day_Night"]).move1.agg(["mean", "size"]).unstack()
        r = r[(r[("size", "Day")] >= 30) & (r[("size", "Night")] >= 30)]
        gaps[site] = r[("mean", "Day")] - r[("mean", "Night")]
        out[site] = {"birds": int(len(r)), "day": float(r[("mean", "Day")].mean()),
                     "night": float(r[("mean", "Night")].mean()), "birds_day_gt_night": int((gaps[site] > 0).sum())}
    p = pd.concat(gaps, axis=1).dropna()
    dd = p.iloc[:, 1] - p.iloc[:, 0]
    est, lo, hi = bird_boot(dd)
    out["paired_diff_in_diff"] = {"birds": int(len(p)), "estimate": est, "ci": [lo, hi],
                                  "birds_positive": int((dd > 0).sum())}
    return out


def individual_consistency(d: pd.DataFrame) -> dict:
    r = d.groupby(["Bird_ID", "Migration_Stage", "season", "Day_Night"]).move1.mean().unstack("Day_Night").dropna()
    n = d.groupby(["Bird_ID", "season"]).size().rename("n")
    r = r.join(n, on=["Bird_ID", "season"])
    r = r[r.n >= 100]
    r["log_ratio"] = np.log((r.Night + 0.01) / (r.Day + 0.01))
    pb = r.groupby(["Bird_ID", "Migration_Stage"]).log_ratio.mean().unstack().dropna()
    rho, p = spearmanr(pb.Winter, pb.Spring)
    out = {"birds_both_stages": int(len(pb)), "spearman_winter_vs_spring": float(rho), "p": float(p)}
    for st in ["Winter", "Spring"]:
        x = r.xs(st, level="Migration_Stage").reset_index()
        k = x.groupby("Bird_ID").size()
        x = x[x.Bird_ID.isin(k[k >= 2].index)]
        m = smf.mixedlm("log_ratio ~ 1", x, groups=x.Bird_ID).fit()
        vb, ve = float(m.cov_re.iloc[0, 0]), float(m.scale)
        out[f"repeatability_{st}"] = {"icc": vb / (vb + ve), "birds": int(x.Bird_ID.nunique()), "seasons": int(len(x))}
    return out


def plot(res: dict, figdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    sites = ["越冬期（鄱阳湖）", "春季停歇期（45–50°N）"]
    bands = ["暗夜", "中等", "亮夜"]
    x = np.arange(len(bands))
    for i, (site, col) in enumerate(zip(sites, [NIGHT, DAY])):
        rates = res["moon"][site]["rates"]
        vals = [rates[b]["rate"] * 100 for b in bands]
        xs = x + (i - 0.5) * 0.36
        ax.bar(xs, vals, width=0.34, color=col, label=site, edgecolor="white", linewidth=1.5)
        for xi, v, b in zip(xs, vals, bands):
            ax.text(xi, v + 0.4, f"{v:.1f}%", ha="center", color=INK2, fontsize=8.5)
    ax.set_xticks(x, ["暗夜（<0.1）", "中等", "亮夜（>0.6）"])
    ax.set_ylabel("夜间≥1 km移动段比例（%）")
    ax.set_title("(a) 月光与夜间移动：越冬期明显，停歇期不明显", loc="left", color=INK)
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="x", visible=False)

    ax = axes[1]
    rows = [(f"越冬 {k[1:]}（{v['birds']}只）", v) for k, v in res["moon"]["winter_by_season"].items()]
    rows.append((f"越冬合计（{res['moon'][sites[0]]['birds']}只）", res["moon"][sites[0]]))
    rows.append((f"春季停歇期（{res['moon'][sites[1]]['birds']}只）", res["moon"][sites[1]]))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, r) in zip(y, rows):
        col = DAY if lab.startswith("春季") else NIGHT
        ax.hlines(yi, r["ci"][0], r["ci"][1], color=col, linewidth=2)
        ax.plot(r["beta"], yi, "o", color=col, markersize=8, markeredgecolor="white", markeredgewidth=1.5)
    ax.axvline(0, color=INK2, linewidth=1, linestyle="--")
    ax.set_yticks(y, [r[0] for r in rows], fontsize=9)
    ax.set_xlabel("月光效应 β（log-odds，满月当空 vs 无月光）")
    ax.set_title("(b) 越冬期效应在各年份均为正", loc="left", color=INK)
    ax.grid(axis="y", visible=False)
    fig.suptitle("图C1  月光对夜间移动的影响随阶段不同（预检，探索性）", x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "C1_moonlight_by_stage.png")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--package", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    d = load(args.package)
    res = {"moon": moon_effects(d), "local_daynight": local_daynight(d),
           "individual_consistency": individual_consistency(d)}
    (args.out / "nocturnal_precheck.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    plot(res, args.out / "figures")
    m = res["moon"]
    print(json.dumps({"winter_beta": m["越冬期（鄱阳湖）"]["beta"], "winter_ci": m["越冬期（鄱阳湖）"]["ci"],
                      "staging_beta": m["春季停歇期（45–50°N）"]["beta"], "staging_ci": m["春季停歇期（45–50°N）"]["ci"],
                      "interaction": m["interaction_staging_vs_winter"], "loo": m["winter_leave_one_bird_out_range"],
                      "bird_equal_winter": m["越冬期（鄱阳湖）"]["bird_equal_bright_minus_dark"],
                      "local": res["local_daynight"], "consistency": res["individual_consistency"]},
                     ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
