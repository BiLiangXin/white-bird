import sys, json
import numpy as np, pandas as pd

P = sys.argv[1]
use = ["Segment_ID","Bird_ID","Tag_ID","Start_Time_UTC","End_Time_UTC","Start_Longitude","Start_Latitude",
       "End_Longitude","End_Latitude","Distance_km","Duration_h","Migration_Stage","Day_Night","Included_Main",
       "Winter_Season","Spring_Year","Start_Sunset_Night_Date","End_Device_Ground_Speed_mps","Start_Apparent_Solar_Hour",
       "Local_Analysis_Date","Contains_Twilight_Endpoint","Main_Status"]
d = pd.read_csv(P, usecols=lambda c: c in use, low_memory=False)
out = {}
out["n_rows"] = len(d)
out["main_status"] = d["Main_Status"].value_counts().to_dict() if "Main_Status" in d else None
out["stage_daynight"] = d.groupby(["Migration_Stage","Day_Night"]).size().to_dict()
d["series"] = d["Bird_ID"].str[0]

def hav(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    a = np.sin((lat2-lat1)/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2
    return 6371*2*np.arcsin(np.sqrt(a))

w = d[d["Migration_Stage"].astype(str).str.contains("Winter", case=False)].copy()
out["winter_rows"] = len(w)
out["winter_lon_range"] = [float(w.Start_Longitude.quantile(q)) for q in (0.01, 0.5, 0.99)]
out["winter_lat_range"] = [float(w.Start_Latitude.quantile(q)) for q in (0.01, 0.5, 0.99)]
out["winter_birds_by_series"] = w.groupby("series")["Bird_ID"].nunique().to_dict()
out["winter_seasons_per_bird"] = w.groupby("Bird_ID")["Winter_Season"].nunique().value_counts().sort_index().to_dict()
out["winter_birds_per_season"] = w.groupby("Winter_Season")["Bird_ID"].nunique().to_dict()
out["winter_segments_per_series_daynight"] = w.groupby(["series","Day_Night"]).size().to_dict()
wn = w[w.Day_Night == "Night"].copy()
wn["night"] = wn["Start_Sunset_Night_Date"].astype(str)
bn = wn.groupby(["Bird_ID","night"]).size()
out["winter_bird_nights"] = int(len(bn))
out["winter_segments_per_night_quantiles"] = bn.quantile([0.1,0.25,0.5,0.75,0.9]).to_dict()
out["winter_bird_nights_per_bird"] = bn.groupby(level=0).size().describe().to_dict()
# GE1km by series at night and day
w["ge1"] = w.Distance_km >= 1
out["winter_ge1_rate_series_daynight"] = w.groupby(["series","Day_Night"])["ge1"].mean().round(4).to_dict()
# per-bird night ge1 rate
out["winter_night_ge1_rate_by_bird"] = wn.assign(ge1=wn.Distance_km>=1).groupby("Bird_ID")["ge1"].mean().round(3).to_dict()

# Night geometry: within each bird-night, ordered segments; positions = start of first + ends
wn = wn.sort_values(["Bird_ID","Start_Time_UTC"])
geo = []
for (b, n), g in wn.groupby(["Bird_ID","night"]):
    if len(g) < 3:
        continue
    lon = np.r_[g.Start_Longitude.iloc[0], g.End_Longitude.values]
    lat = np.r_[g.Start_Latitude.iloc[0], g.End_Latitude.values]
    dist0 = hav(lon[0], lat[0], lon, lat)
    maxd = dist0.max(); last = dist0[-1]
    if maxd < 1:
        cls = "stationary(<1km)"
    elif last < 0.5:
        cls = "out_and_back"
    elif last >= 1:
        cls = "one_way_shift"
    else:
        cls = "partial_return"
    geo.append((b, n, cls, maxd, last))
geo = pd.DataFrame(geo, columns=["bird","night","cls","maxd","last"])
out["winter_night_geometry_n"] = len(geo)
out["winter_night_geometry_counts"] = geo.cls.value_counts().to_dict()
out["winter_night_maxd_quantiles_moving"] = geo.loc[geo.maxd>=1, "maxd"].quantile([0.25,0.5,0.75,0.9]).round(2).to_dict()

# Day-night separation: night first fix vs previous day's daytime centroid (same bird, day date = night date)
wd = w[w.Day_Night == "Day"].copy()
wd["date"] = wd["Local_Analysis_Date"].astype(str)
dc = wd.groupby(["Bird_ID","date"]).agg(lon=("Start_Longitude","median"), lat=("Start_Latitude","median"))
first = wn.groupby(["Bird_ID","night"]).first()[["Start_Longitude","Start_Latitude"]]
first = first.reset_index().rename(columns={"night":"date"})
m = first.merge(dc.reset_index(), on=["Bird_ID","date"], how="inner")
sep = hav(m.Start_Longitude, m.Start_Latitude, m.lon, m.lat)
out["winter_dayCentroid_to_firstNightFix_km_quantiles"] = pd.Series(sep).quantile([0.1,0.25,0.5,0.75,0.9]).round(2).to_dict()
out["winter_dayCentroid_to_firstNightFix_n"] = int(len(m))
# Roost switching: first-night-fix distance between consecutive nights
f2 = wn.groupby(["Bird_ID","night"]).first()[["Start_Longitude","Start_Latitude"]].reset_index()
f2["nd"] = pd.to_datetime(f2.night, errors="coerce")
f2 = f2.sort_values(["Bird_ID","nd"])
f2["plon"] = f2.groupby("Bird_ID").Start_Longitude.shift(); f2["plat"] = f2.groupby("Bird_ID").Start_Latitude.shift()
f2["gap"] = f2.groupby("Bird_ID").nd.diff().dt.days
c = f2[f2.gap == 1]
rs = hav(c.Start_Longitude, c.Start_Latitude, c.plon, c.plat)
out["winter_consecutive_night_roost_shift_km_quantiles"] = pd.Series(rs).quantile([0.25,0.5,0.75,0.9]).round(2).to_dict()
out["winter_consecutive_night_pairs"] = int(len(c))
out["winter_roost_shift_gt1km_share"] = float((rs > 1).mean())

# device ground speed availability
spd = d["End_Device_Ground_Speed_mps"]
out["speed_available_by_series"] = d.assign(has=spd.notna()).groupby("series")["has"].mean().round(3).to_dict()
out["winter_speed_gt3mps_by_daynight"] = w.assign(f=w.End_Device_Ground_Speed_mps > 3).groupby("Day_Night")["f"].sum().to_dict()

# Spring 45-50N night segments
s = d[d["Migration_Stage"].astype(str).str.contains("Spring", case=False)]
sn = s[s.Day_Night == "Night"]
band = sn[(sn.Start_Latitude >= 45) & (sn.Start_Latitude < 50)]
band2 = sn[(sn.Start_Latitude >= 45) & (sn.Start_Latitude <= 50) & (sn.End_Latitude >= 45) & (sn.End_Latitude <= 50)]
out["spring_night_total"] = len(sn)
out["spring_night_start45_50"] = len(band)
out["spring_night_both45_50"] = len(band2)
out["spring_birds_years"] = s.groupby("Bird_ID")["Spring_Year"].nunique().to_dict()
out["spring_months"] = pd.to_datetime(s.Start_Time_UTC).dt.month.value_counts().sort_index().to_dict()
out["winter_months"] = pd.to_datetime(w.Start_Time_UTC).dt.month.value_counts().sort_index().to_dict()
out["winter_hour_counts_night"] = wn.Start_Apparent_Solar_Hour.round().value_counts().sort_index().to_dict() if "Start_Apparent_Solar_Hour" in wn else None
def fix(o):
    if isinstance(o, dict):
        return {str(k): fix(v) for k, v in o.items()}
    return o
print(json.dumps(fix(out), ensure_ascii=False, indent=1, default=str))
