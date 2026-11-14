import sys, json
import numpy as np, pandas as pd
P = sys.argv[1]
use = ["Bird_ID","Start_Time_UTC","Start_Longitude","Start_Latitude","End_Longitude","End_Latitude","Distance_km",
       "Migration_Stage","Day_Night","Start_Sunset_Night_Date","Local_Analysis_Date","Start_Apparent_Solar_Hour","Winter_Season"]
d = pd.read_csv(P, usecols=lambda c: c in use, low_memory=False)
w = d[d.Migration_Stage == "Winter"].copy()
def hav(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    a = np.sin((lat2-lat1)/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2
    return 6371*2*np.arcsin(np.sqrt(a))
wn = w[w.Day_Night == "Night"].sort_values(["Bird_ID","Start_Time_UTC"]).copy()
wn["night"] = pd.to_datetime(wn.Start_Sunset_Night_Date.astype(str), errors="coerce")
wd = w[w.Day_Night == "Day"].copy()
wd["date"] = pd.to_datetime(wd.Local_Analysis_Date.astype(str), errors="coerce")
dc = wd.groupby(["Bird_ID","date"]).agg(lon=("Start_Longitude","median"), lat=("Start_Latitude","median"), n=("Start_Longitude","size")).reset_index()
rows = []
for (b, n), g in wn.groupby(["Bird_ID","night"]):
    if len(g) < 3: continue
    lon = np.r_[g.Start_Longitude.iloc[0], g.End_Longitude.values]
    lat = np.r_[g.Start_Latitude.iloc[0], g.End_Latitude.values]
    d0 = hav(lon[0], lat[0], lon, lat)
    if d0.max() < 1 or d0[-1] < 1: continue
    # hour of first >=1km segment
    k = int(np.argmax(g.Distance_km.values >= 1)) if (g.Distance_km.values >= 1).any() else None
    hr = g.Start_Apparent_Solar_Hour.values[k] if k is not None else np.nan
    nxt = dc[(dc.Bird_ID == b) & (dc.date == n + pd.Timedelta(days=1))]
    prv = dc[(dc.Bird_ID == b) & (dc.date == n)]
    rec = dict(bird=b, night=str(n.date()), shift=d0[-1], hour=hr)
    if len(nxt):
        rec["start_to_nextday"] = float(hav(lon[0], lat[0], nxt.lon.values[0], nxt.lat.values[0]))
        rec["end_to_nextday"] = float(hav(lon[-1], lat[-1], nxt.lon.values[0], nxt.lat.values[0]))
    if len(prv):
        rec["start_to_prevday"] = float(hav(lon[0], lat[0], prv.lon.values[0], prv.lat.values[0]))
        rec["end_to_prevday"] = float(hav(lon[-1], lat[-1], prv.lon.values[0], prv.lat.values[0]))
    rows.append(rec)
r = pd.DataFrame(rows)
out = {"n_one_way_nights": len(r)}
hh = r.hour.dropna()
out["hour_of_first_ge1km_segment_counts"] = hh.round().value_counts().sort_index().to_dict()
x = r.dropna(subset=["start_to_nextday"])
out["n_with_nextday"] = len(x)
out["share_end_closer_to_nextday_centroid"] = float((x.end_to_nextday < x.start_to_nextday).mean())
out["median_end_to_nextday_km"] = float(x.end_to_nextday.median())
out["median_start_to_nextday_km"] = float(x.start_to_nextday.median())
y = r.dropna(subset=["start_to_prevday"])
out["share_end_farther_from_prevday_centroid"] = float((y.end_to_prevday > y.start_to_prevday).mean())
# winter extent per bird: spread of locations (90th pct distance from bird median)
ext = {}
for b, g in w.groupby("Bird_ID"):
    mlon, mlat = g.Start_Longitude.median(), g.Start_Latitude.median()
    ext[b] = round(float(np.quantile(hav(g.Start_Longitude, g.Start_Latitude, mlon, mlat), 0.9)), 1)
out["bird_winter_spread_p90_km"] = ext
out["E_birds_nights"] = wn[wn.Bird_ID.str.startswith("E")].groupby("Bird_ID").night.nunique().to_dict()
out["H_birds_nights"] = wn[wn.Bird_ID.str.startswith("H")].groupby("Bird_ID").night.nunique().to_dict()
out["H_birds_dates"] = {b: [str(g.night.min().date()), str(g.night.max().date())] for b, g in wn[wn.Bird_ID.str.startswith("H")].groupby("Bird_ID")}
out["E_birds_dates"] = {b: [str(g.night.min().date()), str(g.night.max().date())] for b, g in wn[wn.Bird_ID.str.startswith("E")].groupby("Bird_ID")}
print(json.dumps(out, ensure_ascii=False, indent=1, default=str))
