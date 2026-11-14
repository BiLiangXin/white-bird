"""D1 design check: counts of events and informative strata only.
Does NOT estimate any moon-outcome association (keeps the confirmatory test clean)."""
import sys, json, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_body
import astropy.units as u

P = sys.argv[1]
use = ["Bird_ID","Start_Time_UTC","End_Time_UTC","Start_Longitude","Start_Latitude","Distance_km",
       "Migration_Stage","Day_Night","Start_Apparent_Solar_Hour","Start_Sunset_Night_Date"]
d = pd.read_csv(P, usecols=lambda c: c in use, low_memory=False)
w = d[(d.Migration_Stage == "Winter") & (d.Day_Night == "Night") & d.Bird_ID.str.startswith("E")].copy()
t0 = pd.to_datetime(w.Start_Time_UTC, utc=True); t1 = pd.to_datetime(w.End_Time_UTC, utc=True)
mid = t0 + (t1 - t0) / 2
t = Time(mid.dt.tz_convert(None).values.astype("datetime64[ns]"))
loc = EarthLocation(lon=w.Start_Longitude.values*u.deg, lat=w.Start_Latitude.values*u.deg)
malt = get_body("moon", t, loc).transform_to(AltAz(obstime=t, location=loc)).alt.deg
w["moon_up"] = malt > 0
w["y"] = w.Distance_km >= 1
hr = w.Start_Apparent_Solar_Hour.values
w["part"] = np.where((hr >= 17) & (hr < 22), "evening", np.where((hr >= 22) | (hr < 0.5), "midnight", "late"))
w["night"] = w.Bird_ID + "|" + w.Start_Sunset_Night_Date.astype(str)
out = {"segments": len(w), "events_ge1km": int(w.y.sum()), "events_by_part": w.groupby("part").y.sum().astype(int).to_dict()}
g = w.groupby("night").agg(ny=("y", "sum"), n=("y", "size"), mu_min=("moon_up", "min"), mu_max=("moon_up", "max"))
g["y_var"] = (g.ny > 0) & (g.ny < g.n)
g["x_var"] = g.mu_min != g.mu_max
out["bird_nights"] = len(g)
out["nights_with_event"] = int((g.ny > 0).sum())
out["clogit_informative_strata(y_var & x_var)"] = int((g.y_var & g.x_var).sum())
out["events_in_informative_strata"] = int(g.loc[g.y_var & g.x_var, "ny"].sum())
out["events_per_bird"] = w.groupby("Bird_ID").y.sum().astype(int).to_dict()
# first-relocation timing (discrete-time hazard risk set): intervals at risk until first event
w = w.sort_values(["night", "Start_Time_UTC"])
w["cum"] = w.groupby("night").y.cumsum()
w["prev_cum"] = w.groupby("night").cum.shift(fill_value=0)
risk = w[w.prev_cum == 0]
out["hazard_risk_intervals"] = len(risk)
out["hazard_first_events"] = int(risk.y.sum())
print(json.dumps(out, ensure_ascii=False, indent=1))
