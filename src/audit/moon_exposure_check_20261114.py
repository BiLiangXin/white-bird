import sys, json, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_body, get_sun
import astropy.units as u
P = sys.argv[1]
use = ["Bird_ID","Start_Time_UTC","End_Time_UTC","Mid_Longitude","Mid_Latitude","Start_Longitude","Start_Latitude",
       "Migration_Stage","Day_Night","Start_Apparent_Solar_Hour","Start_Sunset_Night_Date","Winter_Season"]
d = pd.read_csv(P, usecols=lambda c: c in use, low_memory=False)
w = d[(d.Migration_Stage == "Winter") & (d.Day_Night == "Night") & d.Bird_ID.str.startswith("E")].copy()
t0 = pd.to_datetime(w.Start_Time_UTC, utc=True); t1 = pd.to_datetime(w.End_Time_UTC, utc=True)
mid = t0 + (t1 - t0) / 2
t = Time(mid.dt.tz_convert(None).values.astype("datetime64[ns]"))
lon = w.Start_Longitude.values; lat = w.Start_Latitude.values
loc = EarthLocation(lon=lon*u.deg, lat=lat*u.deg)
moon = get_body("moon", t, loc); sun = get_sun(t)
aa = AltAz(obstime=t, location=loc)
malt = moon.transform_to(aa).alt.deg
elong = sun.separation(moon).deg
illum = (1 - np.cos(np.radians(elong))) / 2
# waxing if moon is east of sun (ecliptic longitude difference in (0,180))
from astropy.coordinates import GeocentricTrueEcliptic
lm = moon.transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg
ls = sun.transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg
waxing = ((lm - ls) % 360) < 180
w["moon_up"] = malt > 0
w["illum"] = illum
w["bright_up"] = (malt > 0) & (illum > 0.5)
hr = w.Start_Apparent_Solar_Hour.values
w["part"] = np.where((hr >= 17) & (hr < 22), "evening(18-22h)", np.where((hr >= 22) | (hr < 0.5), "midnight(22-24h)", "late(0-4h)"))
w["phase"] = np.where(illum < 0.25, "dark(<0.25)", np.where(waxing, "waxing>=0.25", "waning>=0.25"))
out = {}
out["E_winter_night_segments"] = len(w)
out["share_moon_up"] = round(float(w.moon_up.mean()), 3)
out["segments_part_x_moon_up"] = {f"{a}|{b}": int(n) for (a, b), n in w.groupby(["part","moon_up"]).size().items()}
out["segments_part_x_bright_up"] = {f"{a}|{b}": int(n) for (a, b), n in w.groupby(["part","bright_up"]).size().items()}
out["segments_phase_x_part"] = {f"{a}|{b}": int(n) for (a, b), n in w.groupby(["phase","part"]).size().items()}
# lunar cycles covered (synodic months) per winter
mid_d = mid.dt.tz_convert(None)
out["date_span_by_season"] = {s: [str(g.min().date()), str(g.max().date()), int(g.dt.date.nunique())] for s, g in mid_d.groupby(w.Winter_Season.values)}
# nights where moon state differs within the night (within-night contrast available)
w["night"] = w.Start_Sunset_Night_Date.astype(str)
g = w.groupby(["Bird_ID","night"]).moon_up.agg(["min","max","size"])
out["bird_nights"] = int(len(g))
out["bird_nights_with_within_night_moon_contrast"] = int(((g["min"] == False) & (g["max"] == True)).sum())
out["unique_calendar_nights"] = int(w.night.nunique())
print(json.dumps(out, ensure_ascii=False, indent=1))
