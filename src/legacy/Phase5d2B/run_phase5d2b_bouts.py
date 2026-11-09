import pandas as pd, numpy as np, json
from pathlib import Path
OUT=Path('/mnt/data/Phase5d2B_Manual')
d=pd.read_csv('/mnt/data/Phase5d2A_Manual/Phase5d2A_Relative_Wind_Data.csv',low_memory=False)
s=d[(d.Migration_Stage=='Spring') & (d.Spring_Route_Axis_Valid==True) & (d.Distance_km>=30)].copy()
s['StartTime']=pd.to_datetime(s.Start_Time_UTC,utc=True)
s['EndTime']=pd.to_datetime(s.End_Time_UTC,utc=True)
s=s.sort_values(['Bird_ID','Tag_ID','StartTime','EndTime']).reset_index(drop=True)
g=s.groupby(['Bird_ID','Tag_ID'],sort=False)
prev_id=g.End_Event_ID.shift(1)
prev_time=pd.to_datetime(g.End_Time_UTC.shift(1),utc=True)
continuity=(s.Start_Event_ID==prev_id)&(s.StartTime==prev_time)&np.isclose(s.Start_Longitude,g.End_Longitude.shift(1),atol=1e-8,rtol=0)&np.isclose(s.Start_Latitude,g.End_Latitude.shift(1),atol=1e-8,rtol=0)&(s.Day_Night==g.Day_Night.shift(1))&(s.Spring_Year==g.Spring_Year.shift(1))
s['Is_Bout_Start']=~continuity
s['Bout_Number_Bird']=s.groupby('Bird_ID')['Is_Bout_Start'].cumsum()
s['Bout_ID']=s.Bird_ID.astype(str)+'__'+s.Spring_Year.astype(int).astype(str)+'__'+s.Bout_Number_Bird.astype(str)
assert s.Bout_ID.nunique()==int((~continuity).sum())
B=s.groupby('Bout_ID',sort=False).agg(Bird_ID=('Bird_ID','first'),Spring_Year=('Spring_Year','first'),Day_Night=('Day_Night','first'),N_segments=('Segment_ID','size'),Total_GPS_Sequential_Displacement_km=('Distance_km','sum'),First_UTC=('Start_Time_UTC','first'),Last_UTC=('End_Time_UTC','last'),Duration_sum_h=('Duration_h','sum'),Start_Segment_ID=('Segment_ID','first'),End_Segment_ID=('Segment_ID','last'),Start_Lon=('Start_Longitude','first'),Start_Lat=('Start_Latitude','first'),End_Lon=('End_Longitude','last'),End_Lat=('End_Latitude','last')).reset_index()
# This is sum of successive chords not geodesic net displacement from bout start to end. Do not call true flown distance.
B.to_csv(OUT/'Phase5d2B_Long_Displacement_Sequences.csv',index=False,encoding='utf-8-sig')
summary=B.groupby('Day_Night').agg(N_sequences=('Bout_ID','size'),N_birds=('Bird_ID','nunique'),N_segments=('N_segments','sum'),Median_segments_per_sequence=('N_segments','median'),Median_sequence_duration_h=('Duration_sum_h','median'),Median_sum_chords_km=('Total_GPS_Sequential_Displacement_km','median'),Max_segments_per_sequence=('N_segments','max')).reset_index()
summary.to_csv(OUT/'Phase5d2B_Sequence_Summary.csv',index=False,encoding='utf-8-sig')
report=OUT/'Phase5d2B_Manual_Model_Report.md'
with report.open('a',encoding='utf-8') as f:
 f.write('\n## 11. 核心验收限制：时间依赖与长位移连续序列\n')
 f.write('**本轮仍未通过确认性模型要求：** 以实际共享GPS端点构建的7,653对连续春迁段，其长位移二项模型Pearson残差一阶相关约**+0.582**。这表明长位移事件常连续出现，不能把每个两小时段当作独立迁飞决策。14鸟聚类标准误不能消除这一时间依赖结构的解释与预测问题。\n\n')
 f.write(f'仅在≥30 km且参考迁徙轴有效的540段中，相邻长位移段有{int(continuity.sum())}对，据此按共享事件端点定义**{len(B)}条连续长位移序列**（而非{len(s)}次独立起飞）；其中单段序列{int((B.N_segments==1).sum())}条。这些序列可能被未纳入的短位移/晨昏段拆分，也不能由两小时定位直接确认为连续飞行。\n\n')
 f.write('| 时段 | 连续长位移序列 | 对应原GPS段 | 中位序列时长(h) | 中位连续直线段长度之和(km) |\n|---|---:|---:|---:|---:|\n')
 for _,r in summary.iterrows():
  f.write(f'| {r.Day_Night} | {int(r.N_sequences)} | {int(r.N_segments)} | {r.Median_sequence_duration_h:.1f} | {r.Median_sum_chords_km:.1f} |\n')
 f.write('\n**建议下一阶段转向「长位移连续序列/事件起始」单元，重点是事前风场与事件开始的关系；不要在逐段模型未通过时间相关诊断时继续叠加更多交互项。**\n')
print(summary.to_string(index=False))
print('Bouts',len(B),'contig pairs',int(continuity.sum()))
