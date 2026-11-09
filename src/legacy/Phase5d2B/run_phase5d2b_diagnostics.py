from pathlib import Path
import numpy as np,pandas as pd,scipy.stats as st,statsmodels.formula.api as smf,statsmodels.api as sm,hashlib,zipfile,json
OUT=Path('/mnt/data/Phase5d2B_Manual')
d=pd.read_csv('/mnt/data/Phase5d2A_Manual/Phase5d2A_Relative_Wind_Data.csv',low_memory=False)
s=d[(d.Migration_Stage=='Spring') & (d.Spring_Route_Axis_Valid==True)].copy()
s['bird']=s.Bird_ID.astype(str);s['bird_year']=s.bird+'_'+s.Spring_Year.astype(int).astype(str)
s['date']=pd.to_datetime(s.Start_Time_UTC,utc=True);s['month']=s.date.dt.month.astype('category')
s['night']=(s.Day_Night=='Night').astype(int);s['long30']=(s.Distance_km>=30).astype(int)
s['log_distance']=np.log(s.Distance_km.clip(lower=1e-14))
t=np.deg2rad(s.SpringRoute_Valid_Heading_Deg)
s['along_start']=s.GEE_Start_U10_mps*np.sin(t)+s.GEE_Start_V10_mps*np.cos(t)
s['cross_start']=abs(s.GEE_Start_U10_mps*np.cos(t)-s.GEE_Start_V10_mps*np.sin(t))
for x in ['along_start','cross_start']:
 s[x+'_within']=s[x]-s.groupby('bird_year')[x].transform('mean')
s['Duration_c']=s.Duration_h-2
f='long30 ~ night+along_start_within+cross_start_within+night:along_start_within+C(bird_year)+C(month)'
g=s[s.long30==1].copy()
cf=f.replace('long30 ~','log_distance ~')
fit=smf.glm(f,data=s,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':s.bird})
cfit=smf.ols(cf,data=g).fit(cov_type='cluster',cov_kwds={'groups':g.bird})
fit_dur=smf.glm(f+'+Duration_c',data=s,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':s.bird})
cfit_dur=smf.ols(cf+'+Duration_c',data=g).fit(cov_type='cluster',cov_kwds={'groups':g.bird})

rows=[]
tcrit=st.t.ppf(.975,13)
for nm, ft, unit in [('>=30km_event',fit,'logit'),('>=30km_distance',cfit,'log_km')]:
 names=['along_start_within','along_start_within','cross_start_within','night:along_start_within']
 vecs=[{'along_start_within':1},{'along_start_within':1,'night:along_start_within':1},
       {'cross_start_within':1},{'night:along_start_within':1}]
 for label,mapping in zip(['Day_alongwind','Night_alongwind','Crosswind_both','Night_vs_Day_alongwind'],vecs):
  coef=sum(ft.params[k]*v for k,v in mapping.items())
  vv=sum(ft.cov_params().loc[a,b]*u*v for a,u in mapping.items() for b,v in mapping.items())
  se=np.sqrt(vv)
  rows.append({'Model':nm,'Contrast':label,'Estimate':coef,'SE_cluster':se,'CI95_normal_low':coef-1.96*se,'CI95_normal_high':coef+1.96*se,
               'CI95_t13_reference_low':coef-tcrit*se,'CI95_t13_reference_high':coef+tcrit*se,
               'Exponentiated':np.exp(coef),'Unit':unit,'N_clusters':14,
               'Qualification':'Cluster covariance with only 14 birds; t(13) reference sensitivity, not exact small-sample inference'})
contrasts=pd.DataFrame(rows)
contrasts.to_csv(OUT/'Phase5d2B_Day_Night_Wind_Slopes.csv',index=False,encoding='utf-8-sig')

# contiguous sequence within bird+tag, identical shared event and timestamps; no artificial joins.
q=s.sort_values(['bird','Tag_ID','date']).copy(); prev=q.groupby(['bird','Tag_ID'],sort=False)
q['prev_seg']=prev.Segment_ID.shift(1);q['prev_end_id']=prev.End_Event_ID.shift(1)
q['prev_end_time']=pd.to_datetime(prev.End_Time_UTC.shift(1),utc=True)
q['prev_end_lon']=prev.End_Longitude.shift(1);q['prev_end_lat']=prev.End_Latitude.shift(1)
q['contig']=((q.Start_Event_ID==q.prev_end_id)&(q.date==q.prev_end_time)&
             np.isclose(q.Start_Longitude,q.prev_end_lon,atol=1e-8,rtol=0)&
             np.isclose(q.Start_Latitude,q.prev_end_lat,atol=1e-8,rtol=0))
q['pr_resid']=(q.long30-fit.predict(q))/np.sqrt(np.clip(fit.predict(q)*(1-fit.predict(q)),1e-10,None))
q['prev_pr_resid']=prev.pr_resid.shift(1)
q2=q[q.contig]
lag1=np.corrcoef(q2.pr_resid,q2.prev_pr_resid)[0,1]
# conditional residual for contiguous long movements
r=g.sort_values(['bird','Tag_ID','date']).copy();rr=r.groupby(['bird','Tag_ID'],sort=False)
r['prev_end_id']=rr.End_Event_ID.shift(1);r['prev_end_time']=pd.to_datetime(rr.End_Time_UTC.shift(1),utc=True)
r['prev_end_lon']=rr.End_Longitude.shift(1);r['prev_end_lat']=rr.End_Latitude.shift(1)
r['contig']=((r.Start_Event_ID==r.prev_end_id)&(r.date==r.prev_end_time)&
             np.isclose(r.Start_Longitude,r.prev_end_lon,atol=1e-8,rtol=0)&
             np.isclose(r.Start_Latitude,r.prev_end_lat,atol=1e-8,rtol=0))
r['resid']=cfit.resid;r['prev_resid']=rr.resid.shift(1)
q3=r[r.contig]
lag1cond=np.corrcoef(q3.resid,q3.prev_resid)[0,1] if len(q3)>2 else np.nan

diagnostics={'event_model_resid_lag1_pearson':float(lag1),'N_contiguous_event_pairs':int(len(q2)),
             'conditional_model_resid_lag1_pearson':float(lag1cond),'N_contiguous_conditional_pairs':int(len(q3)),
             'conditional_resid_skew':float(st.skew(cfit.resid)),
             'conditional_resid_excess_kurtosis':float(st.kurtosis(cfit.resid)),
             'duration_sensitivity':{
                model: {k:{'beta':float(ft.params[k]),'CI95_cluster':list(map(float,ft.conf_int().loc[k]))}
                         for k in ['along_start_within','cross_start_within','night:along_start_within','Duration_c']}
                for model,ft in [('event',fit_dur),('distance',cfit_dur)]}}
(OUT/'Phase5d2B_Additional_Diagnostics.json').write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2))

with open(OUT/'Phase5d2B_Manual_Model_Report.md','a',encoding='utf-8') as fp:
 fp.write('\n## 10. 独立补充：昼夜斜率合成、采样时长与残差审查\n')
 fp.write('夜间沿迁徙轴风系数通过「沿轴风主项+夜间交互项」合成，采用模型聚类协方差计算SE；另给t(13)临界值作为仅14个鸟簇时的参考敏感性，不是精确小样本置信区间。\n\n')
 fp.write('| 部分 | 对比 | β | t(13)参考95%区间 |\n|---|---|---:|---:|\n')
 for _,r0 in contrasts.iterrows():
  fp.write(f'| {r0.Model} | {r0.Contrast} | {r0.Estimate:+.4f} | [{r0.CI95_t13_reference_low:+.4f}, {r0.CI95_t13_reference_high:+.4f}] |\n')
 fp.write(f'\n二项模型Pearson残差真实连续邻段相关约{lag1:+.3f}（{len(q2):,}对）；长位移条件回归残差真实连续邻段相关约{lag1cond:+.3f}（{len(q3)}对），后者连续段对数量很少，不能作强结论。条件回归残差偏度{st.skew(cfit.resid):+.3f}、超额峰度{st.kurtosis(cfit.resid):+.3f}。聚类SE只校正不确定性近似，不表示动态残差结构被建模。\n')
 fp.write('另将实际Duration_h减2作为额外控制项拟合的敏感性结果记录在`Phase5d2B_Additional_Diagnostics.json`；不替换预先设定的主模型。\n')

print('Contrasts\n',contrasts[['Model','Contrast','Estimate','CI95_t13_reference_low','CI95_t13_reference_high']].to_string(index=False))
print('Diagnostics',json.dumps(diagnostics,ensure_ascii=False,indent=2))
