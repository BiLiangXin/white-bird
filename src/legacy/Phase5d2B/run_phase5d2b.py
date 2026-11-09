from pathlib import Path
import numpy as np, pandas as pd, scipy.stats as st
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hashlib, json, zipfile, warnings, textwrap
warnings.filterwarnings('ignore', category=FutureWarning)

src=Path('/mnt/data/Phase5d2A_Manual/Phase5d2A_Relative_Wind_Data.csv')
OUT=Path('/mnt/data/Phase5d2B_Manual');OUT.mkdir(exist_ok=True)
D=pd.read_csv(src, low_memory=False)
assert len(D)==17784 and D.Segment_ID.nunique()==17784
S=D[(D.Migration_Stage=='Spring') & (D.Spring_Route_Axis_Valid == True)].copy()
assert len(S)==9779 and S.Bird_ID.nunique()==14
S['bird']=S.Bird_ID.astype(str)
S['year']=S.Spring_Year.astype(int)
S['bird_year']=S.bird+'_'+S.year.astype(str)
S['date']=pd.to_datetime(S.Start_Time_UTC,utc=True)
S['month']=S.date.dt.month.astype('category')
S['night']=(S.Day_Night=='Night').astype(int)
S['long30']=(S.Distance_km>=30).astype(int)
S['long50']=(S.Distance_km>=50).astype(int)
S['log_distance']=np.log(S.Distance_km.clip(lower=1e-14))
S['log1p_distance']=np.log1p(S.Distance_km)
S['theta']=np.deg2rad(S.SpringRoute_Valid_Heading_Deg)
for tag,up,vp in [('start','GEE_Start_U10_mps','GEE_Start_V10_mps'),('path','GEE_U10_mps_Simpson_Mean','GEE_V10_mps_Simpson_Mean')]:
    U=S[up].to_numpy();V=S[vp].to_numpy();theta=S.theta.to_numpy()
    S[f'along_{tag}']=U*np.sin(theta)+V*np.cos(theta)
    S[f'cross_{tag}']=np.abs(U*np.cos(theta)-V*np.sin(theta))
    S[f'along_{tag}_within']=S[f'along_{tag}']-S.groupby('bird_year')[f'along_{tag}'].transform('mean')
    S[f'cross_{tag}_within']=S[f'cross_{tag}']-S.groupby('bird_year')[f'cross_{tag}'].transform('mean')
S['latitude_within_10km']= (S.Start_Latitude - S.groupby('bird_year').Start_Latitude.transform('mean'))/1 # per degree
S['is_tail_start']=S.along_start>0
S['is_head_start']=S.along_start<0
assert S[['along_start','cross_start','along_path','cross_path']].notna().all().all()
assert S.long30.sum()==540

# Use bird-year and month fixed effects to adjust individual year and seasonal timing.
base_F='long30 ~ night + C(bird_year) + C(month)'
main_F='long30 ~ night + along_start_within + cross_start_within + night:along_start_within + C(bird_year) + C(month)'
cond_F='log_distance ~ night + along_start_within + cross_start_within + night:along_start_within + C(bird_year) + C(month)'

def fit_logit(form, df, clustered=True):
    model=smf.glm(form,data=df,family=sm.families.Binomial())
    fit=model.fit(maxiter=150,disp=0)
    if clustered:
        return model.fit(maxiter=150,disp=0,cov_type='cluster',cov_kwds={'groups':df.bird},use_t=False)
    return fit

def fit_ols(form,df,clustered=True):
    model=smf.ols(form,data=df)
    if clustered:
        return model.fit(cov_type='cluster',cov_kwds={'groups':df.bird},use_t=False)
    return model.fit()

log_base=fit_logit(base_F,S)
log_main=fit_logit(main_F,S)
events=S[S.long30==1].copy()
cond_main=fit_ols(cond_F,events)
# Sensitivity: conditioned on start latitude, matched by route phase progression.
log_lat=fit_logit(main_F+' + latitude_within_10km',S)
cond_lat=fit_ols(cond_F+' + latitude_within_10km',events)
# Path-mean weather is knowingly outcome-contaminated as sensitivity only.
path_F=main_F.replace('along_start_within','along_path_within').replace('cross_start_within','cross_path_within')
path_cond_F=cond_F.replace('along_start_within','along_path_within').replace('cross_start_within','cross_path_within')
log_path=fit_logit(path_F,S)
cond_path=fit_ols(path_cond_F,events)
# 50km outcome as pre-existing threshold sensitivity.
log50=fit_logit(main_F.replace('long30','long50'),S)
# Only long events with previous-heading >=1km; compare current seasonal vs prior heading signs.
L=events[events.H1_Previous_Segment_Heading_Deg.notna()].copy()
L['H1_Headwind_or_Tailwind']=np.where(L.H1_AlongTrackWind_mps>0,'Tailwind','Headwind')
L['seasonal_tail']=L.along_start>0
L['h1_tail']=L.H1_AlongTrackWind_mps>0 # path wind! calculate start wind below
Theta_h1=np.deg2rad(L.H1_Previous_Segment_Heading_Deg)
L['h1_along_start']=L.GEE_Start_U10_mps*np.sin(Theta_h1)+L.GEE_Start_V10_mps*np.cos(Theta_h1)
L['h1_tail_start']=L.h1_along_start>0

# Output fitted coefficients for each model with standardized reporting.
fits={
 'Probability_M0_birdYear_month':(log_base,'logit'),
 'Probability_M1_startWind':(log_main,'logit'),
 'Probability_M1_startWind_plus_startLatitude':(log_lat,'logit'),
 'Probability_M1_pathMean_sensitivity':(log_path,'logit'),
 'Probability_long50_startWind_sensitivity':(log50,'logit'),
 'Distance_given_ge30_startWind':(cond_main,'ols_log'),
 'Distance_given_ge30_plus_startLatitude':(cond_lat,'ols_log'),
 'Distance_given_ge30_pathMean_sensitivity':(cond_path,'ols_log')}
coef_rows=[]
for model_name,(fit,typ) in fits.items():
    ci=fit.conf_int()
    for term in fit.params.index:
        if term not in ['Intercept','night','along_start_within','cross_start_within','night:along_start_within',
                        'along_path_within','cross_path_within','night:along_path_within','latitude_within_10km']:
            continue
        estimate=float(fit.params[term]);lo=float(ci.loc[term,0]);hi=float(ci.loc[term,1])
        coef_rows.append({'Model':model_name,'Model_Type':typ,'Term':term,'Estimate':estimate,
                          'Cluster_SE':float(fit.bse[term]),'Cluster_P_approx':float(fit.pvalues[term]),
                          'CI95_Low_approx':lo,'CI95_High_approx':hi,
                          'Exp_Estimate':np.exp(estimate), 'Exp_CI95_Low':np.exp(lo),'Exp_CI95_High':np.exp(hi),
                          'N_segments':int(fit.nobs), 'N_birds':14,
                          'Inference_Qualification':'Exploratory: 14 bird clusters, small-cluster sandwich approximation; not confirmatory'})
coef_df=pd.DataFrame(coef_rows)
coef_df.to_csv(OUT/'Phase5d2B_Model_Coefficients.csv',index=False,encoding='utf-8-sig')

# Bird-balanced descriptive contrast: tailwind vs headwind risk and conditional distance by day/night.
birdrows=[]
for (bird, dn),group in S.groupby(['bird','Day_Night']):
    for windclass,mask in [('Tailwind',group.along_start>0),('Headwind',group.along_start<0)]:
        g=group[mask]
        cg=g[g.long30==1]
        birdrows.append({'Bird_ID':bird,'Day_Night':dn,'Wind_Class':windclass,
                         'N_segments':len(g),'N_long30':int(g.long30.sum()),
                         'Long30_rate':float(g.long30.mean()) if len(g) else np.nan,
                         'N_long30_distance':len(cg),'Median_Long30_Distance_km':float(cg.Distance_km.median()) if len(cg) else np.nan,
                         'Mean_Long30_Distance_km':float(cg.Distance_km.mean()) if len(cg) else np.nan})
bird_df=pd.DataFrame(birdrows)
bird_df.to_csv(OUT/'Phase5d2B_Bird_Level_Wind_Contrasts.csv',index=False,encoding='utf-8-sig')

summ_rows=[]
for dn,group in [('All',S),('Day',S[S.night==0]),('Night',S[S.night==1])]:
    for label,sub in [('Tailwind',group[group.along_start>0]),('Headwind',group[group.along_start<0])]:
        pos=sub[sub.long30==1]
        summ_rows.append({'Day_Night':dn,'Wind_Class':label,'N_segments':len(sub),
                          'N_birds':sub.bird.nunique(),'N_ge30':int(sub.long30.sum()),
                          'Long30_rate_pct':sub.long30.mean()*100,'Median_Distance_all_km':sub.Distance_km.median(),
                          'Median_Distance_when_ge30_km':pos.Distance_km.median() if len(pos) else np.nan,
                          'Mean_Distance_when_ge30_km':pos.Distance_km.mean() if len(pos) else np.nan,
                          'N_birds_with_ge30':pos.bird.nunique(), 'Median_AlongWind_mps':sub.along_start.median()})
summ_df=pd.DataFrame(summ_rows)
summ_df.to_csv(OUT/'Phase5d2B_Descriptive_Summary.csv',index=False,encoding='utf-8-sig')

# Bird-equal differences: only include birds with both wind classes and sufficient samples.
# For threshold event risk, 1+ per wind class is sufficient to calculate per-bird rates,
# but report more conservative 20 observations each for summary.
B=bird_df.pivot(index=['Bird_ID','Day_Night'],columns='Wind_Class',values=['N_segments','N_long30','Long30_rate'])
minn=20
eligible=(B[('N_segments','Tailwind')]>=minn)&(B[('N_segments','Headwind')]>=minn)
bird_delta=(B.loc[eligible,('Long30_rate','Tailwind')]-B.loc[eligible,('Long30_rate','Headwind')]).rename('Bird_equal_rate_difference')
equal_delta=[]
for dn in ['Day','Night']:
    arr=bird_delta.xs(dn,level='Day_Night').dropna().to_numpy()
    equal_delta.append({'Day_Night':dn,'N_eligible_birds':len(arr),
                        'Mean_tail_minus_head_probability':float(np.mean(arr)),
                        'Median_tail_minus_head_probability':float(np.median(arr))})
# Observational bird-level bootstrap (not re-fitting individual segment models).
rng=np.random.default_rng(29522026)
bootrows=[]
for dn in ['Day','Night']:
    arr=bird_delta.xs(dn,level='Day_Night').dropna().to_numpy()
    bs=np.mean(rng.choice(arr,size=(10000,len(arr)),replace=True),axis=1)
    bootrows.append({'Day_Night':dn,'N_eligible_birds':len(arr),'Bird_equal_risk_difference':arr.mean(),
                     'Bird_bootstrap_CI_Low':np.quantile(bs,.025),'Bird_bootstrap_CI_High':np.quantile(bs,.975),
                     'Qualification':'Bird-level descriptive bootstrap, no confounder adjustment, only birds with 20+ records per class'})
boots_df=pd.DataFrame(bootrows)
boots_df.to_csv(OUT/'Phase5d2B_Bird_Equal_Bootstrap.csv',index=False,encoding='utf-8-sig')

# Leave one bird out for primary association, a fast fixed-effects refit without clusters.
loo_rows=[]
for bird in sorted(S.bird.unique()):
    train=S[S.bird != bird];trainC=train[train.long30==1]
    try:
        pfit=fit_logit(main_F,train,clustered=False)
        cfit=fit_ols(cond_F,trainC,clustered=False)
        for mdl,fit in [('event',pfit),('conditional_distance',cfit)]:
            row={'Model':mdl,'Left_Out_Bird':bird,'N_training':int(fit.nobs),'N_birds':train.bird.nunique()}
            for term in ['night','along_start_within','cross_start_within','night:along_start_within']:
                row[term]=float(fit.params.get(term,np.nan))
            loo_rows.append(row)
    except Exception as exc:
        loo_rows.append({'Left_Out_Bird':bird,'Fit_Error':str(exc)})
loo_df=pd.DataFrame(loo_rows)
loo_df.to_csv(OUT/'Phase5d2B_Leave_One_Bird_Out.csv',index=False,encoding='utf-8-sig')

# Leave-one-bird-out predictive comparison, no BirdYear FE because unknown bird-years on holdout.
# Baseline vs start-weather model, both same sample.
cv_base='long30 ~ night + C(month)'
cv_weather='long30 ~ night + along_start + cross_start + night:along_start + C(month)'
cv_rows=[]
for bird in sorted(S.bird.unique()):
    train=S[S.bird != bird]; test=S[S.bird == bird]
    try:
        m0=fit_logit(cv_base,train,clustered=False);m1=fit_logit(cv_weather,train,clustered=False)
        preds={'Baseline':np.clip(m0.predict(test),1e-7,1-1e-7),
               'Wind':np.clip(m1.predict(test),1e-7,1-1e-7)}
        y=test.long30.to_numpy()
        for name,p in preds.items():
            cv_rows.append({'Bird_ID':bird,'Model':name,'N':len(y),'Events':int(y.sum()),
                            'Brier':float(np.mean((p-y)**2)),
                            'LogLoss':float(-np.mean(y*np.log(p)+(1-y)*np.log(1-p))),
                            'Predicted_Events':float(np.sum(p))})
    except Exception as exc:
        cv_rows.append({'Bird_ID':bird,'Model':'Error','Error':str(exc)})
cv_df=pd.DataFrame(cv_rows)
cv_df.to_csv(OUT/'Phase5d2B_CrossBird_Prediction.csv',index=False,encoding='utf-8-sig')

# Diagnostic summary.
check_rows=[]
for nm,df in [('All_route_valid_Spring',S),('Long30_only',events)]:
    check_rows.append({'Audit':nm,'N':len(df),'Birds':df.bird.nunique(),'BirdYears':df.bird_year.nunique(),
                       'Day':int((df.night==0).sum()),'Night':int((df.night==1).sum()),
                       'Long30':int(df.long30.sum()),'Wind_start_missing':int(df.along_start.isna().sum()),
                       'Cross_start_missing':int(df.cross_start.isna().sum())})
pd.DataFrame(check_rows).to_csv(OUT/'Phase5d2B_Coverage_Audit.csv',index=False,encoding='utf-8-sig')

# Curve of predicted event risk as function of wind within-bird-year anomaly:
# G-computation on each observed stratum and month, fix cross within at zero and night at chosen value,
# average equally across birds and month distribution within each bird.
grid=np.linspace(-5,5,81)
curves=[]
for night in [0,1]:
    for along in grid:
        temp=S[['bird','bird_year','month']].copy()
        temp['night']=night;temp['along_start_within']=along;temp['cross_start_within']=0.0
        avg=log_main.predict(temp)
        bird_avgs=pd.Series(np.asarray(avg),index=temp.bird.to_numpy()).groupby(level=0).mean()
        curves.append({'Night':night,'Wind_Anomaly_mps':along,
                       'BirdEqual_Predicted_Prob':bird_avgs.mean(),
                       'Min_BirdMean':bird_avgs.min(),'Max_BirdMean':bird_avgs.max()})
curves_df=pd.DataFrame(curves)
curves_df.to_csv(OUT/'Phase5d2B_Model_Prediction_Curves.csv',index=False,encoding='utf-8-sig')

# Plots with default colors; no fixed style or colors.
fig,ax=plt.subplots(figsize=(9,5.5))
for dn in [0,1]:
    g=curves_df[curves_df.Night==dn]
    ax.plot(g.Wind_Anomaly_mps,100*g.BirdEqual_Predicted_Prob,
            label=('Night' if dn else 'Day'))
ax.set(xlabel='Along-route wind anomaly within bird-year (m/s)',
       ylabel='Estimated probability of >=30 km GPS displacement (%)',
       title='Spring: within-bird-year wind association (exploratory model)')
ax.axvline(0,ls='--');ax.legend();ax.grid(alpha=.25)
fig.tight_layout();fig.savefig(OUT/'01_Long_Displacement_Probability.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(9,5.5))
for night in ['Day','Night']:
    g=summ_df[summ_df.Day_Night==night].set_index('Wind_Class')
    labels=['Headwind','Tailwind'];vals=g.loc[labels,'Long30_rate_pct'].to_numpy()
    shift=-.18 if night=='Day' else .18
    ax.bar(np.arange(2)+shift, vals, width=.35,label=night)
    for i,lab in enumerate(labels):
        ax.text(i+shift,vals[i]+.15,f"{int(g.loc[lab,'N_ge30'])}/{int(g.loc[lab,'N_segments'])}",ha='center',fontsize=9)
ax.set_xticks(range(2),['Against route axis','Along route axis'])
ax.set(ylabel='Observed >=30 km rate (%)',title='Observed long-displacement rate and 10m wind direction')
ax.legend();ax.grid(axis='y',alpha=.25);fig.tight_layout();fig.savefig(OUT/'02_Observed_Long_Displacement_Rates.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(9,5.5))
for night in ['Day','Night']:
    g=events[events.Day_Night==night]
    ax.scatter(g.along_start,g.Distance_km,s=12,alpha=.55,label=night)
ax.axvline(0,ls='--')
ax.set(xlabel='Start-point along-route wind (m/s; positive = assisting)',
       ylabel='GPS net displacement (km)',title='Spring >=30 km events: conditional displacement')
ax.legend();ax.grid(alpha=.25);fig.tight_layout();fig.savefig(OUT/'03_Conditional_Distance.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(9,5.5))
for i,term in enumerate(['along_start_within','night:along_start_within','cross_start_within']):
    a=loo_df[loo_df.Model=='event'][term].dropna().to_numpy()
    ax.scatter([i]*len(a),a,s=20,label=term if i==0 else None)
    ax.hlines(log_main.params[term],i-.2,i+.2,linewidth=3)
ax.axhline(0,ls='--')
ax.set_xticks(range(3),['Along-wind','Night x along-wind','Crosswind'])
ax.set(ylabel='Log-odds slope; each point leaves out one bird',
       title='Leave-one-bird sensitivity of long-displacement model')
ax.grid(axis='y',alpha=.25);fig.tight_layout();fig.savefig(OUT/'04_Leave_One_Bird_Out.png',dpi=170);plt.close(fig)

# Report numbers and model parameter extraction.
def c(name,term):
    row=coef_df[(coef_df.Model==name)&(coef_df.Term==term)]
    return row.iloc[0] if len(row) else None
pmod='Probability_M1_startWind';dmod='Distance_given_ge30_startWind'
rows=[]
for term in ['night','along_start_within','cross_start_within','night:along_start_within']:
    r=c(pmod,term);q=c(dmod,term)
    rows.append((term,r,q))

print('N_spring',len(S),'N_events',len(events),'bird-years',S.bird_year.nunique())
print('Probability model:')
for term,r,q in rows:
    print(term,'logit',round(r.Estimate,4),'CI',round(r.CI95_Low_approx,4),round(r.CI95_High_approx,4),'p',round(r.Cluster_P_approx,4),
          '| cond',round(q.Estimate,4),'CI',round(q.CI95_Low_approx,4),round(q.CI95_High_approx,4),'p',round(q.Cluster_P_approx,4))
print('Descriptions:\n',summ_df.to_string(index=False))
print('Bird equal: \n',boots_df.to_string(index=False))
print('LOO ranges:\n',loo_df.groupby('Model')[['along_start_within','night:along_start_within','cross_start_within']].agg(['min','max']).to_string())
print('CV averaged:\n',cv_df.groupby('Model')[['Brier','LogLoss']].mean().to_string())
print('H1 subset',len(L),'H1/season wind sign agreement', float((L.h1_tail_start == L.seasonal_tail).mean()),'birds',L.bird.nunique())
print('AIC event',log_base.aic,log_main.aic)
print('OLS R2',cond_main.rsquared)

# SAVE compact metrics JSON.
metrics={
 'records':len(S),'birds':int(S.bird.nunique()),'bird_years':int(S.bird_year.nunique()),
 'long30_events':int(S.long30.sum()),'long50_events':int(S.long50.sum()),
 'short_wind_start':{'along_min':float(S.along_start.min()),'along_max':float(S.along_start.max())},
 'primary_main_logit':{term:{'estimate':float(r.Estimate),'lo':float(r.CI95_Low_approx),'hi':float(r.CI95_High_approx)} for term,r,q in rows},
 'primary_distance':{term:{'estimate':float(q.Estimate),'lo':float(q.CI95_Low_approx),'hi':float(q.CI95_High_approx)} for term,r,q in rows},
 'H1_30km_events_with_local_heading':len(L),'H1_versus_route_tailwind_sign_agree':float((L.h1_tail_start==L.seasonal_tail).mean()),
 'CV_mean_by_bird':cv_df.groupby('Model')[['Brier','LogLoss']].mean().round(8).to_dict('index'),
 'AIC_event_base':float(log_base.aic),'AIC_event_wind':float(log_main.aic)}
(OUT/'Phase5d2B_Key_Metrics.json').write_text(json.dumps(metrics,indent=2,ensure_ascii=False),encoding='utf-8')

# Concise Chinese report.
report=[]
report += ['# Phase 5d-2B：春迁迁徙轴风场与GPS净位移的手动探索性模型报告','',
 '**最终判断：探索性关联分析完成；不能视为确认性的风向驱动结论。**',
 '',f'使用已冻结Phase5d2A数据，春迁路线轴有效记录 {len(S):,} 段，{S.bird.nunique()} 只鸟、{S.bird_year.nunique()} 个鸟×春迁年份；其中≥30 km净位移 {len(events)} 段（白天 {int(((S.long30==1)&(S.night==0)).sum())}，夜间 {int(((S.long30==1)&(S.night==1)).sum())}）。',
 '本轮仅使用春迁记录。未重建原GPS、未重新下载气象数据、未删除方向合格样本内的零位移和低位移记录。','',
 '## 1. 分析对象和方法',
 '- **主天气暴露：段起点** ERA5 10 m U/V风分量投影到「排除该段端点后」由鸟×春迁年份其他GPS位置构建的**整个春迁迁徙轴**。它不是实时独立路线，也不等于实际飞行层顺逆风。',
 '- **风分量：**沿季节轴风`Along = U·sinθ + V·cosθ`，>0是支持迁徙轴风；`Crosswind = |U·cosθ−V·sinθ|`。',
 '- **风场组内量：**分别减去鸟×春迁年份的均值，以估计同一只鸟同一春迁年份内的条件关联。',
 '- **第一部分：**在全部方向有效记录中，使用二项Logistic回归检验`Distance≥30 km`的发生；固定效应包括昼夜、沿轴风组内量、侧风组内量、昼夜×沿轴风、鸟×春迁年份截距、UTC月份。',
 '- **第二部分：**在发生`Distance≥30 km`的记录中，使用log(distance) OLS回归检验长位移**发生后的距离大小**，采用相同协变量。',
 '- **不确定性：**以鸟为cluster（14簇）的稳健sandwich标准误，另开展整鸟留一重拟合与鸟等权风向分组自助抽样。14簇的渐近95%区间仅供探索性解读，并非确认性小样本检验。',
 '- **敏感性：**起点纬度组内偏差、路程中点均值风、≥50 km阈值；不按显著性挑选最优模型。',
 '- **不做推断：**这不是风向与真实飞行速度的因果关系，也不能仅通过相邻定位推断鸟当时实际在飞行。','',
 '## 2. 风向分组的直接描述',
 '| 时段 | 风向类别 | GPS段 | ≥30 km段 | 发生率 | 已发生≥30 km段中位距离 |',
 '|---|---|---:|---:|---:|---:|']
for r in summ_rows:
    if r['Day_Night']=='All':continue
    report.append(f"| {r['Day_Night']} | {r['Wind_Class']} | {r['N_segments']:,} | {r['N_ge30']} | {r['Long30_rate_pct']:.2f}% | {r['Median_Distance_when_ge30_km']:.2f} km |")
report += ['', '这些描述仅按起点风沿迁徙轴投影正负划分，未控制鸟和时序，是探索性事实，不能解释为风导致长位移。','',
 '## 3. 第一部分：长距离净位移发生概率',
 '| 模型项 | 估计β(log-odds) | 14鸟聚类近似95%区间 | 条件OR |',
 '|---|---:|---:|---:|']
for term,r,q in rows:
    report.append(f'| `{term}` | {r.Estimate:+.4f} | [{r.CI95_Low_approx:+.4f}, {r.CI95_High_approx:+.4f}] | {np.exp(r.Estimate):.3f} |')
report += ['',
 '注意：沿轴风的主项对应**白天**每+1 m/s；夜间斜率等于主项+交互项，不能只解释交互项而忽略主效应。条件OR不是总体平均的概率比。',
 '', '## 4. 第二部分：达到≥30 km之后的距离大小',
 '| 模型项 | 估计β(log km) | 14鸟聚类近似95%区间 |',
 '|---|---:|---:|']
for term,r,q in rows:
    report.append(f'| `{term}` | {q.Estimate:+.4f} | [{q.CI95_Low_approx:+.4f}, {q.CI95_High_approx:+.4f}] |')
report += ['',
 '该模型只描述**已经发生长位移时**的条件距离；选择`Distance≥30 km`本身受行为和气象影响，选择条件可能造成偏差。两部分估计不能简单合并为普遍的「风使鸟飞更远」。',
 '', '## 5. 个体敏感性与风向口径',
 f'- 在 {len(L)} 个具有有效前一段局地方向的≥30 km春迁记录中，局地方向与季节轴方向判断助推/阻碍的一致率仅为 {100*(L.h1_tail_start==L.seasonal_tail).mean():.1f}%；其余长位移段没有可靠H1方向。不同风向参照不能混称为鸟真正迎风/尾风。',
 f'- 路线轴有效样本包括14只鸟、25个鸟×春迁年份；各鸟轨迹长度和长位移事件数不同，因此留一鸟、鸟等权与未加权记录汇总必须分别报告。',
 '- 模型置信区间以14鸟聚类sandwich估计，存在小样本偏差；连续时间依赖和厚尾残差不会因为聚类SE而消失。',
 '- 起点风是主暴露；沿途平均风依赖移动路径，用于敏感性，不是更严格的因果预测变量。',
 '', '## 6. 鸟等权风向对比及鸟重抽样（描述性）',
 '要求每只鸟在对应昼夜中助推和阻碍风各有≥20段，先计算该鸟长位移率之差，再等权平均鸟；10000次有放回抽鸟仅表示鸟级描述不确定性。',
 '| 昼夜 | 合格鸟数 | 助推-阻碍发生概率差 | 鸟级95%区间 |',
 '|---|---:|---:|---:|']
for r in bootrows:
    report.append(f"| {r['Day_Night']} | {r['N_eligible_birds']} | {100*r['Bird_equal_risk_difference']:+.2f}个百分点 | [{100*r['Bird_bootstrap_CI_Low']:+.2f}, {100*r['Bird_bootstrap_CI_High']:+.2f}]个百分点 |")
report += ['',
 '## 7. 留一鸟预测与稳健性',
 f"- 简化的鸟外推模型（不包含无法预测的新鸟固定截距），鸟等权Brier：Baseline {cv_df.groupby('Model').Brier.mean().get('Baseline',np.nan):.6f}；Wind {cv_df.groupby('Model').Brier.mean().get('Wind',np.nan):.6f}。鸟等权LogLoss：Baseline {cv_df.groupby('Model').LogLoss.mean().get('Baseline',np.nan):.6f}；Wind {cv_df.groupby('Model').LogLoss.mean().get('Wind',np.nan):.6f}。",
 '- 留一鸟的模型系数见CSV；如果符号或量级明显不稳定，则不得推广为普遍风向效应。',
 '- 病态或较大残差没有通过专门的分布与时间序列门槛诊断，不宣称Phase5d-1中的高斯残差问题已解决。',
 '', '## 8. 结果边界与决策',
 '**本轮可确认的是数据关联估计可以完成，而不是生物学机制已验证。** 季节性迁徙轴虽然不使用当前段端点构建，但用了其他时间的轨迹，仍包含未来行为信息；气象风为10m近地风；地理位置和季节进程残留混杂。',
 '若两部分模型方向一致、整鸟留一稳定且鸟外推有所改善，可作为论文探索性补充；否则保留昼夜净位移规律为主线，不应继续叠加新变量寻找显著性。',
 '', '## 9. 交付',
 '- 模型系数、原始描述、鸟级风向对比、鸟等权Bootstrap、留一鸟影响、鸟外推、预测曲线及覆盖审计CSV。',
 '- 4张PNG。',
 '- 全部结果为当前会话手动执行；脚本保存在同目录`run_phase5d2b.py`，并记录输入SHA-256。'
 ]
(OUT/'Phase5d2B_Manual_Model_Report.md').write_text('\n'.join(report)+'\n',encoding='utf-8')

hashes={}
for file in [src,Path('/mnt/data/Phase5a_Frozen_Movement_Data.csv'),Path('/mnt/data/Phase5c_GEE_Full_Matched_Movement_Data.csv')]:
    h=hashlib.sha256();
    with open(file,'rb') as f:
        for part in iter(lambda:f.read(2**20),b''):h.update(part)
    hashes[file.name]=h.hexdigest()
(OUT/'Source_SHA256.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')

print('OUTPUT_DIR',OUT,'COUNT',len(list(OUT.iterdir())))
