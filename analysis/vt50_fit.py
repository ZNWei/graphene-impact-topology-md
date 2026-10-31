"""Re-fit all 150 P2-F analysis observations using the archived estimator.
120 terminal configurations can be re-counted; 30 round-2 observations are
available as processed records only. No missing raw records are fabricated.
"""
import math,json,collections
from common import ROOT,load,write_csv,cell,terminal_counts,ORDER
from model_logistic import logistic_fit

def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n
    center=(p+z*z/(2*n))/den;half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0,center-half),min(1,center+half)
def main():
    rows=[];raw_checks=0
    for stage,filename,key in [('p2f_parent','p2f_parent_80_results.json','valid_results'),('p2f_round1','p2f_round1_40_results.json','new_valid_results'),('p2f_round2','p2f_round2_results.json','new_valid_results')]:
        for r in load('results/source/'+filename)[key]:
            f=r['primary_fraction'];raw=ROOT/'results/raw'/stage/r['run_id']/'final.data'
            status='processed_record_only_raw_unavailable'
            if raw.exists():
                np,nt,nd=terminal_counts(raw);assert (np,nt)==(471,4116)
                assert abs(nd/np-f)<1e-12;assert (nd/np>=.5)==bool(r['primary_event'])
                raw_checks+=1;status='verified_from_terminal_data'
            rows.append({'run_id':r['run_id'],'campaign':stage,'cell':cell(r),'seed':r['seed'],'offset_x_A':r['impact_offset_angstrom'][0],'offset_y_A':r['impact_offset_angstrom'][1],'velocity_mps':r['velocity_mps'],'primary_downstream_fraction':f,'primary_event':int(r['primary_event']),'raw_verification_status':status})
    assert len(rows)==150 and len({r['run_id'] for r in rows})==150 and raw_checks==120
    cfg=load('in/config/p2f_round2_protocol.json');arch=load('results/source/p2f_round2_results.json')['threshold_estimates']
    summaries=[];fits=[];curves=[]
    for c in ORDER:
        g=collections.defaultdict(list)
        for r in rows:
            if r['cell']==c:g[int(r['velocity_mps'])].append(r)
        groups=[]
        for v,rr in sorted(g.items()):
            k=sum(r['primary_event'] for r in rr);n=len(rr);assert n==10
            lo,hi=wilson(k,n);groups.append({'velocity_mps':v,'n':n,'events':k,'event_probability':k/n,'wilson_95_ci':[lo,hi]})
            summaries.append({'cell':c,'velocity_mps':v,'velocity_km_s':v/1000,'n':n,'events':k,'event_probability':k/n,'wilson95_low':lo,'wilson95_high':hi})
        fit=logistic_fit(groups,cfg);assert fit['status']=='estimated'
        old=arch[c[:2]+'_'+c[2:]]
        assert abs(fit['v_t_50_mps']-old['v_t_50_mps'])<1e-6
        assert max(abs(x-y) for x,y in zip(fit['wald_95_ci_mps'],old['wald_95_ci_mps']))<1e-6
        center=cfg['analysis_model']['velocity_center_mps'];scale=cfg['analysis_model']['velocity_scale_mps']
        beta=fit['slope_per_1000_mps'];alpha=(center-fit['v_t_50_mps'])/scale*beta
        fits.append({'cell':c,'n_trajectories':sum(q['n'] for q in groups),'n_velocity_groups':len(groups),'status':fit['status'],'velocity_center_mps':center,'velocity_scale_mps':scale,'intercept_alpha':alpha,'slope_beta':beta,'v_t50_mps':fit['v_t_50_mps'],'v_t50_km_s':fit['v_t_50_mps']/1000,'wald95_low_km_s':fit['wald_95_ci_mps'][0]/1000,'wald95_high_km_s':fit['wald_95_ci_mps'][1]/1000})
        low=min(g);high=max(g)
        for i in range(161):
            v=low+(high-low)*i/160;eta=beta*(v-fit['v_t_50_mps'])/1000
            curves.append({'cell':c,'velocity_km_s':v/1000,'fitted_event_probability':1/(1+math.exp(-eta))})
    assert len(summaries)==15
    write_csv('results/recomputed/Fig2_pooled_150_runs.csv',rows)
    write_csv('results/recomputed/Fig2_grouped_events.csv',summaries)
    write_csv('results/recomputed/Fig2_VT50_fit.csv',fits)
    write_csv('results/recomputed/Fig2_logistic_curves.csv',curves)
    print(json.dumps({'pooled_runs':len(rows),'raw_terminal_fraction_checks':raw_checks,'processed_only_round2_records':30,'grouped_events':15,'VT50_and_CI_match_archive':True},indent=2))
if __name__=='__main__':main()
