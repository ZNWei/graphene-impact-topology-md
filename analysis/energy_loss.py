"""Recompute Fig.1, Fig.3, Table 3 and P3 from supplied raw records.
Run from any directory: python analysis/energy_loss.py
Output: results/recomputed/*.csv. Archived observations are never overwritten.
"""
import statistics,collections,json
from common import ROOT,load,write_csv,cell,track,dump_frames,undercoordination,ORDER

def main():
    specs=load('in/run_specs.json');p1=load('results/source/p1_p4_results_all.json')['results']
    main_rows=[];sensitivity=[];checks=[]
    for r in p1:
        if r['stage'] not in ('mechanism','convergence'):continue
        spec=specs[r['run_id']];raw=ROOT/'results/raw/p1_p4'/r['run_id']
        initial=(raw/'initial_state.txt').read_text(encoding='utf-8').split()
        assert initial[:2]==['INITIAL',r['run_id']]
        end=track(raw/'track.dat',spec['track_columns'])
        loss=float(initial[6])-end['ke_projectile_cm_eV']
        drift=(end['etotal_eV']-float(initial[2]))/max(1,abs(float(initial[2])))
        assert abs(loss-r['projectile_com_ke_loss_eV'])<1e-7
        assert abs(end['ke_target_eV']-r['final']['ke_target_eV'])<1e-7
        checks.append({'run_id':r['run_id'],'energy_loss_matches_archive':True,'target_ke_matches_archive':True})
        row={'run_id':r['run_id'],'cell':cell(r),'seed':r['seed'],'offset_x_A':r['impact_offset_angstrom'][0],'offset_y_A':r['impact_offset_angstrom'][1],'velocity_mps':r['velocity_mps'],'projectile_com_ke_loss_eV':loss,'terminal_target_ke_eV':end['ke_target_eV'],'energy_drift_fraction':drift}
        if r['stage']=='mechanism':
            frames=list(dump_frames(raw/'first_last_frames.lammpstrj'));assert len(frames)==2
            net=undercoordination(frames[0][1],frames[1][1],spec)
            assert net==r['structure']['target_damage_edge_excluded_net']
            row['target_undercoordination_net_atoms']=net;main_rows.append(row)
        else:
            row['axis']=r['variation']['name'];row['setting']=r['variation']['value'];sensitivity.append(row)
    assert len(main_rows)==20 and len(sensitivity)==60
    write_csv('results/recomputed/Fig1_individual_runs.csv',main_rows)
    summaries=[]
    for c in ORDER:
        rr=[r for r in main_rows if r['cell']==c];assert len(rr)==5
        row={'cell':c,'n':len(rr)}
        for k in ['projectile_com_ke_loss_eV','terminal_target_ke_eV','target_undercoordination_net_atoms']:
            row[k+'_mean']=statistics.mean(r[k] for r in rr);row[k+'_sd']=statistics.stdev(r[k] for r in rr)
        summaries.append(row)
    write_csv('results/recomputed/Fig1_summary.csv',summaries)
    g=collections.defaultdict(list)
    for r in sensitivity:g[(r['axis'],str(r['setting']))].append(r['projectile_com_ke_loss_eV'])
    ref=statistics.mean(g[('sheet_size_angstrom','60.0')]);summary=[]
    for (axis,setting),vv in sorted(g.items()):
        m=statistics.mean(vv);sd=statistics.stdev(vv)
        summary.append({'axis':axis,'setting':setting,'n':len(vv),'projectile_com_ke_loss_mean_eV':m,'sample_sd_eV':sd,'reference_mean_eV':ref,'change_percent':100*(m-ref)/ref,'sd_percent_of_reference':100*sd/ref})
    write_csv('results/recomputed/Fig3_individual_runs.csv',sensitivity);write_csv('results/recomputed/Fig3_sensitivity_summary.csv',summary)
    write_csv('results/recomputed/P3_energy_drift.csv',[{'run_id':r['run_id'],'cell':r['cell'],'energy_drift_fraction':r['energy_drift_fraction'],'absolute_fraction':abs(r['energy_drift_fraction'])} for r in main_rows])
    write_csv('results/recomputed/P1_P4_raw_checks.csv',checks)
    print(json.dumps({'mechanism_raw_checks':20,'sensitivity_raw_checks':60,'structure_first_last_checks':20,'maximum_absolute_energy_drift_fraction':max(abs(r['energy_drift_fraction']) for r in main_rows),'maximum_absolute_P4_change_percent':max(abs(r['change_percent']) for r in summary)},indent=2))
if __name__=='__main__':main()
