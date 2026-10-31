"""Verify supplied file hashes and input/structure coverage without running MD."""
import csv,json,re
from common import ROOT,load,sha

def main():
    failures=[];count=0
    for line in (ROOT/'CHECKSUMS.sha256').read_text(encoding='utf-8').splitlines():
        h,name=line.split('  ',1);p=ROOT/name
        if not p.is_file() or sha(p)!=h:failures.append(name)
        count+=1
    specs=load('in/run_specs.json');assert len(specs)==230
    for rid,r in specs.items():
        for p in [ROOT/f"in/impact/{r['package_campaign']}/{rid}.in",ROOT/f"structures/{r['package_campaign']}/data.{r['base_id']}"]:
            if not p.is_file():failures.append(str(p.relative_to(ROOT)))
    print(json.dumps({'files_hashed':count,'formal_runs':len(specs),'initial_structures':len(list((ROOT/'structures').glob('*/data.*'))),'integrity_failures':failures,'round2_raw_outputs_available':False},indent=2))
    if failures:raise SystemExit(1)
if __name__=='__main__':main()
