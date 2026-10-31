"""Prepare and optionally execute one archived LAMMPS run.
Defaults to printing the command; --execute explicitly runs it.
All newly simulated outputs go to results/rerun/, preserving observations.
"""
import argparse,subprocess,shlex,json
from common import ROOT,load

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',required=True);ap.add_argument('--lmp',default='lmp')
    ap.add_argument('--stage',choices=['relax','equil','impact'],default='impact')
    ap.add_argument('--execute',action='store_true');ap.add_argument('--new-base',action='store_true')
    ap.add_argument('--check-only',action='store_true',help='Replace time integration by run 0 in a temporary input to check syntax and force field compatibility.')
    args=ap.parse_args();spec=load('in/run_specs.json')[args.run_id];camp=spec['package_campaign'];base=spec['base_id']
    output=f'results/rerun/{camp}/{args.run_id}';baseout=f'results/rerun/{camp}/structures/data.{base}'
    if args.check_only:output+='/check_'+args.stage
    (ROOT/output).mkdir(parents=True,exist_ok=True);(ROOT/baseout).parent.mkdir(parents=True,exist_ok=True)
    structure=baseout if args.new_base else f'structures/{camp}/data.{base}'
    if args.stage=='impact' and not (ROOT/structure).is_file():raise SystemExit(f'Missing base: {structure}; run --stage equil --execute first')
    wrapper=f'in/in.{args.stage}'
    variables={'relax_input':f'in/relax/{camp}/{base}.in','equil_input':f'in/equil/{camp}/{base}.in','impact_input':f'in/impact/{camp}/{args.run_id}.in','structure_file':structure,'base_output':baseout,'base_status':f'results/rerun/{camp}/structures/{base}.status','output_dir':output}
    if args.check_only:
        if args.stage=='impact':text=(ROOT/variables['impact_input']).read_text(encoding='utf-8')
        else:
            text=(ROOT/variables['relax_input']).read_text(encoding='utf-8')
            if args.stage=='equil':text+='\n'+(ROOT/variables['equil_input']).read_text(encoding='utf-8')
        text='\n'.join('run 0' if line.startswith('run ') else line for line in text.splitlines())+'\n'
        temporary=ROOT/output/'forcefield_check.in';temporary.write_text(text,encoding='utf-8')
        wrapper=str(temporary.relative_to(ROOT)).replace('\\','/')
    cmd=shlex.split(args.lmp,posix=True)+['-in',wrapper,'-log',output+'/log.lammps']
    for key,val in variables.items():cmd+=['-var',key,val]
    print(json.dumps({'working_directory':str(ROOT),'command':cmd,'execute':args.execute,'check_only':args.check_only},indent=2))
    if args.execute:subprocess.run(cmd,cwd=ROOT,check=True)
if __name__=='__main__':main()
