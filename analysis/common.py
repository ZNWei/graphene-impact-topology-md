"""Portable package paths and standard-library data readers."""
from pathlib import Path
import csv,json,math,hashlib

ROOT=Path(__file__).resolve().parents[1]
ORDER=['R0L0','R0L1','R1L0','R1L1']
def load(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))
def write_csv(path, rows):
    rows=list(rows);p=ROOT/path;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def cell(row):return f"R{row['cross_rebo']}L{row['cross_lj']}"
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def track(path,columns):
    last=None
    for l in path.read_text(encoding='utf-8').splitlines():
        if l.strip() and not l.startswith('#'):last=dict(zip(columns,map(float,l.split())))
    if last is None:raise ValueError(f'No track records: {path}')
    return last
def dump_frames(path):
    with path.open(encoding='utf-8') as f:
        while True:
            s=f.readline()
            if not s:return
            assert s.strip()=='ITEM: TIMESTEP'
            step=int(f.readline());assert f.readline().startswith('ITEM: NUMBER OF ATOMS')
            n=int(f.readline());assert f.readline().startswith('ITEM: BOX BOUNDS')
            for _ in range(3):f.readline()
            cols=f.readline().split()[2:];ix={k:i for i,k in enumerate(cols)};atoms=[]
            for _ in range(n):
                v=f.readline().split();atoms.append((int(v[ix['id']]),int(v[ix['type']]),*(float(v[ix[k]]) for k in ('x','y','z'))))
            yield step,atoms
def terminal_counts(path):
    active=False;np=nt=nd=0
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('Atoms'):
            assert '# full' in line;active=True;continue
        if not active or not line.strip():continue
        s=line.split()
        if not s[0].isdigit():break
        typ=int(s[2]);z=float(s[6])
        if typ==4:np+=1;nd+=z<=-500
        elif typ in (1,2,3):nt+=1
    return np,nt,nd
def undercoordination(first,last,spec):
    """Same-layer neighbors within 1.85 A, with edge IDs frozen from frame 0.
    A spatial hash evaluates the same non-periodic Euclidean rule as the
    archived cKDTree implementation. Only first/last frames are needed here.
    """
    cfg=spec['baseline'];cut=cfg['bond_distance_cutoffs_angstrom'][1]
    low=(cfg['box_xy_angstrom']-cfg['sheet_size_angstrom'])/2
    high=(cfg['box_xy_angstrom']+cfg['sheet_size_angstrom'])/2;edge=cfg['edge_width_angstrom']
    edge_ids={i for i,t,x,y,z in first if t in (1,2,3) and (x<=low+edge or x>=high-edge or y<=low+edge or y>=high-edge)}
    def damage(atoms):
        total=0
        for layer in (1,2,3):
            aa=[a for a in atoms if a[1]==layer];grid={}
            for a in aa:
                key=tuple(math.floor(v/cut) for v in a[2:]);grid.setdefault(key,[]).append(a)
            for i,t,x,y,z in aa:
                if i in edge_ids:continue
                key=tuple(math.floor(v/cut) for v in (x,y,z));n=0
                for dx in (-1,0,1):
                    for dy in (-1,0,1):
                        for dz in (-1,0,1):
                            for j,_,xx,yy,zz in grid.get((key[0]+dx,key[1]+dy,key[2]+dz),[]):
                                if i!=j and (x-xx)**2+(y-yy)**2+(z-zz)**2<=cut*cut:n+=1
                total+=n<2.5
        return total
    return damage(last)-damage(first)
