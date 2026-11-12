import argparse, json
from pathlib import Path
import numpy as np
from .q import score
from .propagation import propagate, height_slices, intensity_slices

def load_cfg(path):
    import yaml
    with open(path) as f:return yaml.safe_load(f)
def q_main(c):
    data_path=c.get('source_npz',c.get('data'))
    if not data_path: raise ValueError('Q config requires source_npz (or data)')
    d=np.load(data_path,allow_pickle=False); keys=['F','Fz','Fzz','A','Az','Azz','weights','region','energies']
    missing=[k for k in keys if k not in d]
    if missing: raise ValueError('Q data missing: '+','.join(missing))
    val,detail=score(d['F'],d['Fz'],d['Fzz'],d['A'],d['Az'],d['Azz'],d['weights'],d['region'],d['energies'],c['vacuum_reference_eV'],c['distance_A'],c.get('length_A'))
    out={'Q':float(val),'distance_A':c['distance_A'],'vacuum_reference_eV':c['vacuum_reference_eV'],'source_npz':str(data_path),'components':{k:float(np.sqrt(np.mean(v*v))) for k,v in detail.items() if k in {'source','slope','curvature','residual'} and isinstance(v,np.ndarray)}}
    if c.get('cp2k_outputs') is not None: out['cp2k_outputs']=c['cp2k_outputs']
    Path(c.get('output','q_result.json')).write_text(json.dumps(out,indent=2))
def vp_main(c):
    b=np.load(c['boundary_npy']); u=np.load(c['potential_npy']); z=np.load(c['z_npy'])
    psi,meta=propagate(b,u,c['energy_eV'],c['spacing_A'],c.get('tolerance',1e-8),c.get('maxiter',1000),c.get('allow_allowed',False))
    slices=height_slices(psi,z,c['heights_A'])
    out=c.get('output_npz','vp_result.npz'); np.savez_compressed(out,heights_A=np.asarray(c['heights_A']),psi=slices,intensity=np.abs(slices)**2)
    meta=dict(meta)
    meta.update({'energy_eV': c['energy_eV'], 'spacing_A': c['spacing_A'],
                 'heights_A': list(c['heights_A'])})
    if c.get('cp2k_outputs') is not None: meta['cp2k_outputs']=c['cp2k_outputs']
    Path(c.get('metadata_json','vp_metadata.json')).write_text(json.dumps(meta,indent=2))
def main():
    p=argparse.ArgumentParser();p.add_argument('config');a=p.parse_args();c=load_cfg(a.config)
    if c['mode']=='q':q_main(c)
    elif c['mode']=='vp':vp_main(c)
    else: raise ValueError('mode must be q or vp')
if __name__=='__main__':main()
