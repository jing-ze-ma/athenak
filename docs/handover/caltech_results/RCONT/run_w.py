import json,subprocess,os,glob
from multiprocessing import Pool
import sys
J='/resnick/groups/carnegie_poc/jingze/'; star=sys.argv[1]
o=json.load(open(J+'rcont_1124/out/%s.json'%star)); OP=J+'fresh_prod_1120/tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'
if star=='bsg': OP=J+'fresh_prod_1120/bsg/rosseland_mesakap_gs98_x0.70_z0.008.txt'
def one(x):
    f=x['f']; h=glob.glob(os.path.dirname(os.path.dirname(f))+'/*.user.hst')[0]
    r=subprocess.run([sys.executable,J+'rcont_1124/rph_wedge.py',J+'rsgE_page/vispy_4de8168d',OP,h,f],capture_output=True,text=True)
    return x['tag'],f,(r.stdout.strip() or 'ERR '+r.stderr[-300:])
if __name__=='__main__':
    wf=J+'rcont_1124/out/%s_w.json'%star; old={f:(a,t) for a,f,t in json.load(open(wf))} if os.path.exists(wf) else {}
    todo=[x for x in o if x['f'] not in old or 'ERR' in old[x['f']][1]]
    with Pool(int(os.environ.get('NPROC','4'))) as p: new={f:(a,t) for a,f,t in p.map(one,todo,chunksize=1)}
    old.update(new); out=[[x['tag'],x['f'],old[x['f']][1]] for x in o]
    json.dump(out,open(J+'rcont_1124/out/%s_w.json'%star,'w'),indent=0); print(star,len(out),sum('ERR' in z[2] for z in out),out[0][2][:260],flush=True)
