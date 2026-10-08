import os, sys, tarfile, glob, time
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from species import SP, TNODES
from dace_query.opacity import Molecule, Atom
R = '/orion/ptmp/jinma/rsg_wind_1008/dace/raw/'
jobs = []
for s in SP:
    name, kind, did, ll, ver, tmin, tmax, numax, m = s
    for T in TNODES:
        if tmin <= T <= tmax: jobs.append((s, int(T)))
def one(job):
    (name, kind, did, ll, ver, tmin, tmax, numax, m), T = job
    d = R + name; os.makedirs(d, exist_ok=True)
    expect = int(round(numax/0.01))*4
    fb = glob.glob(f'{d}/Out_*_{T:05d}_n800.bin')
    if fb and os.path.getsize(fb[0]) == expect: return f'{name} {T} skip'
    tf = f'{d}/t{T}.tar'
    for att in range(3):
        try:
            if kind == 'm': Molecule.download(did, ll, ver, (T, T), (-8, -8), output_directory=d, output_filename=f't{T}.tar')
            else:
                a, c = did.split(':'); Atom.download(a, int(c), ll, float(ver), (T, T), (-8, -8), output_directory=d, output_filename=f't{T}.tar')
            with tarfile.open(tf) as t: t.extractall(d)
            os.remove(tf)
            fb = glob.glob(f'{d}/Out_*_{T:05d}_n800.bin')
            sz = os.path.getsize(fb[0]) if fb else -1
            return f'{name} {T} {"OK" if sz == expect else "BADSIZE"} {sz} {expect}'
        except Exception as e:
            err = repr(e); time.sleep(5)
    return f'{name} {T} FAIL {err}'
with ThreadPoolExecutor(4) as ex:
    for r in ex.map(one, jobs): print(r, flush=True)
print('ALLDONE', len(jobs))
