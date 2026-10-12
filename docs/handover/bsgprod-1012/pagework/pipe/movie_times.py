"""movie_times.py PAGEDIR RUN: write PAGEDIR/movie/n.json = {n, t[d] of each f4 frame's dump}."""
import glob
import json
import os
import re
import sys

sp, run = sys.argv[1], sys.argv[2]
fr = sorted(glob.glob(os.path.join(sp, 'movie', 'f4', 'f4_*.png')))
t = []
for f in fr:
    n = int(os.path.basename(f)[3:8])
    with open(os.path.join(run, 'bin', 'bsg3d.hydro_w.%05d.bin' % n), 'rb') as fh:
        h = fh.read(400).decode('latin-1')
    t.append(round(float(re.search(r'time=\s*([-+0-9.eE]+)', h).group(1)) / 86400, 3))
json.dump(dict(n=len(fr), t=t), open(os.path.join(sp, 'movie', 'n.json'), 'w'))
print('n.json', len(fr), t[0], t[-1])
