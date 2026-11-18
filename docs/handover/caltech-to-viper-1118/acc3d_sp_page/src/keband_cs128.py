#!/usr/bin/env python3
"""one-off cache: accsp_1102/tools/keband.py cs on every test128wb hydro_w dump (static + stream, finished run) -> data/cs128_keband.json"""
import glob, json, os, re, subprocess, sys
J = '/resnick/groups/carnegie_poc/jingze'; R = J + '/csacc_1031/test128wb'; OUT = sys.argv[1]
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
fs = sorted(glob.glob(R + '/static/bin/ryper.hydro_w.*.bin')) + [f for f in sorted(glob.glob(R + '/stream/bin/ryper.hydro_w.*.bin')) if not os.path.islink(f)]
for f in fs:
    if f in res: continue
    o = subprocess.run([sys.executable, J + '/accsp_1102/tools/keband.py', 'cs', f], capture_output=True, text=True).stdout
    lat, rad, tot, t = [], [], None, None
    for ln in o.splitlines():
        m = re.match(r'\s*(\d+)-\s*(\d+) deg\s+(.*)\|\s*([0-9.eE+-]+)', ln)
        if m: lat.append([float(x) for x in m.group(3).split()]); rad.append(float(m.group(4)))
        m = re.search(r't ([0-9.]+) KE total ([0-9.eE+-]+) radial ([0-9.eE+-]+) lateral ([0-9.eE+-]+)', ln)
        if m: t = float(m.group(1)); tot = [float(m.group(i)) for i in (2, 3, 4)]
    res[f] = dict(t=t, lat=lat, rad=rad, tot=tot); print('cs128 keband', os.path.basename(f), t, tot, flush=True)
    json.dump(res, open(OUT + '.tmp', 'w')); os.replace(OUT + '.tmp', OUT)
