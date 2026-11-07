#!/usr/bin/env python3
# USER 11-06: re-encode only the segments holding re-rendered stale frames (mk_movies treats an unchanged idx list as nothing to do),
# rejoin the parts they belong to; same _enc/_join ffmpeg settings as mk_movies.py. Also confirms stale vs fresh data differ.
import json, os, shutil, subprocess, numpy as np
PG = '/resnick/groups/carnegie_poc/jingze/hgglob_page/out'; FD = PG + '/frames'; MV = PG + '/publish/movies'; DST = '/resnick/home/jingze/ATHENAK/pages/hggl_glob/movies'
FF = os.environ['FFMPEG']
stale = open(PG + '/stale_1106/list.txt').read().split()
for i in stale:
    a = np.load(f'{PG}/stale_1106/red_{i}.npz'); b = np.load(f'{PG}/data/red_{i}.npz')
    d = np.max(np.abs(a['rho'] / b['rho'] - 1)); print(i, 't', float(a['t']) == float(b['t']), 'max|rho shell ratio-1| %.3e' % d, 'Rph_med', float(a['Rph_med']), float(b['Rph_med']), flush=True)
out = []
for name, pre in (('eq', 'hggl_eq'), ('inner', 'hggl_inner'), ('phot', 'hggl_phot'), ('map', 'hggl_map')):
    SD = f'{FD}/{name}_segs'; pj = SD + '/parts.json'; st = json.load(open(pj))
    shutil.copy(pj, pj + '.pre_stale_1106')
    idx = st['idx']; pos = {x: k for k, x in enumerate(idx)}
    for p in st['parts']:
        hit = False
        for g in p['segs']:
            a, b = pos[g['i0']], pos[g['i1']] + 1
            assert b - a == g['n']
            if not any(g['i0'] <= s <= g['i1'] for s in stale):
                continue
            hit = True
            dst = os.path.join(SD, g['f']); tmp = dst + '.new.mp4'; td = dst + '.in'
            shutil.rmtree(td, ignore_errors=True); os.makedirs(td)
            for k, x in enumerate(idx[a:b]):
                os.symlink(f'{FD}/{name}/{x}.png', f'{td}/{k:06d}.png')
            subprocess.run([FF, '-y', '-loglevel', 'error', '-framerate', str(st['fps']), '-i', td + '/%06d.png', '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
                            '-c:v', 'libx264', '-threads', '4', '-pix_fmt', 'yuv420p', '-crf', str(st['crf']), '-preset', 'slow', '-movflags', '+faststart', tmp], check=True)
            shutil.rmtree(td); os.replace(tmp, dst); g['mb'] = os.path.getsize(dst) / 1e6
            print(name, 're-encoded', g['f'], g['i0'], g['i1'], flush=True)
        if hit:
            f = f"{MV}/{p['file']}"; cl = f + '.join.txt'
            with open(cl, 'w') as fl:
                for g in p['segs']:
                    fl.write("file '%s'\n" % os.path.join(SD, g['f']))
            subprocess.run([FF, '-y', '-loglevel', 'error', '-f', 'concat', '-safe', '0', '-i', cl, '-c', 'copy', '-movflags', '+faststart', f + '.tmp.mp4'], check=True)
            os.replace(f + '.tmp.mp4', f); os.remove(cl); p['mb'] = os.path.getsize(f) / 1e6
            shutil.copy(f, f"{DST}/{p['file']}"); out.append(p['file']); print('part', p['file'], '%.1f MB' % p['mb'], flush=True)
    json.dump(st, open(pj + '.tmp', 'w')); os.replace(pj + '.tmp', pj)
open(PG + '/stale_1106/publish_movies.txt', 'w').write(''.join(f'movies/{f}\n' for f in out)); print('PUBLISH', out)
