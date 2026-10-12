"""Cost of the BSG reproduction run, Raven part (repro_sf_4n, 4 nodes x 4 A100) + viper part (bsg_viper_1006/repro,
8 apu nodes x 2 MI300A; continued there on 10-07 from Raven rst 00027) from their run logs and
sacct (read only).  (Old single-run version: cost_repro.py.bak_viper1007.)

Every link's log (run.log of the running link, run.MMDD_HHMMSS.log of finished links) has
'elapsed=<s> cycle=<n> time=<t> dt=<dt>' lines every ndiag = 50 cycles.  Wall of a part = sum
over its links of the last elapsed (restart overlaps included: they are real cost); sacct
(Raven over the user's ssh socket, skipped if it is not open; viper local) gives the billed
elapsed of every started link incl. start-up and replaces the log wall of a part when
available.  Projection = wall so far + (tlim - t_now) / (model-time rate over the last ~1000
cycles of the newest link): it assumes the current dt and s/cycle persist.
GPU-hours = 16 x wall per part.

usage: python3 cost_repro.py  -> prints JSON, writes cache/cost.json"""
import glob
import json
import os
import re
import subprocess

# RUN = os.environ.get('BSG_REPRO_RUN', '/raven/ptmp/jinma/bsg_raven_1002/repro_sf_4n')  (before 10-07)
PARTS = [
    dict(name='Raven', run='/raven/ptmp/jinma/bsg_raven_1002/repro_sf_4n', gpu='A100',
         ngpu=16, layout='4 Raven nodes x 4 A100', jobs=['30932910'], host='raven'),
    dict(name='viper', run='/viper/ptmp2/jinma/bsg_viper_1006/repro', gpu='MI300A',
         ngpu=16, layout='8 viper apu nodes x 2 MI300A', jobs=['12106576', '12106577', '12106578', '12106579',
                                                     '12106580'], host='local'),
    # 10-08: the Raven twin repro_sf_cont (from the viper rst 00051, 48.16 d) finished the run
    dict(name='Raven (end)', run='/raven/ptmp/jinma/bsg_raven_1002/repro_sf_cont', gpu='A100',
         ngpu=16, layout='4 Raven nodes x 4 A100', jobs=['30975887'], host='raven'),
]
# handover copy (10-09): BSG_COST_PARTS = path of a JSON list of parts (same keys as above;
# host 'local' = sacct here) replaces the viper/Raven PARTS
if os.environ.get('BSG_COST_PARTS'):
    PARTS = json.load(open(os.environ['BSG_COST_PARTS']))
SOCK = os.path.expanduser('~/.ssh/raven.sock')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'cost.json')
TLIM = 4.96e6
DAY = 86400.0
PAT = re.compile(r'^elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)')


def rows(path):
    txt = subprocess.run(['grep', '-a', '^elapsed=', path], capture_output=True,
                         text=True).stdout
    out = []
    for ln in txt.splitlines():
        m = PAT.match(ln)
        if m:
            out.append((float(m.group(1)), int(m.group(2)), float(m.group(3)),
                        float(m.group(4))))
    return out


def sacct(jobs, host):
    """[{job, state, start, elapsed_s}] of the chain's started links, or None (no access)"""
    if not jobs:
        return None
    cmd = ['sacct', '-n', '-P', '-X', '-j', ','.join(jobs), '-o',
           'JobID,State,Start,ElapsedRaw']
    if host == 'raven':
        if not os.path.exists(SOCK):
            return None
        cmd = ['ssh', '-S', SOCK, 'raven01i.mpcdf.mpg.de', ' '.join(cmd)]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=40).stdout
    except Exception:  # noqa: BLE001
        return None
    res = []
    for ln in out.splitlines():
        f = ln.split('|')
        if len(f) == 4 and f[3].isdigit() and int(f[3]) > 0:
            res.append(dict(job=f[0], state=f[1], start=f[2], elapsed_s=int(f[3])))
    return res or None


def part_links(run):
    logs = sorted(glob.glob(os.path.join(run, 'run.*.log')))
    if os.path.exists(os.path.join(run, 'run.log')):
        logs.append(os.path.join(run, 'run.log'))
    links = []
    for f in logs:
        r = rows(f)
        if r:
            links.append(dict(log=f, wall_s=r[-1][0], cycle0=r[0][1],
                              cycle1=r[-1][1], t0=r[0][2], t1=r[-1][2], dt=r[-1][3], rows=r))
    return links


def main():
    parts, all_links = [], []
    for p in PARTS:
        links = part_links(p['run'])
        jobs = sorted(set(p['jobs']) | {os.path.basename(f).split('.')[1] for f in
                                        glob.glob(os.path.join(p['run'], 'dual.*.out'))})
        sa = sacct(jobs, p['host'])
        wall_log = sum(x['wall_s'] for x in links) / 3600
        wall = sum(x['elapsed_s'] for x in sa) / 3600 if sa else wall_log
        parts.append(dict(name=p['name'], run=p['run'], gpu=p['gpu'], ngpu=p['ngpu'],
                          layout=p['layout'], jobs=jobs, n_links=len(links),
                          logs=[x['log'] for x in links], wall_log_h=wall_log,
                          sacct=sa or [], wall_h=wall, wall_src='sacct' if sa else 'log',
                          gpu_h=p['ngpu'] * wall,
                          t0_d=links[0]['t0'] / DAY if links else None,
                          t1_d=links[-1]['t1'] / DAY if links else None,
                          cycle1=links[-1]['cycle1'] if links else None,
                          s_per_cycle_mean=sum(x['wall_s'] for x in links) /
                          max(1, sum(x['cycle1'] - x['cycle0'] for x in links))))
        all_links += links
    cur = all_links[-1]
    r = cur['rows']
    k = max(0, len(r) - 21)                 # last ~1000 cycles (ndiag 50)
    dw, dc, dtm = r[-1][0] - r[k][0], r[-1][1] - r[k][1], r[-1][2] - r[k][2]
    s_cyc = dw / dc
    rate = dtm / dw                          # model s per wall s
    wall = sum(p['wall_h'] for p in parts) * 3600
    t_now = cur['t1']
    rest = max(TLIM - t_now, 0.0) / rate
    ng_now = parts[-1]['ngpu']
    res = dict(parts=parts, run=' + '.join(p['run'] for p in parts),
               jobs=sum((p['jobs'] for p in parts), []), n_links=len(all_links),
               logs=[x['log'] for x in all_links],
               wall_h=wall / 3600, wall_log_h=sum(p['wall_log_h'] for p in parts),
               t_now_d=t_now / DAY, cycle=cur['cycle1'],
               s_per_cycle_recent=s_cyc, dt_recent_s=r[-1][3],
               s_per_cycle_mean=parts[-1]['s_per_cycle_mean'],
               model_d_per_wall_h=rate * 3600 / DAY, gpu_now=parts[-1]['gpu'],
               proj_total_h=(wall + rest) / 3600, gpu_h=sum(p['gpu_h'] for p in parts),
               proj_gpu_h=sum(p['gpu_h'] for p in parts) + ng_now * rest / 3600,
               ngpu=ng_now, tlim_d=TLIM / DAY)
    res['sacct'] = sum((p['sacct'] for p in parts), [])
    json.dump(res, open(OUT, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('parts', 'sacct')},
                     indent=1))
    for p in parts:
        print(p['name'], p['gpu'], '%.2f h (%s)' % (p['wall_h'], p['wall_src']),
              '%.0f GPU-h' % p['gpu_h'], 't %.2f-%.2f d' % (p['t0_d'], p['t1_d']),
              'jobs', p['jobs'])


if __name__ == '__main__':
    main()
