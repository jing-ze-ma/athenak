"""IDORT results page: read the FIN3 final runs (fin3-1012 281dd686 = 1fdbe56e/914e5433 binaries + implicit_predictor_frozen
false; IN's run trees fin3_1012/runs_281dd686 (Cartesian, final-config reruns) and sp_FIN3_281dd686 (sp)) through the shadow
trees idort_page/fin3runs (mklinks_fin3.py: FIN_<x> -> FIN3_<x>, PRODd_<x> -> F2O_<x> = FIN2 control) and cache what the figures need.
(was: final_1012 FIN runs, method-1012 78b2c961; backup figsrc.bak_fin3)
Run inside sbatch (heavy reading).  usage: extract.py [item ...]   (default: all)
Each item writes cache/<item>.pkl (a dict) or cache/<item>.err on failure.  No simulation is run."""
import glob
import json
import math
import os
import pickle
import sys
import traceback

import numpy as np

F = '/raven/ptmp/jinma/suiteF_1012'
RD = '/raven/ptmp/jinma/idort1020/page/p19runs/cart/'    # FIN3 final (names FIN_*)
RDC = RD   # comparison arm: PRODd_* -> F2O (FIN2 control, suiteF runs_js)
CACHE = '/raven/ptmp/jinma/idort1020/page/cache'
os.environ['SPT_RUNS'] = '/raven/ptmp/jinma/idort1020/page/p19runs/sp'
sys.path.insert(0, F + '/ana')
sys.path.insert(0, F + '/sp/ana')
sys.path.insert(0, '/raven/u/jinma/ATHENAK/wt-prodfinal/vis/python')


def log_status(d):
    """rc, Picard mean / NONCONV, FATAL line from a Cartesian run dir"""
    out = {}
    try:
        out['rc'] = open(d + '/rc.txt').read().strip()
    except Exception:
        pass
    try:
        txt = open(d + '/run.log', errors='ignore').read()
        if 'Picard iterations mean=' in txt:
            s = txt.split('Picard iterations mean=')[-1]
            out['picard'] = float(s.split()[0])
            if 'NON-CONVERGED=' in s:
                out['nonconv'] = float(s.split('NON-CONVERGED=')[1].split()[0])
        for ln in txt.splitlines():
            if 'FATAL' in ln or 'DIVERGED' in ln:
                out['fatal'] = ln.strip()[:300]
                break
    except Exception:
        pass
    return out


# ------------------------------------------------------------------ Cartesian
def cb():
    import cb as C
    from lcommon import load
    out = {}
    for run in ('cb/FIN_r4', 'cb/FIN_r12', 'cb/FIN_r4c0.5', 'cb/FIN_r4c32', 'cb/FIN_r4n4x', 'cbh/FIN_r4'):
        d = RD + run
        try:
            v = load(d)
        except Exception:   # idort-1020: a FATAL run (no output) is reported by its status only
            out[run] = dict(st=log_status(d), missing=True); continue
        p = C.readin(d + '/in.athinput')
        ex = C.exact(p, v['x1'], v['x2'])
        m, _ = None, None
        try:
            m = C.metrics(d)
        except Exception:
            m = None
        out[run] = dict(E=v['E'], J=v.get('J'), ex=ex, x1=v['x1'], x2=v['x2'], t=v['t'],
                        rowE=v['E'].sum(axis=0), rowex=ex.sum(axis=0), metrics=m, st=log_status(d),
                        beams=C.beams(p), x1min=float(p['mesh/x1min']), x2min=float(p['mesh/x2min']),
                        x2max=float(p['mesh/x2max']), x1max=float(p['mesh/x1max']))
    return out


def cbg():
    import cb as C
    from lcommon import load
    out = {'PRODd': log_status(RDC + 'cbg/PRODd_g8x32')}
    d = RD + 'cbg/FIN_g8x32'
    v = load(d)
    p = C.readin(d + '/in.athinput')
    ex = C.exact(p, v['x1'], v['x2'])
    out['FIN'] = dict(E=v['E'], J=v.get('J'), ex=ex, x1=v['x1'], x2=v['x2'], t=v['t'],
                      rowE=v['E'].sum(axis=0), rowex=ex.sum(axis=0), st=log_status(d))
    return out


def shd():
    import cb as C
    import shd as S
    from lcommon import load
    out = {}
    for run in ('shd/FIN_r12g2', 'shd/FIN_r36g2', 'shdg/FIN_g8x32'):
        d = RD + run
        v = load(d)
        p = C.readin(d + '/in.athinput')
        ex = S.exact(p, v['x1'], v['x2'])
        out[run] = dict(E=v['E'], J=v.get('J'), ex=ex, x1=v['x1'], x2=v['x2'], t=v['t'], st=log_status(d))
    return out


def hcyl():
    import hcyl as H
    import bin_convert as bc
    out = {}
    for run in ('FIN_g8x32', 'FIN_n512'):
        d = RD + 'hcyl/' + run
        if not os.path.isdir(d):
            continue
        o, f = H.metrics(d)
        r = bc.read_binary_as_athdf(H.lastbin(d, 'm1'))
        out[run] = dict(metrics=o, prof=f, E2=r['m1_e'][0], x1=r['x1v'], x2=r['x2v'], st=log_status(d))
        rv = H.lastbin(d, 'm1_vet')
        if rv:
            q = bc.read_binary_as_athdf(rv)
            if 'm1_d11' in q:
                j0 = np.argmin(np.abs(r['x2v']))
                j1 = j0 - 1 if r['x2v'][j0] > 0 else j0 + 1
                m = r['x1v'] > 0
                out[run]['fxx'] = 0.5*(q['m1_d11'][0][j0] + q['m1_d11'][0][j1])[m]
    ex = np.load(F + '/ana/hcyl_exact_512.npz')
    out['exact512'] = {k: ex[k] for k in ex.files}
    return out


def nlte():
    from s2 import nlte as N
    return {l: N(RD + 'nltek%d/FIN_d' % l) for l in range(1, 9)}


def dif():
    from s2 import dif as D
    out = {'sdif': D(RD + 'sdifk/FIN_d', 'sdif')}
    sys.path.insert(0, F + '/ana/mq')
    import ddts
    for k in ('ddf', 'ddn', 'ddiff'):
        d = RD + '%s/FIN_d' % k
        if not os.path.isdir(d):
            continue
        rows, Dc = ddts.ana(d)
        prof = []
        for f in sorted(glob.glob(d + '/tab/*.m1.*.tab')):
            t, m = ddts.rd(f)
            if t <= 0:
                continue
            x, E = m['x1v'], m['m1_e']
            s = 1 + 160*Dc*t
            xx = (x - ddts.V*t + ddts.L/2) % ddts.L - ddts.L/2
            prof.append(dict(t=t, x=x, E=E, ex=np.exp(-40*xx**2/s)/np.sqrt(s) + ddts.BG))
        out[k] = prof
    for k in ('ddf', 'ddn'):
        if not os.path.isdir(RD + '%s/FIN_d' % k):
            continue
        rows, Dc = ddts.ana(RD + '%s/FIN_d' % k)
        out[k + '_ts'] = dict(rows=rows, D=Dc)
    return out


def eq():
    from s2 import eq_table
    sys.path.insert(0, F + '/ana/mq')
    import eqts
    out = {}
    for c in ('eq1', 'eq2'):
        out[c] = eq_table(RD + '%s/FIN_fine' % c)
        out[c + '_st'] = log_status(RD + '%s/FIN_fine' % c)
    for a in ('FIN', 'PRODd'):
        out['eq3t_' + a] = eqts.series((RD if a == 'FIN' else RDC) + 'eq3t/%s_t' % a)
    return out


def rw():
    rows = {}
    # FIN3 page: the suite's rw cases rerun with FIN3 keys (reruns/fin3/rw); no comparison arm
    for tree, want in (('/raven/ptmp/jinma/idort_page/reruns/fin3/rw', 'FIN'),):
        for f in glob.glob(tree + '/*/results.json'):
            r = json.load(open(f))
            c = r['case']
            if want == 'PRODd' and c['extra'].get('rm_implicit_gas_src') == 'matched':
                continue
            key = (round(3*c['prat'], 6), c['tau'], want)
            rows.setdefault(key, {})[c['N']] = dict(l1_rho=r['l1_rho'], l1_E=r['l1_E'], w_re=r['w_re'], w_im=r['w_im'],
                                                    ref_key=r['ref_key'], refs=r['refs'], x1=r['x1'],
                                                    rho_end=r['rho_end'], e_end=r['e_end'], t=r['t'], absZ=r['absZ'])
    return rows


def shock():
    import shockcmp as sc
    out = {}
    for m in ('1.2', '2', '3'):
        R = sc.ref(m)
        out['ref' + m] = dict(x=R['x'], rho=R['rho'], tg=R['tg'], tr=R['tr'], xs=sc.jump(R['x'], R['rho']))
        for n in ('n512', 'n1024'):
            d = RD + 'rs_m%s/FIN_%s' % (m, n)
            if not os.path.isdir(d):
                continue
            t, x, rho, tg, tr = sc.run(d)
            out['%s_%s' % (m, n)] = dict(t=t, x=x, rho=rho, tg=tg, tr=tr, xs=sc.jump(x, rho), st=log_status(d))
    return out


def cloud():
    """Q2 (NOTE-check.md): FIN with the littests beam (200 a T0^4 per beam, a 2.7x too weak reading of IDORT Fig. 5's
    saturated colour bar) and the Q2 rerun with IDORT's beam I = 1031.3 c a -> 4 pi I w / c = 540 a T0^4 per beam
    (w = 1/24, the in-plane weight of ring 12x2 = IDORT's 24 directions); IDORT's own snapshots (Zenodo shadow_2d)"""
    import bin_convert as bc
    out = {}
    for a, d in (('FIN', RD + 'cloud/FIN_r12g2'), ('FINb540', '/raven/ptmp/jinma/idort1020/runs/xtra/cloud_b540/p')):
        fs = sorted(glob.glob(d + '/bin/*.hydro_w.*.bin'))
        if len(fs) < 9:   # FIN3 page: arm not (yet) run
            continue
        for k in (0, 3, 5, 8):
            v = bc.read_binary_as_athdf(fs[k])
            dd = v['dens'][0]
            m = dd > 3.0
            xc = float((v['x1v'][None, :]*dd*m).sum()/(dd*m).sum()) if m.any() else float('nan')
            j = int(np.argmin(np.abs(v['x2v'])))
            row = dd[j-1:j+1].mean(0)
            x = v['x1v']
            out['%s_%d' % (a, k)] = dict(dens=dd, x1=x, x2=v['x2v'], t=v['Time'], xc=xc,
                                        halo=float(x[np.argmax(row > 1.08)]),
                                        r3=(float(x[m.any(0)].min()), float(x[m.any(0)].max())), rmax=float(dd.max()))
        out['st_' + a] = log_status(d)
    out['st'] = out.get('st_FINb540', out.get('st_FIN'))
    import matplotlib.image as mi
    D = '/raven/ptmp/jinma/suite_1011/ext/idort_mini/x/IDORT_mini/shadow_2d/ds_selected/'
    out['idort_png'] = {k: mi.imread(D + 'irradiation_%03d.png' % k) for k in (3, 5, 8)}
    # the same measures read off IDORT's PNGs along y = 0 (jet, log 1-40 colour bar inverted by nearest colour;
    # pixel calibration of the 991x497 PNG: x = -0.5 / 0.5 at px 132 / 766, y = 0 at px 218, bar 40 / 1 at px 12 / 428)
    meas = {}
    for k, im in out['idort_png'].items():
        im = im[..., :3]
        ys = np.arange(14, 426)
        cb = im[ys, 880]
        vals = 10**(np.log10(40.0)*(1 - (ys - 12)/416.0))
        xs = np.linspace(-0.49, 0.49, 197)
        row = np.array([vals[((cb - im[218, int(132 + (x + 0.5)*634)])**2).sum(1).argmin()] for x in xs])
        meas[k] = dict(halo=float(xs[np.argmax(row > 1.08)]), r3=(float(xs[row > 3].min()), float(xs[row > 3].max())),
                       rmax=float(row.max()))
    out['idort_meas'] = meas
    return out


def bub():
    """Q2: periodic RETARDED free-streaming reference (bubref.py: first-disc-hit, emission inside the light cone
    t - s/c >= 0, periodic images; 2-D in-plane rays E2 for the rings, 3-D directions E3 for the Gauss set), the steady
    reference for comparison, and the metrics of bubeval.py (FIN and PRODd)"""
    import bub as B
    out = {}
    R = np.load(CACHE + '/bubref_2.npz')
    out['ref'] = {k: R[k] for k in R.files}
    out['eval'] = json.load(open(CACHE + '/bubeval_2.json'))
    for run in ('FIN_r100k5', 'FIN_r100k3', 'FIN_g8x32k3'):
        o, fields = B.metrics(RD + 'bub/' + run, ref='E3')
        out[run] = dict(E={t: fields[t][0] for t in fields})
    return out


def marsh():
    ref = np.loadtxt('/raven/u/jinma/ATHENAK/wt-prodfinal/tests_m1/runs_3a/t6_ref_sn.txt')
    out = {'ref': ref}
    for a in ('FIN_c1', 'FIN_c10'):
        d = RD + 'marsh/' + a
        fs = sorted(glob.glob(d + '/tab/*.m1.*.tab'))
        res = []
        for f in fs:
            hdr = [l for l in open(f) if l.startswith('#')][-1][1:].split()
            arr = np.loadtxt(f, ndmin=2)
            hdr = [h for h in hdr if h not in ('k', 'x3v', 'j', 'x2v')] if len(hdr) > arr.shape[1] else hdr
            t = float(open(f).readline().split('time=')[1].split()[0])
            res.append(dict(t=t, x=arr[:, hdr.index('x1v')], E=arr[:, hdr.index('m1_e')]))
        out[a] = dict(tabs=res, st=log_status(d))
    return out


def pipe():
    import pipe as P
    import bin_convert as bc
    out = {}
    for a in ('FIN_c220', 'FIN_c2200'):
        d = RD + 'pipe/' + a
        fs = sorted(glob.glob(d + '/bin/*.hydro_w.*.bin'))
        rows = []
        for f in fs:
            r = bc.read_binary_as_athdf(f)
            T = P.GM1*r['eint'][0]/r['dens'][0]
            x, y = r['x1v'], r['x2v']
            vals = []
            for (px, py) in P.PROBES:
                i = np.argmin(np.abs(x - px))
                j = np.argmin(np.abs(y - py))
                jj = j - 1 if y[j] > py else j + 1
                vals.append(0.5*(T[j, i] + T[jj, i]))
            rows.append([r['Time']] + vals)
        out[a] = dict(probes=np.array(rows), st=log_status(d))
        if a == 'FIN_c220':
            r = bc.read_binary_as_athdf(fs[-1])
            out['map'] = dict(T=P.GM1*r['eint'][0]/r['dens'][0], dens=r['dens'][0], x1=r['x1v'], x2=r['x2v'], t=r['Time'])
            mf = sorted(glob.glob(d + '/bin/*.m1.*.bin'))
            if mf:
                q = bc.read_binary_as_athdf(mf[-1])
                out['map']['E'] = q['m1_e'][0]
    return out


def g5():
    out = {}
    for k in ('12.8', '128', '1280', '12800'):
        for a, run in (('FIN', 'run_fin'),):   # idort-1020: PROD only, no old-binary overlay
            d = F + '/g5/%s/g5_%s_k%s' % (run, a, k)
            if a == 'FIN':   # FIN3 final (IN, fin3_1012/g5, 281dd686)
                d = '/raven/ptmp/jinma/idort1020/runs/xtra/g5_k%s/p' % k
            if os.path.exists(d + '/sigma2.txt'):
                out['%s_%s' % (a, k)] = dict(s=np.loadtxt(d + '/sigma2.txt'), res=open(d + '/result.txt').read())
    return out


# ------------------------------------------------------------------ spherical (vet_gd, runs_fin, arm HFfZ = FIN)
def sp_core():
    import spl
    import s3_eval as S
    out = {}
    for case in ('core_u', 'core_L16_u', 'thin_u', 'thin_L16_u'):
        for n in (64, 128, 256):
            d = S.run_dir(case, 'HFfZ', n)
            t = spl.tab(spl.tabs(d)[-1])
            r, _ = S.rcen(case, n)
            o = dict(r=r, E=t['m1_e'], F=t['m1_f1'], cyc=t['cycle'], st=spl.logstats(d))
            if case.startswith('core'):
                J, H = S.core_exact(r)
                o.update(J=J, H=H)
            out['%s_%d' % (case, n)] = o
    rr = np.exp(np.linspace(np.log(0.2), np.log(7.0), 600))
    J, H = S.core_exact(rr)
    out['core_fine'] = dict(r=rr, J=J, H=H)
    return out


def sp_nlte():
    import spl
    out = {}
    for case, arm in (('nlte2_r1000_tb', 'HFfSgZ'), ('nlte2_tb_cart', 'HCFfZ')):
        d = os.path.join(spl.RUNS, case, '%s_320' % arm)
        t = spl.tab(spl.tabs(d)[-1])
        out[case] = dict(E=t['m1_e'], st=spl.logstats(d), t=t['time'])
    R = np.load(F + '/sp/ref/nlte2_ref.npz')
    out['ref'] = {k: R[k] for k in R.files}
    return out


def sp_shd():
    import spl
    import spt
    import s3_eval as S
    d = S.run_dir('shd', 'HFfZ', 32)
    fs = sorted(f for f in os.listdir(os.path.join(d, 'bin')) if '.m1.' in f)
    b = spl.binall(os.path.join(d, 'bin', fs[-1]))
    E = b['m1_e']
    n3, n2, n1 = E.shape
    c = spt.case_cfg('shd')
    rf = spl.rfaces(c['r'][0], c['r'][1], n1, c['poly'])
    rc = spl.rcent(rf)
    tc = c['th'][0] + (np.arange(n2) + 0.5)*(c['th'][1] - c['th'][0])/n2
    pc = c['ph'][0] + (np.arange(n3) + 0.5)*(c['ph'][1] - c['ph'][0])/n3
    pr = c['prob']
    cr, ct, cp, ca = (float(pr['clump_r']), float(pr['clump_th']), float(pr['clump_ph']), float(pr['clump_a']))
    cen = cr*np.array([math.sin(ct)*math.cos(cp), math.sin(ct)*math.sin(cp), math.cos(ct)])
    k = int(np.argmin(np.abs(pc - cp)))
    ph = pc[k]
    ex = np.zeros((n2, n1))
    free = np.zeros((n2, n1))
    ins = np.zeros((n2, n1), bool)
    for j in range(n2):
        u = np.array([math.sin(tc[j])*math.cos(ph), math.sin(tc[j])*math.sin(ph), math.cos(tc[j])])
        for i in range(n1):
            ex[j, i], free[j, i] = S.shd_exact(rc[i]*u, cen, ca, ns=40)
            # Q2: shd_exact counts the star as visible from INSIDE the clump (its ray test needs the hit in front);
            # the clump is opaque and cold (E ~ 0): mark those cells, the plot masks them
            ins[j, i] = np.linalg.norm(rc[i]*u - cen) < ca
    ex[ins] = 0.0
    ev = S.ev_shd('shd', 'HFfZ', 32)
    # Q2 diagnostic: the same FIN keys with 432 vet_gd directions (18 x 24) instead of 120 (reruns/sp_shd_FIN_Q432)
    dq = '/raven/ptmp/jinma/idort1020/runs/xtra/shd432/p'   # FIN3 deck + nmu 18 nphi 24 (was sp_shd_FIN_Q432, old keys)
    fq = sorted(f for f in os.listdir(os.path.join(dq, 'bin')) if '.m1.' in f) if os.path.isdir(os.path.join(dq, 'bin')) else []
    Eq = spl.binall(os.path.join(dq, 'bin', fq[-1]))['m1_e'] if fq else None
    return dict(E=E[k], E432=None if Eq is None else Eq[k], ex=ex, free=free, ins=ins, rc=rc, rf=rf, tc=tc, th=c['th'], ph=ph, cen=(cr, ct, cp, ca),
                txt=ev['txt'], st=spl.logstats(d), t=b['time'])


def sp_pinc():
    import spl
    import spt
    X = np.load(F + '/sp/ref/pinc_exact_64.npz')
    d = os.path.join(spl.RUNS, 'pinc', 'HFfZ_64')
    fs = sorted(f for f in os.listdir(os.path.join(d, 'bin')) if '.m1.' in f)
    b = spl.binall(os.path.join(d, 'bin', fs[-1]))
    E, F1 = b['m1_e'], b['m1_f1']
    c = spt.case_cfg('pinc')
    cp = float(c['prob']['hole_ph'])
    pc = X['pc']
    k = int(np.argmin(np.abs(pc - cp)))
    tc = X['tc']
    dth = tc[1] - tc[0]
    dph = pc[1] - pc[0]
    dom = (np.cos(tc - dth/2) - np.cos(tc + dth/2))[None, :]*dph*np.ones((len(pc), 1))
    rc = X['rc']
    Lc = np.array([np.sum(F1[:, :, i]*dom)*rc[i]**2 for i in range(len(rc))])
    Lx = np.array([np.sum(X['H'][:, :, i]*dom)*rc[i]**2 for i in range(len(rc))])
    return dict(E=E[k], J=X['J'][k], rc=rc, tc=tc, th=c['th'], Lc=Lc, Lx=Lx, st=spl.logstats(d), t=b['time'],
                rf=spl.rfaces(c['r'][0], c['r'][1], len(rc), c['poly']))


def sp_dif():
    import spl
    import s3_eval as S
    out = {}
    for n in (64, 128, 256):
        d = S.run_dir('dif_r1', 'HFfZ', n)
        r, c = S.rcen('dif_r1', n)
        r0 = float(c['prob']['r0'])
        eo = float(c['prob']['e_out'])
        Dc = float(c['rad']['c_light'])/(3.0*float(c['rad']['kappa_s']))
        prof = []
        for f in spl.tabs(d):
            t = spl.tab(f)
            ex = S.dif_exact(r, t['time'], r0, eo, Dc) if t['time'] > 0 else None
            prof.append(dict(t=t['time'], E=t['m1_e'], ex=ex, sq=math.sqrt(Dc*t['time'])))
        out[n] = dict(r=r, prof=prof, st=spl.logstats(d))
    return out


def sp_relax():
    import spl
    import s3_eval as S
    out = {}
    for case in ('relax1_early_T', 'relax2_early_T'):
        for n in (1, 2, 4, 8, 16):
            try:
                o = S.ev_relax(case, 'HFfZ', n)
                out['%s_%d' % (case, n)] = o['txt']
            except Exception as ex:  # noqa: BLE001
                out['%s_%d' % (case, n)] = 'ERR %s' % ex
        # time series of the dt/1 level, interior mean, vs ODE
        d = S.run_dir(case, 'HFfZ', 1)
        r, c = S.rcen(case, 1)
        T0, E0 = float(c['prob']['gas_temp']), float(c['prob']['e_out'])
        rho = float(c['prob']['gas_rho'])
        cl, kap = float(c['rad']['c_light']), float(c['rad']['kappa_p'])
        fm, fh = spl.tabs(d), spl.tabs(d, 'hydro_w')
        ts, Es, Ts = [], [], []
        for a, b in zip(fm, fh):
            tm, th = spl.tab(a), spl.tab(b)
            m = r < 1.8
            ts.append(tm['time'])
            Es.append(np.mean(tm['m1_e'][m]))
            Ts.append(np.mean((th['eint']*2.0/3.0/th['dens'])[m]))
        tt = np.linspace(0, max(ts), 400)[1:]
        ref = S.relax_ode(T0, E0, rho, 2.0/3.0, 1.0, cl, kap, list(tt))
        out[case + '_ts'] = dict(t=np.array(ts), E=np.array(Es), T=np.array(Ts), tt=tt,
                                 Eref=np.array([ref[x][0] for x in tt]), Tref=np.array([ref[x][1] for x in tt]))
    return out


def sp_marsh():
    import spl
    import s3_eval as S
    out = {}
    for case, arm, levs in (('marsh_r1_u', 'HFfZ', (64, 128, 256)), ('marsh_r1_s', 'HFfZ', (64, 128, 256)),
                            ('marsh_r1000_u', 'HFfZ', (128,)), ('marsh_cart', 'HCFfZ', (128,))):
        for n in levs:
            d = S.run_dir(case, arm, n)
            t = spl.tab(spl.tabs(d)[-1])
            h = spl.tab(spl.tabs(d, 'hydro_w')[-1])
            r, _ = S.rcen(case, n)
            out['%s_%d' % (case, n)] = dict(r=r, E=t['m1_e'], T=h['eint']*2.0/3.0/h['dens'], dens=h['dens'],
                                            t=t['time'], st=spl.logstats(d))
    return out


def sp_wave():
    import spl
    import s3_eval as S
    out = {}
    for case, arm, levs in (('rsw', 'HFfZ', (64, 128, 256)), ('rswv', 'HFfZ', (64, 128, 256)),
                            ('rsw_r1000', 'HFfZ', (128,)), ('rsw_cart', 'HCFfZ', (128,)),
                            ('rswt_r1000', 'HFfZ', (128,)), ('rswt_cart', 'HCFfZ', (128,)),
                            ('pulse_r1', 'HFfZ', (128,)), ('pulse_r1000', 'HFfZ', (128,)), ('pulse_cart', 'HCFfZ', (128,))):
        for n in levs:
            d = S.run_dir(case, arm, n)
            t = spl.tab(spl.tabs(d)[-1])
            o = dict(E=t['m1_e'], F=t['m1_f1'], t=t['time'], st=spl.logstats(d))
            try:
                h = spl.tab(spl.tabs(d, 'hydro_w')[-1])
                o.update(dens=h['dens'], v=h['velx'])
            except Exception:
                pass
            r, _ = S.rcen(case, n)
            o['r'] = r
            o['E0'] = spl.tab(spl.tabs(d)[0])['m1_e']
            out['%s_%d' % (case, n)] = o
    return out


def sp_patm():
    import patm_eval as P
    out = {}
    for case in ('patm_L16_s', 'patmT_L16_s'):
        for n in (64, 128, 256):
            o = P.ev(case, 'HFfZ', n)
            r, rel = o['rel']
            out['%s_%d' % (case, n)] = dict(r=r, rel=rel, tau=P.tau(r), e_ph=o['e_ph'], Fr2=o['Fr2'],
                                            bands={k: o[k] for k in ('tau<0.1', '0.1-1', '1-10', '10-20')})
    return out


ITEMS = ['cb', 'cbg', 'shd', 'hcyl', 'nlte', 'dif', 'eq', 'rw', 'shock', 'cloud', 'bub', 'marsh', 'pipe', 'g5',
         'sp_core', 'sp_nlte', 'sp_shd', 'sp_pinc', 'sp_dif', 'sp_relax', 'sp_marsh', 'sp_wave', 'sp_patm']


def one(name):
    try:
        res = globals()[name]()
        with open(os.path.join(CACHE, name + '.pkl'), 'wb') as fh:
            pickle.dump(res, fh)
        if os.path.exists(os.path.join(CACHE, name + '.err')):
            os.remove(os.path.join(CACHE, name + '.err'))
        return name, 'ok'
    except Exception:
        open(os.path.join(CACHE, name + '.err'), 'w').write(traceback.format_exc())
        return name, 'ERR'


if __name__ == '__main__':
    os.makedirs(CACHE, exist_ok=True)
    names = sys.argv[1:] or ITEMS
    from multiprocessing import Pool
    with Pool(min(len(names), 24)) as p:
        for n, s in p.imap_unordered(one, names):
            print(n, s, flush=True)
