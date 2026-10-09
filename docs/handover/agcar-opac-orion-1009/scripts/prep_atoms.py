"""Parse NIST ASD levels + ionization energies, Kurucz gfall lines and the Verner (1996)
ground-state photoionization fits into one npz per run (atoms.npz).
usage: python prep_atoms.py SRC   SRC = gf08 (Kurucz gfall08oct17, lines between known
levels) or cd1 (Kurucz CD-ROM 1 lowlines+highlines, incl. predicted lines)
-> work/atoms_SRC.npz"""
import csv, glob, os, re, sys
import numpy as np

B = '/orion/ptmp/jinma/agcar_opac_1009/'
ROM = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI']
SYM = ['x', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
       'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
       'Cu', 'Zn']
ZLIST = [1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
         26, 27, 28, 29, 30]
NST = 10            # ionization stages carried: 0..min(Z,NST) (stage index = ion charge)


def clean(s):
    return s.replace('=', '').replace('"', '').strip()


def num(s):
    s = re.sub(r'[\[\]\(\)\?\+x a-zA-Z]', '', clean(s))
    try:
        return float(s)
    except ValueError:
        return None


def read_ie():
    ie = {}
    with open(B + 'dl/nist/ie.csv') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) < 10:
                continue
            Z, q, e = num(row[0]), num(row[2]), num(row[9])
            if Z is None or q is None or e is None:
                continue
            ie[(int(Z), int(q))] = e
    return ie


def read_levels(Z, s, ip):
    f = B + f'dl/nist/lev_{SYM[Z]}_{ROM[s]}.csv'
    g, E, conf = [], [], []
    if os.path.exists(f):
        with open(f) as fh:
            r = csv.reader(fh)
            try:
                next(r)
            except StopIteration:
                r = []
            for row in r:
                if len(row) < 6 or clean(row[3]) == '':
                    continue
                gg, ee = num(row[3]), num(row[5])
                if gg is None or ee is None or ee >= ip:
                    continue
                g.append(gg), E.append(ee), conf.append(clean(row[0]))
    return np.array(g), np.array(E), conf


def read_gfall():
    """returns dict code(Z,q) -> arrays (nu [Hz], gf, Elow [eV], gammaR [s^-1])"""
    out = {}
    cm2eV = 1.239841984e-4
    with open(B + 'dl/kurucz/gfall08oct17.dat') as f:
        for ln in f:
            if len(ln) < 124:
                continue
            try:
                wl = float(ln[0:11])          # nm (vacuum < 200 nm, air above)
                lgf = float(ln[11:18])
                code = float(ln[18:24])
                e1, e2 = abs(float(ln[24:36])), abs(float(ln[52:64]))
                gr = float(ln[80:86])
                hfs = float(ln[109:115]) if ln[109:115].strip() else 0.0
                iso = float(ln[118:124]) if ln[118:124].strip() else 0.0
            except ValueError:
                continue
            Z = int(code)
            q = int(round((code - Z)*100))
            el = min(e1, e2)
            # wavenumber from the level energies (vacuum, exact); fall back on wl
            dE = abs(e2 - e1)
            sig = dE if dE > 0 else 1e7/wl
            out.setdefault((Z, q), []).append((sig*2.99792458e10, 10**(lgf + hfs + iso),
                                               el*cm2eV, 10**gr if gr != 0 else 0.0))
    return {k: np.array(v).T for k, v in out.items()}


def read_cd1():
    """Kurucz CD-ROM 1 packed lists (ASCII): IWL IELION IELO IGFLOG IGR IGS IGW"""
    src = open(B + 'dl/kurucz/readlow.for').read()
    blk = re.findall(r'DATA ELEM[A-H]/(.*?)/', src, re.S)
    codes = []
    for b in blk:
        for ln in b.splitlines():
            codes += [float(x) for x in re.findall(r'\d+\.\d\d', ln[6:])]
    codes = np.array(codes)
    rl = np.log(1.0 + 1.0/2000000.0)
    out = {}
    import gzip
    for fn, op in [('lowlines.asc.gz', gzip.open), ('highlines.asc', open)]:
        import pandas as pd
        a = pd.read_csv(B + 'dl/kurucz/' + fn, sep=r'\s+', header=None, usecols=range(5),
                        dtype=np.int64, engine='c').values
        wl = np.exp(a[:, 0]*rl)                       # vacuum nm
        code = codes[np.abs(a[:, 1])//10 - 1]
        tl = lambda i: 10**((i - 16384)*0.001)
        nu = 2.99792458e17/wl
        gf, el, gr = tl(a[:, 3]), tl(a[:, 2]), tl(a[:, 4])
        el = np.where(a[:, 2] == 16384, 0.0, el)
        Z = np.floor(code + 1e-6).astype(int)
        q = np.round((code - Z)*100).astype(int)
        cm2eV = 1.239841984e-4
        for (z, qq) in set(zip(Z.tolist(), q.tolist())):
            m = (Z == z) & (q == qq)
            v = np.array([nu[m], gf[m], el[m]*cm2eV, gr[m]])
            out[(z, qq)] = np.concatenate([out[(z, qq)], v], 1) if (z, qq) in out else v
        print(fn, len(a), flush=True)
    return out


def read_verner():
    a = np.loadtxt(B + 'work/verner/verner_ground.txt')
    d = {}
    for Z in ZLIST:
        for ne in range(1, Z+1):
            m = (a[:, 0] == Z) & (a[:, 1] == ne)
            d[(Z, Z - ne)] = (a[m, 2], a[m, 3]*1e-18)    # eV, cm^2 (Mb -> cm^2)
    return d


if __name__ == '__main__':
    ie = read_ie()
    SRC = sys.argv[1]
    gfl = read_gfall() if SRC == 'gf08' else read_cd1()
    ver = read_verner()
    sav = {}
    log = []
    for Z in ZLIST:
        for s in range(0, min(Z, NST) + 1):
            key = f'{Z}_{s}'
            if s == Z:                          # bare nucleus
                sav[key + '_g'], sav[key + '_E'] = np.array([1.0]), np.array([0.0])
                sav[key + '_ip'] = np.array([1e9])
                sav[key + '_gc'] = np.array([1], bool)
                continue
            ip = ie[(Z, s)]
            if Z == 1 or (Z == 2 and s == 1):   # hydrogenic: n-levels, filled in the lib
                g, E, conf = np.array([2.0]), np.array([0.0]), ['1s']
            else:
                g, E, conf = read_levels(Z, s, ip)
                if len(g) == 0:
                    log.append(f'{SYM[Z]} {ROM[s]}: no NIST levels -> ground g=1')
                    g, E, conf = np.array([1.0]), np.array([0.0]), ['?']
            o = np.argsort(E)
            g, E = g[o], E[o]
            conf = [conf[i] for i in o]
            gc = np.array([c == conf[0] for c in conf])   # ground-configuration flag
            sav[key + '_g'], sav[key + '_E'], sav[key + '_ip'] = g, E, np.array([ip])
            sav[key + '_gc'] = gc
            if (Z, s) in gfl:
                sav[key + '_lines'] = gfl[(Z, s)]
            if (Z, s) in ver:
                sav[key + '_vE'], sav[key + '_vS'] = ver[(Z, s)]
    nl = sum(v.shape[1] for k, v in sav.items() if k.endswith('_lines'))
    print('lines kept', nl, 'of', sum(v.shape[1] for v in gfl.values()))
    print('\n'.join(log))
    np.savez(B + f'work/atoms_{SRC}.npz', **sav)
