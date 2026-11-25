"""lzf.py: 25t: bsg L_z and mass flux by radial layer from the hydro_w dumps (float32) of FW300 (4a7885d73) and NX300 (79c1cabef): first / 10 % / 90 % / last active radial cell layer.
L_z flux = sum rho v_r v_phi r sin(theta) r^2 dOmega, mass flux = sum rho v_r r^2 dOmega over the wedge (solid angle of the wedge's cells, no normalisation to 4 pi; r, theta at cell centres); rms v_phi and mean rho solid-angle weighted"""
import glob, json, re, sys
import numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
B = '/resnick/groups/carnegie_poc/jingze/b1_1124/fwp/bsg/'
J = '/resnick/groups/carnegie_poc/jingze/'
def asm(d, v):
    lg = np.asarray(d['mb_logical']); b = np.asarray(d['mb_data'][v][0]); n3, n2, n1 = b.shape
    P = np.empty((n3*(lg[:, 2].max() + 1), n2*(lg[:, 1].max() + 1), n1), np.float64)
    for m in range(len(lg)): a, e = lg[m, 1], lg[m, 2]; P[e*n3:(e + 1)*n3, a*n2:(a + 1)*n2] = d['mb_data'][v][m]
    return P
pd = J + 'prod_c1_1124/bsg'; base = 'bsg3d'
_rows = [l.split() for l in open(pd + '/' + base + '.x1grid.txt') if not l.startswith('#')]; _n = int(re.search(r'is = 3 ie = (\d+)', open(pd + '/' + base + '.x1grid.txt').readline()).group(1)) - 2
gr = np.array([x[:4] for x in _rows if len(x) == 4], float); rc = gr[3:3 + _n, 2]
print('nx1 %d, r first %.4e last %.4e; layers (0-based index, r): ' % (_n, rc[0], rc[-1]), end='')
LAY = [0, int(round(0.1*_n)), int(round(0.9*_n)), _n - 1]; print(', '.join('%d (%.4e)' % (i, rc[i]) for i in LAY))
print('| arm | dump | t | layer | r | L_z flux (code/cgs: g cm^2 s^-2) | mass flux (g/s) | mean rho | rms v_phi | mean v_phi | L_z flux / (mass flux x r^2 <sin>^... no) |'.replace(' | L_z flux / (mass flux x r^2 <sin>^... no)', ''))
print('|---|---|---|---|---|---|---|---|---|---|')
R = {}
for arm in ('FW300', 'NX300'):
    for f in sorted(glob.glob(B + arm + '/bin/*.hydro_w.*.bin')):
        h = bc.read_binary(f); h1 = open(f, 'rb').read(300).decode('latin-1'); cy = int(re.search(r'cycle=\s*(\d+)', h1)[1])
        rho, vr, vt, vp = asm(h, 'dens'), asm(h, 'velx'), asm(h, 'vely'), asm(h, 'velz')
        N3, N2, n1 = rho.shape; th = np.linspace(h['x2min'], h['x2max'], N2 + 1); W2 = np.cos(th[:-1]) - np.cos(th[1:]); tc = 0.5*(th[1:] + th[:-1]); dph = (h['x3max'] - h['x3min'])/N3; dO = (W2[None, :, None]*dph)*np.ones((N3, 1, 1))
        for nm, i in zip(('first', '10 %', '90 %', 'last'), LAY):
            q = rho[:, :, i]; w = dO[:, :, 0]; v = vr[:, :, i]; p = vp[:, :, i]; s = np.sin(tc)[None, :]
            lz = (q*v*p*rc[i]*s*rc[i]**2*w).sum(); md = (q*v*rc[i]**2*w).sum()
            R[(arm, cy, nm)] = (lz, md); print('| %s | c%d | %.6e | %s | %.5e | %.5e | %.5e | %.5e | %.5e | %.5e |' % (arm, cy, h['time'], nm, rc[i], lz, md, (q*w).sum()/w.sum(), np.sqrt((p*p*w).sum()/w.sum()), (p*w).sum()/w.sum()))
        del h, rho, vr, vt, vp
ks = sorted(set(k[1] for k in R if k[0] == 'FW300'))
print('\nsolid angle of the wedge cells: %.6f sr (x2 %.5f-%.5f, x3 %.5f-%.5f)' % (float((W2.sum()*dph)*N3), float(th[0]), float(th[-1]), 0.0, float(dph*N3)))
print('\n| dump | layer | L_z flux FW300 | NX300 | FW/NX-1 | mass flux FW300 | NX300 | diff |\n|---|---|---|---|---|---|---|---|')
for c in ks:
    for nm in ('first', '10 %', '90 %', 'last'):
        a, b = R[('FW300', c, nm)], R.get(('NX300', c, nm), (float('nan'),)*2); print('| c%d | %s | %.5e | %.5e | %.3e | %.5e | %.5e | %.3e |' % (c, nm, a[0], b[0], a[0]/b[0] - 1, a[1], b[1], a[1] - b[1]))
