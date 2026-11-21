"""BSGTAB (11-21, USER: BSG with a MESA Rosseland + TOPS Planck table): X 0.70 Z 0.008 GS98.
Rosseland = the tools/ic_star ext2 X0.7 Z0.008 table (Raven cache opac_5b8e2fdfefea, spec bsg.yaml: TOPS hi+lo,
Ferguson log-Z) with MESA r15140 kap defaults spliced in over MESA's whole range -- EXACTLY the recipe of
/raven/ptmp/jinma/acc3d_1017/mk/build_mk.py (rosseland_mesakap_gs98_x0.70_z0.02, join 8.0-8.05), via
tools/ic_star/icstd/mesakap.py (wt-mesascaf 5d93d11c; copy mesakap_5d93d11c.py here).  This script, run on Orion
with kx/kap_data = gs98_z* + lowT_fa05_gs98_z* extracted from MESA kap_data.tar.xz (md5 76857b1d), reproduces the
z0.02 production table bitwise (data rows) from its ext2 base.
Planck = planck_ext2_x0.7_z0.008_5b8e2fdfefea.txt of the same cache (data rows unchanged, header added)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mesakap_5d93d11c as MK
S = sys.argv[1]       # dir with the base tables + kx/kap_data
D = os.path.dirname(os.path.abspath(__file__))
base = S + '/rosseland_ext2_x0.7_z0.008_5b8e2fdfefea.txt'
pbase = S + '/planck_ext2_x0.7_z0.008_5b8e2fdfefea.txt'
out = D + '/rosseland_mesakap_gs98_x0.70_z0.008.txt'
pout = D + '/planck_ext2_gs98_x0.70_z0.008.txt'
RB = '/raven/ptmp/jinma/ic_star_cache/tables/opac_5b8e2fdfefea'
kh = '76857b1d386f411100d07746028ae4be'
mk = dict(X=0.70, Z=0.008, join=[8.0, 8.05])
hdr = ['AthenaK stellar Rosseland opacity table, log10 kappa_R [cm^2/g] on (log10 T[K], log10 rho[g/cm3]); INCLUDES electron scattering',
       'bsgtab-1121 (11-21) MESA-CONSISTENT BSG: %s; MESA kap_data md5 %s; X 0.70 Z 0.008 (lin X/Z interp: MESA Z 0.004 / 0.01 tables, w 2/3 on 0.01; Type2 off)' % (MK.VERSION, kh),
       'JOIN: MESA kappa (gs98 OPAL above log T 3.88, lowT_fa05_gs98 below 3.80, cos blend) at EVERY grid T (join log T 8.0-8.05 > table top 8.0)',
       'wherever MESA log R = log rho - 3 log T + 18 is in [-8, 1]; ext2 base kept for log R > 1 and for log R < -8, blended linearly over',
       'log R -8..-7.5 (0.5 dex); MESA clamps log T < 2.7 to its edge. Base: %s/rosseland_ext2_x0.7_z0.008_5b8e2fdfefea.txt (md5 %s) = tools/ic_star ext2 X0.7 Z0.008 GS98 (TOPS hi+lo, Ferguson05 log-Z, low-rho ext.).' % (RB, MK.md5(base)),
       'Same recipe as rosseland_mesakap_gs98_x0.70_z0.02 (acc3d_1017/mk/build_mk.py). vs L. Ma profile kappa: ratio 0.986-1.002, median 0.995. Planck: planck_ext2_gs98_x0.70_z0.008. Built by %s/build_mk_bsg.py. Read with he_opac_table.' % D]
MK.splice(base, out, S + '/kx/kap_data', mk, print, hdr)
L = open(pbase).read().splitlines()
ph = ['# bsgtab-1121: = %s/planck_ext2_x0.7_z0.008_5b8e2fdfefea.txt (md5 %s; tools/ic_star bsg.yaml), data rows unchanged; TOPS Planck (hi+lo) above log T 4.2,' % (RB, MK.md5(pbase)),
      '# Ferguson05 GS98 Planck (log-Z interp Z 0.004/0.01) below 4.0, blend 4.0-4.2. Pair with rosseland_mesakap_gs98_x0.70_z0.008.txt (same grid).']
with open(pout, 'w') as f:
    f.write(L[0] + '\n' + '\n'.join(ph) + '\n' + '\n'.join(L[1:]) + '\n')
for fn in (out, pout):
    print(fn, MK.md5(fn))
