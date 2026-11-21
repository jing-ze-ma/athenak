"""seam_only.py: tools/ic_star seam fixed point (taper off, startup hopf; = make_star_ic.taper_and_seam)
on the pass-A column, with the bsgsp deck and the 3f15a44c CPU he_star binary.  r136a1_1104 (11-04)."""
import json, os, re, sys
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/wt_prod1102/tools/ic_star')
from icstd import seam as SM
A = '/resnick/groups/carnegie_poc/jingze/r136a1_1104/ic/AB'; B = '/resnick/groups/carnegie_poc/jingze/r136a1_1104/ic/BB'
deck = '/resnick/groups/carnegie_poc/jingze/r136a1_1104/decks/r136a1_128_B_A.athinput'
blk = open(A + '/deck_r136a1B_noseam.athinput').read()
gv = lambda k: re.search(r'^%s\s*=\s*(\S+)' % k, blk, re.M).group(1)
L = 2.60446120e+40; x1max = 6.957e12
ov = ['mesh/nx2=4', 'mesh/nx3=4', 'meshblock/nx2=4', 'meshblock/nx3=4', 'meshblock/nx1=831',
      'problem/he_opac_table=' + gv('he_opac_table'), 'problem/he_planck_table=' + gv('he_planck_table'),
      'problem/he_gm=' + gv('he_gm'), 'problem/he_ic_rad_taua=1', 'problem/he_ic_rad_taub=3',
      'problem/he_nfine=262144', 'problem/he_gm_column=true',
      'rad_m1/implicit_flux_x1min=1.131943656797e15']
cfg = dict(enable=True, deck=deck, binary='/resnick/groups/carnegie_poc/jingze/builds_rstfixx/bin/athena_cpu_he_star_m1_bsgsp_3f15a44c',
           tol_cells=0.1, maxit=8, omp=4, timeout=7200, workdir=B, overrides=ov,
           launcher=['bash', '-lc', 'module purge >/dev/null 2>&1; module load gcc/13.2.0-gcc-13.2.0-w55nxkl '
                     'openmpi/5.0.1-gcc-13.2.0-xpoh5uw >/dev/null 2>&1; exec "$0" "$@"'])
fin, hist = SM.fixed_point(cfg, None, None, L, None, x1max,
                           ['# tools/ic_star r136a1 column, seam fixed point (r136a1_1104/tools/seam_only.py)'],
                           A + '/ic_r136a1B_noseam.txt', lambda x: print(x, flush=True), None, 30.0)
print('SEAM_DONE', fin, json.dumps(hist))
