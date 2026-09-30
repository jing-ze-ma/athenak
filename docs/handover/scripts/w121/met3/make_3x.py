"""w121prod_3x.athinput from w121prod_1x.athinput: only the metallicity-dependent keys
(composition, met/rad_met, 3x ck/CE tables, 3x RCE IC) and what the 3x IC/p_ref move
(ap, grav, x1min, x1max, G_SPARC grid) change.  Writes keys_3x_vs_1x.txt."""
import re
import numpy as np
W = '/viper/ptmp2/jinma/w121prod_0929/'
M = W + 'met3/'
S = np.load(M + 'setup_w121_3x.npz')
CK = '/viper/ptmp2/jinma/ck3x_0930/ckdata3'
TAG = 'w121prod_0930 3x'

c = S['sparc_c']
n = int(S['sparc_nx1'])
kk = [('mesh', 'nx1', '%d' % n), ('meshblock', 'nx1', '%d' % n),
      ('mesh', 'x1min', '%.6e' % S['x1min']), ('mesh', 'x1max', '%.6e' % S['x1max'])]
kk += [('mesh', 'f_stretch_r_c%d' % (i + 1), '%.6f' % v) for i, v in enumerate(c)]
kk += [('hydro', 'rad_met', '0.4771'), ('hydro', 'eos_xh', '0.7188'),
       ('hydro', 'eos_yhe', '0.2420'),
       ('problem', 'met', '0.4771'), ('problem', 'grav', '%.4f' % S['grav']),
       ('problem', 'ap', '%.6e' % S['ap']),
       ('problem', 'ck_table', CK + '/ck/Premixed_3x_g8_11_hiT2.txt'),
       ('problem', 'ck_ce_table', 'CE_tables/FastChem_ck_3x_int_hiT2.txt'),
       ('problem', 'ck_data_dir', CK),
       ('problem', 'ic_profile', M + 'ic_w121_3x.txt')]

blk, out, done, tab = None, [], set(), []
for ln in open(W + 'w121prod_1x.athinput').readlines():
    m = re.match(r'^<(\w+)>', ln)
    if m:
        blk = m.group(1)
        out.append(ln)
        continue
    m = re.match(r'^(\s*)(\w+)(\s*=\s*)(\S*)(.*)$', ln)
    if m and not ln.lstrip().startswith('#'):
        hit = [x for x in kk if x[0] == blk and x[1] == m.group(2)]
        if hit:
            b, k, v = hit[0]
            out.append('%s%s%s%s   # %s: was %s in the 1x input\n'
                       % (m.group(1), k, m.group(3), v, TAG, m.group(4)))
            tab.append('| %s/%s | %s | %s |' % (b, k, m.group(4), v))
            done.add((b, k))
            continue
    out.append(ln)
miss = [x for x in kk if (x[0], x[1]) not in done]
assert not miss, miss
hdr = ('# 3x SOLAR twin of w121prod_1x.athinput (user GO 09-30): identical except the keys marked\n'
       '# "%s" ([M/H] = +0.4771: eos_xh/eos_yhe, met/rad_met, interpolated 3x hiT2 ck/CE\n'
       '# tables ck3x_0930/ckdata3, 3x RCE IC met3/ic_w121_3x.txt, ap/grav/x1min/x1max/G_SPARC\n'
       '# grid from planet_setup.py with the 3x NRS1 p_ref; met3/ and\n'
       '# docs/handover/NOTE-2026-09-30-w121-3x.md).  Solver keys = the 1x values.\n#\n' % TAG)
open(W + 'w121prod_3x.athinput', 'w').write(hdr + ''.join(out))
open(M + 'keys_3x_vs_1x.txt', 'w').write('| key | 1x (w121prod_1x) | 3x (w121prod_3x) |\n'
                                         '|---|---|---|\n' + '\n'.join(tab) + '\n')
print('\n'.join(tab))
