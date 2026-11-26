"""addkeys.py DECK sp|cs : name the keys a restart with -i would otherwise take from the file"""
import sys
f, g = sys.argv[1:3]
K = {'mesh': [('nghost', '3'), ('layout_invariant_mode', 'tile')],
     'time': [('dt_retry_grow', '1.1'), ('dt_retry_ngood', '20')],
     'hydro': [('wellbalance_dynamic', 'true'), ('wb_x1', 'true')]
              + ([('sp_x2_periodic_image', 'true')] if g == 'sp' else [])
              + [('etotgrav', 'true')],
     'rad_vet': [('force_reference', 'wb_arad'), ('clean_init', 'steady'),
                 ('implicit_bc_x1min', 'flux'), ('clean_maxit', '200'), ('opacity', 'table'),
                 ('vet_scatter', 'true'), ('closure', 'sc')]
                + ([('vet_gd_gl_nphi', '18'), ('vet_gd_wall_interp', 'true'),
                    ('vet_gd_wall_sym', 'true'), ('vet_gd_theta_seam_sym', 'true')]
                   if g == 'sp' else [])
                + [('time_scheme', 'hesdirk2')],
     'problem': [('star_sponge_form', 'lawson')]
                + ([('star_ic_rad_atm_maxit', '1000'), ('star_ic_rad_atm_relax', '0.5')]
                   if g == 'cs' else [])}
L = open(f).read().split('\n')
blk, out, named = None, [], set()
for ln in L:
    s = ln.split('#')[0].strip()
    if s.startswith('<'):
        blk = s[1:-1]
    elif '=' in s and blk:
        named.add((blk, s.split('=')[0].strip()))
def flush(b):
    add = [(k, v) for k, v in K.get(b, []) if (b, k) not in named]
    if not add: return []
    w = max(len(k) for k, _ in add)
    return ['# the values the start-up decides (a restart with -i takes them from the deck)'] + \
           ['%s = %s' % (k.ljust(w), v) for k, v in add]
for i, ln in enumerate(L):
    s = ln.strip()
    if s.startswith('<') and blk is not None:
        # insert before trailing blank lines of the previous block
        j = len(out)
        while j > 0 and out[j-1].strip() == '': j -= 1
        out[j:j] = flush(blk)
    if s.startswith('<'):
        blk = s.split('#')[0].strip()[1:-1]
    out.append(ln)
j = len(out)
while j > 0 and out[j-1].strip() == '': j -= 1
out[j:j] = flush(blk)
open(f, 'w').write('\n'.join(out))
