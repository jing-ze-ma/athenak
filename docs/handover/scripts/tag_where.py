"""Where does the floor tag sit?  usage: python3 tag_where.py BIN [va_cgs=2e6] [dfloor=1e-18]
       [gamma=1.4727545958] [BIN_OFF=...]
Reads an AthenaK mhd_w_bcc bin (floor tracer = LAST scalar, as s/rho) and prints the tag
mass sum(rho*r*vol-proxy) by log10(p [bar]) band and by rho/rho_fl band, rho_fl =
max(dfloor, |B|^2/v_A,max^2) (code = cgs units, |B|^2 from bcc; on the cubed sphere bcc is
not orthonormal, so rho_fl is approximate near panel corners).  Pressure is estimated as
(gamma-1) e_int (for banding only).  Cell volumes: x1 face radii^3 x solid angle are not in
the bin, so the 'mass' columns weight cells by rho r^2 dr only within a MeshBlock row
(relative numbers; good for fractions).  With BIN_OFF (same grid, tracer off or tag-only
arm) it also prints, per band, the density difference (BIN - BIN_OFF)/BIN_OFF next to
the tag fraction of BIN_OFF... (if BIN_OFF has the tag) to test whether a shell change
equals the floor mass removed.
"""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/floortrc_1101/suite/vis/python')
import bin_convert as bc  # noqa: E402

args = [a for a in sys.argv[1:] if '=' not in a]
kv = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
va = float(kv.get('va_cgs', 2e6))
dfl = float(kv.get('dfloor', 1e-18))
gam = float(kv.get('gamma', 1.4727545958))


def load(f):
    d = bc.read_binary(f)
    mb = d['mb_data']
    names = d['var_names']
    sn = [n for n in names if n.startswith('s_')]
    out = {}
    for key in ['dens', 'eint', 'bcc1', 'bcc2', 'bcc3'] + sn[-1:]:
        out[key] = np.stack([np.asarray(mb[key][m]) for m in range(len(mb[key]))])
    x1f = np.stack([np.asarray(d['mb_x1f'][m]) for m in range(len(mb['dens']))])
    out['tag'] = out[sn[-1]] if sn else None
    out['x1f'] = x1f
    return out


a = load(args[0])
rho = a['dens']
b2 = a['bcc1']**2 + a['bcc2']**2 + a['bcc3']**2
rfl = np.maximum(dfl, b2/va**2)
x = rho/rfl
pbar = (gam - 1)*a['eint']/1e6
r = a['tag'] if a['tag'] is not None else np.zeros_like(rho)
# radial weight r^2 dr per cell (broadcast over (m, k, j, i))
xf = a['x1f']
w1 = ((xf[:, 1:]**3 - xf[:, :-1]**3)/3.0)[:, None, None, :]
w = np.broadcast_to(w1, rho.shape)
m_all = rho*w
m_tag = rho*r*w
print(f'file {args[0]}: total tag mass fraction {m_tag.sum()/m_all.sum():.3e}')
print('| log p [bar] | tag mass share | tag/mass in band | share at rho<2rho_fl | share 2..4 | '
      'share >4 |')
print('|---|---|---|---|---|---|')
lp = np.log10(np.maximum(pbar, 1e-30))
tot = m_tag.sum()
for lo in range(-12, 1):
    s = (lp >= lo) & (lp < lo + 1)
    if not s.any() or m_tag[s].sum() <= 0:
        continue
    t = m_tag[s]
    xs = x[s]
    print(f'| {lo}..{lo+1} | {t.sum()/tot:.3f} | {t.sum()/m_all[s].sum():.2e} | '
          f'{t[xs < 2].sum()/t.sum():.3f} | {t[(xs >= 2) & (xs < 4)].sum()/t.sum():.3f} | '
          f'{t[xs >= 4].sum()/t.sum():.3f} |')
if len(args) > 1:
    b = load(args[1])
    print(f'\nvs {args[1]}: per band, mass ratio (A/B - 1) and tag/mass of B')
    lpb = np.log10(np.maximum((gam - 1)*b['eint']/1e6, 1e-30))
    mb_all = b['dens']*w
    rb = b['tag'] if b['tag'] is not None else np.zeros_like(rho)
    for lo in range(-9, -2):
        s = (lp >= lo) & (lp < lo + 1)
        sb = (lpb >= lo) & (lpb < lo + 1)
        if not s.any():
            continue
        print(f'| {lo}..{lo+1} | {m_all[s].sum()/mb_all[sb].sum() - 1:+.3e} | '
              f'{(b["dens"]*rb*w)[sb].sum()/mb_all[sb].sum():.3e} |')
