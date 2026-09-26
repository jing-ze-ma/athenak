# g13 sponge-default gate: bin/rst payload identity per pair, hst common columns, restart
import glob, os
R = '/viper/ptmp2/jinma/sprhd_0926/gpu/g15'
def payload(f):
    b = open(f, 'rb').read(); i = b.find(b'<par_end>'); return b[i+9:] if i >= 0 else b
def hst(f):
    return [l.split()[:15] for l in open(f) if not l.startswith('#')]
for a, b in (('Ab5', 'An6'), ('Bb5', 'Bn6')):
    fa = sorted(glob.glob(f'{R}/{a}/bin/*.bin')) + sorted(glob.glob(f'{R}/{a}/rst/*.rst'))
    same = [f for f in fa if payload(f) == payload(f.replace(f'/{a}/', f'/{b}/'))]
    print(a, b, 'files identical %d/%d' % (len(same), len(fa)),
          'hst (13 old columns) identical', hst(f'{R}/{a}/hw.user.hst') == [r[:15] for r in hst(f'{R}/{b}/hw.user.hst')])
rr = sorted(glob.glob(f'{R}/Dn6rr/rst/*.rst'))[-1]
print('restart Dn6rr', os.path.basename(rr), payload(rr) == payload(rr.replace('/Dn6rr/', '/Dn6/')))
ra = [l for l in open(f'{R}/Dn6/hw.user.hst') if not l.startswith('#')]
rb = [l for l in open(f'{R}/Dn6rr/hw.user.hst') if not l.startswith('#')]
print('restart hst rows verbatim %d/%d' % (sum(l in set(ra) for l in rb), len(rb)))
