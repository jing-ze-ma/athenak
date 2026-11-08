"""tanfix-1109 decks (inp2/): arms per (g, lev, case):
 off, on                : twin off (as tanacc-1107)
 off_tw, on_tw          : + the vet_gd twin (fused, det, noesrc)
 on_twj                 : + twin weight j
 on_hd                  : ~400-direction reference (twin off)
 C_off, C_on, C_on_tw   : implicit_flux = central (pure tensor flux; implicit_flux_x1min KEPT)"""
import os, re
import mkdecks as M
D = '/viper/ptmp2/jinma/openitems_1107/tanacc/inp2'
TW = ("vet_gd_twin = false\nvet_gd_twin_fuse = false\n", "vet_gd_twin = true\nvet_gd_twin_fuse = true\nvet_gd_twin_det = true\n")
def twin(s, wj=False):
    assert TW[0] in s
    s = s.replace(TW[0], TW[1] + ("vet_gd_twin_weight = j\n" if wj else ""))
    return s.replace("vet_source_noesrc = false", "vet_source_noesrc = true")
def central(s):
    out = []
    for l in s.split('\n'):
        k = l.split('=')[0].strip()
        if k in ('implicit_flux', 'implicit_flux_faces', 'implicit_flux_beam', 'implicit_hr_face_gd') \
           or k.startswith('implicit_blend') or k.startswith('implicit_thin_'):
            continue
        out.append(l)
    s = '\n'.join(out)
    return s.replace('<rad_m1>\n', '<rad_m1>\nimplicit_flux = central\n', 1)
os.makedirs(D, exist_ok=True)
names = []
for g in ('sp', 'cs'):
    for lev in ('lo', 'hi'):
        for case in ('star', 'shd', 'semi', 'mid', 'hot5', 'hot'):
            on, off = M.deck(g, lev, case, 'on'), M.deck(g, lev, case, 'off')
            arms = {'off': off, 'on': on, 'off_tw': twin(off), 'on_tw': twin(on), 'on_twj': twin(on, True),
                    'on_hd': M.deck(g, lev, case, 'on', 'hd'), 'C_off': central(off), 'C_on': central(on),
                    'C_on_tw': central(twin(on))}
            for a, s in arms.items():
                nm = f'{g}_{lev}_{case}_{a}'
                open(f'{D}/{nm}.athinput', 'w').write(s)
                names.append(nm)
open(D + '/LIST', 'w').write('\n'.join(names) + '\n')
print(len(names))
