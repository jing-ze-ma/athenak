"""atoms_mrg.npz: per ion the larger of the two Kurucz lists (CD-ROM 1 incl. predicted
lines, or gfall08oct17); CD-ROM 1 has no H I and few light-element ions (those are in
Kurucz's separate nltelines), so those come from gfall08oct17.  Prints the choice."""
import numpy as np
B = '/orion/ptmp/jinma/agcar_opac_1009/work/'
a, b = np.load(B + 'atoms_cd1.npz'), np.load(B + 'atoms_gf08.npz')
out = {k: b[k] for k in b.files if not k.endswith('_lines')}
src = []
for k in sorted(set(f for f in a.files + b.files if f.endswith('_lines'))):
    na = a[k].shape[1] if k in a.files else 0
    nb = b[k].shape[1] if k in b.files else 0
    if na > nb:
        out[k] = a[k]
        src.append(f'{k[:-6]}:cd1({na})')
    else:
        out[k] = b[k]
        src.append(f'{k[:-6]}:gf08({nb})')
print(' '.join(src))
np.savez(B + 'atoms_mrg.npz', **out)
