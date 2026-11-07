"""Dark version of Ma+2026 Fig. 4 (raster crop img/paper_fig4.png -> img/paper_fig4_dark.png).

The three image panels and the three colour bars (body + extend triangles) keep their pixels
exactly; everything outside them (white background, black text, axes frames, ticks) is
inverted (255 - v), so the page shows it white on black.  Panel / colour-bar interiors are
found as the large regions enclosed by the figure's dark neutral outlines (axes frames,
colour-bar outlines), i.e. not reachable from the image border through light pixels.

usage: python3 dark_paper_fig4.py [PAGEDIR]   (default /resnick/home/jingze/ATHENAK/pages/bsg_li_caltech)
"""
import sys

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

pd = (sys.argv[1] if len(sys.argv) > 1 else '/resnick/home/jingze/ATHENAK/pages/bsg_li_caltech').rstrip('/')
src, dst = pd + '/img/paper_fig4.png', pd + '/img/paper_fig4_dark.png'
im = Image.open(src).convert('RGB')
a = np.asarray(im).astype(np.int16)
mx, mn = a.max(axis=2), a.min(axis=2)
gray = a.mean(axis=2)
# dark neutral pixels = text, frames, outlines (low saturation, not light)
line = ((mx - mn) < 40) & (gray < 200)
# candidate interiors: everything that is not line; drop the part connected to the image border
lab, nl = ndi.label(~line)
border = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
border = border[border > 0]
cand = (lab > 0) & ~np.isin(lab, border)
# keep only large enclosed regions (panels ~ 1e5 px, colour bars ~ 1e4 px; letter holes < 1e3)
lab2, n2 = ndi.label(cand)
sizes = ndi.sum(cand, lab2, index=np.arange(1, n2 + 1))
big = np.flatnonzero(sizes > 3000) + 1
keep = np.zeros_like(cand)
for k in big:
    keep |= ndi.binary_fill_holes(lab2 == k)   # restore any neutral pixels inside a panel
objs = ndi.find_objects(ndi.label(keep)[0])
print('kept regions (x0, x1, y0, y1, px):')
for s in objs:
    print('  ', s[1].start, s[1].stop, s[0].start, s[0].stop, int(keep[s].sum()))
assert len(objs) == 6, 'expected 3 panels + 3 colour bars, found %d' % len(objs)
out = np.asarray(im).copy()
out[~keep] = 255 - out[~keep]
assert np.array_equal(out[keep], np.asarray(im)[keep])
Image.fromarray(out).save(dst, optimize=True)
print('wrote', dst, out.shape)
