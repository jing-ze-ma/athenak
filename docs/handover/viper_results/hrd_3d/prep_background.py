#!/usr/bin/env python3
"""Prepare the MB24 Fig. 2 background rasters used by make_hrd.py.

Input : HR_fig2_axes.png   (300 dpi crop of figures/HR.pdf, axes frame only)
Output: HR_fig2_axes_grey.png  (light version)
        HR_fig2_axes_dark.png  (dark version, for black slides)

Steps
  1. mask out the f_bin legend text block, pixels x=[8,975], y=[535,1570], with white
  2. desaturate:  g = 255 - (255 - luminance) * 0.80
  3. find the MESA evolutionary tracks and fade them:
     - thin dark features are g below a 41x41 grey closing of g (depth > 25)
     - a connected component is a TRACK if its bounding box spans > 250 px; glyphs
       (region names, mass labels) and the dashed HD-limit segments stay below that
     - light version: g <- g + 0.50 * (closing - g) on track pixels
     - dark version: d = min(255, 2 * (255 - g)), then d <- 0.62 * d on track pixels
"""
import numpy as np
from PIL import Image
from scipy import ndimage

LEG_BOX = (535, 1570, 8, 975)     # y0, y1, x0, x1  of the f_bin legend block
DESAT = 0.80
CLOSE = 41
DEPTH = 25
SPAN = 250                        # px; longer than this -> evolutionary track
FADE_LIGHT = 0.50
FADE_DARK = 0.62

im = np.asarray(Image.open('HR_fig2_axes.png').convert('RGB')).astype(float)
y0, y1, x0, x1 = LEG_BOX
im[y0:y1, x0:x1, :] = 255.0
lum = 0.299*im[:, :, 0] + 0.587*im[:, :, 1] + 0.114*im[:, :, 2]
g = 255.0 - (255.0 - lum)*DESAT

closed = ndimage.grey_closing(g, size=(CLOSE, CLOSE))
lab, n = ndimage.label(closed - g > DEPTH)
track = np.zeros_like(g, dtype=bool)
for i, sl in enumerate(ndimage.find_objects(lab)):
    if sl is None:
        continue
    if max(sl[0].stop - sl[0].start, sl[1].stop - sl[1].start) > SPAN:
        track[sl] |= (lab[sl] == i + 1)     # this component only

light = g.copy()
light[track] = g[track] + FADE_LIGHT*(closed[track] - g[track])
dark = np.minimum(255.0, 2.0*(255.0 - g))
dark[track] *= FADE_DARK

for arr, name in ((light, 'HR_fig2_axes_grey.png'), (dark, 'HR_fig2_axes_dark.png')):
    a = np.clip(arr, 0, 255).astype(np.uint8)
    Image.fromarray(np.repeat(a[:, :, None], 3, axis=2)).save(name)

t = track
print('track pixels: %d of %d components' % (t.sum(), n))
print('light: tracks %d -> %d (median), region shading modes %s, '
      'region-name glyphs %d-%d'
      % (np.median(g[t]), np.median(light[t]),
         sorted(set(np.round(np.unique(g)[-8:]).astype(int)))[:4],
         np.percentile(g[(closed-g > 80)], 5), np.percentile(g[(closed-g > 80)], 95)))
print('dark : tracks %d -> %d (median)'
      % (np.median(np.minimum(255, 2*(255-g))[t]), np.median(dark[t])))
