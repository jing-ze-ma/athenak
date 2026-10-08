"""Digitize Dent+2026 Fig. 5 red curve (T_e vs R/R*) from lit2026/dent2026_fig5.png.
Axes calibrated from the frame: major x ticks 1..5 at px 268.5, 496, 723.5, 951.5, 1179; T=4000 at py 143, 1000 at 823."""
import numpy as np
from PIL import Image
a = np.asarray(Image.open('../lit2026/dent2026_fig5.png').convert('RGB')).astype(int)
r, g, b = a[..., 0], a[..., 1], a[..., 2]
red = (r > 150) & (g < 90) & (b < 90)
red[330:420, 800:1040] = False          # label text
ys, xs = np.nonzero(red)
X = 1 + (xs - 268.5) / 227.6
T = 4000 - (ys - 143) / (823 - 143) * 3000
out = []
for c in np.unique(xs):
    s = xs == c
    out.append((X[s][0], T[s].mean(), T[s].min(), T[s].max()))
out = np.array(out)
# bin to 0.01 R*
xb = np.arange(1.13, 5.0001, 0.01)
tab = []
for x0 in xb:
    s = np.abs(out[:, 0] - x0) < 0.006
    if s.any():
        tab.append((x0, out[s, 1].mean(), out[s, 2].min(), out[s, 3].max()))
tab = np.array(tab)
np.savetxt('out/dent26_fig5_T.txt', tab, fmt='%.3f %.0f %.0f %.0f',
           header='x=R/R*  T_mean  T_min  T_max (column spread) [K]; digitized Fig.5')
for x0 in (1.144, 1.15, 1.2, 1.25, 1.3, 1.5, 1.7, 2, 2.5, 3, 3.5, 4, 4.5, 5):
    i = np.argmin(abs(tab[:, 0] - x0))
    print(x0, tab[i])
