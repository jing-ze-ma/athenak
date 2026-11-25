# hero_rsg: crop viper's mov3d_rsg frames to the star (no title, no bars), pad to 16:9 black, 960x540; loop seam chosen by image difference, 6-frame cross-fade
import numpy as np, glob, os
from PIL import Image
F = sorted(glob.glob('fr/f_*.png'))
def bb(i):
    a = np.asarray(Image.open(F[i]).convert('L'), float); a[:150] = 0; a[1305:] = 0; a[:, 1600:] = 0; y, x = np.where(a > 30); return x.min(), y.min(), x.max(), y.max()
B = np.array([bb(i) for i in range(0, len(F), 3)]); x0, y0 = B[:, 0].min()-12, B[:, 1].min()-12; x1, y1 = B[:, 2].max()+12, B[:, 3].max()+12; print('star box over the movie', x0, y0, x1, y1)
BOX = (int(x0), int(max(y0, 150)), int(x1), int(min(y1, 1305)))
def sm(i): return np.asarray(Image.open(F[i]).convert('L').crop(BOX).resize((255, 275)), float)
I = [sm(i) for i in range(len(F))]
cons = np.median([np.sqrt(((I[i]-I[i+1])**2).mean()) for i in range(len(I)-1)])
best = min((np.sqrt(((I[s]-I[e])**2).mean()), s, e) for s in range(0, len(F)-140) for e in range(s+140, min(len(F), s+186)))
r, s, e = best; print('seam', s, e, 'rms %.1f vs consecutive %.1f' % (r, cons))
K = 6
def full(i):
    c = Image.open(F[i]).convert('RGB').crop(BOX); h = 540; w = int(round(c.width*h/c.height)); c = c.resize((w, h), Image.LANCZOS); o = Image.new('RGB', (960, 540)); o.paste(c, ((960-w)//2, 0)); return o
os.makedirs('seq', exist_ok=True)
for f in glob.glob('seq/*.png'): os.remove(f)
seq = [full(i) for i in range(s+K, e-K+1)] + [Image.blend(full(e-K+k), full(s+k), k/(K+1.0)) for k in range(1, K+1)]
for i, im in enumerate(seq): im.save('seq/h_%05d.png' % i)
full(e-K).save('../out/hero_rsg_poster.png'); print('frames', len(seq), 'source frames', s+K, '-', e, 'of', len(F))
