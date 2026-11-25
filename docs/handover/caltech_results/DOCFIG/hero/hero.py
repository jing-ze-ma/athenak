# hero_bsg: dumps 51-101 plain, then 102-105 cross-faded into 47-50 (seam 105 -> 47: rms difference 19.8 vs 18.6 between consecutive dumps) (the loop seam), 960x540; time stamp taken whole from the dominant frame
from PIL import Image
import os
os.makedirs('heroseq', exist_ok=True)
for f in os.listdir('heroseq'): os.remove('heroseq/' + f)
L = lambda n: Image.open('hero/bsg_%05d.png' % n).convert('RGB')
seq = [L(n) for n in range(51, 102)]
for k, n in enumerate(range(102, 106)):
    a = (k + 1) / 5.0; A, B = L(n), L(n - 55); im = Image.blend(A, B, a); box = (0, 0, 420, 110); im.paste((B if a >= 0.5 else A).crop(box), box[:2]); seq.append(im)
for i, im in enumerate(seq): im.resize((960, 540), Image.LANCZOS).save('heroseq/h_%05d.png' % i)
L(101).resize((960, 540), Image.LANCZOS).save('out/hero_bsg_poster.png'); print(len(seq))
