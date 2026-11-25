import numpy as np
from PIL import Image
I={n:np.asarray(Image.open('hero/bsg_%05d.png'%n).convert('L').resize((480,270)),float) for n in range(30,115)}
for a in I.values(): a[:30,:110]=0
best=sorted(((np.sqrt(((I[s]-I[e])**2).mean()),s,e) for s in range(30,61) for e in range(s+54,115)))[:6]; print(best)
print('consecutive-dump rms (typical motion):',np.median([np.sqrt(((I[n]-I[n+1])**2).mean()) for n in range(30,114)]))
