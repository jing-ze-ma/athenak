#!/bin/bash
# analysis of the GPU jobs (f12, fw2, f4, f8 timing; prof trace).  usage: bash finish.sh
W=/viper/ptmp2/jinma/krylov_0923; R=$W/runs
python3 $W/scripts/tsum.py $R/f12 $R/fw2 $R/f4 $R/f8
python3 $W/scripts/profsum.py $R/prof
echo "== GPU exactness (hst files, repeat a): off vs head, hm vs off (expect identical)"
for p in "f12 head_g1 g1_off" "f12 head_g2 g2_off" "f12 g2_off g2_hm" "fw2 w2_off w2_hm" \
         "f4 g4_off g4_hm" "f4 w4_off w4_hm" "f8 g8_off g8_hm" "f8 w8_off w8_hm"; do
  set -- $p; r=IDENTICAL
  for f in hydro user; do cmp -s $R/$1/${2}_a/m1slab.$f.hst $R/$1/${3}_a/m1slab.$f.hst || r=DIFF; done
  echo "$1 $2 vs $3: $r"
done
echo "== pipe vs off, max |d|/colmax over rows (user.hst, hydro cols 7-10), repeat a"
for p in "f12 g1_off g1_pipe" "f12 g2_off g2_pipe" "f12 g2_off g2_ph" "fw2 w2_off w2_ph" \
         "f4 g4_off g4_ph" "f4 w4_off w4_ph" "f8 g8_off g8_ph" "f8 w8_off w8_ph" \
         "f12 g1_off g2_off"; do
  set -- $p
  python3 - $R/$1/${2}_a $R/$1/${3}_a <<'PY'
import sys, numpy as np
a, b = sys.argv[1], sys.argv[2]
try:
    out = []
    for f, cols in (('user', slice(1, None)), ('hydro', slice(6, 10))):
        x = np.loadtxt(a + '/m1slab.%s.hst' % f)[:, cols]
        y = np.loadtxt(b + '/m1slab.%s.hst' % f)[:, cols]
        n = min(len(x), len(y))
        cm = np.maximum(np.abs(x[:n]).max(0), 1e-300)
        out.append('%s %.1e' % (f, (np.abs(x[:n] - y[:n]) / cm).max()))
    print(a.split('/')[-2], a.split('/')[-1], b.split('/')[-1], ' '.join(out))
except Exception as e:
    print(a, b, 'missing', e)
PY
done
