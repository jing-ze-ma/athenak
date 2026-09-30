import sys, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts'); sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import dhjcs
eos = dhjcs.EOS('/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt')
raw = dhjcs.bin_convert.read_binary('/viper/ptmp2/jinma/w121prod_0929/w1x/bin/dhj.hydro_w.00150.bin')
d = np.asarray(raw['mb_data']['dens'], float); ei = np.asarray(raw['mb_data']['eint'], float)
T, p = eos.invert(d, ei)
for i in (70, 72, 73, 74, 75):
    t = T[..., i]; print(i, 'T min/med/max %.0f %.0f %.0f  frac<=205K %.3f  rho med %.2e' % (np.nanmin(t), np.nanmedian(t), np.nanmax(t), np.mean(t <= 205), np.median(d[..., i])))
