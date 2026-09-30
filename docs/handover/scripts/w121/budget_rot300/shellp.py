"""shell-mean pressure [bar] per x1 index from w1x bin 150 (rot 300), read in place."""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/w121prod_0929/ana_rot300')
sys.argv = ['x', 'none']
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import dhjcs   # noqa
eos = dhjcs.EOS('/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt')
raw = dhjcs.bin_convert.read_binary('/viper/ptmp2/jinma/w121prod_0929/w1x/bin/dhj.hydro_w.00150.bin')
d = np.asarray(raw['mb_data']['dens'], float)
ei = np.asarray(raw['mb_data']['eint'], float)
T, p = eos.invert(d, ei)
lp = np.log10(p/1e6)
out = np.array([10**np.nanmean(lp[..., i]) for i in range(lp.shape[-1])])
np.save('shellp.npy', out)
print(out[[0, 20, 40, 60, 70, 75]])
