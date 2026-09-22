import sys
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert  # noqa: E402
import numpy as np  # noqa: E402
a = bin_convert.read_binary(sys.argv[1])
b = bin_convert.read_binary(sys.argv[2])
out = []
for v in a['var_names']:
    da = np.array(a['mb_data'][v])
    db = np.array(b['mb_data'][v])
    d = np.abs(da - db)
    r = d / np.maximum(np.abs(da), 1e-300)
    out.append('%s maxrel %.2e' % (v, r.max()))
print('  '.join(out))
