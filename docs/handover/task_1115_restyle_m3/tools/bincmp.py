# bincmp.py a b: exit 0 if two files are identical, or (AthenaK .bin / .rst) identical
# after removing the embedded parameter-input block and the header-offset line
import sys, re
a, b = (open(p, 'rb').read() for p in sys.argv[1:3])
if a == b: sys.exit(0)
def norm(x):
    i = x.find(b'<par_end>')
    if i < 0: return x
    j = x.find(b'\n', i) + 1
    pre = x[:x.find(b'\n<')] if b'\n<' in x[:i] else x[:i]
    pre = re.sub(rb'header offset=\d+', b'', pre)
    return pre + b'|' + x[j:]
sys.exit(0 if norm(a) == norm(b) else 1)
