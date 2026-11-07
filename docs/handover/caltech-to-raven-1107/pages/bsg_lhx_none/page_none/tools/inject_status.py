"""inject_status.py PAGEDIR: inline status.json + movie/n.json into PAGEDIR/index.html (as on the
Raven-published page: window.__ST / window.__NJ, so the page renders without the fetches), and set the
movie sliders' static max/value/frame to the newest frame."""
import json, re, sys
pd = sys.argv[1].rstrip('/')
st = json.load(open(pd + '/status.json'))
nj = json.load(open(pd + '/movie/n.json'))
h = open(pd + '/index.html').read()
h = re.sub(r'<script>window\.__ST=.*?</script>\n', '', h, flags=re.S)
js = ('<script>window.__ST=' + json.dumps(st) + ';window.__NJ=' + json.dumps(nj) +
      ';(function(){var d=window.__ST;Object.keys(d).forEach(function(id){var el=document.getElementById(id);'
      'if(el)el.innerHTML=d[id];});})();</script>\n')
key = '<script>\n// Fills run status'
assert h.count(key) == 1
h = h.replace(key, js + key)
n = nj['n'] - 1
for k in ('f3', 'f4', 'fint'):
    h, c1 = re.subn(r'(data-k="%s" min="0" max=")\d+(" value=")\d+(")' % k, r'\g<1>%d\g<2>%d\g<3>' % (n, n), h)
    h, c2 = re.subn(r'(id="mv-time-%s">)t = [^<]*(<)' % k, r'\g<1>t = %.2f d\g<2>' % nj['t'][-1], h)
    h, c3 = re.subn(r'(id="mv-%s" src="movie/%s/%s_)\d{5}(\.png")' % (k, k, k), r'\g<1>%05d\g<2>' % n, h)
    assert c1 == c2 == c3 == 1, (k, c1, c2, c3)
# Caltech 10-14: no 30-57.4 d dumps yet -> the static figure (img/ours_cmp57.png) is replaced by the note
import os
if 'cmp-fig' in st and not os.path.exists(pd + '/img/ours_cmp57.png'):
    h, c4 = re.subn(r'(<div id="cmp-fig">).*?</figure></div>', lambda m: m.group(1) + st['cmp-fig'] + '</div>', h,
                    count=1, flags=re.S)
    assert c4 == 1
open(pd + '/index.html', 'w').write(h)
print('injected', len(st), 'status keys, n =', nj['n'])
