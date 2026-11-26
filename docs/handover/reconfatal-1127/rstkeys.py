"""rstkeys.py EXE RST DECK NLIM [block/key=value ...]

The restart-key check of vet-reconfatal (CLEAN WARNING restart lines), on a production restart file:
  a_deck   -r RST -i DECK         k=v  time/nlim=NLIM   -> the CLEAN WARNING restart lines (keys the deck
                                                            does not name, stored value != code default)
  b_named  -r RST -i DECK.named   k=v  time/nlim=NLIM   -> DECK + every listed key at its stored value: 0 lines
  c_ronly  -r RST                 k=v  time/nlim=NLIM   -> the reference: b_named must equal it bitwise
NLIM is absolute (the rst's cycle + a few).  Pass outputs off as k=v (e.g. output5/dt=1e30); keep the hst
(dcycle 1) and one bin dump at NLIM.  env: ATHENAK = athenak tree (for tests/vet_clean/vett.py), LAUNCH =
launcher prefix (e.g. "srun -n 8 --gpus-per-task=1"), NT = OpenMP threads (default 1).
Writes DECK.named next to DECK and prints the table rows; nothing else is changed."""
import os, re, shlex, subprocess, sys
sys.path.insert(0, os.path.join(os.environ.get('ATHENAK', '.'), 'tests', 'vet_clean'))
import vett
exe, rst, deck, nlim = sys.argv[1:5]; kv = sys.argv[5:] + ['time/nlim=' + nlim]
rst, deck = os.path.abspath(rst), os.path.abspath(deck)
L = shlex.split(os.environ.get('LAUNCH', ''))
def run(sub, args):
    d = os.path.abspath(sub); os.makedirs(d, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get('NT', '1'))
    with open(os.path.join(d, 'run.log'), 'w') as f:
        rc = subprocess.call(L + [exe] + args, cwd=d, stdout=f, stderr=subprocess.STDOUT, env=env)
    log = open(os.path.join(d, 'run.log')).read()
    print(sub, 'rc', rc, 'CLEAN WARNING', log.count('CLEAN WARNING'), flush=True)
    return d, log
a, la = run('a_deck', ['-r', rst, '-i', deck] + kv)
pat = re.compile(r'CLEAN WARNING restart: (\w+)/(\w+) = (\S+) from the restart file \(code default (\S+);')
keys = pat.findall(la)
for b, k, v, dflt in keys:
    print('  stored  %s/%s = %s   (code default %s)' % (b, k, v, dflt))
# DECK.named: each listed key appended to its block at the stored value
txt, out, blk = open(deck).read().split('\n'), [], None
add = {}
for b, k, v, _ in keys: add.setdefault(b, []).append((k, v))
def flush(b):
    if b not in add: return []
    return ['# named at the restart file\'s stored value (rstkeys.py)'] + ['%s = %s' % kv_ for kv_ in add.pop(b)]
for ln in txt:
    s = ln.split('#')[0].strip()
    if s.startswith('<'):
        j = len(out)
        while j and not out[j-1].strip(): j -= 1
        out[j:j] = flush(blk); blk = s[1:-1]
    out.append(ln)
out += flush(blk)
for b in list(add): out += ['', '<%s>' % b] + flush(b)
nd = deck + '.named'; open(nd, 'w').write('\n'.join(out) + '\n')
b_, lb = run('b_named', ['-r', rst, '-i', nd] + kv)
c_, lc = run('c_ronly', ['-r', rst] + kv)
nf, mx, _ = vett.compare_dirs(b_, c_, last_only=True)
print('ROW %s | warnings with DECK %d | with DECK.named %d | b_named vs c_ronly: %d files, maxdiff %g'
      % (os.path.basename(deck), len(keys), lb.count('CLEAN WARNING restart:'), nf, mx))
