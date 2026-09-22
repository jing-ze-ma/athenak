"""Instruction / scratch / global-load counts per rt_chain_ck instantiation.

Usage: python3 count.py <disassembly.txt>

Reads an llvm-objdump -d of the unbundled gfx942 code object (see dis.sh) and prints one
line per correlated-k chain-kernel instantiation, keyed by its compile-time tags: the
radial tier NN and then SPH, BSP, CCH, FRM (and FOP where the build has it).

The tags are read out of the RAW Itanium name, not out of c++filt, because repeated
tag types are substitution-compressed (`SB_`) and c++filt then prints the wrong one.
Only this call's own arguments can appear as substitutions inside the argument list,
and they are registered from the index after the template name's, upward in order of
first appearance; tags() resolves them.
"""
import sys
import re
from collections import defaultdict

ARG = re.compile(r'_clI(.*?)EEDaS2_')
TOK = re.compile(r'St17integral_constantI([ib])L[ib](-?\d+)EE'
                 r'|S[0-9A-Z]*_I([ib])L[ib](-?\d+)EE'
                 r'|S([0-9A-Z]*)_')


def subid(txt):
    return 0 if txt == '' else int(txt, 36) + 1


def tags(name):
    """the five/six tag values of one instantiation, as strings"""
    m = ARG.search(name)
    if not m:
        return None
    out, reg = [], {}
    # the substitution index of std::integral_constant itself; the full tag types are
    # registered from the next index upward, in order of first appearance
    tm = re.search(r'S([0-9A-Z]*)_I[ib]L[ib]', m.group(1))
    base = (subid(tm.group(1)) + 1) if tm else 11
    for t in TOK.finditer(m.group(1)):
        if t.group(2) is not None or t.group(4) is not None:
            v = t.group(2) if t.group(2) is not None else t.group(4)
            k = t.group(1) if t.group(1) is not None else t.group(3)
            val = ('true' if v == '1' else 'false') if k == 'b' else v
            reg[base + len(reg)] = val
            out.append(val)
        else:
            out.append(reg.get(subid(t.group(5)), '?'))
    return out


lines = open(sys.argv[1]).read().split('\n')
cur = None
cnt = defaultdict(int)
scr = defaultdict(int)
gl = defaultdict(int)
for ln in lines:
    m = re.match(r'^[0-9a-f]+ <(.+)>:', ln)
    if m:
        cur = m.group(1)
        continue
    if cur is None:
        continue
    m2 = re.match(r'^\t([a-z][a-z0-9_]*)', ln)
    if m2:
        cnt[cur] += 1
        o = m2.group(1)
        if o.startswith('scratch_'):
            scr[cur] += 1
        if o.startswith('global_load'):
            gl[cur] += 1

rows = []
for k, v in cnt.items():
    if 'T0_T1_T2_T3_' not in k or v < 2000:
        continue
    tg = tags(k)
    if tg:
        rows.append((tg, v, scr[k], gl[k]))
for tg, c, s, g in sorted(rows, key=lambda r: (r[0][1:], int(r[0][0]))):
    print('tier=%-4s SPH=%-5s BSP=%-5s CCH=%-2s FRM=%-2s FOP=%-5s '
          'instr=%-7d scratch=%-6d gload=%d'
          % (tg[0], tg[1], tg[2], tg[3], tg[4],
             (tg[5] if len(tg) > 5 else '-'), c, s, g))
print('# %d chain instantiations' % len(rows))
