"""make_deck.py (dhj_planets_1101): planet deck from the WASP-121b hydro fresh deck
w121hfresh_1101/deck/hfresh.athinput.  Every <block>/key in OVR (a JSON file of
{"block/key": [value, comment]}) replaces the value (the line gets '# <tag>: <comment>
(W121: <old value>)'); keys not in the template are FATAL (a typo must not pass silently),
except those listed under "_add" ({"block/key": [value, comment]}, appended to the block).
usage: make_deck.py OVR.json OUT.athinput TAG "header line" ["header line" ...]
"""
import json
import re
import sys

TEMPLATE = '/viper/ptmp2/jinma/w121hfresh_1101/deck/hfresh.athinput'


def main():
    ovr = json.load(open(sys.argv[1]))
    out, tag, hdr = sys.argv[2], sys.argv[3], sys.argv[4:]
    add = ovr.pop('_add', {})
    lines = open(TEMPLATE).read().split('\n')
    # drop the template's WASP-121b header comments (everything before the first PAR_DUMP)
    lines = lines[next(i for i, s in enumerate(lines) if 'PAR_DUMP' in s):]
    block, done, res = None, set(), []
    for h in hdr:
        res.append('# ' + h)
    res.append('# Template: ' + TEMPLATE + ' (WASP-121b hfresh, fork/fin51-r-fix 182477553).')
    res.append('# Every changed key is marked "# %s:"; all other keys are the template values.'
               % tag)

    def flush_add(b):
        for k, (v, c) in add.items():
            kb, kk = k.split('/')
            if kb == b:
                res.append('%-26s = %-14s # %s: %s (NEW key, not in the template)'
                           % (kk, v, tag, c))
                done.add(k)
    for ln in lines:
        m = re.match(r'^<(\w+)>', ln)
        if m:
            if block is not None:
                flush_add(block)
            block = m.group(1)
            res.append(ln)
            continue
        m = re.match(r'^(\s*)([A-Za-z0-9_]+)(\s*)=(\s*)([^#]*?)(\s*)(#.*)?$', ln)
        if block is not None and m and not ln.lstrip().startswith('#'):
            k = '%s/%s' % (block, m.group(2))
            if k in ovr:
                v, c = ovr[k]
                old = m.group(5).strip()
                res.append('%-26s = %-14s # %s: %s (W121: %s)'
                           % (m.group(2), v, tag, c, old if old else '<empty>'))
                done.add(k)
                continue
        res.append(ln)
    miss = [k for k in list(ovr) + list(add) if k not in done]
    if miss:
        sys.exit('FATAL: keys not found in the template: %s' % miss)
    open(out, 'w').write('\n'.join(res))
    print('wrote', out, len(done), 'keys changed/added')


if __name__ == '__main__':
    main()
