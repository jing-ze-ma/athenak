#!/usr/bin/env python3
"""Export the figures listed in docs/vet_clean/figures/FIGURES.md from the executed notebooks.

    python3 docs/vet_clean/figures/export_figures.py      # from the repository root

Each row's `source` is NOTEBOOK#N, the N-th PNG output (1-based, in cell order). Python
standard library only. Stops if a source has no such output or a PNG here is not listed.
"""
import base64
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))


def rows():
    """[(figure, source, tests, caption)] of the table in FIGURES.md."""
    out = []
    for line in open(os.path.join(HERE, 'FIGURES.md'), encoding='utf-8'):
        c = [x.strip() for x in line.strip().strip('|').replace('\\|', '\x00').split('|')]
        if len(c) == 4 and c[0].endswith('.png'):
            out.append(tuple(x.replace('\x00', '|') for x in c))
    return out


def pngs(nb):
    return [o['data']['image/png'] for c in json.load(open(nb, encoding='utf-8'))['cells']
            for o in c.get('outputs', []) if 'image/png' in o.get('data', {})]


def main():
    cache, listed = {}, set()
    for fig, src, _, _ in rows():
        nb, n = src.rsplit('#', 1)
        if nb not in cache:
            cache[nb] = pngs(os.path.join(ROOT, nb))
        if not 1 <= int(n) <= len(cache[nb]):
            sys.exit('%s: %s has %d PNG outputs' % (fig, nb, len(cache[nb])))
        with open(os.path.join(HERE, fig), 'wb') as f:
            f.write(base64.b64decode(cache[nb][int(n) - 1]))
        listed.add(fig)
        print('wrote', fig)
    extra = sorted(set(f for f in os.listdir(HERE) if f.endswith('.png')) - listed)
    if extra:
        sys.exit('not in FIGURES.md: ' + ', '.join(extra))


if __name__ == '__main__':
    main()
