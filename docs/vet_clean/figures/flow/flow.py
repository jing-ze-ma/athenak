#!/usr/bin/env python3
"""The flow charts of docs/vet_clean: the .mmd files here are the source (Mermaid flowcharts).

  flow.py --check [ROOT]    every function of a box's second line is defined in the file named
                            there, every file exists, at most 20 boxes; the names in backticks of
                            WORKFLOW.md exist.  Exit status 1 and a list if not.
  flow.py --render          X.mmd -> X.svg with graphviz (dot); the site has no Mermaid.
"""
import html
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WIDTH = 1000          # px, the widest chart


def read(f):
    nodes, edges = {}, []
    for line in open(f, encoding='utf-8'):
        s = line.strip()
        if not s or s.startswith('%%') or s.startswith('flowchart'):
            continue
        m = re.match(r'^(\w+)(\[|\{)"(.*)"(\]|\})$', s)
        if m:
            nodes[m.group(1)] = (m.group(3).split('<br/>'), m.group(2) == '{')
            continue
        part = re.split(r'\s*-->(?:\|([^|]*)\|)?\s*', s)
        ids, labs = part[0::2], part[1::2]
        if len(ids) < 2 or any(i not in nodes for i in ids):
            sys.exit('%s: cannot read %r' % (f, s))
        edges += [(a, b, w or '') for a, b, w in zip(ids[:-1], ids[1:], labs)]
    return nodes, edges


def src(root, name):
    for d in ('src/rad_vet/', 'src/'):
        if os.path.exists(os.path.join(root, d, name)):
            return os.path.join(root, d, name)
    return None


def check(root):
    bad = []
    for f in sorted(x for x in os.listdir(HERE) if x.endswith('.mmd')):
        nodes, edges = read(os.path.join(HERE, f))
        if len(nodes) > 20:
            bad.append('%s: %d boxes' % (f, len(nodes)))
        for k, (lab, _) in nodes.items():
            if len(lab) != 2:
                bad.append('%s %s: not two lines' % (f, k))
                continue
            w = [x.strip() for x in lab[1].split(',')]
            path = src(root, w[-1])
            if path is None:
                bad.append('%s %s: no file %s' % (f, k, w[-1]))
                continue
            text = open(path, encoding='utf-8').read()
            for fn in w[:-1]:
                q = fn if '::' in fn else r'RadiationVet::' + fn
                m = re.search(r'^(?![ /]).*\b%s\(' % q, text, re.M)
                if not m:
                    bad.append('%s %s: %s is not defined in %s' % (f, k, fn, w[-1]))
                elif '--list' in sys.argv:
                    print('%-8s %-48s %s:%d' % (f[:-4], lab[0], os.path.relpath(path, root),
                                                text.count('\n', 0, m.start()) + 1))
    page = os.path.join(HERE, '..', '..', 'WORKFLOW.md')
    if os.path.exists(page):
        allsrc = ''
        for d, _, fs in os.walk(os.path.join(root, 'src')):
            allsrc += ''.join(open(os.path.join(d, x), encoding='utf-8', errors='replace').read()
                              for x in fs if x.endswith(('.cpp', '.hpp', '.inc')) and 'kokkos' not in d)
        for name in re.findall(r'`([^`\n]+)`', open(page, encoding='utf-8').read()):
            if re.match(r'^[\w/]+\.(cpp|hpp|inc)$', name):
                if src(root, name) is None:
                    bad.append('WORKFLOW.md: no file %s' % name)
            elif re.match(r'^[A-Z]\w*[a-z]\w*(::\w+)?$', name):
                if not re.search(r'\b%s\b' % re.escape(name), allsrc):
                    bad.append('WORKFLOW.md: no name %s' % name)
            elif re.match(r'^(clean|time2|dt|vet)_\w+$', name):
                if not re.search(r'"%s"' % name, allsrc):
                    bad.append('WORKFLOW.md: no key %s' % name)
    for b in bad:
        print(b)
    return 1 if bad else 0


def render():
    for f in sorted(x for x in os.listdir(HERE) if x.endswith('.mmd')):
        nodes, edges = read(os.path.join(HERE, f))
        g = ['digraph g {', 'graph [bgcolor="white", nodesep=0.28, ranksep=0.34, pad=0.15];',
             'node [fontname="DejaVu Sans", fontsize=12, color="#31586b", fontcolor="#15232b", '
             'style=filled, fillcolor="#eef4f7", margin="0.14,0.07"];',
             'edge [fontname="DejaVu Sans", fontsize=11, color="#4a6674", fontcolor="#15232b", '
             'arrowsize=0.7];']
        for k, (lab, dia) in nodes.items():
            t = '<%s<BR/><FONT POINT-SIZE="10" COLOR="#3d5866">%s</FONT>>' % (
                html.escape(lab[0]), html.escape(lab[1]))
            g.append('%s [shape=%s, %slabel=%s];' % (k, 'box' if not dia else 'octagon',
                     'fillcolor="#fdf3df", color="#8a6a1f", ' if dia else '', t))
        for a, b, w in edges:
            g.append('%s -> %s%s;' % (a, b, ' [label="  %s "]' % w if w else ''))
        g.append('}')
        svg = subprocess.run(['dot', '-Tsvg'], input='\n'.join(g).encode(), check=True,
                             stdout=subprocess.PIPE).stdout.decode()
        svg = svg[svg.index('<svg'):]
        w, h = (float(x) for x in re.search(r'width="([\d.]+)pt" height="([\d.]+)pt"', svg).groups())
        w, h = w*4/3, h*4/3
        if w > WIDTH:
            w, h = WIDTH, h*WIDTH/w
        svg = re.sub(r'width="[\d.]+pt" height="[\d.]+pt"', 'width="%d" height="%d"' % (w, h), svg, 1)
        svg = re.sub(r'<title>[^<]*</title>\n', '', svg).replace('font-family="DejaVu Sans"',
                     'font-family="Helvetica,Arial,sans-serif"')
        svg = re.sub(r'<!--.*?-->\n?', '', svg, flags=re.S)
        open(os.path.join(HERE, f[:-4] + '.svg'), 'w').write(
            '<?xml version="1.0" encoding="UTF-8"?>\n' + svg)
        print('%s.svg %d x %d px, %d boxes' % (f[:-4], w, h, len(nodes)))


if __name__ == '__main__':
    if '--render' in sys.argv:
        render()
    else:
        root = [a for a in sys.argv[1:] if not a.startswith('-')]
        sys.exit(check(root[0] if root else os.path.join(HERE, '../../../..')))
