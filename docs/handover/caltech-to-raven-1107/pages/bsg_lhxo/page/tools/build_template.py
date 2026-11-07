"""Build pipe/index.html for the Caltech FIN3 page: Caltech placeholder template (bsg_1009 pagework)
+ the prose of the Raven-published FIN3 page (10-13, t 20.6 d), with Caltech hardware/footer."""
import sys
cal = open(sys.argv[1]).read().split('\n')     # Caltech template (placeholders)
rav = open(sys.argv[2]).read().split('\n')     # Raven published body (wrapper stripped)
out = sys.argv[3]
# (Caltech 1-based range) -> (Raven 1-based range), from diff
hunks = [((94, 96), (94, 96)), ((114, 114), (113, 113)), ((124, 124), (123, 123)),
         ((129, 129), (128, 128)), ((132, 132), (131, 131)), ((139, 141), (138, 140)),
         ((280, 280), (275, 275)), ((285, 285), (280, 280)), ((309, 309), (305, 305))]
res = list(cal)
for (a, b), (c, d) in sorted(hunks, reverse=True):
    res[a - 1:b] = rav[c - 1:d]
h = '\n'.join(res)
def sub(old, new, n=1):
    global h
    assert h.count(old) == n, (old, h.count(old))
    h = h.replace(old, new)
sub('GPUs (1 Caltech node &times; 2 H200; 16 MeshBlocks, 8 per GPU)',
    'GPUs (1 Caltech node &times; 2 H200; 16 MeshBlocks, 8 per GPU)')
sub('final method FIN3 at 9972a56d, Raven)', 'final method FIN3 at 9972a56d, Caltech)')
sub('16 MeshBlocks of 256 &times; 64 &times; 64, 1 per GPU)', '16 MeshBlocks of 256 &times; 64 &times; 64, 8 per GPU)')
i = h.index('<span>Our run: AthenaK on 16 A100')
j = h.index('</span>', i)
h = h[:i] + ('<span>Our run: AthenaK on 2 H200 (1 Caltech node, chained 11.75-h links; Caltech copy of the '
             'production run bsgP), the paper\'s IC, opacity tables and boundary rules, 10<sup>&minus;6</sup> '
             'temperature seed, no well-balancing, the final implicit two-moment radiation method FIN3 with a '
             'short-characteristics variable Eddington tensor averaged over the last 4 rotated direction sets, '
             'hesdirk2 time integration and the conservative2 gas source (rotavg-rst2-1013 9972a56d, deck '
             'bsgP_seed_fin3.athinput), ideal gas &mu; = 0.62, CFL 0.3. Regenerated on Caltech from the Caltech '
             'run (bsgprod_1012/prod) with viper\'s page pipeline (bsg-files-1009 docs/handover/bsg-hrdet-1009/page) '
             'with Caltech\'s edits (bsgprod-1012/pagework) via <code>bsgprod_1012/pagework_fin3/make_page.sh</code>.'
             ) + h[j:]
for bad in ('Raven', 'A100'):
    for k in range(len(h)):
        pass
open(out, 'w').write(h)
import re
for m in re.finditer(r'[^.<>]{0,80}(Raven|A100)[^.<>]{0,80}', h):
    print('NOTE mentions:', m.group(0)[:170])
