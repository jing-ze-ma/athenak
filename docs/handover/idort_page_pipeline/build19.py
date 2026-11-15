#!/usr/bin/env python3
"""idort-1019 page v2 (USER 10-20: mirror the IDORT paper sections, one large figure per paper figure, no thumbnails).
Reads content.json (summary, conv_table, regress, flags, keys) + fig/*.png; writes a self-contained HTML."""
import json, base64, html, os, sys
D = os.path.dirname(os.path.abspath(__file__)); C = json.load(open(D + '/content.json'))
OUT = sys.argv[1] if len(sys.argv) > 1 else '/raven/ptmp/jinma/idort_page/publish_1020/index.html'
e = lambda s: html.escape('' if s is None else str(s))
def img(f):
    p = os.path.join(D, 'fig', f)
    if not os.path.exists(p): return '<p class="dim">(figure pending)</p>'
    return '<img alt="%s" src="data:image/png;base64,%s">' % (e(f), base64.b64encode(open(p, 'rb').read()).decode())
def badge(v):
    v = (v or '').upper(); k = v.split()[0] if v else ''
    cls = {'PASS': 'ok', 'FAIL': 'bad', 'OPEN': 'warn', 'FLAG': 'warn', 'PENDING': 'pend', 'EXCLUDED': 'excl', 'INFO': 'muted'}.get(k, 'muted')
    return '<span class="badge %s">%s</span>' % (cls, e(v) or 'n/a')
def table(head, rows, raw=()):
    h = ''.join('<th>%s</th>' % e(x) for x in head)
    b = ''.join('<tr>%s</tr>' % ''.join('<td>%s</td>' % (c if i in raw else e(c)) for i, c in enumerate(r)) for r in rows)
    return '<div class="tbl"><table><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>' % (h, b)
SUMKEY = {}
for r in C['summary']: SUMKEY.setdefault(r['test'], r)
def sec(sid, title, ref, figs, text, tests=()):
    rows = [r for r in C['summary'] if any(t in r['test'] or t == r['geom'] for t in tests)]
    t = ''.join('<p class="res">%s <b>%s</b> (%s): %s. <span class="dim">Criterion: %s.</span></p>' % (badge(r['verdict']), e(r['test']), e(r['geom']), e(r['value']), e(r['ref'])) for r in rows)
    f = ''.join('<figure>%s<figcaption>%s</figcaption></figure>' % (img(fn), e(cap)) for fn, cap in figs)
    return '<section id="%s"><h2>%s</h2><p class="ref">%s</p><p>%s</p>%s%s</section>' % (sid, e(title), e(ref), e(text), t, f)
S = C['sections']
css = open(D + '/page.css').read() + '''
section { display:grid; gap:10px }
figure { margin:8px 0 18px } figure img { width:100%; max-width:1100px; height:auto; display:block; border:1px solid var(--line); border-radius:4px; background:#fff }
figcaption { font-size:.9rem; color:var(--dim); max-width:1100px }
.ref { font-family:var(--mono); font-size:.8rem; color:var(--accent) } .res { margin:2px 0 }
.wrap { max-width:1180px }'''
m = C['meta']
summ = table(['Test', 'Geometry', 'Verdict', 'Order (finest pair) / metric'],
             [[r['test'], r['geom'], badge(r['verdict']), r['value']] for r in C['summary']], raw=(2,))
conv = table(['Case', 'Geometry', 'Regime', 'Var', 'Levels', 'Pair orders', 'Status'],
             [[r['case'], r['geom'], r['regime'], r['var'], r['N'], r['pairs'], badge(r['status'])] for r in C['conv_table']], raw=(6,))
cnt = {}
for r in C['summary']: cnt[r['verdict'].split()[0]] = cnt.get(r['verdict'].split()[0], 0) + 1
nav = ''.join('<a href="#%s">%s</a>' % (s['id'], e(s['nav'])) for s in S)
body = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>IDORT Suite PROD</title><style>{css}</style></head><body><div class="wrap">
<div><h1>IDORT Test Suite: production keys</h1><p class="lede">{e(C["lede"])}</p></div>
<div class="summary">{' '.join('<span>%s %d</span>' % (badge(k), v) for k, v in sorted(cnt.items()))}</div>
<div class="card"><p><span class="lbl">Binary</span> <code>{e(m["binary"])}</code></p><p><span class="lbl">Commit</span> {e(m["commit"])}</p>
<p><span class="lbl">Built</span> {e(m["built"])}</p><p><span class="lbl">Runs</span> {e(m["runs"])}</p><p><span class="lbl">PROD arm</span> {e(m["arm"])}</p></div>
<nav class="nav"><a href="#summ">Summary</a>{nav}<a href="#orders">Order table</a><a href="#reg">Flags</a><a href="#diag">Diagnostics</a><a href="#keys">Keys</a></nav>
<h2 id="summ">Summary</h2><p class="dim">{e(C["summary_note"])}</p>{summ}
{''.join(sec(s['id'], s['title'], s['ref'], s['figs'], s['text'], s.get('tests', ())) for s in S)}
<h2 id="orders">Convergence orders (PROD)</h2><p class="dim">{e(C["conv_note"])}</p>{conv}
<h2 id="reg">Flags and diagnoses</h2>
{('<ul>' + ''.join('<li>%s</li>' % e(x) for x in C['regress']) + '</ul>') if C['regress'] else ''}
{''.join('<p><b>%s.</b> %s</p>' % (e(f['title']), e(f['text'])) for f in C['flags'])}
<h2 id="diag">Additional diagnostics (not shown)</h2>
{table(['Item', 'Path'], C['diag'])}
<h2 id="keys">Production method keys (PROD arm)</h2>
<p class="dim">{e(C["keys_note"])}</p>
<div class="tbl"><table><thead><tr><th>Key</th><th>Value</th></tr></thead><tbody>{''.join('<tr><td><code>%s</code></td><td><code>%s</code></td></tr>' % (e(k), e(v)) for k, v in C['keys'])}</tbody></table></div>
<p class="dim">{e(C["drops_note"])}</p>
</div></body></html>'''
os.makedirs(os.path.dirname(OUT), exist_ok=True); open(OUT, 'w').write(body); print('wrote', OUT, len(body))
