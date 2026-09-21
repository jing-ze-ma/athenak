#!/usr/bin/env python3
"""Fetch TOPS/ATOMIC gray (Rosseland AND Planck) mean opacities for a user
mixture from the LANL web service, non-interactively.

The service is a two-step POST (no cookies, no captcha, no registration):
  1) POST https://aphysics2.lanl.gov/submit   -> a summary page
  2) POST https://aphysics2.lanl.gov/results  with output=tabcol -> the numbers
All of the form's text fields must be present in step 1, even the ones the
selected radio buttons do not use (omitting `temps`, `dens` or `energies`
makes the server answer 500).

usage: fetch_tops.py "<mixture>" <tlow_keV> <tup_keV> <rlow> <rup> <nr> <out.txt>
"""
import re
import sys
import html
import time
import urllib.parse
import urllib.request

SUBMIT = 'https://aphysics2.lanl.gov/submit'
RESULTS = 'https://aphysics2.lanl.gov/results'
UA = 'Mozilla/5.0 (research script; contact via LANL opacity@lanl.gov)'


def post(url, fields, referer, timeout=900):
    data = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(url, data=data,
                                 headers={'User-Agent': UA,
                                          'Referer': referer})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')


def query(mixture, tlow, tup, rlow, rup, nr, nthin=1, lib='any',
          fractype='atomic', mixname='m1opac'):
    f = dict(userid='ajn', pageid='formpage', idnumber='19000', numhits='1',
             unique='', lib=lib, fractype=fractype, mixture=mixture,
             mixname=mixname, isotope='noniso',
             tgrid='range', temps='1.', tlow=str(tlow), tup=str(tup),
             nthin=str(nthin),
             rgrid='range', dens='1.', rlow=repr(rlow), rup=repr(rup),
             nr=str(nr), rspace='log',
             egrid='range', energies='.1 1. 10.', egplow='.001',
             egphigh='300.', ngpengs='33', espace='log',
             plasnu='on', datype='gray')
    summary = post(SUBMIT, f, 'https://aphysics2.lanl.gov/apps/')
    m = re.search(r'<form[^>]*action="\./results"[^>]*>(.*?)</form>',
                  summary, re.S)
    if m is None:
        raise RuntimeError('no results form in summary page:\n'
                           + re.sub(r'<[^>]+>', ' ', summary)[:2000])
    g = dict(re.findall(r'name="([^"]*)" value="([^"]*)"', m.group(1)))
    g['output'] = 'tabcol'
    return html.unescape(re.sub(r'<[^>]+>', '', post(RESULTS, g, SUBMIT)))


def parse(text):
    """-> (list of T[keV], list of rho[g/cc], kR[nT][nD], kP[nT][nD])."""
    text = text.replace(u'\xa0', ' ')
    T, rho, kR, kP = [], [], [], []
    cur = None
    for ln in text.splitlines():
        m = re.match(r'\s*Density\s+Ross opa\s+Planck opa.*T=\s*([0-9.eE+-]+)',
                     ln)
        if m:
            T.append(float(m.group(1)))
            cur = ([], [])
            kR.append(cur[0])
            kP.append(cur[1])
            continue
        if cur is not None:
            v = ln.split()
            if len(v) == 5 and all(re.match(r'^[0-9.eE+-]+$', x) for x in v):
                if len(T) == 1:
                    rho.append(float(v[0]))
                cur[0].append(float(v[1]))
                cur[1].append(float(v[2]))
            elif ln.strip() and not ln.strip().startswith('Density'):
                cur = None
    return T, rho, kR, kP


if __name__ == '__main__':
    mix, tl, tu, rl, ru, nr, out = sys.argv[1:8]
    t0 = time.time()
    txt = query(mix, float(tl), float(tu), float(rl), float(ru), int(nr))
    open(out + '.raw', 'w').write(txt)
    T, rho, kR, kP = parse(txt)
    print('%.1f s; nT=%d nrho=%d' % (time.time()-t0, len(T), len(rho)))
    with open(out, 'w') as fo:
        fo.write('# TOPS/ATOMIC gray means, mixture: %s\n' % mix)
        fo.write('# T[keV] rho[g/cc] kappa_R[cm2/g] kappa_P[cm2/g]\n')
        for i, t in enumerate(T):
            for j, d in enumerate(rho):
                fo.write('%.6e %.6e %.6e %.6e\n' % (t, d, kR[i][j], kP[i][j]))
    print('wrote', out)
