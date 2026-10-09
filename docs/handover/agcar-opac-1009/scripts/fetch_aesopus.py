#!/usr/bin/env python3
"""Submit the AESOPUS 2.1 web form (stev.oapd.inaf.it/cgi-bin/aesopus) non-interactively.
usage: fetch_aesopus.py form.html out_prefix lgtmin lgtmax dlgt1 lgtchange dlgt2 lgrmin lgrmax dlgr X Zref
All form fields are posted with their defaults (radio: the checked one) except the overrides;
solmix=3 (GS98), scaled solar (ialpha 0), mass-fraction normalisation irrelevant for scaled solar."""
import re, sys, ssl, urllib.request, uuid
form, out = sys.argv[1], sys.argv[2]
keys = ['lgtmin', 'lgtmax', 'dlgt1', 'lgtchange', 'dlgt2', 'lgrmin', 'lgrmax', 'dlgr', 'xhmin', 'zeta_ref']
ov = dict(zip(keys, sys.argv[3:13]))
ov['solmix'] = '3'
h = open(form).read()
fields = []
for tag in re.findall(r'<input[^>]*>', h, re.I):
    n = re.search(r'name="([^"]*)"', tag); v = re.search(r'value="([^"]*)"', tag)
    t = re.search(r'type="([^"]*)"', tag)
    if not n: continue
    n = n.group(1); v = v.group(1) if v else ''; t = t.group(1).lower() if t else 'text'
    if t == 'submit' and n != 'submit_form': continue
    if t == 'radio':
        if n in ov: continue
        if 'checked' not in tag: continue
    if v == '""': v = ''
    fields.append((n, ov.get(n, v)))
for k in ('solmix',):
    fields.append((k, ov[k]))
b = uuid.uuid4().hex
body = b''
for n, v in fields:
    body += ('--%s\r\nContent-Disposition: form-data; name="%s"\r\n\r\n%s\r\n' % (b, n, v)).encode()
body += ('--%s--\r\n' % b).encode()
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
req = urllib.request.Request('https://stev.oapd.inaf.it/cgi-bin/aesopus_2.1', data=body,
      headers={'Content-Type': 'multipart/form-data; boundary=' + b, 'User-Agent': 'research script'})
with urllib.request.urlopen(req, timeout=3000, context=ctx) as r:
    txt = r.read().decode('utf-8', 'replace')
open(out + '.html', 'w').write(txt)
print(len(txt)); print(re.findall(r'href="([^"]*)"', txt)[:20])
