#!/usr/bin/env python3
"""Print title + lines with mass/radius/separation numbers for a list of text dumps."""
import re, sys, os
d = sys.argv[1]
pat = re.compile(r"(M⊙|R⊙|Porb|P ?\(d\)|P ?\[d\]|yr−1)")
for f in sys.argv[2:]:
    p = os.path.join(d, f + ".txt")
    if not os.path.exists(p):
        print("####", f, "MISSING"); continue
    L = open(p).read().split("\n")
    print("####", f, " ".join(L[:3])[:150])
    n = 0
    for i, l in enumerate(L):
        if pat.search(l):
            print(f"{i+1}: {l[:160]}"); n += 1
            if n >= int(os.environ.get("NMAX", "14")): break
