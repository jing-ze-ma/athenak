"""mkinp.py <base> <out> key=val ... : copy an athinput, set/insert <block>/key = val
(block/key=val), and append hst + rst outputs (every cycle hst, rst at nlim)."""
import sys
src, out = sys.argv[1], sys.argv[2]
sets = {}
for a in sys.argv[3:]:
    k, v = a.split("=", 1)
    b, kk = k.split("/")
    sets.setdefault(b, {})[kk] = v
lines = open(src).read().splitlines()
res, blk, done = [], None, set()


def flush(b):
    for kk, v in sets.get(b, {}).items():
        if (b, kk) not in done:
            res.append("%s = %s" % (kk, v))
            done.add((b, kk))


for ln in lines:
    s = ln.strip()
    if s.startswith("<") and s.endswith(">"):
        if blk:
            flush(blk)
        blk = s[1:-1]
        res.append(ln)
        continue
    key = s.split("=")[0].strip() if "=" in s and not s.startswith("#") else None
    if blk in sets and key in sets[blk]:
        res.append("%s = %s" % (key, sets[blk][key]))
        done.add((blk, key))
        continue
    res.append(ln)
flush(blk)
for b in sets:
    if not any(bb == b for bb, _ in done):
        res.append("<%s>" % b)
        flush(b)
res += ["<output91>", "file_type = hst", "dcycle = 1",
        "<output92>", "file_type = rst", "dcycle = 1000000"]
open(out, "w").write("\n".join(res) + "\n")
