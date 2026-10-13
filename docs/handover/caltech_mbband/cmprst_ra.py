# Caltech: cmprst.py made aware of the rotavg-rst2-1013 (9972a56d) restart change, which is in NEW (mbband-1013)
# but not in OLD (281dd686): the <rad_m1>/vet_rotate_avg header (magic M1ROTAV1) has 11 int64 entries in NEW
# (+ i0) vs 10 in OLD, and NEW stores the history only for i >= i0 (below: zero, never read).
# python3 cmprst_ra.py <rst OLD> <rst NEW> <deck>  ->  first token:
#   BITWISE  data after <par_end> identical once the RA header's extra entry is dropped (first 10 entries equal)
#   RAHIST   all bytes before the vet_rotate_avg history (the last nra slabs; NEW may store fewer, i >= i0) are identical
#            and the RA header agrees except nra (+ i0); i.e. all state is bitwise
#   DIFF     anything else (state differs)
import re, struct, sys
def load(f):
    b = open(f, 'rb').read(); i = b.find(b'<par_end>'); b = b[b.index(b'\n', i) + 1:]
    m = b.find(b'M1ROTAV1')
    if m < 0 or m > 20000: return b, None, b''
    nb = struct.unpack('<Q', b[m + 8:m + 16])[0]
    return b[:m], struct.unpack('<%dq' % (nb // 8), b[m + 16:m + 16 + nb]), b[m + 16 + nb:]
def slab(deck):
    blk, v = None, {}
    for ln in open(deck):
        s = ln.split('#')[0].strip(); mm = re.match(r'<(\w+)>', s)
        if mm: blk = mm.group(1); continue
        if '=' in s: k, x = [t.strip() for t in s.split('=', 1)]; v[(blk, k)] = x
    ng = int(v.get(('mesh', 'nghost'), 2)); n = 1; nmb = 1
    for d in ('nx1', 'nx2', 'nx3'):
        nx = int(v[('mesh', d)]); n *= (nx + 2 * ng) if nx > 1 else 1
        nmb *= nx // int(v.get(('meshblock', d), nx))
    return n * nmb * 8
a, b, deck = sys.argv[1:4]
pa, ha, da = load(a); pb, hb, db = load(b)
if pa == pb and da == db and ((ha is None and hb is None) or (ha is not None and hb is not None and ha[:10] == hb[:10])):
    print('BITWISE', 'ra_hdr', ha, hb); sys.exit()
if pa != pb or ha is None or hb is None or ha[:1] + ha[2:10] != hb[:1] + hb[2:10]:
    print('DIFF', 'pre_eq', pa == pb, 'ra_hdr', ha, hb, 'len', len(da), len(db)); sys.exit()
# state = everything before the RA history (nra slabs at the end; NEW packs only i >= i0, so its nra can be smaller)
sl = slab(deck); ta = len(da) - ha[1] * sl; tb = len(db) - hb[1] * sl
if ta == tb and ta > 8 and da[8:ta] == db[8:tb]:   # bytes 0-7 = data_size (includes nra)
    print('RAHIST', 'state_bytes', ta, 'ra_slabs', ha[1], hb[1], 'len', len(da), len(db), 'ra_hdr', ha, hb); sys.exit()
k = next((i for i in range(8, min(len(da), len(db))) if da[i] != db[i]), None)
print('DIFF', 'first_diff', k, 'state_bytes', ta, tb, 'len', len(da), len(db), 'ra_hdr', ha, hb)
