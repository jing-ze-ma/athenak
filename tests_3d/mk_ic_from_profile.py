#!/usr/bin/env python3
"""Build a red_giant problem/ic_profile file from the LAST record of an rt_profile.bin.

The 1-D relaxation run writes shell-averaged radial profiles (record layout in
src/pgen/red_giant.cpp:494-506); on a 1-D column the shell mean IS the column.  Slots
0 and 6 are rho and eint, and x1v is the stretched grid's own cell centres, so the last
record is exactly the three columns "r rho eint" the ic_profile reader wants.

rt_profile only covers the ACTIVE cells, while the reader demands a file that spans the
mesh AND its radial ghosts.  The tails are therefore taken from the ORIGINAL ic file,
each multiplied by the constant ratio that makes it continuous with the relaxed profile
at the first/last active node, so the extrapolation keeps the original slope.

  usage: mk_ic_from_profile.py rt_profile.bin ic_original.txt ic_out.txt [--record N]
"""
import sys
import struct
import numpy as np


def read_profiles(path):
    recs = []
    with open(path, 'rb') as f:
        blob = f.read()
    off = 0
    while off < len(blob):
        t, = struct.unpack_from('<d', blob, off)
        off += 8
        n1, nv = struct.unpack_from('<ii', blob, off)
        off += 8
        need = 8*n1 + 8*nv*n1
        if off + need > len(blob):
            break
        x1v = np.frombuffer(blob, '<f8', n1, off)
        off += 8*n1
        q = np.frombuffer(blob, '<f8', nv*n1, off).reshape(nv, n1)
        off += 8*nv*n1
        recs.append((t, x1v.copy(), q.copy()))
    return recs


def main():
    prof, ic_in, ic_out = sys.argv[1], sys.argv[2], sys.argv[3]
    irec = -1
    if '--record' in sys.argv:
        irec = int(sys.argv[sys.argv.index('--record')+1])
    recs = read_profiles(prof)
    if not recs:
        sys.exit('no complete record in %s' % prof)
    t, r, q = recs[irec]
    rho, eint = q[0], q[6]
    if not np.all(np.isfinite(rho)) or np.any(rho <= 0.0):
        sys.exit('record %d has non-positive or non-finite rho' % irec)

    old = np.loadtxt(ic_in)
    ro, rhoo, eio = old[:, 0], old[:, 1], old[:, 2]
    lo = ro < r[0]
    hi = ro > r[-1]
    parts_r = []
    for msk, iend in ((lo, 0), (hi, -1)):
        if not msk.any():
            continue
        fd = rho[iend]/np.interp(r[iend], ro, rhoo)
        fe = eint[iend]/np.interp(r[iend], ro, eio)
        parts_r.append((ro[msk], rhoo[msk]*fd, eio[msk]*fe, iend))
    out_r = [p[0] for p in parts_r if p[3] == 0] + [r] + \
            [p[0] for p in parts_r if p[3] == -1]
    out_d = [p[1] for p in parts_r if p[3] == 0] + [rho] + \
            [p[1] for p in parts_r if p[3] == -1]
    out_e = [p[2] for p in parts_r if p[3] == 0] + [eint] + \
            [p[2] for p in parts_r if p[3] == -1]
    R = np.concatenate(out_r)
    D = np.concatenate(out_d)
    E = np.concatenate(out_e)
    assert np.all(np.diff(R) > 0.0), 'r is not strictly ascending'

    with open(ic_out, 'w') as f:
        f.write('# red_giant problem/ic_profile: r[cm] rho[g/cm^3] eint[erg/cm^3]\n')
        f.write('# THE RELAXED 1-D COLUMN, not the MESA structure.  Built by\n')
        f.write('#   tests_3d/mk_ic_from_profile.py %s %s %s\n' % (prof, ic_in, ic_out))
        f.write('# from record %d of that dump, t = %.6e s.  The %d active nodes over\n'
                % (irec, t, len(r)))
        f.write('#   %.6e .. %.6e cm are the run\'s own shell means (slots 0 and 6 of\n'
                % (r[0], r[-1]))
        f.write('# the rt_profile record); outside them the ORIGINAL file %s is used,\n'
                % ic_in)
        f.write('# rescaled by the constant ratio that makes each tail continuous with\n')
        f.write('# the relaxed profile at the first/last active node, so the radial\n')
        f.write('# ghosts keep the original slope.\n')
        for a, b, c in zip(R, D, E):
            f.write('%.10e %.10e %.10e\n' % (a, b, c))
    print('wrote %s: %d nodes, %.6e .. %.6e cm, from t = %.6e s'
          % (ic_out, len(R), R[0], R[-1], t))


if __name__ == '__main__':
    main()
