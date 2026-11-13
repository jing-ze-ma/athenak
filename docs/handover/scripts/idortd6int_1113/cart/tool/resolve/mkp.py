#!/usr/bin/env python3
"""mkp.py <src deck> <args.txt> <out deck> [extra k=v ...]: PROD deck = deck <- args <- PROD method keys (prodkeys_all.txt),
dropping only the keys the test geometry/physics forbids (logged to prod_drop.txt next to the out deck):
  - implicit_marshak_face (spherical-polar only) on non-spherical decks
  - implicit_vimp / implicit_cs2_vimp* / time2_enth_be? no: only the vimp pair when the test runs gas_feedback or coupling = false"""
import sys, os, subprocess
T = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, T); import blk
src, args, out = sys.argv[1:4]; extra = sys.argv[4:]
A = open(args).read().split() if os.path.exists(args) else []
PK = open(T + '/prodkeys_all.txt').read().split()
rad = blk.read(src, 'rad_m1', args if os.path.exists(args) else None)
mesh = blk.read(src, 'mesh', args if os.path.exists(args) else None)
sph = mesh.get('use_spherical_polar', 'false') in ('true', '1') or mesh.get('use_cubed_sphere', 'false') in ('true', '1')
nofb = rad.get('gas_feedback', 'true') in ('false', '0') or rad.get('coupling', 'true') in ('false', '0') or rad.get('dbg_gas_force', '1') in ('false', '0')
per1 = mesh.get('ix1_bc', '') == 'periodic'
drop, keep = [], []
for kv in PK:
    k = kv.split('=')[0].split('/')[-1]
    if (k in ('implicit_marshak_face', 'implicit_face_weight') and not sph) or (nofb and k in ('implicit_vimp', 'implicit_cs2_vimp', 'implicit_cs2_vimp_npass')):
        drop.append(kv)
    elif (k == 'implicit_op_stencil' and per1) or (k in ('vet_gd_gl_nmu', 'vet_gd_gl_nphi') and k in rad):
        drop.append(kv)
    else:
        keep.append(kv)
hyd = blk.read(src, 'hydro')
cs = mesh.get('use_cubed_sphere', 'false') in ('true', '1')
# am_transport is defined on cubed-sphere / multi-D spherical-polar meshes only (hydro_am.cpp:172-177):
# on any other mesh the key is a geometry drop (logged), like implicit_marshak_face
mb = blk.read(src, 'meshblock', args if os.path.exists(args) else None)
spp = mesh.get('use_spherical_polar', 'false') in ('true', '1')
nx2 = int(mb.get('nx2', mesh.get('nx2', '1'))); nx3m = int(mesh.get('nx3', '1'))
amok = cs_ = mesh.get('use_cubed_sphere', 'false') in ('true', '1')
if spp and nx2 >= 2 and not (nx3m < 2 and hyd.get('sp_cart_polar_momentum', 'false') in ('true', '1')):
    amok = True
if hyd and os.environ.get('NOAM', '0') != '1' and not amok:
    drop.append('hydro/am_transport=true')
if hyd and os.environ.get('NOAM', '0') != '1' and amok:
    keep.append('hydro/am_transport=true')          # idort-1020: production <hydro> am_transport
if cs:
    keep.append('mesh/cs_seam_flux=shared')
subprocess.check_call(['python3', T + '/setkeys.py', src, out] + A + keep + extra)
# the PROD canonical key vet_source replaces a test's alias vet_col_source
if any(kv.startswith('rad_m1/vet_source=') for kv in keep):
    L = open(out).read().split('\n'); b = None; o = []
    for ln in L:
        st = ln.strip()
        if st.startswith('<') and st.endswith('>'): b = st[1:-1]
        if b == 'rad_m1' and st.split('=')[0].strip() == 'vet_col_source':
            drop.append('deck:vet_col_source=' + st.split('=', 1)[1].strip()); continue
        o.append(ln)
    open(out, 'w').write('\n'.join(o))
import atexit
atexit.register(lambda: open(os.path.join(os.path.dirname(os.path.abspath(out)), 'prod_drop.txt'), 'w').write(' '.join(drop) + '\n'))
None and open(os.path.join(os.path.dirname(os.path.abspath(out)), 'prod_drop.txt'), 'w').write(' '.join(drop) + '\n')
