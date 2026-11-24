# mkgate.py REV: the Caltech tiny BSG gate deck in rc1 form: H200CHK/bsg_wngl_smoke_dtfac1e6.athinput in the release names, the
#   COMPAT "no effect" keys stripped, + vet_gd_wall_sym / vet_gd_theta_seam_sym (theta-periodic: rc1 FATALs without the latter).
#   Run as before: -i IA123_gate_bsg_smoke.athinput mesh/nx2=32 mesh/nx3=32 meshblock/nx2=32 meshblock/nx3=32 time/nlim=2
#   Reference (Caltech, rc1 da92efe70 with the seam key): clean_init pass 1 3.686e-11, as IA11 / IA12.
import sys, re, hashlib
sys.argv = [sys.argv[0], sys.argv[1]]
src = open('mk.py').read(); ns = {}
exec(src.split('for v in')[0], ns)          # RV / HY / GONE / rename / add from mk.py
REV = sys.argv[1]; F = '/orion/ptmp/jinma/clean_1115/stag2/coord/wt_note/docs/handover/caltech_results/H200CHK/bsg_wngl_smoke_dtfac1e6.athinput'
L = open(F).read().splitlines()
O = [f'# IA12.3-rc1 tiny BSG gate deck ({REV}) = H200CHK/bsg_wngl_smoke_dtfac1e6.athinput, release names, COMPAT no-effect keys stripped, + both seam keys (mkgate.py)']
O += ns['add'](ns['rename'](L), 'gate', rv=ns['RV'][:2], hy=[])
open('IA123_gate_bsg_smoke.athinput', 'w').write('\n'.join(O) + '\n')
print(hashlib.md5(open('IA123_gate_bsg_smoke.athinput', 'rb').read()).hexdigest())
