"""make the gate / timing inputs with every switch spelled out (command-line overrides
of keys that are not in the file are rejected)"""
import re

NEW = ("vet_mb_agg            = true      # sc-scale: one message per neighbour rank\n"
       "vet_mb_mom_fuse       = false     # sc-scale: moments summed in the ray launch\n"
       "vet_mb_agroup         = 1         # sc-scale: hybrid, ranks per angle group\n")
MB = ("vet_mb_halo           = 1         # SC scaling: layers per band exchange\n"
      "vet_mb_overlap        = false     # SC scaling: interior of the next layer under MPI\n"
      "vet_mb_kernel         = ray       # SC scaling: ray | cell\n"
      "vet_mb_angles         = false     # SC scaling: PROTOTYPE angle decomposition\n"
      "vet_mb_mom_batch      = 4         # SC scaling: launches per moment launch\n")


def add_after_nphi(s, extra):
    i = s.index('\nvet_nphi')
    j = s.index('\n', i + 1)
    return s[:j + 1] + extra + s[j + 1:]


C = '/viper/ptmp2/jinma/scscale_0923/cpu/'
s = open(C + 'he_slab_mb.athinput').read()
open(C + 'he_slab_mb_o3.athinput', 'w').write(add_after_nphi(s, MB + NEW))
s = open(C + 'he3d.athinput').read()
open(C + 'he3d_o3.athinput', 'w').write(add_after_nphi(s, MB + NEW))

G = '/viper/ptmp2/jinma/scscale_0923/gpu/'
s = open(G + 'he3d.athinput').read()
s = add_after_nphi(s, NEW)
s = re.sub(r'\nvet_mb_halo *= *1 ', '\nvet_mb_halo           = 3 ', s)
s = re.sub(r'\nclosure *= *eddington', '\nclosure          = vet_sc', s)
FAST = ("implicit_halo_direct  = true      # FAST (README_FAST)\n"
        "implicit_od_cache     = true\n"
        "implicit_krylov_fuse  = 3\n"
        "implicit_op_stencil   = true\n"
        "implicit_precond      = rbgs_fwd\n")
s = add_after_nphi(s, FAST)
open(G + 'he3d_fast.athinput', 'w').write(s)
