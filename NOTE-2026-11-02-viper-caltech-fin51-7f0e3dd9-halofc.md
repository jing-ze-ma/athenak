# NOTE 2026-11-02 viper -> Caltech (USER): fin51-r-fix = 7f0e3dd90: halofc-1102 merged (MI300A bitwise); -3.9 % on W121 MI300A

fin51-r-fix 3957c9562 (TU split) -> 7f0e3dd90 = merge --no-ff of halofc-1102 d999f5eea (clean; bvals.hpp, bvals_fc.cpp only).
MI300A gates vs 3957c9562 (viper 12574584; 2 GPUs, fresh cs MHD starts have cs_seam_fc_nolimit = true by default):
- r1m cs MHD 40 cycles, RKG on 2 ranks / RKG off 2 ranks / RKG on 1 rank: rst whole file + mhd.hst IDENTICAL
- 20 + 20 vs 40: rst data + last 20 mhd.hst rows IDENTICAL
- W121 dhj.00641 (rot 320) +20: rst whole file + mhd.hst IDENTICAL
- cs hydro r1m_viper_hydsc 20 cycles: identical (does not call PackAndSendFC; needs am_hst_split = false, else 25 > 24
  user hst columns FATAL: deck issue)
Timing MI300A, W121 restart +150 cycles, 2 ranks, interleaved x2: 246.8 -> 236.9 ms/cycle = -3.9 % (the pack is a smaller
share of this state than of your r1m H200 case, -18.6 %).
Not covered by any gate (yours or ours): SMR on cs (coarse/fine buffer ranges, coarse seam tables), polar BC + cs, > 2
GPU ranks, blocks < 3 cells across (fallback stencils), CS_SHEAR. The segment table is built once and never
invalidated: correct only while cs refuses AMR / runtime rebalancing (as the existing seam table already assumes).
Please build 7f0e3dd90 (or rebase rkg_strang_fuse / ck_strang on it).
