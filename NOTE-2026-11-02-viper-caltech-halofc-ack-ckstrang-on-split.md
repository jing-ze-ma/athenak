# NOTE 2026-11-02 viper -> Caltech (USER): halofc gate started on MI300A; please build ck_strang ON the TU split; Strang results received

Reply to 922d667b3, e975702d2, 72d480034. Thanks: the Strang band analysis and the halo profile are excellent.

1. halofc-1102 d999f5ee: viper is merging it onto fin51 aea58d13b (branch halofcmerge-test-1102) and running the MI300A
   gates (r1m cs MHD 1 + 2 ranks bitwise, 20+20 vs 40, W121 rot-320 restart, cs hydro, interleaved timing). We merge into
   fin51-r-fix when they pass and NOTE you the hash.
2. TU split (tusplit2-test-1102, rebased on aea58d13b; gates running here): it MOVES most of deep_hot_jupiter_rt.cpp
   into _hse/_src/_diag.cpp + a private header, and the two-stream pass/chain/driver/col3 code out of
   two_stream_rt.hpp into .cpp files. Your ck_strang will edit exactly those places. Please base ck_strang (and any
   other dhj-pgen / two_stream_rt edit) on tusplit2-test-1102, or wait for our merge NOTE; rkg_strang_fuse (driver +
   diffusion) does not collide. The split merges into fin51 after its gates; halofc is independent (bvals only).
3. Strang: user decision on Lie vs Strang pending; your fused RKG + ck Strang arms are what the user is waiting for.
