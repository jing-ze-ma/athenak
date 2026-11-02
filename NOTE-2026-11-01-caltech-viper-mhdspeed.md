# NOTE Caltech -> viper, 11-01 18:10 PST: mhdspeed-1101 merged, then REVERTED on fin51-r-fix (speed regression)

- fin51-r-fix: 0abb64b6 -> 0892e7e9 (merge mhdspeed-1101 4f152fc7) -> aeafae97 (revert). Tree at aeafae97 = 0abb64b6. Please pull; do not build 0892e7e9.
- Bitwise: everything IDENTICAL on H200. ACB vs base: r1m from off_0p50.rst 150 + 300 cycles (rst + mhd.hst), hydro 150, continuation 150+150 = 300.
  Merge 0892e7e9 vs 0abb64b6 binary: 150 + 300 cycles IDENTICAL. CPU 3-rank C8 fresh 30 cycles IDENTICAL, but it ran RKG stage 0, so B is not exercised there.
- Timing (1 H200, developed state, 300 cycles, interleaved x2, on a clean GPU: no other compute app during the window):
  base 81.4, A 121.9, AC 121.4, ACB 124.4 ms/cycle (median; the two replicates agree within 0.1).
  kp profile: SendBuff 1.16 ms/call (base) -> 2.65 ms/call (A). The chunked league redoes the large per-team seam/vertex
  preamble nchunk times, and most of the chunks are empty. C gives nothing measurable. B costs +2.5 % on top of AC.
- Branch mhdspeed-1101 (4f152fc7) is kept for rework. Nothing to do on viper.
