---
name: red-giant-nan-check-hides-origin
description: OPEN 2026-09-09 — T7_nodust_fix (grains off, lid-free, WB fix, clean restart) still died at 1.04e6 with NO precursor: every cell from i~141 up NaN at the 100-cycle NaN check, dt 30.65 to the end. Reading: dt min-reduction ignores NaN, so a NaN born near the top spreads ~2 cells/cycle unseen; nan_check_cycles=100 hides the origin. T8_nancheck (194637) = same run with the check every cycle + rst every 1e5
metadata:
  type: project
---

T6 (8.62e5, span 261..323) and T7 (1.0408e6, span 141..323, 2.44e6 bad cells = ALL cells
above i~141 incl. ghosts) died the same way: no dt collapse, hst smooth, everything NaN above
some radius at the next 100-cycle check. That is what an instantaneous NaN (table miss, log
or sqrt of a negative, 0/0) near the top looks like when dt ignores NaN: ~2 cells/cycle
radially and angularly, and the RT down-sweep carries a top NaN through the thin column in
one call. The span's lower edge is where the spread STOPPED, not where it started.

Recipe: T8_nancheck = T7 (binary ec68df12, restart from T6's t=5e5 rst) + time/nan_check_cycles=1
+ output3/dt=1e5, job 194637, expected to die at ~1.04e6 (CPU runs are bit-reproducible) and
name the first bad cell within one cycle. An Opus agent then reruns the last stretch from the
rst with 1e3-s dumps and an improved checker (active vs ghost counts, neighbour state,
block->panel map). Candidates: EOS-table window (log T 2..7.5) at a cooling top (T7 top
cells 1600-1900 K, draining), rho -> dfloor 1e-18, rt_apply deq with A/E<=0 or ei+de<=0
(rt_newton status in this input?), the rt_top_re column / top slot with p->0, the sponge,
fill_open_out's exp/log before the guard. prod11 (lidded, walled, older binary) never shows
it -> the open top / drained atmosphere is where to look.

**How to apply:** always run red_giant with time/nan_check_cycles=1 while debugging (cheap), and
rst every 1e5 so a death can be replayed. Related: [[red-giant-restart-radial-ke-injection]],
[[red-giant-dust-opacity-kills-lidfree]].
