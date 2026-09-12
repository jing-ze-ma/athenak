---
name: red-giant-seam-floor-eos-garbage
description: 2026-09-12 red-giant GPU death MECHANISM (bench/RG_fofc_long/deathmap/analysis): a convective downflow drains a photospheric column 250x; ONLY on a panel-SEAM row its internal energy then collapses (tau 180 s) to the code floor e/rho = 1.4226e7 erg/g (T 0.07 K, far below the table's 100 K), the across-seam GHOST column hit that floor 2.9e4 s earlier; the next EOS inversion returns a SATURATED GARBAGE state T=2.63e18 K cs=2.66e9 |v|~1e12 in every collapse cell -> dt collapse. FOFC never fires there (it works only at r>3.9e12). 14 corner-ghost cells at 1e7-9e8 K are the only other anomaly
metadata:
  type: project
---

Reproduced deterministically (restart 1.7e6 -> collapse at cycle 57913 = the original). Four collapse cells,
all at r 3.32-3.44e12 (photosphere), 3 of 4 at seam distance 0, one at cube-vertex distance 1; all report the
SAME state T=2.629264e18 K, cs=2.660693e9, |v|~1e12 with only rho differing => a c2p/EOS-inversion artefact,
not heating. Death column time series (gid 15 (9,18,339)): rho 2.5e-10 -> 1.0e-12 over 7e4 s with v_r inflow
-1e6 -> -6e6 (an ordinary downflow lane: 12 equally drained columns exist, seam enrichment neutral), e/rho
flat then -> 1.4226e7 EXACTLY (a floor value) with tau -180 s; the 11 non-seam drained columns stay at
8-9.5e3 K. Every cell with T_proxy < 100 K in every dump is on a seam row (enrichment 8.26 = max).
No dump holds a hot active cell (the runaway is < 780 s); precursor = the seam ghost column floored first.
**Two defects to fix:** (1) the general-EOS inversion must CLAMP to the table's T range and never return
2.6e18 K (cf. [[eos-inversion-nan-trap]], [[eos-table-dump]] table NaN below 71 K); (2) the energy floor that
yields e/rho = 1.4226e7 is inconsistent with tfloor (100 K) at rho ~1e-12 -- floor e_int at rho*e(T_floor),
and the seam resample / ghost fill must not hand the EOS sub-table states. FOFC is irrelevant here.
Corner ghosts at 1e7-9e8 K (14 cells, r 3.4-4.0e12): the sign-preserving limiter allows upward overshoot;
only the conduction cross term reads them ([[cs-wire-fill-wip]]).

**Confirmed clamp-independent (RG_fofc_ctl2, UNCLAMPED binary, 11649849):** dies at cycle 57906 t=1.77264e6 in the SAME
column gid 15 (k,j)=(9,18), i=340, same saturated state T=2.62926e18 cs=2.66069e9 (rho 1.5e-10 this time). The
seam-row runaway + EOS sentinel is robust: 3 of 3 realizations of this trajectory die within 10 cycles of each other.

**FIXED 5fe7be1c (pushed):** e/rho = 1.4226e7 was EnergyFromPressure(rho, pfloor=1e-20) pinned on the root find's
LOWER bracket 10^(eos_logt_min-3) = 0.1 K (pfloor unreachable, no tfloor set); T=2.629264e18 (code units) =
10^(ymax+3) = 3.16e10 K, the UPPER bracket edge returned as a temperature for an out-of-range or non-finite target
(EOSTable::SolveLog, eos_table.hpp ~334-360). Fix: both e->T and p->T inversions clamp to [ymin,ymax], non-finite
targets intercepted, new EventCounters::neos_tclamp + `eos_tclamp` event-log column, efloor gate >= e(rho,T_min).
Bitwise on cs_mhd_prod3 (20 cycles), ck test passes, HIP builds. Residual: a clamped cell has p=p(rho,T_min)
independent of e (e not raised); eos_tclamp makes it visible. Test: bench/RG_fofc_long/eosfix (11652436,
restart 1.7e6 -> 2.0e6 through the 1.7728e6 death).
**VERIFIED 2026-09-12 ~23:30:** eosfix run (11652436, binary 5fe7be1c) went 1.7e6 -> 2.0e6 with dt flat at 30.67 s
through 1.7728e6 (3/3 unfixed realizations died there); eos_tclamp 3e4-3e5 per interval during the episodes.
RG_fofc_long2 (bench/RG_fofc_long2) relaunched from rst 1.7e6 to 3e6 on that binary.
