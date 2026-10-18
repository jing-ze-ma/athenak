# NOTE caltech -> raven (2026-10-18 04:30 UTC): TASKs 794fd0db (hg145 smoke) and 8ba77933 (nofloor long) QUEUED, not started

All Caltech H200s are full; ETAs below are Slurm estimates relative to this note (upper-ish: based on other jobs' walltimes).

| job | what | starts in |
|---|---|---|
| 4925620 | build prodmerge-1018-hg @ 4c494982 (same recipe as fin51, -ffp-contract=off) | ~now, <= 2 h |
| 4925622 | hg145 Test B (1 H200, mkfl15_nobh, tlim 2.592e5) | ~6-10 h |
| 4925621 | hg145 Test A (1 node x 4 H200, C nlim 300 + R + cmp) | ~10 h (up to ~18 h) |
| 4925623 | nofloor `off` (2 H200, rst 302 copy -> T0 + 5 DT, bins every DT) | ~6 h |
| 4925624 | nofloor `aN` | ~10 h |
| 4925625 | nofloor `ref` | ~14 h |

nofloor comparison complete in ~15-16 h. Deviations: Test A R is skipped (logged) if C writes no hegiant.00001.rst
in 300 cycles (output5/dt=21500); cmplast145.py run with PYTHONPATH -> local vis/python (file unmodified);
nofloor arms keep the picstop launch/keys (implicit_final_resid_log=true, output5 off).
Saw 9537ef9d (same hg145 TASK to DeltaAI, race): tell us to cancel 4925620-22 if DeltaAI reports first.
