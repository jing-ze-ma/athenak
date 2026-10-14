# NOTE (XB, Raven) -> DeltaAI: bsgli-val REDIRECT -- build fin4-linhao-bc-8e018272, NOT eb1e032c

TASK-2026-10-14-deltaai-bsgli-val.md (70a37593) stays as written EXCEPT the source:
- branch **fin4-linhao-bc-8e018272** (commit 8e018272) instead of fin4-linhao-bc-eb1e032c.
  eb1e032c ABORTS at program exit when he_bc_linhao is on (Kokkos View freed after finalize; fixed in 37798491);
  8e018272 = 37798491 + the history column L_in = the step-effective wall flux (diagnostic only; physics identical).
- Everything else as the TASK (LI deck + he_bc_linhao = true IN THE DECK -- a command-line key that is not in the deck is fatal;
  C60 vs L30+30 bitwise; Coff gate; md5, s/cycle, memory).
- Optional extra if time allows (<= 15 min): the same deck 20 cycles with output1/dt=1e-6 and report per cycle
  user.hst cols 19 (Etot), 7 (L_in), 6 (L_top) for the energy budget dEtot vs (L_in - L_top) dt.
Results to docs/handover/deltaai_bsgli_val/ as in the TASK.
