# TASK raven -> Caltech: open-ot NPALL, the implicit_cs2_vimp_npass = 3 PROPOSAL on every ACTIVE IDORT run (214 dirs)

From Raven worker OT, via the coordinator. This supersedes the FULL step of CVALL, which tested the rejected key-off variant;
the BIT step of CVALL stays useful. Not for production; the user decides.

## Binary
- 281dd686 + `git am patch/000{1..4}-*.patch` (no kokkos change in these). Plain CPU, same flags. Report the md5.

## Run
`bash run_npall.sh <binary> 8`.
- list.txt: 214 dirs, i.e. the runs where the conservative2 vimp prediction is active (DIFF in the CVALL BIT step):
  11 Cartesian, 134 sp (incl. patm/patmT two-stage), 69 ladders.
- Order of priority if time is short:
  1. sp/*_T and lad/*;
  2. the other sp;
  3. cart.
- Long runs (cloud, pipe, ddf, ddiff) may be skipped; list them.

## Report: docs/handover/NOTE-2026-10-13-caltech-open-ot-npall.md
- md5, rc counts, rc != 0 lines;
- `tar czf npall_runs.tgz --exclude='*.rst' runs`, split if > 90 MB, under docs/handover/caltech_npall/.
