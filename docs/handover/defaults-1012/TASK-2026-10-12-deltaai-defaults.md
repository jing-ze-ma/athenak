# TASK raven -> DeltaAI: method-1012 fresh-run defaults, full AG Car A, 60 cycles (1 job, <= 45 min)

From Raven worker D2 (defaults-1012). The user signed off the RT method; its keys are now the FRESH-run defaults.
This checks on GH200 that a deck naming NO method key gives bitwise the same run as a deck naming them all.
It races Raven gpudev job 31104217 (the same check on A100); first result wins.

## Limits
- 1 node x 4 GH200, 1 rank per GPU, `export MPICH_GPU_SUPPORT_ENABLED=1`, <= 45 min wall.
- Reply SKIP if the queue wait is over about 1 h.

## 1. Build
- `build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1).
- Source: fork branch `method-1012` @ **22f7021d** (the coordinator pushes it). Report the md5.

## 2. Inputs (this directory; G = /work/nvme/bivj/jma20/defaults_1012)
- `agcA_GU.athinput.in`: the current production deck (prodfinal run/agcA/agcar_A_fresh.athinput) with every
  method key AND `vet_gd_rotate_every` removed. Resolve `@FILES@` = docs/handover/agcarA_files -> `$G/agcA_GU.athinput`.
- `agcA_GN.athinput.in`: the same deck with every method key named at the method value. Resolve `@FILES@` the
  same way -> `$G/agcA_GN.athinput.tmpl` (keep `@RV@`; job.sub fills it from GU's log line
  `vet_gd_rotate_dt = auto: 20 x the first dt ... = <RV> (recorded ...)`).
- `cmpout.py`: bitwise compare of two run dirs (bin dump data via bin_convert, hst rows). Copy to `$G/cmpout.py`;
  set `ATHENA_VIS=<your athenak checkout>/vis/python` in the environment (needs numpy).
- `gpuwrap.sh`: copy the one from your consfix_1011 / bsgrec_1011 directory.

## 3. Runs (job.sub; `DRY=1 BP=x bash job.sub` prints the commands)
| run | deck | cycles |
|---|---|---|
| GU | no method key named | 60 |
| GN | every key named, vet_gd_rotate_dt = GU's auto value | 60 |
Both: time/ndiag=1, m1 + hydro_w bin dumps every 2.5e4 s, no rst.

## 4. Report: NOTE-2026-10-12-deltaai-defaults.md
- Per run: rc, FATAL, NON-CONVERGED, the Picard line, s/cycle (steady cycles, excluding cycles 0-2), wall.
- GU's `vet_gd_rotate_dt = auto` line.
- The cmpout.py line (expect `... 0 differ -> BITWISE`); if not bitwise, the first differing file.
- Verdict line: fresh defaults == all keys named on GH200: yes / no.
