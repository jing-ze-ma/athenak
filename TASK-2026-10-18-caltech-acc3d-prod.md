# TASK raven -> Caltech (2026-10-18): 3-D ACCRETOR PRODUCTION copy (race vs Raven), USER GO given

USER GO for the 3-D accretor production ("go. if it queues too slow you can try other machines"). Raven chain is queued
(31274287 + 31274288, est. start ~23 h). Caltech runs the SAME production as a race: **first site to START wins**. Post
the start time in your NOTE as soon as link 1 starts; Raven cancels its chain only after that (never cancel anything
on Raven yourself).

## 1. Build
Branch **accretor-3d-1017 @ 186306cb** (pushed; code 5bc6c3b4). Same cmake recipe as your earlier H200 builds:
`-D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D Kokkos_ENABLE_CUDA_LAMBDA=ON -D Athena_ENABLE_MPI=ON
-D CMAKE_BUILD_TYPE=Release -D PROBLEM=ry_per_accretor -D CMAKE_CXX_FLAGS=-ffp-contract=off` (nofma, as the Raven A100
binary athena_a100_ry3d_v7_nofma md5 210e2fee). Report md5.

## 2. Files: docs/handover/acc3d-prod-1018/ (md5, first 8)
- `acc3d_prod_mk.athinput.in` e9dab743 -> `sed "s#@FILES@#<files dir>#g"`; change no key (3 path keys: env_ic_file,
  env_opac_table, env_planck_table). = Raven prod/acc3d_prod_mk.athinput (6dd6ec53) except paths + one comment.
- `rosseland_mesakap_gs98_x0.70_z0.02.txt` d81b7a46 (MESA r15140 kap spliced into ext2; USER rule: MESA kap),
  `planck_ext2_gs98_x0.70_z0.02.txt` 99d668f6, `column_mesakap_z0.02_p381.txt` 5b813886 (IC column).
  No other input files (general EOS is computed in-code). All files < 5 MB (in git).
- `acc3d_chain.sh` 1ad77aed: the chain logic (run inside each link's allocation; env BN DK PD WT LAUNCH, see header).

## 3. Layout
1 node x 4 H200, **4 ranks, 2 MeshBlocks per GPU** (deck: mesh 420x58x512, MeshBlock 420x58x64 = 8 MBs; do not change),
your CUDA_VISIBLE_DEVICES wrapper, --cpu-bind=cores, UCX_TLS=^xpmem. (12.5 M cells; ran on 2 MI300A and 4 A100 40 GB, so 4 H200 fit).

## 4. Pre-check (short, one job; must pass before the chain)
With the substituted deck DK and `BK="output1/dt=1e-12 output2/dt=1e30 output4/dt=1e30 output3/dt=1e30"`:
- A: `-i $DK time/nlim=10 $BK` ; B: `-r <A>/rst/ryper.00000.rst -i $DK time/nlim=20 $BK` (the rst A writes at its end);
- C: `-i $DK time/nlim=20 $BK` (fresh 20 cycles).
PASS = B vs C final `bin/ryper.hydro_w.*.bin` and `bin/ryper.m1.*.bin` identical except the header line `file_number`
(cmp differs in 1 byte only), hydro.hst rows identical, user.hst columns 1-6,16-20,23 identical. user.hst cumulative
counters (cols 7-15, 21, 22, 24-26) restart from 0 per run (A + B = C; known, by design). Raven/Viper reference:
Viper 12265488 passed exactly this way. Report rc, s/cycle, dt (~1.42e-5), retries (expect 0), GPU peak memory.

## 5. Production chain (2 links for now, the partition's LONGEST walltime each, afterok dependency)
Each link: `BN=<binary> DK=<deck> PD=<run dir> WT=<walltime s - 1200> LAUNCH="<srun ... wrapper>" bash acc3d_chain.sh`;
the link's job script must exit with the chain's rc (exit 1 = STOPPED -> the dependent link must not run; use afterok +
kill-on-invalid-dep). Logic (do not change): static phase `problem/stream=false time/tlim=0.05` (static-start breathing
decays; ~15-30 min) in PD/static, then the deck as is (stream on, **tlim 1e9**, dt_min 1e-8) in PD/stream from the newest
rst until the wall limit; later links restart from the newest stream rst.
STOP rule (in the script): only the end-of-run summary `NON-CONVERGED=[1-9]`, FATAL, "Picard solve DIVERGED", NaN in the
log tail, rc != 0. NEVER per-pass "NON-CONVERGED after" or DT_RETRY lines.
Outputs (deck): hst 2e-4, hydro_w every 0.005 (250 MB), m1 every 0.01 (200 MB), rst every 0.01 code time (**2.9 GB each**);
the script prunes IN-JOB to the 3 newest stream rst (+ the final static rst). Disk ~40 GB per P_orb (0.458) + ~12 GB rst.
Cost/speed (Raven, unmeasured on H200): 4 A100 static 0.39 s/cycle; stream dt ~7e-6 before / ~1.6e-6 after contact (2-D).

## 6. Report: NOTE-2026-10-18-caltech-acc3d-prod.md on bsg-files-1009
Binary md5; pre-check PASS/FAIL (+ numbers above); chain job ids + ETA; **START time of link 1 immediately** (race);
then per link: rc, t reached, s/cycle, dt, retries, summary NON-CONVERGED, chain.txt tail. rsync-able run path.
