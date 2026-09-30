# TASK 2026-09-30 (from Caltech, user 09-29): 3x WASP-121b ck solver A/B on viper

## Why
Caltech 3x production 3640302 (793e03c3, keys = viper's 3x input as shipped) runs clean, but its dt fell
12 s -> 7.5 s between rot 18 and 30 and is still falling (10x settled at ~9.5 s by rot 5 and stayed flat
to rot 300; 1x ~14 s). The limit is the RADIAL CFL in cells 60-65/76: supersonic night-side downflows
(v1 ~ -20 km/s, Mach ~2.4) in 64/6144 columns reaching p ~ 1e-5 bar (10x: 1 column, only at 7e-12 bar;
1x: 11 columns, p < 1.2e-7 bar). T-floor (200 K) cells 1.1 % in 3x, 39 % on the DAY side (1x 0.2 %/6 %,
10x 1.0 %/9 %). Night-side collapse + infall itself follows metallicity (10x strongest) and is not the issue.

The only non-physics input difference 10x vs 3x is the ck solver:

| key | 1x, 3x | 10x |
|---|---|---|
| ck_impl_dtmax | 0.5 | 0.25 |
| ck_impl_maxit | 16 | 24 |
| ck_impl_rsec | off | 20 |
| ck_impl_tol | 1e-8 | 1e-7 |

cknewton_0928 (two_stream_column_ck.hpp ~l.189): with dtmax 0.5 the 10x top cells (p ~1e-10 bar) sit in a
period-2/4 limit cycle (T 1700 <-> 5300 K); dtmax 0.25 + maxit 24 removes it. Hypothesis: 3x hits the same
cycle. Production logs cannot show it (ck_impl_verbose off).

## Caltech run (in progress)
/resnick/scratch/jingze/w121_3x_ckab_0929: restart from 3x prod dhj.00072.rst (rot 36, cycle ~413.9k),
nlim 417000, ck_impl_verbose on; A = 3649979 (3x keys), B = 3649980 (10x keys). Results -> NOTE on Caltech.

## Ask for viper
Repeat the A/B on viper (HIP) once a 3x state at rot >~ 30 exists there (the user copies Caltech
dhj.00072.rst, 141 MB, if needed): same two arms, ~3000 cycles, ck_impl_verbose on. Report
NOT-CONVERGED count, passes per call, top-cell T oscillation, dt. Do NOT change any production keys;
the user decides whether the 3x chain link switches solver keys.

## Snapshot for viper (staged 09-29 19:41 PDT)
Caltech `/resnick/groups/carnegie_poc/jingze/to_viper_0930/` (groups space, not purged):
- `dhj.00072.rst` 141347192 B, md5 `4a73aa7c1de70572e3af0337baca149f` (MD5.txt): 3x prod, rot 36, t 3.9656e6 s, cycle ~413.9k
- `run.athinput` (3x prod input), `COMMIT_binary.txt` (793e03c3), `ab.sub` (Caltech A/B script: nlim 417000, verbose on; arm B adds
  `problem/ck_impl_tol=1.0e-7 problem/ck_impl_maxit=24 problem/ck_impl_rsec=20.0 problem/ck_impl_dtmax=0.25`)
Fetch (user or viper, whichever has the Caltech login): `scp -r jingze@login.hpc.caltech.edu:/resnick/groups/carnegie_poc/jingze/to_viper_0930 .`
then `md5sum -c MD5.txt`. Viper build must be >= 793e03c3 (restart file layout identical; bitwise vs CUDA not expected).

## Viper reply (09-30 ~04:30, user)
User: wait for Caltech's own A/B (3649979/80) first; viper repeats it only if the Caltech result is ambiguous (then the user copies to_viper_0930 to viper). No production key change without the user.
