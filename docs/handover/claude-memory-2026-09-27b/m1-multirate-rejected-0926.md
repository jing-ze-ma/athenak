---
name: m1-multirate-rejected-0926
description: M1 multi-rate (implicit_mr_every 2/4) REJECTED 09-26 even for relaxation - on top of hesdirk2 cfl 0.9 it gives no gain or crashes, and damages the box state permanently
metadata:
  type: project
---
Tested 09-26 (/viper/ptmp2/jinma/mr_0926/RESUME.md). The code accepts mr only with hesdirk2 and even k, so all arms ran hesdirk2 cfl 0.9.
- He box: k2 costs the same as the reference; k4 is -26 % but crashed at t 1623 s (top cells, dt collapse). Both inject KE 1000-4000x the reference.
- Switch-back (k4 to 200 s, then off): still 6-11 % off in plane-mean E at 2000 s, KE 300x. It does NOT relax back.
- Static wedge: k2/k4 are 1.6-1.7x SLOWER (iterations 3.6 -> 39/81). It relaxes back after the switch.
- The old README gains (box -37 %, wedge -46 %) were measured on be and do not hold on hesdirk2 cfl 0.9.
- Caveat: the mr box reference was quiescent, not convecting.
**Verdict:** the failure is the multi-rate sub-stepping, NOT cfl 0.9 (the mr-off reference at hesdirk2 cfl 0.9 ran clean). Do not use mr, not even in relaxation ([[two-phase-relax-then-accurate]]).
Setup pitfalls:
- viper SBATCH_EXPORT=NONE drops env vars;
- a command-line key must exist in the input file;
- the tolsweep box input has two restart blocks (FATAL).
