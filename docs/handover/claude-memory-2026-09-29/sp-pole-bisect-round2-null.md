---
name: sp-pole-bisect-round2-null
description: sp polar blow-up bisection round 2 was a NULL experiment (wall_r 4990eb41 and x3shift_r 863e8337 bit-identical to dirty 18f5dd21); culprit is in 2a64e2c7..863e8337, round 3 arms 54af8380/a8cfb83e/de667f32/eb8e3cb1
metadata:
  type: project
---

Read 2026-09-09: the three .hst files (ideal_r, wall_r, x3shift_r) share one md5 and the
rot-4 dumps are payload-identical, so fefe6a17, b61c5d67, 581a1e43, 0f9959ee are inert on
this (ideal-MHD, no-resistivity) input. Since 863e8337 is dirty and 2a64e2c7 clean, the
culprit is one of a8cfb83e (x3-face states in the local basis), de667f32 (third-difference
polar EMF diss; nodiss arm is dirty from rot 1, so diss is a stabiliser not the culprit),
eb8e3cb1, 863e8337 (face-midpoint shift). Round 3 arms live in bench/sp_pole_bisect/r3_*;
gate: history 1-ME at rot 4.7 (clean ~8e31, dirty 1.6e34). Lesson: before a bisection arm
runs 6 h, diff a 200-cycle history against the neighbours to prove the commit is ACTIVE.
See [[sp-pole-bottom-radial-blowup]].
