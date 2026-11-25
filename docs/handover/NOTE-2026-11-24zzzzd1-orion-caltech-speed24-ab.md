# NOTE 2026-11-24zzzzd1 orion -> Caltech: H200 A/B of the speed branch (branch caltech-speed24-e7b2eea39)

**Order in your queue:** after sdorder5 (zzzzb5) and the cs pattern arms (zzzzc6), one free node; production first. **Trimmed by me** (no arms with little to gain): F on m35mms late, M35ZAMS and t9l62 quiet (not 00015); G on m35mms and M35ZAMS only; H on t9l62 quiet and 00015 only (the wall gives T4 1-2 % on the others). Total about 10 short runs. Nothing of this goes into production without the USER.
**Extra, read-only, no run:** in the production bsg (c1 and earlier bsg runs' logs) grep the MLATDIAG `nfb` column: SPEED24 found 6747 fall-backs in one old bsg arm (BB_bsgcutr); sub-group L_z is dropped in such groups. Report nfb per print of the c1 bsg log.

Text by SPEED24:
Branch vet-speed24 off 95813cf85, tip e7b2eea39: K1 1bd950aa8 (key), T4 e68e5052e (key), **A1+A2 e7b2eea39 (no key)**. Build: CUDA 12.9 + g++13, tiny BSG gate 3.686e-11. With the keys absent the solution is that of 95813cf85 (Orion CPU: 955 run files of the full suite bitwise; m35mms candidate cut, 8 blocks on 8 / 4 / 2 ranks, hst bitwise; not yet seen on a GPU). Idle nodes, copies only, no rst / bin kept, every candidate key named. Same node, back to back with the present 95813cf85 binary.

States: m35mms late (the CAND_1124 m35mms_R rst, c6950), M35ZAMS (m35zams_R rst, c22155), t9l62 quiet (stop rst, PQ), t9l62 00015 (PO). 100 cycles per arm, no outputs.

**F. A1+A2 (binary A/B, keys absent): highest value.** base / tip / base / tip; s/cyc over the last 50.
- hst (user + hydro) must be IDENTICAL between the two binaries; if not, stop and report the first cycle and column. Printed counters (`CLEAN FCONE`, floorcap, fwb) may differ in the last digits.
- Then one nsys (profc.sh, 12 cycles, NCU_ON=0) of the tip on t9l62 quiet and m35mms: us per launch of CleanStagFSolve#6 (was 10962 t9q / 1100 m35), CleanStagFCone#4 (5566 / 196), FloorGasCap#2 (17747 / 8244; kernel names may be renumbered: give the table), gaps > 1 ms per cycle (98 / 119 ms), MPI_Allreduce total and max in stats_r0.

**G. K1** `rad_vet/clean_stag_pcr_factor=true`, tip binary, off / on on all four states: hst identical; s/cyc; GPU memory (+ about 1.1 GB m35mms, 3.4 GB t9l62 per GPU). Ceiling about 4 % on m35mms, 2 % t9l62.

**H. T4** `hydro/dt_courant_sum=0.9` (dt = min(0.9/sum_d r_d, 0.5/max_d r_d); the 0.5 is fixed in the formula, STABAN's PLM guarantee), tip binary, off / on on all four states, dt_report=10. Stop an arm at the first CW / NaN / dt_retry. Report dt per 25 cycles, DTREP lines, passes and its, s/cyc, fofc / pos_gas / floor counts, L_top / L_in, KE, (t9l62 quiet) radial E, T, rho at equal model time. Judge speed as wall per model second.
- On Orion's full suite with the key on, three 1-D tests miss their pins: radwave err_cph 4.6e-3 (bound 3e-3), idort_radshock_m2 Tg +7 %, m5 rho +14 %: in 1-D the cap takes nu from 0.3 to 0.5 and the time error grows with it; in 3-D the sum stays <= 0.9, the worst case of today.
- USER rule cfl 0.3: T4 keeps the 0.3-per-direction worst case, it changes the step of anisotropic cells. Not for production without the USER.

Order: F, then G and H together. PRF_zam (front of the queue) helps read F on M35ZAMS.
-- orion (coordinator), drafted by SPEED24
