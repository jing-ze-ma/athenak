# TASK 2026-11-10 viper -> ORION: raycheap18b -- add the cmem arms (design (a), production candidate) to raycheap18

Follows TASK raycheap18 (e747c5c31). Branch fork/raycheap18-1110 now bb58d80c2 (on dea3e3732): d981e435c
implicit_blend_dtc (diagnostic) + 91fdbe423/9cc04e09d/bb58d80c2 **implicit_hr_cmem** (your 17b design (a); user: blend
keeps only angular info, no dt / opacity in the weight).

## implicit_hr_cmem = true (read only when named; absent = bitwise)
f0x1..3 hold F_c, evolved by the plain face equation F_c' = th (F_c^n - c^2 dt dP/dx ...), so every existing copy (f0x*n,
t2f* saves/restores/avgf, kk slopes / t2i increments, restart) carries F_c with no extra code. The applied flux
(1 - w) F_c' + w F_hr(E') -- exactly what the rows always solved -- is stored in fax1..3 and read by iw F*, the write-back
cell F and the transverse divergence (fp/fm/gp/gm). w unchanged. Inert (echoed) without the half-range beam. FATAL with
implicit_vimp / Linhao bath; with the conservative gas source / rad inertia a WARNING (they take F_c as the old face flux:
exact without gas momentum exchange, off by (1-w)(F_c^n - F_a^n) with feedback -> TEST key there).

## viper measured (CPU login, athena_none_cpu_bb58d80c_rc18; work dir /viper/ptmp2/jinma/rc18_1110/w)
Gates: keys absent bitwise vs dea3e3732 (SP pulse_u T_1/2/4, 17-digit tab); S09+gd with cmem on bitwise (T_1/2); restart
with cmem on bitwise (SP pulse_u T_4: run to t = 1 vs restart at t = 0.5, final tab identical).
| pulse_u arm, E | p (T, l 1..32) |
|---|---|
| SP | -0.25 -0.08 -0.07 0.29 |
| SP + cmem | 0.49 0.75 0.84 1.16 |
| **SP + cmem + x2_fields 13** | **0.78 1.27 1.49 1.91** (F1 0.71 1.13 1.32 1.64) |
| SP + dtc 1e-6 + x2_fields 13 | 0.86 1.27 1.55 2.06 (F1 1.99) |
| SP + cmem + x2 13, C | 0.50 0.15 0.29 (space defect unchanged) |
Space bisection on pulse_u C with dtc + x2 13 + ncore 64 (none restores order): base 0.05 0.09 0.29; thin_corr none
0.07 0.16 0.29; hr_face cell = base; thin_switch none + corr none + beam factor off 0.83 1.19 0.39; rotation off
-0.01 -0.72 1.45; thin_recon dc = base (hcface / js_face / bq not run). C and X agree (X 0.00 0.08 0.30): purely spatial.

## Runs (add to raycheap18)
CM = `implicit_hr_cmem=true time2_vet_x2_fields=13`, CMN = CM + `vet_col_ncore=64`.
(a') pulse_u SP-CM T to l = 64 (does F1 reach 1.8?), and C, X; (c') Part 1b with SP-CMN next to SP-FIXN;
(e') cs Part 2 with CSP-CMN next to CSP-FIXN (cmem also acts on the cs lateral faces).
(f) SPACE: continue (d) of raycheap18 with the untested arms implicit_thin_hcface=mean, implicit_thin_js_face=mean,
implicit_thin_bq=cell, and one diagnostic: the vgd_hr ratios (h+-_r) restricted between n and 2n (do the formal-solution
half-range ratios themselves converge in space?).
(g) Beams: viper ran the blend_1022 IDORT decks (Cartesian vet_sc: cb_FIN3_r4, ws66, shd_FIN3_r12g2) with SP / dtc / cmem
on bb58d80c CPU, but SP itself scores far from the 10-22 tables (cb E_L1 17.5 vs 0.57, ws66 10.0 vs 1.16, shd 97.6 vs
0.013) -> the harness or the deck keys do not carry over to the dea3e3732 base; viper is checking. Not a verdict yet.
Answer as NOTE raycheap18b on bsg-files-1009.
