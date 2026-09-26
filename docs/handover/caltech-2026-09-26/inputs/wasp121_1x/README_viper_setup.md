# wasp121_0925/x1: WASP-121b at 1x solar metallicity (sponge tests)

This is the same setup as `../README.md` (Sing+2024 planet, star and orbit; P = 1.27492504 d;
Bond albedo 0.277; WASP-121 band-flux SED; route B, nq2, sponges; floorbound + kkt_demax).
The 10x ck Newton does not converge on the dayside, so this variant uses 1x solar
metallicity. Only the metallicity-dependent pieces change.

## What changes

- **Metallicity:** `met` and `rad_met` are 0, and `eos_xh` / `eos_yhe` are 0.7381 / 0.2485 (solar,
  as in the old dhj inputs).
- **Tables:** the repo's 1x hiT2 tables, `ck/Premixed_1x_g8_11_hiT2.txt` and
  `CE_tables/FastChem_ck_1x_int_hiT2.txt`, read through `ck_data_dir` = the repo's
  `data/exo_fms_ck`. They equal the base 1x table below 6100 K, the only range the RCE uses.
- **IC:** `ic_profile = ic_w121_1x.txt`, and `ic_profile_mu` is emptied, so there is no
  day/night IC.
  - The profile is the 1x global-mean ck RCE (`rce/`; A_B 0.277; iterated to g 1168.07,
    ap 1.12648e10).
  - It has 269 uniform rows from 9.5e-11 to 2240 bar: isothermal at the top, and the exact
    solar-EOS adiabat below 283 bar. The RCE bottom slope is 0.1471 against nabla_ad 0.1462.

## p_ref

p_ref (`transit_1x.txt`, 1x opacities, NRS1 2.7-3.7 um) is **7.07e-4 bar**.

- The tau_chord = 0.56 level gives 5.7e-4 bar.
- A limb shift of -300 / +300 K gives 6.5e-4 / 1.04e-3 bar.
- TESS gives 7.5e-4 bar.

## planet_setup.py output (`setup_w121_1x.txt`)

- **Geometry:**
  - ap = x1min = 1.126575e10 cm; grav = 1167.875;
  - x1max = 1.631037e10 cm, where p = 1.41e-9 bar on the IC column;
  - p_ref isobar at the pole 1.232081e10, at the equator 1.259071e10, area-equivalent R_p
    1.245391e10.
- **SPARC grid:** nx1 = 76, min 3.28-3.47 cells/H below 1e-6 bar, dt estimate 13.7 s.
- **Production grid:** nx1 = 256, 9.6-10 cells/H (1.92 x plan A), dt estimate 4.8 s.
- **Files:** `grid_w121_1x.env` (G_SPARC, G_PROD, G); `sparc_w121_1x.athinput`,
  `prod_w121_1x.athinput`; `keys_1x_vs_10x.txt` (every key that differs from the 10x input).

## CPU gates (`gate/`, athena_cpu_f637b3be = the source of sparc_0925/athena.gpu e3a8442e)

- Both inputs parse.
- Fresh start (76 x 16 x 16 per panel, ghosts written, T4 + c2 + every 4, 16 cycles):
  - **no floor event in any cycle (0-15)**;
  - ck_implicit: 4 full calls at 3 passes each, 0 NOT-CONVERGED (res 7e-8 at the first
    call, 8e-9 after).
- Deepest cells on the table adiabat to < 2e-5:

  | cell | pole | equator |
  |---|---|---|
  | first active | 218 bar / 4217 K | 685 bar / 4996 K |
  | inner ghosts, up to | | 1091 bar / 5367 K |

  No active cell above 1e-8 bar is more than 2 % off the table.
- r(p_ref) sits at the tool's pole and equator radii to 1.2e-4 and 1.1e-5.

## Scripts

- `sparc_0925/run.sub` and `smoke.sub` take `PLANET=w121x1`. Arms default to
  `/viper/ptmp2/jinma/sparc_w121x1/<arm>`.
- The smoke runs in `/viper/ptmp2/jinma/sparc_w121x1/smoke`, which it clears first (arm dirs
  including STOP files).
- The pre-edit copies are `*.pre_w121x1`.
