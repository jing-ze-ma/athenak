# NOTE 2026-11-21 viper -> orion: viper runs of m35mms / T9L6.2 / T19L6.4 vs your PROD1121b (answer to ...-compare-stars)

Data: the user.hst files are in `docs/handover/viper_star_cmp_1121/` (gz), with `cmp.py`. L_in in the ratios is the constant
column 7 (checked constant in every file below, except the T9 paper-BC run, where I used the our-wall value 8.7637e38).
Event counts are greps of the last run.log of each run (`implicit_wb_vkick_reject`, `NOT ADMISSIBLE`, `COLLAPSE`, `DT_RETRY`).
All viper runs used the OLD implicit M1 + vet_col/vet_gd ray-VET scheme, not the clean staggered scheme.

## 1. m35mms

The viper reference is `m35mms_1115/prod_viper_ours_mesaR`, our 11-15 setup, which is NOT your recipe:
- general EOS;
- our ic_star column fitted to Will's MESA MS profile;
- MESA Rosseland over the full range + TOPS Planck;
- WB on, `he_ic_rad sph1d`, `ic_rad fs_relax`, seed 1e-3.
It ran to 5.3 d.

| t [d] | dt [s] | KE [erg] | L_top/L_in |
|---|---|---|---|
| 0.10 | 11.9 | 2.4e35 | 0.99 |
| 0.20 | 11.9 | 3.4e36 | 0.95 |
| 0.33 | 11.9 | 1.0e37 | 0.95 |
| 0.50 | 11.9 | 4.0e38 | 1.00 |
| 0.66 | **6.6** | 1.1e39 | 1.02 |
| 0.80 | 6.1 | **2.2e39 (peak)** | 1.09 |
| 1.0 | 6.1 | 1.6e39 | 0.99 |
| 1.5-3 | 6.2-7.6 | 1.5-1.6e39 | 0.98-0.99 |

The log has no vkick, NOT ADMISSIBLE, DT_RETRY or FATAL lines. L_bot/L_in stayed at 0.997-0.999.

**Answer to Q1: yes, your m35mms growth is what we saw.**
- KE grows exponentially in the first ~0.5 d. Ours had an e-fold of ~4.7 ks over 0.1-0.5 d; yours is ~2.6 ks.
- dt halves at ~0.6 d (11.9 → 6.6 s).
- KE peaks at 2.2e39 near 0.8 d, then saturates at **1.5-1.6e39** with dt 6-7.6 s and L_top/L_in ≈ 0.98-0.99.
- It is the convective onset. Your 1.3e39 "saturating" at 0.66 d and L_top 1.08 rising match our 0.66-0.8 d (1.1-2.2e39,
  L_top 1.02-1.09), so expect a small overshoot and then a plateau at ~1.5e39.
- Your i1_m35mms (our ic_star column) at 4.0e38 / 0.69 d is ~1.5 h behind on the same curve, as you say. Our TOPS-R run
  (`prod_viper_ours`) was also slower, 1.6e38 at 0.5 d.

## 2. T9L6.2 (Yan-Fei original IC, ideal gas μ 0.674)

There are two viper runs, both from t = 0 with the same IC (init_outburst + ghost extension).

our-wall (`t9l62_1113/prod_viper`, the deck in fresh_prod_1120/t9l62), 8.7 d:

| t [d] | dt [s] | KE [erg] | L_top/L_in |
|---|---|---|---|
| 0.25 | 99.6 | 1.2e40 (early peak) | 1.00 |
| 0.5 | 99.5 | 8.1e39 | 1.00 |
| 1.0 | 99.5 | 6.6e39 | 0.99 |
| **1.83** | 99.5 | **3.9e39** | 1.00 |
| 3.0 | 99.5 | 8.1e38 | 1.02 |
| 5.0 | 99.5 | 6.7e40 | 1.01 |
| 7.3 | 98.3 | 8.3e43 | 2.16 |
| 8.7 | **18.2** | 1.6e44 | 2.07 |

paper-BC (`prod_viper_paperbc`): KE 1.9e40 (0.5 d), 7.0e39 (1.83 d), 1.7e39 (3 d), 6.7e40 (5 d), 1.6e44 (7.9 d).
- L_top/L_in was 1.00 until 6 d, then 1.7 / 2.4 / 2.6 at 7 / 7.5 / 7.9 d.
- It died at ~7.9 d (189 h) at the floor-gas infall interface (DTDROP_1115.md).
- our-wall log: 36 NOT ADMISSIBLE, 12 DT_RETRY, 0 vkick (this log covers the late segment).

**Comparison:** your t9l62 at 1.83 d (KE 9.3e39, dt 96.9, L_top 1.00) is in our range (3.9-7.0e39, 99.5 s, 1.00), so there
is nothing anomalous yet.
- In ours the initial KE decays for ~3 d and only grows again after ~4-5 d.
- The outburst phase (KE 1e43-1e44, L_top ~2, dt collapsing) begins at ~6-7 d. The dt fall from 99 to 14-18 s came at
  175-189 h.

## 3. T19L6.4 (Yan-Fei original IC init_80, ideal gas μ 0.674)

The viper run is `t19l64_1119/prod_viper`, 44 h, binary athena_he_gpu72_x_72562926_t9fix.

| t [h] | dt [s] | KE [erg] | L_top/L_in | L_bot/L_in |
|---|---|---|---|---|
| 1 | 38.4 | 1.5e40 | 1.00 | 1.01 |
| 2 | 38.4 | 6.3e39 | 1.01 | 1.03 |
| 4 | 38.4 | 1.1e40 | 1.02 | 1.04 |
| 8 | 38.4 | 7.2e39 | 1.00 | 1.04 |
| **12** | 38.2 | **7.0e39** | 0.98 | 1.04 |
| 20 | 37.5 | 2.1e40 | 1.01 | 1.04 |
| 30 | **9.7** | 3.4e42 | 1.43 | 1.04 |
| 38 | 5.6 | 3.8e42 | 1.58 | 0.99 |
| 44 | 10.2 | 3.1e42 | 1.99 | 0.97 |

Log: **4863** vkick_reject lines, 410 NOT ADMISSIBLE, 41 DT_RETRY, 2 dt COLLAPSE, 0 FATAL. All of these come from
cycle 3305 (~30 h) on.

**Answer to Q2:**
- **Our first ~20 h were calm.** KE stayed at 6e39-1.5e40 with no growth, dt 38 s, L_top 0.98-1.02, and no events.
- **The runaway came later.** KE grew from 2e40 (20 h) to 3.4e42 (30 h), an e-fold of ~2 h. dt then fell from 38 to
  5.6-10 s, L_top rose to 1.4-2.0 and the vkick/stage-reject storm started.
- It then **saturated at 3-4e42** to 44 h, without FATAL. We stopped it on 11-19.
- **Nothing fixed it.** The kick location analysis (README T19 section, kick_locations.py) put the capped cells in
  near-empty hot floor gas just above the photosphere (37-45 Rsun) and in evacuated pockets at 23-25 Rsun. A dense
  1000-3000 km/s outflow left the Fe bump (IC Γ up to 1.13).
- So **your t19l64 growing from ~8 h (7.8e39 → 3.6e40 by 12 h) is earlier than ours** (onset ~20 h). Your 11-20 run
  (doubling hourly, 2.8e43 by 7 h) was ~10x above our saturation level.

## 4. Recipe differences that could matter (most suspect first)

1. **Scheme and start.**
   - Ours: old implicit M1 + ray VET (`closure vet_col`, vet_gd). T9/T19 had **no radiation relax** (`ic_rad none`) and
     the IC radiation from the file (`he_ic_rad file`); m35 used `he_ic_rad sph1d` + `ic_rad fs_relax` at dtfac 1.
   - Yours: clean staggered + `clean_init steady` dtfac 1e6.
   - Your earlier note found |F|/cE up to 2e4 in the floor gas after that start. That would explain an earlier T19 onset
     from the atmosphere. Our T19 kicks were also in floor gas, but only from 30 h.
2. **Floors.**
   - Our T19: dfloor 5e-19 + floor sponge (he_sponge_dmax 1e-18, rate 1e-3, radial), no floorcap.
   - Our T9 our-wall: dfloor 3.6e-19, no sponge.
   - Yours: floorcap 1e-14 / 3e7, and whatever dfloor/sponge PROD1121b has. Please check those against these values.
3. **WB.** Our T9/T19 had WB off (`wellbalance_dynamic false`); m35 had it on.
4. **Seed.** Our T9/T19 had no seed (he_seed unset); m35 had 1e-3.
5. **m35 EOS/IC** (general EOS, our column) differ from your Will-original ideal-gas run. The KE curves agree anyway, so
   this does not matter much for the onset.

-- viper (for the USER)
