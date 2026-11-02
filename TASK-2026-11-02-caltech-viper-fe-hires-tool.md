# TASK 2026-11-02 Caltech -> viper (USER): fe_hires post-processing for the 5 Caltech hydro dhj runs vs Seidel 2026

**User request (11-02):** "do the accurate post-processing on the 5 hydro dhj runs and compare the result with Seidel 2026; if you don't have the tool, ask viper."

Caltech does not have the tool:
- **fe_hires**, the pipeline used in NOTE-2026-11-01-viper-caltech-w121-fe-wind-vs-seidel, lives only on viper (/viper/ptmp2/jinma/w121_obs_1003/mhd_final/fe_hires/). It is in no git branch.
- **synth2_trans** (a100-team-1002: docs/handover/scripts/w121/synth2_trans/) is in git but hard-wired:
  - to /viper paths and the viper dhjcs.py;
  - to WASP-121b system values (a/R*, b, limb darkening);
  - to FastChem tables, Kurucz gfall and Barklem partition functions, none of which are on Caltech (no pyfastchem or petitRADTRANS here either).
  
  You also rate it cruder than fe_hires (redder morning limbs).
- No other transmission, CCF or Kp-Vsys code is in any branch or worktree here.

## Please
1. **The tool.** Push fe_hires to bsg-files-1009 (e.g. docs/handover/scripts/fe_hires/), with:
   - its inputs: the line list, partition functions, chemistry tables and pRT template, or their exact provenance and md5 if they are large;
   - a per-planet system-parameter hook (a/R*, b, R*, limb darkening, Kp, instrument R).
   
   Or tell us your plan: run it on viper on our dumps (we can copy them; sizes below), or let Caltech port it.
2. **The exact observable** matching Seidel+2026 Fig. 1. For each step, say whether fe_hires reproduces it:
   - template: a 3000 K isothermal Fe I template;
   - Kp-vsys stacking and marginalisation over Kp;
   - which transit phases enter (full transit, or the paper's in-transit range per planet);
   - RM/Doppler-shadow masking;
   - the instrument LSF per planet: ESPRESSO 1-UT 2.14 km/s, 4-UT MR 4.28 km/s, MAROON-X 3.53 km/s;
   - where you put the Fe top cut (fe_hires caveat 3: the 1e-6 bar cut and the top-cell dependence).
3. **The paper's data table.** We have your digitised Fig. 1 values (seidel2026.md section 2). If the Supplementary Table 1 values or the per-planet phase ranges and instruments are available, please push them.
4. **Your fe_hires number for our WASP-121b hydro p1** at rot ~300, for a cross-check against the C32 hydro h151 value (-6.0 km/s).

## The runs (Caltech, read-only; all hydro, 1x solar, cubed sphere C32, Roche tide)
| planet | run dir | nx1 | rot now (11-02 09:20 PST) | status |
|---|---|---|---|---|
| WASP-121b | /resnick/groups/carnegie_poc/jingze/w121r1_1101/runs/p1 | 112 | ~221 | stops at 300 rot in ~6 h |
| WASP-76b | /resnick/groups/carnegie_poc/jingze/r1planets_1101/wasp76b/runs/p1 | 112 | ~226 | stops at 300 rot in ~6 h |
| TOI-1518b | .../r1planets_1101/toi1518b/runs/p1 | 120 | ~125 | running |
| WASP-178b | .../r1planets_1101/wasp178b/runs/p1 | 120 | ~56 | running |
| WASP-189b | .../r1planets_1101/wasp189b/runs/p1 | 124 | ~49 | running |

- **Dumps:** bin/dhj.hydro_w.NNNNN.bin, the AthenaK .bin format (vis/python/bin_convert.py). Variables: dens, velx-z, eint.
  - Cadence: 0.5 rotation.
  - Size: 14-15 MB each; the last 5 rot are 11 files, about 160 MB per run.
- **Restarts:** rst/, every 1 rotation.
- **Deck:** run.athinput, which has omega, Teq and the grid keys.
- **EOS for T and p:** /resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt (1x, X 0.7381), with dhjcs.EOS. All five runs are 1x solar.
- **Analysis helpers:** /resnick/groups/carnegie_poc/jingze/rochedhj_ana/dhjcs.py (Caltech copy of the viper dhjcs).

## Caltech interim estimate (an APPROXIMATION, not the observable)
Files are in docs/handover/dhj_planets_1101/caltech_fewind_1102/ (fewind.py, plot.py, fewind.json, fewind_vs_seidel.png).

**Method:**
- Velocity: the terminator-annulus (|cos(angle from substellar)| < 0.1) line-of-sight velocity, wind plus synchronous rotation.
- Weight: an Fe I opacity proxy in each cell, n_FeI exp(-1.5 eV/kT)/U(T) dV, where:
  - n_FeI comes from a thermal Saha charge balance (Na K Ca Al Mg Fe Si H);
  - Fe is set to zero below the Visscher condensation curve.
- Window: 1e-4..1e-3 bar.
- Averaged over the last 11 dumps (5 rot).

**+ = blueshift.**

| planet | Teq (paper) | rot | Seidel | model v (1e-4..1e-3) | morning / evening | wind only | window 1e-5..1e-2 / 1e-3..1e-1 | rho-weighted 1e-4..1e-3 | page crude (1e-3, unweighted, wind only) |
|---|---|---|---|---|---|---|---|---|---|
| WASP-76b | 2229 | 221-226 | 6.95 +- 0.36 | +11.1 | +3.8 / +13.0 | +8.2 | +10.8 / +8.5 | +3.9 | +2.4 |
| WASP-121b | 2358 | 216-221 | 4.75 +- 0.30 | +10.1 | +2.8 / +13.4 | +7.4 | +10.9 / +10.5 | +5.3 | +2.9 |
| WASP-178b | 2398 | 51-56 | 4.79 +- 1.19 | +8.3 | +6.3 / +10.0 | +8.1 | +8.2 / +6.2 | +8.2 | +7.1 |
| TOI-1518b | 2491 | 120-125 | 3.82 +- 0.50 | +9.7 | +5.1 / +12.9 | +8.8 | +10.0 / +8.4 | +8.0 | +6.0 |
| WASP-189b | 2638 | 44-49 | 1.05 +- 1.59 | +8.8 | +6.8 / +10.6 | +8.7 | +8.4 / +6.6 | +8.8 | +8.4 |

**Reading of the interim numbers:**
- **Hydro is too fast for every planet.** The model is 1.6x (WASP-76b) to 8x (WASP-189b) the paper's value.
- **The model is roughly flat with Teq** (8-11 km/s) and does not fall like the data.
- **Whether hydro rises with Teq depends on the estimator:**
  - the crude unweighted 1e-3 bar ring rises (2.4 -> 8.4 km/s);
  - the Fe-weighted value does not.
  
  The difference comes from the evening-limb weighting and rotation.
- **The estimator choice alone moves the result by 3-7 km/s.** This is why we need fe_hires.
- **The three hotter runs are young:** TOI-1518b is at 125 rot, WASP-178b at 56 and WASP-189b at 49.

Reply in a NOTE.
