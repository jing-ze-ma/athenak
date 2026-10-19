# Disc-forming Algol candidates for the mass-gainer setup (literature scoping, 2026-10-07)

Status: literature only. No code, no jobs.

Conventions:

- **QUOTED**: a verbatim quote, with the source.
- **ESTIMATE**: my own computation. Tools:
  - `cand/ls75.py`: the Lubow & Shu (1975) fits r_min = 0.0488 q^-0.464 a and r_circ = 0.0859 q^-0.426 a, with q = M_donor/M_accretor, plus the exact accretor-L1 distance. Output in `cand/` (`systems.txt`, `systems_m16.txt`).
  - `scripts/roche_stream.py`: the ballistic L1 stream integrator. Output in `cand/stream_runs.txt`.
  - `cand/cost.py`: cost scaled from the RY Per stage-1 measurement. Output in `cand/cost_out.txt`.
  - Kepler's law gives a whenever the paper does not.

Sources are in `lit/`: PDFs plus pypdf text dumps. The line numbers below refer to those `.txt` dumps.

## 1. Thibault Lechien: papers and the systems he cites

arXiv author search (`au:Lechien`, 2023-2026) finds only two astro-ph papers by T. Lechien:

- **arXiv:2505.14780**: Lechien, de Mink, Valli, Rubio, van Son, Klement, Jin, Pols, *"Binary stars take what they get: Evidence for Efficient Mass Transfer from Stripped Stars with Rapidly Rotating Companions"*. This is the mass-transfer paper. Text: `lit/2505.14780.txt`.
- **arXiv:2608.02540**: Lechien et al., *"You're Gonna Need a Bigger Core"*. This one is about calibrating single-star core sizes and is not relevant here.

Two other papers turned up in the search, but Lechien is not an author of either: arXiv:2601.08508 (Wang et al., spin-up of thermal-timescale gainers) and arXiv:2602.06259 (Xing et al., disc-regulated accretion).

### Systems in 2505.14780, with the context of each citation

| System | Context, QUOTED |
|---|---|
| **phi Persei** (the famous one) | "The prototypical Be+sdOB system ϕ Persei provided crucial insights into the evolution of such systems. The inferred mass transfer efficiency was highly conservative, with the Be star in ϕ Persei having accreted at least 70% of the transferred mass (O. R. Pols 2007; A. Schootemeijer et al. 2018). This result challenged theoretical expectations, specifically those resulting from the rotationally limited mass transfer model." (§1, l.122-129). Table 1: "10516 ϕ Per 9.6 ± 0.3 1.2 ± 0.2 126.6982 ± 0.0035 D. Mourard et al. (2015) ** 0.36 0.56 - 1.00" (columns M_Be, M_sdOB, P [d], ref, tier, β_min, best-guess β). Fig. 1 is built on it: "Illustration of the origin of a ϕ Persei–like system". |
| **Regulus** | "At a lower mass and later evolutionary stage, another well-known example is Regulus, which was also found to have accreted at least 70% of the transferred mass (S. Rappaport et al. 2009)." (§1, l.132-136) |
| The rest of the Be+sdOB sample (Table 1) | FY CMa, 59 Cyg, HD 55606, κ Dra, LS Mus, κ Aps, V846 Ara, ι Ara, HR 2142, HD 161306, 7 Vul, 28 Cyg, V2119 Cyg, 60 Cyg, HR 6819. These are all post-transfer systems with no ongoing stream. |
| **Long-period Algols (spin)** | "tides seem to be insufficient to explain the rotation rate of longer-period systems (P > 5 days, A. Dervişoğlu et al. 2010)" (§4, l.586-588). Dervişoğlu et al. 2010 (arXiv:1003.4392, `lit/1003.4392.txt`) Table 1 lists the gainers with q, M, R and v_eq sin i: RY Per, TT Hya, AU Mon, RS Cep, RW Per, S Cnc, U Cep, β Per, and others. I use it in §2. |
| **Plaskett's star (by reference)** | "A strong large-scale magnetic fields has been detected in at least one star that recently gained mass in a binary system (J. H. Grunhut et al. 2013)" (§4, l.589-591). The reference is MNRAS 428, 1686 (l.1425-1427). **My identification** (not quoted): this is the HD 47129 = Plaskett's star magnetic-field paper. |
| Algols in general | "many authors have derived constraints from systems where mass transfer is ongoing, for example, in semi-detached Algols. Here, a wide range of efficiencies are required" (l.635-638). No individual Algol is named. |
| Disc spin regulation (the science link) | "The most promising explanation [...] is that the viscous accretion disk carries away AM efficiently, while still accreting matter on to the star (B. Paczynski 1991; R. Popham & R. Narayan 1991; ...)" (l.597-603) |

**Bottom line for task 1.** The famous system Lechien uses is **phi Per**. It is post-transfer (P = 126.7 d, Be 9.6 + sdO 1.2 Msun), so it has no L1 stream to simulate. Simulating it would mean a **phi-Per-progenitor** in its active Case-B phase. That is a model system: in such a wide orbit the stream misses the accretor by a large margin (r_min/R_acc >> 1, ESTIMATE). It has no observational anchor during the transfer phase. The active systems connected to Lechien's paper are the long-period Algols of Dervişoğlu et al. (2010).

## 2. Candidate parameters

### 2.1 Uniform sources (QUOTED rows)

**Mennickent et al. 2016, arXiv:1510.05628, Table 4 (W Ser stars)** (l.1406-1420). Columns are "system Po Mc Mh Mtot q Rc Rh Tc Th log(Lc) log(Lh) Ṁ/Ṗ", where c = donor and h = gainer. Ṁ is in Msun yr^-1 and Ṗ in s yr^-1.

- "β Lyr 12.94 3.0 13.2 16.2 0.22 15.2 6.0 13200 30200 3.81 4.42 1.6E-5 /18.93"
- "W Cru 198.5 1.2 7.8 9.0 0.16 76 4.0 5500 14000 3.68 2.74 4.4E-8 −1.3E-7/"
- "RX Cas 32.33 1.8 5.8 7.6 0.30 23.5 2.5 4400 - 2.29 - 6E-6 /19.86 K1III +A5eIII"
- "V367 Cyg 18.6 3.3 4.0 7.3 0.82 21.3 2.9 10400 14800 3.68 2.56 5-7E-5 /"
- "SX Cas 36.6 1.5 5.1 6.6 0.29 23.5 3.0 4000 - 2.10 - /-4.8 A6(shell) +K3III"
- "RS Cep 12.42 0.4 2.8 3.2 0.14 7.63 2.65 4610 9400 1.37 1.69 - /- B9.7eV +G8III"

**Table 3 (DPVs)** (l.1346-1353). Columns are "Mc Mh Mtot q Rc Rh Tc Th ... Ṁ".

- "AU Mon 1.2 7.0 8.2 0.17 10.1 5.1 5750 15900 ... 7.6E-6"
- "HD 170582 1.9 9.0 10.9 0.21 15.6 5.5 8000 18000 ... 1.6E-6"
- "V393 Sco 2.0 7.8 9.8 0.25 9.4 4.1 ... 9.5E-9"

**Table 12, disc radii** (l.1782-1799). Columns are "R1 e(R1) Rd e(Rd) a e(a) Rd/a ...".

- "W Cru WSer 4.00 - 126.00 - 299.0"
- "RX Cas WSer 2.50 - 18.80 - 40.0"
- "SX Cas WSer 3.00 0.40 37.00 - 87.1"
- "β Lyrae WSer 6.00 0.20 28.30 0.30 58.5"

The disc data are "from ... Djurašević 1993a [RX Cas], 1993b [SX Cas]". **Flag:** the RX Cas a = 40.0 in Table 12 contradicts Kepler's law for P = 32.33 d and Mtot = 7.6, which gives 84.0 Rsun (ESTIMATE). SX Cas, 87.1, agrees with Kepler. I use 84.0.

**Dervişoğlu et al. 2010, arXiv:1003.4392, Table 1** (l.69-112). Columns are "P/d q M1 M2 R1 R2 i vsyn sin i veq sin i F". The text adds: "veq sin i is the measured projected velocity and F = veq/vsyn".

- "RY Per 6.86 0.271 6.24 1.69 4.06 8.10 83 30 213 7.10"
- "TT Hya 6.95 0.224 2.63 0.59 1.95 5.87 84 15 164 10.9"
- "RS Cep 12.4 0.145 2.83 0.41 2.65 7.63 87 11 170 15.4"
- "RW Per 14.2 0.150 2.56 0.38 2.80 7.30 81 10 161 16.1"
- "S Cnc 9.48 0.090 2.51 0.23 2.15 5.25 83 12 174 14.5"
- "AU Mon 11.1 0.199 5.93 1.18 5.28 10.04 79 24 124 5.17"
- "β Per 2.87 0.217 3.70 0.81 2.74 3.60 82 51 52 1.02"
- "U Cep 2.49 0.550 3.57 1.86 2.41 4.40 88 56 437 7.80"

The same paper on discs: "Algols with P > (4 − 5) d are within the region (Fig. 3b) where we expect stars to have either a permanent or transient disc" (l.191-192).

### 2.2 Per-system anchors (QUOTED)

**β Lyrae (Mourard et al. 2018, arXiv:1807.04789)**

- Masses: "q = mg/md ≃ 4.50 (Harmanec & Scholz 1993) ... (md ≃ 2.9 M⊙, and mg ≃ 13.3 M⊙)" (§1, l.105-108).
- Mass-transfer rate, from the period change: "≈ 2.1 × 10−5 M⊙ yr−1 for a conservative transfer (Harmanec & Scholz 1993), and ≈ 2.9 × 10−5 M⊙ yr−1 for a non-conservative one ... as deduced from the large observed secular change of the orbital period of ≈ 19 s yr−1" (l.125-130).
- Table 7: "P0 (d) 12.913779", "a sin i (R⊙) 58.19", "Rg (R⊙) 6.0" (l.1540-1548).
- The gainer radius is ASSUMED, because the star is hidden: "The radius was set to a value typical for B0.5 IV-V star, Rg = 6 R⊙" (l.895).
- Rotation: the gainer "likely rotates close to its critical velocity" (l.882), an assumption. Mennickent & Djurašević 2013 (arXiv:1303.5812) also assume it: "For β Lyr we assumed critical rotation for the gainer" (l.278).
- Disc: "the radius of the outer rim is 30.0 ± 1.0 R⊙, the semithickness of the disk 6.5 ± 1.0 R⊙, and the binary orbital inclination i = 93.5 ± 1.0 deg" (abstract).
- Interferometry: Zhao et al. 2008 (arXiv:0808.0932), "The images clearly show the mass donor and the thick disk surrounding the mass gainer at all six epochs of observation" (CHARA/MIRC).

**W Crucis (Pavlovski et al. 2006, astro-ph/0603561)**

- "this long-period (198.5 days) eclipsing binary" (abstract).
- "M1 = 8.2 M⊙, and M2 = 1.6 M⊙ ... the separation of the components A = 306 R⊙" (§4, l.563-565).
- Accretor radius: "R1 < 17 R⊙ ... B-type MS star of mass about 8 M⊙ gives R ∼ 4 R⊙" (l.572-577).
- Disc: "The disk is geometrically very extended and its outer radius is about 80 % of the primary's critical lobe" (abstract); "Rd = 124 R⊙ ... zd,out = 17 R⊙" (l.570-572).
- No period change: "period change is not yet found for W Cru" (l.150-151). The disc comes from eclipse mapping, with no interferometry.
- Its fame: "Woolf drew similarity of W Cru to its more famous counterpart β Lyrae" (l.147-148).

**RX Cas, SX Cas**

- Parameters and Ṁ/Ṗ: the Mennickent+16 Table 4 and Table 12 rows above.
- "in RX Cas the long cycle lasts 516.1 days (Kalv 1979)" (1510.05628, l.642-643).
- Mass caveat (Davidge 2023, arXiv:2304.00080, Table 2, l.698-703): "RX Cas 9.2 1.5 3.0 7.6 ± 0.6" and "SX Cas 9.4 1.3 2.6 6.6 ± 0.6". The columns are log age, MSTO mass, total system mass, and the Mennickent+16 mass. The text says: "the initial system masses inferred from the MSTO ... are much smaller than the masses listed by Mennickent et al. (2016)" (l.750-751).
- SX Cas has a period *decrease*, "/-4.8".

**W Ser (Shepard et al. 2025, arXiv:2501.06982; Gies et al. 2025, arXiv:2504.10608)**

- "mass estimates of 2.0 M⊙ and 5.7 M⊙ for the donor and gainer" (2501.06982 abstract); "a R⊙ 48.7", "Rg R⊙ 3.8" (Table, l.825-829). The small gainer is a fit choice: "a small gainer (Rg = 3.8R⊙) surrounded by a large, optically thick disk (Rd = 14.3R⊙)" (l.167).
- Period change and Ṁ: "a period increase rate of 19.7 seconds per year ... implies a mass transfer rate of 1.7 × 10−5M⊙ yr−1" (2504.10608, l.96-97).
- The masses are disputed: "estimates of the total mass range from 2.5M⊙ (Mennickent et al. 2016) to 7.1M⊙ (Erdem & Öztürk 2014)" (2501.06982, l.93-94).

**V367 Cyg (Davidge 2022, arXiv:2209.10449)**

- Table 1 mass estimates span "19 ± 4 12 ± 3" down to "3.3 ± 0.9 4.0 ± 0.5" (l.144-151).
- "estimate a mass transfer rate of 5 − 7 × 10−5 M⊙/year" (l.194).
- Masses this uncertain mean poor anchors. Dropped.

**TT Hya (Miller et al. 2007, astro-ph/0611351, abstract)**

- "Doppler tomography of the observed H-alpha profiles revealed a distinct accretion disk"
- "lower limit to the mass transfer rate of 2e-10 solar masses per year"

**Algol (β Per; Baron et al. 2012, arXiv:1205.0754, Table 4)**

- "RA (R⊙) 2.90 ± 0.04 ... 2.73 ± 0.20", "MA (M⊙) 3.7 ± 0.3 ... 3.17 ± 0.21", "MB (M⊙) 0.81 ± 0.05 ... 0.70 ± 0.08"
- Inner orbit (Table 6): "P (days) 2.867328", "a (mas) 2.3 ± 0.1 2.15 ± 0.05". With "a parallax of 34.7 ± 0.6 mas" this gives a ≈ 13.3 Rsun (ESTIMATE).

**U Cep**: Dervişoğlu Table 1 row above. The same paper: "this system shows transient disc because the eclipse durations vary" (l.156-157).

## 3. Stream geometry (all ESTIMATE)

R_acc is the gainer radius from the sources above. The ballistic r_min comes from roche_stream.py (`cand/stream_runs.txt`); the other columns come from the LS75 fits (`cand/ls75.py`).

| System | P [d] | M_acc / M_don | q | a [Rsun] | R_acc | r_min | r_circ | **r_min/R_acc** | r_circ/R_acc | stream at R_acc |
|---|---|---|---|---|---|---|---|---|---|---|
| W Cru | 198.5 | 7.8 / 1.2 (M16) | 0.154 | 297.9 | 4.0 | 36.5 (ball.) | 56.8 | **8.7-9.1** | 14.2 | misses by far |
| RX Cas | 32.33 | 5.8 / 1.8 | 0.310 | 84.0 | 2.5 | 7.07 (ball.) | 11.9 | **2.83** | 4.75 | misses |
| SX Cas | 36.6 | 5.1 / 1.5 | 0.294 | 87.0 | 3.0 | 7.53 (ball.) | 12.6 | **2.51** | 4.20 | misses |
| S Cnc | 9.48 | 2.51 / 0.23 | 0.092 | 26.4 | 2.15 | 3.90 | 6.27 | **1.81** | 2.92 | misses |
| RS Cep | 12.4 | 2.8 / 0.4 | 0.143 | 33.3 | 2.65 | 4.09 (ball.) | 6.55 | **1.51-1.54** | 2.47 | misses |
| RW Per | 14.2 | 2.56 / 0.38 | 0.148 | 35.4 | 2.80 | 4.18 | 6.84 | **1.49** | 2.44 | misses |
| W Ser (M16 / Shepard) | 14.17 | 1.5/1.0 or 5.7/2.0 | 0.67 / 0.35 | 33.4 / 48.7 | 1.34 / 3.8 | 1.97 / 3.86 | 3.4 / 6.5 | **1.47 / 1.02** | 2.6 / 1.7 | mass-dependent |
| TT Hya | 6.95 | 2.63 / 0.59 | 0.224 | 22.6 | 1.95 | 2.21 | 3.68 | 1.13 | 1.88 | misses, barely |
| HD 170582 | 16.87 | 9.0 / 1.9 | 0.211 | 61.4 | 5.5 | 6.23 (ball.) | 10.2 | 1.13 | 1.86 | misses, barely |
| β Lyr | 12.94 | 13.2 / 3.0 | 0.227 | 58.7 | 6.0 (assumed) | 5.75 (ball.) | 9.48 | **0.96** | 1.58 | grazing: hits at 80 deg from radial, j/j_K = 1.27 |
| AU Mon | 11.11 | 7.0 / 1.2 | 0.171 | 42.3 | 5.1 | 4.72 (ball.) | 7.70 | 0.92 | 1.51 | grazing, 75 deg, j/j_K = 1.24 |
| RY Per (current) | 6.86 | 6.24 / 1.69 | 0.271 | 30.3 | 4.06 | 2.69 (ball.) | 4.54 | 0.66 | 1.12 | impact, 58 deg, j/j_K = 1.07 |
| Algol | 2.87 | 3.70 / 0.81 | 0.219 | 14.0 | 2.74 | 1.39 | 2.30 | 0.51 | 0.84 | direct impact |
| U Cep | 2.49 | 3.57 / 1.86 | 0.521 | 13.6 | 2.41 | 0.90 | 1.54 | 0.37 | 0.64 | direct impact |

Cross-check: Mennickent+16 Table 12 gives "Rd/a" of 0.47 for RX Cas, 0.42 for SX Cas, 0.44 for W Cru and 0.48 for β Lyr. All observed discs are much larger than r_circ.

## 4. Cost (ESTIMATE, `cand/cost.py`)

Basis:

- Same cell shape as RY Per stage 1: dlnr = dphi = 2π/2048, 4 θ cells, r_out = 0.85 d_L1.
- dt scales as r_in·dphi / v_ff(r_in), as in stage 1, which IMPL2.md calls "stream-limited 6.58e-6".
- Cost per orbit = 23 min × (nr/448) × (cycles per orbit / 1.29e5), on one apu node.
- With r_in = R_acc in a disc system, the speed at r_in is about v_K, not v_ff, so these numbers are pessimistic by about √2.
- "Stream cells" = cs_don/Ω divided by r_out·dphi. RY Per has 17, so a system with more can use a coarser grid. Coarsening by a factor f cuts the cost by about f³.

| System | r_in | nr | cycles/orbit | min/orbit (1 node) | stream cells | note |
|---|---|---|---|---|---|---|
| RY Per (ref) | 4.06 = R | 452 | 1.29e5 | 23 | 17 | measured basis |
| β Lyr | 6.0 = R | 548 | 2.0e5 | **44** | 24 | coarsened 1.4x: ~16 |
| RS Cep | 2.65 = R | 650 | 3.1e5 | 81 | 22 | |
| HD 170582 | 5.5 = R | 595 | 2.5e5 | 60 | 23 | |
| RX Cas | 2.5 = R / 4.24 = 0.6 r_min | 936 / 763 | 1.3e6 / 5.9e5 | **496 / 180** | 25 | coarsened 1.45x: ~165 / ~60 |
| SX Cas | 3.0 = R / 4.52 = 0.6 r_min | 890 / 757 | 1.1e6 / 5.7e5 | 380 / 172 | 26 | coarsened 1.5x: ~110 / ~50 |
| W Cru | 4.0 = R / 21.9 = 0.6 r_min | 1227 / 673 | 4.8e6 / 3.5e5 | 2324 / 95 | 44 | coarsened 2.5x: ~150 / ~6 |

Physical clock. One orbit is 13 d for β Lyr, 32 d for RX Cas and 198 d for W Cru. A long period means many inner Keplerian times per orbit, and that is where the cost goes.

With r_in = 0.6 r_min the run sees the disc but not the boundary layer at the star. The spin-up science (the AM flux into the star) needs r_in = R_acc. Even so, an absorbing sphere at R_acc does not have to resolve the photosphere when r_min > R_acc.

## 5. Ranking

| Rank | System | (a) r_min/R | (b) anchors | (c) fame / science | (d) cost (r_in = R) | verdict |
|---|---|---|---|---|---|---|
| 1 | **RX Cas** | 2.8 | Ṗ = 19.86 s/yr gives Ṁ = 6e-6 Msun/yr; disc radius 18.8 Rsun (Djurašević 1993a); 516-d long cycle. No interferometry, no v sin i. Masses disputed (Davidge MSTO total 3.0 vs 7.6). | archetypal W Ser / long-period Algol; disc only | ~165 min/orbit coarsened (~7x RY Per) | best "easy and anchored" choice |
| 2 | **β Lyr** | 0.96, grazing (R_g assumed) | the best: Ṗ ≈ 19 s/yr gives 2-3e-5 Msun/yr; CHARA images of the disc; disc rim 30 Rsun and half-thickness 6.5 Rsun | the most famous; gainer assumed near-critical, jets, hidden B0.5 star | ~16-44 min/orbit, the cheapest | best science, but fails the r_min/R >~ 1.5 criterion. The stream arrives tangentially and super-Keplerian, which is milder than RY Per. |
| 3 | SX Cas | 2.5 | Ṗ = -4.8 s/yr (decreasing, puzzling); disc 37 Rsun; shell star | moderate | ~110 min/orbit coarsened | RX Cas twin; the negative Ṗ is awkward |
| 4 | RS Cep / RW Per / S Cnc | 1.5-1.8 | measured gainer spin: "veq sin i" 170 / 161 / 174 km/s, F = 15.4 / 16.1 / 14.5 (Dervişoğlu Table 1). No Ṁ quoted here. | spin-up link (fast gainers, cited via Lechien → Dervişoğlu); obscure | RS Cep 81 min/orbit | the best spin-up anchor, least famous |
| 5 | W Cru | 8.7 | weak: no Ṗ, model Ṁ 4.4e-8, eclipse-mapped disc | "counterpart of β Lyr" | 2324 min/orbit at R; ~6-95 at 0.6 r_min | cleanest disc, but the star is unreachable |
| - | TT Hya, HD 170582 | 1.13 | Doppler-tomography disc (TT Hya) | | 60 (HD 170582) | margin too small |
| - | W Ser, V367 Cyg | 1.0-1.5 | masses disputed by factors of 3 | | | rejected: geometry undetermined |
| - | AU Mon, RY Per, Algol, U Cep | < 1 | | Algol and U Cep famous | | impact systems, rejected |
| - | phi Per (Lechien) | >> 1 for a progenitor | no active-phase data | the most famous post-transfer system | | possible later "phi-Per progenitor" model run |

## 6. Recommendation

- **RX Cas** is the first easy target. It meets r_min/R_acc ≈ 2.8 with r_circ/R ≈ 4.8, so an absorbing sphere at R_acc = 2.5 Rsun works with no photosphere to resolve. It has a measured Ṗ/Ṁ and a fitted disc radius. Cost is about 1.5-3 h per orbit on one apu node, ESTIMATE.
- **β Lyr** is the famous, high-science follow-up. Its cost per orbit is about RY Per's. It is grazing (r_min/R ≈ 0.96 with the assumed R_g = 6), so it carries some of the RY Per difficulty, but much milder: the stream arrives nearly tangential and super-Keplerian.
- For the spin-up / Plaskett narrative, the Dervişoğlu fast rotators (RS Cep, F = 15) can serve as a cross-check.
- A phi-Per-progenitor run is a model system that ties directly to Lechien's β result.

---

## 7. Follow-up (2026-10-07): massive gainers (>~ 8 Msun) that form a disc

New sources in `lit/`:

- 1803.02379: Schootemeijer+2018, phi Per.
- 2609.25526: Wade+2026, Plaskett's star.
- 1209.6326: Grunhut+2013.
- 2601.08508: Wang+2026.
- Mennickent DPV papers: 2603.20413, 1905.04231, 2009.01352, 2201.13275, 1112.2668, and others.

New scripts in `cand/`:

- `phiper_epochs.py` (output `phiper_epochs.txt`)
- `plaskett_epochs.py` (output `plaskett_epochs.txt`)
- `systems_massive.txt` (run through ls75.py)
- `cost2.py` (output `cost2_out.txt`)

### 7.1 phi Per progenitor

**What the paper gives (QUOTED, Schootemeijer et al. 2018, arXiv:1803.02379)**

- Today, Table 1: "M [M⊙] 1.2 ± 0.2 9.6 ± 0.3" (sdO, Be) and "Porb [d] 126.7 ± 0.001".
- Rotation: "very rapidly rotating (vrot = 0.93 ± 0.08 vcrit) Be star" (l.203).
- Progenitor: "the most probable initial masses of the system are 7.2 ± 0.4 M⊙ for the primary star and 3.8 ± 0.4 M⊙ for the secondary star" (l.529-530).
- "For this initial orbital period, 16 ± 4 days, the system was wide enough to remain detached during the main sequence evolution of the primary star. Mass transfer started after the initial primary started crossing the Hertzsprung gap, so-called case B mass transfer." (l.537-540)
- "The system must have evolved through near-conservative mass transfer." (abstract)
- Accretion is capped at "˙Macc = min(−˙Mdonor, 10Ma/τKH)" (eq. 1).
- On discs: "We thus assume implicitly assume that the accretion disk regulates the angular momentum that is accreted by the mass gainer" (l.281-282).
- On swelling: "the accreting star in such systems swells significantly during the mass transfer phase (Neo et al. 1977) such that it may possibly even fill its Roche lobe too" (l.782-784).
- Earlier work: "Vanbeveren et al. (1998a) and Pols (2007) ... inferred initial masses for the progenitors of around 6 M⊙ for the SES and 5 M⊙ for the Be star, and an initial orbital period in the order of 10 days" (l.543-546).

**What the paper does not give.** The paper does NOT tabulate any state during the transfer: no epoch masses, no gainer radius and no Ṁ. Vinciguerra+2020 is a Be X-ray-binary population study, not a phi Per model. Pols 2007 is not on arXiv. The epochs below are therefore ESTIMATEs.

**Epoch construction (all ESTIMATE, `cand/phiper_epochs.py`)**

- Conservative track, with P = P0 (M1,0 M2,0 / (Md Ma))^3 from 7.2 + 3.8 Msun and 16 d.
- This track ends at 201 d, against the observed 126.7 d. That is consistent with "near"-conservative, i.e. slightly lossy, transfer.
- R_TE = M^0.6 Rsun stands in for a thermal-equilibrium B star. 2 R_TE stands in for an accretion-swollen star.
- The LS75 fits are extrapolated for q > 1.

| Md | Ma | q | P [d] | a [Rsun] | r_min | r_circ | R_TE | r_min/R_TE | r_circ/R_TE | r_min/(2R_TE) | regime |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 6.5 | 4.5 | 1.44 | 13.1 | 52.0 | 2.14 | 3.82 | 2.47 | 0.87 | 1.55 | 0.43 | impact (rapid thermal phase) |
| 5.5 | 5.5 | 1.00 | 11.8 | 48.6 | 2.37 | 4.18 | 2.78 | 0.85 | 1.50 | 0.43 | impact |
| 4.5 | 6.5 | 0.69 | 13.1 | 52.0 | 3.01 | 5.22 | 3.07 | 0.98 | 1.70 | 0.49 | grazing |
| 3.5 | 7.5 | 0.47 | 18.1 | 64.6 | 4.49 | 7.67 | 3.35 | 1.34 | 2.29 | 0.67 | disc only if not swollen |
| **2.5** | **8.5** | 0.29 | **34.2** | 98.5 | 8.48 | 14.3 | 3.61 | **2.35** | 3.95 | 1.17 | **disc** |
| 1.8 | 9.2 | 0.20 | 72.2 | 162 | 16.9 | 27.9 | 3.79 | 4.46 | 7.37 | 2.23 | disc even if swollen |
| 1.2 | 9.8 | 0.12 | 201 | 322 | 41.6 | 67.6 | 3.93 | 10.6 | 17.2 | 5.3 | end of transfer |

Ṁ at these epochs: none is published. My estimate is that the late, slow Case-B phase runs at about 1e-6 to 1e-5 Msun/yr, and the early rapid phase at up to the cap 10 Ma/τKH.

**Recommended epoch: Md ≈ 2.5, Ma ≈ 8.5 Msun, P ≈ 34 d.**

- The disc criterion is met (r_min/R ≈ 2.35) unless the gainer is still swollen by more than a factor 1.6.
- Cost (ESTIMATE, cost2.py): 339 min/orbit at RY Per resolution. Coarsened to the stream width (f = 1.75), about 63 min/orbit on one apu node.
- The Md ≈ 1.8 epoch (Ma 9.2, P 72 d) is robust even if the gainer is swollen. It costs about 88 min/orbit coarsened.
- This is a model system with no observational anchor during the transfer phase.

### 7.2 Plaskett's star (HD 47129) progenitor

**Today (QUOTED, Wade et al. 2026, arXiv:2609.25526)**

- Abstract: "a partially Roche-lobe-filling stripped star with a mass of 5.9 +1.8 −1.4 M⊙, while the secondary is a 40.8 +9.3 −6.7 M⊙, rapidly rotating magnetic accretor".
- Table 1:
  - "P[d] 14.39626"
  - "i[◦] 48±4"
  - "a[R⊙] 78.4 +5.5 −4.5 | 11.4 +1.9 −1.8" (each component's orbit; the sum, 89.8 Rsun, is my addition)
  - "R[R⊙] 12.3±1.2 10.4±1.4"
  - "logL[L⊙] 5.10±0.07 5.12±0.11"
  - "RRL [R⊙] 20.6 | 49.3"
  - "vsini 53.2±2.7 | 350±50"
  - "q=M2/M1 (RVs) 6.9 +1.2 −0.9"
- Rotation: "v eq,2∼430 km s−1, approximately 60% of critical, corresponding to a rotational period of 1.21 days" (l.612-613).
- **Disc vs impact is left open:** "Whether the flow forms an accretion disk or impacts the star directly, angular momentum is deposited in the outer layers of the accretor" (l.733-735).

**Magnetic field (QUOTED, Grunhut et al. 2013, arXiv:1209.6326, abstract)**

- "the rapidly-rotating secondary component is the magnetized star"
- "implying a minimum surface dipole polar strength of 2850 ± 500 G"

**Published evolution model (QUOTED, Wade+2026 §5 and Fig. 3)**

- Initial state: "M1,i = 18.2M⊙ M2,i = 16.8M⊙ Pi = 3.69 [d]"; final state: "M1,f = 6.62M⊙ M2,f = 26.4M⊙ Pf = 14.4 [d]".
- The text gives slightly different initial values: "initial masses of 18M⊙ and 16.2M⊙ (qi = 0.9) with an initial period of 3.81 days".
- Transfer happens in two phases: "an initial phase of mass transfer during the main-sequence (case A) followed by a final stripping phase afterward (case AB)".
- The model misses the observed q: "No solutions were found for q ≥ 5" and "reaches a final mass ratio of q≃4".
- At the start of accretion: "when the secondary starts accreting it has M≃16M⊙ and R≃9R⊙. The accretion rate is ˙M≈10−4M⊙/yr" (l.710-712).
- Their eq. 2 for conservative transfer: "P/Pi = (q/qi)^−3 ((1+q)/(1+qi))^6".

**Epochs (ESTIMATE, `cand/plaskett_epochs.py`)**

- Eq. 2 from 18.2 + 16.8 Msun and 3.69 d.
- Gainer radius R = 9 (M/16)^0.6, anchored on the quoted 9 Rsun at 16 Msun.

| Md | Ma | q_don | P [d] | a | r_min | r_circ | R_acc | r_min/R | r_circ/R |
|---|---|---|---|---|---|---|---|---|---|
| 17.0 | 18.0 | 0.94 | 3.68 | 32.8 | 1.64 | 2.89 | 9.7 | 0.17 | 0.30 |
| 15.0 | 20.0 | 0.75 | 3.91 | 34.1 | 1.90 | 3.32 | 10.3 | 0.19 | 0.32 |
| 12.0 | 23.0 | 0.52 | 5.02 | 40.3 | 2.66 | 4.57 | 11.2 | 0.24 | 0.41 |
| 10.0 | 25.0 | 0.40 | 6.75 | 49.2 | 3.67 | 6.24 | 11.8 | 0.31 | 0.53 |
| 8.0 | 27.0 | 0.30 | 10.5 | 65.9 | 5.65 | 9.50 | 12.3 | 0.46 | 0.77 |
| 6.62 | 28.4 | 0.23 | 15.9 | 87.1 | 8.35 | 13.9 | 12.7 | 0.66 | 1.10 |
| today (5.9 / 40.8, a = 89.8, R = 10.4) | | 0.145 | 14.40 | 89.8 | 10.75 (ballistic 11.02) | 17.6 | 10.4 | **1.03-1.06** | 1.69 |
| today with model q≃4 (6.62 / 26.4) | | 0.25 | 14.40 | 79.9 | 7.41 | 12.4 | 10.4 | 0.71 | 1.19 |

**Verdict.** In this model Plaskett's progenitor is a **direct-impact** system throughout the transfer: r_min/R ≤ 0.66 and r_circ/R ≤ 1.1. The Case-A phase is extreme, with R_acc/a ≈ 0.3. Today's geometry, if transfer were still going on, is grazing: ballistic r_min = 11.0 Rsun against R = 10.4. This confirms the coordinator's rough estimate. Plaskett is therefore an impact-accretion target, harder than RY Per: R/a is larger, and the gainer is a ~40 Msun O star.

**Eddington factor (ESTIMATE).** I use Γ_e = 10^-4.813 (1+X) L/M with today's L = 10^5.12 Lsun and M = 40.8 Msun. This gives Γ_e = 0.086 for X = 0.73 and 0.071 for X = 0.42.

During rapid accretion, Wang et al. 2026 (arXiv:2601.08508, l.1028-1035) write: "The Eddington factor before mass transfer is typically small, for example, less than 0.3 for a 12 M⊙ star, and increases rapidly during mass transfer ... a 12 M⊙ star model with an initial accretion rate of 4.0×10−3 M⊙ yr−1 attains an Eddington factor of about 0.7 at critical rotation, whereas the corresponding model with an initial accretion rate of 1.0×10−3 M⊙ yr−1 reaches only approximately 0.45".

So radiation matters at the gainer surface and in the hot spot during Ṁ ~ 1e-4. For the stream and disc dynamics at Γ_e ≲ 0.1, it is a correction of order 10 %.

### 7.3 Observed massive-gainer systems (gainer >~ 8 Msun)

r_min/R and r_circ/R are ESTIMATEs from `cand/systems_massive.txt` via ls75.py.

| System | QUOTED parameters (source) | r_min/R | r_circ/R | regime |
|---|---|---|---|---|
| OGLE-LMC-DPV-062 | "M1 (M⊙) 10.36±1.14", "M2 (M⊙) 2.73±0.30", "R1 (R⊙) 5.14±0.98", "P o (d) 6.904858", "a orb (R⊙) 35.9±1.4" (arXiv:2603.20413, Table) | 0.63 | 1.06 | impact |
| OGLE-LMC-DPV-065 | "13.8 and 2.81 M⊙ stars of radii 8.8 and 12.6 R⊙", "semimajor axis is 49.9 R⊙", "P [d] 10.0316267"; "the disk has a radius of 25 R⊙" (arXiv:1905.04231) | 0.58 | 0.96 | impact |
| UU Cas | "Mh[M⊙] 17.4±0.2", "Mc[M⊙] 9.0±0.2", "Rh[R⊙] 7.0±0.1", "a stellar separation of 52 R⊙", "P o = 8d.519296(8)", "There is no evidence for orbital period change" (arXiv:2009.01352) | 0.49 | 0.85 | impact |
| V393 Sco (7.8 Msun) | "Mh[M⊙] 7.8±0.2", "Rh[R⊙] 4.4±0.2", "aorb[R⊙] 35.1±0.3" (arXiv:1112.2668) | 0.73 | 1.22 | impact |
| HD 170582 | Mennickent+16 Table 3 (§2.1) | 1.12 | 1.85 | marginal |
| β Lyr | §2.2 | 0.95 | 1.58 | grazing |
| **W Cru** (7.8-8.2 Msun) | "M1 = 8.2 M⊙, and M2 = 1.6 M⊙ ... A = 306 R⊙", "R ∼ 4 R⊙" (astro-ph/0603561, §2.2) | **8.0-8.7** | 13-14 | **disc** |
| **BY Cru** (9.1 Msun) | "BY Cru 106.4 1.7 9.1 10.8 0.19 52 - 11000 -" (Mennickent+16 Table 4: P, Mc, Mh, Mtot, q, Rc, Rh = "-" i.e. not given, Tc) | **4.9 (R = 4.5, ESTIMATE) to 2.5 (R = 9)** | 4-8 | **disc** (gainer radius unmeasured) |

The pattern is clear. Every well-measured massive DPV with P < 17 d is an impact or grazing system. The only observed massive gainers that miss the stream are the long-period W Ser systems W Cru (198.5 d) and BY Cru (106.4 d).

- W Cru's anchors are weak: no period change, an eclipse-mapped disc, no interferometry.
- BY Cru has no published gainer radius. Its parameters come from Daems+1997 via Mennickent+16, and are not on arXiv.

Cost (ESTIMATE, cost2.py, r_in = R_acc, coarsened to the stream width):

| System | min/orbit |
|---|---|
| Plaskett today-geometry | ~12 (impact) |
| HD 170582 | ~26 |
| BY Cru | ~43 |
| phi-Per Md 2.5 | ~63 |
| phi-Per Md 1.8 | ~88 |
| W Cru | ~142 |

### 7.4 Three-way comparison and recommendation

| Line | gainer mass | r_min/R | anchors | cost (min/orbit, coarsened) |
|---|---|---|---|---|
| phi Per progenitor, Md 2.5 / Ma 8.5, P 34 d | 8.5 | 2.35 (1.17 if swollen ×2) | progenitor fit only; famous today; Be at 0.93 v_crit | ~63 |
| Plaskett progenitor (Wade+26 track) | 18 → 28 (40.8 today) | 0.17-0.66; today 1.03 | magnetic field, v_eq 430 km/s, stripped donor | impact: needs the stage-2 envelope physics |
| Observed W Cru / BY Cru | 8.2 / 9.1 | 8.0 / 2.5-4.9 | weak (W Cru) or thin (BY Cru) | ~142 / ~43 |

**Recommendation**

- **The phi-Per-progenitor epoch Md ≈ 2.5, Ma ≈ 8.5 Msun, P ≈ 34 d** as the massive disc case. It is a clean disc (r_min/R ≈ 2.4 with r_circ/R ≈ 4) with a famous end state and Lechien's β result as the science hook. It costs ~1 h per orbit on one node after coarsening.
- If the swelling worry dominates, use Md ≈ 1.8 instead (P 72 d, ~1.5 h per orbit).
- Plaskett belongs to the impact class. It needs the RY Per stage-2 envelope approach, and Γ_e ~ 0.1 (up to ~0.5-0.7 during rapid accretion per Wang+26) argues for radiation.
