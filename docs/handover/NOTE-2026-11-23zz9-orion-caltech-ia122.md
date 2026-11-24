# NOTE 2026-11-23zz9 orion -> Caltech: IA12.2 = caltech-ia122-465aff485 (USER-approved) for m35mms

## Binary
- Branch `caltech-ia122-465aff485` on jing-ze-ma/athenak.
- Build with CUDA 12.9 + g++ 13 (never 12.2). Tiny BSG gate first: pass 1 = 3.686e-11.

## New since IA12.1
- `rad_m1/clean_gas_secant = max_stage2`: HG -15 %, gates pass.
- `hydro/rad_signal_speed_diff = true`: the m35mms dt lever, your D1 / D1h / D2 arms.
- `rad_m1/vet_gd_wall_sym` and `rad_m1/vet_gd_theta_seam_sym`: off unless named.
  - wall_sym is bitwise on your 20 / 60 / 180 deg phi wedges.
  - theta_seam_sym removes a mirror defect at the theta-periodic seam (|F2 + mirror F2| 2.5e-2 -> 1.3e-9 of max|F1|). It changes the numbers in the seam rows.
  - USER decision: name both keys in all NEW theta-periodic runs. All three of your decks are theta-periodic.
- Gate table (all PASS): CPU 42/42 key-off bitwise; release suite 14/14; HG30 L_top 2.76 %; RSG 0.368 %; R136 within time error; bsg / m35mms / t9l62 smokes CW 0; MI300A parity.

## Decks
- In this repo: `docs/handover/ia122_decks/{tile,notile}/IA122_{bsg,m35mms,t9l62}.athinput` (MD5SUMS) plus `IA122_t9l62_restart.args`.
- The decks do NOT name the theta key. Add it on the command line as below.
- The M35ZAMS deck is held for the USER and not included.

## Do (m35mms)
1. Build + gate.
2. On ONE free node, from the newest m35mms production rst: two 300-cycle restart arms on IA12.2. Name every IA12.2 key (rst-stored keys otherwise win):
   - A: the IA122_m35mms keys.
   - B: A + `rad_m1/vet_gd_wall_sym=true rad_m1/vet_gd_theta_seam_sym=true`.
   Report CW / fallbacks / DT_RETRY / s/cyc / dt / L_top/L_in for each.
3. If B is clean (CW 0, no FATAL, s/cyc within 2 % of A): stop the IA11 m35mms production at its next rst and continue it on B. Keep the rst cadence.
4. t9l62: no change. W / WB / SB continue; t9l62 moves with IA12.3 (W's floor fix).
5. bsg: no change for now.

-- orion (coordinator)
