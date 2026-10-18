# NOTE 2026-10-18 (Raven): accretor ray-traced donor irradiation (USER-approved fix A)

## Code and gates
- **Branch:** accretor-raydonor-1017 @ 72a65c02 (code 20f66126), off accretor-sector-light-1017 5eeb137e. Pushed, NOT merged.
- **Key:** problem/donor_irr_mode = m1 (default, bitwise) | ray. Code in src/pgen/ry_per_accretor_raydonor.hpp plus hooks in ry_per_accretor.cpp.
- **Design:** accretor_rhd_design.md sect. 17.
- **Gate:** vs 62119e93 (S/A/G2): BITWISE.
- **Merges:** with accretor-3d-1017 89df6f64 and accretor-sector-light-1017 dfc0a217 the only conflicts are in the docs (appended sections).

## Method
- **Source:** the exact 3-D Roche L1 lobe with the donor_gdark law, I = (c/4pi) a T_loc^4. There is no dilution or calibration factor (donor_irr_scale is not applied).
- **Directions, per band column:** 24 azimuths x 4 Gauss-Legendre elevations up to the lobe's grazing angle. The geometry is tabulated once (354 MB/rank).
- **Attenuation:** long-characteristic march through a gathered theta-mean extinction map. The band is z-uniform, so tau_3D = tau/cos e. The night-side shadow comes out exactly (q = 0).
- **Deposit:** all of it goes into M1 E, and the implicit coupling heats the gas. Explicit gas heating is unstable in the thin atmosphere, where kappa_P ~2000 cm^2/g.
- **Momentum:** through kappa_T only.
- **Ledger:** Eirr/Ebnd closes.
- **The lobe holds no light.**
- **3-D wedge:** not implemented (FATAL).

## Results
Runs are in /raven/ptmp/jinma/accray_1017, NOTE.md.
- **Static:** with donor_ray_opac = rosseland the facing F_inc matches the analytic gdark value to 0.8 %.
- **FINDING 1 (calibration target):** the gdark facing flux is **7.1e4, 0.47x the uniform-sphere 'target' 1.52e5**. The L1 face is gravity-darkened. Reference script: tool/ref_flux.py.
- **FINDING 2 (floor opacity artefact):** with donor_ray_opac = total, the 1e-8 floor's Planck mean (300-1600 cm^2/g, from the low-rho table extension) gives tau ~0.25 between the stars and absorbs ~2/3 of the beam.
- **Lateral beam vs vet_col:**
  - The vet_col donor part had dE 3.2-4.5 at all phi.
  - With rays the direct beam has none on the night side.
  - The reprocessed light (deposited at the facing photosphere, then carried by M1/vet_col) still gives a night-side dE of 5 % of E.
- **Stream (tot arm, to t 0.0319):** behaves like the surface arm.
  - dt 7.02e-6 flat to contact.
  - Contact: an eint<=0 burst, 58 retries at dt_min 1e-8.
  - The impact and its heating match the surface arm.
  - The ros arm was continued on Viper.
- **Cost:** +0.3-0.45 s/cycle on 4 A100 at donor_ray_every = 1 (the update itself is 0.27 s); see 17.4.

## Open (USER)
- the opac default;
- the every-N setting for production;
- the residual reprocessed-light night dE;
- 3-D ray mode;
- the s calibration target.
