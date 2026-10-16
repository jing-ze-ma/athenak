# NOTE 2026-10-16 raven -> caltech: accretor donor light, proposed fix (USER decision) + STOP-rule ack

## Donor light (follows the v3 addendum d6fabde0)
The donor_irr=false diagnostic showed that the donor light, not the stream, drives the thin-atmosphere flow
(v_phi up to 940 km/s at r 9.1-9.6) and the dt drop (7.0e-6 -> 2.4-3.3e-6). In the 2-D band, the irradiation
at the star is ~30x the 3-D expectation (E ~5x at all phi, azimuthal F ~ cE).

USER decision (10-16, via Raven): adopt the addendum's fixes 1 + 2. Caltech owns the code, so please implement:
1. **No irradiation from the back-sink sector.** It doubles Eirr in the stream phase (1.55e7) and is not needed.
2. **A diluted, directed donor source instead of the lobe holding E = a T_d^4.** Raven's preference is the
   physical form: an inward flux of the 3-D value L_d/(4 pi d^2) (d from the donor centre), entering only within
   the donor's solid angle (through r_out or the lobe surface).
   - Simplest acceptable stopgap: scale T_d so that the band flux at the star's sub-donor point equals the 3-D
     value. Here that is ~1/30 of the present flux, i.e. T_d,eff ~ 26 kK x 30^(-1/4) ~ 11 kK.
   - Please say which form you use.

Open but not decided: fix 3 (check vet_col with a lateral beam: vet_col builds the Eddington tensor from radial
columns, so light arriving from the side is not represented) and fix 4 (a force/sponge guard for atmosphere gas).
Raven can run fix 3 as a short test if you want it.

The v3 gpu job 31212975 (impact) is still queued unchanged. Its impact checks will be added to
NOTE-2026-10-16-raven-accretor-bigdom.md.

## 3-D
Raven will wait for the 3-D physics changes (3-D gravity, Coriolis, IC, WB reference, ambient; theta walls;
back-sink row artefact) before any longer 3-D run. The half-resolution benchmark costs ~15 A100 node-h per orbit
(NOTE-2026-10-16-raven-accretor-3d.md).

## Ack d7b8e051 (link STOP rule)
Raven bsgnc_raven.sh / e1chain_raven.sh and Viper link_viper.sh / e1chain_viper.sh already stop only on the summary
count NON-CONVERGED=[1-9] (+ FATAL / "Picard solve DIVERGED" / NaN / rc != 0). DT_RETRY reasons never contain
"DIVERGED", so no change was needed. The rule is recorded on Raven.
