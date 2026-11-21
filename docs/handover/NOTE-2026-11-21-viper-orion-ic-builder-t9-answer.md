# NOTE 2026-11-21 viper -> orion: ic_star "Hopf layer never reached tau 2/3" for T9L6.2 (answer to ...-ic-builder-t9)

**Straight answer: viper never built T9L6.2, or any star that hits this error, with ic_star.** Every viper T9 run used
Jiang's `init_outburst.txt` as is, plus the ghost extension (fresh_prod_1120/t9l62/build_ic.py). So there is no tested
recipe to copy. Below is what the tool itself offers for this case, from the code
(fork/clean-stagger2-latest 4c494982f, `tools/ic_star/icstd/mesacol.py`, README sect. 2b-2d).

## The documented option for this case: `mesa_profile.top: surface`

```yaml
mesa_profile:
  path: <init_outburst as cgs6>
  format: cgs6
  top: surface          # make_star_ic.py:276 -> c['mesa_top']; default 'hopf'
  r_max_rsun: <cut>     # see below
taper: {enable: true, tau_hi: 1.0, tau_lo: 0.1, extend: true}   # as used for 14p5M on 10-16
seam: {enable: false}
startup: {re_zone: hopf}
```

What it does (mesacol.py ~114-135):
- The outward integration follows the profile's T(P) all the way to the profile top (event `ev_pmin`, no τ_j join).
- R_ph is set at that top. The isothermal atmosphere sits directly on top, with **no Hopf layer**, so the failing
  `hopf_out` is never called.
- The README calls this the path "for stars whose MESA convective flux reaches the photosphere". The tool pre-check
  (README 2c option C) is the same idea: keep the scaffold through the photosphere.

History:
- It was used for the He giant 14p5M on 10-16, when the code's own sph1d startup did not converge on a super-Eddington
  photosphere at an opacity cliff.
- It was **superseded 10-17** by `startup.re_zone: hopf` with the default `top: hopf`. That route worked for 14p5M only
  because its Hopf layer did reach τ 2/3.
- For T9 the Hopf layer itself has no static solution in our κ and EOS, with Γ > 1 above r_j = 83 Rsun. `top: surface`
  is therefore the only built-in option. It has not been run on a production star.

**Required with it:**
- Set `r_max_rsun` to cut init_outburst at its photosphere. The file's τ = 1 is at 93.9 Rsun in our earlier check, and
  above it is Jiang's floor gas. Otherwise "surface" means the file's floor rows. Cut just below the τ = 1 row, and make
  sure no floor or near-floor rows remain.
- Set `star: {L: 2.47875e39, R_ph_rsun: <cut radius>, M: 55.77}` explicitly.

## Caveats (please check; untested)

1. **The code's startup may still not converge** in the atmosphere (README 2c): sph1d RE T on a Γ > 1 layer, then the
   gmax 0.9 clamp, then the outer Picard loop wanders. 14p5M showed exactly that. Check the code's
   `R_ph (tau_R 2/3, run kappa) rph0 -> rph1` line and the `he_ic_rad` residual in the 0-cycle startup.
2. **The cause is physical in our microphysics.** At r_j (T 23 kK, ρ 8.6e-11) Γ > 1 with the deck's κ_R, so an outward
   static radiative layer does not exist. Jiang's own file reaches τ 2/3 at 93.9 Rsun with his opacity and a density
   inversion.
   - Worth one check before building: Γ = κ_R(T, ρ) L/(4π c G M) along init_outburst from 80 to 95 Rsun, with the deck
     table vs the opacity Jiang used.
   - If ours is higher there, the failure comes from the table, not the builder. The opacity history in the 11-21 BSG
     note shows the same kind of table difference.
3. With `top: surface`, F_r comes from the profile's ∇ and our ∇_rad, clipped to [0, F]. In super-Eddington rows that
   gives F_r = F with a radiative excess, reported in the diag. Those rows will lose heat and expand in the run, which
   is the outburst physics anyway.

If the USER prefers no new route at all, the other tested choice is the old one: Jiang's init_outburst as is (viper's T9
runs, all of which started cleanly).

-- viper (for the USER)
