# Ma, Bildsten & Jiang 2026 (arXiv 2609.28656): paper vs Linhao's `bsg_star.cpp`

Sources: paper text `/viper/ptmp2/jinma/bsg_1001/2609.28656.txt` (line numbers `txt:N`);
code `/viper/ptmp2/jinma/linhao_1013/bsg_star.cpp` (1034 lines, mtime Oct 13; `cpp:N`);
his files `/viper/ptmp2/jinma/bsg_1001/ma2026_files/BSG/`. No athinput from him exists.
Read-only comparison, 11-04.

## Key differences (most consequential first)

1. **Density seed not in the paper.** Code perturbs the initial density by +-2.5 % white noise in
   every active cell with rho > 0.01 (code units): `if(rho > 0.01 ) amp = 0.05;` (cpp:457) and
   `phydro->u(IDN,k,j,i) = rho * (1.0 + amp * ((double)rand()/(double)RAND_MAX-0.5));` (cpp:463),
   `std::srand(gid);` (cpp:384). Pressure is built from the UNperturbed rho
   (`phydro->u(IEN,k,j,i) = tem * rho/(gamma-1.0);`, cpp:469), so the seed is isobaric (T_gas
   perturbed by -+2.5 %); radiation is not perturbed. The paper only says "We set the initial velocity
   field to zero at every grid point" (txt:329) and mentions no perturbation. Our earlier notes
   (TRUEREPRO.md row 7, MA2026_REPRO_IC.md "no seed") are wrong on this. The cgs threshold is
   0.01 x density_unit (athinput; 1e-11 g/cm^3 only if density_unit = 1e-9). The noise pattern
   depends on gid, i.e. on the meshblock layout.
2. **Atmosphere above the profile: gravity grows outward.** Code: `Real grav_local = grav_rmax *
   pow(x1/ini_profile_max_r,2.0);` and `rho = rho_rmax * exp(-grav_local*(x1-ini_profile_max_r)/(tem_rmax));`
   (cpp:429-430), floored at dfloor (cpp:431). Paper: "rho = rho_photo exp(-(r_grid - r_photo) g / P_photo)
   ... where g is the local gravity at r_grid" (txt:325-327). Differences: code g ~ (r/rmax)^+2, a
   physical local g ~ (r/rmax)^-2; paper formula lacks rho_photo in g/P (typo; code's g/T is the
   correct gas-only rho g/P in code units). Gas pressure only, no radiation pressure, in both.
   Scale at the top (computed from the last profile row, T 63364 K, M 19.99 Msun, mu 0.62):
   H ~ 0.6 Rsun, so the density is ~e^-20 below rho_rmax within ~10 Rsun either way; with
   (r/rmax)^2 the fall-off is steeper (a 1e-20 floor is reached near ~62 Rsun instead of > 70 Rsun).
3. **"Photosphere" = last row of the file, not a photosphere.** Code uses `rmax =
   ini_profile(num_of_points_init-1,0)` (cpp:409) as r_photo. That row is 52.49 Rsun, T = 63364 K,
   rho = 2.64e-11, preceded by ONE interval from 48.32 Rsun (T 70253 K, rho 3.45e-9). The paper says
   "photosphere" / T_photo (txt:310, 320-325) and Teff = 15835 K (txt:153). So the isothermal
   atmosphere is at 63 kK, ~4x Teff. Paper footnote 5 (txt:316-318) admits the integration does not
   always converge to the photosphere.
4. **Profile interpolation is linear in r on linear rho, T, F** (cpp:442-451, also ghosts
   cpp:331-341). Across the 48.32 -> 52.49 Rsun interval rho is interpolated linearly between 3.45e-9
   and 2.64e-11 (e.g. ~1.7e-9 at mid-interval vs ~3e-10 log-linear). Paper: "values interpolated at
   r_grid" (txt:321), method not stated.
5. **Outer hydro BC density differs from the paper.** Code, outflow branch (v_r(ie) >= 0):
   `prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie+i-1) * x1*x1/(x1g*x1g);` (cpp:925), i.e. rho ~ r^-2 in the
   ghosts; inflow branch copies rho and sets v_r = 0 (cpp:921-923). Paper: "gas density, pressure,
   and tangential velocity components in the ghost cells are copied" (txt:434-436). Pressure copy
   (cpp:930) and v_r handling (cpp:923,926) agree with the paper.
6. **Gravity is not in the code.** `// EnrollUserExplicitSourceFunction(GravityPotential);` (cpp:141,
   commented; no such function in the file). Point mass must come from the athinput `GM` key as the
   paper says (txt:489-494). `grav_const_dimensionless` (cpp:138) and the profile's enclosed-mass
   column (19.99 Msun at the top) are used only for the atmosphere g (cpp:412-414).
7. **Planck opacity extrapolation T^-3.5 above the table**, clamping elsewhere (cpp:655-660, 686-702,
   741-743); paper silent. Practically irrelevant at high T (table to log T 6.54, profile max 1.6 MK),
   but at low density (rho < 3.0e-12, table edge log rho -11.52) and low T (< 11600 K) opacities are
   held at the edge value; the atmosphere above ~52 Rsun is in that regime.

## Topic-by-topic table

| Topic | Paper (sec / txt line / quote) | Code (cpp:line, what it does) | Verdict |
|---|---|---|---|
| IC profile source | 2.2.2: radiative-only profile integrated from 16.5 Rsun to the "photosphere" (txt:290-318) | reads `./MESA_profile_combined.txt`, 6 cols r,T,rho,m,kappa,Fr in CODE units (cpp:144-175). We have only the `_cgs` version (16.50-52.49 Rsun, F_r = L/4 pi r^2) | same source; unit conversion via athinput units |
| IC inside r_max | "assigning them the values interpolated at r_grid" (txt:320-321) | linear-in-r interpolation of rho, T, F_r (cpp:434-451); below the first point extrapolates linearly | same; method paper silent (see Key 4) |
| IC above r_max | isothermal at T_photo, rho_photo exp(-(r-r_photo) g/P_photo), g local (txt:322-327) | `tem = tem_rmax`, g = g_rmax (r/rmax)^2, rho floored (cpp:427-431); `radflx = 0.0` (cpp:424) | different (g ~ r^+2); paper formula has a units typo |
| IC velocity | zero (txt:329) | `u(IM1..IM3) = 0.0` (cpp:464-466) | same |
| IC gas energy | "density and energy" mapped (txt:319) | `IEN = tem * rho/(gamma-1.0)` (cpp:469): P_gas = rho T, gas only | same |
| IC radiation | two-zone I+/I- from cE = acT^4, F = L/4 pi r^2, F=0 above photosphere (txt:332-357) | `ir = 0.5*(er/coefa +- radflx/coefb)` with er = T^4, coefa = sum_{mu>0} w, coefb = sum_{mu>0} mu w (cpp:476-501): exact discrete-angle moments | same in intent; paper eq 6 normalisation (I+ + I- = acT^4/4pi, I+ = (acT^4/4+F)/2pi) is not consistent with eqs 4-5 in the text extraction; code is consistent |
| Seed / perturbation | none mentioned | 5 % peak-to-peak density noise where rho > 0.01 code, isobaric, seeded by gid (cpp:384, 457-463) | **different (code only)** |
| Inner BC hydro, rho and P | ghosts = interpolated initial profile rho, P (txt:381-383) | `bd_data` from profile at ghost radii (cpp:320-345); `prim(IDN)=bd_data(0)`, `prim(IEN)=bd_data(0)*bd_data(1)` (cpp:797, 801) | same |
| Inner BC velocity | v_r in the boundary-adjacent ghost layer solved so the HLLC contact speed vanishes; v_t = v_p = 0; v = 0 in other ghosts (txt:412-432) | all ghosts v = 0 (cpp:798-800), then iterate `prim(IVX,k,j,is-1)` up to 100 times, re-reconstructing with the run's xorder, PVRS q-factors (cpp:815-861); warning if not converged to `relative_tolerance = 1.0e-20` (cpp:819, 870-875) | same; note the 1e-20 tolerance is below round-off, so convergence by tolerance is rare and the 100-iteration path / warning is likely the norm |
| Inner BC radiation | ghost T from profile, two-zone with F = L/4 pi r^2 (txt:383-390) | `bd_data(2)` = interpolated F_r column, er = T^4, two-zone, then `ComToLab` with ghost velocities (cpp:984-1024) | same (F column = L/4 pi r^2 to 2e-12); comoving->lab transform paper silent |
| Inner BC B-field | n/a (hydro) | copies face fields if MAGNETIC_FIELDS_ENABLED (cpp:880-902) | unused |
| Outer BC hydro | copy rho, P, v_t; v_r copied if > 0, zero if < 0 (txt:433-438) | rho copied (inflow) or ~ r^-2 (outflow); P, v_t copied; v_r copied / 0 (cpp:916-931) | **different (rho in outflow)** |
| Outer BC radiation | copy outgoing, zero incoming (txt:438-443) | same (cpp:947-957) | same |
| theta, phi BCs | periodic all variables (txt:362-364) | none in code | athinput |
| EOS | ideal gas mu = 0.62, gamma = 5/3 (txt:484-487) | P = rho T in code units (cpp:469, 801); mu from `radiation/molecular_weight` (cpp:133), gamma from athinput (cpp:378) | same form; values athinput |
| Units | not stated | time_unit = l_unit / sqrt(k T_unit/(mu m_H)) (cpp:135-136); G in code units (cpp:138); rho, T, l units from `<radiation>` (cpp:130-132); crat/prat from Athena++ | paper silent |
| Gravity | point mass 20 Msun via `GM`, no self-gravity (txt:488-495) | no gravity source in the code (cpp:141 commented) | code silent -> athinput |
| Opacity source | TOPS, X 0.7, Y 0.292, Z 0.008 GS98 scaled (txt:452-466) | reads 34 x 100 log tables kR, kP, ke on log T 4.06-6.54, log rho -11.52-1.00 (cpp:200-266) | same tables |
| Absorption / scattering split | kappa_a = kR - ke; ke = 0.34 or kR - 1e-5 if kR < 0.34001 (txt:467-479) | `kappa_absorption = kappa_ross - kappa_electron;` (cpp:601); ke table already capped (checked: min ke 0.33545, min kR-ke = 1e-5, none <= 0) | same |
| Opacity use | kappa_a, kappa_e, kappa_Planck given to Athena++ (txt:446-451) | `sigma_s` = ke rho, `sigma_a` = (kR-ke) rho, `sigma_p = sigma_pe` = kP rho (cpp:604-607); T = max(P/rho, tfloor) (cpp:587); also in ghost cells (cpp:566-577) | same; sigma_pe = sigma_p and ghost evaluation paper silent |
| Opacity interpolation | "interpolates gas density and temperature" (txt:480-482) | bilinear in log kappa vs (log rho, log T); clamp to edge in rho and low T; kP x 10^(-3.5 dlogT) above T max (cpp:617-745) | paper silent on method and out-of-table handling |
| RT method | implicit DO, 120 angles (txt:160-164, 247) | nothing (Athena++ core + athinput) | athinput |
| CFL / tolerance | CFL 0.3 hydro, 0.025 "radiation", tol 1e-6, < 10 iterations (txt:496-501) | nothing | athinput |
| Domain | 20-80 Rsun (txt:212-220, "80 M_sun" typo), theta pi/3-2pi/3, phi 0-pi/3 (txt:198) | x1min used only for the ghost check (cpp:191, 324) | athinput |
| Grid / SMR | uniform r, theta, phi; 196 x 128 x 128 base; SMR 35-55 Rsun; 2560 blocks of 16^3 (txt:236-246) | nothing | athinput |
| Riemann solver / reconstruction | HLLC (txt:416-417) | boundary uses `precon->xorder` 1/2/3 (cpp:761-768) and an HLLC contact formula | athinput; BC is only exact if the run uses HLLC |
| Floors | "density floor value" (txt:326), value not given | `tfloor = radiation/tfloor`, `rhofloor = hydro/dfloor` (cpp:126-127); used in IC and opacity T | athinput |
| Run length | t_steady 30 d, run to 57.4 d (txt:557-566) | nothing | athinput |
| User outputs | App. A.2 uses f_rad,r "in Athena++" (txt:1128) | 4 user vars = `phydro->rad_source_term(0..3)` (cpp:547-557); needs Linhao's modified Athena++ (member not in the public tree) | consistent; code needs non-public Athena++ |
| Unused paths | -- | hydro-only fallback `crat = 9166.13; prat = 0.0;` (cpp:391-392); commented opacity output (cpp:521-546) | irrelevant |

## Depends on the athinput we do not have

- `radiation/density_unit`, `T_unit`, `length_unit`, `molecular_weight` (cpp:130-133): fix the code-unit
  profile file, the seed threshold rho > 0.01 (cpp:457) in cgs, and the atmosphere g.
- `hydro/gamma`, `hydro/dfloor`, `radiation/tfloor`, any pfloor; `problem/GM` (point mass) value.
- Mesh: r range, nx1/nx2/nx3, meshblock size, SMR `<refinement*>` blocks, theta/phi periodic BC flags
  (inner/outer x1 must be `user`).
- Radiation: nang/angle scheme (120), nfreq (code assumes 1 in places: `ir_cm(n)`, `mu(0,k,j,is-i,n)`
  without ifr offset, cpp:995-1022), implicit solver tolerance and iterations, cfl_rad, Compton/scattering
  switches, reduced speed of light (crat) if any.
- Hydro: Riemann solver (HLLC?), xorder, integrator, cfl_number, time limit, outputs, restart cadence.
- The code-unit profile file `MESA_profile_combined.txt` itself (we have only `_cgs`).
