---
name: general-eos-project
description: "Status and design decisions for the general (non-ideal) EOS project in AthenaK's Newtonian hydro/MHD"
metadata: 
  node_type: memory
  type: project
  originSessionId: d0224e66-3587-4359-a545-9608682177c3
  modified: 2026-08-13T15:45:20.699Z
---

Ongoing project (started 2026-08-07): replace the ideal-gas assumption in AthenaK's
non-relativistic Hydro and MHD with a general EOS. Target physics is the **optically thick
interior with radiation locked to the gas (LTE), plus H2 dissociation and H/He
ionization** — i.e. `e(rho,T)` = translational + H2 rot/vib + dissociation + ionization +
`aT^4`, and `p = rho k T/(mu(rho,T) m_u) + aT^4/3`.

**Branch:** `general-eos`, off `wb_allmach_cubedsphere`. **Pushed to the fork**; the tip is
`6dddbb4a` as of 2026-08-13 and local == `origin/general-eos`. No PR opened, upstream
untouched. First pushed 2026-08-12, when all 19 commits were author-rewritten to
`Jingze Ma <jingzema0211@gmail.com>` — so the SHAs in the list below are the PRE-rewrite
ones and NO LONGER EXIST. See [[athenak-fork-git-setup]]. The five commits since Stage 3
are real and current:
- `03febafa` Stage 3 — the tabulated analytic EOS ([[general-eos-stage3-table]])
- `4a59a1c9` built-in test pgens through the EOS + `general_eos=table` tests
  ([[general-eos-table-linear-wave]])
- `30cb3cbc` linwave test scripts' stale `vflow` ([[branch-preexisting-breakage]])
- `c94504f7` the last `TODO(stage3)`s ([[general-eos-stage3-loose-ends]])
- `6dddbb4a` `docs/general_eos.md`

Historical commit list (dead SHAs), oldest first:
- `dc36ca57` Stage 1 — general EOS machinery
- `316b42dc` Stage 2 — HLLC/HLLCLM/LHLLC/AUSM+up/HLLD
- `2aec39fb` Well-balanced scheme refuses `eos = general` (SUPERSEDED, see below)
- `5cd44128` Thermal conduction (reads cached `wtemp`; `c_v` from the EOS in NewTimeStep)
- `2d0a8185` Well-balanced scheme (static + dynamic) supports a general EOS
- `9f661639` `deep_hot_jupiter_rt` routed through `getWBq0`
- `1dfd21f1` Source terms
- `7180bc34` Coarse prolongation buffers
- `1e432b81` Fix AMR segfault (null `Mesh::ptree`) — not an EOS bug
- `2a3685fd` Restart inserts blocks into their own panel's tree
- `84b2c9b7` `deep_hot_jupiter_rt` converted to the general EOS (2026-08-09)
- `30d21859` + `4cfdc329` two `deep_hot_jupiter_rt` restart bugs (not EOS-related)
- `fe320594` regression tests, `tst/scripts/{hydro,mhd}/*_general_eos.py`
- `7a9c50dd` ISM cooling number density from `MeanMolecularWeight`
- `6f54cc79` shared `src/pgen/pgen_eos_utils.hpp`; WB sources routed through `getWBq0`
- `b63c09e7` `convection` converted
- `853857e2` solar/cooling convection atmospheres built with the EOS

**`tfloor_kelvin` exists as of 2026-08-16 (commit `c8657f34`).** `tfloor` floors the CODE
temperature, which is `p/d` for an ideal gas (kelvin scale includes `<units>/mu`) but the
EOS's own T under `eos = general` (that scale with mu divided back out) -- the two differ
by exactly mu, so a `tfloor` copied between ideal and general inputs is silently wrong.
`<hydro|mhd>/tfloor_kelvin` states it in kelvin and converts at startup. Additive: `tfloor`
is unchanged. Setting both is fatal, as is using it without `<units>`, under `isothermal`,
or under SR/GR. Test: `tst/scripts/hydro/hydro_tfloor_kelvin.py` (its input sets
`<units>/mu = 2.32408` on purpose so the mu factor does not cancel). Docs in
`docs/general_eos.md`.

**Key design decisions (settled with the user, do not relitigate):**
1. General path added *alongside* the ideal path, not replacing it. `nder = 0` for an ideal
   gas so existing runs allocate nothing extra and take the original code path.
2. Pressure and `Gamma_1` are evaluated ONCE per cell in ConsToPrim into `wder`
   (`IDPR`,`IDG1` in `DerivedIndex`), then reconstructed to interfaces; Riemann solvers
   never call the EOS. `wtemp` caches T as the warm start for the `T(rho,e)` root find.
3. EOS works internally in cgs; `EOS_Data` carries `dens_cgs`/`pres_cgs`/`temp_cgs`.
   Convention is **mu_ref = 1**, so composition lives in the EOS (`p = rho T/mu(rho,e)`),
   not in the unit system. `<units>` block is REQUIRED for `eos = general`.
4. `Gamma_1` must come from the identity `chi_rho + p chi_T^2/(rho T c_v)` — not guessed.
5. Radiation pressure needs no interface change, but makes `T(rho,e)` a quartic root find
   on its own, independent of ionization.

**Done:** C2P/P2C, reconstruction, LLF, HLLE, HLLC, HLLCLM, LHLLC, AUSM+up, HLLD, CFL
timestep, FOFC, reconstruction energy floor (now density dependent), thermal conduction.

**Deliberately refused rather than faked:** the **Roe solver** exits with a FATAL ERROR
under `eos = general` (its Roe average is only defined for an ideal gas; substituting
Gamma_1 breaks `F_R-F_L = A(U_R-U_L)`).

**The well-balanced scheme was ALSO refused at first, and that was wrong** — corrected in
`2d0a8185`. The refusal assumed exact balance requires an analytically integrable
background. It does not: the gravity source is discretised as the background's own
pressure difference across the cell (see `GravitySource` in `hse_atm.cpp`), not `-rho g dx`,
so balance follows from the reconstruction and the source calling the SAME background
routine — now `getWBq0()`, dispatching to `getWBerho` (ideal) or `utils/wb_background.hpp`
(general). Key points, do not relitigate:
- isodensity/isothermal have exact general-EOS forms; isentropic steps `d^(Gamma_1-1)`
  with Gamma_1 of the *coefficient cell*, then projects `e` onto `h + phi = const`.
- Gamma_1 as a LOCAL half-cell coefficient is legitimate (it plays the adiabat constant's
  role); what fails is substituting it into the ideal-gas identity `p = (gamma-1)e`.
- The background walks TWO segments like `getWBerho` — the neighbour supplies the
  coefficient for its own half cell. Integrating everything from the anchor is worse.
- Under a general EOS the solvers read pressure from `wder`, so the pressure channel needs
  deviation reconstruction too, or the whole hydrostatic gradient returns.

## RESUME HERE (last worked 2026-08-13)

**THE WHOLE FEATURE IS IMPLEMENTED, VERIFIED, DOCUMENTED AND PUSHED.** Nothing is
blocking, nothing is half-finished, and the user asked (2026-08-13) to pick the next
direction up in conversation rather than leaving a queued task. So START BY ASKING what
they want to do next; the candidates are in "Next" below, and none is urgent.

**STAGE 3 IS IMPLEMENTED AND VERIFIED** — the tabulated analytic EOS (H2 dissociation,
H/He Saha ionization, optional radiation pressure) is committed as `03febafa`. All ten
`TODO(stage3)` branches in `src/eos/eos.hpp` are filled in. See
[[general-eos-stage3-table]] for the design, the two nasty bugs, and the full
verification.

The composition-model question that used to block this is SETTLED: the user chose a
tabulated `(log rho, log T)` EOS, radiation behind a runtime switch, and a retained
gamma-law mode. `<hydro|mhd>/general_eos` = `gamma` (default) or `table`.

Because the default is `gamma`, everything previously verified is UNCHANGED: all 16
eostest cases still reproduce `eos_compare_baseline.txt` byte for byte and `dhj` is
bitwise identical.

**Items 1 and 2 are DONE (2026-08-13, commit `4a59a1c9`).** `shock_tube` and `linear_wave`
now go through the EOS, `linear_wave`'s hydro AND MHD eigensystems take `Gamma_1` instead
of `gamma`, and two new convergence tests cover `general_eos=table`. Details in
[[general-eos-table-linear-wave]].

**The leftover `TODO(stage3)`s are also DONE** (2026-08-13, commit `c94504f7`) — entropy
floor, the `EnergyFromPressure` re-solve, and `wb_background.hpp`'s five temperature
solves. `grep -rn "TODO(stage3)" src/` now returns nothing. See
[[general-eos-stage3-loose-ends]].

**Documentation is written** (`6dddbb4a`): `docs/general_eos.md` is a wiki-ready parameter
reference — the 14 `eos_*` parameters, the two `general_eos` modes, the required
`<units>` block, the `mu_ref = 1` convention, what is refused and why, the cost model, and
the composition model's validity limits. It is STAGED in `docs/`, not in the wiki: the
project keeps docs in the separate `IAS-Astrophysics/athenak.wiki` repo, which the fork
cannot push to. Move it there and delete the file when a PR is opened.

**Next, in rough priority order:**
1. Nothing is blocking. The remaining work is only for upstreaming: move the doc into the
   wiki, and clean up the pre-existing cpplint over-length lines in the big pgens.
2. Two PRE-EXISTING breakages found while running the harness, neither EOS related — see
   [[branch-preexisting-breakage]]. One is now fixed.

Still true and still worth knowing:
- `solar_convection` and `cooling_convection` HARDCODE `Rgas = 1.38e8` (mu ~ 0.60) in four
  places each. Under `general_eos=table` that value is ignored and the atmosphere changes
  — correct, but it will look alarming.
- `src/pgen/deep_hot_jupiter.cpp` DOES NOT COMPILE at HEAD (`Mesh::use_grid_stretch`,
  since renamed). Nothing to do with the EOS work; dead code sitting in the tree.

**What Stage 3 meant (HISTORICAL — now done, see [[general-eos-stage3-table]]).**
Stage 3 means filling in the `TODO(stage3)` branches in `src/eos/eos.hpp` — now TEN, not
eight: `Pressure`, `Gamma1`, `EnergyFromPressure`, `Temperature`, `EnergyFromTemperature`,
`SpecificHeatCv`, `ChiRho`, `ChiT`, plus `MeanMolecularWeight` and
`DensityFromPressureTemperature` added by `84b2c9b7` (`Enthalpy` is derived and needs
nothing). `ChiRho`/`ChiT` were added for the well-balanced scheme and are also what
`Gamma_1` is built from, so not optional extras. `DensityFromPressureTemperature` is the
(p,T) -> d inversion — the one direction the (d,e) interface cannot evaluate — and is
needed by any pgen that builds a hydrostatic background by marching in (p,T); it is a root
find on d, used at setup only.

Structural point to start from: the physics is naturally a function of `(rho,T)`, but the
whole code calls the interface with `(rho,e)`. So the core is the root find `T(rho,e)`, with
`wtemp` as the warm start (design decision 2), and every other accessor is cheap once `T` is
known. Radiation pressure alone makes that root find quartic (decision 5). A sensible order:
`Temperature` first, then everything else on top of it, then `ChiRho`/`ChiT` by
differentiating the same expressions, then `Gamma1` from the identity in decision 4.

(HISTORICAL, now settled: the composition-model question — Saha per cell vs a tabulated
ionization fraction — was decided in favour of the `(log rho, log T)` table. Do not
reopen it.)

One thing that was expected to bite when the gamma-law placeholder went away: the whole
verification method below would disappear. It did NOT, because the gamma law was retained
as `general_eos = gamma` and is the default, so every ideal-vs-general comparison below is
still live. **The regression tests
now capture it** (`fe320594`) — `tst/scripts/hydro/hydro_general_eos.py` and
`tst/scripts/mhd/mhd_general_eos.py` run the same linear wave with `eos=ideal` and
`eos=general` and require the L1 errors to match, over plm/wenoz x llf/hllc(hlld) x
uniform/SMR, plus an assertion that `roe` aborts. They need their own input files because
`eos=general` requires a `<units>` block and the command line cannot create a block that
does not exist. **RESOLVED: the first option was taken** — the gamma-law mode is
selectable and is the default, so those equality checks stay meaningful, and the table
mode got its own convergence tests instead
(`*_general_eos_table.py`, [[general-eos-table-linear-wave]]).

**`deep_hot_jupiter_rt` is DONE** (`84b2c9b7`) — it was the last pgen still ideal-gas in
its *physics*. Converted: the picket-fence IC (grad_ad now `(Gamma_3-1)/Gamma_1` keeping
the pgen's 0.9 factor), both hydrostatic integrators (now `dln p/dz = rho g/p` closed with
the EOS inversion), the two-stream RT's pressure and Kelvin temperature, the dead implicit
RT solves, and `SetResistivity` in `src/diffusion/` (which had its own ideal-gas T and a
fixed mu). Under a general EOS `<problem>/Rgas` is unused and the pgen warns if it
disagrees with the units block. An ISM-cooling run would still need its number density
rethought (see source terms below) — `MeanMolecularWeight` now exists for exactly that.

**Mesh refinement is done** (`7180bc34`): `ConsToPrimCoarseBndry` in
`bvals/prolong_prims.cpp` called the ideal-gas inversion unconditionally — the floors need
the EOS even though `e = E - KE` does not. There is deliberately no `coarse_wder`:
prolongation produces only primitives, converted straight back to conserved and re-inverted
over the whole block. Verified with STATIC refinement (3 levels, 40 blocks, both
`prolong_primitives` settings, hydro + MHD).

**AMR is verified too** (`1e432b81` first had to fix a null `Mesh::ptree` that made AMR
segfault for an ideal gas — see [[amr-broken-on-branch]]). With the general EOS evaluating a
gamma law, a 3-level adaptive run produces an identical refinement history and identical
final grid to the ideal path, agreeing to 9.1e-15.

**Source terms are done** (`1dfd21f1`): both cooling terms read the cached `wtemp`, and the
code-to-Kelvin factor is `eos_data.temp_cgs`, NOT `punit->temperature_cgs()` — the latter
folds in the fixed `<units>/mu` while the EOS uses `mu_ref = 1`. They agree only when
`mu = 1`. The turbulence driver needed nothing (its non-relativistic path touches conserved
density/momentum only; all `eos.gamma` uses sit inside `flag_relativistic`), and
`ConstantAccel` gates on `is_ideal`, which the general EOS sets true. Left alone: the ISM
cooling number density still uses the fixed `<units>/mu`, so a varying mean molecular weight
is not reflected.

**Dead code, do not be fooled:** `AddIsotropicHeatFluxSpitzerCond` is entirely commented
out upstream (`conduction.cpp` lines 189-339), so the whole Spitzer conduction feature is
non-functional and untestable. Its general-EOS edits are consistent but unexercised.
`NewTimeStep` also has a pre-existing upstream shadowing bug: the `if (spitzer)` block
re-declares `Real temp_unit`/`kappa_unit` locally, leaving the outer ones at 0.0.

**Verification method:** the general EOS currently evaluates a gamma-law, so it must
reproduce the ideal path. **CAUTION — the "agree exactly" claims below are an ARTIFACT of
the test inputs writing `.tab` at `%12.5e`, which hides anything below ~1e-5 relative. At
`%24.16e` the real gap is ~1e-15 absolute, pre-existing and expected (the general path
reconstructs p through PLM, the ideal path recomputes it from the reconstructed e). See
[[general-eos-stage3-perf-cleanup]].** With that caveat: LLF/HLLC/HLLCLM/HLLD agree exactly
on Sod and Orszag-Tang;
LHLLC/AUSM+up agree to ~1e-15 (Gamma_1 is reconstructed, not a compile-time constant);
HLLE differs at truncation level across shocks because the general path uses a
Davis/Einfeldt wave-speed bound instead of a Roe average. **This test disappears once
Stage 3 replaces the gamma-law**, so close remaining leaks before then.

**Test inputs** all live in `/orion/u/jinma/ATHENAK/eostest/` (NOT in the repo, NOT in
[[run-directory-untouchable]]), and each has an `_ideal`/`_gen` pair run as
`athena -i X_gen.athinput -d DIR hydro/eos=general`:
- `wb_iso_*` hydrostatic atmosphere, well-balanced (`hydro/wb_option=` picks the option)
- `wbs_plain_*` / `wbs_perturb_*` static well-balanced
- `cool_*` / `coolm_*` Sod with ISM cooling, hydro / MHD (units scaled so cooling actually
  bites — the obvious unit choice makes it 13 orders too weak to test anything)
- `smr_*` / `smrm_*` static refinement, hydro / MHD; `amr_*` adaptive
- `dhj_*` reduced deep-hot-Jupiter, spherical polar
- `cond_*` thermal conduction
Binary output is compared with `run/bin_convert.py` (`read_binary`, then `mb_data[var]`).
Build is configured with `PROBLEM=built_in_pgens`; the user said not to bother restoring
`cooling_convection`.
See [[athenak-fork-git-setup]] and [[freya-build-procedure]].


## Other problem generators (2026-08-09)

`src/pgen/pgen_eos_utils.hpp` now holds the shared EOS-aware conversions
(`EintFromP`, `PresFromEint`, `TempKelvin`, `DensFromPT`, `GradAd`,
`EintFromDensT`). Each keeps the ORIGINAL ideal-gas arithmetic on the ideal branch so
ideal runs stay bitwise; these are well-balanced problems where the velocity is the
residual of a near cancellation, so last bits show.

| pgen | state |
|---|---|
| `deep_hot_jupiter_rt` | fully converted (`84b2c9b7`), uses the shared header |
| `hse_atm` | converted; verified 1.2e-15 (isothermal) ... 7.3e-7 (isentropic, the documented approximation) |
| `convection` | converted; verified round-off with etotgrav + wellbalance_kappeli |
| `hotbubble`, `cooling_convection`, `solar_convection` | WB source routed through `getWBq0`; the p<->e conversions and the analytic atmospheres are NOT yet done |
| `deep_hot_jupiter` | **SKIPPED at the user's request.** Also does NOT COMPILE at HEAD — references `Mesh::use_grid_stretch`, since renamed to `use_grid_stretch_r` |
| `shallow_hot_jupiter` | **SKIPPED at the user's request** |
| `solid_body_rot`, `sph_mhd_balance` | no EOS dependence, nothing to do |

**The live bug that motivated this:** `hotbubble`, `cooling_convection` and
`solar_convection` called `getWBerho` (the ideal-only background) in their gravity source
while the reconstruction under `eos = general` used the general background. Well balancing
only holds when both come from the SAME background. Fixed in `6f54cc79` — the same
conversion `9f661639` applied to `deep_hot_jupiter_rt`.

**The atmospheres are mostly DONE (`853857e2`), and most needed nothing.** The key
distinction, worth remembering:

- An atmosphere that specifies BOTH `p` and `rho` (closed form, or by integrating
  `dln p = rho g/p dz` with a PRESCRIBED stratification index) satisfies hydrostatic
  balance **independently of the EOS**. Only `e = e(rho,p)` is EOS dependent, and
  `EintFromP` already covers it. This is `hse_atm` (isothermal and nu-polytrope, verified
  algebraically) and `convection` (prescribed `get_gamz`). **Nothing to do.**
- An atmosphere specified in (p,T), or by a thermodynamic adiabat, IS EOS dependent.
  `solar_convection` / `cooling_convection` integrate top down in (p,T): converted, with
  `rho = p/(Rgas T)` -> `DensFromPT` and the Schwarzschild threshold `(gamma-1)/gamma` ->
  `GradAd`. solar_convection also stepped its isentrope as ONE power law from the
  convection onset, only valid for constant nabla_ad; the general path now integrates
  `dlnT = nabla_ad dlnp` step by step (these telescope to the same product for a gamma
  law). cooling_convection's interior gradient is the prescribed `<problem>/grad_int`, a
  specification, so only its threshold needed converting.

**`hotbubble`'s atmosphere: SKIPPED at the user's request** (2026-08-09), along with
`deep_hot_jupiter` and `shallow_hot_jupiter`. Its density is the ideal adiabat
`rho = (p/A)^(1/gamma)` with the bubble carried as a perturbation to `A`; generalising it
would need a decision about how to express an entropy perturbation under an EOS with no
closed-form adiabat (probably `h + Phi = const` plus a chosen entropy parameterisation).
Its well-balanced SOURCE was still fixed in `6f54cc79`, so it is consistent — only the
initial atmosphere remains ideal-gas.

**The pgen sweep is otherwise COMPLETE.**

**Verification trap for solar/cooling convection:** both HARDCODE `Rgas = 1.38e8` (mu ~
0.60, ionized H) in FOUR places each. The gamma-law placeholder has mu = 1, so a direct
ideal-vs-general comparison shows ORDER-UNITY differences — correct behaviour, not a bug,
since under a general EOS composition comes from the EOS and Rgas is unused. To verify,
temporarily set those to `83144621.4563013` (k_B/m_u); the difference then drops to ~5e-7,
which is the documented limiter sensitivity.

**Do not trust per-variable relative metrics on near-zero quantities.** This bit me FOUR
times in one session (bcc3, vely in a static atmosphere, hotbubble velocities at Mach
1e-4, convection velocities). Normalise velocities by the sound speed and fields by
max|B|, and always check the ABSOLUTE difference too.
