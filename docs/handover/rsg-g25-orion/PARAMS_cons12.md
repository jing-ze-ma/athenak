# g25 cons12 -> IDORT module: parameter and Config translation

Reference: old IDORT-dev production cons12 (Fritz job 4325157); Orion stand-in files
`Arepo-IDORT/RSGg_handoff/orion/Config_cons.sh`, `param_g25_orion_cons.txt`. New files: `Config.sh` (VORONOI_CELL),
`Config_classic.sh`, `param.txt`, IC from `convert_ic.py`. Design background: `documentation/modules_idort.md`,
design.md sect. 2.1 (which old options were dropped and why).

## Decisions for the user (not chosen silently)

| # | item | cons12 | this setup | alternatives |
|---|---|---|---|---|
| D1 | face optical depth / closure | `IdortTauWeightMode 1`: w = 1/(1 + tau_f + tau_f^2/30), tau_f = alpha (mean kappa)(mean rho) d, alpha 5 | module defaults first: w = 1/(1 + tau_f), tau_f = alpha min(chi_L, chi_R) d, alpha 5 (`IdortTau0 0`, `IdortTauFace 0`) | cons12 closure (commented out in param.txt): `IdortTau0 30`, `IdortTauFace 1` = the same face tau (kappa = kR + kS of each cell); modules_idort.md: the quadratic form leaves the odd-even mode of thick faces undamped (A5) |
| D2 | thick core | f-hybrid (IRT_STAR_CORE_EOS_RAD) inside 0.85 R* (ramp R1 4.99e13 - R2 5.40e13 cm) + OPTB kappa_eff core | option (b): discrete ordinates everywhere, cons12 kappa_eff below R_c through `IdortOpacityModCoreTable` | f-hybrid on branch idort-fhyb (in progress). README "Q1"/"idort-cprime": (b) drives grid-scale non-radial core motions at <= 20 cells per R |
| D3 | refinement | split/merge, target mass `ReferenceGasPartMass`, REFINEMENT_VOLUME_LIMIT with the CSM/photosphere/density/T-front bands, REFINEMENT_AROUND_DM (forced derefinement commented out) | same (Config refinement block, IdortRefine* of idort-g25ref, REFINEMENT_AROUND_DM_NO_FORCED_DEREF). The reduced-resolution smoke IC must run WITHOUT the block (the mass criterion would split the thinned cells back) | for the reduced IC: scale ReferenceGasPartMass and the band volumes by the thinning factors, or no refinement (done) |
| D4 | time steps | local bins (finest 227 s = top/32 at the start) | local (Config.sh: TIMESTEPS_WAKEUP + MHD + HLLD, B = 0); FORCE_EQUAL_TIMESTEPS as fallback | see the smoke-test result in README.md |
| D5 | solver tolerance | Krylov rtol 1e-2, per-cell tol 1e-3 (IdortCellTol), max 15 outer | `IdortTol 1e-8` (module default, per-cell content norm), Krylov rtol 1e-2 | a looser IdortTol for cost; the budget closes for any outcome |
| D6 | caps/floors of the gas | `GasTmaxCap 1e6` (rho < 1e-9), u floor at 531.2 K, dv caps, limiters | none (design: cons exchange only, no caps). `MinEgySpec 1e10` (base Arepo hydro floor) kept as in cons12 | GAS_TMAX_CAP had a role in the hot thin CSM of cons12; watch T_max of the thin gas |
| D7 | opacity evaluation time | after the previous radiation step (lags one hydro step, `IdortOpacityRefresh 0`) | start of each radiation step (post-hydro state) | not reproducible without code (modules_idort.md, matter tables) |

## Matter (EOS and opacities): checks done

- EOS: `IdortEosMode 2` = the cons12 gas EOS (EOS_OPAL + OPAL_EOS_BLEND_IDEAL + OPAL_IDEAL_STITCH_FIX + EOS_NO_RAD;
  same EOS5_data file). Composition: cons12 passed the cell's X (EOS_OPAL Composition[0]); here `EOS_PASSIVE` carries
  `NuclearComposition` and the module uses Composition[0] as well. X is uniform in the IC (0.6538476353670762, spread 9e-14;
  = the core table's X attribute) and stays uniform (no sources), so the table has ONE X node at `IdortEosX
  0.6538476353670762` (lookups require |X - X0| <= 1e-9: the 8-digit 0.65384764 terminates at the first pressure call).
  Metals: the OPAL table's Z (as cons12).
- Table range (defaults): T 20 K - 1e8 K, rho 1e-26 - 1e4 g/cm^3. IC: rho 1e-20 (background, vacuum-layer cells) -
  3.9e-4 (centre), u 4.84e10 (CSM/background, T ~ 5e2 K) - 4.8e14 erg/g: inside. The vacuum boundary cells keep the
  background state (rho 1e-20, T ~ 500 K): in range; they are frozen and not radiation unknowns (incoming I = 0).
- Table cost: one copy per node in MPI shared memory (163 MB per node, ~1.3 s build on 8 tasks for the default grid;
  all tasks of a node fill it), so 560 ranks = 5 copies; see the run log line of the build.
- Opacity: gs98 Z 0.02 X 0.7 tables, single composition as cons12 (the cell X is not used for kappa, as before).
  `IdortOpacityInterp 1` (bilinear) = cons12 (the module default bicubic differs by rms 1-1.5 %, max 30 % at
  log T < 3.5). Off-table (CSM: logR < table) constant continuation with the cons12 inset 1e-6 and floor 1e-5. The
  modifier stage reproduces cons12 `set_opacity()` to 3e-14 (modules_idort.md).

## param.txt: every cons12 parameter

| cons12 | value | rewrite | note |
|---|---|---|---|
| BoxSize | 1.7515e16 | BoxSize (same) | also written to the IC header (0 in the old IC) |
| InitCondFile | rhoin-..._c1s.ic | g25_full / g25_red | convert_ic.py (BoundaryRole added, BoxSize, Density dropped) |
| OutputDir, SnapshotFileBase, OutputListFilename, OutputListOn 1 | | run-local, OutputListOn 0 | output list was for the crash window |
| TimeLimitCPU 88000, ResubmitOn, ResubmitCommand, CpuTimeBetRestartFile 3600 | | same meaning | |
| ICFormat 3, SnapFormat 3, NumFilesPerSnapshot 1, NumFilesWrittenInParallel 0 | | ICFormat/SnapFormat 3, MaxFilesWithConcurrentIO | parameter renamed in this Arepo |
| ComovingIntegrationOn 0, Omega0, OmegaLambda, OmegaBaryon, HubbleParam | | same | |
| TypeOfTimestepCriterion 0, ErrTolIntAccuracy 0.05, CourantFac 0.3 | | same | |
| MaxSizeTimestep 10000, MinSizeTimestep 1e-3 | | same | top step 7275.58 s via TimeMax/2^k |
| TimeBegin 0, TimeMax 1.5258e10 | | TimeMax 1.5258e10/2^16 (smoke) | same integer timeline (top step 7275.58 s) |
| TimeBetSnapshot 1.526e6, TimeOfFirstSnapshot 0, TimeBetStatistics 1e3 | | per test | |
| PeriodicBoundariesOn 0 | | 0 | GRAVITY_NOT_PERIODIC requires 0; the mesh is periodic, the vacuum layer separates the faces |
| ReferenceGasPartMass, TargetGasMassFactor, RefinementCriterion 1, DerefinementCriterion 1, MaxVolumeDiff 5, MinVolume 0, MaxVolume 1.99e42, RefinementCellsPerSoftening 10 | | same | Config refinement block (D3) |
| RelaxBaseFac 1e30, RelaxTime 0 | | dropped | damping off in cons12 anyway (RELAXOBJECT / DAMP_STAR_CORE not compiled) |
| MultipleDomains 6, TopNodeFactor 8, ActivePartFracForNewDomainDecomp 0.2 | | same | |
| ErrTolTheta 0.5, TypeOfOpeningCriterion 1, ErrTolForceAcc 0.001 | | same | |
| DesNumNgb 32, MaxNumNgbDeviation 1 | | same | |
| InitGasTemp 0, MinGasTemp 0, MinimumDensityOnStartUp 1e-20, LimitUBelow* 0 | | same | |
| MinEgySpec 1e10 | | same | base hydro u floor (D6) |
| MaxMemSize 3000 | | 3000 | |
| Unit* 1, GravityConstantInternal 0 | | same (cgs) | `IdortClight`, `IdortRadConst` -1 = physical |
| Softening* (types 0-5), SofteningTypeOfPartType*, GasSoftFactor 2, MinimumComovingHydroSoftening 8.26e9, AdaptiveHydroSofteningSpacing 1.2 | | same | point mass softening 6.29e11 cm |
| CoolingOn 0, StarformationOn 0 | | same | |
| CellMaxAngleFactor 2.25, CellShapingSpeed 0.5 | | same | |
| EosOpalTable | EOS5_data | IdortEosOpalTable (same file) | + IdortEosX (see Matter) |
| OpacityRossTable, OpacityPlanckTable | gs98 | IdortOpacityRossTable, IdortOpacityPlanckTable | + IdortOpacityMode 2, IdortOpacityInterp 1, IdortOpacityEdge 1 |
| OpacityFloorTmin 3200 | | IdortOpacityModFloorTmin 3200 | + IdortOpacityModifiers 1, KappaFloor 1e-5, EdgeInset 1e-6 (old compiled constants) |
| IrtStarDustMinRadius 1.174e14 | | IdortOpacityModDustRadius | Tlow/Thigh 1200/1500 (old constants) |
| IrtStarPhKapDensity 1e-10 | | IdortOpacityModPhotRho 1e-10, IdortOpacityModPhot 1 | old IRT_STAR_PH_OPACITY_SWITCH_TO_ROSS |
| IrtStarCoreBTable | coreB_side_g25_c1s.hdf5 | IdortOpacityModCoreTable (kappa_eff) | |
| IrtStarCoreRadius 1.761e12, IrtStarCoreRadiusFactor 1 | | IdortOpacityModCoreRadius 1.761e12 | |
| IrtStarCoreLuminosity 6.855e38 | | SphericalSourceLuminosity 6.855e38 | |
| IrtStarCoreHeatRadius 0 (table r_h) | | SphericalSourceRadius 3.522178735e11 (set explicitly, = table attribute r_h = 0.2 R_c) | same cubic-spline kernel, S = sum W V over all cells, active cells over their own step |
| IrtStarCoreRadiusDampFactor 2 | | dropped | damping off |
| IrtStarMin/MaxRadiusCsm, IrtStarVolumeThreshCsm | | IdortRefineCsmRmin / CsmRmax / CsmVmax | idort-g25ref |
| IrtStarMin/MaxRadiusPh, IrtStarVolumeThreshPh | | IdortRefinePhRmin / PhRmax / PhVmax | |
| IrtStarMin/MaxDensityPh, IrtStarVolumeThreshPhDensity | | IdortRefineRhoMin / RhoMax / RhoVmax | |
| IrtStarMin/MaxTempPh, IrtStarVolumeThreshPhTemp, IrtStarMinDensityPhTemp | | IdortRefineTMin / TMax / TVmax / TRhoMin | |
| IrtStarMin/MaxRadiusPhHr 0, IrtStarVolumeThreshPhHr | | dropped | band 0..0 = off in cons12 |
| (IDORT_SHELL_FLUX_LOG radii: 0.05 R_c .. 5 x 4.66e13, 48, + R_c) | | IdortShellFluxEvery 1, N 48, Rmin 8.805e10, Rmax 2.330e14, Rextra 1.761e12 | idort-g25ref; `idort_shell_flux.txt`, same columns as cons12 |
| IrtStarCoreEosRadR1 4.99e13, R2 5.40e13 | | dropped (D2) | f-hybrid, idort-fhyb |
| GravPmEnergyMode 2 | | GRAVITY_PM_ENERGY_CONSERVING (Config) | point-mass energy correction ported from the old code (idort-gpm) |
| IdortAlpha 5 | | IdortAlpha 5 | |
| IdortTauWeightMode 1, IdortTau0 30 | | IdortTau0 0 (D1; cons12: IdortTau0 30 + IdortTauFace 1) | |
| IdortBrulsNum 2 | | IDORT_NUM_ANGLES=24 (Config) | |
| IdortMinIterNum 5, IdortMaxIterNum 10, IdortTolerance 1e-8 | | IdortMaxIter (20), IdortTol 1e-8 (D5) | different solver |
| IdortKrylov* (Mode, Restart 40, MaxIter 400, MaxOuter 15, Tol 1e-2, Stage, LinT, Poly, Prec, Comm, EW, TolMax, PrecF, Check, Orth, FastSweep, CommF, BasisF, CvNum, NegMode) | | IdortKrylovRestart (100), IdortKrylovMaxIter (2000), IdortKrylovRtol (1e-2), IdortPrecSweeps, IdortPrecSweepMode, IdortGmresOrth; IdortDsa 1 | one Newton-FGMRES with block GS + DSA (design 2.1) |
| IdortFailRetry, IdortFailSkip, IdortFailFracMax | | dropped | no fallback chain; conservative skip only (budget column) |
| IdortCellTol 1e-3, IdortCellFracTol 1e-3 | | IdortTol, IdortTolFrac (0) | |
| IdortSafeAlpha, IdortSafeBeta, IdortSafeMargin, IdortSafeNoZero, IdortSafeRadTakeMax, IdortSafeUFloorT 531.2, IdortDvCapCs, IdortDvCapTau, IdortDvCapCsThin, IdortCoupleMode, IdortCoupleGasMax, IdortMinEgySpec | | dropped (D6) | limiters/caps/floors of the old coupling; the module's exchange is conservative by construction (IdortCoupleGas 1) |
| GasTmaxCap 1e6, GasTmaxCapRhoMax 1e-9 | | dropped (D6) | |
| IdortVelMode 5 | | dropped | velocity = start-of-step cell velocity (design 2.1; no half kick needed, README I4) |
| IdortAdFluxTauMin, IdortAdFluxFac, IdortAdFluxMode, IdortFaceFixTauMin | | dropped | ADF / face fix not needed with the conservative face form (design 2.1, README idort-cprime) |
| IdortConsPosFix 4 | | IdortPosRounds 4 (+ IdortPosMode 2, IdortPosClearEvery -1) | positivity policy of the module |
| IdortPostHydroState 1, IdortOpacityRefresh 0 | | by construction / D7 | radiation step after the hydro step |
| IdortRestartOldNumAngles, IDORT_ANGLE_REMAP_RESTART | | dropped | no restart from cons12 files |
| IdortDiagKapPCap 0 | | dropped | (0 = off in cons12) |

## Config.sh: cons12 options

Kept (base Arepo): VORONOI, MHD + RIEMANN_HLLD (B = 0; needed for TIMESTEPS_WAKEUP), VARIABLE_GAMMA,
REGULARIZE_MESH_CM_DRIFT(_USE_SOUNDSPEED), REGULARIZE_MESH_FACE_ANGLE, TREE_BASED_TIMESTEPS, MESHRELAX_DENSITY_IN_INPUT,
SELFGRAVITY, HIERARCHICAL_GRAVITY, CELL_CENTER_GRAVITY, GRAVITY_NOT_PERIODIC, ALLOW_DIRECT_SUMMATION,
DIRECT_SUMMATION_THRESHOLD=1000, EXACT_GRAVITY_FOR_PARTICLE_TYPE=1, NO_SELFGRAVITY_TYPE=1, NO_GRAVITY_TYPE=1,
EXACT_GRAVITY_REACTION, EVALPOTENTIAL, NSOFTTYPES=4, MULTIPLE_NODE_SOFTENING, ADAPTIVE_HYDRO_SOFTENING, CHUNKING,
DMLOWESTTIMEBIN, DMFIXED, DOUBLEPRECISION=1, OUTPUT_IN_DOUBLEPRECISION, VORONOI_DYNAMIC_UPDATE, OUTPUT_VOLUME,
OUTPUTPOTENTIAL, HAVE_HDF5, DEBUG, EOS_NSPECIES=6.

Replaced: EOS_OPAL + OPAL_EOS_BLEND_IDEAL + OPAL_IDEAL_STITCH_FIX + EOS_NO_RAD -> EOS_PASSIVE + IdortEosMode 2;
OPACITY_MESA, OPACITY_PLANCK, OPACITY_TABLE_EDGE_CLAMP, IRT_STAR_PH_OPACITY_SWITCH(_TO_ROSS),
IRT_STAR_DUST_OPACITY_SWITCH, IRT_STAR_CORE_OPTB (opacity part) -> IdortOpacityMode 2 + modifiers; IRT_STAR_CORE_OPTB
(heating) -> SPHERICAL_SOURCE; IDORT_OUTFLOW -> BOUNDARY_CELLS vacuum layer; IDORT, USE_BRULS, IDORT_NUM_ANGLES=24 ->
IDORT, IDORT_NUM_ANGLES=24; IDORT_RELAX_INI(+ITER_NUM_INI) -> IdortInitMode 2; IDORT_COUPLE, IDORT_IT_COEVOL,
IDORT_FULL_IMPLICIT, IDORT_CONSERVATIVE_EXCHANGE -> IdortSolveMode 1 + IdortCoupleGas 1; IDORT_KRYLOV(_FAST) -> module
solver + IDORT_DSA; IDORT_RAD_SOUND_SPEED -> IdortRadSignalSpeed 1; IDORT_FRESH_STATE -> by construction.

Kept/replaced (refinement, D3): REFINEMENT_SPLIT_CELLS, REFINEMENT_MERGE_CELLS, REFINEMENT_VOLUME_LIMIT,
REFINEMENT_AROUND_DM (+ REFINEMENT_AROUND_DM_NO_FORCED_DEREF = the old code's commented-out forced derefinement),
IRT_STAR_REFINE_CSM/_PH/_PH_DENSITY/_PH_TEMP -> IdortRefine* (IRT_STAR_REFINE_PH_HR: off in cons12).

Dropped: MHD_POWELL (B = 0), RELAXOBJECT,
RELAX_RUNTIME(_SPECIFY), DAMP_STAR_CORE(_ALL) (off in cons12), IRT_STAR_CORE, IRT_STAR_CORE_EOS_RAD (D2),
IDORT_TAU_WEIGHT (D1: available as IdortTau0 30 + IdortTauFace 1), IDORT_SAFE_COUPLING, IDORT_SAFE_FIX2, IDORT_COUPLE_CONSISTENT, GAS_TMAX_CAP (D6),
IDORT_VEL_CONSISTENT, IDORT_ADIABAT_FLUX, IDORT_FACE_FIX, IDORT_ANGLE_REMAP_RESTART, IDORT_SNAPSHOT_RESTART_INI,
diagnostics (IDORT_CONV_DIAG, IDORT_DRAG_DIAG, IDORT_SHELL_FLUX_LOG -> IdortShellFlux* (idort-g25ref); analyse.py
also gives L(r) from the snapshot flux), GW_SIGNAL,
COMPUTE_POTENTIAL_ENERGY, DOUBLEPRECISION_FFTW (no PM), output-only options (OUTPUT_*_GRADIENT, OUTPUT_DIVVEL, ...,
OUTPUT_IDORT(_NO_IV), OUTPUT_OPACITY: the module writes IdortTemperature/IdortOpacity/ERad/FRad/PRad),
TETRA_INDEX_IN_FACE, VORONOI_MESH_KEEP_DT_AND_DTC, NGB_TREE_DOUBLEPRECISION, BITS_PER_DIMENSION=42 +
OVERRIDE_PEANOGRID_WARNING, FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES, OUTPUT_COORDINATES_IN_DOUBLEPRECISION,
UPDATE_GRADIENTS_FOR_OUTPUT (not needed for this IC / these outputs), INPUT_IN_DOUBLEPRECISION (not used in this tree:
the IC precision is read from the file).
