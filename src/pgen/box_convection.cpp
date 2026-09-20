//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file box_convection.cpp
//! \brief A Cartesian plane-parallel LOCAL BOX of compressible convection, heated by a
//! fixed radiative flux through the bottom wall and cooled by a Newton-relaxation layer
//! at the top.  Written for the base of a red-giant convection zone, but nothing in it is
//! specific to that star: the state is fixed by (g0, rho_base, t_base) plus the EOS.
//!
//! GEOMETRY.  x1 is the vertical, gravity is the CONSTANT -g0 x1-hat, so the potential is
//! Phi(z) = g0 (z - x1min), increasing upward as the code's well-balanced machinery
//! expects.  x2 and x3 are periodic.  2-D (nx3 = 1) and 3-D both work.
//!
//! INITIAL STATE.  An exactly ISENTROPIC hydrostatic column, integrated with the RUN's
//! OWN EOS by marching in (p, T):
//!     dln p/dz = -rho g0/p,     dln T/dz = grad_ad(p,T) dln p/dz,
//!     rho = DensFromPT(p, T),   e = EintFromDensT(rho, T),
//! with grad_ad and DensFromPT taken from the tabulated EOS (pgen_eos_utils).  The march
//! runs on a fine grid covering the mesh plus its ghosts, in a ONE-THREAD DEVICE KERNEL:
//! a tabulated EOS lives in a DvceArray and cannot be evaluated from host code.  Because
//! the column is isentropic, `<hydro>/wb_option = isentropic` is the EXACT closure for
//! the dynamic well-balanced scheme, and the hydrostatic residual is at round-off.
//!
//! NOTE that a constant-g adiabat has a finite height: T falls linearly and reaches zero
//! at z = H_p(base)/grad_ad, i.e. 2.5 base scale heights for grad_ad = 0.4.  The domain
//! MUST end below that; the start-up report prints where it is.
//!
//! HEATING.  The bottom wall carries a fixed inward radiative flux, supplied by the
//! conduction module: set `<hydro>/rad_flux_inner` to the flux in erg/cm^2/s, with
//! ix1_bc = user (conduction.cpp adds it at the i = is face, and that path is
//! geometry-free -- the Cartesian RKUpdate divides by dx1).  The Rosseland opacity table
//! is read here and installed into the Conduction object, exactly as red_giant.cpp does.
//!
//! COOLING.  In the top `cool_depth` of the box the SPECIFIC internal energy is relaxed
//! toward the initial column's on a timescale `cool_tau`, with a smoothstep ramp:
//!     d(u_IEN)/dt = -ramp(z) rho (eps - eps_0(z))/cool_tau.
//! Relaxing the specific energy rather than the energy density keeps the sink from
//! forcing the density.  In steady state the layer removes exactly what the bottom
//! puts in.
//!
//! FLUX BOOST.  At a star's true flux the convective Mach number here is ~1e-3 and a
//! turnover costs 1e6-1e7 sound-limited steps.  Nothing in this file knows about that:
//! the boost is applied by the INPUT FILE, by raising `<hydro>/rad_flux_inner` and
//! lowering `<hydro>/rad_kappa_fac` by the same factor (so the radiative and convective
//! fluxes keep their ratio).  v ~ F^(1/3) and (nabla - nabla_ad) ~ F^(2/3) scale the
//! answer back.
//!
//! <problem> keys
//!   g0            constant gravity [cm/s^2] (positive; points to x1min)
//!   rho_base      density at x1min [g/cm^3]
//!   t_base        temperature at x1min [K]
//!   dgrad         SUPERADIABATIC EXCESS of the initial column: the march uses
//!                 dln T/dln p = grad_ad + dgrad instead of grad_ad.  0 gives the exact
//!                 adiabat, for which `wb_option = isentropic` is an exact closure and
//!                 the hydrostatic residual is at round-off.  A small positive value
//!                 starts the box CLOSE TO the convective steady state: building an
//!                 excess from an adiabat is limited by the imposed flux and takes a
//!                 Kelvin-Helmholtz time (hundreds of turnovers), whereas SHEDDING one is
//!                 done by the convection itself in a few turnovers.  The WB closure
//!                 stays second-order accurate at any dgrad; at 1e-2 the residual is
//!                 ~1e-4 of the gravity term.
//!   cool_depth    thickness of the cooling layer below x1max [cm]; < 0 -> 0.3 H_p(base)
//!   cool_tau      relaxation time of that layer [s]; < 0 -> 0.1 H_p(base)/v*
//!   bc_mode       the x1 walls.  0: the ghost carries the INITIAL column with the normal
//!                 velocity mirrored (what wb_column.cpp does -- correct only while the
//!                 state stays on that column; in a convecting box the ghost and the
//!                 evolved interior drift apart and the wall runs a steady WIND through
//!                 the box).  1: a plain reflecting mirror, which makes the wall-face
//!                 Riemann problem exactly symmetric and the mass flux exactly zero.
//!                 2 (DEFAULT): the mirror RESCALED by the initial column's own ratio
//!                 across the wall, rho_g = rho_m rho_col(z_g)/rho_col(z_m) and likewise
//!                 for e -- impermeable like 1 up to the stratification over one cell,
//!                 and hydrostatic like 0 at t = 0.  It LEAKS: the rescaled mirror is
//!                 not the hydrostatic continuation the WB background stencil would
//!                 build from the evolved state, so the wall-face Riemann problem keeps
//!                 a one-signed mass flux once the interior departs from the column.
//!                 3: the WB-CONSISTENT wall.  The ghost is the interior cell's own
//!                 (rho,e) walked across the wall with the SAME closure the well
//!                 balanced background uses (utils/wb_background.hpp's WBAdvance), so
//!                 the wall cell's background is an exact hydrostatic continuation of
//!                 what is actually there; on top of that the mass and energy the wall
//!                 face advected are cancelled cell by cell after every stage
//!                 (problem/wall_noflux, red_giant.cpp's treatment), which makes the box
//!                 EXACTLY closed whatever the interior does.
//!                 The walk is GUARDED: a wall cell the EOS table can only clamp (the
//!                 top of a radiation-dominated atmosphere) can send the closure's
//!                 exponential to zero density and, with eos_radiation, to an infinite
//!                 specific energy.  A ghost that is not finite, not positive, or whose
//!                 ratio to the mirror cell is more than problem/wall_walk_maxfac away
//!                 from what the initial column does over the same gap falls back to the
//!                 bc_mode-2 rescaled mirror, and both channels are then floored.
//!   bc_mode_top   the OUTER x1 (top) wall on its own.  -1 (DEFAULT) = whatever
//!                 bc_mode says, i.e. today's behaviour, bitwise.  4 = an OPEN
//!                 (OUTFLOW) top: gas may LEAVE through x1max, nothing comes in.
//!                 Every ghost layer is built from the LAST ACTIVE cell (ie), not
//!                 from the mirror cell, and is
//!                   v1_g = max(v1(ie), 0)          (zero gradient, no inflow)
//!                   v2_g = v2(ie), v3_g = v3(ie)   (zero gradient)
//!                   rho_g = rho(ie) exp(-g dz/(p/rho)(ie))   (ISOTHERMAL
//!                                                    hydrostatic continuation)
//!                   e_g  = eos.EnergyFromTemperature(rho_g, T(ie))
//!                 so the ghost column is in hydrostatic balance with the cell it
//!                 continues -- a plain zero-gradient copy of (rho,e) would put a
//!                 finite pressure gradient with no weight against the top face and
//!                 launch a spurious wind.  Continuing the TEMPERATURE rather than
//!                 the specific energy is what keeps the ghost on the EOS (with
//!                 eos_radiation, e at fixed T is NOT proportional to rho).  The
//!                 exponent is clamped to +-30 and the walk falls back to a plain
//!                 copy if the EOS hands back something that is not finite and
//!                 positive; both channels are floored afterwards.
//!                 wall_noflux is DISABLED ON THE TOP FACE in this mode (the point
//!                 of the mode is that the top face carries a flux); the bottom
//!                 wall, its imposed luminosity and its flux cancellation are
//!                 untouched, and so are the two-stream top (rt_top_vacuum) and the
//!                 vdamp sponge.  The box is then NOT closed: mass and energy leave
//!                 through x1max, which is the intended physics for a lid that is
//!                 pushed out by a near-Eddington flux.
//!   bc_mode_bot   the INNER x1 (bottom) wall on its own.  -1 (DEFAULT) = whatever
//!                 bc_mode says, i.e. today's behaviour, bitwise.  5 = an INFLOW
//!                 bottom, the companion of the open top: with bc_mode_top = 4 the
//!                 box SHEDS mass (the He column loses several per cent of its mass
//!                 per turnover once its vertical pulsation grows), and a closed
//!                 bottom cannot replenish it, so the run can never reach a
//!                 statistical steady state.  The ghost then carries
//!                   rho_g, T_g = the INITIAL (relaxed) COLUMN at the ghost's own
//!                                height -- the same values bc_mode 0 lays down and
//!                                bc_mode 3 falls back to -- with
//!                                e_g = eos.EnergyFromTemperature(rho_g, T_g)
//!                 and the wall itself stays IMPERMEABLE (the velocity is mirrored and
//!                 wall_noflux still cancels the wall-face mass and energy flux): the
//!                 inflow does NOT come through the wall face.  It is added instead as
//!                 an EXACT SOURCE in the cell against the wall -- rho_base v_in of
//!                 mass, rho_base v_in^2 of momentum and rho_base v_in h_base of energy
//!                 per unit area and time, h_base = (e+p)/rho + v_in^2/2 of the base
//!                 state -- so the box mass follows the controller exactly.  Handing
//!                 the bottom face to its own Riemann flux instead was TRIED and is a
//!                 DRAIN, not an inflow: an uncancelled wall face in this stratification
//!                 carries a one-signed acoustic mass flux that took the He column to
//!                 0.32 M0 in 50 s, ~300x the loss the open top makes.
//!                 v_in is set by a MASS-DEFICIT CONTROLLER,
//!                   v_in = min(bc_inflow_vmax,
//!                              max(0, (M0 - M(t))/(rho_base A tau_in))),
//!                 with M(t) the box mass (a global reduction over the active cells
//!                 of u0, done ONCE PER CYCLE), M0 the initial mass, A the horizontal
//!                 area, rho_base the initial column density at x1min and
//!                 tau_in = problem/bc_inflow_time.  The controller is proportional,
//!                 so a box losing mass at a steady rate Mdot settles at a deficit
//!                 Mdot*tau_in rather than at zero deficit; tau_in is the knob.
//!                 M0 is NOT carried in the restart file: it is recomputed at every
//!                 start as the exact sum over the active cells of the same initial
//!                 column interpolation the IC kernel uses, which is deterministic
//!                 (host-side, no MPI) and identical on a restart.
//!                 Everything else at the bottom is untouched: the imposed luminosity
//!                 still enters exactly as it did, through the two-stream's lower
//!                 boundary (rt_bottom_flux) or as <hydro>/rad_flux_inner in that same
//!                 face's energy channel.
//!   bc_inflow_ghost  false (DEFAULT) = the bottom ghost follows bc_mode (the WB
//!                 continuation of the evolved interior) and only the INJECTED mass
//!                 carries the base state; true = the ghost itself is the base state.
//!                 MEASURED: true is a real perturbation, because the wall-face
//!                 MOMENTUM flux is deliberately not cancelled and a frozen-column
//!                 ghost therefore pushes on the wall cell.  On the 1-D B-star column
//!                 (which loses no mass at all through the open top, so the controller
//!                 is idle) the base-state ghost drove E/E0 to 1.40 and Ftop/Fbot to
//!                 1.67 in 3.9e3 s, while the bc_mode ghost reproduced the closed-
//!                 bottom control to 6 digits; on the He column it moved Ftop/Fbot from
//!                 0.9995 to 1.021.
//!   bc_inflow_time  tau_in [code time] of that controller; <= 0 -> one turnover
//!                 H_p(base)/v*, which is what the startup line prints.
//!   bc_inflow_vmax  the cap on v_in [code velocity]; <= 0 -> 0.1 c_s(base).
//!   bc_inflow_print  print the controller state (cycle, t, v_in, M/M0) every N
//!                 cycles; default 100, 0 = off.
//!   wall_walk_maxfac  the slack in that test (default 100).
//!   wall_noflux   cancel the wall-face mass (and energy) flux after each stage.
//!                 Defaults to true under bc_mode 3 and false otherwise.  When a
//!                 diffusive flux (conduction, viscosity) has been added into the same
//!                 face -- at the bottom wall that term IS the imposed luminosity --
//!                 only the energy the cancelled mass carried is removed, as
//!                 red_giant.cpp does; with no diffusive flux the face energy flux is
//!                 cancelled exactly too.
//!   ic_profile    if set, a text file "z rho eint" (cgs, one node per line, increasing
//!                 z, '#' comments) REPLACES the isentropic march.  It must already cover
//!                 the ghosts -- analysis/mkprofile.py writes the horizontally and time
//!                 averaged profile of a finished case and pads both ends
//!                 hydrostatically.  This is how a finer case is started from a coarser
//!                 one's RELAXED stratification: an AMR-style upsampled restart is not
//!                 possible across different meshes, but the 1-D mean profile is exactly
//!                 the part that takes a thermal time to establish.  The flow itself
//!                 still has to grow from the seed, which takes a few turnovers, not a
//!                 Kelvin-Helmholtz time.
//!   vpert         seed amplitude: in units of the LOCAL sound speed for the default
//!                 velocity seed, or the RELATIVE internal-energy amplitude when
//!                 vpert_var = eint
//!   vpert_nk      NUMBER of random horizontal modes (not a wavelength range; default 16)
//!   vpert_seed    RNG seed for those modes (default 1234)
//!   vpert_kmin    smallest integer horizontal wavenumber drawn for k2 and k3, in units
//!   vpert_kmax    of 2 pi / L2 (L3); the draw is uniform on [kmin,kmax] (defaults 1, 4).
//!                 Raise kmin to keep the seed OFF the box-scale surface f-modes of a
//!                 wide box.
//!   vpert_zlo     height range (cm) the seed is applied in; the vertical envelope is
//!   vpert_zhi     sin(pi (z - zlo)/(zhi - zlo)) inside and zero outside.  Defaults are
//!                 the mesh x1min/x1max, i.e. the whole column.
//!   vpert_var     "v1" (default) seeds the vertical VELOCITY; "eint" instead multiplies
//!                 the internal energy (and so the pressure) by
//!                 1 + vpert*envelope*amp at FIXED density and zero velocity, i.e. an
//!                 entropy perturbation.  The linear f-mode seed (seed_fmode_amp) is
//!                 independent of this and always acts on the velocity.
//!   nfine         nodes in the column march (default 8192)
//!   opac_table    Rosseland table, "# nT nD lTmin dlT lDmin dlD" then nT*nD log10 kappa
//!   column_dump   if set, write the initial column to this file
//!   mu            ideal-gas branch only: mean molecular weight
//!   rt_two_stream  run the GREY TWO-STREAM of utils/two_stream_rt.hpp in its
//!                 PLANE-PARALLEL mode instead of relying on the cooling layer alone
//!                 (default false = today's behaviour).  The solver then does the
//!                 cooling: each cell exchanges with the radiation field through the
//!                 same Rosseland opacity the conduction operator uses, the top of the
//!                 domain radiates to space with nothing coming back in
//!                 (problem/rt_top_vacuum, default true), and the deep interior hands
//!                 over to radiative diffusion through the conduction module's tau blend
//!                 (<hydro>/rad_tau_lo, rad_tau_hi) exactly as red_giant.cpp does.  The
//!                 energy still enters at the bottom wall as <hydro>/rad_flux_inner.
//!                 Needs the whole vertical extent in ONE MeshBlock and nx1 + 2*nghost
//!                 <= 520 (the solver's compile-time column tiers).
//!   cool_layer    keep the Newton cooling layer.  Defaults to !rt_two_stream, i.e. on
//!                 by default and OFF as soon as the two-stream runs -- the two are two
//!                 models of the same loss and would double-count.
//!   user_srcs     must be true (gravity and the cooling layer live in the source term)
//!
//! THE EMERGENT-FLUX DIAGNOSTICS (both require problem/rt_two_stream = true).
//!   user_hist     enrol the five-column user history below.  The values are read out of
//!                 the two-stream's OWN face-flux array (two_stream_rt::rt_face_flux(),
//!                 the net longwave flux on the x1 faces), which lives for the run, so
//!                 the history samples the LAST RT call of the cycle -- no extra solve.
//!                 Columns, all box averages over the nx2*nx3 columns of the mesh:
//!                   Ftop   <F(ie+1)>, the net EMERGENT flux, erg/cm^2/s
//!                   Ftop2  <F(ie+1)^2>: the box rms about the mean is
//!                          sqrt(Ftop2 - Ftop^2).  The mean square, not the rms, is what
//!                          is written, because history columns are MPI_SUM-reduced
//!                          across ranks and a root is not a sum.
//!                   Fcut   <F(icut)>, the flux handed IN to the two-stream at the
//!                          deepest face it integrates -- the tau-blend handover.  The
//!                          sweep leaves every face below icut at zero, so the wall face
//!                          i = is carries nothing and is not what to read; the box's
//!                          true bottom flux is the constant <hydro>/rad_flux_inner,
//!                          which is what Ftop must equal in a steady state.
//!                   Ttop   <T(ie)> and
//!                   Ttop2  <T(ie)^2>, the top-cell temperature in K, taken from the
//!                          temperature ConsToPrim already solved (Hydro::wtemp).  The
//!                          tau = 2/3 surface temperature is NOT computed: it needs a
//!                          downward opacity integral per column, which the history
//!                          cadence cannot afford.  Use (Ftop/sigma_SB)^(1/4) instead.
//!                 The history file is written with <outputN>/data_format, which already
//!                 exists and defaults to "%12.5e".  THAT DEFAULT IS TOO COARSE FOR THIS:
//!                 100-s differences of a 4e13 flux fall below one print quantum.  Set
//!                 data_format = %24.16e in the hst output block.
//!   rt_budget_verbose  print the box ENERGY BUDGET term by term every N cycles, each
//!                 term box-integrated over the window and divided by F_bot A sum(dt).
//!                 Diagnostic only; see the block above BoxConvBoxInt for the terms.
//!   rt_surface_dt  cadence, in code time, of a per-column surface dump of F_top(x2,x3).
//!                 <= 0 (default) disables it.  Written from the source term right after
//!                 the RT call, once per cycle, appended to one file.
//!   rt_surface_file  that file (default "rt_surface.bin").  FORMAT: a stream of records,
//!                 each  [int64 ncol][float64 time][ncol x (float64 x2, x3, F_top)],
//!                 little-endian, no padding, no global header.  The columns are in no
//!                 particular order (they arrive rank by rank), which is why each row
//!                 carries its own (x2,x3) cell-centre coordinates.  Reader:
//!                   import numpy as np
//!                   recs = []
//!                   with open("rt_surface.bin","rb") as f:
//!                       while True:
//!                           h = f.read(8)
//!                           if len(h) < 8: break
//!                           n = np.frombuffer(h, "<i8")[0]
//!                           t = np.frombuffer(f.read(8), "<f8")[0]
//!                           a = np.frombuffer(f.read(24*n), "<f8").reshape(n, 3)
//!                           recs.append((t, a))          # a[:,0]=x2 a[:,1]=x3 a[:,2]=F
//!   rt_profile_dt  cadence, in code time, of a HORIZONTALLY AVERAGED x1 PROFILE dump.
//!                 <= 0 (default 0) disables it.  Written from the source term at the
//!                 TOP of the call, before any source touches the state, and only on
//!                 the first RK stage of a cycle (the same pm->time gate the surface
//!                 dump uses), so every record is the state at time = pm->time.  It
//!                 reads the PRIMITIVES w0 (and u0(IEN) for the internal energy, and
//!                 Hydro::wtemp for the temperature) as ConsToPrim last left them,
//!                 i.e. the start-of-cycle state.  Appended to one file.
//!   rt_profile_file  that file (default "rt_profile.bin").  FORMAT: a stream of
//!                 records, each
//!                   [float64 time][int32 nx1][int32 nvar]
//!                   [float64 x1v[nx1]][float64 data[nvar][nx1]]
//!                 little-endian, no padding, no global header.  nvar = 8, and each
//!                 row of data is the mean over the WHOLE horizontal plane (all x2,
//!                 all x3, all MeshBlocks, all ranks) at that x1 index, in code units:
//!                   0 <rho>          1 <v1>            2 <rho v1>
//!                   3 <v1^2>         4 <v2^2 + v3^2>   5 <T>  (x eos.temp_cgs for K)
//!                   6 <e>   internal energy DENSITY, u0(IEN) - rho v^2/2
//!                   7 <v1 (e + p)>   the enthalpy flux, p = w0(IPR)
//!                 Every MeshBlock spans the whole x1 extent (UserProblem fatals
//!                 otherwise), so the local i index IS the global one and the grid is
//!                 uniform in x1: x1v comes straight from the mesh extent.  Reader:
//!                 vis/python/read_rt_profile.py.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "globals.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "driver/driver.hpp"
#include "utils/wb_background.hpp"
#include "diffusion/conduction.hpp"
#include "units/units.hpp"
#include "utils/two_stream_rt.hpp"
#include "utils/rad_taper.hpp"
#include "pgen_eos_utils.hpp"
#include "pgen.hpp"

void BoxConvSrcs(Mesh *pm, Real bdt);
void BoxConvRebuildRadWeights(Mesh *pm, Real bdt);
void BoxConvTransverseApply(Mesh *pm, Real dt);
void BoxConvRTSplit(Mesh *pm, Real bdt);
void BoxConvRTBeforeFlux(Mesh *pm, Real bdt);
void BoxConvRTImEx(Mesh *pm, Driver *pdrive, const int estage);
void BoxConvBC(Mesh *pm);
void BoxConvFinal(ParameterInput *pin, Mesh *pm);
void BoxConvHistory(HistoryData *pdata, Mesh *pm);

namespace {
// the column, on a uniform fine grid covering the mesh plus its ghosts
DvceArray1D<Real> cd_, ce_, cp_, ct_;   // density, eint, pressure, temperature [K]
Real g0_ = 0.0, zlo_ = 0.0, dzf_ = 1.0, zmin_ = 0.0;
Real zcool_ = 0.0, zmax_ = 0.0, tcool_ = 1.0;
int nfine_ = 0, bc_mode_ = 2;
// problem/bc_mode_top (default -1 = the same wall as bc_mode).  The only extra
// value is 4: an OPEN (outflow) top.  See the header block.
int bc_mode_top_ = -1;
// problem/bc_mode_bot (default -1 = the same wall as bc_mode).  The only extra value
// is 5: an INFLOW bottom driven by the mass-deficit controller.  See the header block.
int bc_mode_bot_ = -1;
Real inflow_tau_ = 0.0;      // tau_in
Real inflow_vmax_ = 0.0;     // the cap on v_in
Real inflow_mass0_ = 0.0;    // M0, recomputed from the initial column at every start
Real inflow_area_ = 0.0;     // A, the mesh's horizontal area
Real inflow_rhob_ = 0.0;     // rho_base, the initial column density at x1min
Real inflow_vin_ = 0.0;      // the current v_in
Real inflow_mnow_ = 0.0;     // the current box mass
int inflow_cyc_ = -1;        // the cycle the controller was last updated on
int inflow_print_n_ = 100;   // print cadence in cycles; 0 = off
Real inflow_eb_ = 0.0, inflow_pb_ = 0.0;   // the base state's eint and pressure
bool inflow_ghost_ = false;  // bottom ghost = the base state (false: follow bc_mode)
bool etotgrav_ = false;
bool wall_noflux_ = false;   // cancel the wall-face flux after each stage (bc_mode 3)
Real wall_walk_maxfac_ = 100.0;   // how far the bc_mode-3 walk may depart from the column
// problem/vdamp_top_tau, problem/vdamp_top_time (default 0 = OFF, bitwise inert): a
// RAYLEIGH SPONGE on the VERTICAL velocity in the optically thin lid.  A closed box has
// a vertical acoustic fundamental that convection pumps; its amplitude grows like
// rho^-1/2, so it is the lid cells that evacuate to the density floor and collapse dt.
// The sponge damps m1 (only m1) implicitly per stage,
//     m1 -> m1/(1 + f bdt/vdamp_top_time),
// with the ramp f = 1 where the COLUMN optical depth tau <= vdamp_top_tau/3, f = 0 where
// tau >= vdamp_top_tau, a raised cosine in log tau between (the same shape the transverse
// conduction taper uses, and it reads the same per-plane tau, Conduction::rad_tauf).
// The kinetic energy the sponge removes is subtracted from the TOTAL energy: it is
// REMOVED from the box, not converted to heat.  This is a numerical sponge, not physics.
Real vdamp_tau_ = 0.0, vdamp_time_ = 20.0;
bool vdamp_printed_ = false;
// problem/vdamp_bot_cells, problem/vdamp_bot_time, problem/vdamp_bot_mean_only
// (default 0 = OFF, bitwise inert): the BOTTOM sponge, the companion of the top one.
// The closed bottom wall reflects the box's vertical acoustic fundamental (the He box:
// period ~150 s = twice the sound crossing from the wall to the acoustic cutoff in the
// lid) and the internal gravity waves convection launches downward into the stable
// layer; the star does neither, it lets both propagate away.  The sponge damps m1
// (only m1) in the lowest vdamp_bot_cells ACTIVE cells, i = is .. is+N-1, with the
// raised-cosine weight
//     f(i) = 0.5 (1 + cos(pi (i - is)/N)),
// i.e. f = 1 in the cell against the wall, falling monotonically to 0 in the first
// cell ABOVE the layer (the topmost sponge cell still carries a small weight, so the
// layer has no hard edge).  Per stage, with the exact exponential factor
//     g = 1 - exp(-f bdt/vdamp_bot_time),
//     m1 -> m1 - g m1                      (vdamp_bot_mean_only = false, DEFAULT), or
//     m1 -> m1 - g rho <v1>_h              (vdamp_bot_mean_only = true),
// where <v1>_h is the HORIZONTAL (x2,x3) mean of v1 at that x1 index, taken once per
// stage from the u0 the gravity/cooling kernel just left, with an MPI_Allreduce of the
// plane sums (every rank applies it).  The mean-only form exists for a box whose
// convection reaches the bottom: a plane-coherent v1 IS the radial mode, while a
// convective plume has zero plane mean, so damping the mean alone absorbs the mode and
// leaves the convection untouched.  In the He box the bottom layer is STABLY
// STRATIFIED (the convection zone is the top ~1 Mm), so the full-m1 form is the
// default: it absorbs the g-modes as well.  m2 and m3 are never touched.
// As for the top sponge, the kinetic energy removed is taken OUT of the total energy
// rather than turned into heat.  The sponge is applied BEFORE the bottom inflow source
// (bc_mode_bot = 5) in the same stage, so the momentum the inflow injects into the wall
// cell survives the stage intact.
int vdb_cells_ = 0;
Real vdb_time_ = 20.0;
bool vdb_mean_ = false;
DvceArray1D<Real> vdb_d_;        // (N): the plane sums of v1, then the plane means
HostArray1D<Real> vdb_h_;
bool diff_flux_ = false;     // a diffusive flux shares the wall face's energy channel
bool rt_on_ = false;      // problem/rt_two_stream
// problem/rt_strang (default false, bitwise off): take the grey two-stream OUT of the
// RK stage and apply it Strang-wise around the whole time integrator instead -- half
// the cycle dt before it, half after, through ProblemGenerator::user_split_func.  The
// in-stage call in BoxConvSrcs is then skipped.  It exists to test whether the residual
// dt-INDEPENDENT box-mode forcing that survives the exact column solve
// (rt_implicit_column = 3) is the operator SPLIT rather than the solve.
bool rt_strang_ = false;
// problem/rt_once_per_cycle (default false, bitwise off): like rt_strang, the two-stream
// leaves the RK stage -- but instead of two half-steps it is applied ONCE with the FULL
// cycle dt after the last stage ("after_timeintegrator", followed by ConToPrim).  It
// exists to test whether the SSP-RK stage AVERAGING of an exact per-stage relaxation
// (rt_implicit_column = 3 relaxes the column exactly inside every stage, and stage 2 of
// rk2 then blends that fully relaxed state back with the un-relaxed u^n) is what drives
// the residual dt-independent box mode.  Works for any rt_implicit_column mode.
bool rt_once_ = false;
// problem/rt_col3_once (default false, bitwise off): the SAME once-per-cycle split as
// rt_once_per_cycle, but only legal for the exact block-tridiagonal column solve
// (rt_implicit_column = 3).  In mode 3 the explicit sweep exists ONLY to build the
// column solve's coefficients -- the per-cell apply it would otherwise feed is replaced
// by the solve -- so taking the column out of the RK stage takes the sweep with it, and
// the whole radiation source is then applied ONCE with the FULL cycle dt after the last
// stage.  Halves the cost of the mode-3 source under rk2.
bool rt_col3_once_ = false;
// --- problem/rt_weights_per_stage: RE-CENTRE THE OPACITY/TAU CACHES ON THE STAGE
// STATE.  Conduction::BuildRadWeights (the column optical depth rad_tauf and the
// tau-blend weight rad_w) and Conduction::BuildAngularCoeffs (the transverse
// conductances cap_c2/cap_c3 and the frozen T*, alpha of the ADI operator) are both
// built inside Hydro::Fluxes, from w0 = the STAGE-START primitives.  The two-stream
// column and the radiative force, however, run on u0 AFTER the RK flux update and the
// gravity source.  The caches therefore lag the state they are used with by one
// sub-step, which is an O(dt) error in phase with rho'/T' -- the numerical-kappa
// mechanism the surface f-mode is being tested against.  With this on, the primitives
// are re-formed from the current u0 and both builders are re-run immediately before the
// two-stream call, so everything the column and the force consume is centred on the
// state they act on.  Default false and the code path is untouched when off.
bool rtwps_ = false;
// problem/rt_imex: the mode-3 column solve becomes the IMPLICIT STAGE OPERATOR of the
// ImEx-RK integrator (<time>/integrator = imex2 or imex2+) instead of a source applied
// inside the explicit RK stage.  The splitting error between the hydro update and the
// column -- the O(dt) anti-damping that grows the surface f-mode at gamma = 0.615 dt[s]
// per turnover -- is what this removes: ImEx is jointly second order in the explicit and
// the implicit operator, the first-order split is not.
//
// The RADIATION OPERATOR IS TAKEN AS A WHOLE: the energy exchange AND the radiative
// momentum force (problem/rt_rad_force) and its work are all part of S, so all four
// components (IM1,IM2,IM3,IEN) are stored and recombined.  The force is not stiff, but
// it is evaluated on the same fresh column solve as the energy source, so making it part
// of the same operator costs nothing, needs no hook inside two_stream_rt.hpp, and leaves
// no residual first-order split in the force either.
bool rt_imex_ = false;
// problem/rt_split_transverse (default TRUE whenever rt_strang / rt_once_per_cycle /
// rt_col3_once / rt_imex is on): move the IMPLICIT TRANSVERSE radiative operator
// (<hydro>/rad_implicit_ang, the horizontal ADI) OUT of the RK stage as well, and run it
// inside the same split step as the column, immediately after it and over the SAME bdt.
//
// WHY.  In the lid the column solve is a projection onto radiative equilibrium and the
// horizontal operator is a very fast diffusion (chi_rad ~ 5e17 cm^2/s); in the stage the
// two are applied back to back over the same beta_dt and very nearly cancel.  Taking
// only the column out leaves the horizontal operator acting, at its own stage weight, on
// a state the column has not relaxed -- which is the SAME inconsistency as running with
// the horizontal operator switched off (arm TR), and ARMS 4/5/7 show the same transonic
// lid in all three cases.  Set false to reproduce the ARMS 4/5 behaviour.
bool rt_split_tr_ = false;
// problem/rt_col3_sub (int, default 1 = today's behaviour, bitwise off): SUB-CYCLE the
// whole radiation column operator inside each RK stage.  The hydro step is untouched --
// the stage still advances over beta_dt -- but the radiation operator is applied N times
// with beta_dt/N each.
//
// WHY.  Every reordering of the step (rk2/rk3/imex2/imex2+/Strang, force and heating
// centering, weights, top BC: ARMS 1-9 in plaid/fmode_README.txt) leaves the surface
// f-mode growing at gamma = 0.615 x dt[s] per turnover.  The favoured reading is that the
// lid relaxes FULLY within one application of the operator, so the radiative work
// integral over the mode's cycle is evaluated with an effective thermal time equal to the
// step over which the operator is applied.  If that is right, gamma must follow that step
// and not the hydro dt: N = 2 halves it, N = 4 quarters it.  If gamma is unchanged (or
// doubles with the cost), the reading is refuted.
//
// WHAT IS IN A SUB-STEP: the whole operator, in the "column + transverse together"
// ordering rt_split_transverse established -- the mode-3 column solve (sweep + exact
// block-tridiagonal solve + the radiative force and its work, i.e. all of
// picket_fence_two_stream_RT) followed immediately by the horizontal ADI operator over
// the SAME beta_dt/N.  The in-stage transverse task is therefore made a no-op
// (Conduction::rad_tr_split_out), exactly as the split switches do.
//
// WHAT IS REFRESHED between sub-steps: BoxConvRebuildRadWeights -- one ConsToPrim over
// the full range, then BuildRadWeights (the tau/blend weights) and BuildAngularCoeffs
// (the frozen ADI face conductances, over the SUB-step beta_dt/N).  The column itself
// reads rho and e straight out of u0 (problem/rt_use_cons is required by mode 3), so it
// needs no refresh for its own state; what the ConsToPrim buys is a consistent wtemp --
// the general EOS's temperature GUESS -- and consistent weights.  The rebuild runs before
// EVERY sub-step, the first included, so that all N sub-steps see coefficients formed
// with the same beta_dt/N.
//
// GHOSTS: the column writes active cells only and columns are independent, so it needs no
// halo exchange between sub-steps.  The ADI operator recomputes T and 1/(rho cv) from the
// current u0 on the active cells and exchanges THOSE through its own one-variable halo
// (Conduction::RklConductionUpdate, pbval_tr), so each sub-step's transverse solve sees
// up-to-date transverse neighbours.  What stays as stale as it is today is the mesh u0
// ghost ring itself, which only the frozen face conductances on the outermost faces read.
int rt_col3_sub_ = 1;
// problem/rt_pair_sym (int, default 0 = today's ordering, bitwise off): make the COLUMN
// SOLVE and the HORIZONTAL (transverse ADI) OPERATOR a SYMMETRIC (Strang) PAIR inside
// each radiation call,
//     ADI(bdt/2)  ->  column(bdt) [+ the radiative force, as today]  ->  ADI(bdt/2)
// instead of the Lie ordering column(bdt) -> ADI(bdt) that rt_split_transverse
// established.  WHY: the arm-set-B work budget (problem/work_hist) found that the
// seeded linear f-mode's entire energy budget is the COLUMN HEATING (+1759 in units of
// 2 gamma E) against the HORIZONTAL ADI (-1798), each ~1800x the growth rate and
// cancelling to a net -38 -- i.e. the O(dt) growth IS the imbalance of that pair, which
// is exactly what a Lie split of two operators evaluated at different states produces.
// A symmetric pair makes that imbalance O(dt^2); the prediction is that Wtco+Wtad, and
// with it gamma, drop by a large factor at fixed dt.
//   value 1: the symmetric pair above.  It runs in the PLAIN in-stage path (the one
//     production uses), once per RK stage, and once per rt_col3_sub sub-step when that
//     is also on.  Cost: one extra ADI half-step plus one BoxConvRebuildRadWeights per
//     stage (the frozen ADI face conductances must be formed with the HALF step; the
//     SAME frozen coefficients are then used for both halves, which is what makes the
//     pair adjoint-symmetric).  It implies rt_split_transverse, i.e. it sets
//     Conduction::rad_tr_split_out so the separate in-stage transverse task is silenced
//     and this is the only place the operator runs.
//   value 2 (Picard iteration of the pair) is NOT implemented: it needs a saved copy of
//     the whole start state plus a way to feed the first pass's heating back into the
//     column solve as a lagged source, neither of which exists here (the mode-3 column
//     solve writes u0(IEN) in place from its own converged flux).  See the README.
// work_hist: tags 2/3/4 fire wherever the operators run, so the budget stays complete;
// Wtad is then the SUM of the two half-steps and Wtco+Wtad is still the pair imbalance.
int rt_pair_sym_ = 0;
// problem/rt_before_flux (bool, default false, bitwise off): REVERSE THE LIE ORDER of
// the two operators inside each RK stage.  Today the stage is
//     fluxes + RKUpdate (hydro advection)  ->  BoxConvSrcs (gravity, then the column +
//     radiative force, then the horizontal ADI)
// i.e. a first-order Lie split of the hydro advection against the stiff radiative
// projection.  The leading error of such a split is the COMMUTATOR of the two
// operators: it is proportional to the advecting (convective) velocity -- it vanishes in
// a box at rest -- and to dt, and it CHANGES SIGN when the two operators are exchanged.
// The convecting box's surface f-mode grows at gamma = 1.31 omega^2 dt while the same
// mode without convection is neutral at every dt, and every rearrangement INSIDE the
// radiation step has been null, which leaves the commutator as the candidate.
// With this switch the radiation block (column + force + horizontal ADI, exactly the
// sequence BoxConvSrcs runs today, with the same beta_dt) is applied at the HEAD of the
// stage, on the stage-start state, by Hydro::RTBeforeFlux -- which follows it with a
// full ghost/BC update and a ConToPrim, so the fluxes are built from the radiatively
// updated primitives -- and the radiation part of BoxConvSrcs is skipped for that stage.
// Gravity, cooling, sponge, wall cancellation and inflow stay exactly where they are.
// A sign flip of gamma identifies the commutator; an unchanged gamma refutes it.
// work_hist: the tag-2/3/4 points fire wherever the operators run, so the intervals
// still TILE the stage and the KE closure still holds exactly; only the ORDER changes
// (radiation first, then "Wflx" = the RK combination + flux divergence, then "Wgrv").
bool rt_before_flux_ = false;
// the per-stage implicit sources S^(l), (nimp_stages, nmb, 4, n3, n2, n1), 4 = the
// IM1/IM2/IM3/IEN components in that order.  Allocated on the first call (the Driver,
// which owns nimp_stages, is built after the problem generator).  NOT restarted: every
// slot is written before it is read inside the same cycle (see BoxConvRTImEx).
DvceArray6D<Real> rtimex_src_;
bool cool_on_ = true;     // the Newton cooling layer (off by default once RT is on)
Real rgas_ = 0.0;         // R/mu in code units; the ideal branch's T = p/(Rgas rho)
// --- THE LINEAR f-MODE TEST (problem/seed_fmode_amp, problem/fmode_hist) ------------
// A single horizontal Fourier mode (m,n) of the box, seeded in v1 with the surface-
// gravity-wave depth eigenfunction exp(k_h (z - z_top)), and a pair of history columns
// that project v1 back onto that mode.  Both are inert at their defaults: with
// seed_fmode_amp = 0 the initial state is bit-for-bit the unseeded one, and with
// fmode_hist false the history carries exactly the columns it carried before.  The
// point of the diagnostic is that the (m,n) amplitude of a LINEAR mode is a clean
// exponentially-damped/growing sinusoid, so ln|a| against t measures the scheme's
// numerical damping rate directly.
int fm_m_ = 1, fm_n_ = 1;        // the horizontal mode numbers of the projection
Real fm_kx_ = 0.0, fm_ky_ = 0.0; // 2 pi m / Lx, 2 pi n / Ly
Real fm_kh_ = 0.0;               // sqrt(kx^2 + ky^2)
Real fm_ztop_ = 0.0;             // x1max: the reference depth of the eigenfunction
Real fm_x2min_ = 0.0, fm_x3min_ = 0.0;
bool fm_hist_ = false;           // write the two projection columns
// --- the per-column emergent-flux surface dump (problem/rt_surface_dt) --------------
Real surf_dt_ = 0.0;             // <= 0 disables it
Real surf_next_ = -1.0;          // next dump time; armed at the first source call
char surf_file_[256] = "rt_surface.bin";
HostArray2D<Real> surf_h_;       // (ncol_local, 3): x2, x3, F_top
DvceArray2D<Real> surf_d_;
bool surf_alloc_ = false;

// --- the horizontally averaged x1 profile dump (problem/rt_profile_dt) --------------
constexpr int kNProf = 8;        // see the header block for what each slot is
Real prof_dt_ = 0.0;             // <= 0 disables it
Real prof_next_ = -1.0;          // next dump time; armed at the first source call
char prof_file_[256] = "rt_profile.bin";
HostArray2D<Real> prof_h_;       // (kNProf, nx1): the plane SUMS, then the means
DvceArray2D<Real> prof_d_;
bool prof_alloc_ = false;

// --- problem/rt_budget_verbose: THE BOX ENERGY BUDGET, TERM BY TERM -----------------
// Diagnostic only.  Every term is an ENERGY (erg), box-integrated and accumulated over
// the window between two prints, and is reported divided by W = F_bot * A * sum(bdt),
// so an exactly balanced box reads 1 in and 1 out.  The device slots are
//   0/1  the TOTAL x1-face energy flux at the i = is and i = ie+1 faces, i.e. what
//        Hydro::RKUpdate's flux divergence puts into / takes out of the box.  Both
//        carry the Riemann flux AND whatever AddIsotropicHeatFluxRadiative added to the
//        same channel, which at i = is is the imposed <hydro>/rad_flux_inner.  Every
//        interior face telescopes out of the box integral and x2/x3 are periodic, so
//        these two faces are the whole flux-divergence contribution.
//   2/3  the same faces' MASS flux (the bc_mode-3 leak, before its cancellation)
//   4/5  the bc_mode-3 wall ENERGY removal de, inner / outer wall
//   6/7  the same walls' MASS removal dm
//   10   the radiative force's WORK term v.f (filled by two_stream_rt.hpp)
// and the host-side accumulators carry the operator differences that are measured by
// re-integrating u0(IEN) across a call: 8 = the gravity/WB/cooling kernel, 9 = the
// two-stream call (deposition + v.f), 14 = the wall-correction kernel (the direct
// check on 4 + 5), 11/12 = the two-stream's own <Ftop>/<Fcut> area-integrated, 13 = the
// imposed bottom flux F_bot*A*bdt.
constexpr int kNBud = 16;
int  rtbud_n_ = 0;                 // print every N cycles; 0 = off
DvceArray1D<Real> rtbud_;          // device slots
Real rtbud_h_[kNBud];              // host accumulators
Real rtbud_e0_ = 0.0;              // box energy at the window start
Real rtbud_r0_ = 0.0;              // box LTE radiation energy w aT^4 at the window start
Real rtbud_area_ = 0.0;            // the mesh's horizontal area
Real rtbud_fin_ = 0.0;             // F_bot in code units
int  rtbud_cyc_ = -1;              // the cycle the current stage belongs to
bool rtbud_arm_ = false;

//----------------------------------------------------------------------------------------
//! \fn void BoxConvBoxInt
//! \brief box integrals over the ACTIVE cells: the conserved energy, and the LTE
//! radiation energy w(rho) a T^4 that the EOS taper keeps in it.  Summed over ranks.

void BoxConvBoxInt(Mesh *pm, Real &etot, Real &erad) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int ncell = pmbp->nmb_thispack*nx3*nx2*nx1;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  auto wt = pmbp->phydro->wtemp;
  auto eos = pmbp->phydro->peos->eos_data;
  const bool gen = eos.IsGeneral();
  const bool tap = gen && eos.tbl.rad_taper;
  const Real xlo = eos.tbl.rad_lrho_lo, xhi = eos.tbl.rad_lrho_hi;
  // the same w the EOS used, temperature gate included (rad_taper::WeightGated)
  const Real ylo = eos.tbl.rad_lt_lo, yhi = eos.tbl.rad_lt_hi;
  const bool tg = eos.tbl.rad_tgate;
  const Real arad = eos.tbl.arad;
  const Real tcgs = eos.temp_cgs;
  Real se = 0.0, sr = 0.0;
  Kokkos::parallel_reduce("boxconv_bint",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, ncell),
  KOKKOS_LAMBDA(const int idx, Real &le, Real &lr) {
    const int m = idx/(nx3*nx2*nx1);
    int r = idx - m*(nx3*nx2*nx1);
    const int k = ks + r/(nx2*nx1);
    r -= (r/(nx2*nx1))*(nx2*nx1);
    const int j = js + r/nx1;
    const int i = is + (r - (r/nx1)*nx1);
    const Real dv = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    le += u0(m,IEN,k,j,i)*dv;
    if (tap) {
      const Real tk = wt(m,k,j,i)*tcgs;
      Real wr, dwdx, dwdy;
      rad_taper::WeightGated(log10(w0(m,IDN,k,j,i)), xlo, xhi, log10(tk), ylo, yhi,
                             tg, wr, dwdx, dwdy);
      lr += wr*arad*tk*tk*tk*tk*dv;
    }
  }, se, sr);
  Kokkos::fence();
#if MPI_PARALLEL_ENABLED
  Real snd[2] = {se, sr}, rcv[2];
  MPI_Allreduce(snd, rcv, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  se = rcv[0]; sr = rcv[1];
#endif
  etot = se; erad = sr;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real BoxConvMassInt
//! \brief the box mass over the ACTIVE cells, summed over ranks.  The bottom-inflow
//! controller's M(t); called ONCE PER CYCLE, from the first boundary call of the cycle.

Real BoxConvMassInt(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int ncell = pmbp->nmb_thispack*nx3*nx2*nx1;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  Real sm = 0.0;
  Kokkos::parallel_reduce("boxconv_mint",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, ncell),
  KOKKOS_LAMBDA(const int idx, Real &lm) {
    const int m = idx/(nx3*nx2*nx1);
    int r = idx - m*(nx3*nx2*nx1);
    const int k = ks + r/(nx2*nx1);
    r -= (r/(nx2*nx1))*(nx2*nx1);
    const int j = js + r/nx1;
    const int i = is + (r - (r/nx1)*nx1);
    lm += u0(m,IDN,k,j,i)*size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
  }, sm);
  Kokkos::fence();
#if MPI_PARALLEL_ENABLED
  Real snd = sm, rcv;
  MPI_Allreduce(&snd, &rcv, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  sm = rcv;
#endif
  return sm;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvFtopInt
//! \brief the two-stream's own AREA-INTEGRATED net flux at the top face and at the cut.
//! Same arrays BoxConvHistory reads; valid only after a solver call.

void BoxConvFtopInt(Mesh *pm, Real &ftop, Real &fcut) {
  ftop = 0.0; fcut = 0.0;
  if (!two_stream_rt::rt_face_flux_ready()) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int ncol = pmbp->nmb_thispack*nx3*nx2;
  auto &size = pmbp->pmb->mb_size;
  auto fb = two_stream_rt::rt_face_flux();
  auto icut = two_stream_rt::rt_cut_index();
  const int nblk = two_stream_rt::rt_face_nblk();
  Real st = 0.0, sc = 0.0;
  Kokkos::parallel_reduce("boxconv_ftop",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, ncol),
  KOKKOS_LAMBDA(const int idx, Real &lt, Real &lc) {
    const int m = idx/(nx3*nx2);
    const int kj = idx - m*(nx3*nx2);
    const int k = ks + kj/nx2;
    const int j = js + (kj - (kj/nx2)*nx2);
    const Real da = size.d_view(m).dx2*size.d_view(m).dx3;
    const int ic = icut(m,k,j);
    for (int b=0; b<nblk; ++b) {
      lt += da*fb(m,b,ie+1,k,j);
      lc += da*fb(m,b,ic,k,j);
    }
  }, st, sc);
  Kokkos::fence();
#if MPI_PARALLEL_ENABLED
  Real snd[2] = {st, sc}, rcv[2];
  MPI_Allreduce(snd, rcv, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  st = rcv[0]; sc = rcv[1];
#endif
  ftop = st; fcut = sc;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvBudgetReport
//! \brief close the window: pull the device slots, print every term normalised by
//! W = F_bot A sum(bdt), and re-arm.

void BoxConvBudgetReport(Mesh *pm, const Real enow, const Real rnow) {
  auto hb = Kokkos::create_mirror_view(rtbud_);
  Kokkos::deep_copy(hb, rtbud_);
  Real dev[kNBud];
  for (int n=0; n<kNBud; ++n) dev[n] = hb(n);
#if MPI_PARALLEL_ENABLED
  Real rcv[kNBud];
  MPI_Allreduce(dev, rcv, kNBud, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int n=0; n<kNBud; ++n) dev[n] = rcv[n];
#endif
  const Real w = rtbud_h_[13];
  const Real iw = (w != 0.0) ? 1.0/w : 0.0;
  const Real dE = enow - rtbud_e0_;
  const Real dR = rnow - rtbud_r0_;
  // what the measured operators say the box should have gained
  const Real acc = dev[0] - dev[1] - dev[4] - dev[5] + rtbud_h_[8] + rtbud_h_[9];
  if (global_variable::my_rank == 0) {
    std::printf("### rt_budget ncycle=%d t=%.8e dtsum=%.6e W=%.8e (F_bot*A*dtsum)\n",
                pm->ncycle, pm->time, rtbud_h_[15], w);
    std::printf("###   1 in_bot_face   = %+.6f   (imposed F_bot = %+.6f)\n",
                dev[0]*iw, 1.0);
    std::printf("###   1b in_bot_mass  = %+.6e g\n", dev[2]);
    std::printf("###   2 out_top_face  = %+.6f   out_top_mass = %+.6e g\n",
                dev[1]*iw, dev[3]);
    std::printf("###   2b Ftop_2stream = %+.6f   Fcut_2stream = %+.6f\n",
                rtbud_h_[11]*iw, rtbud_h_[12]*iw);
    std::printf("###   3 wall_de_in    = %+.6f   wall_de_out  = %+.6f"
                "   (dE across the wall kernel = %+.6f)\n",
                -dev[4]*iw, -dev[5]*iw, rtbud_h_[14]*iw);
    std::printf("###   3b wall_dm_in   = %+.6e   wall_dm_out  = %+.6e g\n",
                dev[6], dev[7]);
    std::printf("###   4 rad_force_vf  = %+.6f\n", dev[10]*iw);
    std::printf("###   5 dE_two_stream = %+.6f   (deposition alone = %+.6f)\n",
                rtbud_h_[9]*iw, (rtbud_h_[9] - dev[10])*iw);
    std::printf("###   6 dE_grav_src   = %+.6f\n", rtbud_h_[8]*iw);
    std::printf("###   7 d(w aT^4)     = %+.6f   (box LTE radiation energy)\n", dR*iw);
    std::printf("###   9 dE_box        = %+.6f   E = %.16e\n", dE*iw, enow);
    std::printf("###   R residual      = %+.6f   (dE_box - sum of the terms above;"
                " the implicit x1 + RKL1 operators and the BCs are all that is left)\n",
                (dE - acc)*iw);
    std::fflush(stdout);
  }
  Kokkos::deep_copy(rtbud_, 0.0);
  for (int n=0; n<kNBud; ++n) rtbud_h_[n] = 0.0;
  rtbud_e0_ = enow;
  rtbud_r0_ = rnow;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn ReadOpacityTable
//! \brief the merged Rosseland table: comment lines, one of them
//! "# nT nD lTmin dlT lDmin dlD", then nT*nD values of log10 kappa_R with T slowest.
//! Same format red_giant.cpp reads; the box never goes near the table edges, so none of
//! that file's edge repair is carried over.

void ReadOpacityTable(const std::string &fname, DvceArray2D<Real> &tab,
                      DvceArray1D<Real> &lT, DvceArray1D<Real> &lD, int &nT, int &nD) {
  std::ifstream f(fname);
  if (!f.good()) {
    std::cout << "### FATAL ERROR in box_convection: cannot open problem/opac_table '"
              << fname << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string line;
  Real lt0 = 0.0, dlt = 0.0, ld0 = 0.0, dld = 0.0;
  bool have_grid = false;
  std::vector<Real> vals;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      if (!have_grid) {
        std::istringstream ss(line.substr(1));
        int a, b;
        Real c, d, e, g;
        if (ss >> a >> b >> c >> d >> e >> g) {
          nT = a; nD = b; lt0 = c; dlt = d; ld0 = e; dld = g;
          have_grid = true;
        }
      }
      continue;
    }
    vals.push_back(std::stod(line));
  }
  if (!have_grid || static_cast<int>(vals.size()) != nT*nD) {
    std::cout << "### FATAL ERROR in box_convection: opacity table '" << fname
              << "' has no grid line, or " << vals.size() << " values for "
              << nT << " x " << nD << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Kokkos::realloc(tab, nT, nD);
  Kokkos::realloc(lT, nT);
  Kokkos::realloc(lD, nD);
  auto htab = Kokkos::create_mirror_view(tab);
  auto hlT = Kokkos::create_mirror_view(lT);
  auto hlD = Kokkos::create_mirror_view(lD);
  for (int i=0; i<nT; ++i) {
    hlT(i) = lt0 + i*dlt;
    for (int j=0; j<nD; ++j) htab(i,j) = vals[i*nD + j];
  }
  for (int j=0; j<nD; ++j) hlD(j) = ld0 + j*dld;
  Kokkos::deep_copy(tab, htab);
  Kokkos::deep_copy(lT, hlT);
  Kokkos::deep_copy(lD, hlD);
  return;
}
//----------------------------------------------------------------------------------------
//! \fn par_file_tp
//! \brief p and T of a supplied (rho, e) column, for the start-up report and the dump.
//! A tabulated EOS lives in a DvceArray, so this has to happen on the device.

void par_file_tp(const EOS_Data &eos, DvceArray1D<Real> d, DvceArray1D<Real> e,
                 DvceArray1D<Real> p, DvceArray1D<Real> t, const int n, const Real gm1,
                 const Real rgas) {
  par_for("boxconv_filetp", DevExeSpace(), 0, n-1, KOKKOS_LAMBDA(const int i) {
    Real pp, tk;
    pgen_eos::PresTempFromEint(eos, gm1, rgas, d(i), e(i), -1.0, pp, tk);
    p(i) = pp;
    t(i) = tk;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn BoxConvSurfaceDump
//! \brief append one record of the per-column EMERGENT flux F_top(x2,x3).
//!
//! Called from the source term immediately after the two-stream, so what is written is
//! the flux the solver has just produced, on the state it was handed.  Every MeshBlock
//! holds the whole x1 extent -- UserProblem fatals otherwise, because the column sweep
//! needs it -- so the i = ie+1 face of every block IS the top of the box and every rank
//! owns whole columns.  Each row therefore carries its OWN (x2,x3) cell centre and the
//! reader never has to reconstruct the decomposition.  Format: see the file header.

void BoxConvSurfaceDump(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const int ncol = nmb*nx3*nx2;
  if (!surf_alloc_ || surf_d_.extent_int(0) != ncol) {
    Kokkos::realloc(surf_d_, ncol, 3);
    Kokkos::realloc(surf_h_, ncol, 3);
    surf_alloc_ = true;
  }
  auto fb = two_stream_rt::rt_face_flux();
  const int nblk = two_stream_rt::rt_face_nblk();
  auto &size = pmbp->pmb->mb_size;
  auto sd = surf_d_;
  par_for("boxconv_surf", DevExeSpace(), 0, ncol-1, KOKKOS_LAMBDA(const int idx) {
    const int m = idx/(nx3*nx2);
    const int kj = idx - m*(nx3*nx2);
    const int kk = kj/nx2;
    const int jj = kj - kk*nx2;
    Real ft = 0.0;
    for (int b=0; b<nblk; ++b) ft += fb(m,b,ie+1,ks+kk,js+jj);
    sd(idx,0) = CellCenterX(jj, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    sd(idx,1) = CellCenterX(kk, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    sd(idx,2) = ft;
  });
  Kokkos::deep_copy(surf_h_, surf_d_);

  std::vector<double> mine(3*ncol);
  for (int n=0; n<ncol; ++n) {
    for (int c=0; c<3; ++c) mine[3*n+c] = static_cast<double>(surf_h_(n,c));
  }
  std::vector<double> all;
  int ntot = ncol;
#if MPI_PARALLEL_ENABLED
  const int nr = global_variable::nranks;
  std::vector<int> cnt(nr, 0), disp(nr, 0);
  int mycnt = 3*ncol;
  MPI_Gather(&mycnt, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  int tot = 0;
  for (int r=0; r<nr; ++r) {
    disp[r] = tot;
    tot += cnt[r];
  }
  all.resize((global_variable::my_rank == 0) ? tot : 1);
  MPI_Gatherv(mine.data(), mycnt, MPI_DOUBLE, all.data(), cnt.data(), disp.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ntot = tot/3;
#else
  all.swap(mine);
#endif
  if (global_variable::my_rank != 0) return;
  FILE *pf = std::fopen(surf_file_, "ab");
  if (pf == nullptr) {
    std::cout << "### FATAL ERROR in box_convection: cannot append to "
              << "problem/rt_surface_file '" << surf_file_ << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int64_t n64 = static_cast<int64_t>(ntot);
  const double tnow = static_cast<double>(pm->time);
  std::fwrite(&n64, sizeof(int64_t), 1, pf);
  std::fwrite(&tnow, sizeof(double), 1, pf);
  std::fwrite(all.data(), sizeof(double), 3*ntot, pf);
  std::fclose(pf);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn BoxConvProfileDump
//! \brief append one record of the HORIZONTALLY AVERAGED x1 profile.
//!
//! Diagnostic only: it reads the primitives w0, the conserved energy u0(IEN) and the
//! temperature Hydro::wtemp that ConsToPrim left, and writes nothing back.  Called from
//! the top of BoxConvSrcs, once per cycle, so the record is the start-of-cycle state.
//!
//! Every MeshBlock spans the whole x1 extent (UserProblem fatals otherwise), so the
//! local index i - is IS the global x1 index and no LogicalLocation mapping is needed.
//! One team per x1 index reduces over that plane's (m,k,j); the plane SUMS are gathered
//! with MPI_Reduce and rank 0 divides by the global plane cell count.  Format: see the
//! file header block.

void BoxConvProfileDump(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const int nkj = nmb*nx3*nx2;
  if (!prof_alloc_ || prof_d_.extent_int(1) != nx1) {
    Kokkos::realloc(prof_d_, kNProf, nx1);
    Kokkos::realloc(prof_h_, kNProf, nx1);
    prof_alloc_ = true;
  }
  auto &w0 = pmbp->phydro->w0;
  auto &u0 = pmbp->phydro->u0;
  auto wt = pmbp->phydro->wtemp;
  auto pd = prof_d_;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nx1, Kokkos::AUTO);
  Kokkos::parallel_for("boxconv_prof", policy,
  KOKKOS_LAMBDA(Kokkos::TeamPolicy<>::member_type tmember) {
    const int i = is + tmember.league_rank();
    array_sum::GlobalSum sum;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkj),
    [&](const int idx, array_sum::GlobalSum &ls) {
      const int m = idx/(nx3*nx2);
      const int kj = idx - m*(nx3*nx2);
      const int k = ks + kj/nx2;
      const int j = js + (kj - (kj/nx2)*nx2);
      const Real d = w0(m,IDN,k,j,i);
      const Real v1 = w0(m,IVX,k,j,i);
      const Real v2 = w0(m,IVY,k,j,i);
      const Real v3 = w0(m,IVZ,k,j,i);
      const Real pg = w0(m,IPR,k,j,i);
      const Real ei = u0(m,IEN,k,j,i) - 0.5*d*(v1*v1 + v2*v2 + v3*v3);
      ls.the_array[0] += d;
      ls.the_array[1] += v1;
      ls.the_array[2] += d*v1;
      ls.the_array[3] += v1*v1;
      ls.the_array[4] += v2*v2 + v3*v3;
      ls.the_array[5] += wt(m,k,j,i);
      ls.the_array[6] += ei;
      ls.the_array[7] += v1*(ei + pg);
    }, Kokkos::Sum<array_sum::GlobalSum>(sum));
    Kokkos::single(Kokkos::PerTeam(tmember), [&]() {
      for (int n=0; n<kNProf; ++n) pd(n, i-is) = sum.the_array[n];
    });
  });
  Kokkos::fence();
  Kokkos::deep_copy(prof_h_, prof_d_);

  std::vector<double> buf(kNProf*nx1);
  for (int n=0; n<kNProf; ++n) {
    for (int i=0; i<nx1; ++i) buf[n*nx1+i] = static_cast<double>(prof_h_(n,i));
  }
#if MPI_PARALLEL_ENABLED
  std::vector<double> rbuf((global_variable::my_rank == 0) ? kNProf*nx1 : 1);
  MPI_Reduce(buf.data(), rbuf.data(), kNProf*nx1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  if (global_variable::my_rank == 0) buf.swap(rbuf);
#endif
  if (global_variable::my_rank != 0) return;
  const int nplane = pm->mesh_indcs.nx2*pm->mesh_indcs.nx3;
  const double fnorm = 1.0/static_cast<double>(nplane);
  for (int q=0; q<kNProf*nx1; ++q) buf[q] *= fnorm;
  // the grid is uniform in x1 and every block spans it, so the mesh extent gives x1v
  std::vector<double> x1v(nx1);
  const double x1lo = static_cast<double>(pm->mesh_size.x1min);
  const double x1hi = static_cast<double>(pm->mesh_size.x1max);
  for (int i=0; i<nx1; ++i) {
    x1v[i] = x1lo + (static_cast<double>(i) + 0.5)*(x1hi - x1lo)/nx1;
  }
  FILE *pfp = std::fopen(prof_file_, "ab");
  if (pfp == nullptr) {
    std::cout << "### FATAL ERROR in box_convection: cannot append to "
              << "problem/rt_profile_file '" << prof_file_ << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const double tnow = static_cast<double>(pm->time);
  const int32_t n1 = static_cast<int32_t>(nx1);
  const int32_t nv = static_cast<int32_t>(kNProf);
  std::fwrite(&tnow, sizeof(double), 1, pfp);
  std::fwrite(&n1, sizeof(int32_t), 1, pfp);
  std::fwrite(&nv, sizeof(int32_t), 1, pfp);
  std::fwrite(x1v.data(), sizeof(double), nx1, pfp);
  std::fwrite(buf.data(), sizeof(double), kNProf*nx1, pfp);
  std::fclose(pfp);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvWorkSnap / BoxConvWorkClose / the probes
//! \brief problem/work_hist: THE MODE-PROJECTED, PER-OPERATOR WORK INTEGRALS.
//!
//! Diagnostic only; nothing here writes back to the state, and with problem/work_hist
//! false (the default) not one kernel is launched and not one history column is added.
//!
//! WHAT IS MEASURED.  The seeded (m,n) mode is a LINEAR horizontal Fourier mode, so the
//! box integral of any work rate f.v against it vanishes at first order: the only work
//! that means anything is the work PROJECTED on the mode.  Define, for each x1 row and
//! with W(y,z) = cos(k_x (y-y0)) cos(k_y (z-z0)) the seeded horizontal pattern,
//!    s0 = sum_plane rho            s4 = sum_plane E W
//!    s1 = sum_plane (rho v1) W     s5 = sum_plane p W
//!    s2 = sum_plane (rho v2) W     s6 = sum_plane p
//!    s3 = sum_plane (rho v3) W
//! (plane sums over the WHOLE mesh: every MeshBlock spans x1, so a local i IS a global
//! row, and the plane sums are MPI-reduced onto rank 0, which owns the accumulators.)
//! The mode's momentum amplitude in the row is m_c = 4 s_c/N_h and its mean density
//! R = s0/N_h, with N_h the global plane cell count, so the MODE KINETIC ENERGY is
//!    E_mode = sum_rows (N_h/4) dV m_c^2/(2R)  =  sum_rows 2 dV (s1^2+s2^2+s3^2)/s0 ,
//! which is the "Emod" column.  Everything else is a CHANGE of E_mode across one
//! operator, or a mode-projected heating work.
//!
//! HOW THE TERMS ARE SEPARATED.  A snapshot is taken at five points of the cycle and
//! each interval is charged to the operator that lies inside it.  Between the end of one
//! interval and the start of the next NOTHING else writes u0 in the active cells (the
//! boundary exchange fills ghosts, ConToPrim writes w0), so the intervals TILE the cycle
//! and the KE columns close exactly:
//!    tag 0  entry of BoxConvSrcs        -> "Wflx"  the RK combination + flux divergence
//!                                          (the mode's own pressure and advection work)
//!    tag 1  after the gravity block     -> "Wgrv"  the well-balanced gravity source and
//!                                          the cooling/sponge/wall/inflow sources
//!    tag 2  after the mode-3 column     -> "Wtco"  the COLUMN HEATING's mode work
//!    tag 3  after picket_fence returns  -> "Wfrc"  the RADIATIVE FORCE's work on the
//!                                          mode, and "Wtfr" its v.f heating work
//!    tag 4  end of the transverse ADI   -> "Wtad"  the ADI OPERATOR's mode work
//! "Woth" collects the KE change of the two intervals that should not carry any (the
//! column and the ADI change energy, not momentum): it is the CLOSURE CONTROL, and
//!    Emod(t) - Emod(0)  =  Wflx + Wgrv + Wfrc + Woth
//! must hold to round-off.  What it does NOT contain is anything a floor or a C2P
//! correction puts back into u0, which is exactly what makes it a control.
//!
//! THE HEATING WORK.  An operator that only changes the energy does no work on the mode
//! directly; it drives (or damps) it through the pressure it leaves for the NEXT flux
//! step.  The standard linear measure of that is the pdV work of the heating in phase
//! with the compression, so for a heating that deposits dE per unit volume,
//!    W_th = sum_rows (N_h/4) dV Q (P/Pbar)  =  sum_rows 4 dV dE_proj s5/s6 ,
//! with Q = 4 dE_proj/N_h the mode amplitude of the deposited energy and P/Pbar = 4 s5/s6
//! the mode amplitude of the relative pressure perturbation.  The pressure weight is
//! taken from the state at the START of the interval (w0, i.e. stage-start primitives:
//! the pressure is not re-derived mid-stage).  Its phase error is dt/P ~ 1.6e-3 of a
//! period, far below the effect under test.  NOTE the normalisation: this is the work
//! integral up to the thermodynamic factor (Gamma_3 - 1), which is not applied -- a
//! POSITIVE Wtco/Wtad/Wtfr means heating in phase with compression, i.e. DRIVING, and
//! the terms are compared with each other and with 2 gamma E_mode, not used absolutely.
//!
//! ALL SEVEN "W" COLUMNS ARE CUMULATIVE (erg, summed over every stage since t = 0), so
//! the reader differentiates them and band-passes at the mode frequency.  "Sdsp" is the
//! instantaneous mode amplitude of v1 in the TOP active row (cm/s), 4 s1/s0 there, whose
//! time integral is the (m,n)-projected surface displacement.
//!
//! WHEN THE RADIATION OPERATORS MOVE (rt_strang, rt_imex, rt_col3_sub > 1) the tag-2/3/4
//! points still fire wherever the operator runs, but the intervals no longer tile the
//! stage in the order above; the KE closure still holds, the per-term split is still the
//! operator's own, and only the label "Wflx" becomes "everything since the last tag".
constexpr int kNWk = 7;
bool work_on_ = false;            // problem/work_hist
bool wk_alloc_ = false;
bool wk_armed_ = false;           // the first snapshot has been taken
int wk_nx1_ = 0;
Mesh *wk_pm_ = nullptr;           // for the two_stream_rt probe, which takes no Mesh
// heap-allocated and never freed, like the solver's own scratch: a file-scope Kokkos
// View would be destroyed AFTER Kokkos::finalize and abort the run at exit
DvceArray2D<Real> *wkd_ptr_ = nullptr;   // (kNWk, nx1) the plane sums
HostArray2D<Real> *wkh_ptr_ = nullptr;
std::vector<double> wk_now_, wk_old_;
double wk_dv_ = 0.0;              // the cell volume (uniform mesh)
double wk_acc_[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double wk_emod_ = 0.0, wk_sdsp_ = 0.0;

void BoxConvWorkSnap(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const int nkj = nmb*nx3*nx2;
  if (!wk_alloc_ || wkd_ptr_->extent_int(1) != nx1) {
    if (wkd_ptr_ == nullptr) {
      wkd_ptr_ = new DvceArray2D<Real>("boxconv_wkd", kNWk, nx1);
      wkh_ptr_ = new HostArray2D<Real>("boxconv_wkh", kNWk, nx1);
    } else {
      Kokkos::realloc(*wkd_ptr_, kNWk, nx1);
      Kokkos::realloc(*wkh_ptr_, kNWk, nx1);
    }
    wk_now_.assign(kNWk*nx1, 0.0);
    wk_old_.assign(kNWk*nx1, 0.0);
    wk_nx1_ = nx1;
    const double dx1 = static_cast<double>(pm->mesh_size.x1max - pm->mesh_size.x1min)
                     / static_cast<double>(pm->mesh_indcs.nx1);
    const double dx2 = static_cast<double>(pm->mesh_size.x2max - pm->mesh_size.x2min)
                     / static_cast<double>(pm->mesh_indcs.nx2);
    const double dx3 = static_cast<double>(pm->mesh_size.x3max - pm->mesh_size.x3min)
                     / static_cast<double>(pm->mesh_indcs.nx3);
    wk_dv_ = dx1*dx2*dx3;
    wk_alloc_ = true;
  }
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  auto &size = pmbp->pmb->mb_size;
  auto wd = *wkd_ptr_;
  const Real kx = fm_kx_, ky = fm_ky_;
  const Real x2m = fm_x2min_, x3m = fm_x3min_;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nx1, Kokkos::AUTO);
  Kokkos::parallel_for("boxconv_wksnap", policy,
  KOKKOS_LAMBDA(Kokkos::TeamPolicy<>::member_type tmember) {
    const int i = is + tmember.league_rank();
    array_sum::GlobalSum sum;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkj),
    [&](const int idx, array_sum::GlobalSum &ls) {
      const int m = idx/(nx3*nx2);
      const int kj = idx - m*(nx3*nx2);
      const int k = ks + kj/nx2;
      const int j = js + (kj - (kj/nx2)*nx2);
      const Real x2min = size.d_view(m).x2min, x2max = size.d_view(m).x2max;
      const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;
      const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      const Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      const Real cw = cos(kx*(x2v - x2m))*cos(ky*(x3v - x3m));
      const Real pg = w0(m,IPR,k,j,i);
      ls.the_array[0] += u0(m,IDN,k,j,i);
      ls.the_array[1] += u0(m,IM1,k,j,i)*cw;
      ls.the_array[2] += u0(m,IM2,k,j,i)*cw;
      ls.the_array[3] += u0(m,IM3,k,j,i)*cw;
      ls.the_array[4] += u0(m,IEN,k,j,i)*cw;
      ls.the_array[5] += pg*cw;
      ls.the_array[6] += pg;
    }, Kokkos::Sum<array_sum::GlobalSum>(sum));
    Kokkos::single(Kokkos::PerTeam(tmember), [&]() {
      for (int n=0; n<kNWk; ++n) wd(n, i-is) = sum.the_array[n];
    });
  });
  Kokkos::fence();
  auto wkh_ = *wkh_ptr_;
  Kokkos::deep_copy(wkh_, *wkd_ptr_);
  for (int n=0; n<kNWk; ++n) {
    for (int i=0; i<nx1; ++i) wk_now_[n*nx1+i] = static_cast<double>(wkh_(n,i));
  }
#if MPI_PARALLEL_ENABLED
  std::vector<double> rbuf((global_variable::my_rank == 0) ? kNWk*nx1 : 1);
  MPI_Reduce(wk_now_.data(), rbuf.data(), kNWk*nx1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  if (global_variable::my_rank == 0) wk_now_.swap(rbuf);
#endif
  return;
}

//! the mode kinetic energy of a snapshot, E_mode = sum_rows 2 dV (s1^2+s2^2+s3^2)/s0
double BoxConvWorkKE(const std::vector<double> &s) {
  const int n1 = wk_nx1_;
  double e = 0.0;
  for (int i=0; i<n1; ++i) {
    const double s0 = s[i];
    if (!(s0 > 0.0)) continue;
    const double s1 = s[n1+i], s2 = s[2*n1+i], s3 = s[3*n1+i];
    e += 2.0*wk_dv_*(s1*s1 + s2*s2 + s3*s3)/s0;
  }
  return e;
}

void BoxConvWorkClose(Mesh *pm, const int tag) {
  if (!work_on_) return;
  BoxConvWorkSnap(pm);
  if (global_variable::my_rank == 0) {
    if (wk_armed_) {
      const int n1 = wk_nx1_;
      // the KE bucket: 0 = flux, 1 = gravity, 2 = radiative force, 3 = the control
      const int kb = (tag == 0) ? 0 : ((tag == 1) ? 1 : ((tag == 3) ? 2 : 3));
      wk_acc_[kb] += BoxConvWorkKE(wk_now_) - BoxConvWorkKE(wk_old_);
      if (tag >= 2) {
        // the mode-projected heating work of this interval, with the pressure weight
        // taken at its START (see the header block)
        double wth = 0.0;
        for (int i=0; i<n1; ++i) {
          const double s6 = wk_old_[6*n1+i];
          if (!(s6 > 0.0)) continue;
          wth += 4.0*wk_dv_*(wk_now_[4*n1+i] - wk_old_[4*n1+i])*wk_old_[5*n1+i]/s6;
        }
        wk_acc_[(tag == 2) ? 4 : ((tag == 3) ? 5 : 6)] += wth;
      }
    }
    wk_armed_ = true;
    wk_emod_ = BoxConvWorkKE(wk_now_);
    const int n1 = wk_nx1_;
    const double s0t = wk_now_[n1-1];
    wk_sdsp_ = (s0t > 0.0) ? 4.0*wk_now_[n1 + n1-1]/s0t : 0.0;
  }
  wk_old_.swap(wk_now_);
  return;
}

//! the ProblemGenerator::user_probe_func hook (tag 4, the end of the transverse ADI)
void BoxConvWorkProbe(Mesh *pm, const int tag) {
  BoxConvWorkClose(pm, tag);
  return;
}

//! the two_stream_rt::rt_probe hook (tag 2, the end of the mode-3 column solve)
void BoxConvWorkProbeRT(const int tag) {
  if (wk_pm_ != nullptr) BoxConvWorkClose(wk_pm_, tag);
  return;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_srcs_func = BoxConvSrcs;
  user_bcs_func = BoxConvBC;
  pgen_final_func = BoxConvFinal;
  user_hist_func = BoxConvHistory;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in box_convection: <hydro> is required" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // --- problem/rt_budget_verbose: the term-by-term box energy budget (see the block
  // above BoxConvBoxInt).  Diagnostic only; nothing below changes a source term.
  rtbud_n_ = pin->GetOrAddInteger("problem", "rt_budget_verbose", 0);
  if (rtbud_n_ > 0) {
    rtbud_ = DvceArray1D<Real>("rtbud", kNBud);
    Kokkos::deep_copy(rtbud_, 0.0);
    for (int n=0; n<kNBud; ++n) rtbud_h_[n] = 0.0;
    two_stream_rt::rt_bud_ptr = &rtbud_;
    rtbud_area_ = (pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min)
                 *(pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min);
    rtbud_fin_ = (pmbp->phydro->pcond != nullptr)
               ? pmbp->phydro->pcond->rad_flux_inner : 0.0;
    rtbud_cyc_ = -1;
    rtbud_arm_ = false;
  }
  if (pin->GetOrAddBoolean("mesh", "use_cubed_sphere", false) ||
      pin->GetOrAddBoolean("mesh", "use_spherical_polar", false)) {
    std::cout << "### FATAL ERROR in box_convection: Cartesian meshes only" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;

  const Real g0 = pin->GetReal("problem", "g0");
  const Real rho_b = pin->GetReal("problem", "rho_base");
  const Real t_b = pin->GetReal("problem", "t_base");
  const Real vpert = pin->GetOrAddReal("problem", "vpert", 0.0);
  // --- the linear f-mode seed and its history projection (see the globals above)
  const Real fmamp = pin->GetOrAddReal("problem", "seed_fmode_amp", 0.0);
  fm_m_ = pin->GetOrAddInteger("problem", "seed_fmode_m", 1);
  fm_n_ = pin->GetOrAddInteger("problem", "seed_fmode_n", 1);
  fm_hist_ = pin->GetOrAddBoolean("problem", "fmode_hist", (fmamp > 0.0));
  // --- problem/work_hist: the mode-projected per-operator work integrals (see the
  // block above BoxConvWorkSnap).  Diagnostic only, and off by default.
  work_on_ = pin->GetOrAddBoolean("problem", "work_hist", false);
  if (work_on_) {
    if (!fm_hist_) {
      std::cout << "### FATAL ERROR in box_convection: problem/work_hist needs "
                << "problem/fmode_hist -- the work integrals are projections on the "
                << "SAME (m,n) mode the history columns track" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    wk_pm_ = pmy_mesh_;
    user_probe_func = BoxConvWorkProbe;      // tag 4: the transverse ADI operator
    two_stream_rt::rt_probe = BoxConvWorkProbeRT;   // tag 2: the mode-3 column solve
    if (global_variable::my_rank == 0) {
      std::cout << "### box_convection: problem/work_hist = true, the mode-projected "
                << "work integrals Wflx/Wgrv/Wfrc/Woth/Wtco/Wtfr/Wtad + Emod/Sdsp are "
                << "appended to the user history" << std::endl;
    }
  }
  // re-integrate the supplied ic_profile into EXACT hydrostatic balance under whatever
  // EOS is in force, keeping its T(z); needed when the EOS is changed under a profile
  // that was relaxed with a different one (the radiation-free linear f-mode test)
  const int hse_retune = pin->GetOrAddInteger("problem", "ic_hse_retune", 0);
  if (hse_retune < 0 || hse_retune > 2) {
    std::cout << "### FATAL ERROR in box_convection: problem/ic_hse_retune must be 0 "
              << "(off), 1 (keep T(z), solve p) or 2 (keep rho(z), solve p and T)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int nk = pin->GetOrAddInteger("problem", "vpert_nk", 16);
  const int kseed = pin->GetOrAddInteger("problem", "vpert_seed", 1234);
  // --- the horizontal wavenumber band and the vertical window of the random seed, plus
  // which variable it perturbs.  The defaults reproduce the historical seed exactly.
  const int kmin = pin->GetOrAddInteger("problem", "vpert_kmin", 1);
  const int kmax = pin->GetOrAddInteger("problem", "vpert_kmax", 4);
  if (kmin < 0 || kmax < kmin) {
    std::cout << "### FATAL ERROR in box_convection: need 0 <= problem/vpert_kmin <= "
              << "problem/vpert_kmax" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real pzlo = pin->GetOrAddReal("problem", "vpert_zlo",
                                      pmy_mesh_->mesh_size.x1min);
  const Real pzhi = pin->GetOrAddReal("problem", "vpert_zhi",
                                      pmy_mesh_->mesh_size.x1max);
  if (pzhi <= pzlo) {
    std::cout << "### FATAL ERROR in box_convection: need problem/vpert_zhi > "
              << "problem/vpert_zlo" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const std::string pvar_str = pin->GetOrAddString("problem", "vpert_var", "v1");
  int pvar = 0;
  if (pvar_str.compare("v1") == 0) {
    pvar = 0;
  } else if (pvar_str.compare("eint") == 0) {
    pvar = 1;
  } else {
    std::cout << "### FATAL ERROR in box_convection: problem/vpert_var must be "
              << "\"v1\" or \"eint\"" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int nfine = pin->GetOrAddInteger("problem", "nfine", 8192);
  const Real mu = pin->GetOrAddReal("problem", "mu", 1.3);
  const Real dgrad = pin->GetOrAddReal("problem", "dgrad", 0.0);
  bc_mode_ = pin->GetOrAddInteger("problem", "bc_mode", 2);
  if (bc_mode_ < 0 || bc_mode_ > 3) {
    std::cout << "### FATAL ERROR in box_convection: problem/bc_mode must be 0, 1, 2 or 3"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  bc_mode_top_ = pin->GetOrAddInteger("problem", "bc_mode_top", -1);
  if (bc_mode_top_ < -1 || bc_mode_top_ > 4) {
    std::cout << "### FATAL ERROR in box_convection: problem/bc_mode_top must be -1"
              << " (follow bc_mode) or 0-4" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  bc_mode_bot_ = pin->GetOrAddInteger("problem", "bc_mode_bot", -1);
  if (bc_mode_bot_ < -1 || bc_mode_bot_ > 5 || bc_mode_bot_ == 4) {
    std::cout << "### FATAL ERROR in box_convection: problem/bc_mode_bot must be -1"
              << " (follow bc_mode), 0-3, or 5 (inflow)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  wall_noflux_ = pin->GetOrAddBoolean("problem", "wall_noflux", (bc_mode_ == 3));
  wall_walk_maxfac_ = pin->GetOrAddReal("problem", "wall_walk_maxfac", 100.0);
  diff_flux_ = (pmbp->phydro->pcond != nullptr) || (pmbp->phydro->pvisc != nullptr);
  const std::string dump = pin->GetOrAddString("problem", "column_dump", "");
  const std::string icprof = pin->GetOrAddString("problem", "ic_profile", "");
  g0_ = g0;

  auto &eos = pmbp->phydro->peos->eos_data;
  auto &u0 = pmbp->phydro->u0;
  const Real gamma = eos.gamma;
  const Real igm1 = 1.0/(gamma - 1.0);
  const bool etotgrav = pmbp->phydro->use_etotgrav;
  const bool wbdyn = pmbp->phydro->use_wellbalance_dynamic;
  etotgrav_ = etotgrav;
  // ideal-gas branch only; the general branch carries composition in the table
  Real vunit = 1.0, lunit = 1.0, punit = 1.0;
  if (pmbp->punit != nullptr) {
    vunit = pmbp->punit->velocity_cgs();
    lunit = pmbp->punit->length_cgs();
    punit = pmbp->punit->pressure_cgs();
  }
  const Real rgas = 1.380649e-16/(mu*1.66053906660e-24)/(vunit*vunit);
  rgas_ = rgas;      // BoxConvHistory needs it for the ideal-gas temperature

  // --- the column, on a fine grid over the mesh's x1 extent plus its ghosts
  const Real zmin = pmy_mesh_->mesh_size.x1min;
  const Real zmax = pmy_mesh_->mesh_size.x1max;
  const Real dz = (zmax - zmin)/pmy_mesh_->mesh_indcs.nx1;
  const Real zlo = zmin - (ng + 1)*dz, zhi = zmax + (ng + 1)*dz;
  const Real dzf = (zhi - zlo)/(nfine - 1);
  zlo_ = zlo; dzf_ = dzf; zmin_ = zmin; zmax_ = zmax; nfine_ = nfine;
  DualArray1D<Real> cd("cd", nfine), ce("ce", nfine), cp("cp", nfine), ct("ct", nfine);
  auto cd_d = cd.d_view, ce_d = ce.d_view, cp_d = cp.d_view, ct_d = ct.d_view;
  const int i0 = static_cast<int>((zmin - zlo)/dzf + 0.5);
  // ONE THREAD: the tabulated EOS cannot be evaluated on the host, and the march is
  // sequential anyway.  RK2 (midpoint) in (ln p, T), isentropic by construction.
  par_for("boxconv_col", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(int) {
    // the base pressure from (rho_base, t_base), through the internal energy: the general
    // EOS evaluates e(rho,T) and p(rho,e,T) directly, so no root find is needed here
    const Real e_b = pgen_eos::EintFromDensT(eos, rgas, igm1, rho_b, t_b);
    Real p = pgen_eos::PresFromEint(eos, gamma - 1.0, rho_b, e_b);
    Real tk = t_b;
    cp_d(i0) = p;
    ct_d(i0) = tk;
    // upward
    Real pc = p, tc = tk;
    for (int i=i0+1; i<nfine; ++i) {
      Real d = pgen_eos::DensFromPT(eos, rgas, pc, tc);
      Real gad = pgen_eos::GradAd(eos, gamma, rgas, pc, tc) + dgrad;
      Real dlnp = -d*g0/pc;
      Real pm = pc*exp(0.5*dzf*dlnp);
      Real tm = tc + 0.5*dzf*tc*gad*dlnp;
      d = pgen_eos::DensFromPT(eos, rgas, pm, tm);
      gad = pgen_eos::GradAd(eos, gamma, rgas, pm, tm) + dgrad;
      dlnp = -d*g0/pm;
      pc = pc*exp(dzf*dlnp);
      tc = tc + dzf*tm*gad*dlnp;
      cp_d(i) = pc;
      ct_d(i) = tc;
    }
    // downward
    pc = p; tc = tk;
    for (int i=i0-1; i>=0; --i) {
      Real d = pgen_eos::DensFromPT(eos, rgas, pc, tc);
      Real gad = pgen_eos::GradAd(eos, gamma, rgas, pc, tc) + dgrad;
      Real dlnp = -d*g0/pc;
      Real pm = pc*exp(-0.5*dzf*dlnp);
      Real tm = tc - 0.5*dzf*tc*gad*dlnp;
      d = pgen_eos::DensFromPT(eos, rgas, pm, tm);
      gad = pgen_eos::GradAd(eos, gamma, rgas, pm, tm) + dgrad;
      dlnp = -d*g0/pm;
      pc = pc*exp(-dzf*dlnp);
      tc = tc - dzf*tm*gad*dlnp;
      cp_d(i) = pc;
      ct_d(i) = tc;
    }
  });
  par_for("boxconv_colde", DevExeSpace(), 0, nfine-1, KOKKOS_LAMBDA(const int i) {
    const Real d = pgen_eos::DensFromPT(eos, rgas, cp_d(i), ct_d(i));
    cd_d(i) = d;
    ce_d(i) = pgen_eos::EintFromDensT(eos, rgas, igm1, d, ct_d(i));
  });
  cd.modify_device();  cd.sync_host();
  ce.modify_device();  ce.sync_host();
  cp.modify_device();  cp.sync_host();
  ct.modify_device();  ct.sync_host();
  // --- an externally supplied stratification REPLACES the march
  if (!icprof.empty()) {
    std::ifstream pf(icprof);
    if (!pf.good()) {
      std::cout << "### FATAL ERROR in box_convection: cannot open problem/ic_profile '"
                << icprof << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::vector<Real> zf, df, ef;
    std::string line;
    while (std::getline(pf, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream is(line);
      Real a, b, c;
      if (!(is >> a >> b >> c)) continue;
      zf.push_back(a);
      df.push_back(b);
      ef.push_back(c);
    }
    if (zf.size() < 2 || zf.front() > zlo || zf.back() < zhi) {
      std::cout << "### FATAL ERROR in box_convection: problem/ic_profile has "
                << zf.size() << " nodes spanning ["
                << (zf.empty() ? 0.0 : zf.front()) << ", "
                << (zf.empty() ? 0.0 : zf.back())
                << "], which does not cover the mesh plus its ghosts ["
                << zlo << ", " << zhi << "].  Pad it with analysis/mkprofile.py."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::size_t kk = 0;
    for (int i=0; i<nfine; ++i) {
      const Real z = zlo + i*dzf;
      while (kk + 2 < zf.size() && zf[kk+1] < z) ++kk;
      const Real w = (z - zf[kk])/(zf[kk+1] - zf[kk]);
      // logarithmic in both: the profile spans two decades in density
      cd.h_view(i) = std::exp(std::log(df[kk])*(1.0 - w) + std::log(df[kk+1])*w);
      ce.h_view(i) = std::exp(std::log(ef[kk])*(1.0 - w) + std::log(ef[kk+1])*w);
    }
    cd.modify_host();  cd.sync_device();
    ce.modify_host();  ce.sync_device();
    // T and p of the supplied state, for the report and the dump only
    {
      auto cd_dv = cd.d_view, ce_dv = ce.d_view, cp_dv = cp.d_view, ct_dv = ct.d_view;
      par_file_tp(eos, cd_dv, ce_dv, cp_dv, ct_dv, nfine, gamma - 1.0, rgas);
    }
    cp.modify_device();  cp.sync_host();
    ct.modify_device();  ct.sync_host();
    if (global_variable::my_rank == 0) {
      std::cout << "box_convection: initial stratification READ FROM " << icprof
                << " (" << zf.size() << " nodes); the isentropic march is overridden"
                << std::endl;
    }
  }
  // --- problem/ic_hse_retune: RE-INTEGRATE the column into hydrostatic balance.
  // The march above is hydrostatic by construction, but a profile READ from a file is
  // only hydrostatic under the EOS it was relaxed with.  Change the EOS -- e.g. move the
  // LTE radiation taper so that aT^4 is carried at every density -- and the same (rho,e)
  // column is no longer in balance where the change bites.  This keeps the column's own
  // T(z) (the thermal structure is what the file is FOR) and re-solves dp/dz = -rho g
  // downward from the top node, with rho = rho(p,T) from the live EOS at every step:
  // RK2 midpoint in ln p on the same fine grid the march uses, so the result is
  // hydrostatic to the same order as the march and no better.
  if (hse_retune > 0) {
    auto cd_dv = cd.d_view, ce_dv = ce.d_view, cp_dv = cp.d_view, ct_dv = ct.d_view;
    const int nf = nfine;
    const int hmode = hse_retune;
    par_for("boxconv_hse", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(int) {
      Real pc = cp_dv(nf-1);
      if (hmode == 1) {
        // MODE 1 -- keep T(z), solve for p.  RK2 midpoint in ln p, rho = rho(p,T) from
        // the live EOS.  Fails by construction wherever the supplied T(z) makes the LTE
        // aT^4/3 alone exceed the hydrostatic p: there is then no rho >= 0 that fits.
        for (int i=nf-2; i>=0; --i) {
          const Real tm = 0.5*(ct_dv(i) + ct_dv(i+1));
          Real d = pgen_eos::DensFromPT(eos, rgas, pc, ct_dv(i+1));
          Real dlnp = -d*g0/pc;
          const Real pm = pc*exp(-0.5*dzf*dlnp);
          d = pgen_eos::DensFromPT(eos, rgas, pm, tm);
          dlnp = -d*g0/pm;
          pc = pc*exp(-dzf*dlnp);
          cp_dv(i) = pc;
        }
        for (int i=0; i<nf; ++i) {
          const Real d = pgen_eos::DensFromPT(eos, rgas, cp_dv(i), ct_dv(i));
          cd_dv(i) = d;
          ce_dv(i) = pgen_eos::EintFromDensT(eos, rgas, igm1, d, ct_dv(i));
        }
      } else {
        // MODE 2 -- keep rho(z), solve for p AND T.  This is the one that always has a
        // solution: at fixed rho the EOS pressure is strictly increasing in T from 0 to
        // unbounded, so a bisection in log T always brackets the hydrostatic p.  It also
        // leaves the DEEP column alone: dp/dz = -rho g is unchanged wherever rho is, so
        // p moves only by the constant the lid's missing radiative support contributes,
        // which is negligible against the deep pressure.  The trapezoidal integration is
        // exact for the piecewise-linear rho the fine grid carries.
        for (int i=nf-2; i>=0; --i) {
          pc += 0.5*dzf*g0*(cd_dv(i) + cd_dv(i+1));
          cp_dv(i) = pc;
        }
        for (int i=0; i<nf; ++i) {
          const Real d = cd_dv(i);
          Real tlo = 1.0e-3*ct_dv(i), thi = 1.0e3*ct_dv(i);
          for (int it=0; it<80; ++it) {
            const Real tmid = sqrt(tlo*thi);
            const Real em = pgen_eos::EintFromDensT(eos, rgas, igm1, d, tmid);
            const Real pm = pgen_eos::PresFromEint(eos, gamma - 1.0, d, em);
            if (pm < cp_dv(i)) { tlo = tmid; } else { thi = tmid; }
          }
          const Real tk = sqrt(tlo*thi);
          ct_dv(i) = tk;
          ce_dv(i) = pgen_eos::EintFromDensT(eos, rgas, igm1, d, tk);
          cp_dv(i) = pgen_eos::PresFromEint(eos, gamma - 1.0, d, ce_dv(i));
        }
      }
    });
    ct.modify_device();  ct.sync_host();
    cd.modify_device();  cd.sync_host();
    ce.modify_device();  ce.sync_host();
    cp.modify_device();  cp.sync_host();
    if (global_variable::my_rank == 0) {
      std::cout << "box_convection: problem/ic_hse_retune = true, the column was "
                << "RE-INTEGRATED into hydrostatic balance under the live EOS at fixed "
                << ((hse_retune == 1) ? "T(z)" : "rho(z)") << ", downward from p(top) = "
                << cp.h_view(nfine-1) << "; base (rho,T,p) = " << cd.h_view(i0) << " "
                << ct.h_view(i0) << " " << cp.h_view(i0) << std::endl;
    }
  }
  cd_ = cd.d_view; ce_ = ce.d_view; cp_ = cp.d_view; ct_ = ct.d_view;

  // --- the derived scales of the base state, and the cooling layer
  const Real p_b = cp.h_view(i0);
  const Real hp0 = p_b/(rho_b*g0);
  const Real g1 = pgen_eos::HostGamma1FromP(eos, rho_b, p_b);
  // grad_ad at the base, read back off the column: dlnT/dlnp between the two nodes
  const Real gad_b = std::log(ct.h_view(i0+1)/ct.h_view(i0-1))
                     /std::log(cp.h_view(i0+1)/cp.h_view(i0-1));
  const Real cs0 = std::sqrt(g1*p_b/rho_b);
  Conduction *pc = pmbp->phydro->pcond;
  Real fin = 0.0;
  if (pc != nullptr) fin = pc->rad_flux_inner*punit*vunit;   // back to erg/cm^2/s
  const Real vstar = (fin > 0.0) ? std::cbrt(fin/rho_b) : cs0;
  const Real tturn = hp0/vstar;
  Real cdep = pin->GetOrAddReal("problem", "cool_depth", -1.0);
  // --- THE BOTTOM-INFLOW CONTROLLER'S CONSTANTS (bc_mode_bot = 5) -------------------
  // M0 is recomputed here at EVERY start, restart included: it is the exact sum over
  // the active cells of the very column interpolation the IC kernel evaluates at each
  // cell centre (the x1 grid is uniform and every column is identical, so the sum is
  // one host-side loop over the mesh's nx1), which makes it deterministic and free of
  // any need to carry it in the restart file.
  inflow_tau_ = pin->GetOrAddReal("problem", "bc_inflow_time", -1.0);
  if (inflow_tau_ <= 0.0) inflow_tau_ = tturn;
  inflow_vmax_ = pin->GetOrAddReal("problem", "bc_inflow_vmax", -1.0);
  if (inflow_vmax_ <= 0.0) inflow_vmax_ = 0.1*cs0;
  inflow_print_n_ = pin->GetOrAddInteger("problem", "bc_inflow_print", 100);
  inflow_ghost_ = pin->GetOrAddBoolean("problem", "bc_inflow_ghost", false);
  inflow_area_ = (pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min)
                *(pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min);
  {
    const int nx1m = pmy_mesh_->mesh_indcs.nx1;
    const Real dz1 = (zmax - zmin)/static_cast<Real>(nx1m);
    Real msum = 0.0;
    for (int ic=0; ic<nx1m; ++ic) {
      const Real zc = zmin + (static_cast<Real>(ic) + 0.5)*dz1;
      Real sc = (zc - zlo)/dzf;
      int ii = static_cast<int>(sc);
      ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
      const Real fc = sc - ii;
      msum += cd.h_view(ii)*(1.0 - fc) + cd.h_view(ii+1)*fc;
    }
    inflow_mass0_ = msum*dz1*inflow_area_;
    Real sb = (zmin - zlo)/dzf;
    int ib = static_cast<int>(sb);
    ib = (ib < 0) ? 0 : ((ib > nfine-2) ? nfine-2 : ib);
    const Real fb = sb - ib;
    inflow_rhob_ = cd.h_view(ib)*(1.0 - fb) + cd.h_view(ib+1)*fb;
    inflow_eb_ = ce.h_view(ib)*(1.0 - fb) + ce.h_view(ib+1)*fb;
    inflow_pb_ = cp.h_view(ib)*(1.0 - fb) + cp.h_view(ib+1)*fb;
  }
  if (cdep < 0.0) cdep = 0.3*hp0;
  Real ctau = pin->GetOrAddReal("problem", "cool_tau", -1.0);
  if (ctau < 0.0) ctau = 0.1*tturn;
  zcool_ = zmax - cdep;
  tcool_ = ctau;
  if (zcool_ <= zmin) {
    std::cout << "### FATAL ERROR in box_convection: cooling layer fills the box"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // --- the top vertical-velocity sponge (see the note on vdamp_tau_ above)
  vdamp_tau_ = pin->GetOrAddReal("problem", "vdamp_top_tau", 0.0);
  vdamp_time_ = pin->GetOrAddReal("problem", "vdamp_top_time", 20.0);
  if (vdamp_tau_ > 0.0) {
    if (pc == nullptr || !pc->rad_tau_mode) {
      std::cout << "### FATAL ERROR in box_convection: vdamp_top_tau needs the per-plane "
                << "column optical depth, i.e. the tau blend (<hydro>/rad_tau_hi > 0)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(vdamp_time_ > 0.0)) {
      std::cout << "### FATAL ERROR in box_convection: vdamp_top_time must be > 0"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // --- the bottom vertical-velocity sponge (see the note on vdb_cells_ above)
  vdb_cells_ = pin->GetOrAddInteger("problem", "vdamp_bot_cells", 0);
  vdb_time_ = pin->GetOrAddReal("problem", "vdamp_bot_time", 20.0);
  vdb_mean_ = pin->GetOrAddBoolean("problem", "vdamp_bot_mean_only", false);
  if (vdb_cells_ > 0) {
    if (2*vdb_cells_ >= pmy_mesh_->mesh_indcs.nx1) {
      std::cout << "### FATAL ERROR in box_convection: problem/vdamp_bot_cells = "
                << vdb_cells_ << " must be < nx1/2 = "
                << pmy_mesh_->mesh_indcs.nx1/2 << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(vdb_time_ > 0.0)) {
      std::cout << "### FATAL ERROR in box_convection: vdamp_bot_time must be > 0"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // the sponge absorbs the standing mode of a CLOSED wall; on an open bottom it
    // would instead drive a flow through the boundary.  bc_mode_bot = 5 (inflow) is
    // allowed: that wall is impermeable too, the inflow is a source in the wall cell.
    if (pmy_mesh_->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::user) {
      std::cout << "### FATAL ERROR in box_convection: problem/vdamp_bot_cells needs a "
                << "CLOSED bottom wall, i.e. <mesh>/ix1_bc = user" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // --- the Rosseland table, for the conduction operator
  const std::string opac = pin->GetOrAddString("problem", "opac_table", "");
  if (pc != nullptr && pc->iso_cond_type.compare("radiative") == 0) {
    if (!opac.empty()) {
      if (!(pc->rad_kappa_tab && pc->rad_kappa_rho)) {
        std::cout << "### FATAL ERROR in box_convection: use <hydro>/rad_kappa_src = "
                  << "table_rho, the stellar table is on (T, rho)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      DvceArray2D<Real> ktab;
      DvceArray1D<Real> klT, klD;
      int knT = 0, knD = 0;
      ReadOpacityTable(opac, ktab, klT, klD, knT, knD);
      pc->rad_kr_tab = ktab;
      pc->rad_kr_lT = klT;
      pc->rad_kr_lP = klD;
      pc->rad_kr_nT = knT;
      pc->rad_kr_nP = knD;
    }
    if (pc->rad_kappa_rmax > 0.0) {
      std::cout << "### FATAL ERROR in box_convection: <hydro>/rad_kappa_rmax compares "
                << "against x1v, which a Cartesian mesh never allocates" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // --- the GREY TWO-STREAM (problem/rt_two_stream), utils/two_stream_rt.hpp
  // Everything the solver reads about the geometry and the star goes in here.  The
  // plane-parallel switch is what makes it legal on this mesh: see rt_plane_parallel.
  rt_on_ = pin->GetOrAddBoolean("problem", "rt_two_stream", false);
  // problem/rt_weights_per_stage (see rtwps_ above): rebuild the tau/blend weights and
  // the transverse conductances from the CURRENT state at every stage.  Read
  // unconditionally, so that the parameter always exists and can be set on the command
  // line.
  rtwps_ = pin->GetOrAddBoolean("problem", "rt_weights_per_stage", false);
  if (rtwps_ && rt_on_ && global_variable::my_rank == 0) {
    std::cout << "### box_convection: problem/rt_weights_per_stage = true -- the "
              << "column optical depth, the tau-blend weight and the transverse "
              << "conductances are rebuilt from the stage state" << std::endl;
  }
  cool_on_ = pin->GetOrAddBoolean("problem", "cool_layer", !rt_on_);
  if (rt_on_) {
    namespace ts = two_stream_rt;
    // the column sweep is vertical and private to a thread, so one MeshBlock must hold
    // the whole x1 extent; and the private intensity column is sized at compile time
    if (pmy_mesh_->mb_indcs.nx1 != pmy_mesh_->mesh_indcs.nx1) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_two_stream sweeps a "
                << "whole vertical column inside one MeshBlock, but mesh/nx1 = "
                << pmy_mesh_->mesh_indcs.nx1 << " and meshblock/nx1 = "
                << pmy_mesh_->mb_indcs.nx1 << ". Set meshblock/nx1 = mesh/nx1 and "
                << "decompose in x2/x3 only." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (indcs.nx1 + 2*ng > 520) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_two_stream needs "
                << "nx1 + 2*nghost <= 520 (the grey sweep dispatches on compile-time "
                << "column tiers, the largest of which is 520), but nx1 = "
                << indcs.nx1 << " and nghost = " << ng << " give "
                << (indcs.nx1 + 2*ng) << "." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // the solver is written in CGS throughout (sigma_SB, kappa in cm^2/g, p in bar), so
    // the code units have to BE cgs
    if (pmbp->punit == nullptr ||
        pmbp->punit->length_cgs() != 1.0 || pmbp->punit->density_cgs() != 1.0 ||
        pmbp->punit->velocity_cgs() != 1.0) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_two_stream needs cgs "
                << "code units (<units> with length_cgs = mass_cgs = time_cgs = 1); the "
                << "two-stream solver works in cgs internally." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pc == nullptr || pc->iso_cond_type.compare("radiative") != 0) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_two_stream takes its "
                << "opacity from the conduction module, so <hydro>/isotropic_conduction "
                << "= radiative with rad_kappa_src = table_rho is required." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ts::rt_grey = true;
    ts::rt_split = true;               // implied: the grey kernel lives on the split path
    ts::rt_plane_parallel = true;      // constant g, flat faces, geometry from RegionSize
    ts::rt_top_vacuum = pin->GetOrAddBoolean("problem", "rt_top_vacuum", true);
    ts::rt_top_re = pin->GetOrAddBoolean("problem", "rt_top_re", false);
    ts::rt_de_max = pin->GetOrAddReal("problem", "rt_de_max", 0.5);
    // Unlike red_giant.cpp, which defaults these to the pre-fix solver for backward
    // compatibility, a NEW problem generator defaults to the FIXED one -- the same six
    // opt-ins the production red-giant inputs state explicitly.
    ts::rt_semi_lin = pin->GetOrAddBoolean("problem", "rt_semi_lin", false);
    ts::rt_explicit = pin->GetOrAddBoolean("problem", "rt_explicit", false);
    ts::rt_newton = pin->GetOrAddBoolean("problem", "rt_newton", true);
    ts::rt_rescue_eq = pin->GetOrAddBoolean("problem", "rt_rescue_eq", true);
    ts::rt_ali_diag = pin->GetOrAddBoolean("problem", "rt_ali_diag", true);
    ts::rt_relax_sub = pin->GetOrAddInteger("problem", "rt_relax_sub", 1);
    ts::rt_relax_xcrit = pin->GetOrAddReal("problem", "rt_relax_xcrit", 1.0);
    ts::rt_relax_submax = pin->GetOrAddInteger("problem", "rt_relax_submax", 32);
    ts::rt_src_direct = pin->GetOrAddBoolean("problem", "rt_src_direct", true);
    ts::rt_src_dump = pin->GetOrAddInteger("problem", "rt_src_dump", 0);
    ts::rt_top_clamp = pin->GetOrAddBoolean("problem", "rt_top_clamp", true);
    ts::rt_use_cons = pin->GetOrAddBoolean("problem", "rt_use_cons", true);
    ts::rt_bface = pin->GetOrAddBoolean("problem", "rt_bface", true);
    // the deep-limit gradient in the grey sweep's upward intensity at the cut; the
    // fix is the default, the switch only buys back the old bit pattern
    ts::rt_cut_bc_legacy = pin->GetOrAddBoolean("problem", "rt_cut_bc_legacy", false);
    // the staggered layer source in the sweeps; see rt_layer_legacy.  The fix is
    // the default, the switch only buys back the old bit pattern
    ts::rt_layer_legacy = pin->GetOrAddBoolean("problem", "rt_layer_legacy", false);
    ts::rt_semi_implicit = pin->GetOrAddBoolean("problem", "rt_semi_implicit", true);
    ts::rt_outer_iter = pin->GetOrAddInteger("problem", "rt_outer_iter", 1);
    ts::rt_outer_verbose = pin->GetOrAddBoolean("problem", "rt_outer_verbose",
                                                false);
    // the merged implicit column solve: the two-stream source folded into the radial
    // implicit conduction tridiagonal (see two_stream_rt.hpp, rt_implicit_column)
    ts::rt_implicit_column = pin->GetOrAddInteger("problem", "rt_implicit_column", 0);
    ts::rt_impl_tol = pin->GetOrAddReal("problem", "rt_impl_tol",
                                        (ts::rt_implicit_column == 3) ? 1.0e-8 : 1.0e-6);
    ts::rt_col3_ex_iter = pin->GetOrAddBoolean("problem", "rt_col3_ex_iter",
                                                false);
    // mode 3 (the exact block-tridiagonal column solve) converges in 2-4 Newton steps
    // and is cheap per step, so it gets one more than the linearised modes by default
    ts::rt_impl_maxit = pin->GetOrAddInteger("problem", "rt_impl_maxit",
                                             (ts::rt_implicit_column == 3) ? 8 : 5);
    ts::rt_impl_exjac = pin->GetOrAddBoolean("problem", "rt_impl_exjac", true);
    ts::rt_impl_norm = pin->GetOrAddInteger("problem", "rt_impl_norm", 1);
    ts::rt_impl_norm_eps = pin->GetOrAddReal("problem", "rt_impl_norm_eps", 1.0e-3);
    ts::rt_impl_dstop = pin->GetOrAddBoolean("problem", "rt_impl_dstop", true);
    ts::rt_impl_rescheck = pin->GetOrAddBoolean("problem", "rt_impl_rescheck", true);
    ts::rt_impl_ablate = pin->GetOrAddInteger("problem", "rt_impl_ablate", 0);
    ts::rt_impl_fixit = pin->GetOrAddBoolean("problem", "rt_impl_fixit", false);
    if (ts::rt_impl_ablate != 0 || ts::rt_impl_fixit) {
      std::cout << "### WARNING in box_convection: problem/rt_impl_ablate or "
                << "rt_impl_fixit is set.  These are TIMING INSTRUMENTATION and the "
                << "mode-3 solve they produce is NOT a correct solve." << std::endl;
    }
    ts::rt_impl_cvfreeze = pin->GetOrAddInteger("problem", "rt_impl_cvfreeze", 0);
    // problem/rt_impl_reuse: reuse the block factorisation across Newton passes, with a
    // contraction check (1) or unconditionally (2).  0 (the default) is the old code.
    // problem/rt_impl_mixed: 0 = double (default), 1 = single-precision Newton
    // correction, 2 = and single-precision stored factors.  See two_stream_rt.hpp.
    ts::rt_impl_mixed = pin->GetOrAddInteger("problem", "rt_impl_mixed", 0);
    if (ts::rt_impl_mixed < 0 || ts::rt_impl_mixed > 2) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_impl_mixed must be "
                << "0, 1 or 2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ts::rt_impl_reuse = pin->GetOrAddInteger("problem", "rt_impl_reuse", 0);
    ts::rt_impl_reuse_rho = pin->GetOrAddReal("problem", "rt_impl_reuse_rho", 0.3);
    if (ts::rt_impl_reuse < 0 || ts::rt_impl_reuse > 2) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_impl_reuse must be "
                << "0, 1 or 2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // how often the mode-3 convergence line is printed under rt_outer_verbose
    ts::rt_report_every = pin->GetOrAddInteger("problem", "rt_report_every", 100);
    // how the mode-3 block system is solved: the serial block Thomas (one thread per
    // column) or the team-partitioned solve (see two_stream_column_partition.hpp)
    {
      std::string sv = pin->GetOrAddString("problem", "rt_impl_solver", "thomas");
      if (sv == "thomas") {
        ts::rt_impl_solver = 0;
      } else if (sv == "pcr") {
        ts::rt_impl_solver = 1;
      } else {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_impl_solver = "
                  << sv << " is not one of thomas, pcr" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    ts::rt_impl_redpar = pin->GetOrAddBoolean("problem", "rt_impl_redpar", false);
    ts::rt_col3_hybrid_tau = pin->GetOrAddReal("problem", "rt_col3_hybrid_tau", 0.0);
    ts::rt_col3_split_deep = pin->GetOrAddBoolean("problem", "rt_col3_split_deep",
                                                  false);
    ts::rt_col3_split_w = pin->GetOrAddInteger("problem", "rt_col3_split_w", 8);
    ts::rt_impl_nseg = pin->GetOrAddInteger("problem", "rt_impl_nseg", 64);
    if (ts::rt_impl_nseg < 1) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_impl_nseg must be >= 1"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // problem/rt_impl_warm: warm-start the mode-3 Newton from the previous call's
    // converged Planck function (1) or from a linear extrapolation of the last two (2).
    // 0 (the default) is the old code, bitwise.
    ts::rt_impl_warm = pin->GetOrAddInteger("problem", "rt_impl_warm", 0);
    if (ts::rt_impl_warm < 0 || ts::rt_impl_warm > 2) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_impl_warm must be "
                << "0, 1 or 2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ts::rt_impl_tau_min = pin->GetOrAddReal("problem", "rt_impl_tau_min", 1.0);
    ts::rt_impl_dtmax = pin->GetOrAddReal("problem", "rt_impl_dtmax", 0.25);
    ts::rt_impl_tau_blend = pin->GetOrAddReal("problem", "rt_impl_tau_blend", 1.0);
    rt_strang_ = pin->GetOrAddBoolean("problem", "rt_strang", false);
    rt_once_ = pin->GetOrAddBoolean("problem", "rt_once_per_cycle", false);
    rt_col3_once_ = pin->GetOrAddBoolean("problem", "rt_col3_once", false);
    if (rt_strang_ && rt_once_) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_strang and "
                << "problem/rt_once_per_cycle are mutually exclusive" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (rt_col3_once_ && (rt_strang_ || rt_once_)) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_col3_once is the "
                << "once-per-cycle split of the mode-3 column source and is mutually "
                << "exclusive with problem/rt_strang and problem/rt_once_per_cycle"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (rt_col3_once_ && ts::rt_implicit_column != 3) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_col3_once needs "
                << "problem/rt_implicit_column = 3 (only there does the sweep exist "
                << "solely to feed the column solve).  Got rt_implicit_column = "
                << ts::rt_implicit_column << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // ---- problem/rt_imex: the column solve as the ImEx implicit stage operator ----
    rt_imex_ = pin->GetOrAddBoolean("problem", "rt_imex", false);
    if (rt_imex_) {
      const std::string integ = pin->GetOrAddString("time", "integrator", "rk2");
      if (integ != "imex2" && integ != "imex2+") {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_imex needs "
                  << "<time>/integrator = imex2 or imex2+ (the a_twid/a_impl weights "
                  << "the implicit stages are combined with exist only there).  Got "
                  << "integrator = " << integ << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (ts::rt_implicit_column != 3) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_imex needs "
                  << "problem/rt_implicit_column = 3 -- the ImEx implicit stage must be "
                  << "an exact solve of the whole column, which is what mode 3 is.  Got "
                  << "rt_implicit_column = " << ts::rt_implicit_column << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rt_strang_ || rt_once_ || rt_col3_once_) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_imex is mutually "
                  << "exclusive with rt_strang, rt_once_per_cycle and rt_col3_once -- "
                  << "they are alternative ways of splitting the SAME source"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      user_imex_func = BoxConvRTImEx;
      if (global_variable::my_rank == 0) {
        std::cout << "### box_convection: problem/rt_imex = true, the mode-3 column "
                  << "solve is the IMPLICIT STAGE OPERATOR of " << integ
                  << " and is NOT applied as an in-stage source" << std::endl;
      }
    }
    // enrolled HERE, not next to user_srcs_func: the switch is read only now
    if (rt_strang_ || rt_once_ || rt_col3_once_) {
      user_split_func = BoxConvRTSplit;
      user_split_once = (rt_once_ || rt_col3_once_);
      if (global_variable::my_rank == 0) {
        if (rt_col3_once_) {
          std::cout << "### box_convection: problem/rt_col3_once = true, the mode-3 "
                    << "column source (sweep + exact column solve) is applied ONCE per "
                    << "cycle with the FULL dt after the last RK stage and is NOT "
                    << "applied inside the stages" << std::endl;
        } else if (rt_once_) {
          std::cout << "### box_convection: problem/rt_once_per_cycle = true, the grey "
                    << "two-stream is applied ONCE per cycle with the FULL dt after the "
                    << "last RK stage and is NOT applied inside the stages" << std::endl;
        } else {
          std::cout << "### box_convection: problem/rt_strang = true, the grey two-stream"
                    << " is STRANG-SPLIT around the time integrator (dt/2 before, dt/2 "
                    << "after) and is NOT applied inside the RK stages" << std::endl;
        }
      }
    }
    // ---- problem/rt_col3_sub: sub-cycle the whole column operator in the stage ----
    rt_col3_sub_ = pin->GetOrAddInteger("problem", "rt_col3_sub", 1);
    if (rt_col3_sub_ < 1) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_col3_sub must be >= 1 "
                << "(1 = the un-sub-cycled stage source).  Got " << rt_col3_sub_
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (rt_col3_sub_ > 1) {
      if (!rt_on_) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_col3_sub > 1 needs "
                  << "problem/rt_two_stream = true -- there is no column to sub-cycle"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (ts::rt_implicit_column != 3) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_col3_sub > 1 needs "
                  << "problem/rt_implicit_column = 3 -- a sub-step must be an EXACT "
                  << "solve of the whole column, which is what mode 3 is (the per-cell "
                  << "sub-cycler problem/rt_relax_sub is a different, inert, switch).  "
                  << "Got rt_implicit_column = " << ts::rt_implicit_column << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rt_strang_ || rt_once_ || rt_col3_once_ || rt_imex_) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_col3_sub > 1 "
                  << "sub-cycles the IN-STAGE column source and is mutually exclusive "
                  << "with rt_strang, rt_once_per_cycle, rt_col3_once and rt_imex, which "
                  << "all take that source out of the stage" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "### box_convection: problem/rt_col3_sub = " << rt_col3_sub_
                  << ", the WHOLE radiation operator (mode-3 column solve + radiative "
                  << "force + horizontal ADI) is applied " << rt_col3_sub_
                  << " times per RK stage with beta_dt/" << rt_col3_sub_
                  << " each, at the unchanged hydro dt" << std::endl;
      }
    }
    // ---- problem/rt_pair_sym: the column and the ADI as a symmetric pair -------
    rt_pair_sym_ = pin->GetOrAddInteger("problem", "rt_pair_sym", 0);
    if (rt_pair_sym_ < 0 || rt_pair_sym_ > 1) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_pair_sym must be 0 "
                << "(today's ordering) or 1 (ADI(dt/2) -> column(dt) -> ADI(dt/2)).  "
                << "Value 2 (Picard) is not implemented -- see the declaration.  Got "
                << rt_pair_sym_ << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (rt_pair_sym_ > 0) {
      if (rt_strang_ || rt_once_ || rt_col3_once_ || rt_imex_) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_pair_sym is the "
                  << "IN-STAGE pair ordering and is mutually exclusive with rt_strang, "
                  << "rt_once_per_cycle, rt_col3_once and rt_imex, which move the "
                  << "column out of the stage" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pmbp->phydro == nullptr || pmbp->phydro->pcond == nullptr ||
          !pmbp->phydro->pcond->rad_implicit_ang) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_pair_sym needs the "
                  << "implicit transverse operator (<hydro>/rad_implicit_ang or "
                  << "rad_sts_all): there is no horizontal operator to pair with"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "### box_convection: problem/rt_pair_sym = " << rt_pair_sym_
                  << ", the column solve and the horizontal ADI operator are a "
                  << "SYMMETRIC pair in every radiation call: ADI(bdt/2) -> column(bdt)"
                  << " + force -> ADI(bdt/2), with the face conductances formed at the "
                  << "HALF step" << std::endl;
      }
    }
    // ---- problem/rt_before_flux: the REVERSED Lie order within the stage ---------
    rt_before_flux_ = pin->GetOrAddBoolean("problem", "rt_before_flux", false);
    if (rt_before_flux_) {
      if (rt_strang_ || rt_once_ || rt_col3_once_ || rt_imex_ || rt_pair_sym_ > 0 ||
          rt_col3_sub_ > 1) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_before_flux is the "
                  << "IN-STAGE operator ORDER and is mutually exclusive with rt_strang, "
                  << "rt_once_per_cycle, rt_col3_once, rt_imex, rt_pair_sym and "
                  << "rt_col3_sub > 1, which move or repeat the radiation operator"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!rt_on_) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_before_flux needs "
                  << "the two-stream radiation (<problem>/rt_two_stream)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pmbp->phydro == nullptr || pmbp->phydro->pcond == nullptr ||
          !pmbp->phydro->pcond->rad_implicit_ang) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_before_flux needs "
                  << "the implicit transverse operator (<hydro>/rad_implicit_ang): the "
                  << "whole radiation operator, ADI included, has to move together"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      user_rt_before_flux = BoxConvRTBeforeFlux;
      if (global_variable::my_rank == 0) {
        std::cout << "### box_convection: problem/rt_before_flux = true, the RADIATION "
                  << "OPERATOR (column + force + horizontal ADI) runs at the HEAD of "
                  << "each RK stage, on the stage-start state, BEFORE the hydro flux "
                  << "update -- the REVERSED Lie order (commutator sign test)"
                  << std::endl;
      }
    }
    // problem/rt_split_transverse: with the column out of the stage, take the
    // horizontal ADI operator out with it (see the declaration above).  Default ON
    // whenever a split is active; a no-op otherwise.  rt_col3_sub > 1 keeps the column
    // INSIDE the stage but still has to move the transverse operator into the sub-cycle
    // loop, and uses the same Conduction flag to silence the in-stage task.
    const bool splitout = (rt_strang_ || rt_once_ || rt_col3_once_ || rt_imex_);
    rt_split_tr_ = (splitout || rt_col3_sub_ > 1 || rt_pair_sym_ > 0 ||
                    rt_before_flux_) &&
                   pin->GetOrAddBoolean("problem", "rt_split_transverse", true);
    if (rt_split_tr_ && pmbp->phydro != nullptr && pmbp->phydro->pcond != nullptr &&
        pmbp->phydro->pcond->rad_implicit_ang) {
      pmbp->phydro->pcond->rad_tr_split_out = true;
      if (global_variable::my_rank == 0) {
        std::cout << "### box_convection: problem/rt_split_transverse = true, the "
                  << "IMPLICIT TRANSVERSE radiative operator is run "
                  << ((rt_pair_sym_ > 0) ? "as the SYMMETRIC PAIR around the column"
                     : (rt_before_flux_ ? "with the column at the HEAD of the stage"
                     : ((rt_col3_sub_ > 1) ? "inside each column SUB-STEP"
                                           : "inside the split step with the column")))
                  << " and NOT as its own in-stage task" << std::endl;
      }
    }
    // ---- problem/rt_kappa_frozen: no opacity perturbation for the mode ----------
    // See the note in two_stream_rt.hpp.  The grey opacity kc_g is replaced by its
    // horizontal mean in each x1 row, which removes delta kappa from the column solve
    // AND from the radiative force (both read kc_g) while leaving the mean profile
    // free.  The transverse ADI operator builds its own conductances in Conduction and
    // is NOT covered by this switch.  Off = bitwise unchanged.
    ts::rt_kappa_frozen = pin->GetOrAddBoolean("problem", "rt_kappa_frozen", false);
    if (ts::rt_kappa_frozen) {
      if (!ts::rt_grey) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_kappa_frozen is a "
                  << "GREY-path switch (it freezes kc_g where the grey opacity kernel "
                  << "fills it)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "### box_convection: problem/rt_kappa_frozen = true, the grey "
                  << "opacity is the HORIZONTAL MEAN of each x1 row (delta kappa = 0 "
                  << "for the mode) in the column solve and the radiative force; the "
                  << "transverse ADI conductances are NOT frozen" << std::endl;
      }
    }
    ts::rt_apply_debug = pin->GetOrAddInteger("problem", "rt_apply_debug", 0);
    ts::rt_apply_debug_n = pin->GetOrAddInteger("problem", "rt_apply_debug_n", 8);
    ts::rt_nan_report = pin->GetOrAddBoolean("problem", "nan_report", false);
    ts::rt_dump_file = pin->GetOrAddString("problem", "rt_dump_file", "");
    ts::rt_dump_m = pin->GetOrAddInteger("problem", "rt_dump_m", 0);
    ts::rt_dump_j = pin->GetOrAddInteger("problem", "rt_dump_j", -1);
    ts::rt_dump_k = pin->GetOrAddInteger("problem", "rt_dump_k", -1);
    // ---- the per-column emergent-flux surface dump (see the header block) ----------
    surf_dt_ = pin->GetOrAddReal("problem", "rt_surface_dt", 0.0);
    {
      std::string sf = pin->GetOrAddString("problem", "rt_surface_file",
                                           "rt_surface.bin");
      if (sf.size() >= sizeof(surf_file_)) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_surface_file is "
                  << "longer than 255 characters" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::snprintf(surf_file_, sizeof(surf_file_), "%s", sf.c_str());
    }
    surf_next_ = -1.0;
    surf_alloc_ = false;
    // ---- the horizontally averaged x1 profile dump (see the header block) ---------
    prof_dt_ = pin->GetOrAddReal("problem", "rt_profile_dt", 0.0);
    {
      std::string pf = pin->GetOrAddString("problem", "rt_profile_file",
                                           "rt_profile.bin");
      if (pf.size() >= sizeof(prof_file_)) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_profile_file is "
                  << "longer than 255 characters" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::snprintf(prof_file_, sizeof(prof_file_), "%s", pf.c_str());
    }
    prof_next_ = -1.0;
    prof_alloc_ = false;
    // ---- the thin-region radiative force (see two_stream_rt.hpp, rt_rad_force) ------
    // It is the other half of the EOS's radiation taper: the taper removes (1-w) of the
    // LTE radiation pressure from the gas, and this puts the force that pressure was
    // carrying back as an explicit momentum source.  Neither half is meaningful alone,
    // so both are required together.
    ts::rt_rad_force = pin->GetOrAddBoolean("problem", "rt_rad_force", false);
    ts::rt_force_verbose = pin->GetOrAddInteger("problem", "rt_force_verbose", 0);
    ts::rt_force_grav = g0;
    if (ts::rt_rad_force) {
      if (!pmbp->phydro->peos->eos_data.tbl.rad_taper) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_rad_force is the "
                  << "momentum source that goes with the EOS radiation taper, and is "
                  << "only defined when the taper is on: set <hydro>/eos_rad_rho_hi and "
                  << "eos_rad_rho_lo (with eos_radiation = true), or drop rt_rad_force."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "box_convection: RT radiative momentum source ON -- "
                  << "(1-w) rho kappa F/c + Prad grad w, with w from the EOS taper"
                  << std::endl;
      }
    }
    // problem/rt_force_tau_gate: gate the force above by the optical depth to the top
    // (see two_stream_rt.hpp).  Default OFF and bitwise off.
    ts::rt_force_tau_gate = pin->GetOrAddBoolean("problem", "rt_force_tau_gate", false);
    ts::rt_force_tau_lo = pin->GetOrAddReal("problem", "rt_force_tau_lo", 0.3);
    ts::rt_force_tau_hi = pin->GetOrAddReal("problem", "rt_force_tau_hi", 3.0);
    if (ts::rt_force_tau_gate &&
        (!(ts::rt_force_tau_lo > 0.0) ||
         !(ts::rt_force_tau_hi > ts::rt_force_tau_lo))) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_force_tau_gate needs "
                << "0 < rt_force_tau_lo < rt_force_tau_hi (the ramp is a smoothstep "
                << "in log10 tau)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // ck_nquad: 1 = hemispheric mean (mu = 1/1.66), 2 = two-point Gauss-Legendre
    correlated_k::ck_nq = pin->GetOrAddInteger("problem", "rt_nquad", 2);
    // The internal flux.  It enters the box ONCE, through the bottom wall as
    // <hydro>/rad_flux_inner, so the solver must not inject it a second time at the cut;
    // with the tau blend on, two_stream_rt forces rt_int_at_cut false anyway.
    const Real teff_bot = (fin > 0.0) ? std::pow(fin/5.670374419e-5, 0.25) : 0.0;
    ts::rt_int_at_cut = pin->GetOrAddBoolean("problem", "rt_int_at_cut", false);
    // <problem>/rt_bottom_flux: hand the internal flux to the TWO-STREAM's own lower
    // boundary instead of to the conduction wall face.  Only meaningful when the tau
    // blend is off and the sweep therefore reaches the wall (see ts::rt_bot_flux): with
    // the blend on, the diffusion operator carries the deep flux and must keep the wall.
    const bool rt_botflux_ = pin->GetOrAddBoolean("problem", "rt_bottom_flux", false);
    if (rt_botflux_) {
      if (!pc->rad_tau_mode) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_bottom_flux puts the "
                  << "internal flux on the two-stream's lower boundary and takes it off "
                  << "the conduction wall face, which is only consistent when the sweep "
                  << "reaches that wall -- and the sweep's cut is the deepest face the "
                  << "TAU BLEND still leaves it, so <hydro>/rad_tau_hi must be set."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // ...and the blend weight has to be 0 on every face, i.e. rad_tau_lo deeper than
      // the whole box, or the diffusion operator still carries part of the column and
      // the wall it was fed through has just been taken away from it.  That cannot be
      // checked here (the weights are built per cycle), so it is stated, not enforced.
      if (global_variable::my_rank == 0) {
        std::printf("box_convection: rt_bottom_flux assumes <hydro>/rad_tau_lo = %.3e "
                    "is DEEPER than the bottom of the box, so the blend weight is 0 on "
                    "every face and the cut is the inner wall\n", pc->rad_tau_lo);
      }
      ts::rt_bot_flux = pc->rad_flux_inner;
      pc->rad_flux_inner = 0.0;
      if (global_variable::my_rank == 0) {
        std::printf("box_convection: F_bot = %.5e erg/cm^2/s is carried by the "
                    "TWO-STREAM's lower boundary (rt_bottom_flux); the conduction wall "
                    "face injects nothing\n", fin);
      }
    }
    // ---- problem/rt_col3_skip_sweep: drop the explicit entry sweep under mode 3 ----
    //
    // With the blend weight 0 on every face the column solve consumes NOTHING the sweep
    // makes except the frozen opacity and Planck function (kc, Bb), the cut index and
    // the stellar heating Qb -- all of which the pre-kernels build -- so the sweep is
    // pure cost.  What the sweep ALSO makes, and what the rest of the code reads, is the
    // face flux Fb (rad_f2s, rt_rad_force, the emergent-flux history and the surface
    // dump), the per-cell source Src and the emission rate Em.  Fb is taken over by the
    // column solve, which writes its OWN converged flux (see rt_col3_skip_sweep in
    // two_stream_rt.hpp); Src and Em are left at zero, so every reader of those two has
    // to be off.  Refuse anything that is not served.
    ts::rt_col3_skip_sweep = pin->GetOrAddBoolean("problem", "rt_col3_skip_sweep", false);
    if (ts::rt_col3_skip_sweep) {
      std::string bad;
      if (ts::rt_implicit_column != 3) {
        bad += "\n  problem/rt_implicit_column must be 3 (in every other mode the sweep "
               "IS the solver)";
      }
      if (!ts::rt_col3_ex_iter) {
        bad += "\n  problem/rt_col3_ex_iter must be true (frozen, the handover src_ex is "
               "read out of the sweep's Fb and Src)";
      }
      if (!ts::rt_src_direct) {
        bad += "\n  problem/rt_src_direct must be true (the non-direct handover reads "
               "the sweep's Fb)";
      }
      if (!rt_botflux_) {
        bad += "\n  problem/rt_bottom_flux must be true: it is this code's statement "
               "that <hydro>/rad_tau_lo is deeper than the box, i.e. that the blend "
               "weight is 0 on every face and the handover is identically zero (the "
               "weights themselves are checked once, on the first RT call)";
      }
      if (ts::rt_src_dump > 0) {
        bad += "\n  problem/rt_src_dump reads the sweep's Src";
      }
      if (ts::rt_apply_debug > 0) {
        bad += "\n  problem/rt_apply_debug reads the sweep's Em, Src and Qb";
      }
      if (ts::rt_diag) {
        bad += "\n  problem/rt_diag reads the sweep's Em";
      }
      if (ts::rt_cell_report) {
        bad += "\n  problem/rt_cell_report reads the sweep's per-face intensities";
      }
      if (!bad.empty()) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_col3_skip_sweep "
                  << "cannot be served here:" << bad << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "### box_convection: problem/rt_col3_skip_sweep = true -- the "
                  << "explicit two-stream entry sweep is NOT run; the column solve "
                  << "supplies the face flux Fb from its own converged intensities, and "
                  << "Src/Em stay zero" << std::endl;
        if (ts::rt_outer_verbose) {
          std::cout << "### WARNING in box_convection: with rt_col3_skip_sweep the "
                    << "rt_col3_fsum/rt_col3_flux dumps compare the column against an "
                    << "entry sweep that never ran: their Fb/sw_srcdx columns and the "
                    << "rt_col3 Fbtop slot are the PREVIOUS call's numbers, not this "
                    << "one's" << std::endl;
        }
      }
    }
    // ---- problem/rt_force_center, problem/rt_src_theta: the TIME CENTRING of the
    // two coupling terms the mode-3 column leaves first-order (see two_stream_rt.hpp).
    // Both are default-off and bitwise off; both need the exact column solve.
    ts::rt_force_center = pin->GetOrAddInteger("problem", "rt_force_center", 0);
    if (ts::rt_force_center < 0 || ts::rt_force_center > 2) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_force_center must be "
                << "0 (entry flux), 1 (converged flux) or 2 (the average)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (ts::rt_force_center > 0 && ts::rt_implicit_column != 3) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_force_center needs "
                << "problem/rt_implicit_column = 3 (only the exact column solve has a "
                << "converged flux to centre against)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ts::rt_src_theta = pin->GetOrAddReal("problem", "rt_src_theta", 1.0);
    if (ts::rt_src_theta < 0.0 || ts::rt_src_theta > 1.0) {
      std::cout << "### FATAL ERROR in box_convection: problem/rt_src_theta must be in "
                << "[0, 1]" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (ts::rt_src_theta != 1.0) {
      if (ts::rt_implicit_column != 3) {
        std::cout << "### FATAL ERROR in box_convection: problem/rt_src_theta != 1 needs "
                  << "problem/rt_implicit_column = 3" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (ts::rt_col3_skip_sweep) {
        // the explicit deposit IS the skipped sweep's source, so the sweep has to run:
        // turn the skip off rather than refuse the switch (it costs ~11 % per cycle)
        ts::rt_col3_skip_sweep = false;
        if (global_variable::my_rank == 0) {
          std::cout << "### WARNING in box_convection: problem/rt_src_theta != 1 needs "
                    << "the entry sweep's source, so problem/rt_col3_skip_sweep is "
                    << "turned OFF for this run" << std::endl;
        }
      }
    }
    if (global_variable::my_rank == 0 &&
        (ts::rt_force_center > 0 || ts::rt_src_theta != 1.0)) {
      std::cout << "### box_convection: RT coupling centring -- rt_force_center = "
                << ts::rt_force_center << " (0 entry flux, 1 converged flux, 2 average), "
                << "rt_src_theta = " << ts::rt_src_theta << std::endl;
    }
    ts::rt_tint_override = teff_bot;
    ts::rt_star_teff = 0.0;
    // The star-and-grid carrier the solver reads.  Teq = 0 switches the stellar beam off
    // and grav_point_mass = false makes EffGravAt return the box's constant g.
    hot_jupiter_param.Teq = 0.0;
    hot_jupiter_param.omega = 0.0;
    hot_jupiter_param.grav = g0;
    hot_jupiter_param.ap = 1.0;        // unused: no point mass, no tide, no beam
    hot_jupiter_param.Rgas = rgas;
    hot_jupiter_param.met = pin->GetOrAddReal("problem", "met", 0.0);
    hot_jupiter_param.grav_point_mass = false;
    hot_jupiter_param.stellar_tide = false;
    hot_jupiter_param.rot_potential = false;
    if (global_variable::my_rank == 0) {
      std::printf("box_convection: GREY two-stream ON (plane-parallel), %d-point "
                  "angular quadrature, top %s, T_int = %.5e K\n",
                  correlated_k::ck_nq, ts::rt_top_vacuum ? "VACUUM (no incoming "
                  "intensity)" : "unresolved hydrostatic column", teff_bot);
      std::printf("  deep handover: %s\n", pc->rad_tau_mode
                  ? "tau blend to radiative diffusion (<hydro>/rad_tau_lo, rad_tau_hi)"
                  : "NONE -- the sweep reaches the bottom wall; set rad_tau_lo/hi");
      std::printf("  Newton cooling layer: %s\n", cool_on_ ? "ON" : "off");
    }
  }

  // --- the start-up report: every number the design rests on
  if (global_variable::my_rank == 0) {
    std::printf("box_convection: base rho = %.5e g/cm^3, T = %.5e K, p = %.5e\n",
                rho_b, t_b, p_b);
    std::printf("  g0 = %.5e cm/s^2, H_p = %.5e cm, Gamma_1 = %.5f, c_s = %.5e cm/s\n",
                g0, hp0, g1, cs0);
    std::printf("  domain x1 = [%.5e, %.5e] = %.4f H_p; a constant-g adiabat reaches "
                "T = 0 at z - x1min ~ %.5e cm = %.3f H_p\n", zmin, zmax, (zmax-zmin)/hp0,
                hp0/gad_b, 1.0/gad_b);
    std::printf("  initial column: dln T/dln p = grad_ad + %.5e (dgrad; 0 = exact "
                "adiabat)\n", dgrad);
    std::printf("  column top:  T = %.5e K, rho = %.5e, p = %.5e, ln(p_b/p_top) = %.4f\n",
                ct.h_view(nfine-1-(ng+1)), cd.h_view(nfine-1-(ng+1)),
                cp.h_view(nfine-1-(ng+1)),
                std::log(p_b/cp.h_view(nfine-1-(ng+1))));
    std::printf("  F_bot = %.5e erg/cm^2/s -> v* = %.5e cm/s, Mach = %.4e, "
                "turnover H_p/v* = %.5e s\n", fin, vstar, vstar/cs0, tturn);
    std::printf("  cooling layer: z > %.5e (top %.4f H_p), tau = %.5e s = %.4f"
                " turnover\n", zcool_, cdep/hp0, ctau, ctau/tturn);
    if (vdamp_tau_ > 0.0) {
      std::printf("  top v1 sponge: f = 1 for tau <= %.4e, 0 for tau >= %.4e,"
                  " timescale %.4e s (the cell range is printed at the first stage)\n",
                  vdamp_tau_/3.0, vdamp_tau_, vdamp_time_);
    }
    if (vdb_cells_ > 0) {
      const Real zbot = zmin + vdb_cells_*(zmax - zmin)/pmy_mesh_->mesh_indcs.nx1;
      std::printf("  bottom v1 sponge: i = is .. is+%d, z <= %.5e (%.4f H_p),"
                  " f = 0.5(1+cos(pi (i-is)/N)), timescale %.4e s, %s\n",
                  vdb_cells_-1, zbot, (zbot - zmin)/hp0, vdb_time_,
                  vdb_mean_ ? "the PLANE MEAN of v1 only" : "the FULL v1");
    }
    std::printf("  x1 walls: bc_mode = %d (0 column ghost, 1 mirror, 2 mirror x the"
                " column ratio, 3 WB continuation), wall_noflux = %d%s\n",
                bc_mode_, static_cast<int>(wall_noflux_),
                (wall_noflux_ && diff_flux_) ? " (mass exact, energy via enthalpy)"
                                             : (wall_noflux_ ? " (mass and energy exact)"
                                                             : ""));
    if (bc_mode_top_ >= 0) {
      std::printf("  x1 TOP wall overridden: bc_mode_top = %d%s\n", bc_mode_top_,
                  (bc_mode_top_ == 4)
                  ? " = OPEN (outflow: hydrostatic isothermal ghost, v1 >= 0,"
                    " no top-face flux cancellation)"
                  : "");
    }
    if (bc_mode_bot_ >= 0) {
      std::printf("  x1 BOTTOM wall overridden: bc_mode_bot = %d%s\n", bc_mode_bot_,
                  (bc_mode_bot_ == 5)
                  ? " = INFLOW (base-state ghost under bc_inflow_ghost, wall still"
                    " impermeable, v_in injected as a source in the wall cell)"
                  : "");
    }
    if (bc_mode_bot_ == 5) {
      std::printf("  bottom inflow controller: M0 = %.8e (A = %.5e, rho_base = %.5e),"
                  " tau_in = %.5e s = %.4f turnover, v_in <= %.5e cm/s = %.4f c_s,"
                  " e_base = %.5e, p_base = %.5e, ghost = %s, print every %d cycles\n",
                  inflow_mass0_, inflow_area_, inflow_rhob_,
                  inflow_tau_, inflow_tau_/tturn, inflow_vmax_, inflow_vmax_/cs0,
                  inflow_eb_, inflow_pb_,
                  inflow_ghost_ ? "base state" : "bc_mode", inflow_print_n_);
    }
    if (pc != nullptr) {
      std::printf("  rad_kappa_fac = %.5e (conductivity is 1/rad_kappa_fac x physical)\n",
                  pc->rad_kappa_fac);
    }
    if (ct.h_view(nfine-1) <= 0.0) {
      std::printf("### FATAL: the adiabat runs out of temperature inside the domain;"
                  " shorten x1max\n");
      std::exit(EXIT_FAILURE);
    }
  }
  if (!dump.empty() && global_variable::my_rank == 0) {
    std::ofstream fo(dump);
    fo << "# z rho T p eint\n";
    for (int i=0; i<nfine; ++i) {
      fo << (zlo + i*dzf) << " " << cd.h_view(i) << " " << ct.h_view(i) << " "
         << cp.h_view(i) << " " << ce.h_view(i) << "\n";
    }
  }
  // --- THE GRAVITATIONAL POTENTIAL, WHICH A RESTART DOES NOT CARRY.  restart.cpp
  // writes and reads u0 (and the face fields) only: phicc0 and phi0 are NOT restart
  // state, they are pgen state, and they must be rebuilt on EVERY start.  They used to
  // be filled inside the initial-condition kernel below, i.e. after the `if (restart)`
  // return, so a restarted run ran with a potential of exactly zero -- and this problem
  // needs it twice over.  With <hydro>/etotgrav the conserved energy CARRIES rho*phi, so
  // ConToPrim recovers the internal energy by subtracting d*phicc: with phicc = 0 every
  // cell's e_int came back too large by rho*g0*(z - zmin), which at the top of the box is
  // orders of magnitude above e_int itself.  With <hydro>/wellbalance_dynamic the source
  // term's stencil reads phicc and phi0.x1f directly, so the hydrostatic balance the
  // scheme is built on was being integrated against a flat potential.  Both show up in
  // cycle 0 of the restarted run: the timestep collapses and the two-stream clips on
  // rt_de_max.  red_giant.cpp fills its potential before its own restart return for
  // exactly this reason.
  {
    DvceArray4D<Real> phicc = pmbp->phydro->phicc0;
    DvceArray4D<Real> ph1 = pmbp->phydro->phi0.x1f;
    DvceArray4D<Real> ph2 = pmbp->phydro->phi0.x2f;
    DvceArray4D<Real> ph3 = pmbp->phydro->phi0.x3f;
    const bool have_phi = (etotgrav || wbdyn);
    if (have_phi) {
      par_for("boxconv_phi", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
        const Real z = CellCenterX(i-is, indcs.nx1, x1min, x1max);
        const Real x1l = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
        const Real x1r = LeftEdgeX(i+1-is, indcs.nx1, x1min, x1max);
        const Real phi_c = g0*(z - zmin);
        phicc(m,k,j,i) = phi_c;
        ph1(m,k,j,i) = g0*(x1l - zmin);
        if (i == n1m1) ph1(m,k,j,i+1) = g0*(x1r - zmin);
        ph2(m,k,j,i) = phi_c;
        ph3(m,k,j,i) = phi_c;
        if (j == n2m1) ph2(m,k,j+1,i) = phi_c;
        if (k == n3m1) ph3(m,k+1,j,i) = phi_c;
      });
    }
  }

  if (restart) return;

  // --- the random horizontal modes of the velocity seed
  DualArray2D<Real> md("md", std::max(nk,1), 4);   // k2, k3, amplitude, phase
  {
    std::mt19937 rng(kseed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const Real kminr = static_cast<Real>(kmin);
    const Real nkr = static_cast<Real>(kmax - kmin + 1);
    const bool three_d = (indcs.nx3 > 1);
    Real norm = 0.0;
    for (int n=0; n<nk; ++n) {
      md.h_view(n,0) = kminr + std::floor(nkr*u01(rng));
      md.h_view(n,1) = three_d ? (kminr + std::floor(nkr*u01(rng))) : 0.0;
      md.h_view(n,2) = u01(rng) + 0.25;
      md.h_view(n,3) = 2.0*M_PI*u01(rng);
      norm += md.h_view(n,2)*md.h_view(n,2);
    }
    norm = (norm > 0.0) ? 1.0/std::sqrt(norm) : 1.0;
    for (int n=0; n<nk; ++n) md.h_view(n,2) *= norm;
  }
  md.modify_host();
  md.sync_device();
  auto md_d = md.d_view;

  // --- the state, every cell, ghosts included.  The potentials are already built
  // above, on this path and on the restart path alike; only the etotgrav offset that
  // the INITIAL conserved energy carries is added here.
  const Real x2min_m = pmy_mesh_->mesh_size.x2min, x2max_m = pmy_mesh_->mesh_size.x2max;
  const Real x3min_m = pmy_mesh_->mesh_size.x3min, x3max_m = pmy_mesh_->mesh_size.x3max;
  const Real plz = pzhi - pzlo;
  // --- problem/seed_fmode_amp: the LINEAR SURFACE-GRAVITY-WAVE SEED.  One horizontal
  // Fourier mode (m,n) of the box, in the vertical velocity only, with the f-mode's own
  // depth eigenfunction exp(k_h (z - z_top)).  Density and pressure are left alone: the
  // mode sorts its own thermodynamic part out within the first period, at the cost of
  // shedding half the seed into the counter-propagating branch, which is irrelevant to
  // a growth rate measured over periods 1-4.  The amplitude is in units of the sound
  // speed at the TOP of the box, where the eigenfunction peaks.
  fm_kx_ = 2.0*M_PI*fm_m_/(x2max_m - x2min_m);
  fm_ky_ = (indcs.nx3 > 1) ? (2.0*M_PI*fm_n_/(x3max_m - x3min_m)) : 0.0;
  fm_kh_ = std::sqrt(fm_kx_*fm_kx_ + fm_ky_*fm_ky_);
  fm_ztop_ = zmax;  fm_x2min_ = x2min_m;  fm_x3min_ = x3min_m;
  Real cs_top = 0.0;
  {
    int it = static_cast<int>((zmax - zlo)/dzf);
    it = (it < 0) ? 0 : ((it > nfine-1) ? nfine-1 : it);
    const Real dtop = cd.h_view(it), ptop = cp.h_view(it);
    cs_top = std::sqrt(pgen_eos::HostGamma1FromP(eos, dtop, ptop)*ptop/dtop);
  }
  const Real fm_kx = fm_kx_, fm_ky = fm_ky_, fm_kh = fm_kh_;
  const Real fm_v0 = fmamp*cs_top;
  if (global_variable::my_rank == 0 && fmamp > 0.0) {
    std::printf("box_convection: f-MODE SEED (%d,%d)  k_h = %.5e cm^-1  "
                "omega^2 = g k_h = %.5e s^-2  P = %.4f s  amp = %.3e c_s(top) "
                "= %.5e cm/s\n", fm_m_, fm_n_, fm_kh_, g0*fm_kh_,
                2.0*M_PI/std::sqrt(g0*fm_kh_), fmamp, fm_v0);
  }
  par_for("boxconv_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real x2min = size.d_view(m).x2min, x2max = size.d_view(m).x2max;
    const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;
    const Real z = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    const Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
    const Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
    Real s = (z - zlo)/dzf;
    int ii = static_cast<int>(s);
    ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
    const Real f = s - ii;
    const Real d = cd_d(ii)*(1.0 - f) + cd_d(ii+1)*f;
    const Real e = ce_d(ii)*(1.0 - f) + ce_d(ii+1)*f;
    const Real p = cp_d(ii)*(1.0 - f) + cp_d(ii+1)*f;
    Real v1 = 0.0;
    Real efac = 1.0;
    if (vpert > 0.0 && z > pzlo && z < pzhi) {
      Real amp = 0.0;
      for (int n=0; n<nk; ++n) {
        amp += md_d(n,2)*sin(2.0*M_PI*(md_d(n,0)*(x2v - x2min_m)/(x2max_m - x2min_m)
                                     + md_d(n,1)*(x3v - x3min_m)/(x3max_m - x3min_m))
                             + md_d(n,3));
      }
      const Real env = sin(M_PI*(z - pzlo)/plz);
      if (pvar == 0) {
        const Real cs = sqrt(gamma*p/d);
        v1 = vpert*cs*env*amp;
      } else {
        // entropy seed: internal energy (and so the pressure) scaled at fixed density
        efac = 1.0 + vpert*env*amp;
      }
    }
    if (fm_v0 > 0.0 && z > zmin && z < zmax) {
      v1 += fm_v0*exp(fm_kh*(z - zmax))*cos(fm_kx*(x2v - x2min_m))
                                       *cos(fm_ky*(x3v - x3min_m));
    }
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = e*efac + 0.5*d*v1*v1;
    if (etotgrav) u0(m,IEN,k,j,i) += d*g0*(z - zmin);
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvRebuildRadWeights
//! \brief re-form the primitives from the CURRENT conserved state and rebuild the
//! radiation caches Conduction holds, immediately before the two-stream source.  Serves
//! BOTH problem/rt_weights_per_stage (re-centre the caches on the stage state, once per
//! stage) and a radiation sub-step under problem/rt_col3_sub > 1.  ConsToPrim is run over
//! the full range (ghosts included) exactly as Hydro::ConToPrim runs it, so wtemp/wder --
//! the general EOS's temperature guess and derivatives, which BuildRadWeights, the
//! two-stream and the ADI operator all read -- are consistent with w0.  The ghost cells
//! carry the last exchanged u0, i.e. the same staleness w0 already had there; nothing is
//! made worse.  Cost per call: one ConsToPrim, one BuildRadWeights (one serial x1 sweep
//! per column) and, with the ADI / capped transverse operator on, one BuildAngularCoeffs
//! -- all small next to the column solve itself.

void BoxConvRebuildRadWeights(Mesh *pm, Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->phydro == nullptr) return;
  Conduction *pc = pmbp->phydro->pcond;
  if (pc == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  auto &w0 = pmbp->phydro->w0;
  pmbp->phydro->peos->ConsToPrim(pmbp->phydro->u0, w0, false,
                                 0, n1m1, 0, n2m1, 0, n3m1);
  auto &eosd = pmbp->phydro->peos->eos_data;
  if (pc->rad_tau_mode) pc->BuildRadWeights(w0, eosd);
  // the same gate Conduction::AddIsotropicHeatFluxRadiative applies before it calls this
  // (including its 1-D early return: there are no transverse faces and cap_c2/cap_c3 are
  // sized for a single layer, so the builder must not run)
  if (pm->multi_d && (pc->rad_cap_ang > 0.0 || pc->rad_implicit_ang)
      && !pc->rad_sts_split) {
    pc->stage_beta_dt = bdt;
    pc->BuildAngularCoeffs(w0, eosd, bdt);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvTransverseApply
//! \brief one application of the horizontal (transverse ADI / RKL1) radiative operator
//! over `dt`, on the state as it stands, plus the problem/work_hist tag-4 close.  It is
//! exactly the call the rt_split_transverse paths make inline; it exists so that the
//! symmetric pair (problem/rt_pair_sym) can make it twice with a half step.

void BoxConvTransverseApply(Mesh *pm, Real dt) {
  hydro::Hydro *ph = pm->pmb_pack->phydro;
  if (ph == nullptr || ph->pcond == nullptr || !ph->pcond->rad_implicit_ang) return;
  if (ph->pcond->rad_sts_all) {
    ph->pcond->StsConductionUpdate(ph->u0, ph->peos->eos_data, dt);
  } else {
    ph->pcond->ImplicitTransverseUpdate(ph->u0, ph->peos->eos_data, dt);
  }
  BoxConvWorkClose(pm, 4);   // problem/work_hist: the transverse ADI operator
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvSrcs
//! \brief constant gravity along x1 -- in the well-balanced form under
//! wellbalance_dynamic -- plus the top cooling layer.

void BoxConvSrcs(Mesh *pm, Real bdt) {
  // The horizontally averaged profile dump, on the state as it stands BEFORE any source
  // term of this call.  pm->time is the time at the START of the cycle and does not move
  // between stages, so advancing prof_next_ PAST it here is what makes this fire once
  // per cycle rather than once per stage.  A restart arms it at the restart time.
  if (prof_dt_ > 0.0) {
    if (prof_next_ < 0.0) prof_next_ = pm->time;
    if (pm->time >= prof_next_) {
      BoxConvProfileDump(pm);
      while (prof_next_ <= pm->time) prof_next_ += prof_dt_;
    }
  }
  // problem/work_hist, tag 0: this task runs immediately after RKUpdate and nothing
  // between them writes u0 in the active cells, so the interval that closes here is
  // exactly the RK combination plus the flux divergence.
  BoxConvWorkClose(pm, 0);
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  auto eos = pmbp->phydro->peos->eos_data;
  const bool etotgrav = pmbp->phydro->use_etotgrav;
  const bool wbdyn = pmbp->phydro->use_wellbalance_dynamic;
  const bool wbx1 = pmbp->phydro->use_wb_x1;
  const WBOption wbo = pmbp->phydro->wb_option;
  DvceArray4D<Real> phicc = pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = pmbp->phydro->phi0.x1f;
  DvceArray5D<Real> wbq0 = pmbp->phydro->wbq0;
  const Real g0 = g0_, zlo = zlo_, dzf = dzf_;
  const Real zcool = zcool_, zmax = zmax_, tcool = tcool_;
  const int nfine = nfine_;
  auto cd_d = cd_, ce_d = ce_;
  const bool use_cache = wbx1;
  const Real cwid = (zmax > zcool) ? (zmax - zcool) : 1.0;
  const bool cool_on = cool_on_;
  // the top vertical-velocity sponge: it reads the per-plane COLUMN optical depth that
  // Conduction::BuildRadWeights fills (rad_tauf, x1-face centred; the cell's lower face
  // carries the larger tau, which is the conservative -- weaker -- choice for the ramp)
  // (rad_w_built guards the one stage in which the column tau does not exist yet: an
  // all-zero tau would read as "thin everywhere" and damp the whole box)
  const bool vdamp_on = (vdamp_tau_ > 0.0) && (pmbp->phydro->pcond != nullptr)
                        && pmbp->phydro->pcond->rad_w_built;
  const Real vd_hi = vdamp_tau_, vd_lo = vdamp_tau_/3.0;
  const Real vd_rate = vdamp_on ? bdt/vdamp_time_ : 0.0;
  DvceArray4D<Real> vtauf = vdamp_on ? pmbp->phydro->pcond->rad_tauf
                                     : DvceArray4D<Real>("vdamp_unused", 1, 1, 1, 1);

  // ---- problem/rt_budget_verbose: open/close the window, and take the first of the
  // four energy snapshots this call makes.  The window is closed on the FIRST source
  // call of a cycle, which is stage 1, so E is always sampled at the same phase and the
  // accumulators below cover exactly the interval between two samples.
  const bool bud_on = (rtbud_n_ > 0);
  Real bud_e = 0.0, bud_r = 0.0;
  if (bud_on) {
    if (pm->ncycle != rtbud_cyc_) {
      rtbud_cyc_ = pm->ncycle;
      Real e1, r1;
      BoxConvBoxInt(pm, e1, r1);
      if (!rtbud_arm_) {
        rtbud_e0_ = e1; rtbud_r0_ = r1; rtbud_arm_ = true;
      } else if (pm->ncycle % rtbud_n_ == 0) {
        BoxConvBudgetReport(pm, e1, r1);
      }
    }
    BoxConvBoxInt(pm, bud_e, bud_r);
    rtbud_h_[13] += rtbud_fin_*rtbud_area_*bdt;
    rtbud_h_[15] += bdt;
  }

  par_for("boxconv_srcs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real z = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    const Real dzc = (x1max - x1min)/indcs.nx1;
    const Real d = w0(m,IDN,k,j,i);
    Real src = -bdt*g0*d;
    if (!etotgrav) u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
    if (wbdyn) {
      Real pl, pr, d1, d2, d3;
      if (use_cache) {
        WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
      } else {
        hydro::Hydro::getWBq0(eos, wbo, WBVar::wb_pres,
            w0(m,IDN,k,j,i-1), w0(m,IDN,k,j,i), w0(m,IDN,k,j,i+1),
            w0(m,IEN,k,j,i-1), w0(m,IEN,k,j,i), w0(m,IEN,k,j,i+1),
            phicc(m,k,j,i-1), ph1(m,k,j,i), phicc(m,k,j,i), ph1(m,k,j,i+1),
            phicc(m,k,j,i+1), d1, pl, d2, pr, d3);
      }
      // Cartesian: equal face areas, so the background's own pressure drop over dz
      src = bdt*(pr - pl)/dzc;
    }
    u0(m,IM1,k,j,i) += src;
    // the cooling layer: relax the SPECIFIC internal energy toward the initial column's
    if (cool_on && z > zcool) {
      Real s = (z - zcool)/cwid;
      s = (s > 1.0) ? 1.0 : s;
      const Real ramp = s*s*(3.0 - 2.0*s);
      Real t = (z - zlo)/dzf;
      int ii = static_cast<int>(t);
      ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
      const Real f = t - ii;
      const Real d0 = cd_d(ii)*(1.0 - f) + cd_d(ii+1)*f;
      const Real e0 = ce_d(ii)*(1.0 - f) + ce_d(ii+1)*f;
      u0(m,IEN,k,j,i) -= bdt*ramp*d*(w0(m,IEN,k,j,i)/d - e0/d0)/tcool;
    }
    // the top vertical-velocity sponge.  Implicit in the stage, so it cannot overshoot
    // at any dt.  The kinetic energy it removes is taken OUT of the total energy (not
    // converted to internal energy): this is a numerical sponge on the box's acoustic
    // fundamental, and turning its energy into heat would feed exactly the layer it is
    // meant to quiet.
    if (vdamp_on) {
      const Real f = 1.0 - RadBlendWeight(vtauf(m,k,j,i), vd_lo, vd_hi);
      if (f > 0.0) {
        const Real dc = u0(m,IDN,k,j,i);
        const Real m1o = u0(m,IM1,k,j,i);
        const Real m1n = m1o/(1.0 + f*vd_rate);
        u0(m,IM1,k,j,i) = m1n;
        u0(m,IEN,k,j,i) += 0.5*(SQR(m1n) - SQR(m1o))/dc;
      }
    }
  });
  // ---- the sponge's cell range, once, at the first stage it runs (the column tau only
  // exists once BuildRadWeights has run, so it cannot be known at setup)
  if (vdamp_on && !vdamp_printed_) {
    vdamp_printed_ = true;
    Real zlow = std::numeric_limits<Real>::max();
    Kokkos::parallel_reduce("boxconv_vdamp_range",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmin) {
      if (1.0 - RadBlendWeight(vtauf(m,k,j,i), vd_lo, vd_hi) > 0.0) {
        const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
        const Real z = CellCenterX(i-is, indcs.nx1, x1min, x1max);
        lmin = (z < lmin) ? z : lmin;
      }
    }, Kokkos::Min<Real>(zlow));
#if MPI_PARALLEL_ENABLED
    Real zg = zlow;
    MPI_Allreduce(&zg, &zlow, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      const Real gx1min = pm->mesh_size.x1min, gx1max = pm->mesh_size.x1max;
      const int gnx1 = pm->mesh_indcs.nx1;
      if (zlow < gx1max) {
        const int i0 = static_cast<int>((zlow - gx1min)/((gx1max - gx1min)/gnx1));
        std::printf("  top v1 sponge active: i = %d .. %d (of %d), z >= %.5e,"
                    " f = 1 at tau <= %.4e, timescale %.4e s\n",
                    i0, gnx1-1, gnx1, zlow, vd_lo, vdamp_time_);
      } else {
        std::printf("  top v1 sponge: NO cell has tau < %.4e -- the sponge is inert\n",
                    vd_hi);
      }
    }
  }

  // ---- THE BOTTOM SPONGE (problem/vdamp_bot_cells; see the note on vdb_cells_).  It
  // runs HERE, on the u0 the kernel above left, and BEFORE the bottom-inflow source
  // further down, so the momentum bc_mode_bot = 5 injects into the wall cell is never
  // cancelled by the sponge in the same stage.
  if (vdb_cells_ > 0) {
    const int nb = vdb_cells_;
    const int gnx2 = pm->mesh_indcs.nx2, gnx3 = pm->mesh_indcs.nx3;
    if (vdb_d_.extent_int(0) != nb) {
      Kokkos::realloc(vdb_d_, nb);
      Kokkos::realloc(vdb_h_, nb);
    }
    auto vdb = vdb_d_;
    if (vdb_mean_) {
      // one team per x1 index of the layer sums v1 over that plane's (m,k,j); the plane
      // SUMS are Allreduced (every rank applies the mean) and divided by the global
      // plane cell count, exactly as the profile dump does
      const int lnx2 = indcs.nx2, lnx3 = indcs.nx3;
      const int nkj = (nmb1+1)*lnx3*lnx2;
      Kokkos::TeamPolicy<> vpol(DevExeSpace(), nb, Kokkos::AUTO);
      Kokkos::parallel_for("boxconv_vdbot_mean", vpol,
      KOKKOS_LAMBDA(Kokkos::TeamPolicy<>::member_type tmember) {
        const int i = is + tmember.league_rank();
        Real vs = 0.0;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkj),
        [&](const int idx, Real &ls) {
          const int m = idx/(lnx3*lnx2);
          const int kj = idx - m*(lnx3*lnx2);
          const int k = ks + kj/lnx2;
          const int j = js + (kj - (kj/lnx2)*lnx2);
          ls += u0(m,IM1,k,j,i)/u0(m,IDN,k,j,i);
        }, Kokkos::Sum<Real>(vs));
        Kokkos::single(Kokkos::PerTeam(tmember), [&]() { vdb(i-is) = vs; });
      });
      Kokkos::fence();
      Kokkos::deep_copy(vdb_h_, vdb_d_);
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, vdb_h_.data(), nb, MPI_ATHENA_REAL, MPI_SUM,
                    MPI_COMM_WORLD);
#endif
      const Real fpl = 1.0/static_cast<Real>(gnx2*gnx3);
      for (int q=0; q<nb; ++q) vdb_h_(q) *= fpl;
      Kokkos::deep_copy(vdb_d_, vdb_h_);
    }
    const Real vb_rate = bdt/vdb_time_;
    const bool vb_mean = vdb_mean_;
    par_for("boxconv_vdbot", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, is+nb-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real f = 0.5*(1.0 + std::cos(M_PI*static_cast<Real>(i-is)/nb));
      const Real gg = 1.0 - std::exp(-f*vb_rate);
      const Real dc = u0(m,IDN,k,j,i);
      const Real m1o = u0(m,IM1,k,j,i);
      const Real m1n = m1o - gg*(vb_mean ? dc*vdb(i-is) : m1o);
      u0(m,IM1,k,j,i) = m1n;
      u0(m,IEN,k,j,i) += 0.5*(SQR(m1n) - SQR(m1o))/dc;
    });
  }

  // the gravity/WB/cooling kernel's own contribution, and the two x1 wall faces of the
  // RK flux divergence (every other face telescopes out of the box integral)
  if (bud_on) {
    Real e2, r2;
    BoxConvBoxInt(pm, e2, r2);
    rtbud_h_[8] += e2 - bud_e;
    bud_e = e2;
    auto &flx1b = pmbp->phydro->uflx.x1f;
    auto budf = rtbud_;
    par_for("boxconv_budflx", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real da = size.d_view(m).dx2*size.d_view(m).dx3;
      Kokkos::atomic_add(&budf(0), bdt*da*flx1b(m,IEN,k,j,is));
      Kokkos::atomic_add(&budf(1), bdt*da*flx1b(m,IEN,k,j,ie+1));
      Kokkos::atomic_add(&budf(2), bdt*da*flx1b(m,IDN,k,j,is));
      Kokkos::atomic_add(&budf(3), bdt*da*flx1b(m,IDN,k,j,ie+1));
    });
  }
  // --- THE WALLS ARE IMPERMEABLE.  The bc_mode-3 ghost above is the hydrostatic
  // continuation of the evolved interior, which makes the wall-face mass flux small; but
  // small and one-signed still integrates into a leak over 1e5 stages.  So cancel it
  // exactly, cell by cell: undo the contribution the wall face made to this stage's flux
  // divergence in the mass and energy channels, and leave the MOMENTUM channel alone --
  // its pressure term is the wall force that holds the box up.  This is what
  // red_giant.cpp's problem/wall_noflux does at its inner wall.
  if (wall_noflux_) {
    auto &flx1w = pmbp->phydro->uflx.x1f;
    auto &mb_bcs = pmbp->pmb->mb_bcs;
    const bool difflx = diff_flux_;
    const bool topopen = (bc_mode_top_ == 4);
    const bool bud_w = bud_on;
    auto budw = bud_on ? rtbud_ : DvceArray1D<Real>("budwdummy", kNBud);
    par_for("boxconv_wallflux", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real idz = indcs.nx1/(size.d_view(m).x1max - size.d_view(m).x1min);
      for (int w=0; w<2; ++w) {
        const bool inner = (w == 0);
        const BoundaryFlag bf = inner ? mb_bcs.d_view(m,BoundaryFace::inner_x1)
                                      : mb_bcs.d_view(m,BoundaryFace::outer_x1);
        if (bf != BoundaryFlag::user) continue;
        // an OPEN top (bc_mode_top = 4) is meant to carry a flux: do not cancel it
        if (!inner && topopen) continue;
        const int ic = inner ? is : ie;             // the cell against the wall
        const int ifc = inner ? is : (ie + 1);      // the wall face itself
        // RKUpdate did u0 -= bdt*(flx(ie+1) - flx(is))/dz, so the inner face entered with
        // a + sign and the outer face with a -.
        const Real sgn = inner ? 1.0 : -1.0;
        const Real dm = sgn*bdt*idz*flx1w(m,IDN,k,j,ifc);
        if (dm == 0.0) continue;
        const Real d = u0(m,IDN,k,j,ic);
        if (!(d - dm > 0.0)) continue;
        Real de;
        if (difflx) {
          // conduction (and viscosity) have already been added into this face's ENERGY
          // channel, and at the bottom wall that term IS the imposed luminosity.  Remove
          // only the energy the cancelled mass carried, exactly as red_giant.cpp does.
          Real ei = u0(m,IEN,k,j,ic)
                    - 0.5*(SQR(u0(m,IM1,k,j,ic)) + SQR(u0(m,IM2,k,j,ic))
                           + SQR(u0(m,IM3,k,j,ic)))/d;
          if (etotgrav) ei -= d*phicc(m,k,j,ic);
          // guard the EOS call: a tabulated EOS takes log10(e) and would write a NaN
          // into u0 with no precursor.  Skipping one stage's cancellation is harmless.
          if (!(ei > 0.0)) continue;
          de = dm*(u0(m,IEN,k,j,ic) + eos.Pressure(d, ei))/d;
        } else {
          de = sgn*bdt*idz*flx1w(m,IEN,k,j,ifc);
        }
        u0(m,IDN,k,j,ic) -= dm;
        u0(m,IEN,k,j,ic) -= de;
        if (bud_w) {
          const Real dvw = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
          Kokkos::atomic_add(&budw(inner ? 4 : 5), de*dvw);
          Kokkos::atomic_add(&budw(inner ? 6 : 7), dm*dvw);
        }
      }
    });
  }

  // --- THE BOTTOM INFLOW (bc_mode_bot = 5).  The wall-face cancellation above STAYS:
  // an uncancelled wall face in a stratified box is not a small leak but a drain (the
  // He column ran down to 0.32 M0 in 50 s when the bottom face was left to its own
  // Riemann flux -- ~300x the mass loss the open top makes).  So the wall stays
  // impermeable and the inflow is added here as an EXACT source in the cell against the
  // bottom wall: per unit area and time it is rho_base v_in of mass, rho_base v_in^2 of
  // momentum and rho_base v_in h_base of energy, with h_base the RELAXED BASE STATE's
  // specific enthalpy (e + p)/rho + v_in^2/2 from the initial column at x1min.  The box
  // mass then follows the controller exactly, and the imposed luminosity -- which lives
  // in the same face's ENERGY channel, or under rt_bottom_flux in the two-stream's own
  // lower boundary -- is untouched.
  if (bc_mode_bot_ == 5 && inflow_vin_ > 0.0) {
    auto &mb_bcs_i = pmbp->pmb->mb_bcs;
    const Real vin = inflow_vin_;
    const Real db = inflow_rhob_, eb = inflow_eb_, pb = inflow_pb_;
    const Real hb = (db > 0.0) ? ((eb + pb)/db + 0.5*vin*vin) : 0.0;
    const Real zmn = zmin_;
    par_for("boxconv_inflow", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      if (mb_bcs_i.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
      const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
      const Real idz = indcs.nx1/(x1max - x1min);
      const Real dm = db*vin*bdt*idz;
      Real de = dm*hb;
      if (etotgrav) de += dm*g0*(CellCenterX(0, indcs.nx1, x1min, x1max) - zmn);
      u0(m,IDN,k,j,is) += dm;
      u0(m,IM1,k,j,is) += dm*vin;
      u0(m,IEN,k,j,is) += de;
    });
  }

  if (bud_on) {
    Real e3, r3;
    BoxConvBoxInt(pm, e3, r3);
    rtbud_h_[14] += e3 - bud_e;
    bud_e = e3;
  }

  // problem/work_hist, tag 1: everything above -- the well-balanced gravity source, the
  // cooling layer, the top sponge, the wall-flux cancellation and the bottom inflow.
  BoxConvWorkClose(pm, 1);

  // --- the grey two-stream, after gravity and the cooling layer, exactly where
  // red_giant.cpp calls it: inside the stage, on the state the last ConToPrim left.
  // problem/rt_before_flux: the whole block has already run at the head of this stage,
  // in Hydro::RTBeforeFlux (BoxConvRTBeforeFlux), on the stage-start state.
  if (rt_on_ && !rt_strang_ && !rt_once_ && !rt_col3_once_ && !rt_imex_ &&
      !rt_before_flux_) {
    // problem/rt_col3_sub: N applications of the WHOLE radiation operator per stage with
    // beta_dt/N each, at the unchanged hydro dt.  N = 1 is one pass of exactly the call
    // that was here, with no rebuild and no transverse call -- bitwise the old code.
    const int nsub = rt_col3_sub_;
    const Real sdt = bdt/static_cast<Real>(nsub);
    for (int isub=0; isub<nsub; ++isub) {
      // problem/rt_pair_sym = 1: the frozen ADI face conductances are formed with the
      // HALF step and used for BOTH halves of the symmetric pair, which is what makes
      // the pair adjoint-symmetric; the rebuild also refreshes w0 (and so the tau/blend
      // weights) for the column solve that follows, exactly as the sub-cycle does.
      if (rt_pair_sym_ == 1) {
        BoxConvRebuildRadWeights(pm, 0.5*sdt);
        BoxConvTransverseApply(pm, 0.5*sdt);
      } else if (nsub > 1 || rtwps_) {
        BoxConvRebuildRadWeights(pm, sdt);
      }
      two_stream_rt::picket_fence_two_stream_RT(pm, sdt);
      // problem/work_hist, tag 3: the column solve closed its own interval from inside
      // (two_stream_rt::rt_probe, tag 2), so what closes here is the radiative momentum
      // force and nothing else.
      BoxConvWorkClose(pm, 3);
      // the horizontal ADI operator travels WITH the column, over the same sub-step and
      // on the state the column solve has just relaxed (the ordering rt_split_transverse
      // established).  The in-stage transverse task is a no-op under rad_tr_split_out.
      if (rt_pair_sym_ == 1) {
        // the closing half of the symmetric pair, on the state the column just relaxed
        BoxConvTransverseApply(pm, 0.5*sdt);
      } else if (nsub > 1 && rt_split_tr_) {
        hydro::Hydro *ph = pm->pmb_pack->phydro;
        if (ph != nullptr && ph->pcond != nullptr && ph->pcond->rad_implicit_ang) {
          if (ph->pcond->rad_sts_all) {
            ph->pcond->StsConductionUpdate(ph->u0, ph->peos->eos_data, sdt);
          } else {
            ph->pcond->ImplicitTransverseUpdate(ph->u0, ph->peos->eos_data, sdt);
          }
          BoxConvWorkClose(pm, 4);  // problem/work_hist: the transverse ADI operator
        }
      }
      if (bud_on) {
        Real e4, r4, ft, fc;
        BoxConvBoxInt(pm, e4, r4);
        rtbud_h_[9] += e4 - bud_e;
        bud_e = e4;
        BoxConvFtopInt(pm, ft, fc);
        rtbud_h_[11] += sdt*ft;
        rtbud_h_[12] += sdt*fc;
      }
    }
    // The per-column surface dump, on the flux that call just wrote.  pm->time is the
    // time at the START of the cycle and does not move between stages, so advancing
    // surf_next_ PAST it here is what makes this fire once per cycle rather than once
    // per stage.  A restart arms it at the restart time, so the series simply resumes.
    if (surf_dt_ > 0.0 && two_stream_rt::rt_face_flux_ready()) {
      if (surf_next_ < 0.0) surf_next_ = pm->time;
      if (pm->time >= surf_next_) {
        BoxConvSurfaceDump(pm);
        while (surf_next_ <= pm->time) surf_next_ += surf_dt_;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvRTSplit
//! \brief problem/rt_strang: the grey two-stream as a Strang-split operator around the
//! time integrator.  Hydro::RTStrangSplit calls this twice a cycle with bdt = dt/2, on
//! the state at the start and at the end of the cycle, and runs ConToPrim after each.
//! Nothing else moves: this is exactly the call BoxConvSrcs makes when rt_strang is
//! off, with a different dt and at a different point in the cycle.

void BoxConvRTSplit(Mesh *pm, Real bdt) {
  if (!rt_on_) return;
  // the same budget bookkeeping the in-stage call does, so the energy budget and the
  // Ftop/Fcut integrals stay complete when the source moves out of the stage
  const bool bud_on = (rtbud_n_ > 0);
  Real e0 = 0.0, r0 = 0.0;
  if (bud_on) BoxConvBoxInt(pm, e0, r0);
  if (rtwps_) BoxConvRebuildRadWeights(pm, bdt);
  two_stream_rt::picket_fence_two_stream_RT(pm, bdt);
  BoxConvWorkClose(pm, 3);   // problem/work_hist: the radiative force (see tag 3 above)
  // problem/rt_split_transverse: the horizontal ADI operator moves WITH the column, over
  // the same bdt and on the state the column solve has just relaxed.  Its x2/x3 ghosts
  // are the last exchange's, exactly as they are for the in-stage task it replaces
  // (imptrc runs before SendU), so nothing is more stale than before.
  if (rt_split_tr_) {
    hydro::Hydro *ph = pm->pmb_pack->phydro;
    if (ph != nullptr && ph->pcond != nullptr && ph->pcond->rad_implicit_ang) {
      if (ph->pcond->rad_sts_all) {
        ph->pcond->StsConductionUpdate(ph->u0, ph->peos->eos_data, bdt);
      } else {
        ph->pcond->ImplicitTransverseUpdate(ph->u0, ph->peos->eos_data, bdt);
      }
      BoxConvWorkClose(pm, 4);   // problem/work_hist: the transverse ADI operator
    }
  }
  if (bud_on) {
    Real e4, r4, ft, fc;
    BoxConvBoxInt(pm, e4, r4);
    rtbud_h_[9] += e4 - e0;
    BoxConvFtopInt(pm, ft, fc);
    rtbud_h_[11] += bdt*ft;
    rtbud_h_[12] += bdt*fc;
  }
  if (surf_dt_ > 0.0 && two_stream_rt::rt_face_flux_ready()) {
    if (surf_next_ < 0.0) surf_next_ = pm->time;
    if (pm->time >= surf_next_) {
      BoxConvSurfaceDump(pm);
      while (surf_next_ <= pm->time) surf_next_ += surf_dt_;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvRTBeforeFlux
//! \brief problem/rt_before_flux: the WHOLE radiation operator at the HEAD of the RK
//! stage, on the stage-start state, BEFORE the hydro flux update -- the reversed Lie
//! order (see the declaration of rt_before_flux_).  Enrolled as
//! ProblemGenerator::user_rt_before_flux; Hydro::RTBeforeFlux calls it once per stage
//! with that stage's beta_dt and follows it with a ghost/BC update and a ConToPrim.
//!
//! The sequence is EXACTLY the one BoxConvSrcs runs today at rt_col3_sub = 1,
//! rt_pair_sym = 0 -- the column solve + radiative force, then the horizontal ADI over
//! the same step -- with the same budget bookkeeping and the same per-column surface
//! dump, so nothing but the position in the stage changes.

void BoxConvRTBeforeFlux(Mesh *pm, Real bdt) {
  if (!rt_on_ || !rt_before_flux_) return;
  const bool bud_on = (rtbud_n_ > 0);
  Real e0 = 0.0, r0 = 0.0;
  if (bud_on) BoxConvBoxInt(pm, e0, r0);
  // The tau/blend weights and the frozen ADI face conductances are normally built inside
  // Hydro::Fluxes, from the STAGE-START primitives -- which have not been built yet when
  // this runs (and do not exist at all in the first stage of a run).  Build them here,
  // from the same stage-start state, so the operator sees exactly the coefficients the
  // in-stage ordering gives it: this is the rebuild rt_col3_sub / rt_pair_sym make.
  BoxConvRebuildRadWeights(pm, bdt);
  two_stream_rt::picket_fence_two_stream_RT(pm, bdt);
  // problem/work_hist, tag 3: the column solve closed its own interval from inside
  // (two_stream_rt::rt_probe, tag 2), so what closes here is the radiative force
  BoxConvWorkClose(pm, 3);
  // the horizontal ADI operator travels WITH the column (rt_split_transverse is forced
  // on, so the in-stage transverse task is a no-op); tag 4 closes inside the call
  BoxConvTransverseApply(pm, bdt);
  if (bud_on) {
    Real e4, r4, ft, fc;
    BoxConvBoxInt(pm, e4, r4);
    rtbud_h_[9] += e4 - e0;
    BoxConvFtopInt(pm, ft, fc);
    rtbud_h_[11] += bdt*ft;
    rtbud_h_[12] += bdt*fc;
  }
  // the per-column surface dump, once per cycle (pm->time does not move between stages)
  if (surf_dt_ > 0.0 && two_stream_rt::rt_face_flux_ready()) {
    if (surf_next_ < 0.0) surf_next_ = pm->time;
    if (pm->time >= surf_next_) {
      BoxConvSurfaceDump(pm);
      while (surf_next_ <= pm->time) surf_next_ += surf_dt_;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvRTImEx
//! \brief problem/rt_imex: the mode-3 column solve AS THE IMPLICIT STAGE OPERATOR of the
//! ImEx-RK integrator (<time>/integrator = imex2 or imex2+).  Enrolled as
//! ProblemGenerator::user_imex_func and called by Hydro::RTImEx (once per explicit
//! stage, in the place the in-stage source it replaces occupied) and by
//! Hydro::RTImExFirst (with estage = -1 and 0, at the head of stage 1, for the extra
//! fully implicit stages the tableau adds).
//!
//! It mirrors IonNeutral::ImpRKUpdate exactly -- that is the only true ImEx user in the
//! code, and its conventions are what the Driver's weights are written for:
//!   istage = estage + 2         the implicit stage number, 1..nimp_stages+1
//!   slot   = istage - 1         where this stage's source S^(l) is stored
//!   row    = istage - 2         the a_twid row used to recombine the stored sources
//! and the three steps per stage are
//!   (a) u0 += dt*sum_{l<=row} a_twid[row][l]*S^(l)   the earlier sources, EXPLICITLY
//!   (b) the implicit solve on that state with the effective step a_impl*dt
//!   (c) S^(slot) = (u0_after - u0_before)/(a_impl*dt) = R(U^(istage))
//! (b) and (c) are skipped on the LAST explicit stage (estage == nexp_stages), which is
//! the final combination and does (a) only.  a_twid is NOT the implicit Butcher tableau
//! and the diagonal is NOT a_twid[k][k]: the solve always uses a_impl (0.5 for imex2,
//! 1+1/sqrt(2) for imex2+), exactly as ion-neutral does.
//!
//! Because the column solve is backward Euler over a_impl*dt (mode 3 solves the cell
//! energies implicitly), the increment it returns divided by a_impl*dt IS the source
//! evaluated on the state the solve produced, which is what (c) needs.
//!
//! WHAT IS IN THE OPERATOR: the energy exchange, the radiative momentum force
//! (problem/rt_rad_force) and its work -- the whole of what picket_fence_two_stream_RT
//! writes into u0 -- hence the four stored components IM1,IM2,IM3,IEN.  The force is
//! therefore evaluated on THAT stage's fresh column solve (it is computed in the same
//! kernel that applies the energy source) and recombined with the same weights.  The
//! horizontal ADI operator (<hydro>/rad_implicit_ang) is NOT part of it: it stays the
//! separate operator-split task it is today, right after this one.  It is the obvious
//! next thing to fold in if a residual O(dt) mode survives.
//!
//! RESTARTS need no new state: slot l is written at implicit stage l+1 and first read at
//! stage l+2 of the SAME cycle, so the array is rebuilt from scratch every cycle.

void BoxConvRTImEx(Mesh *pm, Driver *pd, const int estage) {
  if (!rt_on_ || !rt_imex_) return;
  const int istage = estage + 2;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->phydro->u0;
  const Real dt = pm->dt;

  // allocated here, not in UserProblem: the Driver (which owns nimp_stages) is built
  // after the problem generator.  Zeroed, so the coefficient-zero slots imex2+ never
  // writes can still be read.
  if (rtimex_src_.extent(0) == 0) {
    const int nmb = std::max(pmbp->nmb_thispack, pm->nmb_maxperrank);
    const int n1 = indcs.nx1 + 2*indcs.ng;
    const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
    const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
    Kokkos::realloc(rtimex_src_, pd->nimp_stages, nmb, 4, n3, n2, n1);
    Kokkos::deep_copy(rtimex_src_, 0.0);
  }

  // imex2+ (Krapp et al. 2024) has all-zero a_twid rows 0 and 1 and no implicit
  // contribution at its first two implicit stages: they are NO-OPS, exactly as
  // IonNeutral::ImpRKUpdate zeroes its coefficients there.  The practical gain is that
  // imex2+ needs no pre-stage at all, so nothing is ever solved before the stage-1 flux
  // divergence and no ghost zone is ever a stage stale (see the note in RTImExFirst).
  const bool noop = (pd->integrator == "imex2+") && (istage < 3);
  if (noop) return;

  // ---- (a) the earlier stages' sources, re-applied EXPLICITLY --------------------
  if (istage > 1) {
    const int row = istage - 2;
    Real wgt[4];
    for (int l=0; l<4; ++l) wgt[l] = (l <= row) ? (pd->a_twid[row][l])*dt : 0.0;
    const Real w0c = wgt[0], w1c = wgt[1], w2c = wgt[2], w3c = wgt[3];
    const int nl = row;
    auto s_ = rtimex_src_;
    par_for("boxconv_imex_exp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real d[4] = {0.0, 0.0, 0.0, 0.0};
      const Real ww[4] = {w0c, w1c, w2c, w3c};
      for (int l=0; l<=nl; ++l) {
        for (int c=0; c<4; ++c) d[c] += ww[l]*s_(l,m,c,k,j,i);
      }
      u0(m,IM1,k,j,i) += d[0];
      u0(m,IM2,k,j,i) += d[1];
      u0(m,IM3,k,j,i) += d[2];
      u0(m,IEN,k,j,i) += d[3];
    });
  }

  // ---- (b) the implicit solve, and (c) the source it defines ---------------------
  // Skipped on the last explicit stage, which only combines what is stored.
  if (estage < pd->nexp_stages) {
    const int sl = istage - 1;
    const Real adt = (pd->a_impl)*dt;
    auto s_ = rtimex_src_;
    // stash -u0 in the slot: the increment is formed in place, no second array
    par_for("boxconv_imex_pre", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      s_(sl,m,0,k,j,i) = -u0(m,IM1,k,j,i);
      s_(sl,m,1,k,j,i) = -u0(m,IM2,k,j,i);
      s_(sl,m,2,k,j,i) = -u0(m,IM3,k,j,i);
      s_(sl,m,3,k,j,i) = -u0(m,IEN,k,j,i);
    });
    // the solve.  BoxConvRTSplit is the bare call plus the budget/Ftop/surface-dump
    // bookkeeping the in-stage call also does; nothing in it depends on WHERE it is
    // called from, only on the dt it is handed.
    BoxConvRTSplit(pm, adt);
    const Real iadt = 1.0/adt;
    par_for("boxconv_imex_rec", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      s_(sl,m,0,k,j,i) = (s_(sl,m,0,k,j,i) + u0(m,IM1,k,j,i))*iadt;
      s_(sl,m,1,k,j,i) = (s_(sl,m,1,k,j,i) + u0(m,IM2,k,j,i))*iadt;
      s_(sl,m,2,k,j,i) = (s_(sl,m,2,k,j,i) + u0(m,IM3,k,j,i))*iadt;
      s_(sl,m,3,k,j,i) = (s_(sl,m,3,k,j,i) + u0(m,IEN,k,j,i))*iadt;
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvHistory
//! \brief the emergent-flux history columns.  See the header block for what each is and
//! why the MEAN SQUARE, not the rms, is what gets written.
//!
//! Nothing here re-runs the solver: two_stream_rt::rt_face_flux() is the solver's own
//! face-flux array, which lives for the run, so these are the numbers of the LAST RT
//! call of the cycle -- which is what the history is called after.

void BoxConvHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 5 + (fm_hist_ ? 2 : 0) + (work_on_ ? 9 : 0);
  pdata->label[0] = "Ftop";
  pdata->label[1] = "Ftop2";
  pdata->label[2] = "Fcut";
  pdata->label[3] = "Ttop";
  pdata->label[4] = "Ttop2";
  if (fm_hist_) {
    pdata->label[5] = "fmAc";
    pdata->label[6] = "fmAs";
  }
  if (work_on_) {
    pdata->label[7]  = "Wflx";   // cumulative mode-KE change across RKUpdate
    pdata->label[8]  = "Wgrv";   // ... across the gravity/WB + layer sources
    pdata->label[9]  = "Wfrc";   // ... across the radiative momentum force
    pdata->label[10] = "Woth";   // ... across the column + ADI: the closure control
    pdata->label[11] = "Wtco";   // the column heating's mode-projected pdV work
    pdata->label[12] = "Wtfr";   // the radiative force's v.f heating work, same measure
    pdata->label[13] = "Wtad";   // the transverse ADI operator's, same measure
    pdata->label[14] = "Emod";   // the mode kinetic energy (instantaneous)
    pdata->label[15] = "Sdsp";   // mode amplitude of v1 in the top row (cm/s)
  }
  for (int n=0; n<pdata->nhist; ++n) pdata->hdata[n] = 0.0;

  // --- problem/fmode_hist: the (m,n) amplitude of the vertical velocity, projected
  // with the f-mode's own depth weight exp(k_h (z - z_top)).  fmAc is the cos.cos
  // projection (the phase the seed is written in) and fmAs the sin.sin one, which the
  // seed leaves at zero and which is therefore a free control on how much of the signal
  // is the seeded mode and how much is everything else at the same |k|.  Both are box
  // means in cm/s: the 4/N normalisation makes fmAc equal the weighted-mean amplitude of
  // a pure cos.cos mode.  The cells are uniform on this mesh, so a cell count IS the
  // volume weight.
  if (fm_hist_) {
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &indcs = pm->mb_indcs;
    const int is = indcs.is, js = indcs.js, ks = indcs.ks;
    const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
    const int ncell = pmbp->nmb_thispack*nx3*nx2*nx1;
    auto &size = pmbp->pmb->mb_size;
    auto &w0 = pmbp->phydro->w0;
    const Real kx = fm_kx_, ky = fm_ky_, kh = fm_kh_;
    const Real ztop = fm_ztop_, x2m = fm_x2min_, x3m = fm_x3min_;
    const Real inc = 4.0/(static_cast<Real>(pm->mesh_indcs.nx1)*
                          static_cast<Real>(pm->mesh_indcs.nx2)*
                          static_cast<Real>(pm->mesh_indcs.nx3));
    array_sum::GlobalSum sum_fm;
    Kokkos::parallel_reduce("boxconv_fmhist",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, ncell),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &mb_sum) {
      const int m = idx/(nx3*nx2*nx1);
      const int r = idx - m*(nx3*nx2*nx1);
      const int k = ks + r/(nx2*nx1);
      const int r2 = r - (r/(nx2*nx1))*(nx2*nx1);
      const int j = js + r2/nx1;
      const int i = is + r2 - (r2/nx1)*nx1;
      const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
      const Real x2min = size.d_view(m).x2min, x2max = size.d_view(m).x2max;
      const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;
      const Real z = CellCenterX(i-is, nx1, x1min, x1max);
      const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      const Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      const Real wz = exp(kh*(z - ztop))*w0(m,IVX,k,j,i);
      array_sum::GlobalSum hvars;
      for (int n=0; n<NHISTORY_VARIABLES; ++n) hvars.the_array[n] = 0.0;
      hvars.the_array[5] = inc*wz*cos(kx*(x2v - x2m))*cos(ky*(x3v - x3m));
      hvars.the_array[6] = inc*wz*sin(kx*(x2v - x2m))*sin(ky*(x3v - x3m));
      mb_sum += hvars;
    }, Kokkos::Sum<array_sum::GlobalSum>(sum_fm));
    Kokkos::fence();
    pdata->hdata[5] = sum_fm.the_array[5];
    pdata->hdata[6] = sum_fm.the_array[6];
  }
  // problem/work_hist: the accumulators live on rank 0 (the plane sums they are built
  // from are MPI-reduced there), so only rank 0 contributes and the history's own
  // MPI_SUM over ranks lands on the number itself.
  if (work_on_ && global_variable::my_rank == 0) {
    for (int q=0; q<7; ++q) pdata->hdata[7+q] = wk_acc_[q];
    pdata->hdata[14] = wk_emod_;
    pdata->hdata[15] = wk_sdsp_;
  }
  if (!rt_on_ || !two_stream_rt::rt_face_flux_ready()) return;

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int ncol = pmbp->nmb_thispack*nx3*nx2;
  auto fb = two_stream_rt::rt_face_flux();
  auto icut = two_stream_rt::rt_cut_index();
  const int nblk = two_stream_rt::rt_face_nblk();
  auto &w0 = pmbp->phydro->w0;
  auto wt = pmbp->phydro->wtemp;
  auto eos = pmbp->phydro->peos->eos_data;
  const bool gen = eos.IsGeneral();
  const Real tcgs = eos.temp_cgs;
  const Real rgas = rgas_;
  // Each column contributes 1/N of the box mean, so the MPI_SUM the history performs
  // over ranks lands on the mean itself.  N is the whole mesh's column count: one
  // MeshBlock spans the whole x1 extent, so nx2*nx3 of the MESH counts every column once.
  const Real inc = 1.0/(static_cast<Real>(pm->mesh_indcs.nx2)*
                        static_cast<Real>(pm->mesh_indcs.nx3));
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("boxconv_hist",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, ncol),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &mb_sum) {
    const int m = idx/(nx3*nx2);
    const int kj = idx - m*(nx3*nx2);
    const int k = ks + kj/nx2;
    const int j = js + kj - (kj/nx2)*nx2;
    // the sweep starts at icut, the deepest face the two-stream integrates; every face
    // below it is left at zero, so the wall face i = is carries nothing to read
    const int ic = icut(m,k,j);
    Real ft = 0.0, fcut = 0.0;
    for (int b=0; b<nblk; ++b) {
      ft += fb(m,b,ie+1,k,j);
      fcut += fb(m,b,ic,k,j);
    }
    // the temperature ConsToPrim already solved for this very cell; no second inversion
    const Real tk = (gen) ? wt(m,k,j,ie)*tcgs
                          : w0(m,IEN,k,j,ie)/(rgas*w0(m,IDN,k,j,ie));
    array_sum::GlobalSum hvars;
    hvars.the_array[0] = inc*ft;
    hvars.the_array[1] = inc*ft*ft;
    hvars.the_array[2] = inc*fcut;
    hvars.the_array[3] = inc*tk;
    hvars.the_array[4] = inc*tk*tk;
    for (int n=5; n<NHISTORY_VARIABLES; ++n) hvars.the_array[n] = 0.0;
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();
  // slots 0..4 ONLY: the reduction above zeroes everything above 4, and the f-mode and
  // work columns were already written above.  (It used to run to pdata->nhist, which
  // silently wiped fmAc/fmAs in any run that had BOTH radiation and fmode_hist on --
  // arm set L never saw it because it ran with rt_two_stream = false, which returns
  // before this block.)
  for (int n=0; n<5; ++n) pdata->hdata[n] = sum_this_mb.the_array[n];
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvBC
//! \brief the initial column continued into the x1 ghosts with the normal velocity
//! mirrored: a reflecting wall whose ghost state is hydrostatic, so the well-balanced
//! stencil of the wall cell sees a consistent background.  The bottom wall's conductive
//! flux is REPLACED by <hydro>/rad_flux_inner, so the ghost temperature there is inert.

void BoxConvBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  DvceArray4D<Real> phicc = pmbp->phydro->phicc0;
  const Real g0 = g0_, zlo = zlo_, dzf = dzf_, zmin = zmin_;
  const int nfine = nfine_;
  const bool etotgrav = etotgrav_;
  auto cd_d = cd_, ce_d = ce_;
  // the TOP wall may run a different mode; -1 means "the same as bc_mode"
  const int bcm_top = (bc_mode_top_ >= 0) ? bc_mode_top_ : bc_mode_;
  const int bcm_bot = (bc_mode_bot_ == 5) ? (inflow_ghost_ ? 5 : bc_mode_)
                    : ((bc_mode_bot_ >= 0) ? bc_mode_bot_ : bc_mode_);
  auto ct_d = ct_;
  // --- THE BOTTOM-INFLOW CONTROLLER.  One global reduction per CYCLE (not per stage:
  // this function is called once per stage on every rank, and pm->ncycle does not move
  // between the stages of a cycle), so v_in is a constant of the cycle.  Every rank
  // reaches this the same number of times, which is what the Allreduce needs.
  // gated on bc_mode_bot_, NOT on bcm_bot: under bc_inflow_ghost = false the ghost
  // follows bc_mode and bcm_bot is not 5, but the controller must still run
  if (bc_mode_bot_ == 5) {
    if (pm->ncycle != inflow_cyc_) {
      inflow_cyc_ = pm->ncycle;
      inflow_mnow_ = BoxConvMassInt(pm);
      Real vin_new = (inflow_mass0_ - inflow_mnow_)
                     /(inflow_rhob_*inflow_area_*inflow_tau_);
      if (!(vin_new > 0.0)) vin_new = 0.0;
      if (vin_new > inflow_vmax_) vin_new = inflow_vmax_;
      inflow_vin_ = vin_new;
      if (inflow_print_n_ > 0 && global_variable::my_rank == 0
          && (pm->ncycle % inflow_print_n_ == 0)) {
        std::printf("boxconv inflow: cycle %d  t %.6e  v_in %.6e  M/M0 %.8f\n",
                    pm->ncycle, pm->time, inflow_vin_,
                    (inflow_mass0_ > 0.0) ? inflow_mnow_/inflow_mass0_ : 0.0);
        std::fflush(stdout);
      }
    }
  }
  auto eos = pmbp->phydro->peos->eos_data;
  const WBOption wbo = pmbp->phydro->wb_option;
  const Real wfac = wall_walk_maxfac_;
  const Real dfl = eos.dfloor;
  // WHERE THE GHOST FILL READS THE INTERIOR CELL IT CONTINUES: from u0, never from w0.
  // Two reasons, and the second is fatal.
  // (1) ApplyPhysicalBCs runs AFTER RKUpdate, the user source terms (gravity, the
  // cooling layer, the two-stream) and the implicit x1 conduction, and BEFORE ConToPrim
  // -- so w0 here is the PREVIOUS stage's inversion and predates every one of those
  // edits to u0.  At the wall, which is where the imposed luminosity enters, that is a
  // systematic jump across the boundary face rather than an O(dt) smooth-cell error.
  // (2) w0 IS NOT RESTART STATE.  restart.cpp writes u0 only (ghosts included), and
  // Driver::InitBoundaryValuesAndPrimitives calls ApplyPhysicalBCs BEFORE the first
  // ConToPrim, so at the first boundary call of a restarted run w0 is still exactly
  // zero: the bc_mode-3 walk then started from rho = e = 0, the guard (correctly)
  // rejected it, and the fallback rescaled that same zero, so both x1 walls came back as
  // bare floors -- overwriting the ghosts the restart file had restored bitwise with a
  // state the running boundary never produces.  Reading u0 -- which IS restored bitwise,
  // and which the running boundary reads at exactly the same point of the update --
  // makes the fill idempotent, so a restart is a bitwise continuation.
  // This is red_giant.cpp's fix (state_i in RedGiantBC) for the same defect.
  auto state_i = [=] (const int m, const int k, const int j, const int km, const int jm,
                      const int im, Real &d_i, Real &e_i) {
    d_i = u0(m,IDN,km,jm,im);
    const Real di = (d_i > 0.0) ? (1.0/d_i) : 0.0;
    e_i = u0(m,IEN,km,jm,im)
          - 0.5*(SQR(u0(m,IM1,km,jm,im)) + SQR(u0(m,IM2,km,jm,im))
                 + SQR(u0(m,IM3,km,jm,im)))*di;
    if (etotgrav) e_i -= d_i*phicc(m,km,jm,im);
  };
  auto fill = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                            const int km, const int jm, const int im, const int bcmode) {
    // (k,j,i) the ghost cell, (km,jm,im) the active cell it mirrors -- or, under the
    // OPEN top (bcmode 4), the last active cell every ghost layer is continued from
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real zg = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    const Real zm = CellCenterX(im-is, indcs.nx1, x1min, x1max);
    Real d, e;
    if (bcmode == 4) {
      // THE OPEN (OUTFLOW) TOP.  Zero gradient in the velocity with the inflow clamped
      // away (below), and an ISOTHERMAL HYDROSTATIC continuation of the last active
      // cell in (rho,e): with p = rho (p/rho)_a and dp/dz = -rho g,
      //   rho_g = rho_a exp(-g (z_g - z_a) rho_a/p_a),  e_g = e(rho_g, T_a).
      // A plain zero-gradient copy of (rho,e) would leave a pressure gradient with no
      // weight under it at the top face and blow a steady wind out of the box; this
      // ghost carries no residual force, so the face only sees what the interior does.
      // The TEMPERATURE is what is continued, not the specific energy: with
      // eos_radiation e at fixed T is not proportional to rho.
      Real da, ea;
      state_i(m, k, j, km, jm, im, da, ea);
      const Real ta = eos.Temperature(da, ea);
      const Real pa = eos.Pressure(da, ea, ta);
      Real xarg = 0.0;
      if ((da > 0.0) && (pa > 0.0)) xarg = -g0*(zg - zm)*da/pa;
      xarg = (xarg < -30.0) ? -30.0 : ((xarg > 30.0) ? 30.0 : xarg);
      Real dgh = da*Kokkos::exp(xarg);
      Real egh = eos.EnergyFromTemperature(dgh, ta);
      // the same guard the bc_mode-3 walk needs: an EOS-clamped lid can hand back a
      // zero density or, with eos_radiation, an infinite specific energy
      if (!(Kokkos::isfinite(dgh) && (dgh > 0.0)
            && Kokkos::isfinite(egh) && (egh > 0.0))) {
        dgh = da;
        egh = ea;
      }
      d = (dgh > dfl) ? dgh : dfl;
      const Real efl4 = eos.EnergyFloorBound(d);
      e = (egh > efl4) ? egh : efl4;
    } else if (bcmode == 5) {
      // THE INFLOW BOTTOM.  The ghost is the RELAXED BASE STATE: the initial column at
      // the ghost's OWN height -- the same interpolation bc_mode 0 lays down and the
      // bc_mode-3 walk falls back to -- with the energy rebuilt from the column's
      // TEMPERATURE through the EOS, so the ghost sits on the EOS whatever the table's
      // radiation taper does to e(rho,T).  The velocity is set below.
      Real sg5 = (zg - zlo)/dzf;
      int ig5 = static_cast<int>(sg5);
      ig5 = (ig5 < 0) ? 0 : ((ig5 > nfine-2) ? nfine-2 : ig5);
      const Real fg5 = sg5 - ig5;
      const Real dg5 = cd_d(ig5)*(1.0 - fg5) + cd_d(ig5+1)*fg5;
      const Real eg5 = ce_d(ig5)*(1.0 - fg5) + ce_d(ig5+1)*fg5;
      const Real tg5 = ct_d(ig5)*(1.0 - fg5) + ct_d(ig5+1)*fg5;
      d = (dg5 > dfl) ? dg5 : dfl;
      Real e5 = eos.EnergyFromTemperature(d, tg5);
      if (!(Kokkos::isfinite(e5) && (e5 > 0.0))) e5 = eg5;
      const Real efl5 = eos.EnergyFloorBound(d);
      e = (e5 > efl5) ? e5 : efl5;
    } else if (bcmode == 1) {
      state_i(m, k, j, km, jm, im, d, e);
    } else if (bcmode == 3) {
      // THE WB-CONSISTENT WALL.  Walk the mirror cell's OWN (rho,e) across the wall with
      // the very closure the well-balanced background stencil integrates -- for the
      // general EOS that is utils/wb_background.hpp's WBAdvance, which is what
      // WBBackgroundStencil() calls for its own half-cell segments.  The wall cell's
      // background is then the exact hydrostatic continuation of the state that is
      // actually there, not of the initial column, so it carries no residual force and
      // the wall-face Riemann problem stays symmetric however far the interior has
      // drifted.  The closure test is one-sided (the cell one further IN), because the
      // other side of the mirror cell is the ghost being built.
      Real dmm, emm;
      state_i(m, k, j, km, jm, im, dmm, emm);
      // THE FALLBACK, and the scale the walk is judged against: bc_mode 2's rescaled
      // mirror, the initial column's own ratio across this pair.  Always formed, because
      // a walk out of a sick wall cell must not be able to take the ghost with it.
      Real sg = (zg - zlo)/dzf;
      int ig = static_cast<int>(sg);
      ig = (ig < 0) ? 0 : ((ig > nfine-2) ? nfine-2 : ig);
      const Real fg = sg - ig;
      const Real dcg = cd_d(ig)*(1.0 - fg) + cd_d(ig+1)*fg;
      const Real ecg = ce_d(ig)*(1.0 - fg) + ce_d(ig+1)*fg;
      Real sm = (zm - zlo)/dzf;
      int imc = static_cast<int>(sm);
      imc = (imc < 0) ? 0 : ((imc > nfine-2) ? nfine-2 : imc);
      const Real fm = sm - imc;
      const Real dcm = cd_d(imc)*(1.0 - fm) + cd_d(imc+1)*fm;
      const Real ecm = ce_d(imc)*(1.0 - fm) + ce_d(imc+1)*fm;
      const Real rd_col = (dcm > 0.0) ? (dcg/dcm) : 1.0;
      const Real re_col = (ecm > 0.0) ? (ecg/ecm) : 1.0;
      const int in = (im > i) ? (im + 1) : (im - 1);
      Real dnn, enn;
      state_i(m, k, j, km, jm, in, dnn, enn);
      const Real tmm = eos.Temperature(dmm, emm);
      const int wopt = WBOptionNumber(eos, wbo, dmm, emm, dnn, enn, dmm, emm, tmm);
      Real dlntdphi = 0.0;
      if (wopt == 3) {
        const Real zn = CellCenterX(in-is, indcs.nx1, x1min, x1max);
        const Real dphn = g0*(zn - zm);
        if (dphn != 0.0) {
          dlntdphi = (eos.Temperature(dnn, enn, tmm) - tmm)/dphn;
        }
      }
      Real dw = dmm, ew = emm, tw = tmm;
      WBAdvance(eos, wopt, dmm, emm, g0*(zg - zm), dw, ew, tw, tmm, dlntdphi, tmm, tmm);
      // GUARD THE WALK.  A hydrostatic continuation is only meaningful out of a cell that
      // is itself physical.  At the top of a radiation-dominated atmosphere the wall cell
      // can reach a state the EOS table can only clamp (T ~ 1e13 K at rho ~ 2e-9), and
      // then the polytropic branch's stencil gradient a = dT/dPhi is enormous: the
      // segment's exp(k*dphi) underflows the density to ~0, and with eos_radiation the
      // specific energy carries a_rad T^4/rho, so e comes back +inf.  The outermost ghost
      // sees the largest |dphi| and goes first -- exactly the 8-ghost, 0-active death of
      // bench/hestar_fecz/smoke_rt_3msun_bc3_tau100 at cycle 3500.  A test on positivity
      // alone does NOT catch it (+inf > 0), and the file's own comment warns the walk can
      // also come back "FINITE but absurd".  So demand finite, positive, and a ratio to
      // the mirror cell within wall_walk_maxfac of what the initial column does over the
      // same gap; anything else falls back to the rescaled mirror, which is bounded by
      // construction.  Then floor both, so even the fallback cannot hand the Riemann
      // solver a state below what ConsToPrim would accept.
      const Real rdw = dw/dmm, rew = ew/emm;
      const bool walk_ok = Kokkos::isfinite(dw) && Kokkos::isfinite(ew)
                           && (dw > 0.0) && (ew > 0.0)
                           && Kokkos::isfinite(rdw) && Kokkos::isfinite(rew)
                           && (rdw < rd_col*wfac) && (rdw*wfac > rd_col)
                           && (rew < re_col*wfac) && (rew*wfac > re_col);
      d = walk_ok ? dw : dmm*rd_col;
      e = walk_ok ? ew : emm*re_col;
      d = (d > dfl) ? d : dfl;
      const Real efl = eos.EnergyFloorBound(d);
      e = (e > efl) ? e : efl;
    } else {
      Real sg = (zg - zlo)/dzf;
      int ig = static_cast<int>(sg);
      ig = (ig < 0) ? 0 : ((ig > nfine-2) ? nfine-2 : ig);
      const Real fg = sg - ig;
      const Real dg = cd_d(ig)*(1.0 - fg) + cd_d(ig+1)*fg;
      const Real eg = ce_d(ig)*(1.0 - fg) + ce_d(ig+1)*fg;
      if (bcmode == 0) {
        d = dg;
        e = eg;
      } else {
        Real sm = (zm - zlo)/dzf;
        int im2 = static_cast<int>(sm);
        im2 = (im2 < 0) ? 0 : ((im2 > nfine-2) ? nfine-2 : im2);
        const Real fm = sm - im2;
        const Real dm = cd_d(im2)*(1.0 - fm) + cd_d(im2+1)*fm;
        const Real em = ce_d(im2)*(1.0 - fm) + ce_d(im2+1)*fm;
        Real dmm, emm;
        state_i(m, k, j, km, jm, im, dmm, emm);
        d = dmm*(dg/dm);
        e = emm*(eg/em);
      }
    }
    const Real dm_i = u0(m,IDN,km,jm,im);
    const Real idm = (dm_i > 0.0) ? (1.0/dm_i) : 0.0;
    const Real v1r = u0(m,IM1,km,jm,im)*idm;
    // the outflow top copies v1 and refuses inflow; every other mode -- the inflow
    // bottom included, whose mass enters as a source and not through the wall face --
    // mirrors it, which is what keeps the wall-face Riemann problem symmetric
    const Real v1 = (bcmode == 4) ? ((v1r > 0.0) ? v1r : 0.0) : (-v1r);
    const Real v2 = u0(m,IM2,km,jm,im)*idm;
    const Real v3 = u0(m,IM3,km,jm,im)*idm;
    w0(m,IDN,k,j,i) = d;
    w0(m,IEN,k,j,i) = e;
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = d*v2;
    u0(m,IM3,k,j,i) = d*v3;
    Real et = e + 0.5*d*(v1*v1 + v2*v2 + v3*v3);
    if (etotgrav) et += d*g0*(zg - zmin);
    u0(m,IEN,k,j,i) = et;
  };
  par_for("boxconv_bc_x1", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int n) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      fill(m, k, j, is-1-n, k, j, is+n, bcm_bot);
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      // under the open top EVERY ghost layer continues the SAME cell (ie); the other
      // modes mirror layer by layer
      fill(m, k, j, ie+1+n, k, j, (bcm_top == 4) ? ie : (ie-n), bcm_top);
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvFinal
//! \brief releases the file-scope device views before Kokkos::finalize()

void BoxConvFinal(ParameterInput *pin, Mesh *pm) {
  cd_ = DvceArray1D<Real>();
  ce_ = DvceArray1D<Real>();
  cp_ = DvceArray1D<Real>();
  ct_ = DvceArray1D<Real>();
  // the surface-dump buffers are namespace-scope Views, so their destructors run at
  // static-destruction time, which is AFTER Kokkos::finalize(): release them here
  surf_h_ = HostArray2D<Real>();
  surf_d_ = DvceArray2D<Real>();
  surf_alloc_ = false;
  prof_h_ = HostArray2D<Real>();
  prof_d_ = DvceArray2D<Real>();
  prof_alloc_ = false;
  // and the bottom sponge's plane-mean buffers
  vdb_h_ = HostArray1D<Real>();
  vdb_d_ = DvceArray1D<Real>();
  // and the rt_budget_verbose accumulator, for the same reason
  two_stream_rt::rt_bud_ptr = nullptr;
  rtbud_ = DvceArray1D<Real>();
  rtbud_n_ = 0;
  // ...and the ImEx per-stage source store (problem/rt_imex), likewise
  rtimex_src_ = DvceArray6D<Real>();
  return;
}
