#ifndef PGEN_DEEP_HOT_JUPITER_RT_HPP_
#define PGEN_DEEP_HOT_JUPITER_RT_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter_rt.hpp
//! \brief PRIVATE header of the deep_hot_jupiter_rt problem generator
//! (tusplit-1102): the includes, declarations, file-scope state and the inline
//! device/IC helpers shared by deep_hot_jupiter_rt.cpp (UserProblem) and its
//! deep_hot_jupiter_rt_{hse,src,diag}.cpp pieces.  Moved verbatim out of the old
//! single file; what changed: the file-scope variables are `inline` (still one copy
//! in the program) and the two anonymous namespaces are named dhj_rt_local.  See
//! deep_hot_jupiter_rt.cpp for the description of the problem.

#include <sys/stat.h>  // mkdir, for the cyclediag/ subdirectory

// C++ headers
#include <algorithm>  // sort, for problem/ic_profile
#include <cmath>
#include <cstdint>
#include <cstdio>   // snprintf, for the cycle-diagnostic file name
#include <iostream> // cout
#include <fstream>  // ifstream, for the correlated-k table
#include <sstream>  // ostringstream, for the photosphere dump
#include <iomanip>  // setw/setfill, for the cycle-diagnostic file name
#include <string>
#include <type_traits>  // is_same, roche_wb_tan probe
#include <utility>  // pair, for problem/ic_profile
#include <vector>

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/grid_stretch.hpp"  // StretchR / StretchRPoly
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "mhd/mhd.hpp"
#include "srcterms/srcterms.hpp"
#include "utils/random.hpp"
#include "pgen.hpp"
#include "roche_common.hpp"   // roche-dhj-1031: problem/roche
#include "wbtan_common.hpp"   // roche-wbtan-1101: problem/roche_wb_tan
#include "utils/rosseland.hpp"
#include "utils/correlated_k.hpp"
#include "utils/atm_column.hpp"
#include "utils/eint_from_cons.hpp"
#include "utils/two_stream_rt.hpp"
#include "diffusion/conduction.hpp"
#include "units/units.hpp"
#include "pgen_eos_utils.hpp"
#include "diffusion/resistivity.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "outputs/outputs.hpp"
#include "utils/recipe.hpp"   // proddef-1101: fresh-start defaults

#include <Kokkos_Random.hpp>

// The two-stream RT solver, its configuration and scratch arrays, and the small
// geometry/gravity helpers its kernels call, live in utils/two_stream_rt.hpp and
// utils/atm_column.hpp; they were lifted from here verbatim.  These declarations
// keep every use below spelled exactly as it was.
using atm_column::CSCellAngles;
using atm_column::EffGravAt;
using atm_column::GravAccAt;
using atm_column::GravPotAt;
using atm_column::TGuess;
using atm_column::TideAccP;
using atm_column::TideAccR;
using atm_column::TideAccT;
using two_stream_rt::LimitRTSource;
using two_stream_rt::RTSourceLimiterWarn;
using two_stream_rt::get_Tint;
using two_stream_rt::get_albedo;
using two_stream_rt::get_kapr;
using two_stream_rt::get_picket_fence_coeff;
using two_stream_rt::par_reduce_clip3;
using two_stream_rt::par_reduce_clip4;
using two_stream_rt::picket_fence_two_stream_RT;
using two_stream_rt::rt_B_ptr;
using two_stream_rt::rt_Bb_ptr;
using two_stream_rt::rt_Fb_ptr;
using two_stream_rt::rt_Qb_ptr;
using two_stream_rt::rt_Qv_ptr;
using two_stream_rt::rt_T_ptr;
using two_stream_rt::rt_cf_ptr;
using two_stream_rt::ck_beam_sph;
using two_stream_rt::ck_spherical;
using two_stream_rt::rt_cell_report;
using two_stream_rt::rt_report_every;
using two_stream_rt::rt_ck;
using two_stream_rt::rt_ck_pcut;
using two_stream_rt::rt_de_max;
using two_stream_rt::rt_diag;
using two_stream_rt::rt_diag_ptr;
using two_stream_rt::rt_dump_done;
using two_stream_rt::rt_dump_file;
using two_stream_rt::rt_dump_j;
using two_stream_rt::rt_dump_k;
using two_stream_rt::rt_dump_m;
using two_stream_rt::rt_icut_ptr;
using two_stream_rt::rt_int_at_cut;
using two_stream_rt::rt_kc_ptr;
using two_stream_rt::rt_nchain;
using two_stream_rt::rt_pb_ptr;
using two_stream_rt::rt_semi_implicit;
using two_stream_rt::rt_split;
using two_stream_rt::rt_srclim_warned;
using two_stream_rt::rt_star_teff;
using two_stream_rt::rt_tau_ptr;
using two_stream_rt::rt_use_cons;
using two_stream_rt::rt_relax_sub;
using two_stream_rt::rt_relax_xcrit;
using two_stream_rt::rt_relax_submax;
using two_stream_rt::rt_xP_ptr;
using two_stream_rt::rt_xT_ptr;

// EOS-aware conversions shared with the other stratified problem generators
using pgen_eos::EintFromP;
using pgen_eos::PresFromEint;
using pgen_eos::TempKelvin;
using pgen_eos::PresTempFromEint;




//----------------------------------------------------------------------------------------
using pgen_eos::DensFromPT;
using pgen_eos::GradAd;

void HydrostaticEquilibrium(Mesh *pm);
void SourceFunc(Mesh *pm, Real bdt);

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt);
void double_gray_two_stream_RT(Mesh *pm, Real bdt);

KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3],
                                 const Real albedo_set = -1.0);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real &taulim, const Real &A, const Real &B, const Real (&C)[3], const Real (&D)[3], const Real (&E)[3], const Real (&gamv)[3], const Real &tau, Real &T);
template <typename View1D>
void get_picket_fence_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr);
template <typename View1D>
void adjust_ad_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const int &N, View1D Tarr, View1D lgparr);
template <typename View1D>
void read_ic_profile(const std::string &fname, const int &N, View1D Tarr,
                     View1D lgparr);

void DhjPhotosphereDump(ParameterInput *pin, Mesh *pm);
void DhjOlrDump(ParameterInput *pin, Mesh *pm);
void DhjFinalDumps(ParameterInput *pin, Mesh *pm);
// problem/flux_hst: the radiative / energy / mass flux history columns
void DhjFluxHistory(HistoryData *pdata, Mesh *pm);
// problem/flux_hst_wall: four inner-wall in/out split columns (set in UserProblem)
inline bool dhj_flux_hst_wall = false;
// problem/flux_hst_floor: two floor-bookkeeping columns (set in UserProblem)
inline bool dhj_flux_hst_floor = false;
// relaxdrag-1101 problem/relax_drag (default true), problem/relax_drag_rot (default 5):
// the start-up Rayleigh drag of SourceFunc (tau = P/10 for t < 0.4 R rot, log ramp to
// tau = P at R rot, off after; R = relax_drag_rot).  Fresh starts default relax_drag to
// false (recipe::NewDefaults); a restart without the key keeps the drag, bitwise.
inline bool dhj_relax_drag = true;
inline Real dhj_relax_drag_rot = 5.0;
inline Real dhj_floor_t0 = 0.0;   // start of the interval the floor columns average over
// problem/sponge_top_mode (read without recording): "radial" (DEFAULT since 10-01, user
// GO after am_coriolis_1001/ext arm D) makes the top velocity sponge damp the RADIAL
// momentum only (IM1, rhat.m on both grids), so it exerts no axial torque and no drag
// that reaches the jet below it (Shaw & Shepherd 2007): with it the 1e-5..1e-3 bar jet
// matches a domain whose top and sponge are 2 decades higher, at the old dt and cost.
// "all" damps all three components (the behaviour before 10-01, bitwise).
inline bool dhj_sponge_top_radial = true;
// problem/sponge_top_plo_bar, sponge_top_phi_bar, sponge_top_tau (read without recording;
// defaults 1e-7 bar, 1e-6 bar, 1e3 s): the top sponge ramps from 0 at phi to 1 at plo
// (linear in log p) with damping rate fdrag/tau
inline Real dhj_sp_top_plo_bar = 1.0e-7;
inline Real dhj_sp_top_phi_bar = 1.0e-6;
inline Real dhj_sp_top_tau = 1.0e3;
// problem/flux_hst_rkavg (default true; false = the last RK stage, the pre-0929
// columns): the fluid face fluxes of the flux columns (Etot_top, Etot_bot, Mdot_top,
// Mdot_bot) are the RK-weighted sums over the stages of the cycle, i.e. the flux the
// update actually applied (rk1 1; rk2 1/2, 1/2; rk3 1/6, 1/6, 2/3).  Accumulated per
// column by DhjFluxStage at the top of SourceFunc (uflx then holds the stage's fluxes).
inline bool dhj_fhst_avg = false;
inline int dhj_fhst_nst = 0;                     // stages per cycle
inline Real dhj_fhst_w[3] = {0.0, 0.0, 0.0};     // their weights
inline int dhj_fhst_stage = 0;                   // stages accumulated in this cycle
inline int64_t dhj_fhst_cyc = -1;                // the cycle they belong to
// (m, 4, k, j): E_top M_top E_bot M_bot
inline DvceArray4D<Real> *dhj_fhst_ptr = nullptr;
void DhjFluxStage(Mesh *pm);
void DhjCycleDiag(Mesh *pm);

// AXIAL ANGULAR MOMENTUM (branch am-coriolis, 10-01).  In the corotating frame the
// conserved quantity is L_abs = L_rel + Omega I, I = sum M R^2 (R = distance from the
// rotation axis), and the Coriolis torque on L_rel, -2 Omega rho R v_R, is exactly
// -Omega dI/dt by continuity.  The rotating-frame source applies it with the CELL-CENTRE
// rho v, while I changes by the RIEMANN mass fluxes; the two differ (numerical mass
// diffusion moves I without moving the cell-centre momentum), and the difference is an
// angular-momentum source (w121prod_0929/angmom_1001: ~2.5 % of the gross torque, a
// monotonic L_abs loss).  Both corrections below are applied by AmFix, a separate kernel
// after the rotating-frame source, as a force along the zonal direction (-y, x, 0):
// no radial or meridional part, and the total energy is not touched (the Coriolis force
// never entered it), so their small work lands in the internal energy.
//   problem/coriolis_am: take the ZONAL Coriolis torque from the face mass fluxes,
//     -Omega Q_c with
//         Q_c = sum_{faces f of c} Phi_f^out (W_f - W_c),   W = R^2,
//     Phi^out the outward mass flux (flux x area) the update applies and W_f the face
//     centre's R^2.  Q_c -> int rho V.grad(R^2) dV = 2 int rho R v_R dV (consistent), and
//     sum_c Q_c = dI/dt + (outflow of W through the domain boundary) EXACTLY, so the
//     Coriolis torque cancels Omega dI/dt to round-off.
//   problem/am_transport: the transport is not angular-momentum conservative either.
//     Cubed sphere: the covariant tangential momenta move through faces whose basis
//     differs from the cell's, and the angular metric source that makes up for it is
//     evaluated with CELL-CENTRE p, rho v v (and B B) (w121 1x: +1.9e33 g cm^2 s^-2).
//     Spherical polar: the phi-momentum update is flux-based but its effective torque arm
//     is not the cell's R.  The zonal update is corrected by
//         T_cons - T_code,
//     T_code the torque the update applied (AmCell re-evaluates the flux divergences and
//     the angular source of SrcTermsGnomonicEquiangleImpl / SrcTermsSphericalPolar*:
//     keep them in step), T_cons the torque of the same face fluxes taken as
//     ANGULAR-MOMENTUM fluxes with exact face first moments: cs r_f (rhat x F) A on the
//     radial faces and (int_f x dA) x F on the planar xi/eta faces; sp F_phi A Rbar_f,
//     Rbar_f the face's area-mean R.  sum_c T_cons is the boundary flux alone, so L_rel
//     changes only by Coriolis, boundaries and sponges, and a uniform pressure (or a
//     uniform Maxwell stress on the planar faces) exerts no T_cons on a cell.
//   DEFAULTS (user 10-01): both ON whenever the rotating frame is (user_srcs, omega != 0,
//     cubed sphere or spherical polar), OFF otherwise, so non-rotating runs are bitwise
//     unchanged.  A rotating run with a source variant am_transport does not
//     replicate (cs/sp_wellbalanced_src, wellbalance_static, sp_cart_polar_momentum,
//     sp_face_avg) stops with a FATAL naming it, unless problem/am_transport = false is
//     set explicitly; a rotating run with either key false prints a WARNING that L_abs
//     drifts.  Read without recording; the startup log prints the effective values.
//   problem/am_hst (DIAGNOSTIC, default false; read without recording it): eight more
//     history columns after the flux columns, see DhjFluxHistory and AmHstStage.
//   problem/am_hst_split (DIAGNOSTIC, default false, needs am_hst and what am_transport
//     needs): the last two columns become five: the torques of the x1 (with its r-factor
//     source), x2 and x3 flux divergences and of the angular metric source (T_code is
//     their sum), and T_cons; then (condwall-1028) Twall, the inner-wall part of T_cons:
//     the angular-momentum flux IN through the inner wall (the wall torque, Maxwell
//     stress included).  The rest, Ttop = Tcons - Twall, comes in through the top, and
//     the L_z ledger reads dLrel/dt = Tcor + Twall + Ttop.  (Ttop has no column of its
//     own: 9 + 3 floor + 12 = 24 = NHISTORY_VARIABLES; with problem/roche the Ttide
//     column makes it 25, so am_hst_split + flux_hst_floor + roche is a start-up FATAL.)
namespace amh {
inline bool on = false;            // problem/am_hst
inline bool fix = false;           // problem/coriolis_am
inline bool split = false;         // problem/am_hst_split
inline bool trans = false;         // problem/am_transport
inline Real pdeep = 1.0e8;  // problem/am_hst_pdeep [barye]: the Tdif_deep split (100 bar)
inline int nst = 0;                // stages per cycle and their weights (as dhj_fhst_w)
inline Real w[3] = {0.0, 0.0, 0.0};
inline int stage = 0;
inline int64_t cyc = -1;
// Tcc Tcont Tbnd Tdif_deep Tdif_rad | Tx1 Tx2 Tx3 Tsrc Tcons Twall | Tdrag (relaxdrag-1101)
constexpr int NACC = 12;
// since last
inline Real acc[NACC] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
inline Real t0 = 0.0;
inline Real tdrag_rate = 0.0;      // relaxdrag-1101: Tdrag of the last AmHstColumns call
}  // namespace amh
void AmHstStage(Mesh *pm);
Real DhjRelaxDragRate(const Real time, const Real omega);
Real DhjLrelU0(Mesh *pm);
void AmFix(Mesh *pm, const Real bdt);
void AmHstColumns(HistoryData *pdata, Mesh *pm, const int n0);

//----------------------------------------------------------------------------------------
//! \struct AmV
//! \brief everything AmCell reads, captured BY VALUE by the kernels (AmFix, AmHstStage)

struct AmV {
  DvceArray5D<Real> w0, bcc, wder, f1, f2, f3;
  DvceArray4D<Real> a1, a2, a3, vol, xovr, yovr;
  DvceArray3D<Real> ccell, scell;
  DvceArray2D<Real> x1v, x1f, x2v, x2f;
  DualArray1D<RegionSize> size;
  DualArray1D<int> panel;
  EOS_Data eos;
  Real omega;
  int is, ie, js, ks, nx2, nx3;
  bool cs, mhd, gen, trq;   // trq: also the transport torques
};

inline AmV AmMakeV(Mesh *pm, const bool trq) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  AmV v;
  v.mhd = (pmbp->pmhd != nullptr);
  v.w0 = v.mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  v.eos = v.mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  v.gen = v.eos.IsGeneral();
  v.wder = v.mhd ? pmbp->pmhd->wder : pmbp->phydro->wder;
  if (v.mhd) v.bcc = pmbp->pmhd->bcc0;
  DvceFaceFld5D<Real> &uf = v.mhd ? pmbp->pmhd->uflx : pmbp->phydro->uflx;
  v.f1 = uf.x1f;  v.f2 = uf.x2f;  v.f3 = uf.x3f;
  v.a1 = pmbp->pcoord->area.x1f;
  v.a2 = pmbp->pcoord->area.x2f;
  v.a3 = pmbp->pcoord->area.x3f;
  v.vol = pmbp->pcoord->volume;
  v.xovr = pmbp->pcoord->x_ov_rD;
  v.yovr = pmbp->pcoord->y_ov_rC;
  v.ccell = pmbp->pcoord->cos_cell;
  v.scell = pmbp->pcoord->sin_cell;
  v.x1v = pmbp->pcoord->x1v;
  v.x1f = pmbp->pcoord->xx1f;
  v.x2v = pmbp->pcoord->x2v;
  v.x2f = pmbp->pcoord->xx2f;
  v.size = pmbp->pmb->mb_size;
  v.panel = pmbp->pmb->mb_panel;
  v.omega = pm->pgen->hot_jupiter_param.omega;
  auto &ic = pm->mb_indcs;
  v.is = ic.is;  v.ie = ic.ie;  v.js = ic.js;  v.ks = ic.ks;
  v.nx2 = ic.nx2;  v.nx3 = ic.nx3;
  v.cs = pm->use_cubed_sphere;
  v.trq = trq;
  return v;
}

//----------------------------------------------------------------------------------------
//! \struct AmOut
//! \brief AmCell's results for one cell (torques in g cm^2 s^-2, i.e. already x vol).
//! A zonal torque dT is applied as  m2 += dT*c2*bdt,  m3 += dT*c3*bdt.

struct AmOut {
  Real wc;        // R^2 at the cell centre
  Real c2, c3;    // zonal-torque -> momentum-density coefficients (0 when R = 0)
  Real tcc;       // the cell-centre Coriolis torque the rotating-frame source applied
  Real q, qr;     // Q_c and its radial-face part
  Real tccr;      // the radial-velocity part of tcc
  Real tbnd;      // -Omega Phi^out W_f over the domain's x1 faces (is, ie+1)
  Real tc[4];     // T_code: x1 (+ r-factor source), x2, x3, angular metric source
  Real tcons;     // T_cons
  Real twall;     // the inner-wall (is face) part of T_cons: torque IN through the wall
};

//----------------------------------------------------------------------------------------
//! \brief AmCell: the angular-momentum bookkeeping of cell (m,k,j,i), see namespace amh.
//! Cubed sphere: the stored tangential flux components are covariant on the FACE's
//! basis (the cell's for x1 faces) and are raised with that basis' e1.e2; x1 faces have
//! W_f = r_f^2 s2c, xi/eta faces the face centre's r^2 (q_x^2 + q_y^2).  Spherical polar:
//! orthonormal components, R = r sin(theta); W_f on x1 faces r_f^2 sin^2(theta_c), on
//! theta faces r_c^2 sin^2(theta_f), and the phi faces do not change R.

KOKKOS_INLINE_FUNCTION
void AmCell(const AmV &v, const int m, const int k, const int j, const int i, AmOut &o) {
  const Real om = v.omega;
  const Real dv = v.vol(m,k,j,i);
  const Real ia = 1.0/dv;
  const Real r = v.x1v(m,i);
  const Real rl = v.x1f(m,i), rr = v.x1f(m,i+1);
  const Real fac = (rr - rl)/(rr + rl);
  const Real rho = v.w0(m,IDN,k,j,i);
  const Real v1 = v.w0(m,IVX,k,j,i), v2 = v.w0(m,IVY,k,j,i), v3 = v.w0(m,IVZ,k,j,i);
  o.tbnd = 0.0;
  for (int n=0; n<4; ++n) o.tc[n] = 0.0;
  o.tcons = 0.0;
  o.twall = 0.0;
  if (v.cs) {
    const Real x2l = v.size.d_view(m).x2min, x2u = v.size.d_view(m).x2max;
    const Real x3l = v.size.d_view(m).x3min, x3u = v.size.d_view(m).x3max;
    const Real xc = 0.25*M_PI*CellCenterX(j-v.js, v.nx2, x2l, x2u);
    const Real ec = 0.25*M_PI*CellCenterX(k-v.ks, v.nx3, x3l, x3u);
    const Real xl = 0.25*M_PI*LeftEdgeX(j-v.js, v.nx2, x2l, x2u);
    const Real xr = 0.25*M_PI*LeftEdgeX(j+1-v.js, v.nx2, x2l, x2u);
    const Real el = 0.25*M_PI*LeftEdgeX(k-v.ks, v.nx3, x3l, x3u);
    const Real er = 0.25*M_PI*LeftEdgeX(k+1-v.ks, v.nx3, x3l, x3u);
    const int pnl = v.panel.d_view(m);
    Real qc[3], e1[3], e2[3];
    cubed_sphere::PanelToCart(pnl, xc, ec, qc);
    cubed_sphere::PanelTangents(pnl, xc, ec, e1, e2);
    const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
    const Real rh0 = qc[0]*qn, rh1 = qc[1]*qn;
    const Real xx = r*rh0, yy = r*rh1;
    const Real s2c = rh0*rh0 + rh1*rh1;
    const Real wc = xx*xx + yy*yy;
    o.wc = wc;
    o.c2 = (wc > 0.0) ? (-yy*e1[0] + xx*e1[1])/(wc*dv) : 0.0;
    o.c3 = (wc > 0.0) ? (-yy*e2[0] + xx*e2[1])/(wc*dv) : 0.0;
    const Real vcx = v1*rh0 + v2*e1[0] + v3*e2[0];
    const Real vcy = v1*rh1 + v2*e1[1] + v3*e2[1];
    o.tcc = -2.0*om*rho*(xx*vcx + yy*vcy)*dv;
    o.tccr = -2.0*om*rho*v1*r*s2c*dv;
    // Q_c
    const Real r2 = r*r;
    Real q[3];
    o.qr = v.f1(m,IDN,k,j,i+1)*v.a1(m,k,j,i+1)*(rr*rr*s2c - wc)
         - v.f1(m,IDN,k,j,i)*v.a1(m,k,j,i)*(rl*rl*s2c - wc);
    cubed_sphere::PanelToCart(pnl, xr, ec, q);
    Real qh = v.f2(m,IDN,k,j+1,i)*v.a2(m,k,j+1,i)*(r2*(q[0]*q[0] + q[1]*q[1]) - wc);
    cubed_sphere::PanelToCart(pnl, xl, ec, q);
    qh -= v.f2(m,IDN,k,j,i)*v.a2(m,k,j,i)*(r2*(q[0]*q[0] + q[1]*q[1]) - wc);
    cubed_sphere::PanelToCart(pnl, xc, er, q);
    qh += v.f3(m,IDN,k+1,j,i)*v.a3(m,k+1,j,i)*(r2*(q[0]*q[0] + q[1]*q[1]) - wc);
    cubed_sphere::PanelToCart(pnl, xc, el, q);
    qh -= v.f3(m,IDN,k,j,i)*v.a3(m,k,j,i)*(r2*(q[0]*q[0] + q[1]*q[1]) - wc);
    o.q = o.qr + qh;
    if (i == v.ie) o.tbnd -= om*v.f1(m,IDN,k,j,i+1)*v.a1(m,k,j,i+1)*rr*rr*s2c;
    if (i == v.is) o.tbnd += om*v.f1(m,IDN,k,j,i)*v.a1(m,k,j,i)*rl*rl*s2c;
    if (!v.trq) return;
    // T_code: the flux divergences and SrcTermsGnomonicEquiangleImpl's angular source
    const Real c = v.ccell(m,k,j);
    const Real s2 = 1.0 - c*c;
    Real dx1[2], dx2[2], dx3[2];
    for (int n=0; n<2; ++n) {
      const int nm = IM2 + n;
      dx1[n] = -(v.f1(m,nm,k,j,i+1)*v.a1(m,k,j,i+1) - v.f1(m,nm,k,j,i)*v.a1(m,k,j,i))*ia
               - fac*(v.f1(m,nm,k,j,i)*v.a1(m,k,j,i)
                      + v.f1(m,nm,k,j,i+1)*v.a1(m,k,j,i+1))*ia;
      dx2[n] = -(v.f2(m,nm,k,j+1,i)*v.a2(m,k,j+1,i) - v.f2(m,nm,k,j,i)*v.a2(m,k,j,i))*ia;
      dx3[n] = -(v.f3(m,nm,k+1,j,i)*v.a3(m,k+1,j,i) - v.f3(m,nm,k,j,i)*v.a3(m,k,j,i))*ia;
    }
    const Real pr = v.gen ? v.wder(m,IDPR,k,j,i) : v.eos.Pressure(rho, v.w0(m,IEN,k,j,i));
    const Real sine = v.scell(m,k,j);
    const Real sine2 = sine*sine;
    Real ptot = pr, bb2 = 0.0, bb3 = 0.0;
    if (v.mhd) {
      const Real b1 = v.bcc(m,IBX,k,j,i), b2 = v.bcc(m,IBY,k,j,i);
      const Real b3 = v.bcc(m,IBZ,k,j,i);
      bb2 = (sine*b2 - c*b3)*(sine*b2 - c*b3);
      bb3 = b3*b3;
      ptot += 0.5*(b1*b1 + b2*b2 + b3*b3);
    }
    const Real src2 = v.xovr(m,k,j,i)*(ptot + rho*v3*v3*sine2 - bb3);
    const Real src3 = v.yovr(m,k,j,i)*(ptot + rho*v2*v2*sine2 - bb2);
    const Real d2[4] = {dx1[0], dx2[0], dx3[0], src2};
    const Real d3[4] = {dx1[1], dx2[1], dx3[1], src3};
    for (int n=0; n<4; ++n) {
      const Real u2 = (d2[n] - c*d3[n])/s2, u3 = (d3[n] - c*d2[n])/s2;
      o.tc[n] = (xx*(u2*e1[1] + u3*e2[1]) - yy*(u2*e1[0] + u3*e2[0]))*dv;
    }
    // T_cons, radial faces: arm r_f rhat_c, (rhat x F)_z, F on the cell's basis
    Real t = 0.0;
    for (int n=0; n<2; ++n) {
      const int ii = i + n;
      const Real g2 = v.f1(m,IM2,k,j,ii), g3 = v.f1(m,IM3,k,j,ii);
      const Real fa = (g2 - c*g3)/s2, fb = (g3 - c*g2)/s2;
      const Real tr = v.a1(m,k,j,ii)*v.x1f(m,ii)*(rh0*(fa*e1[1] + fb*e2[1])
                                                  - rh1*(fa*e1[0] + fb*e2[0]));
      t += n ? -tr : tr;
      if (n == 0 && i == v.is) o.twall = tr;
    }
    // planar xi / eta faces: (int_f x dA) x F, int_f x dA = (rr^3 - rl^3)/3 tan(psi/2)
    // (q_a + q_b), q_a/q_b the directions of the face's two ends, psi the angle between
    const Real mom = (rr*rr*rr - rl*rl*rl)/3.0;
    for (int dir=0; dir<2; ++dir) {
      for (int n=0; n<2; ++n) {
        Real xf, ef, xa, ea, xb, eb, g2, g3;
        if (dir == 0) {
          xf = n ? xr : xl;  ef = ec;  xa = xf;  ea = el;  xb = xf;  eb = er;
          g2 = v.f2(m,IM2,k,j+n,i);  g3 = v.f2(m,IM3,k,j+n,i);
        } else {
          xf = xc;  ef = n ? er : el;  xa = xl;  ea = ef;  xb = xr;  eb = ef;
          g2 = v.f3(m,IM2,k+n,j,i);  g3 = v.f3(m,IM3,k+n,j,i);
        }
        Real qa[3], qb[3], b1[3], b2[3];
        cubed_sphere::PanelToCart(pnl, xa, ea, qa);
        cubed_sphere::PanelToCart(pnl, xb, eb, qb);
        cubed_sphere::PanelTangents(pnl, xf, ef, b1, b2);
        const Real cf = b1[0]*b2[0] + b1[1]*b2[1] + b1[2]*b2[2];
        const Real sf = 1.0 - cf*cf;
        const Real fa = (g2 - cf*g3)/sf, fb = (g3 - cf*g2)/sf;
        const Real cx = qa[1]*qb[2] - qa[2]*qb[1];
        const Real cy = qa[2]*qb[0] - qa[0]*qb[2];
        const Real cz = qa[0]*qb[1] - qa[1]*qb[0];
        const Real dot = qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2];
        const Real th = sqrt(cx*cx + cy*cy + cz*cz)/(1.0 + dot);   // tan(psi/2)
        const Real mx = mom*th*(qa[0] + qb[0]), my = mom*th*(qa[1] + qb[1]);
        const Real tf = mx*(fa*b1[1] + fb*b2[1]) - my*(fa*b1[0] + fb*b2[0]);
        t += n ? -tf : tf;
      }
    }
    o.tcons = t;
  } else {
    // spherical polar: R = r sin(theta), the phi momentum is IM3 (orthonormal)
    const Real th = v.x2v(m,j);
    const Real thl = v.x2f(m,j), thr = v.x2f(m,j+1);
    const Real sn = sin(th), cn = cos(th);
    const Real rc = r*sn;
    const Real wc = rc*rc;
    o.wc = wc;
    o.c2 = 0.0;
    o.c3 = (rc > 0.0) ? ia/rc : 0.0;
    o.tcc = -2.0*om*rho*(v1*sn + v2*cn)*rc*dv;
    o.tccr = -2.0*om*rho*v1*sn*rc*dv;
    const Real sl = sin(thl), sr = sin(thr);
    o.qr = v.f1(m,IDN,k,j,i+1)*v.a1(m,k,j,i+1)*(rr*rr*sn*sn - wc)
         - v.f1(m,IDN,k,j,i)*v.a1(m,k,j,i)*(rl*rl*sn*sn - wc);
    o.q = o.qr + v.f2(m,IDN,k,j+1,i)*v.a2(m,k,j+1,i)*(r*r*sr*sr - wc)
               - v.f2(m,IDN,k,j,i)*v.a2(m,k,j,i)*(r*r*sl*sl - wc);
    if (i == v.ie) o.tbnd -= om*v.f1(m,IDN,k,j,i+1)*v.a1(m,k,j,i+1)*rr*rr*sn*sn;
    if (i == v.is) o.tbnd += om*v.f1(m,IDN,k,j,i)*v.a1(m,k,j,i)*rl*rl*sn*sn;
    if (!v.trq) return;
    // T_code: R_c vol (flux divergence + SrcTermsSphericalPolar*'s src3)
    const Real fl1 = v.f1(m,IM3,k,j,i)*v.a1(m,k,j,i);
    const Real fr1 = v.f1(m,IM3,k,j,i+1)*v.a1(m,k,j,i+1);
    const Real fl2 = v.f2(m,IM3,k,j,i)*v.a2(m,k,j,i);
    const Real fr2 = v.f2(m,IM3,k,j+1,i)*v.a2(m,k,j+1,i);
    const Real fl3 = v.f3(m,IM3,k,j,i)*v.a3(m,k,j,i);
    const Real fr3 = v.f3(m,IM3,k+1,j,i)*v.a3(m,k+1,j,i);
    o.tc[0] = rc*(-(fr1 - fl1) - fac*(fl1 + fr1));
    o.tc[1] = rc*(-(fr2 - fl2));
    o.tc[2] = rc*(-(fr3 - fl3));
    o.tc[3] = rc*(-v.yovr(m,k,j,i)*(fl2 + fr2));
    // T_cons: F_phi A Rbar_f, Rbar_f the face's area-mean R (exact)
    const Real cl = cos(thl), cr = cos(thr);
    const Real sin2int = 0.5*(thr - thl) - 0.25*(sin(2.0*thr) - sin(2.0*thl));
    const Real dcos = cl - cr;
    const Real r3 = (rr*rr*rr - rl*rl*rl)/3.0, r2h = 0.5*(rr*rr - rl*rl);
    const Real rb1 = (dcos > 0.0) ? sin2int/dcos : 0.0;           // x Rbar/r_f
    const Real rb2 = r3/r2h;                                      // x sin(theta_f)
    const Real rb3 = (thr > thl) ? r3*dcos/(r2h*(thr - thl)) : 0.0;
    o.tcons = -(fr1*rr*rb1 - fl1*rl*rb1) - (fr2*rb2*sr - fl2*rb2*sl) - (fr3 - fl3)*rb3;
    if (i == v.is) o.twall = fl1*rl*rb1;
  }
}
// problem/budget_dt (DIAGNOSTIC, branch dhj-budget-0929; off unless the key exists and
// is > 0): per-shell energy and mass ledger, see EBudCk
namespace ebud {
constexpr int NQ = 20;
inline bool on = false;
inline bool started = false;
inline Real dt_out = 0.0;
inline Real t_next = 0.0;
inline DvceArray2D<Real> *acc = nullptr;    // (NQ, nx1+1)
inline DvceArray2D<Real> *tmp = nullptr;    // (2, nx1+1): shell E and M
inline DvceArray2D<Real> *fhst = nullptr;  // (4, nx1): floors pending for the hst columns
// (nmb, 2, n3, n2, n1): u0(IEN) and e_int before ck
inline DvceArray5D<Real> *cpre = nullptr;
inline int64_t nfull0 = 0;                  // ck_cad_nfull before the call
}  // namespace ebud
// roche-wbtan-1101 problem/roche_wb_tan (default false): the TANGENTIAL (x2/x3) part of
// the potential force (centrifugal of rot_potential + the Roche tide) as the hydro's OWN
// discrete tangential momentum residual of the equipotential equilibrium W_eq, instead of
// the analytic -rho grad_t Phi_tot at the cell centre:
//     S_t = -(rho/rho_eq) R_t(W_eq),   R_t = -(flux divergence) + geometric source,
// evaluated ONCE per run segment (WbTanProbe, first SourceFunc call, fresh start or
// restart) by running the stage's own operators on W_eq: ghost exchange (panel-seam
// resample) + x1 BCs + ConToPrim, Hydro::Fluxes, the seam flux reconciliation, the
// coordinate source and both angular-momentum corrections.  W_eq is the IC formula (the
// column at z_eff(Phi_tot)), so the IC is balanced in x2/x3 to round-off, seams and sp
// poles included.  The work is booked by the etotgrav energy flux (required), as for
// the analytic force before: E_tot stays exactly conservative, the O(h^2) mismatch
// between rho g_t . v and the flux-form potential work lands in e_int, as before.
inline wbtan::WbTan dhj_wbt;   // wbtan_common.hpp
void WbTanProbe(Mesh *pm, const Real bdt);
void EBudStage(Mesh *pm, const Real bdt);
void EBudCk(Mesh *pm, const Real dt, const bool before);
void EBudDrainFloors(Mesh *pm);
// problem/ck_impl_once: the whole correlated-k radiation, once per hydro step,
// with the full cycle dt, after the last RK stage
void DhjCkRtSplit(Mesh *pm, const Real dt);

// PER-CYCLE SINGLE-MESHBLOCK DIAGNOSTIC (problem/diag_gid, default -1 = off).  Which
// global meshblock to dump, and the inclusive cycle window.  File-scope because the
// per-cycle hook takes only the Mesh.
namespace dhj_rt_local {   // tusplit-1102: was an anonymous namespace
inline int diag_gid = -1;
inline int diag_cycle_min = 0;
inline int diag_cycle_max = 2147483647;
// a char buffer, not a std::string: a global string is a style violation, and the
// name is only ever used to build the dump file name
inline char diag_basename[256] = "athena";
inline bool diag_dir_made = false;
}  // namespace dhj_rt_local
using dhj_rt_local::diag_gid;
using dhj_rt_local::diag_cycle_min;
using dhj_rt_local::diag_cycle_max;
using dhj_rt_local::diag_basename;
using dhj_rt_local::diag_dir_made;

KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td);
KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T);
KOKKOS_INLINE_FUNCTION
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T,
                 const int &Ne, const DvceArray1D<Real> &Text, const DvceArray1D<Real> &lgpext);
template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T,
                      const int &Ne, const View1D &Text, const View1D &lgpext);
template <typename View1D>
int ic_adiabat_extension(const EOS_Data &eos, const Real &gamma, const Real &Rgas,
                         const int &N, const View1D &Tarr, const View1D &lgparr,
                         const Real &dex, const int &Ne, View1D Text, View1D lgpext);
template <typename View1D>
void get_wb_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc,
                    const Real &ap, const bool &grav_pmass, const int &N,
                    const Real &zmax, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_wb_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);
template <typename View1D>
void get_init_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc,
                      const Real &ap, const bool &grav_pmass, const View1D &Tarr,
                      const View1D &lgparr, const int &N, const Real &zmax,
                      View1D zarr, View1D logparr, const Real zlow,
                      const int &Ne, const View1D &Text, const View1D &lgpext);
KOKKOS_INLINE_FUNCTION
void get_init_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p,
                  const int &Ne, const DvceArray1D<Real> &Text, const DvceArray1D<Real> &lgpext);

//----------------------------------------------------------------------------------------
//! \fn GravAccAt / GravPotAt
//! \brief the gravity model. `<problem>/grav` is the SURFACE gravity g0 at r = ap.
//!
//! By default gravity is constant, g(r) = g0, and the potential is g0*(r - ap). That is a
//! thin-shell approximation, and this problem is not thin: the production domain reaches
//! r/ap = 1.52, where a point mass with the same g0 gives 0.43*g0. Holding g constant
//! over that range compresses the upper atmosphere -- at 1e-6 bar the level sits ~5 scale
//! heights too low and H is understated by ~2.7x, which propagates straight into
//! transmission radii and feature amplitudes.
//!
//! `<problem>/grav_point_mass = true` switches to a point mass with the SAME g0 at ap:
//!     g(r)   = g0 (ap/r)^2
//!     phi(r) = g0 ap (1 - ap/r)
//! Both reduce to the constant-g forms as r -> ap, so the potential's zero point and the
//! shallow-atmosphere limit are unchanged and the flag is a no-op at the inner boundary.
//!
//! grav_acc is NEGATIVE throughout this file (it is -g0), hence the sign juggling.
//!
//! Every call site carried a commented-out `- 0.5*SQR(omega*r*sin(theta))` after the
//! potential -- the centrifugal term, never enabled. It is recorded here once rather than
//! repeated nine times; omega*r at the top is 6 cm/s^2 against g of 400-940, so it is a
//! sub-percent correction if it is ever wanted.

//----------------------------------------------------------------------------------------
//! \fn ApplyRStretch
//! \brief map a uniform radial coordinate through whichever radial stretch is active.
//!
//! The Coordinates class already stretches x1v/x1f/dx1, so anything reading those arrays
//! needs nothing. But this file also rebuilds radii directly with CellCenterX/LeftEdgeX
//! when it lays down the initial condition and the gravitational potential, and those
//! come out UNSTRETCHED. Without this the coordinates and the initial condition would sit
//! on different grids -- silently, since both are individually self-consistent.
//!
//! StretchR/StretchRPoly are the same file-scope definitions Coordinates uses, so the two
//! cannot drift apart.

// ApplyRStretch now lives in coordinates/grid_stretch.hpp, shared with cs_test.



//! \fn TotPotAt
//! \brief the gravitational plus centrifugal potential, -Omega^2 (r sin theta)^2 / 2 on
//! top of GravPotAt (Phi increases outward here, and rotation lowers it at the equator).
//! om2 = 0 (problem/rot_potential = false) returns the gravitational potential alone.
KOKKOS_INLINE_FUNCTION
Real TotPotAt(const Real grav_acc, const Real ap, const Real r, const Real z,
              const bool pmass, const Real om2, const Real theta) {
  return GravPotAt(grav_acc, ap, r, z, pmass) - 0.5*om2*SQR(r*sin(theta));
}

//! \fn ZEffFromPot
//! \brief the height at which the gravitational potential alone equals `phi`: the
//! equivalent height of an equipotential surface.  The 1D column p(z), rho(z) integrated
//! in the gravitational potential is a function of Phi_grav; sampling it at
//! z_eff(Phi_tot)
//! makes the initial state (and the ghost column, filled by the same kernel) a barotropic
//! equilibrium of the total potential, balanced in r AND theta.
KOKKOS_INLINE_FUNCTION
Real ZEffFromPot(const Real grav_acc, const Real ap, const Real phi, const bool pmass) {
  const Real g = -grav_acc;
  if (!pmass) return phi/g;
  return ap/(1.0 - phi/(g*ap)) - ap;
}

//! \fn TotPotAt (Roche overload, roche-dhj-1031)
//! \brief TotPotAt plus the host star's tidal potential (roche_common.hpp TidePot:
//! companion + indirect term, 0 at the planet centre), with the substellar direction
//! cosine mu = sin(theta) cos(phi) -- the convention of the stellar beam (rt_pre_geom
//! mu0) and of CSCellAngles (phi = longitude from the substellar point).  gms = 0
//! (problem/roche off) returns the plain TotPotAt value BITWISE (no addition at all).
//! Together with om2 = Omega^2 this is the full Roche potential in the planet-centred
//! corotating frame: the centrifugal part lives HERE (om2) and nowhere else -- the
//! explicit source keeps only Coriolis and the tangential projections of grad Phi.
KOKKOS_INLINE_FUNCTION
Real TotPotAt(const Real grav_acc, const Real ap, const Real r, const Real z,
              const bool pmass, const Real om2, const Real theta, const Real phi,
              const Real gms, const Real ra) {
  const Real p0 = TotPotAt(grav_acc, ap, r, z, pmass, om2, theta);
  if (gms == 0.0) return p0;
  return p0 + roche::TidePot(gms, ra, r, sin(theta)*cos(phi));
}


//----------------------------------------------------------------------------------------
//! \fn TideAccR / TideAccT / TideAccP
//! \brief the host star's tidal acceleration, in the planet-centred corotating frame.
//!
//! `<problem>/stellar_tide = true` adds the term the frame has always been missing. The
//! full Hill acceleration in a frame corotating with the orbit is
//!     a = Omega^2 (3x, 0, -z),
//! with x towards the star and z along the spin axis. The block below already applies the
//! centrifugal term about the PLANET's own axis, Omega^2 (x, y, 0) -- see `oor`. What
//! is left to add is the difference,
//!     a_tide = Omega^2 (2x, -y, -z),
//! i.e. the stellar tidal tensor proper: outward at the substellar AND antistellar points
//! (+2 Omega^2 r), inward at the terminators and the poles (-Omega^2 r). A prolate bulge
//! along the star-planet axis.
//!
//! In the spherical basis, with mu = sin(theta) cos(phi) the substellar direction cosine
//! (the same `mu0` the radiative transfer uses for the stellar beam):
//!     a_r     = Omega^2 r (3 mu^2 - 1)
//!     a_theta = 3 Omega^2 r sin(theta) cos(theta) cos^2(phi)
//!     a_phi   = -3 Omega^2 r sin(theta) sin(phi) cos(phi)
//!
//! WHY IT MATTERS HERE: a_tide/g scales as r^3, so it is set by how tall the domain is.
//! At the old constant-g top (r/ap = 1.52) it is 1.5% of g at the terminator; with
//! grav_point_mass and r/ap = 2.18 it is 4.4% at the top and 2.9% at the 1e-6 bar
//! terminator isobar, and 3 Omega^2 r / g reaches 13% at the top. r_L1 = 4.28 ap, so the
//! domain now reaches half way to L1. A hydrostatic estimate against the run's own
//! measured H(r) puts the 1e-6 bar level 3.2 cells lower at the terminator and 9.3 cells
//! higher at the substellar point -- about 0.15 H on the limb observable.
//!
//! Like the centrifugal term, this is applied as an explicit momentum source and is NOT
//! folded into phicc0/phi0_x1f, so the well-balanced reconstruction does not see it. That
//! is deliberate and matches how `oor` has always been treated: the well-balanced path
//! balances the spherically symmetric radial potential, and the tidal potential is not
//! spherically symmetric. So switching this flag on leaves the initial
//! condition (built hydrostatic WITHOUT the tide) no longer exactly balanced; the
//! atmosphere adjusts over the first few dynamical times, which the initial Rayleigh drag
//! damps.
//!
//! WHERE IT IS AND IS NOT APPLIED. Applied at four places: the three momentum components
//! plus the work term in the corotating-frame block, and the radiative transfer's
//! `tau = kappa p / g` closure over the unresolved column above the domain (three call
//! sites, via EffGravAt). Deliberately NOT applied at three others:
//!   * get_wb_eos_arr / get_init_eos_arr, the 1D hydrostatic column integrators behind
//!     the initial condition and the well-balanced background. They have no angle, and
//!     the background must stay spherically symmetric for the well-balanced path.
//!   * the Maxwell outer-BC ghost clamp, which is MHD-only and uses |g| purely as a
//!     reference weight for a dimensionless ratio; the untided value is the right scale
//!     there, and theta/phi are not in scope.




//----------------------------------------------------------------------------------------
//! \fn EffGravAt
//! \brief inward gravity including the tidal correction, for the radiative transfer's
//! `tau = kappa p / g` closure over the unresolved column above the domain.
//!
//! `grav` here is POSITIVE (the RT sites pass g0, not -g0). The tidal radial acceleration
//! is outward-positive, so it is subtracted. Clamped at 10% of the untided value so a
//! domain taken close to L1 cannot divide by zero; at the production top it is
//! -8.8% substellar and +4.4% at the terminator, nowhere near the clamp.



// The correlated-k OPACITY layer now lives in utils/correlated_k.hpp so that other
// problem generators can use the same bands and tables; it was lifted from here
// verbatim.  These declarations keep every use below spelled exactly as it was.
using correlated_k::CK_DIFFUSIVITY;
using correlated_k::CK_NB;
using correlated_k::CK_NG;
using correlated_k::CK_NPF;
using correlated_k::CK_PF_TMIN;
using correlated_k::CK_PF_TMAX;
using correlated_k::CK_NCIA;
using correlated_k::CK_CIA_NTMAX;
using correlated_k::CK_NRAY;
using correlated_k::CK_NCE;
using correlated_k::ck_nq;
using correlated_k::ck_nT;
using correlated_k::ck_nP;
using correlated_k::ck_lT_ptr;
using correlated_k::ck_lP_ptr;
using correlated_k::ck_lk_ptr;
using correlated_k::ck_gw_ptr;
using correlated_k::ck_wl_ptr;
using correlated_k::ck_swf_ptr;
using correlated_k::ck_pf_ptr;
using correlated_k::ck_pf_lTmin;
using correlated_k::ck_pf_idlT;
using correlated_k::ce_nT;
using correlated_k::ce_nP;
using correlated_k::ce_lT_ptr;
using correlated_k::ce_lP_ptr;
using correlated_k::ce_ptr;
using correlated_k::cia_nT_ptr;
using correlated_k::cia_T_ptr;
using correlated_k::cia_k_ptr;
using correlated_k::ray_x_ptr;
using correlated_k::planck_fraction_below;
using correlated_k::build_planck_fractions;
using correlated_k::ck_planck_frac;
using correlated_k::ck_planck_bands;
using correlated_k::read_ck_table;
using correlated_k::read_ck_continuum;
using correlated_k::ck_tp_index;
using correlated_k::ck_kappa;
using correlated_k::ck_continuum;
using correlated_k::ck_selftest;
using correlated_k::ck_rt_selftest;
using correlated_k::ck_build_rosseland_table;


// problem/bc_outer_maxwell: whether the outer-x1 ghost extrapolation carries the
// divergence of the Maxwell stress.
//
// The extrapolation is hydrostatic, and this term was added so the ghost sees the
// magnetic force as well as gravity. It enters as e0 -= e_i*dM1mag/(rho_i*grav_acc), so
// its size relative to the hydrostatic term is (B^2/dr)/(rho g). At the outermost
// shell of the standard setup (rho = 4.0e-10, dr = 5.6e7, g = 942) that is 0.03 at 3 G,
// 0.33 at 10 G, 3.0 at 30 G and 33 at 100 G. Past O(1) it is no longer a correction: it
// swings the ghost density by its own size and the last two active radial shells go.
//
// It is now applied as an effective gravity inside the hydrostatic solve instead, and
// clamped, so it is bounded at any field strength -- see the use site. DEFAULT ON, since
// in that form it is the better physics; set false to drop the magnetic force entirely.
inline bool bc_outer_maxwell = true;
// problem/top_bc (topbc-1101): the outer-x1 boundary.
//   hse      (DEFAULT, bitwise) hydrostatic isothermal extrapolation of rho, e in the
//            full potential phicc0 (gravity + centrifugal + tide), momenta scaled with rho,
//            v_r copied from the last active cell (inflow allowed).  The right top when
//            the WHOLE top shell is inside the planet's Roche lobe (closed equipotentials).
//   outflow1 ONE-WAY OUTFLOW, for a top on OPEN equipotentials (above Phi(L1), where no
//            hydrostatic state exists; the start-up line "top shell on OPEN
//            equipotentials" gives the fraction).  Two parts:
//            (a) ghosts: pure zero-gradient copy of rho, e and the tangential momenta of
//                the last active cell, v_r_ghost = max(v_r_ie, 0); no hydrostatic or
//                equipotential extrapolation anywhere (in the closed substellar cap of
//                WASP-121b the hse ghosts, with a tiny effective g, were a dense
//                reservoir that fed inflow).  The ghost total energy carries the ghost's
//                OWN rho*phicc0 under etotgrav and its own magnetic energy.  MHD face
//                field: tangential B zero-gradient (copied from the last active cell),
//                normal B of the ghost faces from div B = 0 cell by cell outward.
//            (b) the one-way top face (Hydro/MHD diode_ox1): where the interior-side
//                state at the ie+1 face moves inward, the Riemann solver sees its mirror
//                image (a reflecting face: zero mass and advected-energy flux, pressure
//                only); where it moves outward, the ghost state of (a).
//   outflow0 (a) alone -- DIAGNOSTIC ONLY, NOT one-way: the ghost pressure p_ie at the
//            top face exceeds the hydrostatic face value, pushes the top cell inward
//            (v_r ~ -900 m/s over the whole shell) and the Riemann solver refills it from
//            the ghost (C32 WASP-121b smoke, topbc_1101: Mdot_top -7e11 g/s, 12x the hse
//            inflow; -2.4e10 g/s even without the tide).
inline int dhj_top_bc = 0;  // 0 = hse, 1 = outflow1 (ghosts + one-way face), 2 = outflow0

// ---- problem/rt_test_* : TEST HOOKS for the well-posed checks of the implicit
// correlated-k solve (tests_ck_implicit/wellposed/README.md).  All default off, and off
// they change nothing: no hook is enrolled, no array is allocated and the only code on
// the default path is one untaken branch in SourceFunc.
//  rt_test_freeze  the hydro is FROZEN: the RK stages still run, but after the last one
//                  the state is reset to the frozen one and ONLY the radiation is
//                  applied to it, once per cycle (user_split_once), with dt =
//                  rt_test_dt if > 0 (else the hydro dt).  rho and the momenta never
//                  change; the internal energy changes by the radiation alone.
//  rt_test_rho, rt_test_T   (> 0) a horizontally AND radially uniform initial state
//                  [g/cm^3, K] (only meaningful with rt_test_freeze).
//  rt_test_dT, rt_test_nwave   initial perturbation T -> T (1 + dT sin(nwave pi xi)),
//                  xi = (r - x1min)/(x1max - x1min); nwave = 0 is a uniform dT.
//  rt_test_out, rt_test_every  per-cycle diagnostic file (serial runs only).
namespace dhj_rt_local {   // tusplit-1102: was an anonymous namespace
// problem/dtloc_every = N > 0: every N cycles print the cell that set the hydro/MHD dt
// (gid, i, j, k) with its pressure from the ck sweep.  Print only; 0 (default) = off.
inline int dtloc_every = 0;
inline bool rt_test_freeze = false;
inline Real rt_test_dt = 0.0;
inline std::string rt_test_out;
inline int rt_test_every = 1;
inline int rt_test_col[3] = {0, -1, -1};
inline DvceArray5D<Real> *rt_frz_ptr = nullptr;
inline DvceArray4D<Real> *rt_frzT_ptr = nullptr;
inline std::int64_t rt_test_ncall = 0;
inline std::FILE *rt_test_fp = nullptr;
}  // namespace dhj_rt_local
using dhj_rt_local::dtloc_every;
using dhj_rt_local::rt_test_freeze;
using dhj_rt_local::rt_test_dt;
using dhj_rt_local::rt_test_out;
using dhj_rt_local::rt_test_every;
using dhj_rt_local::rt_test_col;
using dhj_rt_local::rt_frz_ptr;
using dhj_rt_local::rt_frzT_ptr;
using dhj_rt_local::rt_test_ncall;
using dhj_rt_local::rt_test_fp;
void DhjRtFreezeSplit(Mesh *pm, const Real dt);

//----------------------------------------------------------------------------------------
//! \fn KOKKOS_INLINE_FUNCTION std::uint64_t SplitMix64Mix()
//! \brief splitmix64 finalizer: an integer hash built from xor-shifts and multiplies.
//! Used to turn (seed, global cell index) into a reproducible pseudo-random number, so
//! the initial-condition perturbation below is a pure function of the GLOBAL cell,
//! independent of the MeshBlock decomposition and of the number of MPI ranks.

KOKKOS_INLINE_FUNCTION
std::uint64_t SplitMix64Mix(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for the shallow hot Jupiter test


//----------------------------------------------------------------------------------------
// tusplit-1102: the device and IC-profile helpers (moved from the end of the old
// single file; every piece that builds an initial or boundary state reads them)

KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td) {
    
    Real bar = 1.0e6;
    Real pl = 1.0e-3*bar;
    Real pt = log10(p/bar);
    Real ptl = log10(pl/bar);
    
    Real fn[13];
    Real fd[14];
    
    fn[0] = 1388.77348;
    fn[1] = 279.575848;
    fn[2] = -213.835822;
    fn[3] = 21.0010475;
    fn[4] = 100.938036;
    fn[5] = 12.7972336;
    fn[6] = -13.9266925;
    fn[7] = -3.70783272;
    fn[8] = 0.522370269;
    fn[9] = 0.320837882;
    fn[10]= 0.0451831612;
    fn[11]= 2.18195583e-3;
    fn[12]= 3.98938097e-6;

    fd[0] = 2152.06036;
    fd[1] = 29.3485512;
    fd[2] = -183.318696;
    fd[3] = 46.3893130;
    fd[4] = 19.8116485;
    fd[5] = -28.5473177;
    fd[6] = -2.52726545;
    fd[7] = 8.43627538;
    fd[8] = 2.62945375;
    fd[9] = -0.297098168;
    fd[10]= -0.286871487;
    fd[11]= -0.0590629443;
    fd[12]= -5.38679474e-3;
    fd[13]= -1.89972415e-4;
    
    Real Tnstar = 0.0;
    Real Tdstar = 0.0;
    Real Tnstar_pl = 0.0;
    Real Tdstar_pl = 0.0;
    Real ptn = 1.0;
    Real ptln = 1.0;
    for(int ilogp=0; ilogp<13; ilogp++) {
        Tnstar += fn[ilogp]*ptn;
        Tnstar_pl += fn[ilogp]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }
    ptn = 1.0;
    ptln = 1.0;
    for(int ilogp=0; ilogp<14; ilogp++) {
        Tdstar += fd[ilogp]*ptn;
        Tdstar_pl += fd[ilogp]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }

    Tn = Tnstar;
    Td = Tdstar;
//    if (p < pl) {
//        Tn = Tnstar_pl * exp(0.1*log10(p/pl));
//        if (Tn < 250) {
//            Tn = 250;
//        }
//        Td = Tdstar_pl * exp(0.015*log10(p/pl));
//        if (Td < 1000) {
//            Td = 1000;
//        }
//    }
//    if (Td < Tn) {
//        Td = Tn;
//    }
    Real x = log10(p / pl);

    Real Tn_new = Tnstar_pl;// * exp(0.1   * x);
    Real Td_new = Tdstar_pl;// * exp(0.015 * x);

    Tn_new = fmax(Tn_new, 250.0);
    Td_new = fmax(Td_new, 1000.0);

    // blend
    Tn = (p < pl) ? Tn_new : Tn;
    Td = (p < pl) ? Td_new : Td;

    // enforce Td >= Tn
    Td = fmax(Td, Tn);

    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T,
                 const int &Ne, const DvceArray1D<Real> &Text, const DvceArray1D<Real> &lgpext) {
    // problem/ic_adiabat_ext: below the table bottom, the adiabat (ic_adiabat_extension);
    // Ne = 0 or inside the table: the original lookup below, unchanged
    if (Ne > 1 && log10(p) > lgparr(N-1)) {
      const Real lge = log10(p);
      const Real dle = (lgpext(Ne-1) - lgpext(0))/(Ne-1);
      int k = static_cast<int>(std::floor((lge - lgpext(0))/dle));
      k = (k < 0) ? 0 : ((k > Ne-2) ? Ne-2 : k);
      T = Text(k) + (Text(k+1) - Text(k))/(lgpext(k+1) - lgpext(k))*(lge - lgpext(k));
      return;
    }
    
    Real lgp = log10(p);
    Real dlgp = (lgparr(N-1)-lgparr(0))/(N-1);
    int Nt = std::floor((lgp - lgparr(0))/dlgp);
    int NN = (Nt < 0) ? 0 : ((Nt > N-2) ? N-2 : Nt);   // never past the table end
//    for (int it=Nt-2; it<Nt+3; ++it)
//    {
//        if (lgp < lgparr(it) && lgp >= lgparr(it-1)) {
//            NN = it-1;
//            break;
//        }
//    }
    T = Tarr(NN) + (Tarr(NN+1)-Tarr(NN))/(lgparr(NN+1)-lgparr(NN))*(lgp-lgparr(NN));
    T = (Nt < 0) ? Tarr(0) : T;
    
//    Real Tn, Td;
//    get_daynight_Tp(p, Tn, Td);
//
////    T = Tn;
//
//    Real Tn2 = Tn*Tn;
//    Real Tn4 = Tn2*Tn2;
//    Real Td2 = Td*Td;
//    Real Td4 = Td2*Td2;
//    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;
//
//    Real Tmid = sqrt(sqrt(Tmid4));
//    T = Tmid;
//
////    Real bar = 1.0e6;
////    Real Teq = 1469.0;
//////    Real Tirr = Teq*sqrt(2);
//////    Real Tint = 100.0;
//////    Real g = 942.0;
//////    Real mus = cos(50.0/90.0*M_PI/2.0);
//////    Real fH = 0.5;
//////    Real fK = 1.0/3.0;
//////    Real kap_v = 4.0e-3;
//////    Real kap_ir = 1.0e-2;
//////    Real gam = kap_v/kap_ir;
//////    Real tau = p/(g/kap_ir);
//////    tau = (tau < 0.0) ? 0.0 : tau;
//////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK) + 0.25*SQR(SQR(Tirr))*(mus/fH+SQR(mus)/gam/fK+(gam-SQR(mus)/gam/fK)*exp(-gam*tau/mus));
////////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK);
//////    T = sqrt(sqrt(T4)); //2581.5574;
////    T = Teq;
    
    
    return;
}

template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T,
                 const int &Ne, const View1D &Text, const View1D &lgpext) {
    // problem/ic_adiabat_ext: below the table bottom, the adiabat (ic_adiabat_extension);
    // Ne = 0 or inside the table: the original lookup below, unchanged
    if (Ne > 1 && log10(p) > lgparr(N-1)) {
      const Real lge = log10(p);
      const Real dle = (lgpext(Ne-1) - lgpext(0))/(Ne-1);
      int k = static_cast<int>(std::floor((lge - lgpext(0))/dle));
      k = (k < 0) ? 0 : ((k > Ne-2) ? Ne-2 : k);
      T = Text(k) + (Text(k+1) - Text(k))/(lgpext(k+1) - lgpext(k))*(lge - lgpext(k));
      return;
    }
    
    Real lgp = log10(p);
    Real dlgp = (lgparr(N-1)-lgparr(0))/(N-1);
    int Nt = std::floor((lgp - lgparr(0))/dlgp);
    int NN = (Nt < 0) ? 0 : ((Nt > N-2) ? N-2 : Nt);   // never past the table end
//    for (int it=Nt-2; it<Nt+3; ++it)
//    {
//        if (lgp < lgparr(it) && lgp >= lgparr(it-1)) {
//            NN = it-1;
//            break;
//        }
//    }
    T = Tarr(NN) + (Tarr(NN+1)-Tarr(NN))/(lgparr(NN+1)-lgparr(NN))*(lgp-lgparr(NN));
    T = (Nt < 0) ? Tarr(0) : T;
    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T) {
    
    Real Tn, Td;
    get_daynight_Tp(p, Tn, Td);
    
//    T = Tn;
    
    Real Tn2 = Tn*Tn;
    Real Tn4 = Tn2*Tn2;
    Real Td2 = Td*Td;
    Real Td4 = Td2*Td2;
    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;

    Real Tmid = sqrt(sqrt(Tmid4));
    T = Tmid;
    
//    Real bar = 1.0e6;
//    Real Teq = 1469.0;
////    Real Tirr = Teq*sqrt(2);
////    Real Tint = 100.0;
////    Real g = 942.0;
////    Real mus = cos(50.0/90.0*M_PI/2.0);
////    Real fH = 0.5;
////    Real fK = 1.0/3.0;
////    Real kap_v = 4.0e-3;
////    Real kap_ir = 1.0e-2;
////    Real gam = kap_v/kap_ir;
////    Real tau = p/(g/kap_ir);
////    tau = (tau < 0.0) ? 0.0 : tau;
////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK) + 0.25*SQR(SQR(Tirr))*(mus/fH+SQR(mus)/gam/fK+(gam-SQR(mus)/gam/fK)*exp(-gam*tau/mus));
//////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK);
////    T = sqrt(sqrt(T4)); //2581.5574;
//    T = Teq;
    
    return;
}


template <typename View1D>
void get_wb_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc,
                    const Real &ap, const bool &grav_pmass, const int &N,
                    const Real &zmax, View1D zarr, View1D logparr) {
    
//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;
        // gravity at THIS height, so the integration is consistent with the source term
        const Real gz = GravAccAt(grav_acc, ap, ap + zarr(n), grav_pmass);
        const Real fac = gz/Rgas*dz;

        get_wb_Tp(p,T);

        // n+1 is guarded: logparr holds N entries, so the last pass of this loop used to
        // integrate one step PAST THE END of it. A one-element write off the end of a
        // Kokkos allocation is silent in a Release build -- it lands in whatever follows
        // -- but it is undefined behaviour, a Kokkos Debug build aborts on it, and the
        // value written was never read. zarr(n) above still has to be set for every n.
        if (n+1 < N) {
          if (eos.IsGeneral()) {
            // dln p/dz = rho g/p, closed with the EOS's (p,T) -> rho inversion. The ideal
            // branch is the same thing with rho = p/(Rgas T), kept in its original form.
            Real rho = DensFromPT(eos, Rgas, p, T);
            logparr(n+1) = logparr(n) + gz*dz*rho/p;
          } else {
            logparr(n+1) = logparr(n) + fac/T;
          }
        }
    }

    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {
    
//    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        // clamp to the last segment (linear extrapolation beyond it): cells above the
        // column top (z_eff > 1.1 (r1-r0) under rot_potential) read past the end of
        // logparr, which faulted on the GPU (ck-hitemp2-gpu). In range: unchanged.
        const int Nmax = zarr.extent_int(0) - 2;
        if (Nt > Nmax) Nt = Nmax;
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_wb_Tp(p,T);
        rho = DensFromPT(eos, Rgas, p, T);
    } else {
        Real T0;
        get_wb_Tp(p0,T0);
        Real rho0 = DensFromPT(eos, Rgas, p0, T0);
//        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}

template <typename View1D>
void get_init_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc,
                      const Real &ap, const bool &grav_pmass, const View1D &Tarr,
                      const View1D &lgparr, const int &N, const Real &zmax,
                      View1D zarr, View1D logparr, const Real zlow,
                      const int &Ne, const View1D &Text, const View1D &lgpext) {
    // zlow < 0 (problem/ic_profile only): the column also covers z in [zlow, 0) -- the
    // cells and inner ghosts below the 250 bar anchor, deepest on the equator under
    // rot_potential -- integrated DOWN from p0 through the same table T(p), so they sit
    // on the table's adiabat instead of get_init_eos's isothermal extension.  The node
    // n0 is exactly z = 0.  zlow = 0 (default) is the original loop below, unchanged.
    if (zlow < 0.0) {
      const Real bar_ = 1.0e6;
      const Real dzd = (zmax - zlow)/N;
      const int n0 = static_cast<int>(std::ceil(-zlow/dzd));
      for (int n=0; n<N; n++) zarr(n) = (n - n0)*dzd;
      logparr(n0) = std::log(250.0*bar_);
      for (int n=n0; n<N-1; n++) {           // up, as the original loop
        Real T;
        const Real p = exp(logparr(n));
        const Real gz = GravAccAt(grav_acc, ap, ap + zarr(n), grav_pmass);
        get_init_Tp_host(N, Tarr, lgparr, p, T, Ne, Text, lgpext);
        logparr(n+1) = logparr(n) + gz*dzd*DensFromPT(eos, Rgas, p, T)/p;
      }
      for (int n=n0; n>0; n--) {              // down
        Real T;
        const Real p = exp(logparr(n));
        const Real gz = GravAccAt(grav_acc, ap, ap + zarr(n), grav_pmass);
        get_init_Tp_host(N, Tarr, lgparr, p, T, Ne, Text, lgpext);
        logparr(n-1) = logparr(n) - gz*dzd*DensFromPT(eos, Rgas, p, T)/p;
      }
      return;
    }

//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;
        // gravity at THIS height, so the integration is consistent with the source term
        const Real gz = GravAccAt(grav_acc, ap, ap + zarr(n), grav_pmass);
        const Real fac = gz/Rgas*dz;

        get_init_Tp_host(N, Tarr, lgparr, p, T, Ne, Text, lgpext);

        // n+1 is guarded: logparr holds N entries, so the last pass of this loop used to
        // integrate one step PAST THE END of it. A one-element write off the end of a
        // Kokkos allocation is silent in a Release build -- it lands in whatever follows
        // -- but it is undefined behaviour, a Kokkos Debug build aborts on it, and the
        // value written was never read. zarr(n) above still has to be set for every n.
        if (n+1 < N) {
          if (eos.IsGeneral()) {
            // dln p/dz = rho g/p, closed with the EOS's (p,T) -> rho inversion. The ideal
            // branch is the same thing with rho = p/(Rgas T), kept in its original form.
            Real rho = DensFromPT(eos, Rgas, p, T);
            logparr(n+1) = logparr(n) + gz*dz*rho/p;
          } else {
            logparr(n+1) = logparr(n) + fac/T;
          }
        }
    }


    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p,
                  const int &Ne, const DvceArray1D<Real> &Text, const DvceArray1D<Real> &lgpext) {

//    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= zarr(0)) {   // zarr(0) = 0 unless the ic_profile column reaches below z = 0
        int Nt = std::floor((z - zarr(0))/dz);
        if (Nt > N-2) Nt = N-2;   // as in get_wb_eos: never read past the column top
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_init_Tp(N, Tarr, lgparr, p,T, Ne, Text, lgpext);
        rho = DensFromPT(eos, Rgas, p, T);
    } else {
        Real T0;
        get_init_Tp(N, Tarr, lgparr, p0,T0, Ne, Text, lgpext);
        Real rho0 = DensFromPT(eos, Rgas, p0, T0);
//        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}





KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3], const Real albedo_set) {
    
    Real Tirr4 = SQR(SQR(Tirr));
    Real Tint4 = SQR(SQR(Tint));
    Real Teq = Tirr/sqrt(2.0);
    
    Real Teff0 = sqrt(sqrt(Tint4+Tirr4/sqrt(3.0)));
    Real albedo;
    get_albedo(Teff0,grav,albedo);
    if (albedo_set >= 0.0) albedo = albedo_set;   // problem/albedo (observed Bond albedo)
    
    Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
    Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
    get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);
    
    Real R = gamir1/gamir2;
    Real gamp = gamir1 + gamir2 - SQR(gamir2)*R;
    taulim = sqrt(R/3.0)*sqrt(beta*SQR(R-1.0)-SQR(beta*(R-1.0))+R)/SQR(gamir1);
    Real At1 = SQR(gamir1)*log(1.0+1.0/(taulim*gamir1));
    Real At2 = SQR(gamir2)*log(1.0+1.0/(taulim*gamir2));
    
    Real a0 = 1.0/gamir1 + 1.0/gamir2;
    Real a1 = -1.0/(3.0*SQR(taulim))*(gamp/(1.0-gamp)*(gamir1+gamir2-2.0)/(gamir1+gamir2) + (gamir1+gamir2)*taulim - (At1+At2)*SQR(taulim));
    Real b0 = 1.0/(gamir1*gamir2/(gamir1-gamir2)*(At1-At2)/3.0 - SQR(gamir1*gamir2)/sqrt(3.0*gamp) - SQR(gamir1*gamir2)*gamir1*gamir2/(1.0-gamir1)/(1.0-gamir2)/(gamir1+gamir2));
    A = 1.0/3.0*(a0+a1*b0);
    B = -1.0/3.0*SQR(gamir1*gamir2)/gamp*b0;
    
//    Real T4 = 3.0/4.0*Tint4*(tau + A + B*exp(-tau/taulim));
    
    for (int iv=0; iv<3; ++iv) {
        Real gamv;
        if (iv==0) gamv = gamv1;
        if (iv==1) gamv = gamv2;
        if (iv==2) gamv = gamv3;
        gamvv[iv] = gamv;
        
        Real longf = (3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*(gamir1+gamir2) - 3.0*gamv*(6.0*SQR(gamir1*gamir2)-SQR(gamv)*(SQR(gamir1)+SQR(gamir2)));
        Real a2 = SQR(taulim)/(gamp*SQR(gamv)) * longf/(1.0-SQR(gamv*taulim));
        Real Av1 = SQR(gamir1)*log(1.0+gamv/gamir1);
        Real Av2 = SQR(gamir2)*log(1.0+gamv/gamir2);
        Real a3 = -SQR(taulim) * (3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*(Av1+Av2) / (gamp*SQR(gamv)*gamv*(1.0-SQR(gamv*taulim)));
        Real b1 = gamir1*gamir2*(3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*SQR(taulim) / (gamp*SQR(gamv)*(SQR(gamv*taulim)-1.0));
        Real b2 = 3.0*(gamir1+gamir2)*SQR(gamv)*gamv / ((3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv)));
        Real b3 = (Av2-Av1) / (gamv*(gamir1-gamir2));
        
        C[iv] = -1.0/3.0*(b0*b1*(1.0+b2+b3)*a1 + a2 + a3);
        D[iv] = -B*b1*(1.0+b2+b3);
        E[iv] = (3.0-SQR(gamv/gamir1))*(3.0-SQR(gamv/gamir2)) / (9.0*gamv*(SQR(gamv*taulim)-1.0));
//        T4 += 3.0/4.0*1.0/3.0*Tirr4*mus*(C + D*exp(-tau/taulim) + E*exp(-gamv*tau));
    }
    
//    T = sqrt(sqrt(T4));
    return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real taulim, const Real &A, const Real &B, const Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamv)[3], const Real &tau, Real &T) {
    
    Real Tirr4 = SQR(SQR(Tirr));
    Real Tint4 = SQR(SQR(Tint));
    Real T4 = 3.0/4.0*Tint4*(tau + A + B*exp(-tau/taulim));
    for (int iv=0; iv<3; ++iv) {
        T4 += 3.0/4.0*1.0/3.0*Tirr4*mus*(C[iv] + D[iv]*exp(-tau/taulim) + E[iv]*exp(-gamv[iv]*tau));
    }
    T = sqrt(sqrt(T4));
    return;
}

template <typename View1D>
void get_picket_fence_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr) {
    Real bar = 1.0e6;
    Real tautop = 1.0e-6;
    Real Ttop;
    Real taulim, A, B, C[3], D[3], E[3], gamv[3];
    get_picket_fence_Ttau_coeff(Tint, Tirr, met, grav, mus, taulim, A, B, C, D, E, gamv,
                                two_stream_rt::rt_albedo_set);
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tautop, Ttop);
    Real kapr, ptop;
    ptop = 2.0; // 2e-6 bar
    get_kapr(Ttop, ptop, met, kapr);
    tautop = kapr/grav*ptop;
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tautop, Ttop);
    
    Real lgptop = log10(ptop);
    Real pbot = 300.0*bar;
    Real lgpbot = log10(pbot);
    Real dlgp = (lgpbot-lgptop)/(N-1);
    Real lgp = lgptop;
    Real p = pow(10.0,lgp);
    Real tau = tautop;
    Real lgtau = log10(tau);
    Real T;
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tau, T);
    get_kapr(Ttop, p, met, kapr);
    for (int ip=0; ip<N; ++ ip) {
        Tarr(ip) = T;
        lgparr(ip) = lgp;
        Real K = p*kapr/grav*log(10.0);
        tau += K*dlgp;
        lgp += dlgp;
        p = pow(10.0,lgp);
//        lgp += dlgp;
//        Real pp = pow(10.0,lgp);
//        Real dp = pp-p;
//        Real K = kapr/grav;
//        tau += K*dp;
//        p = pp;
        get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tau, T);
        get_kapr(T, p, met, kapr);
    }
    
    adjust_ad_pT_arr(eos, Rgas, gamma, N, Tarr, lgparr);

    return;
}

//----------------------------------------------------------------------------------------
//! \fn void read_ic_profile
//! \brief problem/ic_profile: read a text table of "log10 p [BARYE]  T [K]" (one pair per
//! line, '#' comments, any order in p) and resample it onto the N-point array the
//! hydrostatic integrator reads, UNIFORM in log10 p (get_init_Tp_host assumes uniform
//! spacing), linear in log10 p.  The table must reach 250 bar (get_init_eos_arr starts
//! there and get_init_Tp_host has no clamp below the bottom of the array).  Everything
//! downstream -- the hydrostatic z(log p), the probini fill -- is unchanged.

template <typename View1D>
void read_ic_profile(const std::string &fname, const int &N, View1D Tarr, View1D lgparr) {
  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "problem/ic_profile: cannot open '" << fname << "'"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::vector<std::pair<Real, Real>> tab;
  std::string line;
  while (std::getline(fin, line)) {
    const std::size_t c = line.find('#');
    if (c != std::string::npos) line.erase(c);
    std::istringstream ls(line);
    Real lp, t;
    if (ls >> lp >> t) tab.emplace_back(lp, t);
  }
  std::sort(tab.begin(), tab.end());
  const int nt = static_cast<int>(tab.size());
  bool ok = (nt >= 2);
  for (int i=0; ok && i<nt; ++i) {
    if (!(tab[i].second > 0.0) || (i > 0 && !(tab[i].first > tab[i-1].first))) ok = false;
  }
  if (!ok || tab[nt-1].first < std::log10(250.0e6)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "problem/ic_profile '" << fname << "': need >= 2 rows of "
              << "strictly increasing log10 p [barye] with T > 0, reaching 250 bar "
              << "(log10 p >= 8.398)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real lg0 = tab[0].first, lg1 = tab[nt-1].first;
  const Real dlgp = (lg1 - lg0)/(N - 1);
  int it = 0;
  for (int ip=0; ip<N; ++ip) {
    const Real lgp = (ip == N-1) ? lg1 : lg0 + ip*dlgp;
    while (it < nt-2 && tab[it+1].first < lgp) ++it;
    const Real f = (lgp - tab[it].first)/(tab[it+1].first - tab[it].first);
    lgparr(ip) = lgp;
    Tarr(ip) = tab[it].second + f*(tab[it+1].second - tab[it].second);
  }
  std::cout << "deep_hot_jupiter_rt: initial T(p) from problem/ic_profile '" << fname
            << "' (" << nt << " rows, log10 p = " << lg0 << " .. " << lg1
            << " barye)" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn int ic_adiabat_extension
//! \brief problem/ic_adiabat_ext: T(log10 p) BELOW the resampled IC table (Tarr, lgparr,
//! N points), along the adiabat through its bottom row: dlog10 T/dlog10 p = GradAd(p, T)
//! of the deck's EOS (pgen_eos_utils.hpp; (gamma-1)/gamma for an ideal gas), classical
//! RK4 on Ne points uniform in log10 p from lgparr(N-1) (exactly; Text(0) = Tarr(N-1))
//! over `dex` decades.  Host, setup only (GradAd root finds the density).  Returns Ne
//! (0 if dex <= 0).  T is kept below 10^5.9 K (inside the general EOS table).

template <typename View1D>
int ic_adiabat_extension(const EOS_Data &eos, const Real &gamma, const Real &Rgas,
                         const int &N, const View1D &Tarr, const View1D &lgparr,
                         const Real &dex, const int &Ne, View1D Text, View1D lgpext) {
  if (!(dex > 0.0) || Ne < 2) return 0;
  const Real lg0 = lgparr(N-1), h = dex/(Ne - 1);
  const Real tmax = std::pow(10.0, eos.IsGeneral() ? 5.9 : 30.0);
  auto slope = [&](const Real lgp, const Real lt) {
    const Real tk = std::min(std::pow(10.0, lt), tmax);
    return GradAd(eos, gamma, Rgas, std::pow(10.0, lgp), tk);
  };
  lgpext(0) = lg0;
  Text(0) = Tarr(N-1);
  Real lt = std::log10(Tarr(N-1));
  for (int k=0; k<Ne-1; ++k) {
    const Real lgp = lg0 + k*h;
    const Real k1 = slope(lgp, lt);
    const Real k2 = slope(lgp + 0.5*h, lt + 0.5*h*k1);
    const Real k3 = slope(lgp + 0.5*h, lt + 0.5*h*k2);
    const Real k4 = slope(lgp + h, lt + h*k3);
    lt += h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    lgpext(k+1) = lg0 + (k+1)*h;
    Text(k+1) = std::min(std::pow(10.0, lt), tmax);
  }
  if (global_variable::my_rank == 0) {
    const Real gb = GradAd(eos, gamma, Rgas, std::pow(10.0, lg0), Tarr(N-1));
    const Real gt = (Tarr(N-1) > Tarr(N-2)) ?
        std::log10(Tarr(N-1)/Tarr(N-2))/(lgparr(N-1) - lgparr(N-2)) : 0.0;
    const int k2x = static_cast<int>(std::lround((Ne - 1)*std::log10(2.0)/dex));
    std::cout << "deep_hot_jupiter_rt: problem/ic_adiabat_ext: below the IC table bottom ("
              << std::pow(10.0, lg0)/1.0e6 << " bar, " << Tarr(N-1) << " K; table slope "
              << gt << ", nabla_ad " << gb << ") T follows the adiabat: "
              << std::pow(10.0, lgpext(k2x))/1.0e6 << " bar " << Text(k2x) << " K .. "
              << std::pow(10.0, lgpext(Ne-1))/1.0e6 << " bar " << Text(Ne-1) << " K"
              << std::endl;
  }
  return Ne;
}

template <typename View1D>
void adjust_ad_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const int &N, View1D Tarr, View1D lgparr) {

  // --- find convective boundary (search from bottom) ---
  int ic = -1;

  for (int ip = 0; ip < N-1; ++ip) {
    int i_inv = N - 1 - ip;
    int i_inv_p1 = N - 1 - (ip + 1);

    Real lgT  = log10(Tarr(i_inv));
    Real lgT1 = log10(Tarr(i_inv_p1));

    Real lgp  = lgparr(i_inv);
    Real lgp1 = lgparr(i_inv_p1);

    Real nabla = (lgT - lgT1) / (lgp - lgp1);

    Real T = Tarr(i_inv);
    Real nabla_ad = 0.9*GradAd(eos, gamma, Rgas, pow(10.0,lgp), T);

    if (nabla < nabla_ad) {
      ic = i_inv_p1;  // map back to original indexing
      break;
    }
  }

  // fallback if no crossing found
  if (ic > 0) {
    // --- enforce adiabat downward ---
    for (int ip = ic; ip < N-1; ++ip) {
      Real T = Tarr(ip);
      Real nabla_ad = 0.9*GradAd(eos, gamma, Rgas, pow(10.0,lgparr(ip)), T);

      Tarr(ip+1) = pow(10.0,
        nabla_ad * (lgparr(ip+1) - lgparr(ip)) + log10(T)
      );
    }
  }
    
  return;
}

#endif  // PGEN_DEEP_HOT_JUPITER_RT_HPP_
