#ifndef RAD_M1_RAD_M1_HPP_
#define RAD_M1_RAD_M1_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1.hpp
//! \brief definitions for the grey photon two-moment (M1) radiation module, <rad_m1>.
//! See docs/dev/rad_m1_design.md.  MILESTONE 1a: explicit transport only -- the closure,
//! the HLL flux with thick_flux = none, the PD-ARS time integration with a null source,
//! sub-cycling and the cell-centred plumbing.  The implicit matter coupling (design
//! sect. 4) is milestone 1b and enters at the marked Coupling task in both stages.
//!
//! MILESTONE 1b adds: analytic (constant / power-law / user-hook) opacities stored as
//! rho*kappa per cell, the two thick-limit fluxes (ap_hll, scaled), and the implicit
//! local matter coupling to <hydro>, called inside BOTH PD-ARS stages.  Still absent
//! (milestone 1c): the advective enthalpy-flux split of the E equation for moving media.

#include <map>
#include <memory>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"

// forward declarations
class Driver;

namespace radm1 {

// indices of the evolved moments in u0(m,n,k,j,i)
constexpr int M1_E  = 0;
constexpr int M1_F1 = 1;
constexpr int M1_F2 = 2;
constexpr int M1_F3 = 3;
constexpr int M1_NVAR = 4;

// number of explicit stages of the PD-ARS (Chu et al. 2019) integrator used per substep.
// With S == 0 it reduces to Heun / SSP-RK2.
constexpr int M1_NSTAGE = 2;

// components of the per-cell opacity array, all stored as rho*kappa (inverse length)
constexpr int M1_OP_P = 0;   // rho kappa_P, Planck (emission) mean
constexpr int M1_OP_E = 1;   // rho kappa_E, energy (absorption) mean
constexpr int M1_OP_T = 2;   // rho (kappa_F + kappa_s), the TRANSPORT opacity: what the
                             // flux source relaxes F with, and what tau_face is built of
constexpr int M1_NOPAC = 3;

// <rad_m1>/opacity
constexpr int M1_OPAC_CONST    = 0;
constexpr int M1_OPAC_POWERLAW = 1;
constexpr int M1_OPAC_USER     = 2;

// entries of the implicit-solve diagnostic counter (see RadiationM1::cnt)
constexpr int M1_CNT_NSOLVE = 0;   // cells handed to the bracketed Newton
constexpr int M1_CNT_ITSUM  = 1;   // total iterations over those cells
constexpr int M1_CNT_ITMAX  = 2;   // largest iteration count of any one cell
constexpr int M1_CNT_NFAIL  = 3;   // cells that hit maxit without reaching the tolerance
constexpr int M1_CNT_NBRAK  = 4;   // cells for which no bracket [T_lo, T_hi] was found
constexpr int M1_NCNT = 5;

// bracketed-Newton controls (design sect. 4a)
constexpr int  M1_MAXIT = 100;
constexpr Real M1_RTOL  = 1.0e-12;

//----------------------------------------------------------------------------------------
//! \struct RadiationM1TaskIDs
//! \brief container to hold TaskIDs of all rad_m1 tasks

struct RadiationM1TaskIDs {
  TaskID irecv;
  TaskID copyu;
  TaskID closure;
  TaskID flux;
  TaskID sendf;
  TaskID recvf;
  TaskID opac;
  TaskID update;
  TaskID coupl;
  TaskID c2p;
  TaskID restu;
  TaskID sendu;
  TaskID recvu;
  TaskID bcs;
  TaskID prol;
  TaskID csend;
  TaskID crecv;
};

//----------------------------------------------------------------------------------------
//! \class RadiationM1

class RadiationM1 {
 public:
  RadiationM1(MeshBlockPack *ppack, ParameterInput *pin);
  ~RadiationM1();

  // parameters read from <rad_m1>
  Real c_light;        // speed of light in code units (required)
  Real chat_over_c;    // reduced speed of light factor (default 1)
  Real chat;           // = chat_over_c * c_light
  Real cfl_rad;        // CFL number for the radiation substep (default 0.4)
  Real e_floor;        // floor on E
  bool subcycle;       // sub-cycle the module inside the mesh timestep
  bool eddington;      // closure = eddington (chi = 1/3); default is the M1 closure
  std::string thick_flux_str;   // none | ap_hll | scaled
  int thick_flux;               // M1_THICK_*
  Real scaled_pref;             // thick_flux = scaled: tau_c = scaled_pref*tau_face

  // opacities, per unit mass, in code units (design sect. 4).  kappa = kappa0 for
  // opacity = const, kappa0 (rho/rho_ref)^opac_a (T/t_ref)^opac_b for powerlaw.
  int opacity_type;             // M1_OPAC_*
  Real kappa_p, kappa_e, kappa_f, kappa_s;
  Real opac_rho_ref, opac_t_ref, opac_a, opac_b;
  bool opac_zero;               // all four opacities are identically zero

  // matter coupling (design sect. 4)
  bool coupling;       // run the implicit local solve at all
  bool gas_feedback;   // write the solve's reaction back into hydro u0.  With this off
                       // the medium is a PRESCRIBED static background: the radiation
                       // still feels the opacity, the gas does not move or heat.  That
                       // is what the static-medium tests T3/T3b/T3c mean by "no hydro".
  Real arad;           // radiation constant in CODE units: the equilibrium energy
                       // density is arad*T^4 with T the EOS's code temperature

  // true when this module's dtnew must be folded into Mesh::NewTimeStep, i.e. when
  // sub-cycling is off OR when no other module sets the mesh timestep
  bool sets_mesh_dt;

  // evolved variables
  DvceArray5D<Real> u0;         // (E, F1, F2, F3)
  DvceArray5D<Real> u1;         // state at the start of the substep
  DvceArray5D<Real> coarse_u0;  // on 2x coarser grid (for SMR/AMR)
  DvceFaceFld5D<Real> uflx;     // fluxes on cell faces
  DvceArray5D<Real> opac;       // rho*kappa per cell, M1_NOPAC components
  DvceArray5D<Real> ugas1;      // hydro (IEN, IM1, IM2, IM3) at the start of the substep

  // implicit-solve diagnostics, accumulated on the device and printed at the end of the
  // run (M1_CNT_*).  Reals rather than ints so that the sums cannot wrap on a long run.
  DualArray1D<Real> cnt;

  // boundary communication
  MeshBoundaryValuesCC *pbval_u;

  // timestep bookkeeping
  Real dtnew;     // cfl_rad*min(dx)/chat, over this pack
  Real dt_sub;    // the dt of the substep currently being taken
  int nsub;       // number of substeps in the current mesh step

  // PD-ARS explicit-stage weights, u0 = gam0*u0 + gam1*u1 - beta*dt_sub*div(F)
  Real gam0[M1_NSTAGE], gam1[M1_NSTAGE], beta[M1_NSTAGE];

  // reconstruction method (PLM only in milestone 1a)
  ReconstructionMethod recon_method;

  RadiationM1TaskIDs id;

  // functions
  void AssembleRadM1Tasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  //! set dt_sub and nsub for a mesh step of length dt_mesh; returns nsub
  int SetSubsteps(Real dt_mesh);
  //! print the implicit-solve counters (called from the destructor, rank 0)
  void ReportCounters();

  // ...in "m1_before_stagen"
  TaskStatus InitRecv(Driver *d, int stage);
  // ...in "m1_stagen"
  TaskStatus CopyCons(Driver *d, int stage);
  TaskStatus ApplyClosureLimits(Driver *d, int stage);
  TaskStatus Opacity(Driver *d, int stage);
  TaskStatus CalculateFluxes(Driver *d, int stage);
  TaskStatus SendFlux(Driver *d, int stage);
  TaskStatus RecvFlux(Driver *d, int stage);
  TaskStatus Update(Driver *d, int stage);
  TaskStatus Coupling(Driver *d, int stage);
  TaskStatus HydroConToPrim(Driver *d, int stage);
  TaskStatus RestrictU(Driver *d, int stage);
  TaskStatus SendU(Driver *d, int stage);
  TaskStatus RecvU(Driver *d, int stage);
  TaskStatus ApplyPhysicalBCs(Driver *d, int stage);
  TaskStatus Prolongate(Driver *d, int stage);
  TaskStatus NewTimeStep(Driver *d, int stage);
  // ...in "m1_after_stagen"
  TaskStatus ClearSend(Driver *d, int stage);
  TaskStatus ClearRecv(Driver *d, int stage);

 private:
  MeshBlockPack *pmy_pack;
};

} // namespace radm1
#endif // RAD_M1_RAD_M1_HPP_
