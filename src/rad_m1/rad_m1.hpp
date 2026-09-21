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
  TaskID update;
  TaskID coupl;
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
  std::string thick_flux_str;   // none | ap_hll | scaled; only "none" is implemented

  // true when this module's dtnew must be folded into Mesh::NewTimeStep, i.e. when
  // sub-cycling is off OR when no other module sets the mesh timestep
  bool sets_mesh_dt;

  // evolved variables
  DvceArray5D<Real> u0;         // (E, F1, F2, F3)
  DvceArray5D<Real> u1;         // state at the start of the substep
  DvceArray5D<Real> coarse_u0;  // on 2x coarser grid (for SMR/AMR)
  DvceFaceFld5D<Real> uflx;     // fluxes on cell faces

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

  // ...in "m1_before_stagen"
  TaskStatus InitRecv(Driver *d, int stage);
  // ...in "m1_stagen"
  TaskStatus CopyCons(Driver *d, int stage);
  TaskStatus ApplyClosureLimits(Driver *d, int stage);
  TaskStatus CalculateFluxes(Driver *d, int stage);
  TaskStatus SendFlux(Driver *d, int stage);
  TaskStatus RecvFlux(Driver *d, int stage);
  TaskStatus Update(Driver *d, int stage);
  TaskStatus Coupling(Driver *d, int stage);
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
