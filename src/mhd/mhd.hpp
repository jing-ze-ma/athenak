#ifndef MHD_MHD_HPP_
#define MHD_MHD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd.hpp
//  \brief definitions for MHD class

#include <map>
#include <memory>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/cell_locations.hpp"

// forward declarations
class EquationOfState;
class Coordinates;
class Viscosity;
class Resistivity;
class Conduction;
class SourceTerms;
class OrbitalAdvectionCC;
class OrbitalAdvectionFC;
class ShearingBoxCC;
class ShearingBoxFC;
class Driver;

// function ptr for user-defined MHD boundary functions enrolled in problem generator
namespace mhd {
using MHDBoundaryFnPtr = void (*)(int m, Mesh* pm, MHD* pmhd, DvceArray5D<Real> &u);
}

// constants that enumerate MHD Riemann Solver options
enum class MHD_RSolver {advect, llf, hlle, hlld, roe,   // non-relativistic
                        llf_sr, hlle_sr,                // SR
                        llf_gr, hlle_gr};                       // GR

//----------------------------------------------------------------------------------------
//! \struct MHDTaskIDs
//  \brief container to hold TaskIDs of all mhd tasks

struct MHDTaskIDs {
  TaskID savest;
  TaskID irecv;
  TaskID copyu;
  TaskID flux;
  TaskID sendf;
  TaskID recvf;
  TaskID rkupdt;
  TaskID srctrms;
  TaskID sendu_oa;
  TaskID recvu_oa;
  TaskID restu;
  TaskID sendu;
  TaskID recvu;
  TaskID sendu_shr;
  TaskID recvu_shr;
  TaskID efld;
  TaskID sende;
  TaskID recve;
  TaskID ct;
  TaskID sendb_oa;
  TaskID recvb_oa;
  TaskID restb;
  TaskID sendb;
  TaskID recvb;
  TaskID sendb_shr;
  TaskID recvb_shr;
  TaskID bcs;
  TaskID prol;
  TaskID c2p;
  TaskID newdt;
  TaskID csend;
  TaskID crecv;
};

namespace mhd {

//----------------------------------------------------------------------------------------
//! \class MHD

class MHD {
 public:
  MHD(MeshBlockPack *ppack, ParameterInput *pin);
  ~MHD();

  // data
  ReconstructionMethod recon_method;
  MHD_RSolver rsolver_method;
  EquationOfState *peos;   // chosen EOS

  int nmhd;                // number of mhd variables (5/4 for ideal/isothermal EOS)
  int nscalars;            // number of passive scalars
  DvceArray5D<Real> u0;    // conserved variables
  DvceArray5D<Real> w0;    // primitive variables
  DvceFaceFld4D<Real> b0;  // face-centered magnetic fields
  DvceArray5D<Real> bcc0;  // cell-centered magnetic fields

  DvceArray5D<Real> coarse_u0;    // conserved variables on 2x coarser grid (for SMR/AMR)
  DvceArray5D<Real> coarse_w0;    // primitive variables on 2x coarser grid (for SMR/AMR)
  DvceFaceFld4D<Real> coarse_b0;  // face-centered B-field on 2x coarser grid

  // Objects containing boundary communication buffers and routines for u and b
  MeshBoundaryValuesCC *pbval_u;
  MeshBoundaryValuesFC *pbval_b;
  MHDBoundaryFnPtr MHDBoundaryFunc[6];

  // Orbital advection and shearing box BCs
  OrbitalAdvectionCC *porb_u = nullptr;
  OrbitalAdvectionFC *porb_b = nullptr;
  ShearingBoxCC *psbox_u = nullptr;
  ShearingBoxFC *psbox_b = nullptr;

  // Object(s) for extra physics (viscosity, resistivity, thermal conduction, srcterms)
  Viscosity *pvisc = nullptr;
  Resistivity *presist = nullptr;
  Conduction *pcond = nullptr;
  SourceTerms *psrc = nullptr;

  // following only used for time-evolving flow
  DvceArray5D<Real> u1;       // conserved variables, second register
  DvceFaceFld4D<Real> b1;     // face-centered magnetic fields, second register
  DvceFaceFld5D<Real> uflx;   // fluxes of conserved quantities on cell faces
  DvceEdgeFld4D<Real> efld;   // edge-centered electric fields (fluxes of B)
  // temporary variables used to store face-centered electric fields returned by RS
  DvceArray4D<Real> e3x1, e2x1;
  DvceArray4D<Real> e1x2, e3x2;
  DvceArray4D<Real> e2x3, e1x3;
  Real dtnew;

  // following used for time derivatives in computation of jcon
  bool wbcc_saved = false;
  DvceArray5D<Real> wsaved;
  DvceArray5D<Real> bccsaved;

  // following used for FOFC algorithm
  DvceArray4D<bool> fofc;  // flag for each cell to indicate if FOFC is needed
  bool use_fofc = false;   // flag to enable FOFC
    
  // following only used for including time-independent gravity in the conserved energy equation
  bool use_etotgrav = false;   // flag to enable etotgrav
  DvceFaceFld4D<Real> phi0;     // face-centered gravitational potential energy
  DvceArray4D<Real> phicc0;     // cell-centered gravitational potential energy

  // container to hold names of TaskIDs
  MHDTaskIDs id;

  // functions...
  void SetSaveWBcc();
  void AssembleMHDTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  // ...in "before_timeintegrator" task list
  TaskStatus SaveMHDState(Driver *d, int stage);
  // ...in "before_stagen_tl" task list
  TaskStatus InitRecv(Driver *d, int stage);
  // ...in "stagen_tl" task list
  TaskStatus CopyCons(Driver *d, int stage);
  TaskStatus Fluxes(Driver *d, int stage);
  TaskStatus SendFlux(Driver *d, int stage);
  TaskStatus RecvFlux(Driver *d, int stage);
  TaskStatus RKUpdate(Driver *d, int stage);
  TaskStatus MHDSrcTerms(Driver *d, int stage);
  TaskStatus SendU_OA(Driver *d, int stage);
  TaskStatus RecvU_OA(Driver *d, int stage);
  TaskStatus RestrictU(Driver *d, int stage);
  TaskStatus SendU(Driver *d, int stage);
  TaskStatus RecvU(Driver *d, int stage);
  TaskStatus SendU_Shr(Driver *d, int stage);
  TaskStatus RecvU_Shr(Driver *d, int stage);
  TaskStatus CornerE(Driver *d, int stage);
  TaskStatus EField(Driver *d, int stage);
  TaskStatus SendE(Driver *d, int stage);
  TaskStatus RecvE(Driver *d, int stage);
  TaskStatus CT(Driver *d, int stage);
  TaskStatus SendB_OA(Driver *d, int stage);
  TaskStatus RecvB_OA(Driver *d, int stage);
  TaskStatus RestrictB(Driver *d, int stage);
  TaskStatus SendB(Driver *d, int stage);
  TaskStatus RecvB(Driver *d, int stage);
  TaskStatus SendB_Shr(Driver *d, int stage);
  TaskStatus RecvB_Shr(Driver *d, int stage);
  TaskStatus ApplyPhysicalBCs(Driver* pdrive, int stage);
  TaskStatus Prolongate(Driver* pdrive, int stage);
  TaskStatus ConToPrim(Driver *d, int stage);
  TaskStatus NewTimeStep(Driver *d, int stage);
  // ...in "after_stagen_tl" task list
  TaskStatus ClearSend(Driver *d, int stage);
  TaskStatus ClearRecv(Driver *d, int stage);  // also in Driver::Initialize

  // CalculateFluxes function templated over Riemann Solvers
  template <MHD_RSolver T>
  void CalculateFluxes(Driver *d, int stage);

  // first-order flux correction
  void FOFC(Driver *d, int stage);
    
  void AddGravFlux(const DvceFaceFld4D<Real> &phi0, DvceFaceFld5D<Real> &flx);
  void AddGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons, const int il, const int iu, const int jl, const int ju, const int kl, const int ku);
  void RemoveGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons, const int il, const int iu, const int jl, const int ju, const int kl, const int ku);

  DvceArray5D<Real> utest, bcctest;  // scratch arrays for FOFC
    
    KOKKOS_INLINE_FUNCTION
    void PLM_nonuniform(const Real &q_im1, const Real &q_i, const Real &q_ip1,
                        const Real &dxL, const Real &dxR, const Real &dxLh, const Real &dxRh,
                        Real &ql_ip1, Real &qr_i) {
      Real cL = dxL/dxLh;
      Real cR = dxR/dxRh;

      // Left/right slopes (properly scaled)
      Real sL = (q_i   - q_im1) / dxL;
      Real sR = (q_ip1 - q_i  ) / dxR;

      // modified van Leer (Mignone 2014)
      Real slope = sL*sR*(cR*sL+cL*sR) / (SQR(sL)+(cR+cL-2.0)*sL*sR+SQR(sR));
      slope = (sL * sR > 0.0) ? slope : 0.0;

      // Reconstruct to faces
      ql_ip1 = q_i + slope * dxRh;
      qr_i   = q_i - slope * dxLh;
      return;
    };
    
    KOKKOS_INLINE_FUNCTION
    void GridPiecewiseLinearX1(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
        const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf,
        ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
      int nvar = q.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real x_im1  = xv(m,i-1);
          Real x_imh  = xf(m,i);
          Real x_i    = xv(m,i);
          Real x_iph  = xf(m,i+1);
          Real x_ip1  = xv(m,i+1);
          Real dxLh = x_i-x_imh;
          Real dxRh = x_iph-x_i;
          Real dxL = x_i-x_im1;
          Real dxR = x_ip1-x_i;
                
          PLM_nonuniform(q(m,n,k,j,i-1), q(m,n,k,j,i), q(m,n,k,j,i+1), dxL, dxR, dxLh, dxRh, ql(n,i+1), qr(n,i));
        });
      }
      return;
    }
    
    KOKKOS_INLINE_FUNCTION
    void GridPiecewiseLinearX2(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
        const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf,
        ScrArray2D<Real> &ql_jp1, ScrArray2D<Real> &qr_j) {
      int nvar = q.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real x_jm1  = xv(m,j-1);
          Real x_jmh  = xf(m,j);
          Real x_j    = xv(m,j);
          Real x_jph  = xf(m,j+1);
          Real x_jp1  = xv(m,j+1);
          Real dxLh = x_j-x_jmh;
          Real dxRh = x_jph-x_j;
          Real dxL = x_j-x_jm1;
          Real dxR = x_jp1-x_j;
                
          PLM_nonuniform(q(m,n,k,j-1,i), q(m,n,k,j,i), q(m,n,k,j+1,i), dxL, dxR, dxLh, dxRh, ql_jp1(n,i), qr_j(n,i));
        });
      }
      return;
    }
    
    KOKKOS_INLINE_FUNCTION
    void GridPiecewiseLinearX3(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
        const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf,
        ScrArray2D<Real> &ql_kp1, ScrArray2D<Real> &qr_k) {
      int nvar = q.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real x_km1  = xv(m,k-1);
          Real x_kmh  = xf(m,k);
          Real x_k    = xv(m,k);
          Real x_kph  = xf(m,k+1);
          Real x_kp1  = xv(m,k+1);
          Real dxLh = x_k-x_kmh;
          Real dxRh = x_kph-x_k;
          Real dxL = x_k-x_km1;
          Real dxR = x_kp1-x_k;
                
          PLM_nonuniform(q(m,n,k-1,j,i), q(m,n,k,j,i), q(m,n,k+1,j,i), dxL, dxR, dxLh, dxRh, ql_kp1(n,i), qr_k(n,i));
        });
      }
      return;
    }

 private:
  MeshBlockPack* pmy_pack;   // ptr to MeshBlockPack containing this MHD
  // temporary variables used to store face-centered electric fields returned by RS
  DvceArray4D<Real> e1_cc, e2_cc, e3_cc;
};

} // namespace mhd
#endif // MHD_MHD_HPP_
