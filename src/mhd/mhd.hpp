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
#include "reconstruct/plm.hpp"
#include "eos/eos.hpp"
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
    
  // following used for well-balanced scheme
  bool use_wellbalance_static = false;    // flag to enable static wellbalance
  bool use_wellbalance_static_reconst_perturb = false;    // flag to enable reconstructing perturbed primitive variables (less robust against large deviations)
  bool use_wellbalance_dynamic = false;    // flag to enable dynamical well-balanced method by Kappeli & Mishra 2014, 2016
  bool use_wb_x1 = false;    // flag for directions
  bool use_wb_x2 = false;    // flag for directions
  bool use_wb_x3 = false;    // flag for directions
  bool use_wb_rho = false;   // flag to enable local well-balanced method also in reconstructing rho
  WBOption wb_option;
  DvceArray5D<Real> u0wb;   // background conserved variables
  DvceArray5D<Real> w0wb;   // background primitive variables
  DvceFaceFld5D<Real> w0facewb;   // face-centered background primitive variables

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

  DvceArray5D<Real> utest, bcctest;  // scratch arrays for FOFC
    
    void PolarAzimuthalAverageEr(void);
    
    void AddGravFlux(const DvceFaceFld4D<Real> &phi0, DvceFaceFld5D<Real> &flx);
    void AddGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons, const int il, const int iu, const int jl, const int ju, const int kl, const int ku);
    void RemoveGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons, const int il, const int iu, const int jl, const int ju, const int kl, const int ku);
    
    void RemoveWbFlux(const DvceFaceFld5D<Real> &w0facewb, DvceFaceFld5D<Real> &flx);
    void AddWbVar(const DvceArray5D<Real> &varwb, DvceArray5D<Real> &var);
    void RemoveWbVar(const DvceArray5D<Real> &varwb, DvceArray5D<Real> &var);

      //----------------------------------------------------------------------------------------
      //! \fn AddWbPrimFace()
      //! \brief Adds background face-centered variables onto ql(i+1) and qr(i).

      KOKKOS_INLINE_FUNCTION
      void AddWbPrimFace(const Real &qlwb_ip1, const Real &qrwb_i,
               Real &ql_ip1, Real &qr_i) {
        ql_ip1 += qlwb_ip1;
        qr_i   += qrwb_i;
        return;
      };

      //----------------------------------------------------------------------------------------
      //! \fn AddWbPrimFaceX1()
      //! \brief Adds background states onto face-centered primitive variables in x1-direction.
      //! This function should be called over [is-1,ie+1] to get BOTH L/R states over [is,ie]

      KOKKOS_INLINE_FUNCTION
      void AddWbPrimFaceX1(TeamMember_t const &member, const int m, const int k, const int j,
           const int il, const int iu, const DvceArray5D<Real> &q,
           ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
        int nvar = q.extent_int(1);
        for (int n=0; n<nvar; ++n) {
          par_for_inner(member, il, iu, [&](const int i) {
            AddWbPrimFace(q(m,n,k,j,i+1), q(m,n,k,j,i), ql(n,i+1), qr(n,i));
          });
        }
        return;
      };

      //----------------------------------------------------------------------------------------
      //! \fn AddWbPrimFaceX2()
      //! \brief Adds background states onto face-centered primitive variables in x2-direction.
      //! This function should be called over [js-1,je+1] to get BOTH L/R states over [js,je]

      KOKKOS_INLINE_FUNCTION
      void AddWbPrimFaceX2(TeamMember_t const &member, const int m, const int k, const int j,
           const int il, const int iu, const DvceArray5D<Real> &q,
           ScrArray2D<Real> &ql_jp1, ScrArray2D<Real> &qr_j) {
        int nvar = q.extent_int(1);
        for (int n=0; n<nvar; ++n) {
          par_for_inner(member, il, iu, [&](const int i) {
            AddWbPrimFace(q(m,n,k,j+1,i), q(m,n,k,j,i), ql_jp1(n,i), qr_j(n,i));
          });
        }
        return;
      };

      //----------------------------------------------------------------------------------------
      //! \fn AddWbPrimFaceX3()
      //! \brief Adds background states onto face-centered primitive variables in x3-direction.
      //! This function should be called over [ks-1,ke+1] to get BOTH L/R states over [ks,ke]

      KOKKOS_INLINE_FUNCTION
      void AddWbPrimFaceX3(TeamMember_t const &member, const int m, const int k, const int j,
           const int il, const int iu, const DvceArray5D<Real> &q,
           ScrArray2D<Real> &ql_kp1, ScrArray2D<Real> &qr_k) {
        int nvar = q.extent_int(1);
        for (int n=0; n<nvar; ++n) {
          par_for_inner(member, il, iu, [&](const int i) {
            AddWbPrimFace(q(m,n,k+1,j,i), q(m,n,k,j,i), ql_kp1(n,i), qr_k(n,i));
          });
        }
        return;
      };

      KOKKOS_INLINE_FUNCTION
      void WbLocalPiecewiseLinearX1(TeamMember_t const &member, const EOS_Data &eos, const int m, const int k, const int j,
           const int il, const int iu, const DvceArray5D<Real> &q, const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
           ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
          
        const auto wb_option_ = wb_option;
          
        Real igm1 = 1.0/(eos.gamma-1.0);
        int nvar = q.extent_int(1);
        for (int n=0; n<nvar; ++n) {
          if (n == (IEN) || (n == (IDN) && use_wb_rho)) {
            par_for_inner(member, il, iu, [&](const int i) {
              Real q0_im1, q0_ip1, q0_imh, q0_iph, q0_i;
              getWBerho(n, eos.gamma,
                       q(m,IDN,k,j,i-1),q(m,IDN,k,j,i),q(m,IDN,k,j,i+1),
                       q(m,IEN,k,j,i-1),q(m,IEN,k,j,i),q(m,IEN,k,j,i+1),
                       phicc(m,k,j,i-1),phi(m,k,j,i),phicc(m,k,j,i),phi(m,k,j,i+1),phicc(m,k,j,i+1),
                       q0_im1, q0_imh, q0_i, q0_iph, q0_ip1);
                
              Real q1_im1 = q(m,n,k,j,i-1) - q0_im1;
              Real q1_i = q(m,n,k,j,i) - q0_i;
              Real q1_ip1 = q(m,n,k,j,i+1) - q0_ip1;
                
              PLM(q1_im1, q1_i, q1_ip1, ql(n,i+1), qr(n,i));
                
              ql(n,i+1) += q0_iph;
              qr(n,i) += q0_imh;
                
                // Left/right slopes (properly scaled)
                Real sL = q(m,n,k,j,i)   - q(m,n,k,j,i-1);
                Real sR = q(m,n,k,j,i+1) - q(m,n,k,j,i);
                
                if (ql(n,i+1) < 0.0 || qr(n,i) < 0.0 || sL * sR <= 0.0) {
                    PLM(q(m,n,k,j,i-1), q(m,n,k,j,i), q(m,n,k,j,i+1), ql(n,i+1), qr(n,i));
                }
            });
          } else {
            par_for_inner(member, il, iu, [&](const int i) {
              PLM(q(m,n,k,j,i-1), q(m,n,k,j,i), q(m,n,k,j,i+1), ql(n,i+1), qr(n,i));
            });
          }
        }
        return;
      };
      
      KOKKOS_INLINE_FUNCTION
      void WbLocalPiecewiseLinearX2(TeamMember_t const &member, const EOS_Data &eos, const int m, const int k, const int j,
           const int il, const int iu, const DvceArray5D<Real> &q, const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
           ScrArray2D<Real> &ql_jp1, ScrArray2D<Real> &qr_j) {
        Real igm1 = 1.0/(eos.gamma-1.0);
        int nvar = q.extent_int(1);
        for (int n=0; n<nvar; ++n) {
          if (n == (IEN) || (n == (IDN) && use_wb_rho)) {
            par_for_inner(member, il, iu, [&](const int i) {
              Real q0_jm1, q0_jp1, q0_jmh, q0_jph, q0_j;
              getWBerho(n, eos.gamma,
                       q(m,IDN,k,j-1,i),q(m,IDN,k,j,i),q(m,IDN,k,j+1,i),
                       q(m,IEN,k,j-1,i),q(m,IEN,k,j,i),q(m,IEN,k,j+1,i),
                       phicc(m,k,j-1,i),phi(m,k,j,i),phicc(m,k,j,i),phi(m,k,j+1,i),phicc(m,k,j+1,i),
                       q0_jm1, q0_jmh, q0_j, q0_jph, q0_jp1);
                
              Real q1_jm1 = q(m,n,k,j-1,i) - q0_jm1;
              Real q1_j = q(m,n,k,j,i) - q0_j;
              Real q1_jp1 = q(m,n,k,j+1,i) - q0_jp1;
                
              PLM(q1_jm1, q1_j, q1_jp1, ql_jp1(n,i), qr_j(n,i));
                
              ql_jp1(n,i) += q0_jph;
              qr_j(n,i) += q0_jmh;
                
                // Left/right slopes (properly scaled)
                Real sL = q(m,n,k,j,i)   - q(m,n,k,j-1,i);
                Real sR = q(m,n,k,j+1,i) - q(m,n,k,j,i);
                
                if (ql_jp1(n,i) < 0.0 || qr_j(n,i) < 0.0 || sL * sR <= 0.0) {
                    PLM(q(m,n,k,j-1,i), q(m,n,k,j,i), q(m,n,k,j+1,i), ql_jp1(n,i), qr_j(n,i));
                }
            });
          } else {
            par_for_inner(member, il, iu, [&](const int i) {
              PLM(q(m,n,k,j-1,i), q(m,n,k,j,i), q(m,n,k,j+1,i), ql_jp1(n,i), qr_j(n,i));
            });
          }
        }
        return;
      };
      
      KOKKOS_INLINE_FUNCTION
      void WbLocalPiecewiseLinearX3(TeamMember_t const &member, const EOS_Data &eos, const int m, const int k, const int j,
           const int il, const int iu, const DvceArray5D<Real> &q, const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
           ScrArray2D<Real> &ql_kp1, ScrArray2D<Real> &qr_k) {
        Real igm1 = 1.0/(eos.gamma-1.0);
        int nvar = q.extent_int(1);
        for (int n=0; n<nvar; ++n) {
          if (n == (IEN) || (n == (IDN) && use_wb_rho)) {
            par_for_inner(member, il, iu, [&](const int i) {
              Real q0_km1, q0_kp1, q0_kmh, q0_kph, q0_k;
              getWBerho(n, eos.gamma,
                       q(m,IDN,k-1,j,i),q(m,IDN,k,j,i),q(m,IDN,k+1,j,i),
                       q(m,IEN,k-1,j,i),q(m,IEN,k,j,i),q(m,IEN,k+1,j,i),
                       phicc(m,k-1,j,i),phi(m,k,j,i),phicc(m,k,j,i),phi(m,k+1,j,i),phicc(m,k+1,j,i),
                       q0_km1, q0_kmh, q0_k, q0_kph, q0_kp1);
                
              Real q1_km1 = q(m,n,k-1,j,i) - q0_km1;
              Real q1_k = q(m,n,k,j,i) - q0_k;
              Real q1_kp1 = q(m,n,k+1,j,i) - q0_kp1;
                
              PLM(q1_km1, q1_k, q1_kp1, ql_kp1(n,i), qr_k(n,i));
                
              ql_kp1(n,i) += q0_kph;
              qr_k(n,i) += q0_kmh;
                
                // Left/right slopes (properly scaled)
                Real sL = q(m,n,k,j,i)   - q(m,n,k-1,j,i);
                Real sR = q(m,n,k+1,j,i) - q(m,n,k,j,i);
                
                if (ql_kp1(n,i) < 0.0 || qr_k(n,i) < 0.0 || sL * sR <= 0.0) {
                    PLM(q(m,n,k-1,j,i), q(m,n,k,j,i), q(m,n,k+1,j,i), ql_kp1(n,i), qr_k(n,i));
                }
            });
          } else {
            par_for_inner(member, il, iu, [&](const int i) {
              PLM(q(m,n,k-1,j,i), q(m,n,k,j,i), q(m,n,k+1,j,i), ql_kp1(n,i), qr_k(n,i));
            });
          }
        }
        return;
      };
    
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
    void GridPiecewiseLinearX1(TeamMember_t const &member, const EOS_Data &eos, const int m, const int k, const int j,
         const int il, const int iu, const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf, const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi, const bool hyd,
        ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
      auto &size = pmy_pack->pmb->mb_size;
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int is = indcs.is;
      Real gamma = eos.gamma;
      int nvar = q.extent_int(1);
      for (int n=0; n<nvar; ++n) {
          if ((n == (IEN) || (n == (IDN) && use_wb_rho)) && use_wellbalance_dynamic && use_wb_x1 && hyd) {
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
                
              Real q0_im1, q0_ip1, q0_imh, q0_iph, q0_i;
              getWBerho(n, eos.gamma,
                     q(m,IDN,k,j,i-1),q(m,IDN,k,j,i),q(m,IDN,k,j,i+1),
                     q(m,IEN,k,j,i-1),q(m,IEN,k,j,i),q(m,IEN,k,j,i+1),
                     phicc(m,k,j,i-1),phi(m,k,j,i),phicc(m,k,j,i),phi(m,k,j,i+1),phicc(m,k,j,i+1),
                     q0_im1, q0_imh, q0_i, q0_iph, q0_ip1);
                
              Real q1_im1 = q(m,n,k,j,i-1) - q0_im1;
              Real q1_i = q(m,n,k,j,i) - q0_i;
              Real q1_ip1 = q(m,n,k,j,i+1) - q0_ip1;
                  
              PLM_nonuniform(q1_im1, q1_i, q1_ip1, dxL, dxR, dxLh, dxRh, ql(n,i+1), qr(n,i));
                  
              ql(n,i+1) += q0_iph;
              qr(n,i) += q0_imh;
                
                // Left/right slopes (properly scaled)
                Real sL = (q(m,n,k,j,i)   - q(m,n,k,j,i-1)) / dxL;
                Real sR = (q(m,n,k,j,i+1) - q(m,n,k,j,i)  ) / dxR;
                
                if (ql(n,i+1) < 0.0 || qr(n,i) < 0.0 || sL * sR <= 0.0) {
                    PLM_nonuniform(q(m,n,k,j,i-1), q(m,n,k,j,i), q(m,n,k,j,i+1), dxL, dxR, dxLh, dxRh, ql(n,i+1), qr(n,i));
                }
            });
          } else {
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
      }
      return;
    };
    
    KOKKOS_INLINE_FUNCTION
    void GridPiecewiseLinearX2(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
        const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf,
        ScrArray2D<Real> &ql_jp1, ScrArray2D<Real> &qr_j) {
      Real x_jm1  = xv(m,j-1);
      Real x_jmh  = xf(m,j);
      Real x_j    = xv(m,j);
      Real x_jph  = xf(m,j+1);
      Real x_jp1  = xv(m,j+1);
      Real dxLh = x_j-x_jmh;
      Real dxRh = x_jph-x_j;
      Real dxL = x_j-x_jm1;
      Real dxR = x_jp1-x_j;
      int nvar = q.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        par_for_inner(member, il, iu, [&](const int i) {
          PLM_nonuniform(q(m,n,k,j-1,i), q(m,n,k,j,i), q(m,n,k,j+1,i), dxL, dxR, dxLh, dxRh, ql_jp1(n,i), qr_j(n,i));
        });
      }
      return;
    };
    
    KOKKOS_INLINE_FUNCTION
    void GridPiecewiseLinearX3(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
        const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf,
        ScrArray2D<Real> &ql_kp1, ScrArray2D<Real> &qr_k) {
      Real x_km1  = xv(m,k-1);
      Real x_kmh  = xf(m,k);
      Real x_k    = xv(m,k);
      Real x_kph  = xf(m,k+1);
      Real x_kp1  = xv(m,k+1);
      Real dxLh = x_k-x_kmh;
      Real dxRh = x_kph-x_k;
      Real dxL = x_k-x_km1;
      Real dxR = x_kp1-x_k;
      int nvar = q.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        par_for_inner(member, il, iu, [&](const int i) {
          PLM_nonuniform(q(m,n,k-1,j,i), q(m,n,k,j,i), q(m,n,k+1,j,i), dxL, dxR, dxLh, dxRh, ql_kp1(n,i), qr_k(n,i));
        });
      }
      return;
    };
    
    KOKKOS_INLINE_FUNCTION
    void getWBerho(const int &n, const Real &gamma,
                const Real &rho_im1, const Real &rho_i, const Real &rho_ip1,
                const Real &e_im1, const Real &e_i, const Real &e_ip1,
                const Real &phi_im1, const Real &phi_imh, const Real &phi_i, const Real &phi_iph, const Real &phi_ip1,
                Real &q0_im1, Real &q0_imh, Real &q0_i, Real &q0_iph, Real &q0_ip1) {
      const auto wb_option_ = wb_option;
      int wb_option_num;
      Real igm1 = 1.0/(gamma-1.0);
      Real gigm1 = gamma*igm1;
      Real gm1ig = (gamma-1.0)/gamma;
      Real ig = 1.0/gamma;
        
      switch (wb_option_) {
        case WBOption::isodensity:
          {
            wb_option_num = 0;
          }
          break;
        case WBOption::isothermal:
          {
            wb_option_num = 1;
          }
          break;
        case WBOption::isentropic:
          {
            wb_option_num = 2;
          }
          break;
        case WBOption::adaptive:
          {
            Real dTdivT = fabs((e_ip1/rho_ip1 - e_im1/rho_im1) / (e_i/rho_i));
            Real dsdivs = fabs((e_ip1/pow(rho_ip1,gamma) - e_im1/pow(rho_im1,gamma)) / (e_i/pow(rho_i,gamma)));
            if (dTdivT > dsdivs) {
              wb_option_num = 2;
            } else {
              wb_option_num = 1;
            }
          }
          break;
        default:
          break;
      }
        
        Real factor_im1, factor_i, factor_ip1;
        if (n == (IEN)) {    // reconstructing e
          if (wb_option_num == 0)
              {
                factor_im1 = rho_im1*igm1;
                factor_i = rho_i*igm1;
                factor_ip1 = rho_ip1*igm1;
                q0_i = e_i;
              }
          else if (wb_option_num == 1)
              {
                factor_im1 = rho_im1/e_im1*igm1;
                factor_i = rho_i/e_i*igm1;
                factor_ip1 = rho_ip1/e_ip1*igm1;
                q0_i = log(e_i);
              }
          else
              {
                factor_im1 = rho_im1/pow(e_im1,ig)*ig;
                factor_i = rho_i/pow(e_i,ig)*ig;
                factor_ip1 = rho_ip1/pow(e_ip1,ig)*ig;
                q0_i = pow(e_i,gm1ig);
              }
        } else {    // reconstructing rho
            if (wb_option_num == 0)
              {
                q0_im1 = 0.0;
                q0_ip1 = 0.0;
                q0_imh = 0.0;
                q0_iph = 0.0;
                q0_i = 0.0;
                return;
              }
            else if (wb_option_num == 1)
              {
                factor_im1 = rho_im1/e_im1*igm1;
                factor_i = rho_i/e_i*igm1;
                factor_ip1 = rho_ip1/e_ip1*igm1;
                q0_i = log(rho_i);
              }
            else
              {
                factor_im1 = pow(rho_im1,gamma)/e_im1*ig;
                factor_i = pow(rho_i,gamma)/e_i*ig;
                factor_ip1 = pow(rho_ip1,gamma)/e_ip1*ig;
                q0_i = pow(rho_i,gamma-1.0);
              }
        }
          
        Real dphi_im1 = phi_imh-phi_im1;
        Real dphi_imh = phi_i-phi_imh;
        Real dphi_i = phi_iph-phi_i;
        Real dphi_iph = phi_ip1-phi_iph;
        
        q0_im1 = q0_i + factor_im1 * dphi_im1 + factor_i * dphi_imh;
        q0_ip1 = q0_i - factor_i * dphi_i - factor_ip1 * dphi_iph;
        q0_imh = q0_i + factor_i * dphi_imh;
        q0_iph = q0_i - factor_i * dphi_i;
          
        if (n == (IEN)) {    // reconstructing e
            q0_i = e_i;
            if (wb_option_num == 1)
                {
                  q0_im1 = exp(q0_im1);
                  q0_ip1 = exp(q0_ip1);
                  q0_imh = exp(q0_imh);
                  q0_iph = exp(q0_iph);
                }
            if (wb_option_num == 2)
                {
                  q0_im1 = pow(q0_im1,gigm1);
                  q0_ip1 = pow(q0_ip1,gigm1);
                  q0_imh = pow(q0_imh,gigm1);
                  q0_iph = pow(q0_iph,gigm1);
                }
        } else {    // reconstructing rho
            q0_i = rho_i;
            if (wb_option_num == 1)
                {
                  q0_im1 = exp(q0_im1);
                  q0_ip1 = exp(q0_ip1);
                  q0_imh = exp(q0_imh);
                  q0_iph = exp(q0_iph);
                }
            if (wb_option_num == 2)
                {
                  q0_im1 = pow(q0_im1,igm1);
                  q0_ip1 = pow(q0_ip1,igm1);
                  q0_imh = pow(q0_imh,igm1);
                  q0_iph = pow(q0_iph,igm1);
                }
        }
      
      return;
    };

 private:
  MeshBlockPack* pmy_pack;   // ptr to MeshBlockPack containing this MHD
  // temporary variables used to store face-centered electric fields returned by RS
  DvceArray4D<Real> e1_cc, e2_cc, e3_cc;
};

} // namespace mhd
#endif // MHD_MHD_HPP_
