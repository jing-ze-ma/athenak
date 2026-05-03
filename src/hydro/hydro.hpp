#ifndef HYDRO_HYDRO_HPP_
#define HYDRO_HYDRO_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro.hpp
//  \brief definitions for Hydro class

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
class Conduction;
class SourceTerms;
class OrbitalAdvectionCC;
class ShearingBoxCC;
class Driver;

// constants that enumerate Hydro Riemann Solver options
enum class Hydro_RSolver {advect, llf, hlle, hllc, lhllc, hllclm, ausmpup, roe,    // non-relativistic
                          llf_sr, hlle_sr, hllc_sr,        // SR
                          llf_gr, hlle_gr};                // GR
// constants that enumerate local well-balanced scheme options
enum WB_Option {isodensity, isothermal, isentropic};

//----------------------------------------------------------------------------------------
//! \struct HydroTaskIDs
//  \brief container to hold TaskIDs of all hydro tasks

struct HydroTaskIDs {
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
  TaskID bcs;
  TaskID prol;
  TaskID c2p;
  TaskID newdt;
  TaskID csend;
  TaskID crecv;
};

namespace hydro {

//----------------------------------------------------------------------------------------
//! \class Hydro

class Hydro {
 public:
  Hydro(MeshBlockPack *ppack, ParameterInput *pin);
  ~Hydro();

  // data
  ReconstructionMethod recon_method;
  Hydro_RSolver rsolver_method;
  EquationOfState *peos;  // chosen EOS

  int nhydro;             // number of hydro variables (5/4 for ideal/isothermal EOS)
  int nscalars;           // number of passive scalars
  DvceArray5D<Real> u0;   // conserved variables
  DvceArray5D<Real> w0;   // primitive variables

  DvceArray5D<Real> coarse_u0;  // conserved variables on 2x coarser grid (for SMR/AMR)
  DvceArray5D<Real> coarse_w0;  // primitive variables on 2x coarser grid (for SMR/AMR)

  // Boundary communication buffers and functions for u
  MeshBoundaryValuesCC *pbval_u;

  // Orbital advection and shearing box BCs
  OrbitalAdvectionCC *porb_u = nullptr;
  ShearingBoxCC *psbox_u = nullptr;

  // Object(s) for extra physics (viscosity, thermal conduction, srcterms)
  Viscosity *pvisc = nullptr;
  Conduction *pcond = nullptr;
  SourceTerms *psrc = nullptr;

  // following only used for time-evolving flow
  DvceArray5D<Real> u1;       // conserved variables at intermediate step
  DvceFaceFld5D<Real> uflx;   // fluxes of conserved quantities on cell faces
  Real dtnew;

  // following used for FOFC
  DvceArray4D<bool> fofc;  // flag for each cell to indicate if FOFC is needed
  bool use_fofc = false;   // flag to enable FOFC
  DvceArray5D<Real> utest;  // scratch array for FOFC
    
  // following only used for including time-independent gravity in the conserved energy equation
  bool use_etotgrav = false;   // flag to enable etotgrav
  DvceFaceFld4D<Real> phi0;     // face-centered gravitational potential energy
  DvceArray4D<Real> phicc0;     // cell-centered gravitational potential energy
    
  // following used for well-balanced scheme
  bool use_wellbalance = false;    // flag to enable wellbalance
  bool use_wellbalance_reconst_perturb = false;    // flag to enable reconstructing perturbed primitive variables (less robust against large deviations)
  bool use_wellbalance_local = false;    // flag to enable local well-balanced method by Kappeli & Mishra 2014, 2016
  bool use_wb_x1 = false;    // flag for directions
  bool use_wb_x2 = false;    // flag for directions
  bool use_wb_x3 = false;    // flag for directions
  bool use_wb_rho = false;   // flag to enable local well-balanced method also in reconstructing rho
  WB_Option wb_option;
  DvceArray5D<Real> u0wb;   // background conserved variables
  DvceArray5D<Real> w0wb;   // background primitive variables
  DvceFaceFld5D<Real> w0facewb;   // face-centered background primitive variables
    
  bool use_reconst_logp_x1 = false;
  bool use_reconst_logp_x2 = false;
  bool use_reconst_logp_x3 = false;

  // container to hold names of TaskIDs
  HydroTaskIDs id;

  // functions...
  void AssembleHydroTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  // ...in "before_stagen_tl" list
  TaskStatus InitRecv(Driver *d, int stage);
  // ...in "stagen_tl" list
  TaskStatus CopyCons(Driver *d, int stage);
  TaskStatus Fluxes(Driver *d, int stage);
  TaskStatus SendFlux(Driver *d, int stage);
  TaskStatus RecvFlux(Driver *d, int stage);
  TaskStatus RKUpdate(Driver *d, int stage);
  TaskStatus HydroSrcTerms(Driver *d, int stage);
  TaskStatus SendU_OA(Driver *d, int stage);
  TaskStatus RecvU_OA(Driver *d, int stage);
  TaskStatus RestrictU(Driver *d, int stage);
  TaskStatus SendU(Driver *d, int stage);
  TaskStatus RecvU(Driver *d, int stage);
  TaskStatus SendU_Shr(Driver *d, int stage);
  TaskStatus RecvU_Shr(Driver *d, int stage);
  TaskStatus ApplyPhysicalBCs(Driver* pdrive, int stage);
  TaskStatus Prolongate(Driver* pdrive, int stage);
  TaskStatus ConToPrim(Driver *d, int stage);
  TaskStatus NewTimeStep(Driver *d, int stage);
  // ...in "after_stagen_tl" list
  TaskStatus ClearSend(Driver *d, int stage);
  TaskStatus ClearRecv(Driver *d, int stage);  // also in Driver::Initialize

  // CalculateFluxes function templated over Riemann Solvers
  template <Hydro_RSolver T>
  void CalculateFluxes(Driver *d, int stage);

  // first-order flux correction
  void FOFC(Driver *d, int stage);
    
  void AddGravFlux(const DvceFaceFld4D<Real> &phi0, DvceFaceFld5D<Real> &flx);
  void AddGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons, const int il, const int iu, const int jl, const int ju, const int kl, const int ku);
  void RemoveGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons, const int il, const int iu, const int jl, const int ju, const int kl, const int ku);
    
  void RemoveWbFlux(const DvceFaceFld5D<Real> &w0facewb, DvceFaceFld5D<Real> &flx);
  void AddWbVar(const DvceArray5D<Real> &varwb, DvceArray5D<Real> &var);
  void RemoveWbVar(const DvceArray5D<Real> &varwb, DvceArray5D<Real> &var);
////  KOKKOS_INLINE_FUNCTION
//  void AddWbPrimFace(const Real &qlwb_ip1, const Real &qrwb_i, Real &ql_ip1, Real &qr_i);
////  KOKKOS_INLINE_FUNCTION
//  void AddWbPrimFaceX1(TeamMember_t const &member, const int m, const int k, const int j,
//         const int il, const int iu, const DvceArray5D<Real> &q,
//                       ScrArray2D<Real> &ql, ScrArray2D<Real> &qr);
////  KOKKOS_INLINE_FUNCTION
//  void AddWbPrimFaceX2(TeamMember_t const &member, const int m, const int k, const int j,
//         const int il, const int iu, const DvceArray5D<Real> &q,
//                         ScrArray2D<Real> &ql_jp1, ScrArray2D<Real> &qr_j);
////  KOKKOS_INLINE_FUNCTION
//  void AddWbPrimFaceX3(TeamMember_t const &member, const int m, const int k, const int j,
//         const int il, const int iu, const DvceArray5D<Real> &q,
//                       ScrArray2D<Real> &ql_kp1, ScrArray2D<Real> &qr_k);
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
    }
    
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

      // Left/right slopes (properly scaled)
      Real sL = (q_i   - q_im1) / dxL;
      Real sR = (q_ip1 - q_i  ) / dxR;

      // Monotonized central slope (harmonic mean limiter)
      Real slope = (2.0 * sL * sR) / (sL + sR);  // harmonic mean
      slope = (sL * sR > 0.0) ? slope : 0.0;

      // Reconstruct to faces
      ql_ip1 = q_i + slope * dxRh;
      qr_i   = q_i - slope * dxLh;
      return;
    };
    
    KOKKOS_INLINE_FUNCTION
    void GridPiecewiseLinearX1(TeamMember_t const &member, const EOS_Data &eos, const int m, const int k, const int j,
         const int il, const int iu, const DvceArray5D<Real> &q, const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
        ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
      auto &size = pmy_pack->pmb->mb_size;
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int is = indcs.is;
      Real gamma = eos.gamma;
      int nvar = q.extent_int(1);
      for (int n=0; n<nvar; ++n) {
          if ((n == (IEN) || (n == (IDN) && use_wb_rho)) && use_wellbalance_local && use_wb_x1) {
            par_for_inner(member, il, iu, [&](const int i) {
                Real x1v_im1  = CellCenterX(i-1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_imh  = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_i  = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_iph  = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_ip1  = CellCenterX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_0 = pmy_pack->pmesh->mesh_size.x1min;
                Real x1v_1 = pmy_pack->pmesh->mesh_size.x1max;
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_im1);
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_imh);
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_i);
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_iph);
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_ip1);
                Real dxLh = x1v_i-x1v_imh;
                Real dxRh = x1v_iph-x1v_i;
                Real dxL = x1v_i-x1v_im1;
                Real dxR = x1v_ip1-x1v_i;
                
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
                Real x1v_im1  = CellCenterX(i-1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_imh  = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_i  = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_iph  = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_ip1  = CellCenterX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
                Real x1v_0 = pmy_pack->pmesh->mesh_size.x1min;
                Real x1v_1 = pmy_pack->pmesh->mesh_size.x1max;
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_im1);
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_imh);
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_i);
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_iph);
                pmy_pack->pcoord->StretchR(x1v_0,x1v_1,x1v_ip1);
                Real dxLh = x1v_i-x1v_imh;
                Real dxRh = x1v_iph-x1v_i;
                Real dxL = x1v_i-x1v_im1;
                Real dxR = x1v_ip1-x1v_i;
                
                Real qm1 = q(m,n,k,j,i-1);
                Real q0 = q(m,n,k,j,i);
                Real qp1 = q(m,n,k,j,i+1);
                Real qlp1, qr0;
                if (n == (IEN) && use_reconst_logp_x1) {
                  qm1 = log(qm1);
                  q0 = log(q0);
                  qp1 = log(qp1);
                }
              PLM_nonuniform(qm1, q0, qp1, dxL, dxR, dxLh, dxRh, qlp1, qr0);
                ql(n,i+1) = qlp1;
                qr(n,i) = qr0;
                if (n == (IEN) && use_reconst_logp_x1) {
                  ql(n,i+1) = exp(ql(n,i+1));
                  qr(n,i) = exp(qr(n,i));
                }
            });
          }
      }
      return;
    }
    
    KOKKOS_INLINE_FUNCTION
    void getWBerho(const int &n, const Real &gamma,
                const Real &rho_im1, const Real &rho_i, const Real &rho_ip1,
                const Real &e_im1, const Real &e_i, const Real &e_ip1,
                const Real &phi_im1, const Real &phi_imh, const Real &phi_i, const Real &phi_iph, const Real &phi_ip1,
                Real &q0_im1, Real &q0_imh, Real &q0_i, Real &q0_iph, Real &q0_ip1) {
      const auto wb_option_ = wb_option;
      Real igm1 = 1.0/(gamma-1.0);
      Real gigm1 = gamma*igm1;
      Real gm1ig = (gamma-1.0)/gamma;
      Real ig = 1.0/gamma;
        
        Real factor_im1, factor_i, factor_ip1;
        if (n == (IEN)) {    // reconstructing e
          switch (wb_option_) {
            case WB_Option::isodensity:
              {
                factor_im1 = rho_im1*igm1;
                factor_i = rho_i*igm1;
                factor_ip1 = rho_ip1*igm1;
                q0_i = e_i;
              }
              break;
            case WB_Option::isothermal:
              {
                factor_im1 = rho_im1/e_im1*igm1;
                factor_i = rho_i/e_i*igm1;
                factor_ip1 = rho_ip1/e_ip1*igm1;
                q0_i = log(e_i);
              }
              break;
            case WB_Option::isentropic:
              {
                factor_im1 = rho_im1/pow(e_im1,ig)*ig;
                factor_i = rho_i/pow(e_i,ig)*ig;
                factor_ip1 = rho_ip1/pow(e_ip1,ig)*ig;
                q0_i = pow(e_i,gm1ig);
              }
              break;
            default:
              break;
          }
        } else {    // reconstructing rho
          switch (wb_option_) {
            case WB_Option::isodensity:
              {
                q0_im1 = 0.0;
                q0_ip1 = 0.0;
                q0_imh = 0.0;
                q0_iph = 0.0;
                q0_i = 0.0;
                return;
              }
              break;
            case WB_Option::isothermal:
              {
                factor_im1 = rho_im1/e_im1*igm1;
                factor_i = rho_i/e_i*igm1;
                factor_ip1 = rho_ip1/e_ip1*igm1;
                q0_i = log(rho_i);
              }
              break;
            case WB_Option::isentropic:
              {
                factor_im1 = pow(rho_im1,gamma)/e_im1*ig;
                factor_i = pow(rho_i,gamma)/e_i*ig;
                factor_ip1 = pow(rho_ip1,gamma)/e_ip1*ig;
                q0_i = pow(rho_i,gamma-1.0);
              }
              break;
            default:
              break;
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
            switch (wb_option_) {
              case WB_Option::isothermal:
                {
                  q0_im1 = exp(q0_im1);
                  q0_ip1 = exp(q0_ip1);
                  q0_imh = exp(q0_imh);
                  q0_iph = exp(q0_iph);
                }
                break;
              case WB_Option::isentropic:
                {
                  q0_im1 = pow(q0_im1,gigm1);
                  q0_ip1 = pow(q0_ip1,gigm1);
                  q0_imh = pow(q0_imh,gigm1);
                  q0_iph = pow(q0_iph,gigm1);
                }
                break;
              default:
                break;
            }
        } else {    // reconstructing rho
            q0_i = rho_i;
            switch (wb_option_) {
              case WB_Option::isothermal:
                {
                  q0_im1 = exp(q0_im1);
                  q0_ip1 = exp(q0_ip1);
                  q0_imh = exp(q0_imh);
                  q0_iph = exp(q0_iph);
                }
                break;
              case WB_Option::isentropic:
                {
                  q0_im1 = pow(q0_im1,igm1);
                  q0_ip1 = pow(q0_ip1,igm1);
                  q0_imh = pow(q0_imh,igm1);
                  q0_iph = pow(q0_iph,igm1);
                }
                break;
              default:
                break;
            }
        }
      
      return;
    };

 private:
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Hydro
};

} // namespace hydro
#endif // HYDRO_HYDRO_HPP_
