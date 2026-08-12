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
#include "utils/wb_background.hpp"

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

  // Derived thermodynamic variables (pressure, Gamma_1) for a general EOS. Evaluated once
  // per cell in ConsToPrim and reconstructed to interfaces, so the Riemann solvers never
  // call the EOS. Only allocated when the EOS is general; empty for an ideal gas.
  DvceArray5D<Real> wder;
  // Cached temperature, used as the warm start for the T(d,e) root find. Retaining the
  // previous step's value reduces the inversion to 1-2 Newton iterations. Only allocated
  // when the EOS is general.
  DvceArray4D<Real> wtemp;

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

  // Background gas pressure, evaluated ONCE from the static background (w0wb, w0facewb)
  // rather than per cell per stage. Under a general EOS p(d,e) is a root find, so asking
  // the EOS for the background pressure inside the reconstruction, the flux correction
  // and the coordinate source terms -- roughly twenty calls per cell per stage -- would
  // cost more than the rest of the well-balanced scheme put together. The background does
  // not evolve, so none of that work is ever different. See SetWbBackgroundPressure().
  DvceArray4D<Real> pwb;        // cell-centered background gas pressure
  DvceFaceFld4D<Real> pfacewb;  // face-centered background gas pressure

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
    
  // evaluate the background pressure from the static background state
  void SetWbBackgroundPressure();
  void RemoveWbFlux(const DvceFaceFld4D<Real> &pfacewb, DvceFaceFld5D<Real> &flx);
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
            getWBq0(eos, (n == (IEN)) ? WBVar::wb_eint : WBVar::wb_dens,
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
            getWBq0(eos, (n == (IEN)) ? WBVar::wb_eint : WBVar::wb_dens,
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
            getWBq0(eos, (n == (IEN)) ? WBVar::wb_eint : WBVar::wb_dens,
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
    
    //------------------------------------------------------------------------------------
    //! \fn WbStaticPiecewiseLinearDerX1()
    //! \brief static well-balanced reconstruction of the derived thermodynamic variables.
    //!
    //! The static scheme reconstructs the perturbation of the primitives about a
    //! problem-generator-supplied background and adds the background back at the faces.
    //! Under a general EOS the pressure the Riemann solver sees comes from wder, not from
    //! the reconstructed primitives, so it needs the same treatment or the stratified
    //! part
    //! of the pressure goes through plain PLM and the deviation reconstruction buys
    //! nothing. The background pressure is not written as (gamma-1)*e_bg -- it depends on
    //! the background density too -- but neither is it asked of the EOS here: the
    //! background is static, so its pressure is evaluated once into pwb/pfwb at startup
    //! by SetWbBackgroundPressure(). Under a general EOS each of these five reads would
    //! otherwise be a temperature root find, in every cell of every stage.

    KOKKOS_INLINE_FUNCTION
    void WbStaticPiecewiseLinearDerX1(TeamMember_t const &member,
         const int m, const int k, const int j, const int il, const int iu,
         const DvceArray4D<Real> &pwb, const DvceArray4D<Real> &pfwb,
         const DvceArray5D<Real> &qd, ScrArray2D<Real> &dl, ScrArray2D<Real> &dr) {
      int nvar = qd.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        if (n == (IDPR)) {
          par_for_inner(member, il, iu, [&](const int i) {
            Real q0_im1 = pwb(m,k,j,i-1);
            Real q0_i   = pwb(m,k,j,i);
            Real q0_ip1 = pwb(m,k,j,i+1);

            PLM(qd(m,n,k,j,i-1) - q0_im1, qd(m,n,k,j,i) - q0_i,
                qd(m,n,k,j,i+1) - q0_ip1, dl(n,i+1), dr(n,i));

            dl(n,i+1) += pfwb(m,k,j,i+1);
            dr(n,i) += pfwb(m,k,j,i);
          });
        } else {
          par_for_inner(member, il, iu, [&](const int i) {
            PLM(qd(m,n,k,j,i-1), qd(m,n,k,j,i), qd(m,n,k,j,i+1), dl(n,i+1), dr(n,i));
          });
        }
      }
      return;
    };

    //------------------------------------------------------------------------------------
    //! \fn WbStaticPiecewiseLinearDerX2()
    //! \brief static well-balanced reconstruction of the derived thermodynamic variables.
    //!
    //! The static scheme reconstructs the perturbation of the primitives about a
    //! problem-generator-supplied background and adds the background back at the faces.
    //! Under a general EOS the pressure the Riemann solver sees comes from wder, not from
    //! the reconstructed primitives, so it needs the same treatment or the stratified
    //! part
    //! of the pressure goes through plain PLM and the deviation reconstruction buys
    //! nothing. The background pressure is not written as (gamma-1)*e_bg -- it depends on
    //! the background density too -- but neither is it asked of the EOS here: the
    //! background is static, so its pressure is evaluated once into pwb/pfwb at startup
    //! by SetWbBackgroundPressure(). Under a general EOS each of these five reads would
    //! otherwise be a temperature root find, in every cell of every stage.

    KOKKOS_INLINE_FUNCTION
    void WbStaticPiecewiseLinearDerX2(TeamMember_t const &member,
         const int m, const int k, const int j, const int il, const int iu,
         const DvceArray4D<Real> &pwb, const DvceArray4D<Real> &pfwb,
         const DvceArray5D<Real> &qd, ScrArray2D<Real> &dl, ScrArray2D<Real> &dr) {
      int nvar = qd.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        if (n == (IDPR)) {
          par_for_inner(member, il, iu, [&](const int i) {
            Real q0_jm1 = pwb(m,k,j-1,i);
            Real q0_j   = pwb(m,k,j,i);
            Real q0_jp1 = pwb(m,k,j+1,i);

            PLM(qd(m,n,k,j-1,i) - q0_jm1, qd(m,n,k,j,i) - q0_j,
                qd(m,n,k,j+1,i) - q0_jp1, dl(n,i), dr(n,i));

            dl(n,i) += pfwb(m,k,j+1,i);
            dr(n,i) += pfwb(m,k,j,i);
          });
        } else {
          par_for_inner(member, il, iu, [&](const int i) {
            PLM(qd(m,n,k,j-1,i), qd(m,n,k,j,i), qd(m,n,k,j+1,i), dl(n,i), dr(n,i));
          });
        }
      }
      return;
    };

    //------------------------------------------------------------------------------------
    //! \fn WbStaticPiecewiseLinearDerX3()
    //! \brief static well-balanced reconstruction of the derived thermodynamic variables.
    //!
    //! The static scheme reconstructs the perturbation of the primitives about a
    //! problem-generator-supplied background and adds the background back at the faces.
    //! Under a general EOS the pressure the Riemann solver sees comes from wder, not from
    //! the reconstructed primitives, so it needs the same treatment or the stratified
    //! part
    //! of the pressure goes through plain PLM and the deviation reconstruction buys
    //! nothing. The background pressure is not written as (gamma-1)*e_bg -- it depends on
    //! the background density too -- but neither is it asked of the EOS here: the
    //! background is static, so its pressure is evaluated once into pwb/pfwb at startup
    //! by SetWbBackgroundPressure(). Under a general EOS each of these five reads would
    //! otherwise be a temperature root find, in every cell of every stage.

    KOKKOS_INLINE_FUNCTION
    void WbStaticPiecewiseLinearDerX3(TeamMember_t const &member,
         const int m, const int k, const int j, const int il, const int iu,
         const DvceArray4D<Real> &pwb, const DvceArray4D<Real> &pfwb,
         const DvceArray5D<Real> &qd, ScrArray2D<Real> &dl, ScrArray2D<Real> &dr) {
      int nvar = qd.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        if (n == (IDPR)) {
          par_for_inner(member, il, iu, [&](const int i) {
            Real q0_km1 = pwb(m,k-1,j,i);
            Real q0_k   = pwb(m,k,j,i);
            Real q0_kp1 = pwb(m,k+1,j,i);

            PLM(qd(m,n,k-1,j,i) - q0_km1, qd(m,n,k,j,i) - q0_k,
                qd(m,n,k+1,j,i) - q0_kp1, dl(n,i), dr(n,i));

            dl(n,i) += pfwb(m,k+1,j,i);
            dr(n,i) += pfwb(m,k,j,i);
          });
        } else {
          par_for_inner(member, il, iu, [&](const int i) {
            PLM(qd(m,n,k-1,j,i), qd(m,n,k,j,i), qd(m,n,k+1,j,i), dl(n,i), dr(n,i));
          });
        }
      }
      return;
    };

    //------------------------------------------------------------------------------------
    //! \fn WbPiecewiseLinearDerX1/X2/X3()
    //! \brief well-balanced reconstruction of the DERIVED thermodynamic variables.
    //!
    //! Under a general EOS the Riemann solvers do not recompute pressure from the
    //! reconstructed (d,e); they read it from wder, which is reconstructed to the
    //! interfaces independently. Reconstructing d and e in deviation form while leaving
    //! pressure to plain PLM would leave the whole hydrostatic gradient in the pressure
    //! the solver actually sees, and the scheme would not be balanced at all. So the
    //! pressure channel gets the same deviation treatment, against the pressure of the
    //! same background. Gamma_1 is left alone: it is a smooth O(1) quantity that is not
    //! stratified over many scale heights, so plain PLM is appropriate.
    //!
    //! These are only ever called for a general EOS; an ideal gas allocates no wder.

    KOKKOS_INLINE_FUNCTION
    void WbPiecewiseLinearDerX1(TeamMember_t const &member, const EOS_Data &eos,
         const int m, const int k, const int j, const int il, const int iu,
         const DvceArray5D<Real> &q, const DvceArray5D<Real> &qd,
         const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
         ScrArray2D<Real> &dl, ScrArray2D<Real> &dr) {
      int nvar = qd.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        if (n == (IDPR)) {
          par_for_inner(member, il, iu, [&](const int i) {
            Real q0_im1, q0_imh, q0_i, q0_iph, q0_ip1;
            getWBq0(eos, WBVar::wb_pres,
                    q(m,IDN,k,j,i-1),q(m,IDN,k,j,i),q(m,IDN,k,j,i+1),
                    q(m,IEN,k,j,i-1),q(m,IEN,k,j,i),q(m,IEN,k,j,i+1),
                    phicc(m,k,j,i-1),phi(m,k,j,i),phicc(m,k,j,i),
                    phi(m,k,j,i+1),phicc(m,k,j,i+1),
                    q0_im1, q0_imh, q0_i, q0_iph, q0_ip1);

            Real q1_im1 = qd(m,n,k,j,i-1) - q0_im1;
            Real q1_i   = qd(m,n,k,j,i)   - q0_i;
            Real q1_ip1 = qd(m,n,k,j,i+1) - q0_ip1;

            PLM(q1_im1, q1_i, q1_ip1, dl(n,i+1), dr(n,i));

            dl(n,i+1) += q0_iph;
            dr(n,i)   += q0_imh;

            // fall back to plain reconstruction where the deviation form misbehaves,
            // matching what the primitive channels do
            Real sL = qd(m,n,k,j,i)   - qd(m,n,k,j,i-1);
            Real sR = qd(m,n,k,j,i+1) - qd(m,n,k,j,i);
            if (dl(n,i+1) < 0.0 || dr(n,i) < 0.0 || sL * sR <= 0.0) {
              PLM(qd(m,n,k,j,i-1), qd(m,n,k,j,i), qd(m,n,k,j,i+1), dl(n,i+1), dr(n,i));
            }
          });
        } else {
          par_for_inner(member, il, iu, [&](const int i) {
            PLM(qd(m,n,k,j,i-1), qd(m,n,k,j,i), qd(m,n,k,j,i+1), dl(n,i+1), dr(n,i));
          });
        }
      }
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void WbPiecewiseLinearDerX2(TeamMember_t const &member, const EOS_Data &eos,
         const int m, const int k, const int j, const int il, const int iu,
         const DvceArray5D<Real> &q, const DvceArray5D<Real> &qd,
         const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
         ScrArray2D<Real> &dl_jp1, ScrArray2D<Real> &dr_j) {
      int nvar = qd.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        if (n == (IDPR)) {
          par_for_inner(member, il, iu, [&](const int i) {
            Real q0_jm1, q0_jmh, q0_j, q0_jph, q0_jp1;
            getWBq0(eos, WBVar::wb_pres,
                    q(m,IDN,k,j-1,i),q(m,IDN,k,j,i),q(m,IDN,k,j+1,i),
                    q(m,IEN,k,j-1,i),q(m,IEN,k,j,i),q(m,IEN,k,j+1,i),
                    phicc(m,k,j-1,i),phi(m,k,j,i),phicc(m,k,j,i),
                    phi(m,k,j+1,i),phicc(m,k,j+1,i),
                    q0_jm1, q0_jmh, q0_j, q0_jph, q0_jp1);

            Real q1_jm1 = qd(m,n,k,j-1,i) - q0_jm1;
            Real q1_j   = qd(m,n,k,j,i)   - q0_j;
            Real q1_jp1 = qd(m,n,k,j+1,i) - q0_jp1;

            PLM(q1_jm1, q1_j, q1_jp1, dl_jp1(n,i), dr_j(n,i));

            dl_jp1(n,i) += q0_jph;
            dr_j(n,i)   += q0_jmh;

            Real sL = qd(m,n,k,j,i)   - qd(m,n,k,j-1,i);
            Real sR = qd(m,n,k,j+1,i) - qd(m,n,k,j,i);
            if (dl_jp1(n,i) < 0.0 || dr_j(n,i) < 0.0 || sL * sR <= 0.0) {
              PLM(qd(m,n,k,j-1,i), qd(m,n,k,j,i), qd(m,n,k,j+1,i),
                  dl_jp1(n,i), dr_j(n,i));
            }
          });
        } else {
          par_for_inner(member, il, iu, [&](const int i) {
            PLM(qd(m,n,k,j-1,i), qd(m,n,k,j,i), qd(m,n,k,j+1,i), dl_jp1(n,i), dr_j(n,i));
          });
        }
      }
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void WbPiecewiseLinearDerX3(TeamMember_t const &member, const EOS_Data &eos,
         const int m, const int k, const int j, const int il, const int iu,
         const DvceArray5D<Real> &q, const DvceArray5D<Real> &qd,
         const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
         ScrArray2D<Real> &dl_kp1, ScrArray2D<Real> &dr_k) {
      int nvar = qd.extent_int(1);
      for (int n=0; n<nvar; ++n) {
        if (n == (IDPR)) {
          par_for_inner(member, il, iu, [&](const int i) {
            Real q0_km1, q0_kmh, q0_k, q0_kph, q0_kp1;
            getWBq0(eos, WBVar::wb_pres,
                    q(m,IDN,k-1,j,i),q(m,IDN,k,j,i),q(m,IDN,k+1,j,i),
                    q(m,IEN,k-1,j,i),q(m,IEN,k,j,i),q(m,IEN,k+1,j,i),
                    phicc(m,k-1,j,i),phi(m,k,j,i),phicc(m,k,j,i),
                    phi(m,k+1,j,i),phicc(m,k+1,j,i),
                    q0_km1, q0_kmh, q0_k, q0_kph, q0_kp1);

            Real q1_km1 = qd(m,n,k-1,j,i) - q0_km1;
            Real q1_k   = qd(m,n,k,j,i)   - q0_k;
            Real q1_kp1 = qd(m,n,k+1,j,i) - q0_kp1;

            PLM(q1_km1, q1_k, q1_kp1, dl_kp1(n,i), dr_k(n,i));

            dl_kp1(n,i) += q0_kph;
            dr_k(n,i)   += q0_kmh;

            Real sL = qd(m,n,k,j,i)   - qd(m,n,k-1,j,i);
            Real sR = qd(m,n,k+1,j,i) - qd(m,n,k,j,i);
            if (dl_kp1(n,i) < 0.0 || dr_k(n,i) < 0.0 || sL * sR <= 0.0) {
              PLM(qd(m,n,k-1,j,i), qd(m,n,k,j,i), qd(m,n,k+1,j,i),
                  dl_kp1(n,i), dr_k(n,i));
            }
          });
        } else {
          par_for_inner(member, il, iu, [&](const int i) {
            PLM(qd(m,n,k-1,j,i), qd(m,n,k,j,i), qd(m,n,k+1,j,i), dl_kp1(n,i), dr_k(n,i));
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
         const int il, const int iu, const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf, const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
        ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
      auto &size = pmy_pack->pmb->mb_size;
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int is = indcs.is;
      Real gamma = eos.gamma;
      int nvar = q.extent_int(1);
      for (int n=0; n<nvar; ++n) {
          if ((n == (IEN) || (n == (IDN) && use_wb_rho)) && use_wellbalance_dynamic && use_wb_x1) {
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
              getWBq0(eos, (n == (IEN)) ? WBVar::wb_eint : WBVar::wb_dens,
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

    //------------------------------------------------------------------------------------
    //! \fn GridPiecewiseLinearDerX1()
    //! \brief non-uniform-grid counterpart of WbPiecewiseLinearDerX1, for the derived
    //! thermodynamic variables on a spherical-polar mesh. GridPiecewiseLinearX2/X3 are
    //! already generic in the array they reconstruct and do no well-balancing, so wder
    //! can
    //! be passed to them directly; only x1 needs its own version.

    KOKKOS_INLINE_FUNCTION
    void GridPiecewiseLinearDerX1(TeamMember_t const &member, const EOS_Data &eos,
         const int m, const int k, const int j, const int il, const int iu,
         const DvceArray5D<Real> &q, const DvceArray5D<Real> &qd,
         const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf,
         const DvceArray4D<Real> &phicc, const DvceArray4D<Real> &phi,
         ScrArray2D<Real> &dl, ScrArray2D<Real> &dr) {
      int nvar = qd.extent_int(1);
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

          if (n == (IDPR) && use_wellbalance_dynamic && use_wb_x1) {
            Real q0_im1, q0_imh, q0_i, q0_iph, q0_ip1;
            getWBq0(eos, WBVar::wb_pres,
                    q(m,IDN,k,j,i-1),q(m,IDN,k,j,i),q(m,IDN,k,j,i+1),
                    q(m,IEN,k,j,i-1),q(m,IEN,k,j,i),q(m,IEN,k,j,i+1),
                    phicc(m,k,j,i-1),phi(m,k,j,i),phicc(m,k,j,i),
                    phi(m,k,j,i+1),phicc(m,k,j,i+1),
                    q0_im1, q0_imh, q0_i, q0_iph, q0_ip1);

            Real q1_im1 = qd(m,n,k,j,i-1) - q0_im1;
            Real q1_i   = qd(m,n,k,j,i)   - q0_i;
            Real q1_ip1 = qd(m,n,k,j,i+1) - q0_ip1;

            PLM_nonuniform(q1_im1, q1_i, q1_ip1, dxL, dxR, dxLh, dxRh,
                           dl(n,i+1), dr(n,i));

            dl(n,i+1) += q0_iph;
            dr(n,i)   += q0_imh;

            Real sL = (qd(m,n,k,j,i)   - qd(m,n,k,j,i-1)) / dxL;
            Real sR = (qd(m,n,k,j,i+1) - qd(m,n,k,j,i)  ) / dxR;
            if (dl(n,i+1) < 0.0 || dr(n,i) < 0.0 || sL * sR <= 0.0) {
              PLM_nonuniform(qd(m,n,k,j,i-1), qd(m,n,k,j,i), qd(m,n,k,j,i+1),
                             dxL, dxR, dxLh, dxRh, dl(n,i+1), dr(n,i));
            }
          } else {
            PLM_nonuniform(qd(m,n,k,j,i-1), qd(m,n,k,j,i), qd(m,n,k,j,i+1),
                           dxL, dxR, dxLh, dxRh, dl(n,i+1), dr(n,i));
          }
        });
      }
      return;
    };

    //------------------------------------------------------------------------------------
    //! \fn getWBq0()
    //! \brief background values of ONE reconstruction channel at the five stencil points
    //! {i-1, i-1/2, i, i+1/2, i+1}.  Dispatches to the ideal-gas closed forms in
    //! getWBerho() below, or, for a general EOS, to the integrated background in
    //! utils/wb_background.hpp.  Every well-balanced reconstruction and every
    //! gravitational source term must go through here, so that the background subtracted
    //! before reconstruction and the background differenced by the source term are the
    //! same object -- that consistency, not the background's accuracy, is what makes the
    //! scheme well balanced.

    KOKKOS_INLINE_FUNCTION
    void getWBq0(const EOS_Data &eos, const int var,
                const Real &rho_im1, const Real &rho_i, const Real &rho_ip1,
                const Real &e_im1, const Real &e_i, const Real &e_ip1,
                const Real &phi_im1, const Real &phi_imh, const Real &phi_i,
                const Real &phi_iph, const Real &phi_ip1,
                Real &q0_im1, Real &q0_imh, Real &q0_i, Real &q0_iph, Real &q0_ip1) {
      if (eos.IsGeneral()) {
        WBState s_im1, s_imh, s_i, s_iph, s_ip1;
        WBBackgroundStencil(eos, wb_option, rho_im1, rho_i, rho_ip1, e_im1, e_i, e_ip1,
                            phi_im1, phi_imh, phi_i, phi_iph, phi_ip1,
                            s_im1, s_imh, s_i, s_iph, s_ip1);
        if (var == WBVar::wb_dens) {
          q0_im1 = s_im1.d; q0_imh = s_imh.d; q0_i = s_i.d;
          q0_iph = s_iph.d; q0_ip1 = s_ip1.d;
        } else if (var == WBVar::wb_eint) {
          q0_im1 = s_im1.e; q0_imh = s_imh.e; q0_i = s_i.e;
          q0_iph = s_iph.e; q0_ip1 = s_ip1.e;
        } else {
          q0_im1 = s_im1.p; q0_imh = s_imh.p; q0_i = s_i.p;
          q0_iph = s_iph.p; q0_ip1 = s_ip1.p;
        }
      } else {
        // ideal gas: unchanged closed forms.  The pressure channel is never reconstructed
        // for an ideal gas (nder is 0), but p = (gamma-1)e keeps this total anyway.
        int n = (var == WBVar::wb_dens) ? IDN : IEN;
        getWBerho(n, eos.gamma, rho_im1, rho_i, rho_ip1, e_im1, e_i, e_ip1,
                  phi_im1, phi_imh, phi_i, phi_iph, phi_ip1,
                  q0_im1, q0_imh, q0_i, q0_iph, q0_ip1);
        if (var == WBVar::wb_pres) {
          Real gm1 = eos.gamma - 1.0;
          q0_im1 *= gm1; q0_imh *= gm1; q0_i *= gm1; q0_iph *= gm1; q0_ip1 *= gm1;
        }
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
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Hydro
};

} // namespace hydro
#endif // HYDRO_HYDRO_HPP_
