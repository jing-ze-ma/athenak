#ifndef DIFFUSION_RESISTIVITY_HPP_
#define DIFFUSION_RESISTIVITY_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity.hpp
//  \brief Contains data and functions that implement various non-ideal MHD (resistive)
//  processes: Ohmic diffusion and ambipolar diffusion (<mhd>/ambipolar). The Hall effect
//  is not implemented.

#include <map>
#include <memory>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/meshblock.hpp"

struct ResistivityTaskIDs {
  TaskID res_totstage;
  TaskID res_coeff;
  TaskID mhd_irecv;
  TaskID res_flux;
  TaskID mhd_sendf;
  TaskID mhd_recvf;
  TaskID res_copyiniuf;
  TaskID res_rkupdt;
  TaskID mhd_restu;
  TaskID mhd_sendu;
  TaskID mhd_recvu;
  TaskID res_efld;
  TaskID mhd_sende;
  TaskID mhd_recve;
  TaskID res_copyinie;
  TaskID res_ct;
  TaskID mhd_restb;
  TaskID mhd_sendb;
  TaskID mhd_recvb;
  TaskID mhd_bcs;
  TaskID mhd_prol;
  TaskID res_copyu;
  TaskID mhd_c2p;
  TaskID mhd_newdt;
  TaskID mhd_csend;
  TaskID mhd_crecv;
  TaskID res_eta;
  TaskID mhd_finalnewdt;
};

//----------------------------------------------------------------------------------------
//! \class Resistivity
//  \brief data and functions that implement various resistive physics

class Resistivity {
 public:
  Resistivity(MeshBlockPack *pp, ParameterInput *pin);
  ~Resistivity();
  template <class V> void DtRetryVisit(V &v);   // resistivity_dtretry_visit.hpp

  // data
  Real dtnew;
  std::string iso_resist_type;  // "constant", "perna" or "eos"
  Real eta_ohm_const;
  DvceArray4D<Real> eta_b;      // total resistivity of non-ideal MHD
  Real max_eta;
  bool use_rkg_sts;             // RKG super-time-stepping
  // <mhd>/rkg_dt_safety (0 = off): the dt_diff that sets the RKG stage count is safety x
  // the forward-Euler limit 2/lambda of the operator, from the per-cell bound
  //   lambda_c = 4 eta_r (1/dx1)(1/dc- + 1/dc+)/2 + 4 eta_h (1/dx2^2 + 1/dx3^2)/s^2
  // (dc = centroid spacings; eta_r, eta_h = max of eta over the 3x3x3 block; s = sine of
  // the cubed-sphere cell angle, 1 elsewhere), not multiplied by cfl_number. Off: the
  // usual cfl x (1/6) min_d dx_d^2/eta. See NewTimeStepRKGFE.
  Real rkg_dt_safety = 0.0;
  bool rkg_dt_fe = false;
  Real LambdaBoundMax();
  void NewTimeStepRKGFE();
  // conducting walls: the resistive edge EMF tangential to the reflect faces
  // (<mhd>/wall_conducting = resistive|full) and/or to the is face of a user inner x1
  // wall is 0, and so is the resistive Poynting flux through those faces
  bool wall_cond_reflect = false;
  bool wall_cond_ix1 = false;
  // track_x1_bflx (set by the problem generator): RKUpdate follows the energy flux
  // through the is and ie+1 faces along the RKG recurrence, for ledgers that book the x1
  // face fluxes of the RK stages only: rkg_bflx(m, 0|1, k, j) = F_0 (is | ie+1),
  // (2|3) = E_{j-2}, (4|5) = E_{j-1}. After the last super-stage (4|5) is the flux the
  // whole super-step applied per unit dt and area.
  bool track_x1_bflx = false;
  DvceArray4D<Real> rkg_bflx;
  // <mhd>/rkg_lean: the intermediate super-stages (stage < s) skip what only the last
  // stage needs (conserved-variable halo exchange, EOS inversion, timestep) and copy
  // only the energy slot of the conserved registers. Uniform cubed-sphere meshes only.
  bool rkg_lean = false;
  // <mhd>/rkg_strang: the super-step is split in two halves of dt/2, one before and one
  // after the RK stages (Driver::Execute), so the resistive operator is Strang-split.
  // Each half takes its stage count from dt/2 and the dt_diff of the cycle start; eta
  // (and ad_coef) is refreshed from w0 before the second half (RefreshEta). The rkg_bflx
  // flux of the cycle is the mean of the two halves (slots 6|7 hold the first).
  bool rkg_strang = false;
  // <mhd>/rkg_strang_fuse (needs rkg_strang): the second half of cycle n and the first
  // half of cycle n+1 are one super-step over (dt_n + dt_{n+1})/2 at the start of cycle
  // n+1 (see Driver::Execute). The super-step that leads a cycle is rkg_half = 3 (eta
  // refreshed first, rkg_bflx kept in slots 6|7); the closing half of a sync cycle is
  // rkg_half = 2.
  bool rkg_strang_fuse = false;
  int rkg_half = 0;        // set by the driver: 0 = one full-dt super-step, 1|2 = half,
                           // 3 = the leading super-step of a fused cycle
  Real rkg_dtfac = 1.0;    // set by the driver: fraction of pmesh->dt it spans
                           // (fused: (dt_n/2 + dt_{n+1}/2)/dt_{n+1})
  void RefreshEta();
  // <mhd>/cs_resist_x1_centred: on a cubed sphere with a stretched radial grid, take the
  // radial derivative of the Stokes loop at the x1 face instead of at the midpoint of
  // the two cell centres (resistivity_gnomonic.cpp)
  bool cs_resist_x1_centred;
  // <mhd>/sp_resist_x1_centred: the same face-centred radial derivative for the
  // spherical-polar current (current_density.hpp)
  bool sp_resist_x1_centred = false;
  Real tau, alpha, w1, mu, nu, mut, gat;
  int s = 0;   // RKG stage count; set by the first TotStage()
  Real bjm2, bjm1, bj;
  DvceArray5D<Real> u_ideal;        // conserved variables, third register
  DvceFaceFld4D<Real> b_ideal;      // face-centered magnetic fields, third register
  DvceFaceFld5D<Real> uflx_ideal;   // fluxes of conserved variables, third register
  DvceArray5D<Real> u2;             // conserved variables, fourth register
  DvceFaceFld4D<Real> b2;           // face-centered magnetic fields, fourth register

  DvceEdgeFld4D<Real> efld_resist;  // edge-centered non-ideal electric field (eta J)
  DvceEdgeFld4D<Real> jnorm;        // cubed sphere: J in the face-normal frame, the
                                    // intermediate of the two-pass curl
  DvceEdgeFld4D<Real> efld_ideal;   // edge-centered electric fields

  // ambipolar diffusion, <mhd>/ambipolar = none (default) | constant | eos. The edge EMF
  // gains eta_A (J - (J.B) B/B^2), eta_A = ad_coef |B|^2 per cell, capped at max_eta_ad:
  //   constant: ad_coef = <mhd>/eta_ad (Athena++ convention eta_A = eta_ad B^2)
  //   eos:      ad_coef = pres_cgs/(gamma_in rho_i rho_n) (cgs), rho_i = x_e n m_ion,
  //             gamma_in = <mhd>/ad_sigv_in/(m_ion + m_n); see SetAmbipolar()
  int ad_type = 0;              // 0 none, 1 constant, 2 eos
  int xe_idx = -1;              // index of the advected x_e scalar in w, -1 = local Saha
  Real xe_nhfac = 0.0;          // n_H = xe_nhfac rho_cgs (the scalar holds n_e/n_H)
  Real eta_ad = 0.0;            // <mhd>/eta_ad (constant)
  Real ad_sigv_in = 1.9e-9;     // <mhd>/ad_sigv_in [cm^3/s], Langevin rate
  Real ad_mion_amu = 30.0;      // <mhd>/ad_mion_amu, mean ion mass [amu]
  Real max_eta_ad = 1.0e300;    // <mhd>/max_eta_ad, cap on eta_A alone
  DvceArray4D<Real> ad_coef;    // eta_A/|B|^2 per cell (code units)
  DvceArray4D<Real> ad_eta_c;   // min(ad_coef |B|^2, cap) per cell, scratch
  DvceEdgeFld4D<Real> jcur;     // edge current density J (ambipolar only)
  bool OhmicNeedsUpdate() const {
    return (iso_resist_type.compare("perna") == 0 || iso_resist_type.compare("eos") == 0);
  }
  bool AmbipolarNeedsUpdate() const { return (ad_type == 2); }
  void AddEMFAmbipolar(const DvceFaceFld4D<Real> &b0, const DvceArray5D<Real> &bc,
                       DvceEdgeFld4D<Real> &efld);
  void AddEMFAmbipolarGnomonic(const DvceFaceFld4D<Real> &b0, const DvceArray5D<Real> &bc,
                               DvceEdgeFld4D<Real> &efld);
  void SetAmbipolar(const DvceArray5D<Real> &w, const EOS_Data &eos, const Real &Rgas,
                    const int il, const int iu, const int jl, const int ju,
                    const int kl, const int ku);

  // container to hold names of TaskIDs
  ResistivityTaskIDs id;
  void AssembleResistRKGTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  TaskStatus TotStage(Driver *d, int stage);
  TaskStatus RKGCoeff(Driver *d, int stage);
  TaskStatus Fluxes(Driver *d, int stage);
  TaskStatus CopyIniConsAndFluxes(Driver *d, int stage);
  TaskStatus RKUpdate(Driver *d, int stage);
  TaskStatus EField(Driver *d, int stage);
  TaskStatus CopyIniE(Driver *d, int stage);
  TaskStatus CT(Driver *d, int stage);
  TaskStatus CopyCons(Driver *d, int stage);
  TaskStatus SendFluxRKG(Driver *d, int stage);
  TaskStatus RecvFluxRKG(Driver *d, int stage);
  TaskStatus UpdateResistivity(Driver *d, int stage);
  TaskStatus ConToPrimRKG(Driver *d, int stage);
  TaskStatus NewTimeStepRKG(Driver *d, int stage);

  // functions to add resistive E-Field and energy flux
  void AddResistiveEMFs(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void AddResistiveFluxes(const DvceFaceFld4D<Real> &b0, const DvceArray5D<Real> &bc,
                          DvceFaceFld5D<Real> &flx);
  void AddEMFGeneralResist(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void AddEMFGnomonicResist(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  // <mesh>/cs_vertex_3cell: J (jnorm) at the cube-vertex edges without the corner block:
  // the vertex radial J (j1 = true, radial cells il..iu) and the corner cell's two
  // transverse edges (j1 = false, faces il..iu)
  void CsVertexJFix(const DvceFaceFld4D<Real> &b0, const int il, const int iu,
                    const bool j1);
  void AddEMFDirect(const DvceEdgeFld4D<Real> &efld_resist, DvceEdgeFld4D<Real> &efld);
  void AddFluxGeneralResist(const DvceFaceFld4D<Real> &b, const DvceArray5D<Real> &bc,
                            DvceFaceFld5D<Real> &flx);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos_data);
  void NewTimeStepGeneralResist(const DvceArray5D<Real> &w, const EOS_Data &eos_data);
  void ClearResistiveEMFs(DvceEdgeFld4D<Real> &efld);
  void SetResistivity(const DvceArray5D<Real> &w, const EOS_Data &eos, const Real &Rgas,
                      DvceArray4D<Real> &eta_b, const int il, const int iu, const int jl,
                      const int ju, const int kl, const int ku);

  //! \fn void ResistivityEOS
  //! \brief Ohmic diffusivity from the electron fraction x_e of the EOS. Two collision
  //! channels, added because the collision frequencies add:
  //!   eta_en = 230 sqrt(T)/x_e   electrons off neutrals (Blaes & Balbus 1994)
  //!   eta_ei = 5.2e11 lnL/T^1.5  electrons off ions (Spitzer, lnL = 20)
  //! Static, with max_eta passed in, because it is called inside a KOKKOS_LAMBDA (a
  //! non-static member function would capture `this`).
  KOKKOS_INLINE_FUNCTION
  static void ResistivityEOS(const Real &xe, const Real &T, const Real &max_eta,
                             Real &eta) {
    Real xemin = 230.0*sqrt(T)/max_eta;
    Real xu = (xe > xemin) ? xe : xemin;
    Real eta_en = 230.0*sqrt(T)/xu;
    Real eta_ei = 5.2e11*20.0/(T*sqrt(T));
    eta = eta_en + eta_ei;
    eta = (eta < max_eta) ? eta : max_eta;
    return;
  }

  //! \fn void ResistivityPerna
  //! \brief Ohmic diffusivity of a weakly ionized gas from the thermal ionization of
  //! potassium (Perna et al. 2010, after Balbus & Hawley 2000). Static, as above.
  KOKKOS_INLINE_FUNCTION
  static void ResistivityPerna(const Real &nn, const Real &T, const Real &max_eta,
                               Real &eta) {
    Real ak = 1.0e-7;
    Real td25 = sqrt(sqrt(T/1.0e3));
    Real lgxe = log10(6.47e-13/(1.15e-11)*sqrt(ak/1.0e-7)*td25*td25*td25)
                + 0.5*log10(2.4e15) - 0.5*log10(nn) - 25188.0/T/log(10.0);
    Real min_xe = (230.0*sqrt(T))/max_eta;
    lgxe = (lgxe < log10(min_xe)) ? log10(min_xe) : lgxe;
    eta = 230.0*sqrt(T)/pow(10.0,lgxe);
    return;
  }

 private:
  MeshBlockPack* pmy_pack;
};

#endif // DIFFUSION_RESISTIVITY_HPP_
