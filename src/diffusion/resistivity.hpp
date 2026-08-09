#ifndef DIFFUSION_RESISTIVITY_HPP_
#define DIFFUSION_RESISTIVITY_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity.hpp
//  \brief Contains data and functions that implement various non-ideal MHD (resistive)
//  processes, such as Ohmic diffusion. TODO(@user): add ambipolar diffusion, Hall effect

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

  // data
  Real dtnew;
  std::string iso_resist_type;  // "constant" or "perna"
  Real eta_ohm_const;
  DvceArray4D<Real> eta_b; // total resistivity of non-ideal MHD
//  Real min_xe;
  Real max_eta;
  bool use_rkg_sts; // RKG super-stepping (Mattia+2026)
  Real tau, alpha, w1, mu, nu, mut, gat;
  int s;
  Real bjm2, bjm1, bj;
  DvceArray5D<Real> u_ideal;       // conserved variables, third register
  DvceFaceFld4D<Real> b_ideal;     // face-centered magnetic fields, third register
  DvceFaceFld5D<Real> uflx_ideal;   // fluxes of conserved quantities on cell faces, third register
  DvceArray5D<Real> u2;       // conserved variables, fourth register
  DvceFaceFld4D<Real> b2;     // face-centered magnetic fields, fourth register
    
  DvceEdgeFld4D<Real> efld_resist;   // edge-centered electric fields due to non-ideal effects (E_{resistive} = \eta J)
  DvceEdgeFld4D<Real> efld_ideal;   // edge-centered electric fields
    
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
    TaskStatus UpdateResistivity(Driver *d, int stage);

  // functions to add resistive E-Field and energy flux
  void AddResistiveEMFs(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void AddResistiveFluxes(const DvceFaceFld4D<Real> &b0, const DvceArray5D<Real> &bc, DvceFaceFld5D<Real> &flx);
  void AddEMFConstantResist(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void AddEMFGeneralResist(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void AddEMFDirect(const DvceEdgeFld4D<Real> &efld_resist, DvceEdgeFld4D<Real> &efld);
  void AddFluxConstantGridResist(const DvceFaceFld4D<Real> &b, DvceFaceFld5D<Real> &flx);
  void AddFluxGeneralResist(const DvceFaceFld4D<Real> &b, const DvceArray5D<Real> &bc, DvceFaceFld5D<Real> &flx);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos_data);
  void NewTimeStepConstantGridResist(const DvceArray5D<Real> &w, const EOS_Data &eos_data);
  void NewTimeStepGeneralResist(const DvceArray5D<Real> &w, const EOS_Data &eos_data);
  void ClearResistiveEMFs(DvceEdgeFld4D<Real> &efld);
  void SetResistivity(const DvceArray5D<Real> &w, const EOS_Data &eos, const Real &Rgas, DvceArray4D<Real> &eta_b, const int il, const int iu, const int jl, const int ju, const int kl, const int ku);

  KOKKOS_INLINE_FUNCTION
  void ResistivityPerna(const Real &nn, const Real &T, Real &eta) {
      // Perna+2010 <- Balbus & Hawley 2000
      Real ak = 1.0e-7;
      Real td25 = sqrt(sqrt(T/1.0e3));
      Real x = 25188.0/T;
      Real lgxe = log10(6.47e-13/(1.15e-11)*sqrt(ak/1.0e-7)*td25*td25*td25)+0.5*log10(2.4e15)-0.5*log10(nn)-25188.0/T/log(10.0);
      Real min_xe = (230.0*sqrt(T))/max_eta;
      lgxe = (lgxe < log10(min_xe)) ? log10(min_xe) : lgxe;
      eta = 230.0*sqrt(T)/pow(10.0,lgxe);
      return;
  };

    KOKKOS_INLINE_FUNCTION
    void ResistivityKumar(const Real &lgrho,
                          const Real &lgT,
                          Real &eta)
    {
        const int DEG_RHO = 3;
        const int DEG_T   = 10;
        
        const Real coeff[44] = {
            17756.83765359371,     -173325.8092641056,      168767.0872366536,       -5723.276575903336,
           -37326.19648478865,      790498.1743027489,     -840928.4180831918,        43325.59798308437,
           -49248.27029547974,    -1478318.436231536,     1764541.005654571,        -126982.0747314041,
           234130.2791272859,      1446036.229095392,    -2021980.922938012,        194005.3996851212,
          -310074.1307858078,      -773268.0679357855,    1371509.834590079,       -170580.4558143873,
           207008.8037814973,       205890.886646132,     -557875.2688617986,        88339.33302144373,
           -75659.75112439394,     -13471.97879598774,     131593.3843760589,       -26241.94440105469,
            14728.0225312468,       -4914.572538661793,    -16574.12962411476,        4194.810984646928,
            -1361.682380805944,        915.8663770325294,      969.022257281818,       -373.4843052758166,
               45.32716559295613,       -41.45292106666464,      -18.95627191853248,        46.76985780921621,
               -0.1374597443224571,       0.1358831125439432,       -2.447290454239667,        -6.506619810083309
       };

        Real rho_power[DEG_RHO+1];
        rho_power[0] = 1.0;

        for(int m=1; m<=DEG_RHO; m++) {
            rho_power[m] = rho_power[m-1]*lgrho;
        }


        Real T_power[DEG_T+1];
        T_power[0] = 1.0;

        for(int k=1; k<=DEG_T; k++) {
            T_power[k] = T_power[k-1]*lgT;
        }


        Real result = 0.0;

        for(int k=0; k<=DEG_T; k++)
        {
            Real a = 0.0;

            for(int m=0; m<=DEG_RHO; m++) {
                a += coeff[k*4+m]*rho_power[m];
            }

            result += a*T_power[DEG_T-k];
        }


        eta = 1.0e4 /
              (4.0*M_PI*1.0e-7) /
               pow(10.0, result);
        
        return;
    };

 private:
  MeshBlockPack* pmy_pack;
};

#endif // DIFFUSION_RESISTIVITY_HPP_
