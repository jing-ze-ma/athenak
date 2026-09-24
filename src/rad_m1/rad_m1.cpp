//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1.cpp
//! \brief constructor, parameter reading and allocation for the RadiationM1 class

#include <float.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

RadiationM1::RadiationM1(MeshBlockPack *ppack, ParameterInput *pin) :
    u0("m1_u0",1,1,1,1,1),
    u1("m1_u1",1,1,1,1,1),
    coarse_u0("m1_coarse_u0",1,1,1,1,1),
    uflx("m1_uflx",1,1,1,1,1),
    pmy_pack(ppack) {
  // (1) parameters
  c_light = pin->GetReal("rad_m1","c_light");
  if (!(c_light > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/c_light must be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  chat_over_c = pin->GetOrAddReal("rad_m1","chat_over_c",1.0);
  if (!(chat_over_c > 0.0) || chat_over_c > 1.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/chat_over_c must be in (0,1]" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  chat = chat_over_c*c_light;
  // ---- MILESTONE 3a hook: which transport scheme runs.  Everything the implicit one
  // needs is read and allocated by ImplicitInit at the end of this constructor.
  {std::string tr = pin->GetOrAddString("rad_m1","transport","explicit");
  if (tr.compare("explicit") == 0) {
    transport = M1_TRANSPORT_EXPLICIT;
  } else if (tr.compare("implicit_x1") == 0) {
    transport = M1_TRANSPORT_IMPLICIT_X1;
  } else if (tr.compare("implicit") == 0) {
    transport = M1_TRANSPORT_IMPLICIT;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = '" << tr << "' is not a valid choice "
      << "(explicit | implicit_x1 | implicit)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // GEOMETRY GUARD (docs/dev/rad_m1_curvilinear_design.md sect. 0).  Every M1 kernel,
  // explicit and implicit, is CARTESIAN with a UNIFORM cell: divergences divide by
  // mb_size.dx1..3 (rad_m1_update.cpp, rad_m1_fluxes.cpp, rad_m1_implicit.cpp), the
  // flux vector is Cartesian, and no metric enters.  On a spherical-polar or cubed-sphere
  // mesh the run would proceed silently with that arithmetic (on sp dx2 is Delta theta
  // in radians: the transverse transport is wrong by r and r sin(theta)), so refuse it.
  // The radial stretches are refused too although they only act inside the sp/cs
  // geometry (on a plain Cartesian mesh the hydro ignores them and so would M1), so
  // that no stretched input can reach the uniform-dx kernels; and relativistic
  // coordinates, since the M1 here is the O(v/c) flat-space system.
  //
  // STAGE S1 (tests_m1/runs_5a_sp_s1) lifts the refusal for ONE configuration: a
  // spherical-polar WEDGE that does not touch a pole, with the radial stretches allowed,
  // transport = implicit, closure = eddington and the restrictions checked by
  // SphericalS1Check after ImplicitInit; STAGE S2 (tests_m1/runs_5b_sp_s2) adds the
  // chi(f) closures m1 / minerbo / kershaw there.  Everything else on sp (the poles,
  // explicit transport, vet_sc, tau) and everything on the cubed sphere stays refused.
  sph_geom = false;
  {
    Mesh *pm_ = ppack->pmesh;
    std::string why;
    const bool sp_ = pm_->use_spherical_polar && !pm_->use_cubed_sphere;
    if (pm_->use_spherical_polar && !sp_) {why += " mesh/use_spherical_polar";}
    if (pm_->use_cubed_sphere) {why += " mesh/use_cubed_sphere";}
    if (pm_->use_polar_boundary) {
      why += " mesh/use_polar_boundary (the poles are not supported yet: run a wedge "
             "whose theta range stays clear of theta = 0 and pi)";
    } else if (sp_ && (pm_->mesh_size.x2min <= 1.0e-10 ||
                       pm_->mesh_size.x2max >= M_PI - 1.0e-10)) {
      why += " a theta range that reaches a pole (the poles are not supported yet: run "
             "a wedge whose theta range stays clear of theta = 0 and pi)";
    }
    if (pm_->use_grid_stretch_r && !sp_) {why += " mesh/use_grid_stretch_r";}
    if (pm_->use_grid_stretch_r_poly && !sp_) {why += " mesh/use_grid_stretch_r_poly";}
    if (pm_->use_grid_stretch_theta) {why += " mesh/use_grid_stretch_theta";}
    if (ppack->pcoord != nullptr &&
        (ppack->pcoord->is_special_relativistic ||
         ppack->pcoord->is_general_relativistic ||
         ppack->pcoord->is_dynamical_relativistic)) {
      why += " relativistic coordinates";
    }
    if (!why.empty()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1> (explicit and implicit transport) supports a "
        << "uniform Cartesian mesh and, for the implicit solve only (stages S1, S2: "
        << "closure eddington | m1 | minerbo | kershaw), "
        << "a spherical-polar wedge clear of the poles; this input sets:" << why << "."
        << std::endl
        << "The M1 kernels would run with Cartesian uniform-dx arithmetic on it.  See "
        << "docs/dev/rad_m1_curvilinear_design.md (branch m1-curv-design) for the "
        << "staged plan that lifts this." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    sph_geom = sp_;
  }
  nstage = M1_NSTAGE;
  impl_cfl = -1.0;
  impl_tol = 1.0e-8;
  impl_maxit = 30;
  impl_opac_update = false;
  impl_allow_multid = false;
  marshak_q = 0.5;
  ibc_x1min = M1_IBC_MARSHAK;
  ibc_x1max = M1_IBC_MARSHAK;
  iflux_x1min = 0.0;
  iflux_x1max = 0.0;
  iebath_x1min = 0.0;
  iebath_x1max = 0.0;
  impl_nstep = 0.0;
  impl_itsum = 0.0;
  impl_itmax = 0.0;
  impl_nfail = 0.0;
  // ---- milestone 3b phase B
  impl_solver = M1_ISOLV_LINE_JACOBI;
  impl_linsum = 0.0;
  impl_linmax = 0.0;
  impl_lin_tol = 1.0e-10;
  trans_on = false;
  trans_x3 = false;
  pbval_th = nullptr;
  pbval_tq = nullptr;
  halo_shell = false;
  // ---- milestone 3b phase C
  bicg_on = false;
  impl_lin_maxit = 200;
  pbval_kr = nullptr;
  bcg_nsolve = 0.0;
  bcg_itsum = 0.0;
  bcg_itmax = 0.0;
  bcg_nbreak = 0.0;
  bcg_nfall = 0.0;
  bcg_nred = 0.0;
  impl_bcg_sync = 0;
  impl_plog = 0;
  bcg_r0rel = 0.0;
  impl_lres_test = true;
  impl_conv_est = false;
  impl_ew_max = 0.0;
  impl_ew_gam = 0.9;
  impl_lin_cnorm = 0.0;
  ew_fprev = 0.0;
  ew_etaprev = 0.0;
  impl_pred = false;
  pred_ok = false;
  pred_dt = 0.0;
  // ---- milestone 3b phase D
  impl_offdiag = M1_OD_LAGGED;
  od_now = M1_OD_LAGGED;
  od_nfall = 0.0;
  od_emin = 1.0e300;
  impl_vimp = false;
  vimp_now = false;
  iw_vimp = -1;
  vimp_nfall = 0.0;
  vimp_emin = 1.0e300;
  pbval_vm = nullptr;
  time_scheme = 0;
  t2_ok = false;
  t2_afmode = 0;
  t2_solve = 0;
  t2_fail = false;
  t2_dbg_fail = -1;
  t2_nstep = 0.0;
  t2_nbe = 0.0;
  t2_nfall = 0.0;
  t2_dtprev = 0.0;
  t2_vprev = false;
  t2_vext = true;
  t2_nclip = 0.0;
  pred2_ok = false;
  pred2_dt = 0.0;
  impl_crelax = 1.0;
  impl_crelax_thin = false;
  impl_clag_step = false;
  impl_tlim = M1_TLIM_NONE;
  impl_tfmax = 1.0;
  // ---- milestone 3e: Anderson acceleration of the Picard map
  impl_accel = M1_IACC_NONE;
  impl_and_m = 5;
  impl_and_beta = 1.0;
  impl_and_start = 1;
  aa_nc = 0;
  aa_nh = 0;
  aa_head = 0;
  aa_hasp = false;
  aa_fnp = -1.0;
  aa_nacc = 0.0;
  aa_nrst = 0.0;
  // ---- milestone 3g: the gas-radiation energy coupling
  impl_gas_newton = false;
  impl_eos_cache = false;
  impl_line_solver = 0;
  impl_pcr_team = 0;
  impl_pcr_check = false;
  impl_halo_direct = false;
  halo_direct_on = false;
  impl_odc = false;
  impl_stencil = false;
  st_edges = true;
  impl_prec_float = false;
  impl_kfuse = 0;
  impl_prec = 0;
  pcr_chk_max = 0.0;
  pcr_chk_n = 0.0;
  impl_ecnt = 2;
  impl_eccheck = true;
  iw_gas = -1;
  impl_nec = 0;
  newt_nfb = 0.0;
  gas_ncell = 0.0;
  ec_nmiss = 0.0;
  ec_emax = 0.0;
  ec_tmax = 0.0;
  // ---- end of the 3a hook
  cfl_rad = pin->GetOrAddReal("rad_m1","cfl_rad",0.4);
  e_floor = pin->GetOrAddReal("rad_m1","e_floor",(FLT_MIN));
  subcycle = pin->GetOrAddBoolean("rad_m1","subcycle",true);

  // thick-limit flux (design sect. 3, and the 1b correction to it recorded in
  // rad_m1_closure.hpp: ap_hll is the BLEND form, which is what survives PLM).
  thick_flux_str = pin->GetOrAddString("rad_m1","thick_flux","none");
  if (thick_flux_str.compare("none") == 0) {
    thick_flux = M1_THICK_NONE;
  } else if (thick_flux_str.compare("ap_hll") == 0) {
    thick_flux = M1_THICK_APHLL;
  } else if (thick_flux_str.compare("scaled") == 0) {
    thick_flux = M1_THICK_SCALED;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/thick_flux = '" << thick_flux_str
      << "' is not a valid choice (none | ap_hll | scaled)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // which algebraic form of the ap_hll E-flux (design sect. 10; milestone 1c).
  // `unified` needs no reconstruct-dependent switch, but it is MEASURABLY worse on
  // every 1b gate that discriminates (T3 at tau_cell = 10, the resolution order, the
  // Nyquist decay rate, T6), so the 1b pair stays the default: see rad_m1_closure.hpp.
  ap_form_str = pin->GetOrAddString("rad_m1","ap_form","alpha2");
  if (ap_form_str.compare("unified") == 0) {
    ap_form = M1_APFORM_UNIFIED;
  } else if (ap_form_str.compare("alpha2") == 0) {
    ap_form = M1_APFORM_ALPHA2;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/ap_form = '" << ap_form_str
      << "' is not a valid choice (unified | alpha2)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  scaled_pref = pin->GetOrAddReal("rad_m1","scaled_prefactor",20.0);
  if (!(scaled_pref > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/scaled_prefactor must be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (1b) opacities, per unit mass, in code units.  kappa_e defaults to kappa_p and
  // kappa_f to kappa_p as well, so a single kappa_p gives a grey absorbing medium.
  {std::string op = pin->GetOrAddString("rad_m1","opacity","const");
  if (op.compare("const") == 0) {
    opacity_type = M1_OPAC_CONST;
  } else if (op.compare("powerlaw") == 0) {
    opacity_type = M1_OPAC_POWERLAW;
  } else if (op.compare("user") == 0) {
    opacity_type = M1_OPAC_USER;
  } else if (op.compare("table") == 0) {
    opacity_type = M1_OPAC_TABLE;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/opacity = '" << op << "' is not a valid choice "
      << "(const | powerlaw | user | table)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // the unit conversions of the TABLE lookup, which is in kelvin and g/cm^3 (design
  // sect. 4 / the stage-2 plan sect. 4).  T[K] = T_code*temp_unit_kelvin: the same
  // direction Conduction uses with Units::temperature_cgs(), which for the He box's
  // <units> 1/1/1 with mu = 1 is mu m_u/k_B = 1.202724e-8 K per code unit.
  otab.tunit = pin->GetOrAddReal("rad_m1","temp_unit_kelvin",1.0);
  otab.dunit = pin->GetOrAddReal("rad_m1","rho_unit_cgs",1.0);
  otab.kunit = pin->GetOrAddReal("rad_m1","kappa_unit",1.0);
  if (!(otab.tunit > 0.0) || !(otab.dunit > 0.0) || !(otab.kunit > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/temp_unit_kelvin, rho_unit_cgs and kappa_unit must all "
      << "be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  kappa_p = pin->GetOrAddReal("rad_m1","kappa_p",0.0);
  kappa_e = pin->GetOrAddReal("rad_m1","kappa_e",kappa_p);
  kappa_f = pin->GetOrAddReal("rad_m1","kappa_f",kappa_p);
  kappa_s = pin->GetOrAddReal("rad_m1","kappa_s",0.0);
  opac_rho_ref = pin->GetOrAddReal("rad_m1","rho_ref",1.0);
  opac_t_ref = pin->GetOrAddReal("rad_m1","t_ref",1.0);
  opac_a = pin->GetOrAddReal("rad_m1","opac_a",0.0);
  opac_b = pin->GetOrAddReal("rad_m1","opac_b",0.0);
  if (!(opac_rho_ref > 0.0) || !(opac_t_ref > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/rho_ref and t_ref must be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (kappa_p < 0.0 || kappa_e < 0.0 || kappa_f < 0.0 || kappa_s < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> opacities must be non-negative" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  opac_zero = (kappa_p == 0.0 && kappa_e == 0.0 && kappa_f == 0.0 && kappa_s == 0.0 &&
               opacity_type != M1_OPAC_USER && opacity_type != M1_OPAC_TABLE);

  // (1c) matter coupling.  It reads rho, v and the gas energy from hydro's CONSERVED
  // u0 and writes u0(IEN) and u0(IM1..3) back, so it needs a <hydro> block unless every
  // opacity vanishes (the pure-transport tests of milestone 1a).
  bool have_hydro = pin->DoesBlockExist("hydro");
  // (1c) the advective enthalpy-flux split of the E equation (design sect. 3 "Moving
  // fluid").  It is what carries (4/3) E v where alpha -> 0 switches the transport flux
  // off, so it is on by default exactly where alpha exists and there is a medium; with
  // thick_flux = scaled | none it is available for comparison but off by default.
  advect_split = pin->GetOrAddBoolean("rad_m1","advect_split",
                                      have_hydro && (thick_flux == M1_THICK_APHLL));
  // (1c-B) which velocity builds the enthalpy flux at the face.  `cell` is the 1c-A
  // behaviour and is first order there; `recon` reconstructs v with the same method as
  // (E, f_i), which is what makes A second order in a sheared flow.
  {std::string sv = pin->GetOrAddString("rad_m1","split_vel","recon");
  if (sv.compare("recon") == 0) {
    split_vel_recon = true;
  } else if (sv.compare("cell") == 0) {
    split_vel_recon = false;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/split_vel = '" << sv << "' is not a valid choice "
      << "(recon | cell)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // the O(v/c) CONTROL of T4/T4b.  Never a production setting: it is here to exhibit the
  // spurious heating ~ (4/3) c rho kappa beta^2 E that the full source form cancels.
  {std::string sf = pin->GetOrAddString("rad_m1","source_form","full");
  if (sf.compare("full") == 0) {
    source_ovc = false;
  } else if (sf.compare("ovc") == 0) {
    source_ovc = true;
    std::cout << "### WARNING: <rad_m1>/source_form = ovc truncates the matter coupling "
      << "to O(v/c); this is the FAILING CONTROL of T4/T4b, not a physical option"
      << std::endl;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/source_form = '" << sf << "' is not a valid choice "
      << "(full | ovc)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // (1d) the interface (well-balanced) form of the implicit F source, design sect. 13.
  // Default `cell` reproduces milestone 1c bit for bit.
  {std::string fs = pin->GetOrAddString("rad_m1","f_source","cell");
  if (fs.compare("cell") == 0) {
    f_source_wb = false;
  } else if (fs.compare("wb") == 0) {
    f_source_wb = true;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/f_source = '" << fs << "' is not a valid choice "
      << "(cell | wb)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  coupling = pin->GetOrAddBoolean("rad_m1","coupling",!opac_zero);
  gas_feedback = pin->GetOrAddBoolean("rad_m1","gas_feedback",true);
  // debug switches (see rad_m1.hpp); every default is the production value
  opac_freeze = pin->GetOrAddBoolean("rad_m1","opac_freeze",false);
  opac_frozen = false;
  dbg_gas_force = pin->GetOrAddBoolean("rad_m1","dbg_gas_force",true);
  dbg_gas_heat = pin->GetOrAddBoolean("rad_m1","dbg_gas_heat",true);
  dbg_gas_force_trans = pin->GetOrAddBoolean("rad_m1","dbg_gas_force_trans",true);
  dbg_trans_memory = pin->GetOrAddReal("rad_m1","dbg_trans_memory",1.0);
  // DIAGNOSTIC, VET scaffolding (tests_m1/runs_3i_tensor): the Eddington tensor of the
  // multi-D implicit solve is NOT taken from the cell's own flux.
  //   frozen  computed once from the initial state (selected closure), never updated
  //   tilt    frozen, with the tensor axis rotated by tilt*sin(2 pi x2/L2) about x3
  //   tau     rebuilt EVERY step from the evolving gas: chi = exact grey K/J at
  //           the column optical depth, axis along -grad tau; never reads the local flux
  {std::string st = pin->GetOrAddString("rad_m1","dbg_tensor","none");
  dbg_tensor_init = false;
  dbg_tensor_tilt = pin->GetOrAddReal("rad_m1","dbg_tensor_tilt",0.3);
  if (st.compare("none") == 0) {
    dbg_tensor = 0;
  } else if (st.compare("frozen") == 0) {
    dbg_tensor = 1;
  } else if (st.compare("tilt") == 0) {
    dbg_tensor = 2;
  } else if (st.compare("tau") == 0) {
    dbg_tensor = 3;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/dbg_tensor = '" << st << "' not implemented "
      << "(none | frozen | tilt | tau)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (dbg_tensor != 0) {
    std::cout << "<rad_m1> DIAGNOSTIC dbg_tensor = " << st << " (tilt amplitude "
              << dbg_tensor_tilt << ")" << std::endl;
  }
  }
  if ((opac_freeze || !dbg_gas_force || !dbg_gas_heat || !dbg_gas_force_trans) &&
      global_variable::my_rank == 0) {
    std::cout << "### WARNING: <rad_m1> DEBUG switches are active: opac_freeze="
              << (opac_freeze ? "true" : "false") << " dbg_gas_force="
              << (dbg_gas_force ? "true" : "false") << " dbg_gas_heat="
              << (dbg_gas_heat ? "true" : "false") << " dbg_gas_force_trans="
              << (dbg_gas_force_trans ? "true" : "false")
              << " -- this is not a physical configuration" << std::endl;
  }
  if (coupling && !have_hydro) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> matter coupling requires a <hydro> block" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!opac_zero && !coupling) {
    std::cout << "### WARNING: <rad_m1> has non-zero opacities but coupling = false: "
      << "the opacity enters the thick-limit FLUX only" << std::endl;
  }
  // an isothermal gas has no energy equation to couple to
  if (have_hydro && (coupling || !opac_zero)) {
    std::string heos = pin->GetOrAddString("hydro","eos","ideal");
    if (heos.compare("isothermal") == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1> needs a hydro EOS with an energy equation, not "
        << "'isothermal'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // the design's gas-only EOS requirement (sect. 1)
  if (have_hydro && pin->GetOrAddBoolean("hydro","eos_radiation",false)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> requires a GAS-ONLY EOS: <hydro>/eos_radiation must be "
      << "false (the radiation energy and pressure are evolved, not in the EOS)"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // radiation constant, code units: the equilibrium energy density is arad*T^4 with T
  // the EOS's own code temperature.  Required whenever there is emission/absorption.
  if (coupling && (kappa_p > 0.0 || kappa_e > 0.0 ||
                   opacity_type == M1_OPAC_USER ||
                   opacity_type == M1_OPAC_TABLE)) {
    arad = pin->GetReal("rad_m1","arad");
    if (!(arad > 0.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/arad must be positive" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    arad = pin->GetOrAddReal("rad_m1","arad",0.0);
  }

  // (2a) <rad_m1>/force_reference.  `wb_arad` means: an EXTERNAL well-balanced scheme
  // already delivers rho*arad_ref(z) to the gas momentum (box_convection's effective
  // potential), so the coupling must hand the gas only the RESIDUAL radiative force.
  // The reference array itself comes from the problem generator (SetForceReference);
  // without it the option is a fatal, not a silent no-op.
  {std::string fr = pin->GetOrAddString("rad_m1","force_reference","none");
  if (fr.compare("none") == 0) {
    force_ref = M1_FREF_NONE;
  } else if (fr.compare("wb_arad") == 0) {
    force_ref = M1_FREF_WB_ARAD;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/force_reference = '" << fr << "' is not a valid choice "
      << "(none | wb_arad)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  if (force_ref != M1_FREF_NONE && !(coupling && gas_feedback)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/force_reference needs coupling = true and "
      << "gas_feedback = true: it modifies the momentum handed to the gas" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (2a) the RSLA start-up check of design sect. 2
  rsla_warn = pin->GetOrAddReal("rad_m1","rsla_warn",0.1);
  rsla_vmax = pin->GetOrAddReal("rad_m1","rsla_vmax",-1.0);
  rsla_force = pin->GetOrAddBoolean("rad_m1","rsla_force",false);
  rsla_done = false;

  // closure
  vet_sc = false;
  vet_full = false;
  vet_mbs = nullptr;
  vet_eig_min = 0.0;
  vet_nguard = 0.0;
  vet_ncell = 0.0;
  vet_guard_max = 0.0;
  vet_d11_min = 1.0;
  // read only when given, so that an input without it writes a byte-identical restart
  dbg_opac_patch = 1.0;
  dbg_opac_x1lo = dbg_opac_x1hi = dbg_opac_x2lo = dbg_opac_x2hi = 0.0;
  if (pin->DoesParameterExist("rad_m1","dbg_opac_patch")) {
    dbg_opac_patch = pin->GetReal("rad_m1","dbg_opac_patch");
    dbg_opac_x1lo = pin->GetReal("rad_m1","dbg_opac_x1lo");
    dbg_opac_x1hi = pin->GetReal("rad_m1","dbg_opac_x1hi");
    dbg_opac_x2lo = pin->GetReal("rad_m1","dbg_opac_x2lo");
    dbg_opac_x2hi = pin->GetReal("rad_m1","dbg_opac_x2hi");
  }
  vet_time = 0.0;
  vet_itime = 0.0;
  vet_ncall = 0.0;
  tau_closure = false;
  tau_ready = false;
  tau_time = 0.0;
  tau_ncall = 0.0;
  vet_col = false;
  vcol_sph = false;
  vcol_axis_flux = false;
  vcol_nc = vcol_np = vcol_nmu = vcol_every = vcol_nray = vcol_dump_every = 0;
  vcol_built = false;
  vcol_time = vcol_ncall = vcol_nskip = 0.0;
  vcol_sq = false;
  vcol_team = false;
  vcol_qmin = vcol_qmax = 0.0;
  vcol_lc = 1;
  vcol_ts = vcol_lcin = 0;
  {std::string cl = pin->GetOrAddString("rad_m1","closure","m1");
  chi_kind = M1_CHI_LEVERMORE;
  if (cl.compare("m1") == 0) {
    eddington = false;
  } else if (cl.compare("minerbo") == 0) {
    eddington = false;
    chi_kind = M1_CHI_MINERBO;
  } else if (cl.compare("kershaw") == 0) {
    eddington = false;
    chi_kind = M1_CHI_KERSHAW;
  } else if (cl.compare("eddington") == 0) {
    eddington = true;
  } else if (cl.compare("vet_sc") == 0) {
    // the implicit multi-D solve reads (chi, n) from the short-characteristics formal
    // solution (rad_m1_vet.cpp); every other use of chi (explicit wave speeds, the
    // 1-D branch) keeps Levermore's
    eddington = false;
    vet_sc = true;
  } else if (cl.compare("tau") == 0) {
    // the multi-D implicit solve reads (chi, n) from the column optical depth
    // (rad_m1_tau.cpp), in the uniaxial form; every other use of chi (explicit wave
    // speeds, the 1-D branch, the explicit coupling) keeps Levermore's, as closure = m1
    eddington = false;
    tau_closure = true;
  } else if (cl.compare("vet_col") == 0) {
    // STAGE S5 (rad_m1_vetcol.cpp): the tensor of a per-column 1-D formal solution,
    // carried by the tau closure's machinery (tau_ten, built once per step)
    eddington = false;
    tau_closure = true;
    vet_col = true;
    vcol_nc = pin->GetOrAddInteger("rad_m1","vet_col_ncore",8);
    vcol_np = pin->GetOrAddInteger("rad_m1","vet_col_nsub",1);
    vcol_nmu = pin->GetOrAddInteger("rad_m1","vet_col_nmu",4);
    vcol_every = pin->GetOrAddInteger("rad_m1","vet_col_every",1);
    vcol_sq = pin->GetOrAddBoolean("rad_m1","vet_col_surface_q",false);
    vcol_qmin = pin->GetOrAddReal("rad_m1","vet_col_surface_qmin",1.0e-3);
    vcol_qmax = pin->GetOrAddReal("rad_m1","vet_col_surface_qmax",1.0);
    vcol_team = pin->GetOrAddBoolean("rad_m1","vet_col_team",true);
    vcol_ts = pin->GetOrAddInteger("rad_m1","vet_col_team_size",0);
    vcol_lcin = pin->GetOrAddInteger("rad_m1","vet_col_chunk",0);
    std::string ax = pin->GetOrAddString("rad_m1","vet_col_axis","radial");
    if (ax.compare("flux") == 0) {
      vcol_axis_flux = true;
    } else if (ax.compare("radial") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/vet_col_axis = '" << ax << "' not implemented "
        << "(radial | flux)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    vcol_dump = "";
    if (pin->DoesParameterExist("rad_m1","vet_col_dump")) {
      vcol_dump = pin->GetString("rad_m1","vet_col_dump");
      vcol_dump_every = pin->GetOrAddInteger("rad_m1","vet_col_dump_every",1);
      if (vcol_dump.compare("none") == 0 || vcol_dump_every <= 0) {vcol_dump = "";}
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/closure = '" << cl << "' not implemented "
      << "(m1 | minerbo | kershaw | eddington | vet_sc | tau | vet_col)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (vet_sc && transport != M1_TRANSPORT_IMPLICIT) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/closure = vet_sc needs <rad_m1>/transport = implicit "
      << "on a multi-D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (tau_closure && transport != M1_TRANSPORT_IMPLICIT) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/closure = tau needs <rad_m1>/transport = implicit "
      << "on a multi-D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (chi_kind != M1_CHI_LEVERMORE && transport == M1_TRANSPORT_EXPLICIT) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/closure = '" << cl << "' needs an implicit transport: "
      << "the explicit HLL wave speeds are Levermore's" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }

  // reconstruction of (E, f_i).  dc was added in 1b so that the thick-limit blend can
  // be compared with the first-order scheme it reduces to (Berthon & Turpault / Bloch
  // et al.); the three high-order kernels of src/reconstruct/ were added in 1c-B.
  recon_str = pin->GetOrAddString("rad_m1","reconstruct","plm");
  if (recon_str.compare("plm") == 0) {
    recon_method = ReconstructionMethod::plm;
    recon_code = M1_RECON_PLM;
  } else if (recon_str.compare("dc") == 0) {
    recon_method = ReconstructionMethod::dc;
    recon_code = M1_RECON_DC;
  } else if (recon_str.compare("ppm4") == 0) {
    recon_method = ReconstructionMethod::ppm4;
    recon_code = M1_RECON_PPM4;
  } else if (recon_str.compare("ppmx") == 0) {
    recon_method = ReconstructionMethod::ppmx;
    recon_code = M1_RECON_PPMX;
  } else if (recon_str.compare("wenoz") == 0) {
    recon_method = ReconstructionMethod::wenoz;
    recon_code = M1_RECON_WENOZ;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/reconstruct = '" << recon_str << "' not implemented "
      << "(dc | plm | ppm4 | ppmx | wenoz)" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // PLM of (E, f_i) reads two cells beyond the face, so two ghost zones are the
  // minimum; the 5-cell kernels need the face stencil i-3..i+2, hence three.
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int ng_need = (recon_code >= M1_RECON_PPM4) ? 3 : 2;
  if (indcs.ng < ng_need) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/reconstruct = " << recon_str << " requires at least "
      << ng_need << " ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (2) who owns the mesh timestep.  With sub-cycling on, the radiation CFL is NOT
  // folded into Mesh::NewTimeStep -- unless nothing else sets a timestep at all, in
  // which case the mesh step would be unbounded.
  bool other_sets_dt = (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd") ||
                        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("particles"));
  sets_mesh_dt = (!subcycle) || (!other_sets_dt);

  // (3) PD-ARS explicit weights.  With S == 0 the tableau (design sect. 5) is
  //   U1   = U^n + dt L(U^n)
  //   U^n1 = U^n + dt/2 [L(U^n) + L(U1)]
  // i.e. Heun / SSP-RK2 written in the AthenaK gam0/gam1/beta form.
  gam0[0] = 1.0;  gam1[0] = 0.0;  beta[0] = 1.0;
  gam0[1] = 0.5;  gam1[1] = 0.5;  beta[1] = 0.5;

  dtnew = (FLT_MAX);
  dt_sub = 0.0;
  nsub = 1;

  // (4) allocate arrays
  int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(u0, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(u1, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(uflx.x1f, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(uflx.x2f, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(uflx.x3f, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  // rho*kappa is needed in the GHOST cells too: the face opacity at the first and last
  // active face is the arithmetic mean over a ghost and an active cell.
  Kokkos::realloc(opac, nmb, M1_NOPAC, ncells3, ncells2, ncells1);
  Kokkos::deep_copy(opac, 0.0);
  if (coupling) {
    Kokkos::realloc(ugas1, nmb, 4, ncells3, ncells2, ncells1);
  }
  // force_reference: allocated (and zeroed) here so that the kernel always has a valid
  // View to capture; the pgen REPLACES it with its own through SetForceReference.
  if (force_ref != M1_FREF_NONE) {
    Kokkos::realloc(arad_ref, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(arad_ref, 0.0);
  }
  Kokkos::realloc(cnt, M1_NCNT);
  for (int n=0; n<M1_NCNT; ++n) {cnt.h_view(n) = 0.0;}
  cnt.modify_host();
  cnt.sync_device();

  if (ppack->pmesh->multilevel) {
    int nccells1 = indcs.cnx1 + 2*(indcs.ng);
    int nccells2 = (indcs.cnx2 > 1)? (indcs.cnx2 + 2*(indcs.ng)) : 1;
    int nccells3 = (indcs.cnx3 > 1)? (indcs.cnx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(coarse_u0, nmb, M1_NVAR, nccells3, nccells2, nccells1);
  }

  // (4b) MILESTONE 3a hook: the implicit solver's own parameters, checks and arrays.
  // Returns immediately with transport = explicit.
  ImplicitInit(pin);
  if (sph_geom) {SphericalS1Check(pin);}

  // (5) boundary buffers
  pbval_u = new MeshBoundaryValuesCC(ppack, pin, false);
  pbval_u->InitializeBuffers(M1_NVAR);
  // (F2, F3) is the tangential pair (what the default (IVY, IVZ) rule gives by luck);
  // E and F1 (radial) are scalars.  Stated explicitly, see SetVectorPairs.
  pbval_u->SetVectorPairs(M1_NVAR, {{M1_F2, M1_F3}});

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: c=" << c_light << " chat/c=" << chat_over_c
              << " cfl_rad=" << cfl_rad << " thick_flux=" << thick_flux_str
              << " ap_form=" << ap_form_str
              << " recon=" << recon_str
              << " advect_split=" << (advect_split ? "true" : "false")
              << " split_vel=" << (split_vel_recon ? "recon" : "cell")
              << " source_form=" << (source_ovc ? "ovc" : "full")
              << " closure=" << (eddington ? "eddington" : "m1")
              << " subcycle=" << (subcycle ? "true" : "false") << std::endl;
    std::cout << "         kappa_p=" << kappa_p << " kappa_e=" << kappa_e
              << " kappa_f=" << kappa_f << " kappa_s=" << kappa_s
              << " coupling=" << (coupling ? "true" : "false")
              << " gas_feedback=" << (gas_feedback ? "true" : "false")
              << " arad=" << arad << std::endl;
    std::cout << "         opacity=" << (opacity_type == M1_OPAC_TABLE ? "table" : "")
              << " force_reference="
              << (force_ref == M1_FREF_WB_ARAD ? "wb_arad" : "none")
              << " temp_unit_kelvin=" << otab.tunit
              << " rho_unit_cgs=" << otab.dunit << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::SetOpacityTables

void RadiationM1::SetOpacityTables(const DvceArray2D<Real> &kr,
                                   const DvceArray2D<Real> &kp,
                                   const DvceArray1D<Real> &lT,
                                   const DvceArray1D<Real> &lD,
                                   const int nT, const int nD) {
  otab.kr = kr;
  otab.kp = kp;
  otab.lT = lT;
  otab.lD = lD;
  otab.nT = nT;
  otab.nD = nD;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::SetForceReference

void RadiationM1::SetForceReference(const DvceArray4D<Real> &a) {
  arad_ref = a;
}

//----------------------------------------------------------------------------------------
// destructor

RadiationM1::~RadiationM1() {
  if (tau_closure) {TauClosureReport();}
  ReportCounters();
  ImplicitReport();   // milestone 3a; a no-op in transport = explicit
  Time2Report();      // time_scheme = hesdirk2; a no-op otherwise
  if (vet_sc) {VetReport();}
  VetFree();
  delete pbval_u;
  if (pbval_th != nullptr) {delete pbval_th;}
  if (pbval_tq != nullptr) {delete pbval_tq;}
  if (pbval_kr != nullptr) {delete pbval_kr;}
  if (pbval_vm != nullptr) {delete pbval_vm;}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ReportCounters
//! \brief one line at the end of the run with the cost and the health of the implicit
//! energy-exchange solve.  A non-zero failure count means cells left the bracketed
//! Newton at M1_MAXIT without reaching M1_RTOL and is a result to distrust.

void RadiationM1::ReportCounters() {
  if (!coupling) return;
  cnt.modify_device();
  cnt.sync_host();
  Real nsolve = cnt.h_view(M1_CNT_NSOLVE);
#if MPI_PARALLEL_ENABLED
  {
    Real loc[M1_NCNT], glob[M1_NCNT];
    for (int n=0; n<M1_NCNT; ++n) {loc[n] = cnt.h_view(n);}
    MPI_Reduce(loc, glob, M1_NCNT, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    Real lmax = cnt.h_view(M1_CNT_ITMAX), gmax;
    MPI_Reduce(&lmax, &gmax, 1, MPI_ATHENA_REAL, MPI_MAX, 0, MPI_COMM_WORLD);
    for (int n=0; n<M1_NCNT; ++n) {cnt.h_view(n) = glob[n];}
    cnt.h_view(M1_CNT_ITMAX) = gmax;
    nsolve = cnt.h_view(M1_CNT_NSOLVE);
  }
#endif
  if (global_variable::my_rank != 0) return;
  Real mean = (nsolve > 0.0) ? (cnt.h_view(M1_CNT_ITSUM)/nsolve) : 0.0;
  std::cout << "<rad_m1> implicit energy solve: cells=" << nsolve
            << " iterations mean=" << mean << " max=" << cnt.h_view(M1_CNT_ITMAX)
            << " FAILURES=" << cnt.h_view(M1_CNT_NFAIL)
            << " (no bracket: " << cnt.h_view(M1_CNT_NBRAK) << ")" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::SetSubsteps
//! \brief Set the substep dt and the number of substeps for a mesh step dt_mesh.
//! N_sub = ceil(dt_mesh/dt_rad) with dt_rad = cfl_rad*min(dx)/chat (= dtnew), and every
//! substep uses dt_mesh/N_sub.  When this module owns the mesh timestep there is nothing
//! to sub-cycle and N_sub = 1.

int RadiationM1::SetSubsteps(Real dt_mesh) {
  if (sets_mesh_dt || !subcycle) {
    nsub = 1;
  } else {
    nsub = 1;
    if (dtnew > 0.0 && dt_mesh > dtnew) {
      nsub = static_cast<int>(std::ceil(dt_mesh/dtnew));
      if (nsub < 1) {nsub = 1;}
    }
  }
  dt_sub = dt_mesh/static_cast<Real>(nsub);
  return nsub;
}

} // namespace radm1
