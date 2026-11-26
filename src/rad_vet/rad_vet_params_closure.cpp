//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_params_closure.cpp
//! \brief keys of the closure: rays, rotation, boundary emission.
//!
//!  called from: rad_vet.cpp, rad_vet_closure_faces.cpp, rays/sc_cart_setup.cpp,
//!    rays/sc_cs_sweep.cpp, rays/sc_sph_sweep.cpp
//!  calls into: rays/sc_directions.cpp

#include <math.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "utils/recipe.hpp"

namespace radvet {

extern bool g_vra_adv;   // vet_rotate_advance (rays/sc_directions.cpp)

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetBcInit
//! \brief x1 boundary intensity of the formal solutions, one rule for vet_sc and vet_gd.
//!   vet_bc_inner = auto | diff | emit (alias vet_gd_inner): rays entering at a Marshak
//!     or flux inner face; diff = the bottom cell's diffusion intensity (vet_bc_diff),
//!     emit = the face's own emission (bath implicit_ebath_x1min, or black sphere 4 F/c),
//!     auto = emit at Marshak, diff at flux. Restart files without it: vet_sc auto if
//!     vet_bc_bath, else diff; vet_gd diff.
//!   vet_bc_reflect = auto | mirror | open at a reflecting x1 end: mirror = specular
//!     (H = 0 at the face); open = diffusion at x1min, vacuum at x1max; auto = diffusion
//!     at x1min, mirror at x1max on vet_sc. vet_gd cannot mirror (its global set has no
//!     specular images). Restart files without it: open.

void RadiationVet::VetBcInit(ParameterInput *pin, const bool gd) {
  const bool rst = global_variable::restart_run;
  const std::string kin = gd ? VetKeyName(pin, "vet_bc_inner", {"vet_gd_inner"})
                             : std::string("vet_bc_inner");
  const bool bnamed = pin->DoesParameterExist("rad_vet", "vet_bc_bath");
  const bool bath = bnamed ? pin->GetBoolean("rad_vet", "vet_bc_bath") : !rst;
  std::string vi;
  if (pin->DoesParameterExist("rad_vet", kin)) {
    vi = pin->GetString("rad_vet", kin);
  } else if (!gd && (rst || (bnamed && !bath))) {
    vi = bath ? "auto" : "diff";     // vet_sc: the vet_bc_bath rule
  } else if (gd && rst) {
    vi = "diff";                     // vet_gd: vet_gd_inner unset
  } else {
    vi = pin->GetOrAddString("rad_vet", "vet_bc_inner", "auto");
  }
  if (vi.compare("diff") == 0) {
    vet_bc_inner = 0;
  } else if (vi.compare("emit") == 0) {
    vet_bc_inner = 1;
  } else if (vi.compare("auto") == 0) {
    vet_bc_inner = 2;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_vet>/" << kin << " = '" << vi << "' is not a choice "
              << "(auto | diff | emit)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string vr;
  if (pin->DoesParameterExist("rad_vet", "vet_bc_reflect")) {
    vr = pin->GetString("rad_vet", "vet_bc_reflect");
  } else if (rst) {
    vr = "open";
  } else {
    vr = pin->GetOrAddString("rad_vet", "vet_bc_reflect", "auto");
  }
  if (vr.compare("auto") == 0) {
    // inner end: diffusion intensity (the isotropic part of H = 0); outer end: mirror
    vet_bc_mirror = false;
    vet_bc_mirror_hi = !gd;
  } else if (vr.compare("mirror") == 0) {
    if (gd) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_bc_reflect = mirror is not available on "
                << "vet_gd (its direction set is global: the specular image of a "
                << "direction is not in it); use auto or open" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    vet_bc_mirror = true;
    vet_bc_mirror_hi = true;
  } else if (vr.compare("open") == 0) {
    vet_bc_mirror = false;
    vet_bc_mirror_hi = false;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_vet>/vet_bc_reflect = '" << vr << "' is not a choice "
              << "(auto | mirror | open)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetBcEmit
//! \brief vet_bc_inner resolves to emit at the current inner x1 face

bool RadiationVet::VetBcEmit() const {
  if (ibc_x1min == VET_IBC_OUTFLOW) {return true;}   // clean: nothing enters
  if (vet_bc_inner == 1) {
    return (ibc_x1min == VET_IBC_MARSHAK) || (ibc_x1min == VET_IBC_FLUX);
  }
  if (vet_bc_inner == 2) {return (ibc_x1min == VET_IBC_MARSHAK);}
  return false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotAvgParse
//! \brief <rad_vet>/vet_rotate_avg (+ _weight) and vet_rotate_avg_ema*, after the
//! closure's rotation keys

void RadiationVet::RotAvgParse(ParameterInput *pin) {
  // vet_rotate_avg = K (read only when named; < 2 off): closure averaged over the last
  // K rotated sets (rad_vet_rotavg.cpp)
  if (pin->DoesParameterExist("rad_vet","vet_rotate_avg")) {
    vra_k = pin->GetInteger("rad_vet","vet_rotate_avg");
    if (vra_k > VET_RA_NF) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_rotate_avg must be <= " << VET_RA_NF
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    vra_on = (vra_k >= 2) && (vgd_on || vet_sc);
    const bool rot = vgd_on ? VgdRotOn() : vsc_rot;
    if (vra_on && !rot && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_rotate_avg = " << vra_k << ": the direction set does "
                << "not rotate (vet_gd_rotate_dt / _every 0, or pinned): inert"
                << std::endl;
    }
    if (vra_on && rot && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_rotate_avg = " << vra_k << ": the "
                << (vgd_on ? "vet_gd" : "vet_sc") << " closure is the mean over the "
                << "last " << vra_k << " rotated sets" << std::endl;
    }
  }
  // vet_rotate_avg_ema = N (default 128 for a fresh run whose set rotates; 0 off): the
  // exponential moving average over the rotated sets
  const bool nw = !global_variable::restart_run;
  const bool erot = nw && (vgd_on || vet_sc) && (vgd_on ? VgdRotOn() : vsc_rot);
  recipe::Old(pin, "rad_vet", "vet_rotate_avg_ema", "0");
  if (erot || pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema")) {
    vra_ema_n = pin->GetOrAddReal("rad_vet","vet_rotate_avg_ema",128.0);
  }
  if (vra_ema_n >= 1.0 && (vgd_on || vet_sc)) {
    vra_ema = true;
    vra_k = 3;   // slot 0 the EMA, slot 1 the change-test aux, slot 2 the Delta
    {
      // vet_rotate_avg_ema_change_form = linear | step
      const std::string fm = pin->GetOrAddString("rad_vet",
                                                 "vet_rotate_avg_ema_change_form",
                                                 "linear");
      if (fm.compare("step") == 0) {
        vra_ema_chstep = true;
      } else if (fm.compare("linear") != 0) {
        std::cout << "### FATAL ERROR <rad_vet>/vet_rotate_avg_ema_change_form = '" << fm
                  << "' is not a choice (linear | step)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    vra_ema_wb = vgd_on;
    recipe::Rule(pin, "rad_vet", "vet_rotate_avg_ema_wb", vgd_on ? "true" : "false");
    vra_on = true;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_rotate_avg_ema = " << vra_ema_n << ": "
                << "the closure is an exponential moving average over the rotated sets "
                << "(a = 1/N per block; change test eps " << kEmaEps << ")"
                << "; vet_rotate_avg is not used"
                << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ReadClosureKeys
//! \brief the closure: rays, rotation, boundary emission

void RadiationVet::ReadClosureKeys(ParameterInput *pin) {
  const bool nw = !global_variable::restart_run;
  // the fresh default of a scheme key, on a restart without it too (fix)
  auto fix = [&](const char *k, const char *v, const char *o) {
    return nw || recipe::Fix(pin, "rad_vet", k, v, o);
  };
  MeshBlockPack *ppack = pmy_pack;
  const bool fresh = !global_variable::restart_run;
  // closure
  vet_sc = false;
  vet_full = false;
  vet_mbs = nullptr;
  vet_nguard = 0.0;
  vet_ncell = 0.0;
  vet_guard_max = 0.0;
  vet_d11_min = 1.0;
  vet_time = 0.0;
  vet_itime = 0.0;
  vet_ncall = 0.0;
  vray_ready = false;
  // sc: short characteristics, on a Cartesian mesh and (vgd_on) on spherical polar and
  // the cubed sphere.  An older name is read as sc only where it ran the same closure:
  // vet_sc on a Cartesian mesh, vet_col with vet_gd = true on a curved one; with an
  // analytic clean_closure there also vet_col and, in a restart file, m1, minerbo, kershaw
  // and eddington (a file without the key ran m1)
  {std::string cl = pin->GetOrAddString("rad_vet","closure",
                                        global_variable::restart_run ? "m1" : "sc");
  if (cl.compare("sc") != 0) {
    const bool an = sph_geom && VetAnalytic(pin);
    const bool gd = sph_geom && pin->DoesParameterExist("rad_vet","vet_gd") &&
                    pin->GetBoolean("rad_vet","vet_gd");
    const bool oldcl = (cl.compare("m1") == 0) || (cl.compare("minerbo") == 0) ||
                       (cl.compare("kershaw") == 0) || (cl.compare("eddington") == 0);
    if (!((cl.compare("vet_sc") == 0 && !sph_geom) ||
          (cl.compare("vet_col") == 0 && (gd || an)) ||
          (oldcl && an && global_variable::restart_run))) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/closure = '" << cl << "' was removed (only sc: the "
        << "rays, on every mesh; vet_sc is read as sc on a Cartesian mesh, vet_col with "
        << "vet_gd = true on spherical polar and the cubed sphere; the analytic closures "
        << "are clean_closure = eddington | sphere).  Name closure = sc to change to the "
        << "rays" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "### NOTE: <rad_vet>/closure = " << cl << " is read as sc (the same "
                << "closure under its new name)" << std::endl;
    }
    pin->SetString("rad_vet","closure","sc");
  }
  if (!sph_geom) {
    // a Cartesian mesh (rad_vet_vet.cpp)
    vet_sc = true;
    VetRotInit(pin);   // the vet_gd rotation keys and defaults
    // vet_gd_rotate_every has no vet_sc clock: inert
    if (pin->DoesParameterExist("rad_vet","vet_gd_rotate_every") &&
        pin->GetInteger("rad_vet","vet_gd_rotate_every") > 0 &&
        global_variable::my_rank == 0) {
      std::cout << "### NOTE: <rad_vet>/vet_gd_rotate_every has no effect on a "
                << "Cartesian mesh (the set does not rotate); name vet_gd_rotate_dt"
                << std::endl;
    }
    RotAvgParse(pin);  // vet_rotate_avg (shared key)
    // clean_rot_ring = live: the rotating set's cross-fade on the clean ring
    vgd_xf = vsc_rot && cln_ring == 1 && VgdRotOn();
    if (VgdRotOn()) {
      recipe::Rule(pin, "rad_vet", "vet_gd_rotate_xfade", vgd_xf ? "true" : "false");
    }
    // vet_source_noesrc (rad_vet_vet.cpp)
    vsrc_noes = pin->GetOrAddBoolean("rad_vet","vet_source_noesrc",
                                     fix("vet_source_noesrc", "true", "false"));
  } else {
    // spherical polar, cubed sphere (rad_vet_vetray.cpp, rad_vet_vetgd.cpp)
    vgd_on = VetCurvedRays(pin);
    if (pin->DoesParameterExist("rad_vet","vet_source_noesrc")) {
      vsrc_noes = pin->GetBoolean("rad_vet","vet_source_noesrc");
    }
    // vet_bc_diff = face | cell: inner diffusion intensity of the formal solutions at the
    // face (S + mu(-dS/dtau) with S of the face, 1.5 S_a - 0.5 S_b) or from the first
    // cell
    {
      const char *bdd = fix("vet_bc_diff", "face", "cell") ? "face" : "cell";
      std::string bd = pin->GetOrAddString("rad_vet","vet_bc_diff",bdd);
      if (bd.compare("face") == 0) {
        vet_bc_face = true;
      } else if (bd.compare("cell") == 0) {
        vet_bc_face = false;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/vet_bc_diff must be cell or face"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // Athena++ x2rat / x3rat (not AthenaK keys) noted for VetRayInit
    for (const char *rk : {"x2rat", "x3rat"}) {
      if (pin->DoesParameterExist("mesh", rk) && pin->GetReal("mesh", rk) != 1.0) {
        vray_rat_keys += std::string(" mesh/") + rk + " = " + pin->GetString("mesh", rk);
      }
    }
    // default of a fresh start on spherical polar: a fixed set (vet_gd_rotate_dt = 0)
    if (vgd_on && ppack->pmesh->use_spherical_polar && !global_variable::restart_run &&
        recipe::NewDefaults(pin) && !pin->DoesParameterExist("rad_vet","vet_gd_rotate_dt")
        && !pin->DoesParameterExist("rad_vet","vet_gd_rotate_every")) {
      pin->GetOrAddString("rad_vet","vet_gd_rotate_dt","0");
      recipe::Log("<rad_vet>/vet_gd_rotate_dt = 0 (default on spherical polar)");
    }
    // default set of a fresh start that names no set key: 3-D spherical polar the
    // staggered GL 7 x 18 (126 directions); 2-D axisymmetric with a phi extent below 2 pi
    // the staggered GL 9 x nphi, nphi nearest 12 with 2 pi/nphi dividing the extent
    // (fatal if none); else the geodesic vet_gd_geo = 3 (92). A restart naming none
    // keeps HEALPix nside 3.
    if (vgd_on && !pin->DoesParameterExist("rad_vet","vet_gd_gl_nmu") &&
        !pin->DoesParameterExist("rad_vet","vet_gd_ls") &&
        !pin->DoesParameterExist("rad_vet","vet_gd_geo")) {
      const auto &ms = ppack->pmesh->mesh_size;
      const double per = ms.x3max - ms.x3min;
      const bool ax2 = ppack->pmesh->use_spherical_polar && !ppack->pmesh->three_d &&
          std::abs(per - 2.0*M_PI) >= 1.0e-8;
      const bool fresh = recipe::NewDefaults(pin) && !global_variable::restart_run;
      const bool sp3 = fresh && ppack->pmesh->use_spherical_polar && !ax2;
      int np = sp3 ? 18 : 0;
      for (int q = 4; q <= 64 && fresh && ax2; ++q) {
        const double m = per*q/(2.0*M_PI);
        if (std::round(m) >= 1.0 && std::fabs(m - std::round(m)) < 1.0e-9 &&
            (np == 0 || std::abs(q - 12) < std::abs(np - 12))) {np = q;}
      }
      if (fresh && ax2 && np == 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_vet>/vet_gd: 2-D axisymmetric, no 2 pi/nphi (nphi 4..64) "
          << "divides the phi extent " << per << ": name the set (vet_gd_gl_nmu and "
          << "vet_gd_gl_nphi)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (fresh && !ax2 && !sp3) {
        pin->GetOrAddInteger("rad_vet","vet_gd_geo",3);
      }
      if (fresh && (ax2 || sp3)) {
        pin->GetOrAddInteger("rad_vet","vet_gd_gl_nmu",sp3 ? 7 : 9);
        pin->GetOrAddInteger("rad_vet","vet_gd_gl_nphi",np);
        pin->GetOrAddBoolean("rad_vet","vet_gd_gl_stagger",true);
      }
      recipe::Log(!fresh ? std::string("<rad_vet>/vet_gd: HEALPix nside 3 (restart input "
                                       "has no key)")
          : (ax2 ? "<rad_vet>/vet_gd_gl_nmu = 9, vet_gd_gl_nphi = " + std::to_string(np) +
                   ", stagger (default on 2-D axisymmetric sp: exact phi walls)"
           : sp3 ? std::string("<rad_vet>/vet_gd_gl_nmu = 7, vet_gd_gl_nphi = 18, "
                               "stagger (default on 3-D spherical polar)")
                 : std::string("<rad_vet>/vet_gd_geo = 3 (default)")));
    }
    // a fresh vet_gd run whose set rotates and names no vet_rotate_avg averages 4 rotated
    // sets (with the time clock: the cross-fade)
    if (vgd_on && recipe::NewDefaults(pin) && !global_variable::restart_run &&
        !pin->DoesParameterExist("rad_vet","vet_rotate_avg")) {
      const bool rdt = !pin->DoesParameterExist("rad_vet","vet_gd_rotate_dt") ||
          pin->GetString("rad_vet","vet_gd_rotate_dt") == "auto" ||
          pin->GetReal("rad_vet","vet_gd_rotate_dt") > 0.0;
      const bool rev = pin->DoesParameterExist("rad_vet","vet_gd_rotate_every") &&
          pin->GetInteger("rad_vet","vet_gd_rotate_every") > 0;
      if (rdt || rev) {
        pin->GetOrAddInteger("rad_vet","vet_rotate_avg",4);
        recipe::Log("<rad_vet>/vet_rotate_avg = 4 (default with a rotating vet_gd set)");
      }
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_gl_nmu")) {
      vgd_gl_nmu = pin->GetInteger("rad_vet","vet_gd_gl_nmu");
    }
    if (vgd_on && vgd_gl_nmu > 0) {
      vgd_gl_nphi = pin->GetOrAddInteger("rad_vet","vet_gd_gl_nphi",vgd_gl_nphi);
      vgd_gl_stag = pin->GetOrAddBoolean("rad_vet","vet_gd_gl_stagger",vgd_gl_stag);
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_ls")) {
      vgd_ls = pin->GetInteger("rad_vet","vet_gd_ls");
    }
    // vet_gd_geo = n: geodesic set (vgd_ls = -n)
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_geo")) {
      const int ng_ = pin->GetInteger("rad_vet","vet_gd_geo");
      if (ng_ > 0) {vgd_ls = -ng_;}
    }
    // one direction set: the product set takes precedence over vet_gd_ls / vet_gd_geo, so
    // naming both would silently ignore the latter (rule 12: FATAL, name one; 0 = off)
    if (vgd_on && vgd_gl_nmu > 0 && vgd_ls != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
                << "<rad_vet>: both a product set (vet_gd_gl_nmu = " << vgd_gl_nmu
                << ") and " << ((vgd_ls > 0) ? "vet_gd_ls" : "vet_gd_geo") << " are named; "
                << "the product set would be used and the other ignored. Name one set "
                << "(vet_gd_ls = 0 / vet_gd_geo = 0 for the product set)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_rotate_every")) {
      vgd_rot = pin->GetInteger("rad_vet","vet_gd_rotate_every");
    }
    // rotation keys (VetRotInit): vet_gd_rotate_dt = <dt_rot> | auto (20 x the first
    // step's dt, recorded) on a physical-time clock; default auto on a fresh vet_gd run
    // naming neither _dt nor _every and no vet_gd_beam. Sub-keys _t0, _align, _stretch,
    // _slope.
    VetRotInit(pin);   // shared with vet_sc (rad_vet_vetgd.cpp)
    // clean sp / cs, fresh, rotating on the time clock: vet_rotate_avg = 4 (recorded)
    const bool clk = vgd_on && VgdRotOn() && (vgd_rot_dt > 0.0 || vgd_rot_auto);
    if (clk && nw) pin->GetOrAddInteger("rad_vet","vet_rotate_avg",4);
    RotAvgParse(pin);   // also after vet_sc's VetRotInit
    // the mean of sets rotating in time is cross-faded from the current to the next ray
    // set over each block (two sweeps per build, fused where VetGdXffOk allows), except
    // under clean_rot_ring = jump
    vgd_xf = vgd_xff = clk && vra_on && !(cln_ring_named && cln_ring != 1);
    if (VgdRotOn()) {
      recipe::Rule(pin, "rad_vet", "vet_gd_rotate_xfade", vgd_xf ? "true" : "false");
    }
    // the cubed sphere with clean_rays = sc: vet_gd's own sweep (intensities to cln_fvi)
    cln_cssc = pmy_pack->pmesh->use_cubed_sphere && vgd_on;
    // (vet_rotate_avg_ema_wb: the ring carries the cells' half ranges, CLN_HO, so that a
    // set change keeps them in a static state, as the moments)
    cln_rfh = cln_cssc && vra_ema_wb;
    // clean_rot_ring unnamed: live where the vet_gd rays rotate in time with a cross-fade
    // (decided before a restart file's channels are counted)
    if (!cln_ring_named && vgd_on && vra_on && g_vra_adv &&
        vgd_xf && vra_k <= VET_RA_NF) {
      cln_ring = 1;
    }
    // vet_bc_bath: the top shell's inward rays start from the incident bath
    // implicit_ebath_x1max of a Marshak outer x1
    if (vgd_on) {
      vgd_tbath = pin->GetOrAddBoolean("rad_vet", "vet_bc_bath",
                                       fix("vet_bc_bath", "true", "false"));
    }
    // vet_bc_inner / vet_bc_reflect shared with vet_sc (VetBcInit; alias vet_gd_inner).
    if (vgd_on) {VetBcInit(pin, true);}
    if (vgd_on) vgd_allow_clamp =
        pin->GetOrAddBoolean("rad_vet","vet_gd_allow_clamp",vgd_allow_clamp);
    // an sp full sphere: the global band
    {
      Mesh *pm = ppack->pmesh;
      vgd_gband = vgd_on && pm->use_spherical_polar && pm->use_polar_boundary &&
      pin->GetOrAddString("mesh", "ix2_bc", "") == "polar" &&
      pin->GetOrAddString("mesh", "ox2_bc", "") == "polar" &&
      pin->GetOrAddString("mesh", "ix3_bc", "") == "periodic" &&
      std::abs(pm->mesh_size.x3max - pm->mesh_size.x3min - 2.0*M_PI) < 1.0e-8;
      if (vgd_on) {
        recipe::Rule(pin, "rad_vet", "vet_gd_band_global", vgd_gband ? "true" : "false");
      }
    }
    // vet_gd_band_ring (read only when named, VetGdInit): band intensities lent per shell
    // from a few slabs; unnamed: on where it applies
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_band_ring")) {
      vgd_ring_par = pin->GetBoolean("rad_vet","vet_gd_band_ring");
    }
    vgd_hcomp = vgd_on ? 2 : 0;   // (VetGdInit switches it off where it does not apply)
    // defaults of a fresh 3-D spherical-polar run: wall_interp; with a fixed product set
    // wall_sym, and theta_seam_sym on a theta-periodic wedge centred on the equator with
    // a set that is its equatorial mirror (odd nmu, or no stagger)
    const auto &msz = ppack->pmesh->mesh_size;
    const bool w3 = nw && vgd_on && ppack->pmesh->use_spherical_polar &&
                    ppack->pmesh->three_d && !ppack->pmesh->use_polar_boundary;
    const bool wgl = w3 && vgd_gl_nmu > 0 && !VgdRotOn();
    if (vgd_on) vgd_wint = pin->GetOrAddBoolean("rad_vet","vet_gd_wall_interp",w3);
    if (vgd_on) vgd_wsym = pin->GetOrAddBoolean("rad_vet","vet_gd_wall_sym",
        wgl && vgd_wint);
    if (vgd_on) vgd_wst = pin->GetOrAddBoolean("rad_vet","vet_gd_theta_seam_sym",
        wgl && vgd_wint && (vgd_gl_nmu%2 == 1 || !vgd_gl_stag) &&
        ppack->pmesh->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic &&
        std::fabs(msz.x2min + msz.x2max - M_PI) <= 1.0e-12);
    // vet_gd_theta_wall: periodic | mirror (specular) | vacuum (no incoming intensity at
    // a physical theta wall) | extrap (zero gradient)
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_theta_wall")) {
      const std::string tw = pin->GetString("rad_vet","vet_gd_theta_wall");
      if (tw.compare("mirror") == 0) {
        vgd_thmir = true;
      } else if (tw.compare("vacuum") == 0) {
        vgd_thvac = 1;
      } else if (tw.compare("extrap") == 0) {
        vgd_thvac = 2;
      } else if (tw.compare("periodic") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/vet_gd_theta_wall must be periodic, mirror,"
                  << " vacuum or extrap" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (vgd_on && pmy_pack->pmesh->use_spherical_polar &&
               !pmy_pack->pmesh->use_polar_boundary) {
      // unnamed: a reflecting theta wall mirrors, a periodic one wraps
      vgd_thmir = (pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x2] ==
                   BoundaryFlag::reflect);
    }
    // exact halo lists where they apply (VetGdInit switches them off elsewhere)
    vgd_hx_on = true;
    // vet_gd_halo_ipc: default on CUDA, HIP without a theta wall (VetGdSetup keeps it on
    // the raw 3-D sweep only; off-node partners take MPI)
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    vgd_ipc_on = vgd_on && !vgd_thmir && vgd_thvac == 0;
#endif
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_halo_ipc")) {
      vgd_ipc_on = pin->GetBoolean("rad_vet","vet_gd_halo_ipc");
    }
    // vet_col_lat_cut_thin (CleanRayCut)
  }
  }
}

}  // namespace radvet
