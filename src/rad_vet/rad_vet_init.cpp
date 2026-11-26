//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_init.cpp
//! \brief arrays of the solve, the device memory arena shared with the rays, and
//! CleanArgs, what the kernels read.
//!
//!  called from: rad_vet_closure.cpp, rad_vet_closure_faces.cpp,
//!    rad_vet_closure_rays.cpp, rad_vet_closure_ring.cpp, rad_vet_fct.cpp,
//!    rad_vet_flux.cpp, rad_vet_flux_bounds.cpp, rad_vet_flux_schur.cpp,
//!    rad_vet_flux_seam.cpp, rad_vet_gas.cpp, rad_vet_ledger.cpp, rad_vet_newton.cpp,
//!    rad_vet_oldvector.cpp, rad_vet_params_implicit.cpp, rad_vet_rhs.cpp,
//!    rad_vet_solver.cpp, rad_vet_steady.cpp, rad_vet_tasks.cpp, rad_vet_update.cpp,
//!    rays/sc_setup.cpp, driver/dt_retry.cpp, hydro/hydro_tasks.cpp
//!  calls into: rad_vet_geometry.cpp

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_args.hpp"

namespace radvet {

extern bool g_vra_adv;   // vet_rotate_advance (rays/sc_directions.cpp)

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanInit
//! \brief allocates the arrays of the clean solve and checks that the setup is supported

void RadiationVet::CleanInit() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb = pmy_pack->nmb_thispack;
  // clean_closure = iterate: not the closure a restart file holds (Time2RstSet)
  if (!cln_cip) {
    Kokkos::realloc(cln_s, nmb, CleanNs(pmy_pack->pmesh->use_cubed_sphere, cln_rfh),
                    n3, n2, n1);
  }
  if (cln_phic.extent_int(0) != nmb) {
    Kokkos::realloc(cln_phic, nmb, 3, n3, n2, n1);
  }
  Kokkos::realloc(cln_v, nmb, 3, n3, n2, n1);
  // staggered on the cubed sphere: cell terms and faces of CleanStagFSolve
  if (pmy_pack->pmesh->use_cubed_sphere) Kokkos::realloc(cln_vf, nmb, 18, n3, n2, n1);
  // mem_arena: the shared slots follow cln_w.o (CleanWoBind), unless the closure is built
  // inside the solve (clean_closure = iterate | extrap*)
  cln_w.s = DvceArray5D<Real>();
  if (!marena_on || cln_closit || cln_extrap > 0) {
    Kokkos::realloc(cln_w.s, nmb, CLW_NS, n3, n2, n1);
  }
  // mem_arena: cln_w.o is made at the first solve, over the intensities (CleanWoBind)
  cln_wo_nc = CLW_NX - CLW_NS;
  cln_w.o = DvceArray5D<Real>();
  if (!marena_on) Kokkos::realloc(cln_w.o, nmb, cln_wo_nc, n3, n2, n1);
  cln_w.ns = CLW_NS;
  // opacities at the iterate's T where they depend on it (tables, powerlaw T exponent)
  cln_opt = !opac_zero && (opacity_type == VET_OPAC_TABLE ||
            (opacity_type != VET_OPAC_CONST && opac_b != 0.0));
  // staggered bracket's lagged flux DG: one more slot (CleanStagBracket)
  cln_fdg = cln_consist ? CLF_NF : -1;
  cln_xqs = CLF_NF + (cln_fdg >= 0);   // clean_stag_xdonor's q
  const int nfx = CLF_NF + (cln_fdg >= 0) + 3;
  Kokkos::deep_copy(cln_w.s, 0.0);
  Kokkos::deep_copy(cln_w.o, 0.0);
  for (CleanF *f : {&cln_fc1, &cln_fc2, &cln_fc3}) {
    Kokkos::realloc(f->s, nmb, CLF_NS, n3, n2, n1);
    Kokkos::realloc(f->o, nmb, nfx - CLF_NS, n3, n2, n1);
    Kokkos::deep_copy(f->s, 0.0);
    Kokkos::deep_copy(f->o, 0.0);
  }
  for (int d = 0; d < 3 && cln_rstd.extent_int(0) == nmb && cln_fsolve; ++d) {
    auto fc = (d == 0) ? cln_fc1 : ((d == 1) ? cln_fc2 : cln_fc3);   // from a restart
    auto AL = Kokkos::ALL;
    Kokkos::deep_copy(Kokkos::subview(fc.s, AL, CLF_N + 3, AL, AL, AL),
                      Kokkos::subview(cln_rstd, AL, d, AL, AL, AL));
    Kokkos::deep_copy(Kokkos::subview(fc.s, AL, CLF_KD, AL, AL, AL),
                      Kokkos::subview(cln_rstd, AL, 3 + d, AL, AL, AL));
  }
  if ((force_ref == VET_FREF_CLEAN || force_ref == VET_FREF_WB_MD) &&
      cln_fref.extent_int(0) != nmb) {   // unless a restart has read it
    Kokkos::realloc(cln_fref, nmb, (force_ref == VET_FREF_WB_MD) ? 9 : 3 + cln_lzt, n3,
                    n2,
                    n1);
  }
  Kokkos::realloc(cln_st, nmb, CLN_ST, indcs.nx3, indcs.nx2, indcs.nx1);
  if (pmy_pack->pmesh->use_cubed_sphere) {
    Kokkos::realloc(cln_csy, nmb, CLN_CSY, indcs.nx3, indcs.nx2, indcs.nx1);
  }
  if (cln_ysum) Kokkos::realloc(cln_yl, nmb, 9, n3, n2, n1);
  // multigrid levels of the preconditioner: halve x2 and x3 where the tile's are even
  // (>= 4), clean_stag_mg_levels coarse levels at most (x1 is never halved)
  std::array<int, 3> nc = {indcs.nx1, indcs.nx2, indcs.nx3};
  if (cln_mgt[0] == 0) cln_mgt = nc;   // layout_tile: the tiles' cells
  std::array<int, 3> tc = cln_mgt;
  while (static_cast<int>(cln_mgn.size()) < cln_mgnl) {
    bool coarser = false;
    for (int d = 1; d < 3; ++d) {
      if (tc[d] >= 4 && tc[d] % 2 == 0) {
        nc[d] /= 2;
        tc[d] /= 2;
        coarser = true;
      }
    }
    if (!coarser) break;
    cln_mgn.push_back(nc);
    cln_mg.emplace_back("cln_mg", nmb, CLS_N, nc[2], nc[1], nc[0]);
  }
  if (cln_extrap > 0) {
    // a restart may have read t^n and the end of the last step
    const int ns = cln_s.extent_int(1), c3 = indcs.nx3, c2 = indcs.nx2, c1 = indcs.nx1;
    if (cln_mn.extent_int(0) != nmb) Kokkos::realloc(cln_mn, nmb, ns, c3, c2, c1);
    if (cln_mc.extent_int(0) != nmb) Kokkos::realloc(cln_mc, nmb, ns, c3, c2, c1);
    if (!CleanExtrapSx()) Kokkos::realloc(cln_mp, nmb, ns, c3, c2, c1);
    if (cln_extrap == 2) Kokkos::realloc(cln_m1, nmb, ns, c3, c2, c1);
  }
  Mesh *pm = pmy_pack->pmesh;
  if (pm->multilevel || chat != c_light) {
    CleanFatal("<rad_vet>/scheme = clean currently needs a uniform mesh and chat = c");
  }
  // boundaries: open (bath, beam; outflow without a bath), reflect, periodic and an
  // imposed x1 flux.  The imposed flux is the face's whole energy flux, taken with no E
  // weights (CleanFaceWeights, CleanFaceSource), kept as the face F by the F solve; the
  // inner half cell absorbs its momentum (CleanMomentumChange), as at open faces.
  const int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  for (int s = 0; s < 6; ++s) {
    const bool flux = bc[s] == VET_IBC_FLUX || (bc[s] == VET_IBC_X2FLUX && ifx2_on);
    if (!(bc[s] == VET_IBC_MARSHAK || bc[s] == VET_IBC_OUTFLOW || bc[s] == VET_IBC_REFLECT
          ||
          bc[s] == VET_IBC_PERIODIC || (flux && s < 2))) {
      CleanFatal("<rad_vet>/implicit_bc_* with scheme = clean: auto | marshak | outflow "
                 "| "
                 "reflect | periodic, flux (x1); not efix");
    }
  }
  if (iebath_hi_col_on) {
    CleanFatal("scheme = clean has no per-column outer bath");
  }
  if ((force_ref == VET_FREF_WB_ARAD || force_ref == VET_FREF_WB_MD) && !fref_wsplit) {
    CleanFatal("<rad_vet>/force_reference = wb_arad(_md) with scheme = clean needs the "
               "reference work given by the problem generator (force_reference_work = "
               "split)");
  }
  if (force_ref == VET_FREF_WB_MD && !cln_steady && !cln_fref_set) {
    CleanFatal("<rad_vet>/force_reference = wb_arad_md needs clean_init = steady (its "
               "reference is the steady start's force)");
  }
  if ((vgd_on || (vsc_rot && cln_ring)) && VgdRotOn()) {
    // rotating rays (CleanFvSets): the closure is the J-weighted mean of the sets, made
    // continuous across blocks by the cross-fade, no slope correction.
    // vet_rotate_avg_ema: only with clean_rot_ring = live (CleanRingEma)
    if ((vra_ema && cln_ring != 1) ||
        (vgd_xf && !(vgd_rot_dt > 0.0 || vgd_rot_auto)) ||
        (vgd_rot3 && !pm->three_d)) {
      CleanFatal("<rad_vet>/scheme = clean with a rotating vet_gd set: vet_rotate_avg = "
                 "K "
                 "(J-weighted; _ema only with clean_rot_ring = live), vet_gd_rotate_xfade with vet_gd_rotate_dt, and on a "
                 "2-D axisymmetric mesh the z axis");
    }
  }
  const bool ring_ok = (vgd_on || vsc_rot) && vra_on && g_vra_adv &&
                       (vgd_rot_dt > 0.0 || vgd_rot_auto) &&
                       vra_k <= VET_RA_NF;
  if (cln_ring && !(ring_ok && vgd_xf == (cln_ring == 1))) {
    CleanFatal("<rad_vet>/clean_rot_ring needs vet_gd or vet_sc rays rotating in time: "
               "vet_rotate_avg = K and vet_gd_rotate_dt");
  }
  // sc rays sweep one set at a time: a mean or cross-fade of sets needs the ring
  if (vgd_on && VgdRotOn() && (vra_on || vgd_xf) && !cln_ring) {
    CleanFatal("<rad_vet>/clean_rays = sc with vet_rotate_avg or vet_gd_rotate_xfade "
               "needs clean_rot_ring");
  }
  if (pm->use_cubed_sphere) {
    // rays: vet_gd's (clean_rays = sc: its band, VetGdCsSweep) or the analytic sphere
    // (eddington: infinite radius); an even panel resolution keeps cell centres off the z
    // axis (CleanCsFrame)
    if (!(cln_sphere || cln_edd || vgd_on) || (cln_consist && !vgd_on) ||
        pm->mesh_indcs.nx2 % 2 != 0) {
      CleanFatal("<rad_vet>/scheme = clean on the cubed sphere needs "
                 "clean_closure = sphere | eddington (no clean_consistency) or the "
                 "rays (closure = sc) and an even <mesh>/nx2");
    }
    if (!vgd_on) {
      if (cln_edd) cln_sph_r = std::numeric_limits<Real>::max();
      cln_sphere = true;
    }
    CleanCsGeometry();
    cln_init = true;
    return;
  }
  if (pm->use_spherical_polar) {
    // rays: vet_gd's fixed direction set swept from the first shell, or the analytic
    // sphere (clean_closure = sphere | eddington)
    if (!(cln_sphere || vgd_on || cln_edd) ||
        (cln_consist && !vgd_on)) {
      CleanFatal("<rad_vet>/scheme = clean on spherical polar needs clean_closure = "
                 "sphere | eddington (no clean_consistency) or the "
                 "rays (closure = sc)");
    }
    // theta walls: vet_gd_theta_wall or the mesh's (rad_vet.cpp); the sweep starts at the
    // first shell or the ray cut (vray_icut).  No clean solve reads the half-range ratios
    // vgd_hr, so implicit_flux_beam = halfrange has no effect.
    impl_beam_hr = false;
    cln_init = true;
    return;
  }
  if (cln_sphere || !vet_sc || !vet_full) {
    CleanFatal("<rad_vet>/scheme = clean on a Cartesian mesh needs "
               "vet_tensor = full (clean_closure = sphere is for spherical polar)");
  }
  // the formal solution provides the half-range sums and C1
  impl_beam_hr = true;
  thsw_c1 = true;
  cln_init = true;
}

//----------------------------------------------------------------------------------------
//! \fn CleanArgs RadiationVet::CleanMakeArgs
//! \brief collects the arrays and constants the device functions need

CleanArgs RadiationVet::CleanMakeArgs() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  CleanArgs x;
  x.cs = cln_s;
  x.w = cln_w;
  x.fc[0] = cln_fc1;
  x.fc[1] = cln_fc2;
  x.fc[2] = cln_fc3;
  x.beam = ibm_f;
  x.bstress = ibm_s;
  x.h = dt_sub;
  x.c = c_light;
  for (int d = 0; d < 3; ++d) x.idx[d] = 1.0/MeshBlock0Dx(d);
  const int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  // an outflow face is the open one without an entering bath
  for (int s = 0; s < 6; ++s) {
    x.bc[s] = (bc[s] == VET_IBC_X2FLUX) ? VET_IBC_FLUX :
              ((bc[s] == VET_IBC_OUTFLOW) ? VET_IBC_MARSHAK : bc[s]);
  }
  x.iflux[0] = iflux_x1min;
  x.iflux[1] = iflux_x1max;
  x.ifx2[0] = ifx2_lo;
  x.ifx2[1] = ifx2_hi;
  x.mbc = pmy_pack->pmb->mb_bcs.d_view;
  for (int n = 0; n < 6; ++n) x.ebath[n] = 0.0;
  x.ebath[0] = iebath_x1min;
  x.ebath[1] = iebath_x1max;
  for (int s = 0; s < 2; ++s) x.ebath[s] *= (bc[s] != VET_IBC_OUTFLOW);
  x.use_beam = (ibm_n > 0);
  x.feedback = cln_feedback;
  x.vzero = cln_vzero;
  x.hi = false;
  x.mcl = false;
  x.fa = false;
  x.fwc = false;
  x.fs = false;
  x.fd = cln_fsolve && cln_fkin;
  x.gf = cln_gflux;
  x.fre = cln_freal;
  x.fdg = cln_fdg;
  x.xq = cln_xd_on ? cln_xqs : -1;
  x.drag = cln_feedback;
  x.upw = cln_upw;
  x.sar = !cln_gflux;
  x.op = opac;
  x.nf = pmy_pack->pmesh->three_d ? 3 : (pmy_pack->pmesh->multi_d ? 2 : 1);
  x.is = indcs.is;
  x.ie = indcs.ie;
  x.js = indcs.js;
  x.je = indcs.je;
  x.ks = indcs.ks;
  x.ke = indcs.ke;
  x.multi_d = pmy_pack->pmesh->multi_d;
  x.three_d = pmy_pack->pmesh->three_d;
  x.sph = pmy_pack->pmesh->use_spherical_polar;
  auto &ms = pmy_pack->pmesh->mesh_size;
  x.cseam = (pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic)
            ? 0.5/tan(ms.x2min) + 0.5/tan(ms.x2max) : 0.0;   // (read on that wrap only)
  x.cub = pmy_pack->pmesh->use_cubed_sphere;
  x.ys = !x.cub && cln_ysum;
  x.lwb = x.cub && cln_latwb;
  x.yl = cln_yl;
  if (x.sph) {
    auto *pc = pmy_pack->pcoord;
    x.vol = pc->volume;
    x.area[0] = pc->area.x1f;
    x.area[1] = pc->area.x2f;
    x.area[2] = pc->area.x3f;
    x.dxf[0] = pc->dxface.x1f;
    x.dxf[1] = pc->dxface.x2f;
    x.dxf[2] = pc->dxface.x3f;
    x.x1v = pc->x1v;
    x.x1f = pc->xx1f;
    x.x2v = pc->x2v;
    x.x2f = pc->xx2f;
    const int nmb = pmy_pack->nmb_thispack, n2 = pc->x2v.extent_int(1);
    if (x.sph && (cln_tg.extent_int(0) != nmb || cln_tg.extent_int(2) != n2 + 1)) {
      Kokkos::realloc(cln_tg, nmb, 4, n2 + 1);   // (once: the trig of theta, per j)
      auto tg = cln_tg;
      auto x2v = pc->x2v, x2f = pc->xx2f;
      par_for("cln_tg", DevExeSpace(), 0, nmb - 1, 0, n2,
      KOKKOS_LAMBDA(const int m, const int j) {
        if (j < n2) {
          tg(m,0,j) = tan(x2v(m,j));
          tg(m,2,j) = sin(x2v(m,j));
        }
        tg(m,1,j) = tan(x2f(m,j));
        tg(m,3,j) = sin(x2f(m,j));
      });
    }
    x.tg = cln_tg;
    if (vgd_wst) {
      if (cln_sgv.extent_int(0) != nmb) {CleanSeamImage();}
      x.vol = cln_sgv;
      x.area[0] = cln_sga0;
      x.area[1] = cln_sga1;
      x.dxf[1] = cln_sgd1;
      x.dxf[2] = cln_sgd2;
      x.x2v = cln_sgx;
      x.tg = cln_sgtg;
    }
  }
  if (x.cub) {
    x.vol = pmy_pack->pcoord->volume;
    x.ae0 = pmy_pack->pcoord->area.x1f;
    for (int d = 0; d < 3; ++d) {
      x.area[d] = cln_ga[d];
      x.dxf[d] = cln_gd[d];
    }
    x.g = cln_gn;
    x.x1v = pmy_pack->pcoord->x1v;
    x.x1f = pmy_pack->pcoord->xx1f;
  }
  if (!x.sph && !x.cub && cln_lzsp) {   // (clean_lz_spin, Cartesian: z x r of cells)
    const int nmb = pmy_pack->nmb_thispack, ng = indcs.ng;
    const int n1 = indcs.nx1 + 2*ng, n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
    if (cln_cx[0].extent_int(0) != nmb) {
      Kokkos::realloc(cln_cx[0], nmb, n1);
      Kokkos::realloc(cln_cx[1], nmb, n2);
      auto cx = cln_cx[0], cy = cln_cx[1];
      auto size = pmy_pack->pmb->mb_size.d_view;
      const int is = indcs.is, js = indcs.js, nx1 = indcs.nx1, nx2 = indcs.nx2;
      par_for("cln_cx", DevExeSpace(), 0, nmb - 1, 0, std::max(n1, n2) - 1,
      KOKKOS_LAMBDA(const int m, const int i) {
        if (i < n1) cx(m,i) = CellCenterX(i - is, nx1, size(m).x1min, size(m).x1max);
        if (i < n2) cy(m,i) = CellCenterX(i - js, nx2, size(m).x2min, size(m).x2max);
      });
    }
    x.x1v = cln_cx[0];
    x.x2v = cln_cx[1];
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  x.etg = fl.on && fl.use_etotgrav;
  x.mhd = fl.on && fl_mhd;
  x.phi = fl.phicc0;
  x.emag = emag0;
  return x;
}

//----------------------------------------------------------------------------------------
//! \fn const void *RadiationVet::CleanDevCopy
//! \brief device copy (slot s < 8, n <= 8 KB) of an argument struct for non-inlined
//! device calls: read in place, not copied to every thread's stack as a capture is

const void *RadiationVet::CleanDevCopy(const void *p, const size_t n, const int s) {
  // (host builds too: the copy outlives the struct)
  constexpr size_t slot = 8192;
  if (n > slot || s < 0 || s >= 8) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "CleanDevCopy: " << n << " bytes in slot " << s
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (cln_xdev.extent(0) == 0) Kokkos::realloc(cln_xdev, 8*slot);
  Kokkos::View<const char *, HostMemSpace, Kokkos::MemoryUnmanaged>
      h(static_cast<const char *>(p), n);
  auto d = Kokkos::subview(cln_xdev, std::make_pair(s*slot, s*slot + n));
#if defined(KOKKOS_ENABLE_CUDA)
  // (CUDA: in stream order, no fence; a small pageable copy is staged at the call)
  Kokkos::deep_copy(DevExeSpace(), d, h);
#else
  Kokkos::deep_copy(d, h);
#endif
  return cln_xdev.data() + s*slot;
}

const CleanArgs *RadiationVet::CleanDevArgs(const CleanArgs &x, const int s) {
  return CleanDev(x, s);
}

//========================================================================================
// device memory arena (vet_gd intensities, clean work vectors over them)

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaInit
//! \brief enables the arena with vet_gd; allocates iw (not in the arena: some of its
//! channels live across a VetGdBuild)

void RadiationVet::MemArenaInit(ParameterInput *pin, const int nmb, const int niw,
                                const int n3, const int n2, const int n1) {
  marena_on = vgd_on;
  Kokkos::realloc(iw, nmb, niw, n3, n2, n1);
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::MemArenaVgd
//! \brief binds a.d (n Reals) into the arena after the previous intensity array;
//! false: not bound, the caller allocates a.d

bool RadiationVet::MemArenaVgd(VgdRag &a, const int64_t n) {
  if (!marena_on || marena_sealed) return false;
  const int64_t off = marena_sw;
  const int64_t need = off + n;
  MemArenaGrow(need, false);   // (the earlier array holds zeros, as a new arena)
  a.d = DvceArray1D<Real>(marena.data() + off, n);
  marena_sw = need + ((64 - (need % 64)) % 64);   // 64-Real alignment
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaGrow
//! \brief grows the arena to need Reals and rebinds the views in it (same offsets);
//! keep: copy the old contents, else free the old allocation first, nothing copied

void RadiationVet::MemArenaGrow(const int64_t need, const bool keep) {
  const int64_t have = static_cast<int64_t>(marena.extent(0));
  if (need <= have) return;
  Kokkos::fence();
  const Real *o0 = marena.data(), *o1 = marena.data() + have;
  if (!keep) marena = DvceArray1D<Real>();
  DvceArray1D<Real> nb("vet_arena", need);
  if (have > 0 && keep) {
    Kokkos::deep_copy(Kokkos::subview(nb, std::make_pair(int64_t(0), have)), marena);
  }
  VgdRag *bound[2] = {&vgd_i, &vgd_i_xb};
  for (VgdRag *b : bound) {
    if (have > 0 && b->d.data() != nullptr && b->d.data() >= o0 && b->d.data() < o1) {
      b->d = DvceArray1D<Real>(nb.data() + (b->d.data() - o0), b->d.extent(0));
    }
  }
  for (DvceArray5D<Real> *w : {&cln_w.o, &cln_w.s}) {   // bound by CleanWoBind
    if (have > 0 && w->data() != nullptr && w->data() >= o0 && w->data() < o1) {
      *w = DvceArray5D<Real>(nb.data() + (w->data() - o0), w->extent(0), w->extent(1),
                             w->extent(2), w->extent(3), w->extent(4));
    }
  }
  marena = nb;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaSeal
//! \brief called by every DtRetry::Save before it copies; the first call after the
//! intensities exist moves cln_w.o into the arena (copied, rebound), then seals it

void RadiationVet::MemArenaSeal() {
  // each checkpoint takes the handles anew, so a retry after the move restores them
  if (marena_sealed || !marena_on || marena.data() == nullptr) return;
  marena_sealed = true;
  DvceArray5D<Real> *cand[] = {&cln_w.o, &ost, &fmm_t, &odc, &d3w, &aa_xc, &aa_fc,
                               &aa_xp, &aa_fp};
  const char *nm[] = {"cln_wo", "ost", "fmm_t", "odc", "d3w", "aa_xc", "aa_fc", "aa_xp",
                      "aa_fp"};
  // grow to hold cln_w.o where it outgrows the intensities
  MemArenaGrow(static_cast<int64_t>(cln_w.o.span()));
  const int64_t have = static_cast<int64_t>(marena.extent(0));
  int64_t off = 0;
  std::string moved;
  Kokkos::fence();
  for (int q = 0; q < static_cast<int>(sizeof(cand)/sizeof(cand[0])); ++q) {
    DvceArray5D<Real> &a = *cand[q];
    const int64_t n = static_cast<int64_t>(a.span());
    // only cln_w.o (q = 0) moves
    if (q > 0 || a.data() == nullptr || n == 0 ||
        !a.span_is_contiguous() || off + n > have || a.data() == marena.data()) {
      // already in the arena (CleanWoBind)
      if (a.data() == marena.data() && n > 0) off += n + ((64 - (n % 64)) % 64);
      continue;
    }
    DvceArray5D<Real> b(marena.data() + off, a.extent(0), a.extent(1), a.extent(2),
                        a.extent(3), a.extent(4));
    Kokkos::deep_copy(b, a);
    a = b;
    off += n + ((64 - (n % 64)) % 64);
    moved += std::string(" ") + nm[q];
  }
  MemArenaBind();
  Kokkos::fence();
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> mem_arena: intensities " << marena_sw*sizeof(Real)/1.0e9
              << " GB in one allocation; solve scratch moved over them (" << off*
              sizeof(Real)/1.0e9 << " GB on rank 0):" << moved << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaBind
//! \brief once cln_w.o lies in the arena: moves the scratch arrays that are dead across a
//! closure build, written before they are read in each use, over the intensities too:
//! after cln_w the solver's (matrix rows, lateral sums, multigrid levels, write-back
//! scratch), and over both, from the start, the hydro's flux-stage arrays (dead across
//! the radiation stage: x2, x3 fluxes, FOFC's trial state, the well-balanced cache, the
//! seam, block-end and dfloor_flux_limit records).  The history reads uflx.x1f: kept.
//! Behind both: DualEnergy's pre-conversion state, used within Hydro::ConToPrim, and K2
//! (t2k2, t2pf2)

void RadiationVet::MemArenaBind() {
  const Real *a0 = marena.data();
  const int64_t have = static_cast<int64_t>(marena.extent(0));
  if (marena_bound || a0 == nullptr || cln_w.o.data() != a0) return;
  marena_bound = true;
  Kokkos::fence();
  auto align = [](const int64_t n) {return n + ((64 - (n % 64)) % 64);};
  int64_t off = 0;
  for (auto *w : {&cln_w.o, &cln_w.s}) {
    const int64_t o = w->data() - a0;
    if (w->data() != nullptr && o >= 0 && o < have) {
      off = std::max(off, align(o + static_cast<int64_t>(w->span())));
    }
  }
  int64_t moved = 0;
  auto bind = [&](DvceArray5D<Real> &a) {
    const int64_t n = static_cast<int64_t>(a.span());
    if (a.data() == nullptr || n <= 1 || off + n > have) return;
    a = DvceArray5D<Real>(marena.data() + off, a.extent(0), a.extent(1), a.extent(2),
                          a.extent(3), a.extent(4));
    off = align(off + n);
    moved += n;
  };
  // (cubed sphere without the live ring: cln_st holds the build's face-layer chunks,
  // CleanGdFaceSums)
  if (!pmy_pack->pmesh->use_cubed_sphere || cln_ring == 1) bind(cln_st);
  bind(cln_yl);
  bind(cln_vf);
  if (CleanExtrapSx()) bind(cln_csy);   // (else CleanRingMoments' scratch)
  for (auto &a : cln_mg) bind(a);
  // (the closure built once per stage: no pass keeps the faces' coefficients over one)
  if (!cln_closit) {
    for (CleanF *f : {&cln_fc1, &cln_fc2, &cln_fc3}) bind(f->o);
  }
  // (the write-back's: their first use makes them, here already over the arena)
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int n1 = id.nx1 + 2*id.ng, n2 = (id.nx2 > 1) ? id.nx2 + 2*id.ng : 1;
  const int n3 = (id.nx3 > 1) ? id.nx3 + 2*id.ng : 1, nmb = pmy_pack->nmb_thispack;
  for (const auto &s : {std::make_pair("cln_wbfp", 4), std::make_pair("cln_xf", 3)}) {
    const int64_t n = static_cast<int64_t>(nmb)*s.second*n3*n2*n1;
    if (off + n > have) continue;
    ScratchBind(s.first, DvceArray5D<Real>(marena.data() + off, nmb, s.second, n3, n2,
                                           n1));
    off = align(off + n);
    moved += n;
  }
  hydro::Hydro *ph = pmy_pack->phydro;
  if (ph != nullptr) {
    const int64_t os = off;
    off = 0;
    if (ph->rsolver_method != Hydro_RSolver::advect) {   // (it leaves IEN's fluxes zero)
      bind(ph->uflx.x2f);
      bind(ph->uflx.x3f);
    }
    marena_hyd = off;
    bind(ph->wbq0);
    bind(ph->utest);
    auto record = [&](DvceFaceFld5D<Real> *r) {
      if (r == nullptr) return;
      bind(r->x2f);
      bind(r->x3f);
    };
    record(ph->pseamfs);
    if (ph->pmlfs != ph->pseamfs) record(ph->pmlfs);
    if (ph->use_dfl) ph->DfloorFluxInit();
    record(ph->pdfl_fs);
    // (behind both: the conserved state before each conversion to primitives)
    off = std::max(off, os);
    bind(ph->dual.pre);
  }
  // K2 and its face fluxes: from the write-back of stage 1 to the end of the step; with
  // extrap1 no closure is built between them (sdirk2: the stage-2 solve writes K1)
  if (cln_extrap == 3 && t2_tab == 2) {
    bind(t2k2);
    bind(t2pf2);
    marena_k2 = true;
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> mem_arena: " << moved*sizeof(Real)/1.0e9
              << " GB of solver and hydro scratch over the intensities" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaHydro
//! \brief the hydro's fluxes over the intensities, zero at each flux stage as their own
//! allocations are where nothing writes them (FOFC's trial reads faces beyond the
//! sweeps), and the well-balanced cache's outer radial ghosts

void RadiationVet::MemArenaHydro() {
  if (marena_hyd == 0) return;
  auto h = Kokkos::subview(marena, std::make_pair(int64_t(0), marena_hyd));
  Kokkos::deep_copy(DevExeSpace(), h, 0.0);
  // (the x1 sweep under use_fofc reads the well-balanced cache one cell beyond the cells
  // BuildWBCache fills, is-1..ie+1: a zero background there, plain PLM)
  auto c = pmy_pack->phydro->wbq0;
  if (c.data() == nullptr || c.extent(1) < 15 || c.data() < marena.data() ||
      c.data() >= marena.data() + marena.extent(0)) return;
  const int is = pmy_pack->pmesh->mb_indcs.is, ie = pmy_pack->pmesh->mb_indcs.ie;
  const int nx = is - 1;   // (cells per side)
  par_for("vet_arena_wbq0", DevExeSpace(), 0, c.extent_int(0)-1, 0, 14, 0,
          c.extent_int(2)-1, 0, c.extent_int(3)-1, 0, 2*nx-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int q) {
    c(m,n,k,j,(q < nx) ? q : ie + 2 + q - nx) = 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanWoBind
//! \brief allocates cln_w.o at the first solve: over the intensities in the arena when
//! it exists (dead outside a build), else its own allocation.  cln_w.s, if CleanInit
//! left it unmade, follows cln_w.o (after it in the arena, or its own allocation)

void RadiationVet::CleanWoBind() {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int n1 = id.nx1 + 2*id.ng, n2 = (id.nx2 > 1) ? id.nx2 + 2*id.ng : 1;
  const int n3 = (id.nx3 > 1) ? id.nx3 + 2*id.ng : 1, nmb = pmy_pack->nmb_thispack;
  if (cln_w.o.data() == nullptr) {
    if (marena.data() == nullptr) {
      Kokkos::realloc(cln_w.o, nmb, cln_wo_nc, n3, n2, n1);
    } else {
      MemArenaGrow(static_cast<int64_t>(nmb)*cln_wo_nc*n3*n2*n1, false);
      cln_w.o = DvceArray5D<Real>(marena.data(), nmb, cln_wo_nc, n3, n2, n1);
    }
  }
  if (cln_w.s.data() != nullptr) {
    MemArenaBind();
    return;
  }
  if (cln_w.o.data() != marena.data()) {
    Kokkos::realloc(cln_w.s, nmb, CLW_NS, n3, n2, n1);
    return;
  }
  const int64_t no = static_cast<int64_t>(cln_w.o.span());
  const int64_t off = no + ((64 - (no % 64)) % 64);   // 64-Real alignment
  MemArenaGrow(off + static_cast<int64_t>(nmb)*CLW_NS*n3*n2*n1, false);
  cln_w.s = DvceArray5D<Real>(marena.data() + off, nmb, CLW_NS, n3, n2, n1);
  MemArenaBind();
}

}  // namespace radvet
