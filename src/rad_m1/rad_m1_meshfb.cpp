//========================================================================================
// meshfb-1106: the mesh-wide positivity fallback of the implicit solve, measured and made
// local.  Both keys are read only when named; without them nothing here runs.
//
// <rad_m1>/implicit_fb_diag = N (DIAGNOSTIC): the first N passes that solved to E <= 0
//   write every such cell (cap 5000 per rank and pass) to fbdiag.rank<r>.txt in the run
//   directory and one line with the global min E cell to stdout; solve and pass
//   statistics of ALL solves (dropped or not) print every 2000 solves and at the end.
//   It adds a reduction per pass and moves no state: the run itself is unchanged.
// <rad_m1>/implicit_fb_local = true (default rounds 2, radius 2): a pass that solves to
//   E <= 0 drops D_r,lat and the implicit_vimp kick only in the L1 ball of the radius
//   implicit_fb_local_radius around the bad cells; the next pass assembles without them
//   there.  Rounds per solve: implicit_fb_local_rounds; once they are used, or if the
//   solve runs implicit_offdiag = operator (no local form), the global path fires.
//
// WHY THE L1 BALL OF RADIUS 2.  A row of the 19-point stencil reads cells with
// |di|+|dj|+|dk| <= 2 and D_r,lat at those (M1CsLat0: q's lateral neighbours at the x1
// faces; M1CsLatF: the radial neighbours at the lateral faces), and the vimp row reaches
// 2 along each axis.  With D_r,lat = 0 and no kick in that ball the row of a bad cell is
// the M-matrix row of `none`.  It is grown by repeated ONE-cell axis dilations with a
// ghost exchange between them, so it reads face ghosts only: no corner ghost (a cube
// vertex has none), and a block edge sees what a block interior sees (layout invariant).
//
// WHY A COPY OF tau_ten.  The solve reads D_r,lat through VlatTT(): tau_ten, or, after a
// local trigger, vlat_tt_loc = tau_ten with LAT1..3 zero in the mask and its ghosts
// exchanged (VetLatTTGhosts).  So the operator, its stencil, the assembly and the face
// update of one rank see the same field as their neighbours across ranks and seams, and
// tau_ten itself (the time2 extrapolation, the restart state, the dumps) never changes.
// The D_rr fold in slot 0 stays, as on the global path: it is a diagonal term.
//========================================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include <array>
#include <vector>
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "utils/recipe.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/rad_m1_parfor.hpp"
#include "rad_m1/m1_fluid.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::ImplicitFbKeys: the keys (from ImplicitInit)
void RadiationM1::ImplicitFbKeys(ParameterInput *pin) {
  if (pin->DoesParameterExist("rad_m1", "implicit_fb_diag")) {
    fbd_n = pin->GetInteger("rad_m1", "implicit_fb_diag");
  }
  if (pin->DoesParameterExist("rad_m1", "implicit_seam_dump")) {   // seamE-1106
    sdm_n = pin->GetInteger("rad_m1", "implicit_seam_dump");
  }
  if (pin->DoesParameterExist("rad_m1", "implicit_seam_dump_stage")) {   // seamE-1106
    sdm_stage = pin->GetBoolean("rad_m1", "implicit_seam_dump_stage");
  }
  // prod-1108 (USER 11-08 "make some defaults", r8 key set): implicit_fb_local DEFAULT
  // true with implicit_fb_local_rounds DEFAULT 4 (were false / 2; recipe none keeps them)
  if (recipe::GetMode(pin) != recipe::Mode::none) {
    fbl_on = true;
    fbl_rounds = 4;
  }
  if (pin->DoesParameterExist("rad_m1", "implicit_fb_local")) {
    fbl_on = pin->GetBoolean("rad_m1", "implicit_fb_local");
  }
  if (fbl_on) {
    if (pin->DoesParameterExist("rad_m1", "implicit_fb_local_rounds")) {
      fbl_rounds = pin->GetInteger("rad_m1", "implicit_fb_local_rounds");
    }
    if (pin->DoesParameterExist("rad_m1", "implicit_fb_local_radius")) {
      fbl_rad = pin->GetInteger("rad_m1", "implicit_fb_local_radius");
    }
    if (fbl_rad < 1 || fbl_rad > 8) {
      std::cout << "### FATAL ERROR <rad_m1>/implicit_fb_local_radius must be 1..8"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    vimp_kill = true;   // the kick's mask, ImplicitVimpBuild
  }
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::ImplicitFbMin: the global min E of this pass and the count <= 0
Real RadiationM1::ImplicitFbMin(Real *nbad) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  Real emin = 1.0e300, nb = 0.0;
  par_reduce_flat("m1_fb_min", 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lm) {
    const Real r = iw_(m,M1_IW_S2,k,j,i);
    lm = (r < lm) ? r : lm;
  }, Kokkos::Min<Real>(emin));
  par_reduce_flat("m1_fb_nb", 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a) {
    if (!(iw_(m,M1_IW_S2,k,j,i) > 0.0)) {a += 1.0;}
  }, Kokkos::Sum<Real>(nb));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &emin, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &nb, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  *nbad = nb;
  return emin;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::VetLatFbFill: vlat_tt_loc = tau_ten without LAT1..3 in the mask
void RadiationM1::VetLatFbFill() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  if (vlat_tt_loc.extent_int(0) == 0) {
    Kokkos::realloc(vlat_tt_loc, tau_ten.extent(0), tau_ten.extent(1), tau_ten.extent(2),
                    tau_ten.extent(3), tau_ten.extent(4));
  }
  Kokkos::deep_copy(DevExeSpace(), vlat_tt_loc, tau_ten);
  auto tl_ = vlat_tt_loc;
  auto vk_ = vkill_;
  const bool c3 = cs_geom;
  const bool cst = vgd_cstan;
  par_for("m1_fb_ttz", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (vk_(m,k,j,i) > 0.0) {
      tl_(m,M1_TT_LAT0+1,k,j,i) = 0.0;
      tl_(m,M1_TT_LAT0+2,k,j,i) = 0.0;
      if (c3) {tl_(m,M1_TT_LAT0+3,k,j,i) = 0.0;}
      if (cst) {   // cs-tan-1106: Q off in the mask too
        tl_(m,M1_TT_LAT0+4,k,j,i) = 0.0;
        tl_(m,M1_TT_LAT0+5,k,j,i) = 0.0;
      }
    }
  });
  VetLatTTGhosts(&vlat_tt_loc);
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::ImplicitFbLocal
//! \brief true = this pass's E <= 0 is handled by the local mask (the caller skips the
//! global fallback and takes another pass); false = no trigger, or the global path
bool RadiationM1::ImplicitFbLocal(int it, int odm) {
  if (!fbl_on) {return false;}
  Real nb = 0.0;
  const Real emin = ImplicitFbMin(&nb);
  if (emin > 0.0) {return false;}
  // csmetric-1112 vet_col_lat_posfix = lag: the first trigger of a step takes the global
  // path, which moves the lateral terms to their lagged form (nothing masked); a later
  // trigger masks as before
  if (vlat_posfix_lag && VlatOp()) {return false;}
  fbl_ntrig += 1.0;
  if (!vlat_on || !vlat_ready || (odm == M1_OD_OPERATOR) || (fbl_used >= fbl_rounds)) {
    fbl_nglobal += 1.0;
    return false;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto vk_ = vkill_;
  auto lx_ = vlat_lx;
  const int n3 = static_cast<int>(vk_.extent(1)) - 1;
  const int n2 = static_cast<int>(vk_.extent(2)) - 1;
  const int n1 = static_cast<int>(vk_.extent(3)) - 1;
  par_for("m1_fb_seed", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (!(iw_(m,M1_IW_S2,k,j,i) > 0.0)) {vk_(m,k,j,i) = 1.0;}
  });
  const bool thrd = trans_x3;
  const bool twod = pmy_pack->pmesh->multi_d;
  for (int r = 0; r < fbl_rad; ++r) {
    // the mask's face ghosts through the scratch of VetLatTTGhosts (slot 0)
    par_for("m1_fb_lxo", DevExeSpace(), 0, nmb1, 0, n3, 0, n2, 0, n1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      lx_(m,0,k,j,i) = vk_(m,k,j,i);
    });
    VetLatExchange(vlat_lx, vlat_tt_c, pbval_vt);
    // one cell along each axis, into lx slot 1 (out of place), then back
    par_for("m1_fb_dil", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      bool a = (lx_(m,0,k,j,i) != 0.0) || (lx_(m,0,k,j,i-1) != 0.0) ||
               (lx_(m,0,k,j,i+1) != 0.0);
      if (twod) {a = a || (lx_(m,0,k,j-1,i) != 0.0) || (lx_(m,0,k,j+1,i) != 0.0);}
      if (thrd) {a = a || (lx_(m,0,k-1,j,i) != 0.0) || (lx_(m,0,k+1,j,i) != 0.0);}
      vk_(m,k,j,i) = a ? 1.0 : 0.0;
    });
  }
  Real nc = 0.0;
  par_reduce_flat("m1_fb_nc", 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a) {
    if (vk_(m,k,j,i) > 0.0) {a += 1.0;}
  }, Kokkos::Sum<Real>(nc));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &nc, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  vkill_any = true;
  fbl_any = true;
  fbl_used += 1;
  fbl_nlocal += 1.0;
  fbl_ncell += nc;
  fbl_nbad += nb;
  VetLatFbFill();
  if (vimp_now) {vimp_emin = std::min(vimp_emin, emin);}
  od_emin = std::min(od_emin, emin);
  (void) it;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::ImplicitFbDiagPass: after the linear solve of a pass, BEFORE any
//! fallback acts: the pass's statistics, and the records of a pass with E <= 0
void RadiationM1::ImplicitFbDiagPass(int it) {
  const int kd = t2_inprobe ? 3 : ((t2_solve == M1_T2S_STAGE1) ? 1
                 : ((t2_solve == M1_T2S_STAGE2) ? 2 : 0));
  fbd_kind = kd;
  const bool lat = VlatOp(), vim = vimp_now;
  fbd_done = (tau_ten.extent_int(1) < M1_TT_LAT0 + 4) ? fbd_n : fbd_done;   // no LAT slots
  Real nb = 0.0;
  const Real emin = ImplicitFbMin(&nb);
  Real *s = fbd_s[kd];
  s[1] += 1.0;
  if (lat) {s[2] += 1.0;}
  if (vim) {s[3] += 1.0;}
  if (lat && vim && !fbl_any) {s[4] += 1.0;}
  fbd_last_full = lat && vim && !fbl_any;
  if (!(emin > 0.0)) {
    if (lat || vim) {
      s[5] += 1.0;
    } else {
      s[6] += 1.0;   // E <= 0 with both terms gone: they were not the cause
    }
    s[9] += nb;
    if (fbd_first < 0) {
      fbd_first = it;
      s[7] += 1.0;
      s[8] += static_cast<Real>(it);
    }
  }
  if (!(emin > 0.0) && sdm_done < sdm_n && (!sdm_stage || kd == 1 || kd == 2)) {   // seamE-1106
    sdm_done += 1;
    ImplicitSeamRowDump(it, kd, nb);
  }
  if (emin > 0.0 || fbd_done >= fbd_n) {return;}
  fbd_done += 1;
  // the records: host copies of what they read, on a recorded trigger only
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  // (a subview of a component is strided: no device -> host deep_copy; go through a
  // contiguous device array)
  auto &iwv = iw;
  auto &opv = opac;
  const int nh3 = iw.extent_int(2), nh2 = iw.extent_int(3), nh1 = iw.extent_int(4);
  DvceArray4D<Real> dt3("m1_fb_d3", nmb, nh3, nh2, nh1);
  auto fetch = [&](const DvceArray5D<Real> &src, const int c) {
    auto sv = src;
    auto dv = dt3;
    par_for("m1_fb_fetch", DevExeSpace(), 0, nmb - 1, 0, nh3 - 1, 0, nh2 - 1, 0, nh1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      dv(m,k,j,i) = sv(m,c,k,j,i);
    });
    Kokkos::fence();
    auto hv = Kokkos::create_mirror(HostMemSpace(), dv);   // rsgdiag-1107: a new array also on CPU (mirror_view aliased dt3: every fetch read the last one)
    Kokkos::deep_copy(hv, dv);
    return hv;
  };
  auto s2 = fetch(iwv, M1_IW_S2);
  auto ep = fetch(iwv, M1_IW_EP);
  auto chi = fetch(opv, M1_OP_T);
  auto tth = Kokkos::create_mirror_view_and_copy(HostMemSpace(), tau_ten);
  auto &size = pmy_pack->pmb->mb_size;
  auto x1fv = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->xx1f);
  auto lloc = pmy_pack->pmesh->lloc_eachmb;
  const int gs = pmy_pack->gids;
  char fn[64];
  std::snprintf(fn, sizeof(fn), "fbdiag.rank%05d.txt", global_variable::my_rank);
  FILE *fp = std::fopen(fn, "a");
  // the rank's min E record: E, panel, ig, jg, kg, r, tau_top, r/R_ph
  double mn[8] = {1.0e300, -1, -1, -1, -1, 0, 0, 0};
  const int cyc = pmy_pack->pmesh->ncycle;
  int nw = 0;
  (void) size;
  for (int m = 0; m < nmb; ++m) {
    const LogicalLocation &ll = lloc[gs + m];
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        bool any = false;
        for (int i = is; i <= ie; ++i) {any = any || !(s2(m,k,j,i) > 0.0);}
        if (!any) {continue;}
        // the column in this block: tau_top(i) and the radius of tau_top = 1 (a radially
        // split layout would stop at the block's top: all production decks have 1 block)
        std::vector<double> tt(ie - is + 2, 0.0);
        double rph = 0.0;
        for (int i = ie; i >= is; --i) {
          const double dr = x1fv(m,i+1) - x1fv(m,i);
          tt[i-is] = tt[i-is+1] + chi(m,k,j,i)*dr;
          if (rph == 0.0 && tt[i-is] >= 1.0) {rph = 0.5*(x1fv(m,i) + x1fv(m,i+1));}
        }
        for (int i = is; i <= ie; ++i) {
          if (s2(m,k,j,i) > 0.0) {continue;}
          const double r = 0.5*(x1fv(m,i) + x1fv(m,i+1));
          const double dr = x1fv(m,i+1) - x1fv(m,i);
          const double lmax = std::max(std::max(std::fabs(tth(m,M1_TT_LAT0+1,k,j,i)),
                                                std::fabs(tth(m,M1_TT_LAT0+2,k,j,i))),
                                       std::fabs(tth(m,M1_TT_LAT0+3,k,j,i)));
          const int ig = ll.lx1*indcs.nx1 + i - is, jg = ll.lx2*indcs.nx2 + j - js;
          const int kg = ll.lx3*indcs.nx3 + k - ks;
          const double rr = (rph > 0.0) ? r/rph : 0.0;
          if (s2(m,k,j,i) < mn[0]) {
            mn[0] = s2(m,k,j,i); mn[1] = ll.panel; mn[2] = ig; mn[3] = jg; mn[4] = kg;
            mn[5] = r; mn[6] = tt[i-is]; mn[7] = rr;
          }
          if (fp != nullptr && nw < 5000) {
            nw += 1;
            std::fprintf(fp, "%d %d %d %.0f %d %d %d %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e "
                         "%.6e %d %d\n", cyc, it, kd, nb, ll.panel, ig, jg, kg, r, rr,
                         tt[i-is], chi(m,k,j,i)*dr, static_cast<double>(s2(m,k,j,i)),
                         static_cast<double>(ep(m,k,j,i)), lmax,
                         static_cast<double>(tth(m,0,k,j,i)), lat ? 1 : 0, vim ? 1 : 0);
          }
        }
      }
    }
  }
  if (fp != nullptr) {std::fclose(fp);}
#if MPI_PARALLEL_ENABLED
  {
    std::vector<double> all(8*global_variable::nranks);
    MPI_Allgather(mn, 8, MPI_DOUBLE, all.data(), 8, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int q = 0; q < global_variable::nranks; ++q) {
      if (all[8*q] < mn[0]) {std::copy(&all[8*q], &all[8*q] + 8, mn);}
    }
  }
#endif
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> FBDIAG cycle " << cyc << " pass " << it << " kind "
              << ((kd == 1) ? "stage1" : ((kd == 2) ? "stage2" : ((kd == 3) ? "probe" : "be")))
              << " nE<=0 " << nb << " min E " << mn[0] << " at panel " << mn[1]
              << " i,j,k " << mn[2] << " " << mn[3] << " " << mn[4] << " r " << mn[5]
              << " tau_top " << mn[6] << " r/R_ph " << mn[7] << " vlat " << lat
              << " vimp " << vim << " masked " << fbl_any << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::ImplicitSeamRowDump (seamE-1106, DIAGNOSTIC, key implicit_seam_dump = N)
//! The rank's min E cell of the first N passes with E <= 0: its assembled row, split into
//! the implicit part and what is on the right-hand side, the four lateral faces (two-point
//! implicit flux vs the flux the row was given; the difference is explicit), and the
//! neighbours.  Writes seamrow.rank<r>.txt.  Moves no state.
void RadiationM1::ImplicitSeamRowDump(int it, int kd, Real nb) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const int nh3 = iw.extent_int(2), nh2 = iw.extent_int(3), nh1 = iw.extent_int(4);
  DvceArray4D<Real> dt3("m1_sd_d3", nmb, nh3, nh2, nh1);
  auto fetch = [&](const DvceArray5D<Real> &src, const int c) {
    auto sv = src;
    auto dv = dt3;
    par_for("m1_sd_fetch", DevExeSpace(), 0, nmb - 1, 0, nh3 - 1, 0, nh2 - 1, 0, nh1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      dv(m,k,j,i) = sv(m,c,k,j,i);
    });
    Kokkos::fence();
    auto hv = Kokkos::create_mirror(HostMemSpace(), dv);   // rsgdiag-1107: a new array also on CPU (mirror_view aliased dt3: every fetch read the last one)
    Kokkos::deep_copy(hv, dv);
    return hv;
  };
  auto s2 = fetch(iw, M1_IW_S2);
  // the rank's min E cell
  int bm = -1, bk = 0, bj = 0, bi = 0;
  double bv = 1.0e300;
  for (int m = 0; m < nmb; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (s2(m,k,j,i) < bv) {bv = s2(m,k,j,i); bm = m; bk = k; bj = j; bi = i;}
        }
      }
    }
  }
  if (!(bv <= 0.0)) {return;}
  auto ep = fetch(iw, M1_IW_EP);
  auto en = fetch(iw, M1_IW_EN);
  auto kt = fetch(iw, M1_IW_KT);
  auto wchi = fetch(iw, M1_IW_WCHI);
  auto v1 = fetch(iw, M1_IW_V1);
  auto v2 = fetch(iw, M1_IW_V2);
  auto v3 = fetch(iw, M1_IW_V3);
  auto a2 = fetch(iw, M1_IW_A2);
  auto a3 = fetch(iw, M1_IW_A3);
  auto g0 = fetch(iw, M1_IW_G0);
  auto n1 = fetch(iw, M1_IW_N1);
  auto srcb = fetch(iw, M1_IW_SRCB);
  auto srcr = fetch(iw, M1_IW_SRCR);
  auto tdia = fetch(iw, M1_IW_TDIA);
  auto trhs = fetch(iw, M1_IW_TRHS);
  auto ta = fetch(iw, M1_IW_TA);
  auto tb = fetch(iw, M1_IW_TB);
  auto tc = fetch(iw, M1_IW_TC);
  auto tr = fetch(iw, M1_IW_TR);
  auto cjm = fetch(iw, M1_IW_CJM);
  auto cjp = fetch(iw, M1_IW_CJP);
  auto ckm = fetch(iw, M1_IW_CKM);
  auto ckp = fetch(iw, M1_IW_CKP);
  auto kb = fetch(iw, M1_IW_KB);
  auto ke_ = fetch(opac, M1_OP_E);
  auto egn = fetch(iw, M1_IW_EGN);
  auto f0x2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), f0x2);
  auto f0x3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), f0x3);
  auto f2nh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), f0x2n);
  auto f3nh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), f0x3n);
  auto f1nh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), f0x1n);
  auto volh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->volume);
  auto cdx = M1Dxf();
  auto dx2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cdx.x2f);
  auto dx3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cdx.x3f);
  auto ar2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->area.x2f);
  auto ar3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->area.x3f);
  auto sn2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                  pmy_pack->pcoord->sin_face_xi);
  auto sn3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                  pmy_pack->pcoord->sin_face_eta);
  auto x1vh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x1v);
  auto x1fh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->xx1f);
  auto g2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), csg2);
  auto g3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), csg3);
  auto csh = cs_seam.h_view;
  DvceArray5D<Real> ddv("m1_sd_dd", nmb, 7, nh3, nh2, nh1);
  {
    auto dd = ddv; auto iw_ = iw; auto vd_ = vet_cell; const bool dfull = vet_full;
    par_for("m1_sd_dd", DevExeSpace(), 0, nmb - 1, 0, nh3 - 1, 0, nh2 - 1, 0, nh1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      dd(m,0,k,j,i) = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
      dd(m,1,k,j,i) = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
    });
    // the lateral term of the x1 face equation at active cells, with and without the seam
    // foot -> midpoint shift
    if (tau_ten.extent_int(1) >= M1_TT_LAT0 + 4) {
      auto vlt_ = VlatTT();
      auto cvol = pmy_pack->pcoord->volume;
      auto carf = cs_area_eff;
      auto cfn2 = cs_fn2;
      auto cfn3 = cs_fn3;
      auto cse_ = cs_seam.d_view;
      auto cg2_ = csg2;
      auto cg3_ = csg3;
      const int c0l = M1_TT_LAT0;
      const bool thrd = trans_x3;
      par_for("m1_sd_lat", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        dd(m,2,k,j,i) = M1CsLat0(iw_,vlt_,c0l,cvol,carf,cfn2,cfn3,m,k,j,i,thrd,M1_IW_EP,true,
                                 cse_,cg2_,cg3_,js,je,ks,ke);
        dd(m,3,k,j,i) = M1CsLat0(iw_,vlt_,c0l,cvol,carf,cfn2,cfn3,m,k,j,i,thrd,M1_IW_EP,false,
                                 cse_,cg2_,cg3_,js,je,ks,ke);
      });
    }
    Kokkos::fence();
  }
  {
    // the gas state of the row: T_k, T(EGN), e(T_k), rho (the stage old vector), uh(IEN)
    FluidRef fl = FluidRef::Get(pmy_pack);
    if (fl.on) {
      auto dd = ddv; auto iw_ = iw; auto uh = fl.u0; auto eos = fl.eos;
      par_for("m1_sd_gas", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real r = uh(m,IDN,k,j,i), tk = iw_(m,M1_IW_TP,k,j,i);
        const Real eg = iw_(m,M1_IW_EGN,k,j,i);
        Real te, pp, cr, ct, tcv;
        eos.ThermoAt(r, tk, te, pp, cr, ct, tcv);
        dd(m,4,k,j,i) = te;
        dd(m,5,k,j,i) = (eg > 0.0) ? eos.Temperature(r, eg) : -1.0;
        dd(m,6,k,j,i) = r;
      });
      Kokkos::fence();
    }
  }
  auto ddh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ddv);
  auto lloc = pmy_pack->pmesh->lloc_eachmb;
  char fn[64];
  std::snprintf(fn, sizeof(fn), "seamrow.rank%05d.txt", global_variable::my_rank);
  FILE *fp = std::fopen(fn, "a");
  if (fp == nullptr) {return;}
  // rsgdiag-1107: every cell with E <= 0 of this rank (cap 40), not only the min
  std::vector<std::array<int,4>> cand;
  for (int mm = 0; mm < nmb; ++mm) {
    for (int kk = ks; kk <= ke; ++kk) {
      for (int jj = js; jj <= je; ++jj) {
        for (int ii = is; ii <= ie; ++ii) {
          if (!(s2(mm,kk,jj,ii) > 0.0) && cand.size() < 40) {cand.push_back({mm,kk,jj,ii});}
        }
      }
    }
  }
  for (const auto &cnd : cand) {
  const int m = cnd[0], k = cnd[1], j = cnd[2], i = cnd[3];
  const LogicalLocation &ll = lloc[pmy_pack->gids + m];
  const int ig = ll.lx1*indcs.nx1 + i - is, jg = ll.lx2*indcs.nx2 + j - js;
  const int kg = ll.lx3*indcs.nx3 + k - ks;
  const double cl = c_light, ch = chat, dt = dt_sub, cr = ch/cl;
  const double iv = dt/volh(m,k,j,i);
  std::fprintf(fp, "TRIG cyc %d it %d kind %d nb %.0f panel %d i %d j %d k %d seam %d %d %d %d "
               "vlat %d vimp %d dt %.5e ch %.5e cl %.5e r %.6e\n",
               static_cast<int>(pmy_pack->pmesh->ncycle), it, kd, nb, ll.panel, ig, jg, kg,
               csh(m,0), csh(m,1), csh(m,2), csh(m,3), VlatOp() ? 1 : 0, vimp_now ? 1 : 0, dt, ch,
               cl, static_cast<double>(x1vh(m,i)));
  const double rhs_lag = kb(m,k,j,i) - tr(m,k,j,i) + trhs(m,k,j,i);
  // the row, S2 reading ghosts when it has them (see row_res: only meaningful if S2's
  // ghosts are valid at this point)
  const double lhs = ta(m,k,j,i)*s2(m,k,j,i-1) + tb(m,k,j,i)*s2(m,k,j,i)
                     + tc(m,k,j,i)*s2(m,k,j,i+1) - cjm(m,k,j,i)*0.0;
  std::fprintf(fp, "ROW TA %.6e TB %.6e TC %.6e CJM %.6e CJP %.6e CKM %.6e CKP %.6e TR %.6e "
               "KB %.6e EN %.6e EP %.6e E %.6e SRCB %.6e SRCR %.6e TDIA %.6e TRHS %.6e "
               "latlag %.6e radrest %.6e radlhs %.6e\n",
               ta(m,k,j,i), tb(m,k,j,i), tc(m,k,j,i), cjm(m,k,j,i), cjp(m,k,j,i),
               ckm(m,k,j,i), ckp(m,k,j,i), tr(m,k,j,i), kb(m,k,j,i), en(m,k,j,i), ep(m,k,j,i),
               s2(m,k,j,i), srcb(m,k,j,i), srcr(m,k,j,i), tdia(m,k,j,i), trhs(m,k,j,i), rhs_lag,
               tr(m,k,j,i) - en(m,k,j,i) - srcr(m,k,j,i) - trhs(m,k,j,i), lhs);
  // the lateral faces: 2 = x2, 3 = x3; side 0 = m, 1 = p.  Fd: the two-point flux the
  // row holds implicitly; lag = dt A/V cr (F - Fd) its explicit rest, signed as it enters
  // the right-hand side.  Fmem: theta F^n, the old flux in F.
  auto face = [&](const int d, const int sd) {
    const int jn = (d == 2) ? (j + 2*sd - 1) : j, kn = (d == 3) ? (k + 2*sd - 1) : k;
    const int jf = (d == 2) ? (j + sd) : j, kf = (d == 3) ? (k + sd) : k;
    const bool seam = (d == 2) ? ((sd == 0 && j == js && csh(m,0)) ||
                                  (sd == 1 && j == je && csh(m,1)))
                               : ((sd == 0 && k == ks && csh(m,2)) ||
                                  (sd == 1 && k == ke && csh(m,3)));
    const double dsh = 0.5*(x1fh(m,i+1)*x1fh(m,i+1) - x1fh(m,i)*x1fh(m,i));
    double A, L;
    if (d == 2) {
      A = seam ? dsh*g2h(m,0,k,sd) : ar2h(m,kf,jf,i);
      L = seam ? x1vh(m,i)*g2h(m,1,k,sd) : dx2h(m,kf,jf,i)*sn2h(m,kf,jf);
    } else {
      A = seam ? dsh*g3h(m,0,j,sd) : ar3h(m,kf,jf,i);
      L = seam ? x1vh(m,i)*g3h(m,1,j,sd) : dx3h(m,kf,jf,i)*sn3h(m,kf,jf);
    }
    const int dq = d - 2;
    const double th = 1.0/(1.0 + ch*dt*0.5*(kt(m,k,j,i) + kt(m,kn,jn,i)));
    const double dc = ddh(m,dq,k,j,i), dn = ddh(m,dq,kn,jn,i);
    const double sgn = (sd == 1) ? 1.0 : -1.0;
    const double fd = -sgn*th*ch*cl*dt*(dn*ep(m,kn,jn,i) - dc*ep(m,k,j,i))/L;
    const double f = (d == 2) ? f0x2h(m,kf,jf,i) : f0x3h(m,kf,jf,i);
    const double fn_ = (d == 2) ? f2nh(m,kf,jf,i) : f3nh(m,kf,jf,i);
    const double lag = -sgn*iv*A*cr*(f - fd);
    std::fprintf(fp, "FACE x%d %s seam %d A/V*dt*cr %.5e L %.5e th %.5e g=c2dt2A/(VL) %.5e Dc "
                 "%.5e Dn %.5e En %.6e kt_n %.5e F %.6e Fd %.6e Fmem %.6e lag %.6e\n",
                 d, sd ? "p" : "m", seam ? 1 : 0, iv*A*cr, L, th, iv*A*th*ch*ch*dt/L, dc, dn,
                 ep(m,kn,jn,i), kt(m,kn,jn,i), f, fd, th*fn_, lag);
  };
  face(2,0); face(2,1); face(3,0); face(3,1);
  {
    auto ar1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->area.x1f);
    auto dx1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cdx.x1f);
    const double nup = ar1h(m,k,j,i+1)*iv, num = ar1h(m,k,j,i)*iv;
    const double thp = 1.0/(1.0 + ch*dt*0.5*(kt(m,k,j,i) + kt(m,k,j,i+1)));
    const double thm = 1.0/(1.0 + ch*dt*0.5*(kt(m,k,j,i) + kt(m,k,j,i-1)));
    std::fprintf(fp, "FACEX1 nup %.5e num %.5e thp %.5e thm %.5e dxp %.5e dxm %.5e gp %.5e gm "
                 "%.5e F0n_p %.6e F0n_m %.6e mem_rhs %.6e latp %.6e latm %.6e lat_rhs %.6e "
                 "latnc_rhs %.6e\n", nup, num, thp, thm, dx1h(m,k,j,i+1), dx1h(m,k,j,i),
                 nup*thp*ch*ch*dt/dx1h(m,k,j,i+1), num*thm*ch*ch*dt/dx1h(m,k,j,i),
                 f1nh(m,k,j,i+1), f1nh(m,k,j,i),
                 -nup*cr*thp*f1nh(m,k,j,i+1) + num*cr*thm*f1nh(m,k,j,i),
                 0.5*(ddh(m,2,k,j,i) + ddh(m,2,k,j,i+1)),
                 0.5*(ddh(m,2,k,j,i-1) + ddh(m,2,k,j,i)),
                 -nup*cr*thp*ch*cl*dt*0.5*(ddh(m,2,k,j,i) + ddh(m,2,k,j,i+1))
                 + num*cr*thm*ch*cl*dt*0.5*(ddh(m,2,k,j,i-1) + ddh(m,2,k,j,i)),
                 -nup*cr*thp*ch*cl*dt*0.5*(ddh(m,3,k,j,i) + ddh(m,3,k,j,i+1))
                 + num*cr*thm*ch*cl*dt*0.5*(ddh(m,3,k,j,i-1) + ddh(m,3,k,j,i)));
    if (vimp_now && iw_vimp >= 0) {
      auto jd = fetch(iw, iw_vimp + M1_IV_JD);
      auto j1m = fetch(iw, iw_vimp + M1_IV_J1M);
      auto j1p = fetch(iw, iw_vimp + M1_IV_J1P);
      auto jr = fetch(iw, iw_vimp + M1_IV_JRHS);
      std::fprintf(fp, "VIMP JD %.6e J1M %.6e J1P %.6e JRHS %.6e\n", jd(m,k,j,i), j1m(m,k,j,i),
                   j1p(m,k,j,i), jr(m,k,j,i));
    }
    {
      auto ifh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ifw);
      std::fprintf(fp, "IFW1 DGp %.6e DGm %.6e HCLp %.5e HCRp %.5e HCLm %.5e HCRm %.5e ALp %.4e ALm "
                   "%.4e dg_rhs %.6e\n", ifh(m,M1_IFW_DG,k,j,i+1), ifh(m,M1_IFW_DG,k,j,i),
                   ifh(m,M1_IFW_HCL,k,j,i+1), ifh(m,M1_IFW_HCR,k,j,i+1), ifh(m,M1_IFW_HCL,k,j,i),
                   ifh(m,M1_IFW_HCR,k,j,i), ifh(m,M1_IFW_AL,k,j,i+1), ifh(m,M1_IFW_AL,k,j,i),
                   -nup*ifh(m,M1_IFW_DG,k,j,i+1) + num*ifh(m,M1_IFW_DG,k,j,i));
    }
    if (kd == 1 || kd == 2) {
      auto ti = Kokkos::create_mirror_view_and_copy(HostMemSpace(), t2inc);
      auto k1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), t2k1);
      auto k2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), t2k2);
      auto tph = fetch(iw, M1_IW_TP);
      auto opp = fetch(opac, M1_OP_P);
      auto de0h = fetch(iw, M1_IW_DE0);
      std::fprintf(fp, "GAS Tk %.5e TEGN %.5e ek %.5e egn %.5e rho %.5e kP %.5e kE %.5e de0 %.5e "
                   "emisT4 %.5e\n", tph(m,k,j,i), ddh(m,5,k,j,i), ddh(m,4,k,j,i), egn(m,k,j,i),
                   ddh(m,6,k,j,i), opp(m,k,j,i), ke_(m,k,j,i), de0h(m,k,j,i),
                   dt*ch*opp(m,k,j,i)*arad*std::pow(tph(m,k,j,i), 4));
      for (const int q : {0, -1, 1}) {
        std::fprintf(fp, "T2 di %d incE %.6e incEN %.6e k1E %.6e k1EN %.6e k2E %.6e k2EN %.6e "
                     "incM1 %.5e\n", q, ti(m,M1_T2_E,k,j,i+q), ti(m,M1_T2_EN,k,j,i+q),
                     k1h(m,M1_T2_E,k,j,i+q), k1h(m,M1_T2_EN,k,j,i+q), k2h(m,M1_T2_E,k,j,i+q),
                     k2h(m,M1_T2_EN,k,j,i+q), ti(m,M1_T2_M1,k,j,i+q));
      }
    }
    auto tth = Kokkos::create_mirror_view_and_copy(HostMemSpace(), tau_ten);
    std::fprintf(fp, "LATC D0 %.5e L1 %.5e L2 %.5e L3 %.5e\n", tth(m,0,k,j,i),
                 tth(m,M1_TT_LAT0+1,k,j,i), tth(m,M1_TT_LAT0+2,k,j,i),
                 tth(m,M1_TT_LAT0+3,k,j,i));
  }
  // the stencil
  const int off[][3] = {{0,0,0},{-2,0,0},{-1,0,0},{1,0,0},{2,0,0},{0,-2,0},{0,-1,0},{0,1,0},
                        {0,2,0},{0,0,-2},{0,0,-1},{0,0,1},{0,0,2}};
  for (const auto &o : off) {
    const int ii = i + o[0], jj = j + o[1], kk = k + o[2];
    std::fprintf(fp, "CELL %d %d %d EN %.6e EP %.6e S2 %.6e KT %.5e KE %.5e Dx2 %.4e Dx3 %.4e "
                 "chi %.4e n1 %.4e v1 %.4e v2 %.4e v3 %.4e a2 %.4e a3 %.4e g0 %.4e egn %.5e "
                 "SRCB %.4e SRCR %.5e\n", o[0], o[1], o[2], en(m,kk,jj,ii), ep(m,kk,jj,ii),
                 s2(m,kk,jj,ii), kt(m,kk,jj,ii), ke_(m,kk,jj,ii), ddh(m,0,kk,jj,ii),
                 ddh(m,1,kk,jj,ii), wchi(m,kk,jj,ii), n1(m,kk,jj,ii), v1(m,kk,jj,ii),
                 v2(m,kk,jj,ii), v3(m,kk,jj,ii), a2(m,kk,jj,ii), a3(m,kk,jj,ii), g0(m,kk,jj,ii),
                 egn(m,kk,jj,ii), srcb(m,kk,jj,ii), srcr(m,kk,jj,ii));
  }
  }   // rsgdiag-1107: cand
  std::fclose(fp);
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::ImplicitFbDiagSolve: after the pass loop of a solve
void RadiationM1::ImplicitFbDiagSolve(int npass) {
  Real *s = fbd_s[fbd_kind];
  s[0] += 1.0;
  if (fbd_last_full) {s[10] += 1.0;}   // the accepted pass ran with both terms everywhere
  if (npass <= 1) {s[11] += 1.0;}
  fbd_first = -1;
  fbd_nsolve += 1;
  if ((fbd_nsolve % 2000) == 0) {ImplicitFbSummary();}
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::ImplicitFbSummary (rank 0, stdout)
void RadiationM1::ImplicitFbSummary() {
  if (global_variable::my_rank != 0) {return;}
  if (fbd_n > 0) {
    const char *nm[4] = {"be", "stage1", "stage2", "probe"};
    for (int q = 0; q < 4; ++q) {
      const Real *s = fbd_s[q];
      if (!(s[0] > 0.0)) {continue;}
      std::cout << "<rad_m1> FBDIAG_SUM " << nm[q] << " solves " << s[0] << " 1-pass "
                << s[11] << " passes " << s[1] << " vlat-active " << s[2]
                << " vimp-active " << s[3] << " both-unmasked " << s[4]
                << " | solves with E<=0 " << s[7] << " mean first pass "
                << ((s[7] > 0.0) ? s[8]/s[7] : 0.0) << " | E<=0 passes: terms on " << s[5]
                << " terms off " << s[6] << " bad cell-passes " << s[9]
                << " | solves whose last pass ran both terms everywhere " << s[10]
                << std::endl;
    }
  }
  if (fbl_on) {
    std::cout << "<rad_m1> FBLOCAL triggers " << fbl_ntrig << " local " << fbl_nlocal
              << " global " << fbl_nglobal << " bad cells " << fbl_nbad
              << " masked cells " << fbl_ncell << std::endl;
  }
}

} // namespace radm1
