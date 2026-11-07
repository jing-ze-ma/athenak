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
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/rad_m1_parfor.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::ImplicitFbKeys: the keys (from ImplicitInit)
void RadiationM1::ImplicitFbKeys(ParameterInput *pin) {
  if (pin->DoesParameterExist("rad_m1", "implicit_fb_diag")) {
    fbd_n = pin->GetInteger("rad_m1", "implicit_fb_diag");
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
  par_for("m1_fb_ttz", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (vk_(m,k,j,i) > 0.0) {
      tl_(m,M1_TT_LAT0+1,k,j,i) = 0.0;
      tl_(m,M1_TT_LAT0+2,k,j,i) = 0.0;
      if (c3) {tl_(m,M1_TT_LAT0+3,k,j,i) = 0.0;}
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
  if (emin > 0.0 || fbd_done >= fbd_n) {return;}
  fbd_done += 1;
  // the records: host copies of what they read, on a recorded trigger only
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  auto s2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
             Kokkos::subview(iw, Kokkos::ALL, static_cast<int>(M1_IW_S2), Kokkos::ALL,
                             Kokkos::ALL, Kokkos::ALL));
  auto ep = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
             Kokkos::subview(iw, Kokkos::ALL, static_cast<int>(M1_IW_EP), Kokkos::ALL,
                             Kokkos::ALL, Kokkos::ALL));
  auto chi = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
             Kokkos::subview(opac, Kokkos::ALL, static_cast<int>(M1_OP_T), Kokkos::ALL,
                             Kokkos::ALL, Kokkos::ALL));
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
