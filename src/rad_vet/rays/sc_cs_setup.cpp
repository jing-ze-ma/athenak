//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cs_setup.cpp
//! \brief cubed-sphere rays: set-up, source band, tables.
//!
//!  called from: rays/sc_build.cpp, rays/sc_setup.cpp, rays/sc_tables.cpp
//!  calls into: rays/sc_cs_halo.cpp

#include <math.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_cs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

// build generation per compact-list cache entry; not a member, so the dt-retry
// checkpoint keeps it and DtRetryAfterRestore drops entries rebuilt by the failed try
int64_t g_vgd_csgen[16] = {0};        // VET_VGD_CSCL_MAX

int64_t g_vgd_csgen_save[16] = {0};

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsInit
//! \brief cs band setup (from VetGdInit): depth vgd_w, panel tables, unit vector of
//! every band cell, exact exchange lists

void RadiationVet::VetGdCsInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int g0 = pmy_pack->gids;
  const int me = global_variable::my_rank;
  const int nc = pm->mesh_indcs.nx2;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  if (nc < 2 || pm->mesh_indcs.nx3 != nc || (nc % nx2) != 0 || (nc % nx3) != 0) {
    VgdCsFatal("needs mesh nx2 = nx3 = N >= 2 per panel, divided by the MeshBlock sizes");
  }
  if (pm->multilevel) {VgdCsFatal("SMR/AMR not supported");}
  vgd_csn = nc;
  const double h = 0.5*M_PI/nc;
  const double q4 = 0.25*M_PI;
  const int nbj = nc/nx2, nbk = nc/nx3;
  const int ngb = pm->nmb_total;
  // blocks: panel and global first cell; a block's band comes from its own x1 level
  std::vector<int> bp(ngb), bj(ngb), bk(ngb), gof(vgd_nb1*6*nbj*nbk, -1);
  auto gofi = [&](const int g, const int p, const int b3, const int b2) -> int & {
    return gof[((pm->lloc_eachmb[g].lx1*6 + p)*nbk + b3)*nbj + b2];
  };
  for (int g = 0; g < ngb; ++g) {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    if (l.panel < 0 || l.panel > 5 || l.lx2 < 0 || l.lx2 >= nbj || l.lx3 < 0 ||
        l.lx3 >= nbk) {
      VgdCsFatal("unexpected MeshBlock location (panel / lx2 / lx3)");
    }
    bp[g] = l.panel;
    bj[g] = l.lx2*nx2;
    bk[g] = l.lx3*nx3;
    gofi(g, l.panel, l.lx3, l.lx2) = g;
  }
  for (size_t t = 0; t < gof.size(); ++t) {
    if (gof[t] < 0) {VgdCsFatal("incomplete MeshBlock grid");}
  }
  {
    // the local blocks' sizes agree with the logical locations
    auto &mbs = pmy_pack->pmb->mb_size;
    for (int m = 0; m < nmb; ++m) {
      const double x2 = -1.0 + 2.0*bj[g0 + m]/nc, x3 = -1.0 + 2.0*bk[g0 + m]/nc;
      if (std::fabs(mbs.h_view(m).x2min - x2) > 1.0e-12 ||
          std::fabs(mbs.h_view(m).x3min - x3) > 1.0e-12) {
        VgdCsFatal("MeshBlock size does not match its logical location");
      }
    }
  }
  // panel frames and the neighbour across each side (0 -x2, 1 +x2, 2 -x3, 3 +x3)
  double fa[6][3], fb[6][3], fn[6][3];
  for (int p = 0; p < 6; ++p) {
    Real a[3], b[3], n[3];
    cubed_sphere::PanelFrame(p, a, b, n);
    for (int e = 0; e < 3; ++e) {fa[p][e] = a[e]; fb[p][e] = b[e]; fn[p][e] = n[e];}
  }
  int nbr[6][4];
  for (int p = 0; p < 6; ++p) {
    for (int s = 0; s < 4; ++s) {
      double t[3];
      for (int e = 0; e < 3; ++e) {
        const double v = (s < 2) ? fa[p][e] : fb[p][e];
        t[e] = (s % 2 == 0) ? -v : v;
      }
      nbr[p][s] = -1;
      for (int q = 0; q < 6; ++q) {
        if (t[0]*fn[q][0] + t[1]*fn[q][1] + t[2]*fn[q][2] > 0.5) {nbr[p][s] = q;}
      }
      if (nbr[p][s] < 0) {VgdCsFatal("panel neighbour table");}
    }
  }
  auto cen = [&](const int p, const int a, const int t, double c[3]) {
    Real q[3];
    cubed_sphere::PanelToCart(p, -q4 + (a + 0.5)*h, -q4 + (t + 0.5)*h, q);
    for (int e = 0; e < 3; ++e) {c[e] = q[e];}
  };
  // unfolded position (J, K) of panel p -> the global cell (q, a, t) it holds
  auto unfold = [&](const int p, const int J, const int K, int &q, int &a,
                    int &t) -> bool {
    const bool ij = (J >= 0 && J < nc), ik = (K >= 0 && K < nc);
    if (ij && ik) {q = p; a = J; t = K; return true;}
    if (!ij && !ik) {return false;}
    int s, jm = J, km = K;
    if (!ij) {
      s = (J < 0) ? 0 : 1;
      jm = (J < 0) ? (-1 - J) : (2*nc - 1 - J);
    } else {
      s = (K < 0) ? 2 : 3;
      km = (K < 0) ? (-1 - K) : (2*nc - 1 - K);
    }
    if (jm < 0 || jm >= nc || km < 0 || km >= nc) {return false;}
    q = nbr[p][s];
    double c[3], e[3];
    cen(p, jm, km, c);
    for (int u = 0; u < 3; ++u) {e[u] = (fn[p][u] - fn[q][u])/std::sqrt(2.0);}
    const double dc = c[0]*e[0] + c[1]*e[1] + c[2]*e[2];
    Real cr[3];
    for (int u = 0; u < 3; ++u) {cr[u] = c[u] - 2.0*dc*e[u];}
    Real xi, eta;
    cubed_sphere::CartToPanel(q, cr, xi, eta);
    const double ga = (xi + q4)/h - 0.5, gt = (eta + q4)/h - 0.5;
    a = static_cast<int>(std::lround(ga));
    t = static_cast<int>(std::lround(gt));
    if (std::fabs(ga - a) > 1.0e-6 || std::fabs(gt - t) > 1.0e-6 || a < 0 || a >= nc ||
        t < 0 || t >= nc) {
      VgdCsFatal("unfold: the mirror cell is off the neighbour's lattice");
    }
    return true;
  };
  // inverse maps (q's (a, t) -> p's unfolded (J, K)): signed permutation + shift
  std::vector<int> cmh(6*6*7, 0);
  auto CM = [&](const int p, const int q, const int c) -> int & {
    return cmh[(p*6 + q)*7 + c];
  };
  for (int p = 0; p < 6; ++p) {
    CM(p,p,0) = 1; CM(p,p,1) = 1; CM(p,p,5) = 1;
    for (int s = 0; s < 4; ++s) {
      const int jb = (s == 0) ? -1 : ((s == 1) ? nc : 0);
      const int kb = (s == 2) ? -1 : ((s == 3) ? nc : 0);
      const int sj = (s == 0) ? -1 : 1, sk = (s == 2) ? -1 : 1;
      int q0, a0, t0, q1, a1, t1, q2, a2, t2;
      const bool o0 = unfold(p, jb, kb, q0, a0, t0);
      const bool o1 = unfold(p, jb + ((s < 2) ? sj : 1), kb, q1, a1, t1);
      const bool o2 = unfold(p, jb, kb + ((s < 2) ? 1 : sk), q2, a2, t2);
      if (!o0 || !o1 || !o2 || q1 != q0 || q2 != q0) {VgdCsFatal("unfold strip");}
      const int stj = (s < 2) ? sj : 1, stk = (s < 2) ? 1 : sk;
      const int uj0 = (a1 - a0)*stj, uj1 = (t1 - t0)*stj;   // (a, t) per unit J
      const int uk0 = (a2 - a0)*stk, uk1 = (t2 - t0)*stk;   // (a, t) per unit K
      if (std::abs(uj0) + std::abs(uj1) != 1 || std::abs(uk0) + std::abs(uk1) != 1 ||
          uj0*uk0 + uj1*uk1 != 0) {
        VgdCsFatal("unfold: not a signed permutation");
      }
      CM(p,q0,0) = 1;
      CM(p,q0,1) = uj0;
      CM(p,q0,2) = uj1;
      CM(p,q0,3) = jb - uj0*a0 - uj1*t0;
      CM(p,q0,4) = uk0;
      CM(p,q0,5) = uk1;
      CM(p,q0,6) = kb - uk0*a0 - uk1*t0;
    }
  }
  // check: forward unfold and inverse map agree on every strip cell
  for (int p = 0; p < 6; ++p) {
    for (int s = 0; s < 4; ++s) {
      for (int d = 0; d < nc; ++d) {
        for (int g = 0; g < nc; ++g) {
          const int J = (s == 0) ? (-1 - d) : ((s == 1) ? (nc + d) : g);
          const int K = (s == 2) ? (-1 - d) : ((s == 3) ? (nc + d) : g);
          int q, a, t;
          if (!unfold(p, J, K, q, a, t)) {VgdCsFatal("unfold check");}
          const int jj = CM(p,q,1)*a + CM(p,q,2)*t + CM(p,q,3);
          const int kk = CM(p,q,4)*a + CM(p,q,5)*t + CM(p,q,6);
          if (jj != J || kk != K) {VgdCsFatal("unfold: inverse map mismatch");}
        }
      }
    }
  }
  // band depth: sp formula with the narrowest angular spacing (+3: one more than sp for
  // the index kink across a seam); capped at N
  double dmin = 1.0e300;
  for (int a = 0; a < nc; ++a) {
    for (int t = 0; t < nc; ++t) {
      double c0[3], c1[3];
      cen(0, a, t, c0);
      if (a + 1 < nc) {
        cen(0, a + 1, t, c1);
        dmin = std::min(dmin, std::acos(std::min(1.0, c0[0]*c1[0] + c0[1]*c1[1] +
                                                         c0[2]*c1[2])));
      }
      if (t + 1 < nc) {
        cen(0, a, t + 1, c1);
        dmin = std::min(dmin, std::acos(std::min(1.0, c0[0]*c1[0] + c0[1]*c1[1] +
                                                         c0[2]*c1[2])));
      }
    }
  }
  {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->x1v);
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->xx1f);
    vgd_wsh.assign(indcs.nx1 + 2*indcs.ng, 0);
    vgd_wsi.assign(indcs.nx1 + 2*indcs.ng, 0);
    vgd_wso.assign(indcs.nx1 + 2*indcs.ng, 0);
    int wneed = 0;
    auto xf_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_x1f);
    for (int m = 0; m < nmb; ++m) {
      // shells next to another level: radii from the x1 ghosts
      const bool lo = (xf_h(m) & 1), hi = (xf_h(m) & 2);
      for (int l = indcs.is; l <= indcs.ie; ++l) {
        const double r = x1v_h(m,l);
        const double rlo = (l > indcs.is || !lo) ? x1v_h(m,l-1) : x1f_h(m,l);
        double a = 2.0*std::sqrt(std::max(r*r - rlo*rlo, 0.0));
        if (l < indcs.ie || !hi) {
          a = std::max(a, std::sqrt(x1v_h(m,l+1)*x1v_h(m,l+1) - r*r));
        }
        wneed = std::max(wneed, static_cast<int>(std::ceil(a/(rlo*dmin))) + 3);
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &wneed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    vgd_w = std::min(wneed, nc);
    vgd_capped = (vgd_w < wneed);
    for (int l = indcs.is; l <= indcs.ie; ++l) {
      vgd_wsh[l] = vgd_w;
      vgd_wsi[l] = vgd_w;
      vgd_wso[l] = vgd_w;
    }
    if (me == 0) {
      std::cout << "<rad_vet> vet_gd (cubed sphere): N = " << nc << " per panel, "
                << "narrowest "
                << "cell spacing " << dmin << " rad, dense lateral band " << vgd_w
                << " cells (needed " << wneed << ")"
                << (vgd_capped ? ": capped at N, near-tangent reads beyond are clamped "
                                 "(counted)" : "") << std::endl;
    }
  }
  const int w = vgd_w;
  const int c2w = nx2 + 2*w, c3w = nx3 + 2*w;
  // device tables
  Kokkos::realloc(vgd_cscm, 6, 6, 7);
  Kokkos::realloc(vgd_csnb, 6, 4);
  Kokkos::realloc(vgd_csblk, nmb, 3);
  {
    auto cm_h = Kokkos::create_mirror_view(vgd_cscm);
    auto nb_h = Kokkos::create_mirror_view(vgd_csnb);
    auto bl_h = Kokkos::create_mirror_view(vgd_csblk);
    for (int p = 0; p < 6; ++p) {
      for (int q = 0; q < 6; ++q) {
        for (int c = 0; c < 7; ++c) {cm_h(p,q,c) = CM(p,q,c);}
      }
      for (int s = 0; s < 4; ++s) {nb_h(p,s) = nbr[p][s];}
    }
    for (int m = 0; m < nmb; ++m) {
      bl_h(m,0) = bp[g0 + m];
      bl_h(m,1) = bj[g0 + m];
      bl_h(m,2) = bk[g0 + m];
    }
    Kokkos::deep_copy(vgd_cscm, cm_h);
    Kokkos::deep_copy(vgd_csnb, nb_h);
    Kokkos::deep_copy(vgd_csblk, bl_h);
  }
  // unit vector of the cell at every band position (0: none)
  Kokkos::realloc(vgd_csr, nmb, c3w, c2w, 3);
  vgd_csr_h.assign(static_cast<size_t>(nmb)*c3w*c2w*3, 0.0);
  {
    auto cr_h = Kokkos::create_mirror_view(vgd_csr);
    for (int m = 0; m < nmb; ++m) {
      const int g = g0 + m;
      for (int kk = 0; kk < c3w; ++kk) {
        for (int jj = 0; jj < c2w; ++jj) {
          int q, a, t;
          double c[3] = {0.0, 0.0, 0.0};
          if (unfold(bp[g], bj[g] + jj - w, bk[g] + kk - w, q, a, t)) {cen(q, a, t, c);}
          for (int e = 0; e < 3; ++e) {
            cr_h(m,kk,jj,e) = c[e];
            vgd_csr_h[((static_cast<size_t>(m)*c3w + kk)*c2w + jj)*3 + e] = c[e];
          }
        }
      }
    }
    Kokkos::deep_copy(vgd_csr, cr_h);
  }
  // exchange lists: receiver blocks in global order, band positions in (k, j) order;
  // the same iteration on every rank, so send and receive lists match in order
  std::vector<int> loc;
  std::map<int, std::vector<int>> snd, rcv;
  for (int gr = 0; gr < ngb; ++gr) {
    const int rr = pm->rank_eachmb[gr];
    for (int kk = 0; kk < c3w; ++kk) {
      for (int jj = 0; jj < c2w; ++jj) {
        if (jj >= w && jj < w + nx2 && kk >= w && kk < w + nx3) {continue;}
        int q, a, t;
        if (!unfold(bp[gr], bj[gr] + jj - w, bk[gr] + kk - w, q, a, t)) {continue;}
        const int go = gofi(gr, q, t/nx3, a/nx2);
        const int ro = pm->rank_eachmb[go];
        if (rr != me && ro != me) {continue;}
        const int ko = t - bk[go] + w, jo = a - bj[go] + w;
        if (rr == me && ro == me) {
          loc.insert(loc.end(), {gr - g0, kk, jj, go - g0, ko, jo});
        } else if (rr == me) {
          rcv[ro].insert(rcv[ro].end(), {gr - g0, kk, jj});
        } else {
          snd[rr].insert(snd[rr].end(), {go - pm->gids_eachrank[me], ko, jo});
        }
      }
    }
  }
  vgd_csnloc = static_cast<int>(loc.size()/6);
  // owner of every band position for direct on-rank reads; remote / none: itself
  Kokkos::realloc(vgd_csown, nmb, c3w, c2w, 3);
  {
    auto o_h = Kokkos::create_mirror_view(vgd_csown);
    for (int m = 0; m < nmb; ++m) {
      for (int kk = 0; kk < c3w; ++kk) {
        for (int jj = 0; jj < c2w; ++jj) {
          o_h(m,kk,jj,0) = m;
          o_h(m,kk,jj,1) = kk;
          o_h(m,kk,jj,2) = jj;
        }
      }
    }
    for (int e = 0; e < vgd_csnloc; ++e) {
      const int *x = &loc[6*e];
      o_h(x[0],x[1],x[2],0) = x[3];
      o_h(x[0],x[1],x[2],1) = x[4];
      o_h(x[0],x[1],x[2],2) = x[5];
    }
    Kokkos::deep_copy(vgd_csown, o_h);
  }
  Kokkos::realloc(vgd_csloc, std::max(vgd_csnloc, 1), 6);
  {
    auto l_h = Kokkos::create_mirror_view(vgd_csloc);
    for (int e = 0; e < vgd_csnloc; ++e) {
      for (int c = 0; c < 6; ++c) {l_h(e,c) = loc[6*e + c];}
    }
    Kokkos::deep_copy(vgd_csloc, l_h);
  }
  vgd_csprk.clear();
  for (const auto &x : snd) {vgd_csprk.push_back(x.first);}
  for (const auto &x : rcv) {
    if (snd.find(x.first) == snd.end()) {vgd_csprk.push_back(x.first);}
  }
  std::sort(vgd_csprk.begin(), vgd_csprk.end());
  const int np = static_cast<int>(vgd_csprk.size());
  vgd_csps.assign(np, 0);
  vgd_cspn.assign(np, 0);
  vgd_cspr.assign(np, 0);
  vgd_cspm.assign(np, 0);
  std::vector<int> sl, rl;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    vgd_csps[p] = static_cast<int>(sl.size()/3);
    vgd_cspr[p] = static_cast<int>(rl.size()/3);
    if (snd.count(rk)) {sl.insert(sl.end(), snd[rk].begin(), snd[rk].end());}
    if (rcv.count(rk)) {rl.insert(rl.end(), rcv[rk].begin(), rcv[rk].end());}
    vgd_cspn[p] = static_cast<int>(sl.size()/3) - vgd_csps[p];
    vgd_cspm[p] = static_cast<int>(rl.size()/3) - vgd_cspr[p];
  }
  vgd_csnsnd = static_cast<int>(sl.size()/3);
  vgd_csnrcv = static_cast<int>(rl.size()/3);
  Kokkos::realloc(vgd_cssl, std::max(vgd_csnsnd, 1), 3);
  Kokkos::realloc(vgd_csrl, std::max(vgd_csnrcv, 1), 3);
  {
    auto s_h = Kokkos::create_mirror_view(vgd_cssl);
    auto r_h = Kokkos::create_mirror_view(vgd_csrl);
    for (int e = 0; e < vgd_csnsnd; ++e) {
      for (int c = 0; c < 3; ++c) {s_h(e,c) = sl[3*e + c];}
    }
    for (int e = 0; e < vgd_csnrcv; ++e) {
      for (int c = 0; c < 3; ++c) {r_h(e,c) = rl[3*e + c];}
    }
    Kokkos::deep_copy(vgd_cssl, s_h);
    Kokkos::deep_copy(vgd_csrl, r_h);
  }
  // message buffers: one shell of every direction (source exchange chunked)
  const int nent = std::max(std::max(vgd_csnsnd, vgd_csnrcv), 1);
  const int nvb = std::max(vgd_n, 2*8);
  Kokkos::realloc(vgd_cssb, static_cast<size_t>(nent)*nvb);
  Kokkos::realloc(vgd_csrb, static_cast<size_t>(nent)*nvb);
  vgd_hmpi = (np > 0);
  // every side of every block carries the band: no direct neighbour reads
  Kokkos::realloc(vgd_hloc, 8*nmb);
  Kokkos::deep_copy(vgd_hloc, -1);
  vgd_hrank.assign(8*nmb, me);
  vgd_hlid.assign(8*nmb, -1);
  // compact lists (vet_gd_halo_compact > 0, MPI band): receive entry of every band
  // position (-1: own interior, on-rank owner or no cell)
  vgd_cscomp = (vgd_hcomp > 0) && vgd_hmpi;
  Kokkos::realloc(vgd_csrix, nmb, c3w, c2w);
  {
    auto x_h = Kokkos::create_mirror_view(vgd_csrix);
    for (size_t q = 0; q < x_h.size(); ++q) {x_h.data()[q] = -1;}
    for (int e = 0; e < vgd_csnrcv; ++e) {x_h(rl[3*e], rl[3*e+1], rl[3*e+2]) = e;}
    Kokkos::deep_copy(vgd_csrix, x_h);
  }
  Kokkos::realloc(vgd_csmiss, 1);
  // one cache entry per set of vet_rotate_avg (the K sets cycle); entries beyond
  // vgd_csncl stay unallocated
  vgd_csncl = 2;
  if (vgd_cscomp) {
    vgd_csncl = std::min(std::max(2, vra_k), static_cast<int>(VET_VGD_CSCL_MAX));
  }
  vgd_csper = static_cast<int64_t>(vgd_n)*indcs.nx1;
  vgd_csper = 16*((vgd_csper + 15)/16);   // padded for the packed-word lists
  for (int c = 0; c < vgd_csncl; ++c) {
    vgd_cscl[c].alpha = -9.0;
    vgd_cscl[c].use = 0;
    vgd_cscl[c].msk = decltype(vgd_cscl[c].msk)();   // masks live in vgd_csmb
  }
  if (me == 0) {
    std::cout << "<rad_vet> vet_gd (cubed sphere): rank 0 band copies on-rank "
              << vgd_csnloc << ", sent " << vgd_csnsnd << ", received " << vgd_csnrcv
              << " (cells per shell), partners " << np << "; MPI band "
              << (vgd_cscomp ? "COMPACT (vet_gd_halo_compact > 0: per shell and pass "
                               "only the read (cell, direction) values)"
                             : (vgd_hmpi ? "dense (every remote band cell, every "
                                           "direction, per shell and pass)"
                                         : "none"))
              << "; list cache " << vgd_csncl << " sets" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsSource
//! \brief band of ln chi, ln S (vgd_cs) for every shell of the sweep

void RadiationVet::VetGdCsSource() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ilo = indcs.is, ie = indcs.ie;
  VetGdCsHalo(vgd_cs, 2, ilo, ie, true);
  Kokkos::fence();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsTab
//! \brief host tables of direction set d on cs: no walls, identity map, branch n . r_hat
//! of every band cell

void RadiationVet::VetGdCsTab(const std::vector<double> &d) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c2 = indcs.nx2 + 2*vgd_w, c3 = indcs.nx3 + 2*vgd_w;
  const int n = vgd_n;
  VgdTabHost &H = vgd_tabh;
  for (int m = 0; m < nmb; ++m) {
    for (int k = 0; k < c3; ++k) {
      for (int j = 0; j < c2; ++j) {
        H.w(m,k,j) = 0;
        const double *rh = &vgd_csr_h[((static_cast<size_t>(m)*c3 + k)*c2 + j)*3];
        for (int q = 0; q < n; ++q) {
          H.mp(m,k,j,q) = q;
          H.mr(m,k,j,q) = d[4*q]*rh[0] + d[4*q+1]*rh[1] + d[4*q+2]*rh[2];
        }
      }
    }
  }
  H.nwall = 0;
  H.nexact = 0;
  H.nmap = 0;
}

}  // namespace radvet
