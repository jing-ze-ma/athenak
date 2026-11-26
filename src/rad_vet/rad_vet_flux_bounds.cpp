//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_flux_bounds.cpp
//! \brief bounds of the face flux: |F| <= c E (clean_stag_fcone), the projection of the
//! steady balance's residual (clean_stag_fwb).
//!
//!  called from: rad_vet_flux.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_flux_seam.cpp, rad_vet_init.cpp

#include <algorithm>
#include <array>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "rad_vet/rad_vet_flux.hpp"
#include "utils/p_ledger.hpp"
#include "utils/repro_sum.hpp"
#include "utils/scratch_view.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluxCone
//! \brief clean_stag_fcone: the cell flux vector (CellVector, F_i = sum_f F_f
//! u_f) within the realisable cone |F_i| <= c E_i (E of work vector e), by construction.
//! Each face is scaled towards 0 by l_f = min(l_L, l_R), l_i = min(1, c E_i/S_i), S_i^2 =
//! sum_fg max(0, F_f u_f . F_g u_g), so |F_i| <= l_i S_i <= c E_i.  One value per face
//! (cubed sphere: the lower gid's) keeps flux form.  The gas takes the absorption of the
//! scaled F (CLF_N + 1), G + d is kept, the momentum cut is booked and printed (CLEAN
//! FCONE).

void RadiationVet::FluxCone(const int e, const int kd) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  SolveArgs x = MakeArgs();
  const bool cub = x.cub, dd = cln_fkin && e == CLW_X && cln_dok, wbc = e == CLW_X;
  // (clean_stag_fcone_lz: l_i of the faces without L_z arm (fa), CLW_CP >= 0 with the
  // others' 1, or -l_i of the others with fa's 0, where those alone exceed c E;
  // clean_lz_spin: the others' cut withheld, its L_z X the gas's rigid rotation; the
  // cubed sphere's gas phi momentum is lz_tel's torque, which holds the cut already)
  const bool lz = cln_fclz, spn = cln_lzsp, sph = x.sph;
  auto fa = KOKKOS_LAMBDA(const int q) {return !LzFamily(sph, cub, q);};
  auto lf = KOKKOS_LAMBDA(const Real cp, const int q) {
    return !lz ? cp : (fa(q) ? fmax(cp, 0.0) : ((cp < 0.0) ? -cp : 1.0));
  };
  const bool emx = cln_fcemx;   // (clean_stag_fcone_emax)
  if (cln_fconc.extent_int(0) != 5) Kokkos::realloc(cln_fconc, 5);
  if (spn && cln_fcx.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(cln_fcx, nmb1 + 1, 3, cln_w.s.extent_int(2) + 1,
                    cln_w.s.extent_int(3) + 1, cln_w.s.extent_int(4) + 1);
  }
  auto fcx = cln_fcx;
  Kokkos::deep_copy(cln_fconc, 0.0);
  auto fcc = cln_fconc;
  SolveFor("cln_fcone_l", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // l_i = min(1, c E_i/S_i), S_i^2 = sum_fg max(0, v_f.v_g) of the faces' vectors v_f =
    // F_f u_f: |sum_f l_f v_f|^2 <= l_i^2 S_i^2 for every l_f in [0, l_i], S_i = |F_i| of
    // a beam on an orthogonal mesh.  (sa, sb, sc: the faces with L_z arm, the two kinds'
    // products, the others, clean_stag_fcone_lz)
    Real v[6][3], sa = 0.0, sb = 0.0, sc = 0.0;
    for (int d = 0; d < 2*x.nf; ++d) {
      const int q = d/2, h = d%2;
      Real c[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      c[d] = x.fc[q](m,CLF_N,k+h*(q == 2),j+h*(q == 1),i+h*(q == 0));
      CellVector(x, m, k, j, i, c, v[d]);
      for (int g = 0; g <= d; ++g) {
        const Real p = ((g < d) ? 2.0 : 1.0)*
            fmax(v[d][0]*v[g][0] + v[d][1]*v[g][1] + v[d][2]*v[g][2], 0.0);
        const int na = (lz && !fa(q)) + (lz && !fa(g/2));
        if (na == 2) {
          sa += p;
        } else if (na == 1) {
          sb += 0.5*p;
        } else {
          sc += p;
        }
      }
    }
    const Real ce = x.c*fmax(x.w(m,e,k,j,i), 0.0), r = SQR(ce) - sa;
    x.w(m,CLW_CP,k,j,i) = (r < 0.0) ? -ce/sqrt(sa) :
        ((2.0*sb + sc > r) ? (sqrt(SQR(sb) + sc*r) - sb)/sc : 1.0);
  });
  GhostFill(cln_w, CLW_CP, 1);
  for (int d = 0; d < x.nf; ++d) {
    const int dk = (d == 2), dj = (d == 1), di = (d == 0);
    auto fc = x.fc[d];
    SolveFor("cln_fcone_f", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // a physical boundary face: its inner cell's l, after the cubed sphere's tie; an
      // opaque face (tau_f = chi_f dx_f) keeps its share tau/(1 + tau) of F (a diffusion
      // front's face mean is no free-streaming excess)
      const bool bd = IsBoundaryFace(x, d, m, k, j, i);
      const Real tau = 0.5*(x.op(m,VET_OP_T,k,j,i) + x.op(m,VET_OP_T,k-dk,j-dj,i-di))*
                       CellFaceDistance(x, d, m, k, j, i);
      const Real l = fmin(lf(x.w(m,CLW_CP,k-dk,j-dj,i-di), d),
                          lf(x.w(m,CLW_CP,k,j,i), d));
      Real fs = (bd ? 1.0 : 1.0 - (1.0 - l)/(1.0 + tau))*fc(m,CLF_N,k,j,i);
      if (emx && !bd) {   // (clean_stag_fcone_emax: |F_f| <= c max(E_L, E_R))
        const Real ce = x.c*fmax(fmax(x.w(m,e,k,j,i), x.w(m,e,k-dk,j-dj,i-di)), 0.0);
        fs = fmin(fmax(fs, -ce), ce);
      }
      fc(m,CLF_N+2,k,j,i) = fs;
    });
  }
  if (cub) {   // a face two MeshBlocks share: its lower gid's value
    auto v = cln_vf;
    SolveFor("cln_fcone_v", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int q = 0; q < 3; ++q) {
        v(m,9+2*q,k,j,i) = x.fc[q](m,CLF_N+2,k,j,i);
        v(m,9+2*q+1,k,j,i) = x.fc[q](m,CLF_N+2,k+(q == 2),j+(q == 1),i+(q == 0));
      }
    });
    const int f0 = (pmy_pack->pmesh->nmb_rootx1 == 1 && !pmy_pack->pmesh->multilevel)
                   ? 2 : 0;
    GhostFill(cln_vf, 9 + f0, 6 - f0, true, 0, 1, FluxMask(f0));
    FluxCsShared();
  }
  auto sl = ScratchView<DvceArray1D<Real>>("cln_csl", 6*CLN_NSL);
  if (wbc) Kokkos::deep_copy(sl, 0.0);
  for (int d = 0; d < x.nf; ++d) {
    const int dk = (d == 2), dj = (d == 1), di = (d == 0);
    auto fc = x.fc[d];
    SolveFor("cln_fcone_b", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // a physical boundary face has no F equation: the gas takes its half cell's
      // absorption of fs (FaceMomentumChange)
      const bool bd = IsBoundaryFace(x, d, m, k, j, i);
      const bool lin = (d == 0) ? (i > is) : ((d == 1) ? (j > js) : (k > ks));
      const Real f = fc(m,CLF_N,k,j,i);
      const Real fs = !bd ? fc(m,CLF_N+2,k,j,i) :
                      lf(x.w(m,CLW_CP,k-lin*dk,j-lin*dj,i-lin*di), d)*f;
      fc(m,CLF_N+2,k,j,i) = fs;
      if (spn) fcx(m,d,k,j,i) = 0.0;
      if (fs == f) return;
      const Real it = 1.0/fc(m,CLF_TH,k,j,i);
      // (clean_lz_spin: the arm faces' cut withheld from the gas per face, WriteBack)
      if (spn && !bd && !fa(d) && !cub) fcx(m,d,k,j,i) = (f - fs)*it/(x.c*x.c);
      const bool wh = !bd && !(lz && !fa(d) && !(spn && !cub));   // (withheld)
      if (wh) fc(m,CLF_N+1,k,j,i) -= (f - fs)*it/(x.h*x.c*x.c);
      if (dd && !bd) {   // (the energy flux's memory G + d kept)
        fc(m,CLF_N+3,k,j,i) += f - fs;
        if (kd > 0) fc(m,kd,k,j,i) += (f - fs)/x.h;
      }
      fc(m,CLF_N,k,j,i) = fs;
      if (!wbc) return;   // counters: the write-back's, as printed
      SlotAdd(sl, 0, m, k, j, i, 1.0);
      const Real dp = (f - fs)*it*0.5*(CellVolume(x, m, k, j, i) +
                      CellVolume(x, m, k-dk, j-dj, i-di))/(x.c*x.c);
      SlotAdd(sl, 1, m, k, j, i, fabs(dp));
      // (the momentum withheld from the gas along each face family, signed, its own
      // cells' halves: fcone_lz, p_ledger)
      const bool rin = (d == 0) ? (i <= ie) : ((d == 1) ? (j <= je) : (k <= ke));
      if (wh) SlotAdd(sl, 2 + d, m, k, j, i, (f - fs)*it*0.5*((rin ? CellVolume(x, m,
                      k, j, i) : 0.0) + (lin ? CellVolume(x, m, k-dk, j-dj, i-di) :
                      0.0))/(x.c*x.c));
    });
  }
  if (wbc) SlotSum(sl, fcc, 5);
  if (e == CLW_X) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_fconc);
    const Real n = GlobalReduce(h(0), false), p = GlobalReduce(h(1), false);
    if (global_variable::my_rank == 0 && cln_verbose > 0 && n > 0.0) {
      std::printf("CLEAN FCONE cycle %d t2s %d faces %.0f momentum cut sum |dF| V/c^2 "
                  "%.6e\n", pmy_pack->pmesh->ncycle, t2_solve, n, p);
    }
    Real q[3];   // (every rank: the sums are collective)
    for (int d = 0; d < 3; ++d) q[d] = GlobalReduce(h(2+d), false);
    if (cln_feedback && !cln_vzero) {
      p_ledger::Withheld(0, q, PSlope(t2_solve), dt_sub, !t2_inprobe);
    }
    if (n > 0.0 && cln_fclz) {
      if (global_variable::my_rank == 0) {
        std::printf("CLEAN FCONE_LZ cycle %d t2s %d momentum removed (signed, V/c^2) x1 "
                    "%.6e x2 %.6e x3 %.6e\n", pmy_pack->pmesh->ncycle, t2_solve, q[0],
                    q[1], q[2]);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RingMatesInit
//! \brief the MeshBlocks of each local block's line along d, position by position
//! (FluxBalanceProject): src(m, p) >= 0 the local block, else -1 - slot of a received
//! one; snd / rcv the blocks exchanged (block or slot, rank, gid), one message each

void RadiationVet::RingMatesInit(const int d, RingMates &rm) {
  Mesh *pm = pmy_pack->pmesh;
  const int nr[3] = {pm->nmb_rootx1, pm->nmb_rootx2, pm->nmb_rootx3};
  const int nmb = pmy_pack->nmb_thispack, me = global_variable::my_rank;
  rm.nmb = nmb;
  rm.nt = pm->nmb_total;
  rm.snd.clear();
  rm.rcv.clear();
  std::map<int, int> slot;   // (received gid -> slot)
  std::set<std::array<int, 3>> sent;
  Kokkos::realloc(rm.src, nmb, nr[d]);
  auto h = Kokkos::create_mirror_view(rm.src);
  auto loc = [&](const int g) {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    return std::array<int, 4>{static_cast<int>(l.lx1), static_cast<int>(l.lx2),
                              static_cast<int>(l.lx3), l.level};
  };
  for (int m = 0; m < nmb; ++m) {
    const int gm = pmy_pack->pmb->mb_gid.h_view(m);
    std::array<int, 4> w = loc(gm);
    const int p0 = w[d];
    for (int p = 0; p < nr[d]; ++p) {
      w[d] = p;
      int g = 0;
      while (g < pm->nmb_total && loc(g) != w) ++g;
      if (g == pm->nmb_total) RadVetFatal("fwb: no MeshBlock at a line's position");
      const int rk = pm->rank_eachmb[g];
      if (rk == me) {
        h(m,p) = g - pm->gids_eachrank[me];
        continue;
      }
      if (slot.find(g) == slot.end()) {
        const int k = static_cast<int>(slot.size());
        slot[g] = k;
        rm.rcv.push_back({k, rk, g});
      }
      h(m,p) = -1 - slot[g];
      if (p != p0) sent.insert({m, rk, gm});
    }
  }
  rm.snd.assign(sent.begin(), sent.end());
  rm.nslot = static_cast<int>(slot.size());
  Kokkos::deep_copy(rm.src, h);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluxBalanceProject
//! \brief clean_stag_fwb_ring (d = 2): fwb's residual s of
//! the faces along d (rs, its weight w/(1 + tau) rw, FluxSolve) less mu rw,
//! mu = sum V_f s/sum V_f rw over the whole line through the mesh (a radial line, a phi
//! ring), added to face slot slot.  A radial line's faces share their direction, a
//! ring's their lever arm R sin(theta): the net s V of each line is 0, so fwb gives no
//! momentum along it and no L_z.  A line across MeshBlocks: each block's sums at its
//! position along d in one mesh-wide array, summed over the ranks (every slot from one
//! block, the others 0: exact), then over the positions in order, so mu does not depend
//! on the rank count (one MeshBlock along d: the block's own sums).  Tile mode: the
//! sums of each tile along a lateral d, then the tiles in order, so mu does not depend on
//! the MeshBlocks either

void RadiationVet::FluxBalanceProject(const int d, const int slot,
                                      const DvceArray4D<Real> &rs,
                                      const DvceArray4D<Real> &rw, const bool cnt) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks, nmb = pmy_pack->nmb_thispack;
  Mesh *pm = pmy_pack->pmesh;
  const int nx[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int nr[3] = {pm->nmb_rootx1, pm->nmb_rootx2, pm->nmb_rootx3};
  const int da = (d == 0) ? 1 : 0, db = (d == 2) ? 1 : 2;   // transverse dirs
  const int n = nx[d], na = nx[da], nb = nx[db];
  const int dk = (d == 2), dj = (d == 1), di = (d == 0);
  const bool wrap = nr[d] == 1 && pm->mesh_bcs[2*d] == BoundaryFlag::periodic;
  SolveArgs x = MakeArgs();
  auto fc = x.fc[d];
  auto fwc = cln_fwbc;
  const int tn = (d > 0 && repro::tile_n > 0) ? repro::tile_n : n, nc = n/tn;
  auto mu = ScratchView<DvceArray4D<Real>>("cln_fwmu" + std::to_string(d), nmb, nb, na,
                                           2*nc);
  auto face = KOKKOS_LAMBDA(const int t, const int a, const int b, int &k, int &j,
                            int &i) {
    i = is + ((d == 0) ? t : ((da == 0) ? a : b));
    j = js + ((d == 1) ? t : ((da == 1) ? a : b));
    k = ks + ((d == 2) ? t : b);
  };
  par_for("cln_fwbp", DevExeSpace(), 0, nmb - 1, 0, nb - 1, 0, na - 1,
  KOKKOS_LAMBDA(const int m, const int b, const int a) {
    for (int c = 0; c < nc; ++c) {
      Real ss = 0.0, sw = 0.0;
      for (int t = c*tn; t < (c + 1)*tn; ++t) {   // (V_f: the face's weight, own faces)
        int k, j, i;
        face(t, a, b, k, j, i);
        const Real vf = CellVolume(x, m, k, j, i) + CellVolume(x, m, k-dk, j-dj, i-di);
        ss += vf*rs(m,k,j,i);
        sw += vf*rw(m,k,j,i);
      }
      mu(m,b,a,2*c) = ss;
      mu(m,b,a,2*c+1) = sw;
    }
  });
  int nq = nc;
  auto mq = mu;
  if (nr[d] > 1) {   // (the other MeshBlocks' sums, then all positions in order)
    RingMates &rm = cln_fwr[d];
    if (rm.nmb != nmb || rm.nt != pm->nmb_total) RingMatesInit(d, rm);
    auto rb = ScratchView<DvceArray4D<Real>>("cln_fwrb" + std::to_string(d),
                                             std::max(rm.nslot, 1), nb, na, 2*nc);
    mq = ScratchView<DvceArray4D<Real>>("cln_fwmq" + std::to_string(d), nmb, nb, na, 2);
#if MPI_PARALLEL_ENABLED
    const int nw = nb*na*2*nc;
    std::vector<MPI_Request> rq(rm.rcv.size() + rm.snd.size());
    DevExeSpace().fence();
    size_t r = 0;
    for (auto &v : rm.rcv) {
      MPI_Irecv(rb.data() + static_cast<size_t>(v[0])*nw, nw, MPI_ATHENA_REAL, v[1],
                4250 + v[2], MPI_COMM_WORLD, &rq[r++]);
    }
    for (auto &v : rm.snd) {
      MPI_Isend(mu.data() + static_cast<size_t>(v[0])*nw, nw, MPI_ATHENA_REAL, v[1],
                4250 + v[2], MPI_COMM_WORLD, &rq[r++]);
    }
    MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
#endif
    const int np = nr[d], ncc = nc;
    auto src = rm.src;
    auto mo = mq;
    nq = 1;
    par_for("cln_fwbr", DevExeSpace(), 0, nmb - 1, 0, nb - 1, 0, na - 1, 0, 1,
    KOKKOS_LAMBDA(const int m, const int b, const int a, const int c) {
      Real s = 0.0;
      for (int p = 0; p < np; ++p) {
        const int sl = src(m,p);
        for (int t = 0; t < ncc; ++t) {
          s += (sl >= 0) ? mu(sl,b,a,2*t+c) : rb(-1-sl,b,a,2*t+c);
        }
      }
      mo(m,b,a,c) = s;
    });
  }
  par_for("cln_fwbq", DevExeSpace(), 0, nmb - 1, 0, nb - 1, 0, na - 1,
  KOKKOS_LAMBDA(const int m, const int b, const int a) {
    Real ss = 0.0, sw = 0.0;
    for (int c = 0; c < nq; ++c) {
      ss += mq(m,b,a,2*c);
      sw += mq(m,b,a,2*c+1);
    }
    const Real q = (sw > 0.0) ? ss/sw : 0.0;
    Real cn = 0.0, cp = 0.0, c1 = 0.0;   // (the line's counters: one atomic each)
    for (int t = 0; t <= n; ++t) {   // (one MeshBlock round a periodic wrap: low face)
      int k, j, i, kr, jr, ir;
      face(t, a, b, k, j, i);
      face((wrap && t == n) ? 0 : t, a, b, kr, jr, ir);
      const Real s = rs(m,kr,jr,ir) - q*rw(m,kr,jr,ir);
      fc(m,slot,k,j,i) += s/(x.h*x.c*x.c);
      if (!cnt || t == n) continue;   // counters: the printed s, the write-back's
      const Real p = fc(m,CLF_TH,k,j,i)*s*0.5*(CellVolume(x, m, k, j, i) +
                     CellVolume(x, m, k-dk, j-dj, i-di))/(x.c*x.c);
      cn += 1.0;
      cp += fabs(p);
      c1 += p;
    }
    if (!cnt) return;
    Kokkos::atomic_add(&fwc(0), cn);
    Kokkos::atomic_add(&fwc(1), cp);
    if (d == 0) Kokkos::atomic_add(&fwc(2), c1);
  });
}

}  // namespace radvet
