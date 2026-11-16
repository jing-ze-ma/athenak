//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean_solver.cpp
//! \brief Linear algebra of the clean implicit radiation solve: face coefficients, the
//! operator A E = E (1 + h c kE (1-r)) + h div F_lin(E), an x1-line tridiagonal
//! preconditioner and right-preconditioned BiCGStab.

#include <cmath>
#include <limits>

#include "athena.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_clean.hpp"
#include "rad_m1/rad_m1_clean_cell.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn CleanArgs RadiationM1::CleanMakeArgs
//! \brief collects the arrays and constants the device functions need

CleanArgs RadiationM1::CleanMakeArgs() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  CleanArgs x;
  x.cs = cln_s;
  x.w = cln_w;
  x.fc[0] = cln_fc1;
  x.fc[1] = cln_fc2;
  x.fc[2] = cln_fc3;
  x.beam = ibm_f;
  x.bstress = ibm_s;
  x.u = cln_x;
  x.h = dt_sub;
  x.c = c_light;
  x.idx[0] = 1.0/size.h_view(0).dx1;
  x.idx[1] = 1.0/size.h_view(0).dx2;
  x.idx[2] = 1.0/size.h_view(0).dx3;
  const int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  // an outflow face is the open one without an entering bath
  for (int s = 0; s < 6; ++s) {
    x.bc[s] = (bc[s] == M1_IBC_X2FLUX) ? M1_IBC_FLUX :
              ((bc[s] == M1_IBC_OUTFLOW) ? M1_IBC_MARSHAK : bc[s]);
  }
  x.iflux[0] = iflux_x1min;
  x.iflux[1] = iflux_x1max;
  x.ifx2[0] = ifx2_lo;
  x.ifx2[1] = ifx2_hi;
  x.mbc = pmy_pack->pmb->mb_bcs.d_view;
  for (int n = 0; n < 6; ++n) x.ebath[n] = 0.0;
  x.ebath[0] = iebath_x1min;
  x.ebath[1] = iebath_x1max;
  for (int s = 0; s < 2; ++s) x.ebath[s] *= (bc[s] != M1_IBC_OUTFLOW);
  x.marshak_q = marshak_q;
  x.marshak = (cln_bdry_form == 1);
  x.use_beam = (ibm_n > 0);
  x.feedback = cln_feedback;
  x.vzero = cln_vzero;
  x.recon = cln_recon;
  x.all_low = cln_all_low;
  x.consist = cln_consist;
  x.milne = cln_milne;
  x.fvrays = (cln_rfv > 0);
  x.hi = false;
  x.fs = false;
  x.fd = cln_fsolve && cln_fkin;
  x.gf = cln_gflux;
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
  x.cub = pmy_pack->pmesh->use_cubed_sphere;
  if (x.sph || (!x.cub && cln_rfv > 0)) {
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
  }
  if (x.cub) {
    x.vol = pmy_pack->pcoord->volume;
    for (int d = 0; d < 3; ++d) {
      x.area[d] = cln_ga[d];
      x.dxf[d] = cln_gd[d];
    }
    x.g = cln_gn;
    x.x1v = pmy_pack->pcoord->x1v;
    x.x1f = pmy_pack->pcoord->xx1f;
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  x.etg = fl.on && fl.use_etotgrav;
  x.mhd = fl.on && fl_mhd;
  x.phi = fl.phicc0;
  x.emag = emag0;
  return x;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanGhostFill
//! \brief ghost cells of components [c0, c0+nc) of a work array: the neighbours' cells
//! (other MeshBlocks, periodic wrap; across a cubed-sphere seam the neighbour's own
//! cells, without the along-seam resample) by the cell-centred boundary exchange on a
//! scratch array of nc components (only the ghosts and the ng layers the exchange reads
//! are copied in, only the ghosts back), then (phys) a copy of the edge cell on physical
//! boundaries.  part = 1 only starts the exchange (nc <= 32), part = 2 finishes it: work
//! on the cells CleanInner by ng can run in between.  The closure's off-diagonal D12,
//! D13, D23 take the edge cell's with their parity: odd across a reflecting radiation
//! wall normal to one of their axes, and D_r,theta, D_r,phi across the pole (e_theta,
//! e_phi flip), so that d_b(D E) puts no shear traction on the wall (as ce91f387d)

void RadiationM1::CleanGhostFill(DvceArray5D<Real> a, const int c0, const int nc,
                                 const bool phys, const int part) {
  if (nc > 32) {   // the exchange moves at most 32 variables at once
    CleanGhostFill(a, c0, 32, phys);
    CleanGhostFill(a, c0 + 32, nc - 32, phys);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ng = indcs.ng;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = a.extent_int(4), n2 = a.extent_int(3), n3 = a.extent_int(2);
  CleanHalo &hl = cln_halo[nc];
  if (!hl.pb) {
    // (the exchanges run one at a time: their scratch arrays share one allocation)
    if (cln_hbuf.extent_int(0) == 0) Kokkos::realloc(cln_hbuf, nmb1 + 1, 32, n3, n2, n1);
    hl.a = DvceArray5D<Real>(cln_hbuf.data(), nmb1 + 1, nc, n3, n2, n1);
    Kokkos::realloc(hl.ca, nmb1 + 1, nc, 1, 1, 1);
    hl.pb.reset(new MeshBoundaryValuesCC(pmy_pack, nullptr, false));
    hl.pb->InitializeBuffers(nc);
    hl.pb->SetVectorPairs(nc, {});
    hl.pb->cs_noresample = true;
  }
  auto h = hl.a;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  MeshBoundaryValuesCC *pb = hl.pb.get();
  if (part != 2) {
    par_for("cln_hcopy", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      if (!CleanInner(k, j, i, ng, is, ie, js, je, ks, ke, md, td)) {
        h(m,n,k,j,i) = a(m,c0+n,k,j,i);
      }
    });
    while (pb->InitRecv(nc) == TaskStatus::incomplete) {}
    while (pb->PackAndSendCC(hl.a, hl.ca) == TaskStatus::incomplete) {}
    if (part == 1) return;
  }
  while (pb->RecvAndUnpackCC(hl.a, hl.ca) == TaskStatus::incomplete) {}
  while (pb->ClearSend() == TaskStatus::incomplete) {}
  while (pb->ClearRecv() == TaskStatus::incomplete) {}
  par_for("cln_hback", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    if (!CleanInner(k, j, i, 0, is, ie, js, je, ks, ke, md, td)) {
      a(m,c0+n,k,j,i) = h(m,n,k,j,i);
    }
  });
  if (!phys) return;

  auto bc = pmy_pack->pmb->mb_bcs.d_view;
  auto pf = [] KOKKOS_FUNCTION (const BoundaryFlag f) {
    return f != BoundaryFlag::block && f != BoundaryFlag::periodic &&
           f != BoundaryFlag::panel;
  };
  // (the cubed sphere's D is Cartesian: a radial wall mirrors it about n, not done here)
  const bool ts = (a.data() == cln_s.data()) && !pmy_pack->pmesh->use_cubed_sphere;
  const int ib[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  bool rf[6];
  for (int s = 0; s < 6; ++s) rf[s] = (ib[s] == M1_IBC_REFLECT);
  auto sg = [=] KOKKOS_FUNCTION (const int c, const int s, const int m) {
    const int p = c - CLN_D11 - 3;   // (D12, D13, D23: axes 01, 02, 12; p != 2 - d)
    if (!ts || p < 0 || p > 2) return 1.0;
    if (bc(m,s) == BoundaryFlag::polar) return (p < 2) ? -1.0 : 1.0;
    return (rf[s] && p != 2 - s/2) ? -1.0 : 1.0;
  };
  par_for("cln_ghost1", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, n2-1, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int g) {
    const Real sl = sg(c0+n, 0, m), sh = sg(c0+n, 1, m);
    if (pf(bc(m,BoundaryFace::inner_x1))) a(m,c0+n,k,j,is-1-g) = sl*a(m,c0+n,k,j,is);
    if (pf(bc(m,BoundaryFace::outer_x1))) a(m,c0+n,k,j,ie+1+g) = sh*a(m,c0+n,k,j,ie);
  });
  if (pmy_pack->pmesh->multi_d) {
    par_for("cln_ghost2", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, ng-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int g, const int i) {
      const Real sl = sg(c0+n, 2, m), sh = sg(c0+n, 3, m);
      if (pf(bc(m,BoundaryFace::inner_x2))) a(m,c0+n,k,js-1-g,i) = sl*a(m,c0+n,k,js,i);
      if (pf(bc(m,BoundaryFace::outer_x2))) a(m,c0+n,k,je+1+g,i) = sh*a(m,c0+n,k,je,i);
    });
  }
  if (pmy_pack->pmesh->three_d) {
    par_for("cln_ghost3", DevExeSpace(), 0, nmb1, 0, nc-1, 0, ng-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int g, const int j, const int i) {
      const Real sl = sg(c0+n, 4, m), sh = sg(c0+n, 5, m);
      if (pf(bc(m,BoundaryFace::inner_x3))) a(m,c0+n,ks-1-g,j,i) = sl*a(m,c0+n,ks,j,i);
      if (pf(bc(m,BoundaryFace::outer_x3))) a(m,c0+n,ke+1+g,j,i) = sh*a(m,c0+n,ke,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanReduce
//! \brief maximum (max) or sum of v over all ranks

Real RadiationM1::CleanReduce(Real v, const bool max) {
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_ATHENA_REAL, max ? MPI_MAX : MPI_SUM,
                MPI_COMM_WORLD);
#endif
  return v;
}

//----------------------------------------------------------------------------------------
//! \fn std::array<Real, 2> RadiationM1::CleanDots
//! \brief the n <= 2 sums over active cells and the nc components of w(ab[2q])
//! w(ab[2q+1]) (q < n): one pass, one reduction over the ranks

std::array<Real, 2> RadiationM1::CleanDots(DvceArray5D<Real> w, const int nc, const int n,
                                           const std::array<int, 4> ab) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  std::array<Real, 2> d = {0.0, 0.0};
  const int a0 = ab[0], a1 = ab[1], a2 = ab[2], a3 = ab[3];   // (std::array is host-only)
  Kokkos::parallel_reduce("cln_dot", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &s0, Real &s1) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    for (int c = 0; c < nc; ++c) {
      s0 += w(m,a0+c,k+ks,j+js,i+is)*w(m,a1+c,k+ks,j+js,i+is);
      if (n > 1) s1 += w(m,a2+c,k+ks,j+js,i+is)*w(m,a3+c,k+ks,j+js,i+is);
    }
  }, Kokkos::Sum<Real>(d[0]), Kokkos::Sum<Real>(d[1]));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, d.data(), n, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFaceCoefficients
//! \brief per face: the blend weight w(C1), theta, the old-vector flux G, and the
//! half-range factors h+-_f.  h+_f = H+_f/J_L and h-_f = H-_f/J_R, with H+-_f the rays'
//! half-range fluxes reconstructed to the face from the upwind side (minmod-limited);
//! |h+-_f| <= 1 keeps the operator an M-matrix.

void RadiationM1::CleanFaceCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real h = dt_sub, c = c_light;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  const bool stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  auto cs = cln_s;
  auto op = opac;
  auto inc = t2inc;
  const bool fsk = cln_fsolve, gfx = cln_gflux;

  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = M1_T2_F1 + dir;

    par_for("cln_facecoef", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      const Real chi = 0.5*(op(m,M1_OP_T,k,j,i) + op(m,M1_OP_T,kl,jl,il));
      const Real c1 = fmax(cs(m,CLN_C1,k,j,i), cs(m,CLN_C1,kl,jl,il));
      fc(m,CLF_W,k,j,i) = CleanSmoothStep((c1 - c1lo)/(c1hi - c1lo));
      fc(m,CLF_TH,k,j,i) = 1.0/(1.0 + h*c*chi);
      fc(m,CLF_G,k,j,i) = f0(m,k,j,i) + (stage ? inc(m,tf,k,j,i) : 0.0);
      if (fsk) fc(m,CLF_N+2,k,j,i) = fc(m,CLF_G,k,j,i);   // (fsolve: start of its F)
      if (gfx) fc(m,CLF_SLO,k,j,i) = 0.0;   // (gflux: every face high order)

      // rays' half-range fluxes H+ = h+ J (from the left) and H- = h- J (from the right)
      auto hp = [&](const int kk, const int jj, const int ii) {
        return cs(m,CLN_HP+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);
      };
      auto hm = [&](const int kk, const int jj, const int ii) {
        return cs(m,CLN_HM+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);
      };
      auto minmod = [](const Real a, const Real b) {
        return (a*b <= 0.0) ? 0.0 : ((fabs(a) < fabs(b)) ? a : b);
      };
      const Real hpl = hp(kl,jl,il);
      const Real hpf = hpl + 0.5*minmod(hpl - hp(kl-dk,jl-dj,il-di), hp(k,j,i) - hpl);
      const Real hmr = hm(k,j,i);
      const Real hmf = hmr - 0.5*minmod(hmr - hm(kl,jl,il), hm(k+dk,j+dj,i+di) - hmr);
      const Real jl_ = cs(m,CLN_J,kl,jl,il), jr_ = cs(m,CLN_J,k,j,i);
      fc(m,CLF_HP,k,j,i) = (jl_ > 0.0) ? fmin(fmax(hpf/jl_, 0.0), 1.0)
                                       : cs(m,CLN_HP+dir,kl,jl,il);
      fc(m,CLF_HM,k,j,i) = (jr_ > 0.0) ? fmax(fmin(hmf/jr_, 0.0), -1.0)
                                       : cs(m,CLN_HM+dir,k,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanApply
//! \brief w(dst) = A w(src) over the active cells (the linear part of the operator)

void RadiationM1::CleanApply(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanGhostFill(cln_w, src, 1);
  CleanArgs x = CleanMakeArgs();
  auto op = opac;

  par_for("cln_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
    x.w(m,dst,k,j,i) = x.w(m,src,k,j,i)*(1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*qr)
                       + CleanFluxDivergence(x, src, m, k, j, i, false);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPreconditionerSetup
//! \brief coefficients of the preconditioner: the x1-line tridiagonal part of A plus the
//! diagonal of its x2 and x3 parts (the off-diagonal Eddington terms are left to the
//! Krylov iteration)

void RadiationM1::CleanPreconditionerSetup() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto op = opac;

  par_for("cln_pcsetup", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real al_p, ar_p, al_m, ar_m;
    CleanFaceWeights(x, 0, m, k, j, i+1, al_p, ar_p);
    CleanFaceWeights(x, 0, m, k, j, i, al_m, ar_m);
    const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
    Real a0 = 1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*qr + x.h*x.idx[0]*(al_p - ar_m);
    x.w(m,CLW_AM,k,j,i) = -x.h*x.idx[0]*al_m;
    x.w(m,CLW_AP,k,j,i) = x.h*x.idx[0]*ar_p;
    if (x.sph || x.cub) {   // (A_f/V weights, CleanDivW)
      const Real wp = CleanDivW(x, 0, m, k, j, i, 1), wm = CleanDivW(x, 0, m, k, j, i, 0);
      a0 = 1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*qr + x.h*(wp*al_p - wm*ar_m);
      x.w(m,CLW_AM,k,j,i) = -x.h*wm*al_m;
      x.w(m,CLW_AP,k,j,i) = x.h*wp*ar_p;
    }
    if (x.multi_d) {
      CleanFaceWeights(x, 1, m, k, j+1, i, al_p, ar_p);
      CleanFaceWeights(x, 1, m, k, j, i, al_m, ar_m);
      a0 += (x.sph || x.cub) ? x.h*(CleanDivW(x, 1, m, k, j, i, 1)*al_p -
                         CleanDivW(x, 1, m, k, j, i, 0)*ar_m) : x.h*x.idx[1]*(al_p - ar_m);
    }
    if (x.three_d) {
      CleanFaceWeights(x, 2, m, k+1, j, i, al_p, ar_p);
      CleanFaceWeights(x, 2, m, k, j, i, al_m, ar_m);
      a0 += x.cub ? x.h*(CleanDivW(x, 2, m, k, j, i, 1)*al_p -
                         CleanDivW(x, 2, m, k, j, i, 0)*ar_m)
                  : x.h*x.idx[2]*(al_p - ar_m);
    }
    x.w(m,CLW_A0,k,j,i) = a0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPrecondition
//! \brief w(dst) = M^{-1} w(src): tridiagonal (Thomas) solve along every x1 line

void RadiationM1::CleanPrecondition(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;

  par_for("cln_thomas", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    Real cprev = 0.0, dprev = 0.0;
    for (int i = is; i <= ie; ++i) {
      const Real am = (i > is) ? w(m,CLW_AM,k,j,i) : 0.0;
      const Real ap = (i < ie) ? w(m,CLW_AP,k,j,i) : 0.0;
      const Real den = w(m,CLW_A0,k,j,i) - am*cprev;
      cprev = ap/den;
      dprev = (w(m,src,k,j,i) - am*dprev)/den;
      w(m,CLW_CP,k,j,i) = cprev;
      w(m,dst,k,j,i) = dprev;
    }
    for (int i = ie-1; i >= is; --i) {
      w(m,dst,k,j,i) -= w(m,CLW_CP,k,j,i)*w(m,dst,k,j,i+1);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanBiCGStab
//! \brief right-preconditioned BiCGStab for A X = B, starting from the X it is given:
//! scalar E vectors of cln_w (staggered layout) or 4-component vectors of cln_x (cell
//! layout).  Returns the number of iterations, negative if not converged.

int RadiationM1::CleanBiCGStab() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool cell = cln_cell;
  auto w = cell ? cln_x : cln_w;
  const int nc = cell ? 4 : 1;
  const int sx = cell ? CLC_X : CLW_X, sb = cell ? CLC_B : CLW_B;
  const int sr = cell ? CLC_R : CLW_R, sr0 = cell ? CLC_R0 : CLW_R0;
  const int sp = cell ? CLC_P : CLW_P, sv = cell ? CLC_V : CLW_V;
  const int ss = cell ? CLC_S : CLW_S, sq = cell ? CLC_Q : CLW_Q;
  const int sy = cell ? CLC_Y : CLW_Y, sz = cell ? CLC_Z : CLW_Z;

  // w(dst) = sa w(a) + sb w(b)
  auto axpby = [&](const int dst, const int a, const Real sa, const int b,
                   const Real sb_) {
    par_for("cln_axpby", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,dst+c,k,j,i) = sa*w(m,a+c,k,j,i) + sb_*w(m,b+c,k,j,i);
      }
    });
  };
  auto apply = [&](const int src, const int dst) {
    if (cell) {
      CleanCellApply(src, dst);
    } else {
      CleanApply(src, dst);
    }
  };
  auto precondition = [&](const int src, const int dst) {
    if (cell) {
      CleanCellPrecondition(src, dst);
    } else {
      CleanPrecondition(src, dst);
    }
  };
  auto dots = [&](const int n, const std::array<int, 4> ab) {
    return CleanDots(w, nc, n, ab);
  };

  // (re)start from the true residual of X with the shadow residual r0 = r.  A breakdown
  // (r0.r at the rounding level of |r0| |r|, or r0.v or omega zero) restarts; a second
  // breakdown right after a restart fails the solve.  The dots of an iteration take four
  // reductions: r0.v, s.s, (q.q, q.s), and (r.r, r0.r) for the next one.
  Real rho = 1.0, alpha = 1.0, omega = 1.0, rr = 0.0, r0r0 = 0.0, rho_new = 0.0, bb = 0.0;
  int it0 = 0;
  auto restart = [&](const int it) {
    if (it == it0 + 1) return false;
    apply(sx, sr);
    axpby(sr, sb, 1.0, sr, -1.0);
    axpby(sr0, sr, 1.0, sr, 0.0);
    axpby(sp, sr, 0.0, sr, 0.0);
    axpby(sv, sr, 0.0, sr, 0.0);
    rho = alpha = omega = 1.0;
    const auto d = dots((it < 0) ? 2 : 1, {sr, sr, sb, sb});
    rr = r0r0 = rho_new = d[0];
    if (it < 0) bb = d[1];
    it0 = it;
    return true;
  };
  restart(-1);
  const Real tol = cln_tol*fmax(std::sqrt(bb), 1.0e-300);
  if (std::sqrt(rr) <= tol) return 0;

  const Real eps = std::numeric_limits<Real>::epsilon();
  for (int it = 1; it <= cln_maxit; ++it) {
    if (!(fabs(rho_new) > eps*std::sqrt(r0r0*rr))) {
      if (restart(it)) continue;
      return -it;
    }
    const Real beta = (rho_new/rho)*(alpha/omega);
    rho = rho_new;
    par_for("cln_pupdate", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,sp+c,k,j,i) = w(m,sr+c,k,j,i)
                          + beta*(w(m,sp+c,k,j,i) - omega*w(m,sv+c,k,j,i));
      }
    });
    precondition(sp, sy);
    apply(sy, sv);
    const Real r0v = dots(1, {sr0, sv, 0, 0})[0];
    if (!std::isfinite(r0v) || r0v == 0.0) {
      if (restart(it)) continue;
      return -it;
    }
    alpha = rho/r0v;
    axpby(ss, sr, 1.0, sv, -alpha);
    if (std::sqrt(dots(1, {ss, ss, 0, 0})[0]) <= tol) {
      axpby(sx, sx, 1.0, sy, alpha);
      return it;
    }
    precondition(ss, sz);
    apply(sz, sq);
    const auto d = dots(2, {sq, sq, sq, ss});
    omega = (d[0] > 0.0) ? d[1]/d[0] : 0.0;
    par_for("cln_xupdate", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,sx+c,k,j,i) += alpha*w(m,sy+c,k,j,i) + omega*w(m,sz+c,k,j,i);
        w(m,sr+c,k,j,i) = w(m,ss+c,k,j,i) - omega*w(m,sq+c,k,j,i);
      }
    });
    const auto e = dots(2, {sr, sr, sr0, sr});
    rr = e[0];
    rho_new = e[1];
    if (std::sqrt(rr) <= tol) return it;
    if (omega == 0.0 && !restart(it)) return -it;
  }
  return -cln_maxit;
}

}  // namespace radm1
