//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_ringfourier.cpp
//! \brief ringsmooth-1102 <mesh>/polar_ring_avg_filter = fourier: the conservative Fourier
//! truncation of the polar rings that replaces the piecewise-constant group average of
//! ringavg-1029 (see hydro_ringavg.hpp for the keys and the energy rule).
//!
//! WHY.  The group average leaves a phi STAIRCASE in every averaged row, while the
//! well-balanced tangential force g_t (roche_wb_tan: the fine operator's own residual of
//! the smooth tidal IC) and the radiation keep fine-grid smooth fields: the staircase
//! state is never the discrete equilibrium, and the polar rows of the sp accretor grow a
//! v_phi pattern at the tidal m = 1, 2 (accsp_1102 test: lateral KE +41 % at t 0.005,
//! confined to the averaged rows).  The truncation leaves the low modes exactly alone.
//!
//! HOW.  Rings are global: (pole, row from the pole, global i).  For every ring and
//! conserved variable (IEN replaced by the internal energy e = E - KE [- rho Phi]),
//! a_m = sum_k u_k cos(2 pi m k/nphi), b_m = sum_k u_k sin(...), m = 0..m_c: each rank
//! sums its own phi blocks in a fixed order (deterministic), then MPI_Allreduce (SUM)
//! on the host.  u_f,k = (a_0 + 2 sum_{m>=1} (a_m cos + b_m sin))/nphi (m_c < nphi/2).
//! Then per ring: sum E, sum E_f (E_f = e_f + KE(rho_f, m_f) [+ rho_f Phi]) and
//! min(rho_f - thr rho), min(e_f - thr e) (second small Allreduce, SUM and MIN); a ring
//! passes when both mins (the second with the uniform share dE/nphi added when dE is gas
//! heat) are >= 0, and gets e_f + dE/nphi; else the old group average on that ring.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro_ringavg.hpp"

namespace hydro {

RingAvg::~RingAvg() {
  if (fourier && global_variable::my_rank == 0 && ncall > 0.0) {
    std::cout << "RINGAVG fourier: " << ncall << " filter calls, " << nfall
              << " ring-i fallbacks to the group average" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn RingAvg::FourierInit
//! \brief rings, cutoffs, slots (vjp: global row jg of every entry)

void RingAvg::FourierInit(ParameterInput *pin, const std::vector<int> &vjg) {
  Mesh *pm = pmy_pack->pmesh;
  auto &ic = pm->mb_indcs;
  const int nrow = pm->mesh_indcs.nx2;
  nphi = pm->mesh_indcs.nx3;
  nx1g = pm->mesh_indcs.nx1;
  const Real mfac = pin->GetOrAddReal("mesh", "polar_ring_avg_mfac", 1.0);
  fthr = pin->GetOrAddReal("mesh", "polar_ring_avg_thr", 0.5);
  // the group size of every global row (identical on all ranks: same rule, same x2v);
  // rings are the rows with N > 1, from the global nfac table rebuilt here
  std::vector<int> ng(nrow, 1);
  {
    auto hfac = Kokkos::create_mirror_view_and_copy(HostMemSpace(), nfac);
    std::vector<int> loc(nrow, 0);
    for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
      const int lx2 = static_cast<int>(pm->lloc_eachmb[pmy_pack->gids + m].lx2);
      for (int j = ic.js; j <= ic.je; ++j) {
        loc[lx2*ic.nx2 + (j - ic.js)] = static_cast<int>(hfac(m,j));
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(loc.data(), ng.data(), nrow, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
    ng = loc;
#endif
  }
  std::vector<int> rid(nrow, -1), hmc, hoff;
  nring = 0;
  for (int jg = 0; jg < nrow; ++jg) {
    if (ng[jg] > 1) {
      rid[jg] = nring++;
      const Real n = static_cast<Real>(ng[jg]);
      int mc = static_cast<int>(std::floor(mfac*(nphi/M_PI)*std::asin(1.0/n)));
      mc = std::max(0, std::min(mc, nphi/2 - 1));
      hmc.push_back(mc);
    }
  }
  // per entry: ring, global k0, global i0
  std::vector<int> her(nent), hk0(nent), hi0(nent);
  auto hm = Kokkos::create_mirror_view_and_copy(HostMemSpace(), em);
  for (int e = 0; e < nent; ++e) {
    const int gid = pmy_pack->gids + hm(e);
    her[e] = rid[vjg[e]];
    hk0[e] = static_cast<int>(pm->lloc_eachmb[gid].lx3)*ic.nx3;
    hi0[e] = static_cast<int>(pm->lloc_eachmb[gid].lx1)*ic.nx1;
  }
  // slots: (ring, i0) with the local entries in entry order
  std::vector<int> hsr, hsi, hsb, hsl;
  std::vector<int> done(nent, 0);
  for (int e = 0; e < nent; ++e) {
    if (done[e]) continue;
    hsr.push_back(her[e]); hsi.push_back(hi0[e]); hsb.push_back(static_cast<int>(hsl.size()));
    for (int f = e; f < nent; ++f) {
      if (!done[f] && her[f] == her[e] && hi0[f] == hi0[e]) {hsl.push_back(f); done[f] = 1;}
    }
  }
  nslot = static_cast<int>(hsr.size());
  hsb.push_back(static_cast<int>(hsl.size()));
  auto up = [](DvceArray1D<int> &d, const std::vector<int> &h, const char *nm) {
    d = DvceArray1D<int>(std::string(nm), static_cast<int>(std::max<size_t>(h.size(), 1)));
    auto hh = Kokkos::create_mirror_view(d);
    for (size_t q = 0; q < h.size(); ++q) hh(q) = h[q];
    Kokkos::deep_copy(d, hh);
  };
  up(ering, her, "rf_ering"); up(ek0, hk0, "rf_ek0"); up(eig0, hi0, "rf_eig0");
  up(rmc, hmc, "rf_rmc");
  up(sring, hsr, "rf_sring"); up(sig0, hsi, "rf_sig0"); up(sbeg, hsb, "rf_sbeg");
  up(slist, hsl, "rf_slist");
  tcos = DvceArray1D<Real>("rf_tcos", nphi);
  tsin = DvceArray1D<Real>("rf_tsin", nphi);
  auto hc = Kokkos::create_mirror_view(tcos);
  auto hs = Kokkos::create_mirror_view(tsin);
  for (int q = 0; q < nphi; ++q) {
    hc(q) = std::cos(2.0*M_PI*q/nphi);
    hs(q) = std::sin(2.0*M_PI*q/nphi);
  }
  fmcmax = 0;
  for (int r = 0; r < nring; ++r) fmcmax = std::max(fmcmax, hmc[r]);
  Kokkos::deep_copy(tcos, hc);
  Kokkos::deep_copy(tsin, hs);
  if (global_variable::my_rank == 0) {
    std::cout << "RINGAVG polar_ring_avg_filter = fourier (ringsmooth-1102): " << nring
              << " rings, cutoff m_c per row from the pole:";
    for (int r = 0; r < nring/2 + (nring % 2); ++r) std::cout << " " << hmc[r];
    std::cout << " (nphi " << nphi << ", mfac " << mfac << ", thr " << fthr
              << "); energy defect -> uniform heat" << std::endl;
  }
}

void RingAvg::FourierAlloc(const int nvar) {
  fnvar = nvar;
  auto hmc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rmc);
  auto ho = std::vector<int>(nring);
  int off = 0;
  for (int r = 0; r < nring; ++r) {ho[r] = off; off += nx1g*nvar*(hmc(r) + 1)*2;}
  nfoff = off;
  roff = DvceArray1D<int>("rf_roff", std::max(nring, 1));
  auto hr = Kokkos::create_mirror_view(roff);
  for (int r = 0; r < nring; ++r) hr(r) = ho[r];
  Kokkos::deep_copy(roff, hr);
  coef = DualArray1D<Real>("rf_coef", std::max(nfoff, 1));
  rsum = DualArray1D<Real>("rf_rsum", std::max(nring*nx1g*2, 1));
  rmin = DualArray1D<Real>("rf_rmin", std::max(nring*nx1g*2, 1));
  auto &ic = pmy_pack->pmesh->mb_indcs;
  scr = DvceArray4D<Real>("rf_scr", std::max(nent*nvar, 1), ic.nx3, 1, ic.nx1);
  part = DvceArray3D<Real>("rf_part", std::max(nent, 1), ic.nx1, 4);
}

//----------------------------------------------------------------------------------------
//! \fn RingAvg::Filter
//! \brief the conservative ring Fourier truncation (file header).  Collective: every rank
//! calls it (rings span ranks).

void RingAvg::Filter(DvceArray5D<Real> &u, const int nvar, const DvceArray4D<Real> &phi,
                     const bool grav, DvceArray5D<Real> *urad, const int irad,
                     const Real hc) {
  if (nring == 0) return;
  if (fnvar != nvar) FourierAlloc(nvar);
  auto &ic = pmy_pack->pmesh->mb_indcs;
  const int is = ic.is, ks = ic.ks, nx1 = ic.nx1, nx3 = ic.nx3;
  const int np = nphi, nxg = nx1g;
  auto u_ = u;
  auto phi_ = phi;
  auto em_ = em, ej_ = ej, en_ = en;
  auto er_ = ering, ek0_ = ek0, ei0_ = eig0, mc_ = rmc, ro_ = roff;
  auto sr_ = sring, si_ = sig0, sb_ = sbeg, sl_ = slist;
  auto tc_ = tcos, ts_ = tsin;
  auto cf = coef.d_view;
  const bool gr = grav;

  // (1) partial DFT sums of this rank's blocks, per (slot, i, v, m)
  Kokkos::deep_copy(cf, 0.0);
  const int mcmax = fmcmax;
  if (nslot > 0) {
    par_for("ringf_dft", DevExeSpace(), 0, nslot-1, 0, nx1-1, 0, nvar-1, 0, mcmax,
    KOKKOS_LAMBDA(const int s, const int i, const int v, const int mm) {
      const int r = sr_(s);
      const int mc = mc_(r);
      if (mm > mc) return;
      Real a = 0.0, b = 0.0;
      for (int q = sb_(s); q < sb_(s+1); ++q) {
        const int e = sl_(q);
        const int m = em_(e), j = ej_(e), k0 = ek0_(e);
        for (int kk = 0; kk < nx3; ++kk) {
          const int k = ks + kk, ii = is + i;
          Real x = u_(m,v,k,j,ii);
          if (v == IEN) {
            const Real d = u_(m,IDN,k,j,ii);
            const Real m1 = u_(m,IM1,k,j,ii), m2 = u_(m,IM2,k,j,ii), m3 = u_(m,IM3,k,j,ii);
            x -= 0.5*(m1*m1 + m2*m2 + m3*m3)/d;
            if (gr) x -= d*phi_(m,k,j,ii);
          }
          const int t = static_cast<int>((static_cast<long>(mm)*(k0 + kk)) % np);
          a += x*tc_(t);
          b += x*ts_(t);
        }
      }
      const int o = ro_(r) + (((si_(s) + i)*nvar + v)*(mc + 1) + mm)*2;
      cf(o) = a;
      cf(o+1) = b;
    });
  }
  coef.modify_device();
  coef.sync_host();
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, coef.h_view.data(), nfoff, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  coef.modify_host();
  coef.sync_device();

  // (2) the filtered state into scr, per-(entry, i) sums and mins
  auto sc = scr;
  auto pt = part;
  const Real thr = fthr;
  if (nent > 0) {
    par_for("ringf_syn", DevExeSpace(), 0, nent-1, 0, nx1-1,
    KOKKOS_LAMBDA(const int e, const int i) {
      const int m = em_(e), j = ej_(e), k0 = ek0_(e), r = er_(e), ig = ei0_(e) + i;
      const int mc = mc_(r), ii = is + i;
      const Real rn = 1.0/static_cast<Real>(np);
      Real se = 0.0, sf = 0.0, md = 1.0e300, me = 1.0e300;
      for (int kk = 0; kk < nx3; ++kk) {
        const int k = ks + kk;
        for (int v = 0; v < nvar; ++v) {
          const int o = ro_(r) + ((ig*nvar + v)*(mc + 1))*2;
          Real x = cf(o);
          for (int mm = 1; mm <= mc; ++mm) {
            const int t = static_cast<int>((static_cast<long>(mm)*(k0 + kk)) % np);
            x += 2.0*(cf(o + 2*mm)*tc_(t) + cf(o + 2*mm + 1)*ts_(t));
          }
          sc(e*nvar + v, kk, 0, i) = x*rn;
        }
        const Real d = u_(m,IDN,k,j,ii);
        const Real m1 = u_(m,IM1,k,j,ii), m2 = u_(m,IM2,k,j,ii), m3 = u_(m,IM3,k,j,ii);
        const Real ph = gr ? phi_(m,k,j,ii) : 0.0;
        const Real eo = u_(m,IEN,k,j,ii) - 0.5*(m1*m1 + m2*m2 + m3*m3)/d - d*ph;
        const Real df = sc(e*nvar + IDN, kk, 0, i);
        const Real f1 = sc(e*nvar + IM1, kk, 0, i), f2 = sc(e*nvar + IM2, kk, 0, i);
        const Real f3 = sc(e*nvar + IM3, kk, 0, i), ef = sc(e*nvar + IEN, kk, 0, i);
        se += u_(m,IEN,k,j,ii);
        sf += ef + ((df > 0.0) ? 0.5*(f1*f1 + f2*f2 + f3*f3)/df : 0.0) + df*ph;
        md = fmin(md, df - thr*d);
        me = fmin(me, ef - thr*eo);
      }
      pt(e,i,0) = se; pt(e,i,1) = sf; pt(e,i,2) = md; pt(e,i,3) = me;
    });
  }
  // slot reduction (fixed entry order), then the ring Allreduce (SUM and MIN)
  auto rs = rsum.d_view, rmn = rmin.d_view;
  const int nr = nring;
  Kokkos::deep_copy(rs, 0.0);
  Kokkos::deep_copy(rmn, 1.0e300);
  if (nslot > 0) {
    par_for("ringf_red", DevExeSpace(), 0, nslot-1, 0, nx1-1,
    KOKKOS_LAMBDA(const int s, const int i) {
      const int r = sr_(s), ig = si_(s) + i;
      Real a = 0.0, b = 0.0, c = 1.0e300, d = 1.0e300;
      for (int q = sb_(s); q < sb_(s+1); ++q) {
        const int e = sl_(q);
        a += pt(e,i,0); b += pt(e,i,1);
        c = fmin(c, pt(e,i,2)); d = fmin(d, pt(e,i,3));
      }
      rs((r*nxg + ig)*2) = a; rs((r*nxg + ig)*2 + 1) = b;
      rmn((r*nxg + ig)*2) = c; rmn((r*nxg + ig)*2 + 1) = d;
    });
  }
  rsum.modify_device(); rsum.sync_host();
  rmin.modify_device(); rmin.sync_host();
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, rsum.h_view.data(), nr*nxg*2, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, rmin.h_view.data(), nr*nxg*2, MPI_ATHENA_REAL, MPI_MIN,
                MPI_COMM_WORLD);
#endif
  const bool torad = (urad != nullptr);
  {
    double nf = 0.0;
    for (int q = 0; q < nr*nxg; ++q) {
      const Real de = (rsum.h_view(2*q) - rsum.h_view(2*q+1))/static_cast<Real>(np);
      const bool ok = rmin.h_view(2*q) >= 0.0 &&
                      rmin.h_view(2*q+1) + (torad ? 0.0 : de) >= 0.0;
      if (!ok) nf += 1.0;
    }
    nfall += nf;
    ncall += 1.0;
  }
  rsum.modify_host(); rsum.sync_device();
  rmin.modify_host(); rmin.sync_device();

  // (3) write back: the filtered state + dE/nphi (gas heat or radiation), or the group
  // average where the ring failed
  if (nent == 0) return;
  DvceArray5D<Real> ur;
  if (torad) ur = *urad;
  const int ir = irad;
  const Real hcr = hc;
  par_for("ringf_put", DevExeSpace(), 0, nent-1, 0, nx1-1,
  KOKKOS_LAMBDA(const int e, const int i) {
    const int m = em_(e), j = ej_(e), r = er_(e), ig = ei0_(e) + i, ii = is + i;
    const int q = r*nxg + ig;
    const Real de = (rs(2*q) - rs(2*q+1))/static_cast<Real>(np);
    const bool ok = rmn(2*q) >= 0.0 && rmn(2*q+1) + (torad ? 0.0 : de) >= 0.0;
    if (ok) {
      for (int kk = 0; kk < nx3; ++kk) {
        const int k = ks + kk;
        for (int v = 0; v < nvar; ++v) {
          if (v != IEN) u_(m,v,k,j,ii) = sc(e*nvar + v, kk, 0, i);
        }
        const Real df = u_(m,IDN,k,j,ii);
        const Real f1 = u_(m,IM1,k,j,ii), f2 = u_(m,IM2,k,j,ii), f3 = u_(m,IM3,k,j,ii);
        Real en = sc(e*nvar + IEN, kk, 0, i) + 0.5*(f1*f1 + f2*f2 + f3*f3)/df;
        if (gr) en += df*phi_(m,k,j,ii);
        if (torad) {
          ur(m,ir,k,j,ii) += hcr*de;
        } else {
          en += de;
        }
        u_(m,IEN,k,j,ii) = en;
      }
    } else {
      // the ringavg-1029 group average (hydro::RingAvg::Average), etotgrav included
      const int n = en_(e);
      const Real rn = 1.0/static_cast<Real>(n);
      for (int k = ks; k < ks + nx3; k += n) {
        for (int v = 0; v < nvar; ++v) {
          Real s = 0.0;
          for (int kk = k; kk < k + n; ++kk) {s += u_(m,v,kk,j,ii);}
          s *= rn;
          for (int kk = k; kk < k + n; ++kk) {u_(m,v,kk,j,ii) = s;}
        }
        if (gr) {
          Real pb = 0.0;
          for (int kk = k; kk < k + n; ++kk) {pb += phi_(m,kk,j,ii);}
          pb *= rn;
          const Real db = u_(m,IDN,k,j,ii);
          for (int kk = k; kk < k + n; ++kk) {
            u_(m,IEN,kk,j,ii) += db*(phi_(m,kk,j,ii) - pb);
          }
        }
      }
    }
  });
}

}  // namespace hydro
