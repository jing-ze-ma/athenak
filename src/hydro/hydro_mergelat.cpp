//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_mergelat.cpp
//! \brief mrglat-1019: merged lateral cells near r_in (see hydro_mergelat.hpp).

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_mergelat.hpp"

namespace hydro {

namespace {
std::vector<std::string> MlatSplit(const std::string &s) {
  std::vector<std::string> out;
  std::string tok;
  // ',' or ':' (sbatch --export splits on commas)
  std::string s2 = s;
  for (auto &c : s2) {if (c == ':') c = ',';}
  std::stringstream ss2(s2);
  while (std::getline(ss2, tok, ',')) {
    const size_t a = tok.find_first_not_of(" \t");
    const size_t b = tok.find_last_not_of(" \t");
    if (a != std::string::npos) out.push_back(tok.substr(a, b - a + 1));
  }
  return out;
}
[[noreturn]] void MlatFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<hydro>/merge_lat: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
}  // namespace

bool MergeLat::Named(ParameterInput *pin) {
  return pin->DoesParameterExist("hydro", "merge_lat_rc") ||
         pin->DoesParameterExist("hydro", "merge_lat_diag");
}

MergeLat::MergeLat(MeshBlockPack *ppack, ParameterInput *pin) : pmy_pack(ppack) {
  for (int n = 0; n < kMlatMaxShell; ++n) {rc[n] = 0.0; fs[n] = 1;}
  if (pin->DoesParameterExist("hydro", "merge_lat_diag")) {
    diag_every = pin->GetInteger("hydro", "merge_lat_diag");
  }
  if (!pin->DoesParameterExist("hydro", "merge_lat_rc")) return;
  on = true;
  auto vr = MlatSplit(pin->GetString("hydro", "merge_lat_rc"));
  auto vf = MlatSplit(pin->GetString("hydro", "merge_lat_f"));
  if (vr.empty() || vr.size() != vf.size() ||
      static_cast<int>(vr.size()) > kMlatMaxShell) {
    MlatFatal("merge_lat_rc and merge_lat_f need the same number (1..8) of entries");
  }
  nsh = static_cast<int>(vr.size());
  Mesh *pm = ppack->pmesh;
  auto &ic = pm->mb_indcs;
  for (int n = 0; n < nsh; ++n) {
    rc[n] = std::stod(vr[n]);
    fs[n] = std::stoi(vf[n]);
    if (n > 0 && !(rc[n] > rc[n-1])) MlatFatal("merge_lat_rc must increase");
    if (fs[n] < 1 || fs[n] > 32) MlatFatal("merge_lat_f must be in 1..32");
    if ((ic.nx2 % fs[n]) != 0 || (ic.nx3 % fs[n]) != 0) {
      MlatFatal("MeshBlock nx2 = " + std::to_string(ic.nx2) + ", nx3 = " +
                std::to_string(ic.nx3) + " not divisible by merge_lat_f = " +
                std::to_string(fs[n]) + " (groups may not straddle a MeshBlock)");
    }
    fbig = std::max(fbig, fs[n]);
  }
  if (!(pm->use_spherical_polar || pm->use_cubed_sphere)) {
    MlatFatal("needs use_spherical_polar or use_cubed_sphere");
  }
  if (!pm->three_d) MlatFatal("needs a 3-D mesh");
  if (pm->multilevel) MlatFatal("SMR/AMR not supported");
  const std::string mm = pin->GetOrAddString("hydro", "merge_lat_mom", "rigid");
  if (mm == "rigid") {
    rigid = true;
  } else if (mm == "cart") {
    rigid = false;
  } else {
    MlatFatal("merge_lat_mom = rigid | cart");
  }
  rad = pin->GetOrAddBoolean("hydro", "merge_lat_rad", false);
  if (pin->DoesParameterExist("hydro", "merge_lat_check")) {
    check = pin->GetBoolean("hydro", "merge_lat_check");
  }
  if (global_variable::my_rank == 0) {
    std::cout << "MLAT merged lateral cells: " << nsh << " shell(s):";
    for (int n = 0; n < nsh; ++n) {
      std::cout << " x1v < " << rc[n] << " f " << fs[n] << ";";
    }
    std::cout << " mom " << (rigid ? "rigid" : "cart") << ", rad "
              << (rad ? "merged" : "fine") << std::endl;
  }
}

MergeLatGeo MergeLat::Geo() {
  MergeLatGeo g;
  Coordinates *pco = pmy_pack->pcoord;
  auto &ic = pmy_pack->pmesh->mb_indcs;
  g.vol = pco->volume;
  g.x1v = pco->x1v;  g.x2v = pco->x2v;  g.x3v = pco->x3v;
  g.ccell = pco->cos_cell;
  g.size = pmy_pack->pmb->mb_size;
  g.panel = pmy_pack->pmb->mb_panel;
  g.js = ic.js;  g.ks = ic.ks;  g.nx2 = ic.nx2;  g.nx3 = ic.nx3;
  g.cs = pmy_pack->pmesh->use_cubed_sphere;
  return g;
}

void MergeLat::Build() {
  if (built) return;
  built = true;
  auto &ic = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int n1 = ic.nx1 + 2*ic.ng;
  Kokkos::realloc(fi, nmb, n1);
  auto fi_ = fi;
  auto x1v_ = pmy_pack->pcoord->x1v;
  const int nsh_ = on ? nsh : 0;
  Real rc0 = rc[0], rc1 = rc[1], rc2 = rc[2], rc3 = rc[3];
  Real rc4 = rc[4], rc5 = rc[5], rc6 = rc[6], rc7 = rc[7];
  int f0 = fs[0], f1 = fs[1], f2 = fs[2], f3 = fs[3];
  int f4 = fs[4], f5 = fs[5], f6 = fs[6], f7 = fs[7];
  par_for("mlat_fi", DevExeSpace(), 0, nmb-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int i) {
    const Real r = x1v_(m,i);
    const Real rr[8] = {rc0, rc1, rc2, rc3, rc4, rc5, rc6, rc7};
    const int ff[8] = {f0, f1, f2, f3, f4, f5, f6, f7};
    int f = 1;
    for (int n = 0; n < nsh_; ++n) {
      if (r < rr[n]) {f = ff[n]; break;}
    }
    fi_(m,i) = f;
  });
  // the largest merged interior i on this rank
  auto hfi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fi);
  imax_loc = -1;
  for (int m = 0; m < nmb; ++m) {
    for (int i = ic.is; i <= ic.ie; ++i) {
      if (hfi(m,i) > 1) imax_loc = std::max(imax_loc, i);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn MergeLat::Average
//! \brief replace every merged group's conserved state by its group value (see the file
//! header of hydro_mergelat.hpp).  One thread per group (the group's first cell), sums in
//! a fixed order: deterministic.

void MergeLat::Average(DvceArray5D<Real> &u, const int nhyd, const int nvar,
                       const int kind, const DvceArray4D<Real> &phi, const bool grav) {
  if (!on) return;
  if (nvar - nhyd > 4) MlatFatal("at most 4 passive scalars");
  Build();
  if (imax_loc < 0) return;
  ++ncall;
  auto &ic = pmy_pack->pmesh->mb_indcs;
  const int is = ic.is, js = ic.js, je = ic.je, ks = ic.ks, ke = ic.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const MergeLatGeo g = Geo();
  auto fi_ = fi;
  auto u_ = u;
  const bool rigid_ = rigid;
  const int nhyd_ = nhyd, nvar_ = nvar;
  auto phi_ = phi;
  const bool grav_ = grav;
  DvceArray1D<int> nfb("mlat_nfb", 1);
  const bool chk_ = check && (kind == 0);
  DvceArray1D<Real> cerr("mlat_cerr", 4);
  par_for("mlat_avg", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, imax_loc,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int f = fi_(m,i);
    if (f <= 1) return;
    if (((j - js) % f) != 0 || ((k - ks) % f) != 0) return;
    if (kind == 1) {
      // radiation: E volume average, F the Cartesian volume average
      Real vs = 0.0, es = 0.0, fc[3] = {0.0, 0.0, 0.0};
      for (int kk = k; kk < k + f; ++kk) {
        for (int jj = j; jj < j + f; ++jj) {
          Real rh[3], b1[3], b2[3], cc, v[3];
          MlatBasis(g, m, kk, jj, rh, b1, b2, cc);
          const Real dv = g.vol(m,kk,jj,i);
          MlatToCart(rh, b1, b2, cc, u_(m,1,kk,jj,i), u_(m,2,kk,jj,i), u_(m,3,kk,jj,i), v);
          vs += dv;
          es += dv*u_(m,0,kk,jj,i);
          for (int q = 0; q < 3; ++q) {fc[q] += dv*v[q];}
        }
      }
      const Real iv = 1.0/vs;
      const Real eb = es*iv;
      for (int q = 0; q < 3; ++q) {fc[q] *= iv;}
      for (int kk = k; kk < k + f; ++kk) {
        for (int jj = j; jj < j + f; ++jj) {
          Real rh[3], b1[3], b2[3], cc;
          MlatBasis(g, m, kk, jj, rh, b1, b2, cc);
          u_(m,0,kk,jj,i) = eb;
          MlatFromCart(rh, b1, b2, fc, u_(m,1,kk,jj,i), u_(m,2,kk,jj,i), u_(m,3,kk,jj,i));
        }
      }
      return;
    }
    // ---- hydro
    const Real r = g.x1v(m,i);
    Real vs = 0.0, ms = 0.0, es = 0.0, ss[4] = {0.0, 0.0, 0.0, 0.0}, vp = 0.0;
    Real P[3] = {0.0, 0.0, 0.0}, X[3] = {0.0, 0.0, 0.0}, Lo[3] = {0.0, 0.0, 0.0};
    const int nsc = nvar_ - nhyd_;
    for (int kk = k; kk < k + f; ++kk) {
      for (int jj = j; jj < j + f; ++jj) {
        Real rh[3], b1[3], b2[3], cc, v[3];
        MlatBasis(g, m, kk, jj, rh, b1, b2, cc);
        const Real dv = g.vol(m,kk,jj,i);
        MlatToCart(rh, b1, b2, cc, u_(m,IM1,kk,jj,i), u_(m,IM2,kk,jj,i),
                   u_(m,IM3,kk,jj,i), v);
        vs += dv;
        ms += dv*u_(m,IDN,kk,jj,i);
        es += dv*u_(m,IEN,kk,jj,i);
        if (grav_) vp += dv*phi_(m,kk,jj,i);
        for (int n = 0; n < nsc && n < 4; ++n) {ss[n] += dv*u_(m,nhyd_+n,kk,jj,i);}
        for (int q = 0; q < 3; ++q) {P[q] += dv*v[q]; X[q] += dv*r*rh[q];}
        if (chk_) {
          Lo[0] += dv*r*(rh[1]*v[2] - rh[2]*v[1]);
          Lo[1] += dv*r*(rh[2]*v[0] - rh[0]*v[2]);
          Lo[2] += dv*r*(rh[0]*v[1] - rh[1]*v[0]);
        }
      }
    }
    const Real iv = 1.0/vs;
    const Real db = ms*iv;
    Real xc[3], mb[3];
    for (int q = 0; q < 3; ++q) {xc[q] = X[q]*iv; mb[q] = P[q]*iv;}
    Real om[3] = {0.0, 0.0, 0.0};
    Real ekin = 0.0;   // group kinetic energy (x V) of the new momentum field
    // etotgrav: IEN carries rho Phi; the group potential energy of the new state
    const Real epot = db*vp;
    if (db > 0.0) {
      const Real pp = P[0]*P[0] + P[1]*P[1] + P[2]*P[2];
      ekin = 0.5*pp*iv/db;
      if (rigid_) {
        // angular momentum about the centroid and the inertia tensor of the cells
        Real L[3] = {0.0, 0.0, 0.0};
        Real J[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        for (int kk = k; kk < k + f; ++kk) {
          for (int jj = j; jj < j + f; ++jj) {
            Real rh[3], b1[3], b2[3], cc, v[3], d[3];
            MlatBasis(g, m, kk, jj, rh, b1, b2, cc);
            const Real dv = g.vol(m,kk,jj,i);
            MlatToCart(rh, b1, b2, cc, u_(m,IM1,kk,jj,i), u_(m,IM2,kk,jj,i),
                       u_(m,IM3,kk,jj,i), v);
            for (int q = 0; q < 3; ++q) {d[q] = r*rh[q] - xc[q];}
            L[0] += dv*(d[1]*v[2] - d[2]*v[1]);
            L[1] += dv*(d[2]*v[0] - d[0]*v[2]);
            L[2] += dv*(d[0]*v[1] - d[1]*v[0]);
            const Real d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            for (int a = 0; a < 3; ++a) {
              for (int b = 0; b < 3; ++b) {
                J[a][b] += dv*((a == b ? d2 : 0.0) - d[a]*d[b]);
              }
            }
          }
        }
        // Omega = J^-1 L (cofactors; J symmetric positive definite for f >= 2)
        const Real c00 = J[1][1]*J[2][2] - J[1][2]*J[2][1];
        const Real c01 = J[1][2]*J[2][0] - J[1][0]*J[2][2];
        const Real c02 = J[1][0]*J[2][1] - J[1][1]*J[2][0];
        const Real det = J[0][0]*c00 + J[0][1]*c01 + J[0][2]*c02;
        const Real tr = J[0][0] + J[1][1] + J[2][2];
        if (det > 1.0e-12*tr*tr*tr) {
          const Real c11 = J[0][0]*J[2][2] - J[0][2]*J[2][0];
          const Real c12 = J[0][1]*J[2][0] - J[0][0]*J[2][1];
          const Real c22 = J[0][0]*J[1][1] - J[0][1]*J[1][0];
          const Real id = 1.0/det;
          // inverse = cofactor^T / det (symmetric)
          om[0] = (c00*L[0] + c01*L[1] + c02*L[2])*id;
          om[1] = (c01*L[0] + c11*L[1] + c12*L[2])*id;
          om[2] = (c02*L[0] + c12*L[1] + c22*L[2])*id;
          const Real erot = 0.5*(om[0]*L[0] + om[1]*L[1] + om[2]*L[2])/db;
          // keep at least half of the group's thermal(+potential) energy
          if (erot >= 0.0 && erot <= 0.5*(es - ekin - epot)) {
            ekin += erot;
          } else {
            om[0] = om[1] = om[2] = 0.0;
            Kokkos::atomic_add(&nfb(0), 1);
          }
        }
      }
    }
    const Real eb = (es - ekin - epot)*iv;   // thermal energy density
    Real sb[4];
    for (int n = 0; n < 4; ++n) {sb[n] = ss[n]*iv;}
    for (int kk = k; kk < k + f; ++kk) {
      for (int jj = j; jj < j + f; ++jj) {
        Real rh[3], b1[3], b2[3], cc, v[3], d[3];
        MlatBasis(g, m, kk, jj, rh, b1, b2, cc);
        for (int q = 0; q < 3; ++q) {d[q] = r*rh[q] - xc[q];}
        v[0] = mb[0] + (om[1]*d[2] - om[2]*d[1]);
        v[1] = mb[1] + (om[2]*d[0] - om[0]*d[2]);
        v[2] = mb[2] + (om[0]*d[1] - om[1]*d[0]);
        u_(m,IDN,kk,jj,i) = db;
        MlatFromCart(rh, b1, b2, v, u_(m,IM1,kk,jj,i), u_(m,IM2,kk,jj,i),
                     u_(m,IM3,kk,jj,i));
        const Real ke = (db > 0.0) ? 0.5*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/db : 0.0;
        u_(m,IEN,kk,jj,i) = eb + ke + (grav_ ? db*phi_(m,kk,jj,i) : 0.0);
        for (int n = 0; n < nsc && n < 4; ++n) {u_(m,nhyd_+n,kk,jj,i) = sb[n];}
      }
    }
    if (chk_) {
      // <hydro>/merge_lat_check: the group sums after the replacement vs before
      Real m2 = 0.0, e2 = 0.0, p2[3] = {0.0, 0.0, 0.0}, l2[3] = {0.0, 0.0, 0.0};
      for (int kk = k; kk < k + f; ++kk) {
        for (int jj = j; jj < j + f; ++jj) {
          Real rh[3], b1[3], b2[3], cc, v[3];
          MlatBasis(g, m, kk, jj, rh, b1, b2, cc);
          const Real dv = g.vol(m,kk,jj,i);
          MlatToCart(rh, b1, b2, cc, u_(m,IM1,kk,jj,i), u_(m,IM2,kk,jj,i),
                     u_(m,IM3,kk,jj,i), v);
          m2 += dv*u_(m,IDN,kk,jj,i);
          e2 += dv*u_(m,IEN,kk,jj,i);
          for (int q = 0; q < 3; ++q) {p2[q] += dv*v[q];}
          l2[0] += dv*r*(rh[1]*v[2] - rh[2]*v[1]);
          l2[1] += dv*r*(rh[2]*v[0] - rh[0]*v[2]);
          l2[2] += dv*r*(rh[0]*v[1] - rh[1]*v[0]);
        }
      }
      // momentum and L errors scaled by the group's |P| + M c-free scale: |m| sum
      Real pn = 0.0, ln = 0.0, dp = 0.0, dl = 0.0;
      for (int q = 0; q < 3; ++q) {
        pn += P[q]*P[q];  ln += Lo[q]*Lo[q];
        dp += (p2[q] - P[q])*(p2[q] - P[q]);  dl += (l2[q] - Lo[q])*(l2[q] - Lo[q]);
      }
      const Real pr = sqrt(dp/fmax(pn, 1.0e-300));
      const Real lr = sqrt(dl)/fmax(sqrt(ln) + r*sqrt(pn), 1.0e-300);
      Kokkos::atomic_max(&cerr(0), fabs(m2 - ms)/fmax(fabs(ms), 1.0e-300));
      Kokkos::atomic_max(&cerr(1), fabs(e2 - es)/fmax(fabs(es), 1.0e-300));
      Kokkos::atomic_max(&cerr(2), pr);
      Kokkos::atomic_max(&cerr(3), lr);
    }
  });
  if (chk_) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cerr);
    for (int n = 0; n < 4; ++n) {chkerr[n] = std::max(chkerr[n], static_cast<double>(h(n)));}
  }
  if (rigid && kind == 0) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), nfb);
    nfallback += h(0);
  }
}

int MergeLat::CountNonUniform(const DvceArray5D<Real> &u, const int n0) {
  if (!on) return 0;
  Build();
  int cnt = 0;
  if (imax_loc >= 0) {
    auto &ic = pmy_pack->pmesh->mb_indcs;
    const int is = ic.is, js = ic.js, je = ic.je, ks = ic.ks, ke = ic.ke;
    const int nmb = pmy_pack->nmb_thispack;
    const int ni = imax_loc - is + 1, nj = je - js + 1, nk = ke - ks + 1;
    const int nkji = nk*nj*ni, nji = nj*ni;
    auto fi_ = fi;
    auto u_ = u;
    Kokkos::parallel_reduce("mlat_uni", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            nmb*nkji),
    KOKKOS_LAMBDA(const int idx, int &c) {
      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji + ks;
      const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
      const int i = (idx - m*nkji - (k - ks)*nji - (j - js)*ni) + is;
      const int f = fi_(m,i);
      if (f <= 1 || ((j - js) % f) != 0 || ((k - ks) % f) != 0) return;
      const Real q0 = u_(m,n0,k,j,i);
      bool same = true;
      for (int kk = k; kk < k + f; ++kk) {
        for (int jj = j; jj < j + f; ++jj) {
          if (u_(m,n0,kk,jj,i) != q0) same = false;
        }
      }
      if (!same) c += 1;
    }, Kokkos::Sum<int>(cnt));
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  return cnt;
}

//----------------------------------------------------------------------------------------
//! \fn MergeLat::Diag
//! \brief global gas sums: M, E(IEN), P (Cartesian), L = sum V x x m (Cartesian, cell
//! centre arm).  Rank 0 prints one MLATDIAG line.

void MergeLat::Diag(const DvceArray5D<Real> &u, const char *tag) {
  Mesh *pm = pmy_pack->pmesh;
  auto &ic = pm->mb_indcs;
  const int is = ic.is, js = ic.js, ks = ic.ks;
  const int nx1 = ic.nx1, nx2 = ic.nx2, nx3 = ic.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const MergeLatGeo g = Geo();
  auto u_ = u;
  array_sum::GlobalSum s;
  Kokkos::parallel_reduce("mlat_diag", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &a) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = (idx - m*nkji - (k - ks)*nji - (j - js)*nx1) + is;
    Real rh[3], b1[3], b2[3], cc, v[3];
    MlatBasis(g, m, k, j, rh, b1, b2, cc);
    const Real dv = g.vol(m,k,j,i);
    const Real r = g.x1v(m,i);
    MlatToCart(rh, b1, b2, cc, u_(m,IM1,k,j,i), u_(m,IM2,k,j,i), u_(m,IM3,k,j,i), v);
    a.the_array[0] += dv*u_(m,IDN,k,j,i);
    a.the_array[1] += dv*u_(m,IEN,k,j,i);
    for (int q = 0; q < 3; ++q) {a.the_array[2+q] += dv*v[q];}
    const Real x[3] = {r*rh[0], r*rh[1], r*rh[2]};
    a.the_array[5] += dv*(x[1]*v[2] - x[2]*v[1]);
    a.the_array[6] += dv*(x[2]*v[0] - x[0]*v[2]);
    a.the_array[7] += dv*(x[0]*v[1] - x[1]*v[0]);
  }, Kokkos::Sum<array_sum::GlobalSum>(s));
  double h[8];
  for (int n = 0; n < 8; ++n) h[n] = s.the_array[n];
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, h, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  long nf = nfallback;
  MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  double ce[4] = {chkerr[0], chkerr[1], chkerr[2], chkerr[3]};
  MPI_Allreduce(MPI_IN_PLACE, ce, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  long nf = nfallback;
  double ce[4] = {chkerr[0], chkerr[1], chkerr[2], chkerr[3]};
#endif
  if (global_variable::my_rank == 0) {
    std::printf("MLATDIAG %s cyc %d t %.10e dt %.6e M %.16e E %.16e P %.10e %.10e %.10e "
                "L %.16e %.16e %.16e nfb %ld chk %.2e %.2e %.2e %.2e\n", tag, pm->ncycle,
                pm->time, pm->dt, h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], nf,
                ce[0], ce[1], ce[2], ce[3]);
    std::fflush(stdout);
  }
}

}  // namespace hydro
