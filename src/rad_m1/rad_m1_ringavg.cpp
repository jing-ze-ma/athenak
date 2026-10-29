//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_ringavg.cpp
//! \brief ringavg-1029 <rad_m1>/polar_ring_avg = true: the implicit M1 transport of the
//! polar rows on the COARSE phi groups of <mesh>/polar_ring_avg (hydro_ringavg.hpp, whose
//! group table it reads), i.e. each group of N_j cells is ONE radiation cell, consistent
//! with the averaged hydro.  Opt-in; default false: nothing below runs (bitwise).
//!
//! WHY.  Next to a pole the phi face is a small fraction of the equatorial one, so the
//! phi couplings of the implicit operator are ~(dr/(varpi dphi))^2 stiffer than the
//! radial ones there; they grow with dt (which the hydro averaging doubles) and the
//! x1-line preconditioners do not see them (S1 full sphere: BiCGStab up to lin_maxit).
//!
//! HOW (a Galerkin restriction of the SAME operator; no kernel of the assembly changes):
//!  (1) the solve is restricted to the group-uniform vectors: x = P y, the rows summed
//!      over each group, R A P y = R b (R = P^T; the cells of a row have one volume, so
//!      R/N is the group average).  BiCGStab runs on that system as is: b and x0 are
//!      averaged (RavKrylovStart), and the output of every operator and preconditioner
//!      application (ImplicitApplyOp, ImplicitPrecond(X)) is averaged, so every Krylov
//!      vector stays group-uniform and the inner products are those of the restricted
//!      system.  The rows are conservative (face flux in one row = face flux out of the
//!      neighbour row), so the faces INSIDE a group cancel exactly in R A P: the
//!      converged E is the coarse-cell solution and the energy is conserved by the sum.
//!  (2) the group geometry on the faces BETWEEN groups: the centre distance the phi
//!      gradient is taken over (Coordinates dxface.x3f, read by every implicit kernel
//!      through M1Dxf()) is N_j times the fine one there (the coarse cells' centre
//!      distance N_j varpi dphi), so the phi transport is the coarse grid's (diffusion
//!      coupling down by N_j^2 overall: N_j from the distance, N_j from R/N).  The
//!      operator, the right-hand side (lagged and explicit face terms), the theta
//!      limiter and the face-flux update all read the same dxface: one operator.
//!  (3) the preconditioner (mg, rbgs, pcr: the 7-point rows TA..CKP) keeps the FINE rows
//!      and its output is averaged: P R/N M^{-1} is a fixed linear map, a valid right
//!      preconditioner of the restricted system.  Study key polar_ring_avg_lump = true
//!      lumps the intra-group phi couplings onto the diagonal first (RavLump: TB += CKM,
//!      CKM = 0 on an interior face; exact on group-uniform vectors, so the system is
//!      unchanged); measured WORSE (S1 t 0.015-0.017: inner iterations mean 44.9 lumped,
//!      23.8 fine rows, 38.8 without the key): default false.
//!  (4) the persistent face state before the write-back (RavFaces): the radial faces
//!      f0x1 and the theta faces f0x2 averaged over the phi groups (theta face j over
//!      the smaller group of its two rows; block x2 edges are left alone); the phi faces
//!      INSIDE a group set to the linear interpolation of the group's two boundary
//!      faces, i.e. the coarse cell's flux (every fine cell then has the coarse cell's
//!      flux divergence, and the cell F_phi = (f_L + f_R)/2 of the coarse cell on
//!      average).  The write-back then forms the cell F and the gas kick from these.
//!  (5) after the solve (RadiationM1::HydroConToPrim): u0 (E, F1..F3, componentwise
//!      like the hydro momentum) and the gas the coupling wrote (hydro u0, as
//!      hydro::RingAvg::Average, etotgrav included) replaced by their group averages,
//!      so gas + radiation stay uniform per group; sums kept: conservative.
//! The closure (vet_col) and opacities are built from these group-uniform states (the
//! hydro user sources after its averaging can leave a small phi spread, as in the hydro).
//! The implicit M1 has no phi CFL (rad_m1_newdt.cpp), nothing to change there.
//!
//! Needs <mesh>/polar_ring_avg = true, transport = implicit, 3-D spherical polar,
//! implicit_solver = bicgstab, implicit_op_stencil = true and implicit_krylov_fuse <= 1
//! (the fused operator-reduction kernels of krylov_fuse 2 | 3 reduce before the
//! average).

#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_ringavg.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

namespace {
[[noreturn]] void RavM1Fatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_m1>/polar_ring_avg: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RavInit
//! \brief read <rad_m1>/polar_ring_avg, check the configuration, build the coarse-group
//! dxface (x3f times N_j on the faces between groups; x1f, x2f shared with Coordinates)

void RadiationM1::RavInit(ParameterInput *pin) {
  rav_on = pin->GetOrAddBoolean("rad_m1", "polar_ring_avg", false);
  if (!rav_on) return;
  hydro::Hydro *ph = pmy_pack->phydro;
  if (ph == nullptr || ph->pravg == nullptr) {
    RavM1Fatal("needs <mesh>/polar_ring_avg = true (the hydro group table)");
  }
  if (transport != M1_TRANSPORT_IMPLICIT || !trans_x3) {
    RavM1Fatal("needs transport = implicit on a 3-D mesh");
  }
  if (!sph_geom || cs_geom) RavM1Fatal("needs the spherical-polar mesh");
  if (impl_solver != M1_ISOLV_BICGSTAB) RavM1Fatal("needs implicit_solver = bicgstab");
  if (!impl_stencil) RavM1Fatal("needs implicit_op_stencil = true");
  if (impl_kfuse > 1) RavM1Fatal("needs implicit_krylov_fuse = 0 or 1");
  prav = ph->pravg;
  // study keys: (3) the lumped preconditioner rows, (2) the coarse phi distance
  rav_lump = pin->GetOrAddBoolean("rad_m1", "polar_ring_avg_lump", false);
  rav_dist = pin->GetOrAddBoolean("rad_m1", "polar_ring_avg_dist", true);

  auto *pc = pmy_pack->pcoord;
  rav_dxf = pc->dxface;   // shallow: x1f, x2f are Coordinates' own
  rav_dxf.x3f = DvceArray4D<Real>("m1ravdxf3", pc->dxface.x3f.extent(0),
                                  pc->dxface.x3f.extent(1), pc->dxface.x3f.extent(2),
                                  pc->dxface.x3f.extent(3));
  Kokkos::deep_copy(rav_dxf.x3f, pc->dxface.x3f);
  auto d3 = rav_dxf.x3f;
  auto nf = prav->nfac;
  const int ks = pmy_pack->pmesh->mb_indcs.ks;
  const int n3f = static_cast<int>(d3.extent(1)), n2 = static_cast<int>(d3.extent(2));
  const int n1 = static_cast<int>(d3.extent(3));
  if (rav_dist) {
  par_for("m1_rav_dxf", DevExeSpace(), 0, pmy_pack->nmb_thispack-1, 0, n3f-1, 0, n2-1,
          0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int n = static_cast<int>(nf(m,j));
    if (n > 1 && (((k - ks) % n) + n) % n == 0) {d3(m,k,j,i) *= static_cast<Real>(n);}
  });
  }
  if (global_variable::my_rank == 0) {
    std::cout << "RINGAVG <rad_m1>/polar_ring_avg: the implicit M1 on the coarse phi "
              << "groups (restricted Krylov system, coarse phi distance, lumped "
              << "preconditioner rows, group-averaged E, F and gas)" << std::endl;
  }
}

DvceFaceFld4D<Real> RadiationM1::M1Dxf() const {
  return rav_on ? rav_dxf : pmy_pack->pcoord->dxface;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RavAvg
//! \brief components v0..v0+nv-1 of a (m, n, k, j, i) array to their phi-group averages
//! over the active cells of the averaged rows (the loop of hydro::RingAvg::Average)

void RadiationM1::RavAvg(const DvceArray5D<Real> &a, const int v0, const int nv) {
  if (!rav_on || prav->nent == 0) return;
  auto &ic = pmy_pack->pmesh->mb_indcs;
  const int is = ic.is, ie = ic.ie, ks = ic.ks, ke = ic.ke;
  auto a_ = a;
  auto em_ = prav->em, ej_ = prav->ej, en_ = prav->en;
  par_for("m1_rav_avg", DevExeSpace(), 0, prav->nent-1, ks, ke, is, ie,
  KOKKOS_LAMBDA(const int e, const int k, const int i) {
    const int n = en_(e);
    if (((k - ks) % n) != 0) return;
    const int m = em_(e), j = ej_(e);
    const Real rn = 1.0/static_cast<Real>(n);
    for (int v = v0; v < v0 + nv; ++v) {
      Real s = 0.0;
      for (int kk = k; kk < k + n; ++kk) {s += a_(m,v,kk,j,i);}
      s *= rn;
      for (int kk = k; kk < k + n; ++kk) {a_(m,v,kk,j,i) = s;}
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RavLump
//! \brief (3) of the file header, polar_ring_avg_lump = true only: the phi couplings
//! inside a group onto the diagonal of the 7-point rows, before the pass's stencil /
//! preconditioner build.  Idempotent; the assembly rewrites TB and CKM/CKP every pass.

void RadiationM1::RavLump() {
  if (!rav_on || !rav_lump || prav->nent == 0) return;
  auto &ic = pmy_pack->pmesh->mb_indcs;
  const int is = ic.is, ie = ic.ie, ks = ic.ks, ke = ic.ke;
  auto iw_ = iw;
  auto em_ = prav->em, ej_ = prav->ej, en_ = prav->en;
  par_for("m1_rav_lump", DevExeSpace(), 0, prav->nent-1, ks, ke, is, ie,
  KOKKOS_LAMBDA(const int e, const int k, const int i) {
    const int n = en_(e), m = em_(e), j = ej_(e);
    const int q = (k - ks) % n;
    Real tb = iw_(m,M1_IW_TB,k,j,i);
    if (q > 0) {
      tb += iw_(m,M1_IW_CKM,k,j,i);
      iw_(m,M1_IW_CKM,k,j,i) = 0.0;
    }
    if (q < n - 1) {
      tb += iw_(m,M1_IW_CKP,k,j,i);
      iw_(m,M1_IW_CKP,k,j,i) = 0.0;
    }
    iw_(m,M1_IW_TB,k,j,i) = tb;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RavKrylovStart
//! \brief (1): the right-hand side and the first iterate of a linear solve to the
//! group-uniform space (called right after x0 = E^k)

void RadiationM1::RavKrylovStart() {
  RavAvg(iw, M1_IW_KB, 1);
  RavAvg(iw, M1_IW_KX, 1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RavFaces
//! \brief (4) of the file header: the stored face fluxes of the converged solve, before
//! the write-back

void RadiationM1::RavFaces() {
  if (!rav_on || prav->nent == 0) return;
  auto &ic = pmy_pack->pmesh->mb_indcs;
  const int is = ic.is, ie = ic.ie, ks = ic.ks, ke = ic.ke;
  auto em_ = prav->em, ej_ = prav->ej, en_ = prav->en;
  auto nf_ = prav->nfac;
  auto f1_ = f0x1;
  par_for("m1_rav_f1", DevExeSpace(), 0, prav->nent-1, ks, ke, is, ie+1,
  KOKKOS_LAMBDA(const int e, const int k, const int i) {
    const int n = en_(e);
    if (((k - ks) % n) != 0) return;
    const int m = em_(e), j = ej_(e);
    Real s = 0.0;
    for (int kk = k; kk < k + n; ++kk) {s += f1_(m,kk,j,i);}
    s /= static_cast<Real>(n);
    for (int kk = k; kk < k + n; ++kk) {f1_(m,kk,j,i) = s;}
  });
  if (!trans_on) return;
  auto f2_ = f0x2;
  // the lower theta face of each averaged row, over the smaller group of its two rows
  // (an x2 ghost row counts 1: the faces on a MeshBlock x2 edge are left alone)
  par_for("m1_rav_f2", DevExeSpace(), 0, prav->nent-1, ks, ke, is, ie,
  KOKKOS_LAMBDA(const int e, const int k, const int i) {
    const int m = em_(e), j = ej_(e);
    const int nb = static_cast<int>(nf_(m,j-1));
    const int n = (en_(e) < nb) ? en_(e) : nb;
    if (n < 2 || ((k - ks) % n) != 0) return;
    Real s = 0.0;
    for (int kk = k; kk < k + n; ++kk) {s += f2_(m,kk,j,i);}
    s /= static_cast<Real>(n);
    for (int kk = k; kk < k + n; ++kk) {f2_(m,kk,j,i) = s;}
  });
  if (!trans_x3) return;
  auto f3_ = f0x3;
  par_for("m1_rav_f3", DevExeSpace(), 0, prav->nent-1, ks, ke, is, ie,
  KOKKOS_LAMBDA(const int e, const int k, const int i) {
    const int n = en_(e);
    if (((k - ks) % n) != 0) return;
    const int m = em_(e), j = ej_(e);
    const Real fl = f3_(m,k,j,i), fr = f3_(m,k+n,j,i);
    const Real rn = 1.0/static_cast<Real>(n);
    for (int q = 1; q < n; ++q) {
      f3_(m,k+q,j,i) = fl + (static_cast<Real>(q)*rn)*(fr - fl);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RavState
//! \brief (5) of the file header: the radiation u0 and (gas = true) the hydro u0 to their
//! group averages after the solve

void RadiationM1::RavState(const bool gas) {
  if (!rav_on) return;
  RavAvg(u0, 0, M1_NVAR);
  if (!gas) return;
  hydro::Hydro *ph = pmy_pack->phydro;
  prav->Average(ph->u0, ph->nhydro + ph->nscalars, ph->phicc0, ph->use_etotgrav);
}

}  // namespace radm1
