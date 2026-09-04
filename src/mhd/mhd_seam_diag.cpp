//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_seam_diag.cpp
//! \brief CUBED SPHERE, MEASUREMENT ONLY: is a seam ghost cell's INTERNAL ENERGY
//! inconsistent with its own ghost field?
//!
//! THE HYPOTHESIS.  The cell-centred halo carries the CONSERVED variables -- SendU packs
//! u0, whose IEN is the TOTAL energy with 0.5|B_src|^2 already inside it -- while the
//! face field b0 crosses the seam through a SEPARATE exchange with its own co-location,
//! gnomonic transform, along-seam resample and monotone clamp.  ConsToPrim then forms
//!
//!     e_int = E_ghost - KE_ghost - 0.5|B_ghost|^2
//!
//! from two halves that were resampled by DIFFERENT operators: the energy as a scalar,
//! the field componentwise and then squared.  Resample-then-square is not
//! square-then-resample, and the whole residual lands in e_int, amplified by 1/beta.  If
//! that is large it is a floor trigger IN THE HALO -- boundary-associated, worst at low
//! beta, and invisible to every improvement made to the interior scheme.
//!
//! HOW IT IS MEASURED, without needing a second communication.  A seam is normal to x2 or
//! x3, i.e. to a PANEL-TANGENTIAL direction, along which the solution is smooth (the
//! atmosphere's steep gradient is radial).  So a quadratic extrapolation from the three
//! active cells nearest the boundary predicts what the first ghost cell should hold, to
//! O(h^3).  The diagnostic reports |q_ghost - q_extrap|/|q_extrap| for three quantities:
//!
//!     dens   -- a plainly communicated scalar.  THE CONTROL: it shares the resample
//!               but not the E - 0.5B^2 subtraction, so it carries the resample alone.
//!     eint   -- the suspect.
//!     magE   -- 0.5|bcc|^2, to show how big the magnetic energy is that gets subtracted.
//!
//! and splits every count two ways that make the answer falsifiable:
//!
//!     SEAM        vs  same-panel BLOCK FACE.  A block face uses the same halo machinery
//!                     with NO transform and NO resample -- a plain copy -- so it is the
//!                     null control.  If eint's anomaly is no worse at a seam than at a
//!                     block face, this mechanism is dead.
//!     low beta    vs  high beta.  The hypothesis is specifically that the residual is
//!                     amplified by 1/beta, so the seam/low-beta cell must be the worst
//!                     of the four by a wide margin.  If eint tracks dens everywhere,
//!                     the subtraction is innocent and only the resample is at work.
//!
//! It changes no answer: it reads w0 and bcc0 after ConsToPrim has filled them over the
//! ghost zones, and writes only its own counters.

#include <iostream>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd.hpp"
#include "eos/eos.hpp"

namespace mhd {

// rows: 0 = seam & low beta, 1 = seam & high beta, 2 = block face & low beta,
//       3 = block face & high beta.  cols: 3*q + {sum, count, max} for q in
//       {dens, eint, magE}.
static constexpr int NROW = 4;
static constexpr int NCOL = 9;


//----------------------------------------------------------------------------------------
//! helpers.  Plain KOKKOS_INLINE_FUNCTIONs rather than lambdas: an extended device lambda
//! captured inside another device lambda is the construct hipcc and nvcc restrict.

//! quadratic extrapolation one cell beyond the last active cell, L(-1) from points 0,1,2
KOKKOS_INLINE_FUNCTION
static Real Extrap3(const Real q0, const Real q1, const Real q2) {
  return 3.0*q0 - 3.0*q1 + q2;
}

KOKKOS_INLINE_FUNCTION
static Real MagE(const DvceArray5D<Real> &b, const int m, const int k, const int j,
                 const int i) {
  return 0.5*(SQR(b(m,IBX,k,j,i)) + SQR(b(m,IBY,k,j,i)) + SQR(b(m,IBZ,k,j,i)));
}

//! accumulate one ghost cell into the (seam/block) x (low/high beta) row.
//! beta is only a BIN BOUNDARY here, so the proxy 2*eint/B^2 is used rather than a real
//! EOS pressure: it differs from 2p/B^2 by the Gamma_1 - 1 factor, well inside one bin.
KOKKOS_INLINE_FUNCTION
static void Tally(const DvceArray2D<Real> &a, const int seam,
                  const Real dg, const Real eg, const Real bg,
                  const Real de, const Real ee, const Real be) {
  const Real beta = (bg > 0.0) ? (eg/bg) : 1.0e30;
  const int row = 2*(seam ? 0 : 1) + ((beta < 1.0) ? 0 : 1);
  const Real rq[3] = {(de != 0.0) ? fabs(dg - de)/fabs(de) : 0.0,
                      (ee != 0.0) ? fabs(eg - ee)/fabs(ee) : 0.0,
                      (be != 0.0) ? fabs(bg - be)/fabs(be) : 0.0};
  for (int q=0; q<3; ++q) {
    Kokkos::atomic_add(&a(row, 3*q+0), rq[q]);
    Kokkos::atomic_add(&a(row, 3*q+1), 1.0);
    Kokkos::atomic_max(&a(row, 3*q+2), rq[q]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::SeamHaloDiag
//! \brief accumulate and print the halo consistency measurement described above.

void MHD::SeamHaloDiag() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w0_ = w0;
  auto bcc_ = bcc0;
  auto &size_ = pmy_pack->pmb->mb_size;
  // A panel spans [-1,1] in both tangential coordinates, so a block whose x2min sits on
  // the mesh x2min has its x2-lo face ON a panel edge, i.e. a SEAM.  Everything else is
  // an ordinary same-panel block face.  Comparing the block's own extent against the mesh
  // extent needs no neighbour lookup and no panel table.
  const Real mx2min = pmy_pack->pmesh->mesh_size.x2min;
  const Real mx2max = pmy_pack->pmesh->mesh_size.x2max;
  const Real mx3min = pmy_pack->pmesh->mesh_size.x3min;
  const Real mx3max = pmy_pack->pmesh->mesh_size.x3max;
  const Real tol = 1.0e-10*(mx2max - mx2min);

  DvceArray2D<Real> acc("seam_diag", NROW, NCOL);
  Kokkos::deep_copy(acc, 0.0);

  // ---- the two x2 boundaries
  par_for("seamdiag_x2", DevExeSpace(), 0, nmb1, ks, ke, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    for (int side=0; side<2; ++side) {
      const bool lo = (side == 0);
      const int jg = lo ? (js-1) : (je+1);
      const int j0 = lo ? js : je;
      const int j1 = lo ? (js+1) : (je-1);
      const int j2 = lo ? (js+2) : (je-2);
      const Real edge = lo ? size_.d_view(m).x2min : size_.d_view(m).x2max;
      const Real panel = lo ? mx2min : mx2max;
      const int seam = (fabs(edge - panel) < tol) ? 1 : 0;
      Tally(acc, seam,
            w0_(m,IDN,k,jg,i), w0_(m,IEN,k,jg,i), MagE(bcc_,m,k,jg,i),
            Extrap3(w0_(m,IDN,k,j0,i), w0_(m,IDN,k,j1,i), w0_(m,IDN,k,j2,i)),
            Extrap3(w0_(m,IEN,k,j0,i), w0_(m,IEN,k,j1,i), w0_(m,IEN,k,j2,i)),
            Extrap3(MagE(bcc_,m,k,j0,i), MagE(bcc_,m,k,j1,i), MagE(bcc_,m,k,j2,i)));
    }
  });

  // ---- the two x3 boundaries
  par_for("seamdiag_x3", DevExeSpace(), 0, nmb1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int j, const int i) {
    for (int side=0; side<2; ++side) {
      const bool lo = (side == 0);
      const int kg = lo ? (ks-1) : (ke+1);
      const int k0 = lo ? ks : ke;
      const int k1 = lo ? (ks+1) : (ke-1);
      const int k2 = lo ? (ks+2) : (ke-2);
      const Real edge = lo ? size_.d_view(m).x3min : size_.d_view(m).x3max;
      const Real panel = lo ? mx3min : mx3max;
      const int seam = (fabs(edge - panel) < tol) ? 1 : 0;
      Tally(acc, seam,
            w0_(m,IDN,kg,j,i), w0_(m,IEN,kg,j,i), MagE(bcc_,m,kg,j,i),
            Extrap3(w0_(m,IDN,k0,j,i), w0_(m,IDN,k1,j,i), w0_(m,IDN,k2,j,i)),
            Extrap3(w0_(m,IEN,k0,j,i), w0_(m,IEN,k1,j,i), w0_(m,IEN,k2,j,i)),
            Extrap3(MagE(bcc_,m,k0,j,i), MagE(bcc_,m,k1,j,i), MagE(bcc_,m,k2,j,i)));
    }
  });

  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), acc);
  std::vector<double> v(NROW*NCOL);
  for (int r=0; r<NROW; ++r) { for (int c=0; c<NCOL; ++c) { v[r*NCOL+c] = h(r,c); } }
#if MPI_PARALLEL_ENABLED
  // sums and counts add; the maxima have to be reduced as maxima, so they go separately
  std::vector<double> mx(NROW*3);
  for (int r=0; r<NROW; ++r) {
    for (int q=0; q<3; ++q) { mx[r*3+q] = v[r*NCOL+3*q+2]; v[r*NCOL+3*q+2] = 0.0; }
  }
  MPI_Allreduce(MPI_IN_PLACE, v.data(), NROW*NCOL, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, mx.data(), NROW*3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  for (int r=0; r<NROW; ++r) {
    for (int q=0; q<3; ++q) { v[r*NCOL+3*q+2] = mx[r*3+q]; }
  }
#endif
  if (global_variable::my_rank != 0) return;

  const char *rname[NROW] = {"SEAM  beta<1", "SEAM  beta>1",
                             "block beta<1", "block beta>1"};
  std::cout << "## cs_seam_diag: |ghost - quadratic extrapolation| / |extrapolation|"
            << std::endl;
  std::cout << "##   region          ncell        dens mean/max         "
            << "eint mean/max         magE mean/max" << std::endl;
  for (int r=0; r<NROW; ++r) {
    const double n = v[r*NCOL+1];
    std::cout << "##  " << rname[r] << std::setw(11)
              << static_cast<std::int64_t>(n);
    for (int q=0; q<3; ++q) {
      const double mean = (n > 0.0) ? v[r*NCOL+3*q+0]/n : 0.0;
      std::cout << "   " << std::scientific << std::setprecision(2) << mean
                << " " << v[r*NCOL+3*q+2];
    }
    std::cout << std::endl;
  }
  std::cout.unsetf(std::ios_base::floatfield);
}

}  // namespace mhd
