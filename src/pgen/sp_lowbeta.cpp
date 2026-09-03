//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sp_lowbeta.cpp
//! \brief SPHERICAL-POLAR twin of the cubed sphere's low-beta stability test.
//!
//! WHY THIS EXISTS.  On the cubed sphere, a uniform medium AT REST carrying a uniform
//! Cartesian (hence force-free) field is linearly UNSTABLE below beta ~ 0.05, with the
//! threshold beta rising as the grid refines -- see cs_test iprob = 8 and
//! inputs/tests/cubed_sphere_mhd_strat.athinput.  A radial monopole at the same beta is
//! stable, so the feedback needs TANGENTIAL field.
//!
//! The question this file answers is whether that is a property of the GNOMONIC
//! machinery -- the non-orthogonal tangent basis, its face rotations and its geometric
//! source term -- or of AthenaK's curvilinear MHD source terms generally.  Spherical
//! polar discriminates: SrcTermsSphericalPolarMHD has the same centre-versus-face
//! structure (a cell-centred bcc in the source against reconstructed face states in the
//! flux) but an ORTHOGONAL basis and no seam.  Same exact solution, same beta, same
//! reconstruction and Riemann solver, same radial extent.
//!
//! THE EXACT SOLUTION IS THE INITIAL CONDITION, FOREVER, and v is exactly zero in it, so
//! the kinetic energy in the history file IS the error and the gate is whether it grows
//! exponentially.
//!
//! iprob = 1  uniform Cartesian field B = b0*bhat, built from A = (1/2) B x r so that
//!            div B is zero to round-off.  Tangential components everywhere.
//! iprob = 2  radial monopole B_r = b0*(r0/r)^2, the control: force-free too, but with no
//!            tangential component at all.
//!
//! Boundaries.  x3 (phi) spans the full 2 pi and is periodic.  x1 (r) and, if the domain
//! is a theta WEDGE, x2 both take the EXACT analytic state in their ghosts -- cells and
//! face fields alike -- so that nothing at a boundary can seed or drive what is being
//! measured.  Set ix1_bc/ox1_bc (and ix2_bc/ox2_bc) to `user` and <problem>/user_bcs.
//! A wedge is the default because it keeps the polar axis, which the cubed sphere has no
//! analogue of, out of the comparison.

#include <cmath>
#include <iostream>
#include <sstream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

namespace {
// shared between the generator and the boundary function
int  spl_iprob = 1;
Real spl_d0 = 1.0, spl_p0 = 0.1, spl_gm1 = 2.0/3.0;
Real spl_bx = 0.0, spl_by = 0.0, spl_bz = 1.0;   // the field VECTOR, not a direction
Real spl_b0 = 1.0, spl_r0 = 1.0;                 // monopole amplitude, reference radius

//! \fn SPLowBetaField
//! \brief the exact field at (r, theta, phi), in the SPHERICAL-POLAR orthonormal triple.
KOKKOS_INLINE_FUNCTION
void SPLowBetaField(const int iprob, const Real bx, const Real by, const Real bz,
                    const Real b0, const Real r0,
                    const Real r, const Real th, const Real ph,
                    Real &br, Real &bt, Real &bp) {
  if (iprob == 2) {
    br = b0*SQR(r0/r);
    bt = 0.0;
    bp = 0.0;
    return;
  }
  const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
  br =  bx*st*cp + by*st*sp + bz*ct;
  bt =  bx*ct*cp + by*ct*sp - bz*st;
  bp = -bx*sp    + by*cp;
}

//! \fn SPLowBetaBsq
//! \brief |B|^2 at radius r.  ANALYTIC, deliberately: taking it from the face-averaged
//! bcc instead would put the initial pressure out of balance by exactly the averaging
//! error, which is a perturbation of the same order as the thing being measured.
KOKKOS_INLINE_FUNCTION
Real SPLowBetaBsq(const int iprob, const Real bsq_uniform, const Real b0, const Real r0,
                  const Real r) {
  return (iprob == 2) ? SQR(b0*SQR(r0/r)) : bsq_uniform;
}
}  // namespace

void SPLowBetaBC(Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "sp_lowbeta requires <mhd>" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!pmy_mesh_->use_spherical_polar) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "sp_lowbeta requires <mesh>/use_spherical_polar"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  spl_iprob = pin->GetOrAddInteger("problem", "iprob", 1);
  spl_d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  spl_p0 = pin->GetOrAddReal("problem", "p0", 0.1);
  spl_b0 = pin->GetOrAddReal("problem", "b0c", 1.0);
  const Real bhx = pin->GetOrAddReal("problem", "bhx", 0.37);
  const Real bhy = pin->GetOrAddReal("problem", "bhy", 0.61);
  const Real bhz = pin->GetOrAddReal("problem", "bhz", 0.70);
  const Real bn = std::sqrt(bhx*bhx + bhy*bhy + bhz*bhz);
  spl_bx = spl_b0*bhx/bn;
  spl_by = spl_b0*bhy/bn;
  spl_bz = spl_b0*bhz/bn;
  spl_r0 = pmy_mesh_->mesh_size.x1min;
  spl_gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;
  user_bcs_func = SPLowBetaBC;
  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc = pmbp->pmhd->bcc0;
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x3v_ = pmbp->pcoord->x3v;
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x2f_ = pmbp->pcoord->xx2f;
  auto &x3f_ = pmbp->pcoord->xx3f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &area2 = pmbp->pcoord->area.x2f;
  auto &area3 = pmbp->pcoord->area.x3f;
  auto &dxe2 = pmbp->pcoord->dxedge.x2e;
  auto &dxe3 = pmbp->pcoord->dxedge.x3e;
  const int iprob = spl_iprob;
  const Real d0 = spl_d0, p0 = spl_p0, gm1 = spl_gm1;
  const Real bvx = spl_bx, bvy = spl_by, bvz = spl_bz;
  const Real b0m = spl_b0, r0m = spl_r0;
  const Real bsq_u = SQR(spl_b0);
  const int nmhd = pmbp->pmhd->nmhd;
  const int nscal = pmbp->pmhd->nscalars;

  // --- the face fields, from the vector potential ------------------------------------
  // A = (1/2) B x r for the uniform field, whose two nonzero components on the
  // orthonormal triple are  A.that = (r/2) B.phat  and  A.phat = -(r/2) B.that
  // -- note A.rhat vanishes identically, which is what lets the two TANGENTIAL faces
  // below take the circulation over their tangential edges alone.  A.thetahat happens not
  // to depend on theta at all, so the theta edge needs only (r, phi).
  par_for("splb_b", DevExeSpace(), 0,(pmbp->nmb_thispack-1), ks,ke+1, js,je+1, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (iprob == 2) {
      if (j <= je && k <= ke) b0.x1f(m,k,j,i) = b0m*SQR(r0m/x1f_(m,i));
      if (i <= ie && k <= ke) b0.x2f(m,k,j,i) = 0.0;
      if (i <= ie && j <= je) b0.x3f(m,k,j,i) = 0.0;
      return;
    }
    // A on the THETA edge (r face, theta centre, phi face) and on the PHI edge
    // (r face, theta face, phi centre), each projected on that edge's own unit tangent.
    auto Ath = [&](const int ii, const int kk) {
      const Real rf = x1f_(m,ii), ph = x3f_(m,kk);
      return 0.5*rf*(-bvx*sin(ph) + bvy*cos(ph));
    };
    auto Aph = [&](const int ii, const int jj, const int kk) {
      const Real rf = x1f_(m,ii), th = x2f_(m,jj), ph = x3v_(m,kk);
      return -0.5*rf*(bvx*cos(th)*cos(ph) + bvy*cos(th)*sin(ph) - bvz*sin(th));
    };
    // B.n = (1/area) * circulation of A around the face, exactly as mhd_ct.cpp does it
    if (j <= je && k <= ke) {
      b0.x1f(m,k,j,i) =
          (dxe3(m,k,j+1,i)*Aph(i,j+1,k) - dxe3(m,k,j,i)*Aph(i,j,k)
         - dxe2(m,k+1,j,i)*Ath(i,k+1)   + dxe2(m,k,j,i)*Ath(i,k))/area1(m,k,j,i);
    }
    if (i <= ie && k <= ke) {
      b0.x2f(m,k,j,i) =
          -(dxe3(m,k,j,i+1)*Aph(i+1,j,k) - dxe3(m,k,j,i)*Aph(i,j,k))/area2(m,k,j,i);
    }
    if (i <= ie && j <= je) {
      b0.x3f(m,k,j,i) =
          (dxe2(m,k,j,i+1)*Ath(i+1,k) - dxe2(m,k,j,i)*Ath(i,k))/area3(m,k,j,i);
    }
  });

  // --- the hydro state: uniform, at rest ----------------------------------------------
  par_for("splb_w", DevExeSpace(), 0,(pmbp->nmb_thispack-1), ks,ke, js,je, is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    w0(m,IDN,k,j,i) = d0;
    w0(m,IEN,k,j,i) = p0/gm1;
    w0(m,IVX,k,j,i) = 0.0;
    w0(m,IVY,k,j,i) = 0.0;
    w0(m,IVZ,k,j,i) = 0.0;
    bcc(m,IBX,k,j,i) = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
    bcc(m,IBY,k,j,i) = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
    bcc(m,IBZ,k,j,i) = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
    u0(m,IDN,k,j,i) = d0;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = p0/gm1
                    + 0.5*SPLowBetaBsq(iprob, bsq_u, b0m, r0m, x1v_(m,i));
    for (int n=nmhd; n<nmhd+nscal; ++n) {
      w0(m,n,k,j,i) = 1.0;
      u0(m,n,k,j,i) = d0;
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn SPLowBetaBC
//! \brief the EXACT state in every ghost of a boundary flagged `user`.
//!
//! Cells and face fields both.  A face is left alone when it is one the CT update evolves
//! -- that set is is <= i <= ie+1 (x1f), js <= j <= je+1 (x2f), ks <= k <= ke+1 (x3f),
//! with the other two indices in their ACTIVE range -- and every other face beyond a user
//! boundary is overwritten.  Faces in the ghosts of an INTERNAL block boundary are not
//! touched: those belong to the exchange, and filling them analytically would quietly
//! hand the run a perfect halo.

void SPLowBetaBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = indcs.nx2 + 2*ng;
  const int n3 = indcs.nx3 + 2*ng;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x3v_ = pmbp->pcoord->x3v;
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x2f_ = pmbp->pcoord->xx2f;
  auto &x3f_ = pmbp->pcoord->xx3f;
  const int iprob = spl_iprob;
  const Real d0 = spl_d0, p0 = spl_p0, gm1 = spl_gm1;
  const Real bvx = spl_bx, bvy = spl_by, bvz = spl_bz;
  const Real b0m = spl_b0, r0m = spl_r0;
  const Real bsq_u = SQR(spl_b0);
  const int nmhd = pmbp->pmhd->nmhd;
  const int nscal = pmbp->pmhd->nscalars;

  // cells
  par_for("splb_bc_c", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,(n3-1), 0,(n2-1),
          0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const bool g = (i < is && mb_bcs.d_view(m,BoundaryFace::inner_x1)
                              == BoundaryFlag::user)
                || (i > ie && mb_bcs.d_view(m,BoundaryFace::outer_x1)
                              == BoundaryFlag::user)
                || (j < js && mb_bcs.d_view(m,BoundaryFace::inner_x2)
                              == BoundaryFlag::user)
                || (j > je && mb_bcs.d_view(m,BoundaryFace::outer_x2)
                              == BoundaryFlag::user);
    if (!g) return;
    w0(m,IDN,k,j,i) = d0;
    w0(m,IEN,k,j,i) = p0/gm1;
    w0(m,IVX,k,j,i) = 0.0;
    w0(m,IVY,k,j,i) = 0.0;
    w0(m,IVZ,k,j,i) = 0.0;
    u0(m,IDN,k,j,i) = d0;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = p0/gm1
                    + 0.5*SPLowBetaBsq(iprob, bsq_u, b0m, r0m, x1v_(m,i));
    for (int n=nmhd; n<nmhd+nscal; ++n) {
      w0(m,n,k,j,i) = 1.0;
      u0(m,n,k,j,i) = d0;
    }
  });

  // face fields
  par_for("splb_bc_b", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,n3, 0,n2, 0,n1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const bool ix1u = (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user);
    const bool ox1u = (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user);
    const bool ix2u = (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user);
    const bool ox2u = (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user);
    Real br, bt, bp;
    // x1 face: normal is rhat, at the (j,k) cell centre
    if (j < n2 && k < n3 &&
        ((i <= is-1 && ix1u) || (i >= ie+2 && ox1u)
         || (j < js && ix2u) || (j > je && ox2u))) {
      SPLowBetaField(iprob, bvx, bvy, bvz, b0m, r0m,
                     x1f_(m,i), x2v_(m,j), x3v_(m,k), br, bt, bp);
      b0.x1f(m,k,j,i) = br;
    }
    // x2 face: normal is thetahat, at the theta FACE
    if (i < n1 && k < n3 &&
        ((j <= js-1 && ix2u) || (j >= je+2 && ox2u)
         || (i < is && ix1u) || (i > ie && ox1u))) {
      SPLowBetaField(iprob, bvx, bvy, bvz, b0m, r0m,
                     x1v_(m,i), x2f_(m,j), x3v_(m,k), br, bt, bp);
      b0.x2f(m,k,j,i) = bt;
    }
    // x3 face: normal is phihat, at the phi FACE.  x3 is periodic, so a phi face is only
    // ever rewritten because it sits beyond an x1 or x2 user boundary.
    if (i < n1 && j < n2 &&
        ((i < is && ix1u) || (i > ie && ox1u) || (j < js && ix2u) || (j > je && ox2u))) {
      SPLowBetaField(iprob, bvx, bvy, bvz, b0m, r0m,
                     x1v_(m,i), x2v_(m,j), x3f_(m,k), br, bt, bp);
      b0.x3f(m,k,j,i) = bp;
    }
  });
  return;
}
