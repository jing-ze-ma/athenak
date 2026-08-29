#ifndef MESH_PROLONGATION_HPP_
#define MESH_PROLONGATION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file prolongation.hpp
//! \brief prolongation operators for cell-centered and face-centered variables,
//! implemented as inline functions so they can be used both in Bval and AMR functions.

#include "z4c/z4c.hpp"

//----------------------------------------------------------------------------------------
//! \fn ProlongCC()
//! \brief 2nd-order (piecewise-linear) prolongation operator for cell-centered variables

KOKKOS_INLINE_FUNCTION
void ProlongCC(const int m, const int v, const int k, const int j, const int i,
               const int fk, const int fj, const int fi,
               const bool multi_d, const bool three_d,
               const DvceArray5D<Real> &ca, const DvceArray5D<Real> &a) {
  // calculate x1-gradient using the min-mod limiter
  Real dl = ca(m,v,k,j,i  ) - ca(m,v,k,j,i-1);
  Real dr = ca(m,v,k,j,i+1) - ca(m,v,k,j,i  );
  Real dvar1 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));

  // calculate x2-gradient using the min-mod limiter
  Real dvar2 = 0.0;
  if (multi_d) {
    dl = ca(m,v,k,j  ,i) - ca(m,v,k,j-1,i);
    dr = ca(m,v,k,j+1,i) - ca(m,v,k,j  ,i);
    dvar2 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  }

  // calculate x1-gradient using the min-mod limiter
  Real dvar3 = 0.0;
  if (three_d) {
    dl = ca(m,v,k  ,j,i) - ca(m,v,k-1,j,i);
    dr = ca(m,v,k+1,j,i) - ca(m,v,k  ,j,i);
    dvar3 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  }

  // interpolate to the finer grid
  a(m,v,fk,fj,fi  ) = ca(m,v,k,j,i) - dvar1 - dvar2 - dvar3;
  a(m,v,fk,fj,fi+1) = ca(m,v,k,j,i) + dvar1 - dvar2 - dvar3;
  if (multi_d) {
    a(m,v,fk,fj+1,fi  ) = ca(m,v,k,j,i) - dvar1 + dvar2 - dvar3;
    a(m,v,fk,fj+1,fi+1) = ca(m,v,k,j,i) + dvar1 + dvar2 - dvar3;
  }
  if (three_d) {
    a(m,v,fk+1,fj  ,fi  ) = ca(m,v,k,j,i) - dvar1 - dvar2 + dvar3;
    a(m,v,fk+1,fj  ,fi+1) = ca(m,v,k,j,i) + dvar1 - dvar2 + dvar3;
    a(m,v,fk+1,fj+1,fi  ) = ca(m,v,k,j,i) - dvar1 + dvar2 + dvar3;
    a(m,v,fk+1,fj+1,fi+1) = ca(m,v,k,j,i) + dvar1 + dvar2 + dvar3;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn ProlongFCSharedX1Face()
//! \brief 2nd-order (piecewise-linear) prolongation operator for face-centered variables
//! on shared X1-faces between fine and coarse cells

KOKKOS_INLINE_FUNCTION
void ProlongFCSharedX1Face(const int m, const int k, const int j, const int i,
                   const int fk, const int fj, const int fi,
                   const bool multi_d, const bool three_d,
                   const DvceArray4D<Real> &cbx1f, const DvceArray4D<Real> &bx1f) {
  // Prolongate b.x1f (v=0) by interpolating in x2/x3
  Real dvar2 = 0.0;
  if (multi_d) {
    Real dl = cbx1f(m,k,j  ,i) - cbx1f(m,k,j-1,i);
    Real dr = cbx1f(m,k,j+1,i) - cbx1f(m,k,j  ,i);
    dvar2 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  }

  Real dvar3 = 0.0;
  if (three_d) {
    Real dl = cbx1f(m,k  ,j,i) - cbx1f(m,k-1,j,i);
    Real dr = cbx1f(m,k+1,j,i) - cbx1f(m,k  ,j,i);
    dvar3 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  }

  bx1f(m,fk,fj,fi) = cbx1f(m,k,j,i) - dvar2 - dvar3;
  if (multi_d) {
    bx1f(m,fk,fj+1,fi) = cbx1f(m,k,j,i) + dvar2 - dvar3;
  }
  if (three_d) {
    bx1f(m,fk+1,fj  ,fi) = cbx1f(m,k,j,i) - dvar2 + dvar3;
    bx1f(m,fk+1,fj+1,fi) = cbx1f(m,k,j,i) + dvar2 + dvar3;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn ProlongFCSharedX2Face()
//! \brief 2nd-order (piecewise-linear) prolongation operator for face-centered variables
//! on shared X2-faces between fine and coarse cells

KOKKOS_INLINE_FUNCTION
void ProlongFCSharedX2Face(const int m, const int k, const int j, const int i,
                   const int fk, const int fj, const int fi,
                   const bool three_d,
                   const DvceArray4D<Real> &cbx2f, const DvceArray4D<Real> &bx2f) {
  // Prolongate b.x2f (v=1) by interpolating in x1/x3
  Real dl = cbx2f(m,k,j,i  ) - cbx2f(m,k,j,i-1);
  Real dr = cbx2f(m,k,j,i+1) - cbx2f(m,k,j,i  );
  Real dvar1 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));

  Real dvar3 = 0.0;
  if (three_d) {
    dl = cbx2f(m,k  ,j,i) - cbx2f(m,k-1,j,i);
    dr = cbx2f(m,k+1,j,i) - cbx2f(m,k  ,j,i);
    dvar3 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  }

  bx2f(m,fk  ,fj,fi  ) = cbx2f(m,k,j,i) - dvar1 - dvar3;
  bx2f(m,fk  ,fj,fi+1) = cbx2f(m,k,j,i) + dvar1 - dvar3;
  if (three_d) {
    bx2f(m,fk+1,fj,fi  ) = cbx2f(m,k,j,i) - dvar1 + dvar3;
    bx2f(m,fk+1,fj,fi+1) = cbx2f(m,k,j,i) + dvar1 + dvar3;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn ProlongFCSharedX3Face()
//! \brief 2nd-order (piecewise-linear) prolongation operator for face-centered variables
//! on shared X3-faces between fine and coarse cells

KOKKOS_INLINE_FUNCTION
void ProlongFCSharedX3Face(const int m, const int k, const int j, const int i,
                   const int fk, const int fj, const int fi,
                   const bool multi_d,
                   const DvceArray4D<Real> &cbx3f, const DvceArray4D<Real> &bx3f) {
  // Prolongate b.x3f (v=2) by interpolating in x1/x2
  Real dl = cbx3f(m,k,j,i  ) - cbx3f(m,k,j,i-1);
  Real dr = cbx3f(m,k,j,i+1) - cbx3f(m,k,j,i  );
  Real dvar1 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));

  Real dvar2 = 0.0;
  if (multi_d) {
    dl = cbx3f(m,k,j  ,i) - cbx3f(m,k,j-1,i);
    dr = cbx3f(m,k,j+1,i) - cbx3f(m,k,j  ,i);
    dvar2 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  }

  bx3f(m,fk,fj  ,fi  ) = cbx3f(m,k,j,i) - dvar1 - dvar2;
  bx3f(m,fk,fj  ,fi+1) = cbx3f(m,k,j,i) + dvar1 - dvar2;
  if (multi_d) {
    bx3f(m,fk,fj+1,fi  ) = cbx3f(m,k,j,i) - dvar1 + dvar2;
    bx3f(m,fk,fj+1,fi+1) = cbx3f(m,k,j,i) + dvar1 + dvar2;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn ProlongFCSharedX1FaceCurvi(), X2, X3
//! \brief The shared-face FC prolongation on a CURVILINEAR grid.
//!
//! The Cartesian operators above write B_c +/- dvar2 +/- dvar3, whose PLAIN mean is B_c.
//! That is the right thing only when the four fine faces have equal areas: what the
//! divergence actually sums is the FLUX A*B, so the property to preserve is
//! sum(A_f B_f) = A_c B_c with A_c the total area of the coarse face. Taking A_c to be
//! sum(A_f) -- which is exact by construction, unlike evaluating the gnomonic area
//! formula on the coarse mesh -- that reduces to sum(A_f d_f) = 0 on the deviations.
//! The Cartesian operator only gives sum(d_f) = 0, so subtract the AREA-WEIGHTED mean of
//! the deviations. On a uniform grid the correction is identically zero and these agree
//! with the operators above term for term.

KOKKOS_INLINE_FUNCTION
void ProlongFCSharedX1FaceCurvi(const int m, const int k, const int j, const int i,
                   const int fk, const int fj, const int fi,
                   const bool multi_d, const bool three_d,
                   const DvceArray4D<Real> &cbx1f, const DvceArray4D<Real> &bx1f,
                   const DvceArray4D<Real> &ax1f) {
  Real dvar2 = 0.0;
  if (multi_d) {
    Real dl = cbx1f(m,k,j  ,i) - cbx1f(m,k,j-1,i);
    Real dr = cbx1f(m,k,j+1,i) - cbx1f(m,k,j  ,i);
    dvar2 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  }
  Real dvar3 = 0.0;
  if (three_d) {
    Real dl = cbx1f(m,k  ,j,i) - cbx1f(m,k-1,j,i);
    Real dr = cbx1f(m,k+1,j,i) - cbx1f(m,k  ,j,i);
    dvar3 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  }
  const Real b0 = cbx1f(m,k,j,i);
  if (three_d) {
    const Real a00 = ax1f(m,fk  ,fj  ,fi), a10 = ax1f(m,fk  ,fj+1,fi);
    const Real a01 = ax1f(m,fk+1,fj  ,fi), a11 = ax1f(m,fk+1,fj+1,fi);
    const Real d00 = -dvar2 - dvar3, d10 =  dvar2 - dvar3;
    const Real d01 = -dvar2 + dvar3, d11 =  dvar2 + dvar3;
    const Real w = (a00*d00 + a10*d10 + a01*d01 + a11*d11)/(a00 + a10 + a01 + a11);
    bx1f(m,fk  ,fj  ,fi) = b0 + d00 - w;
    bx1f(m,fk  ,fj+1,fi) = b0 + d10 - w;
    bx1f(m,fk+1,fj  ,fi) = b0 + d01 - w;
    bx1f(m,fk+1,fj+1,fi) = b0 + d11 - w;
  } else if (multi_d) {
    const Real a0 = ax1f(m,fk,fj,fi), a1 = ax1f(m,fk,fj+1,fi);
    const Real w = (a1 - a0)*dvar2/(a0 + a1);
    bx1f(m,fk,fj  ,fi) = b0 - dvar2 - w;
    bx1f(m,fk,fj+1,fi) = b0 + dvar2 - w;
  } else {
    bx1f(m,fk,fj,fi) = b0;
  }
  return;
}

KOKKOS_INLINE_FUNCTION
void ProlongFCSharedX2FaceCurvi(const int m, const int k, const int j, const int i,
                   const int fk, const int fj, const int fi, const bool three_d,
                   const DvceArray4D<Real> &cbx2f, const DvceArray4D<Real> &bx2f,
                   const DvceArray4D<Real> &ax2f) {
  Real dl = cbx2f(m,k,j,i  ) - cbx2f(m,k,j,i-1);
  Real dr = cbx2f(m,k,j,i+1) - cbx2f(m,k,j,i  );
  Real dvar1 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  Real dvar3 = 0.0;
  if (three_d) {
    Real el = cbx2f(m,k  ,j,i) - cbx2f(m,k-1,j,i);
    Real er = cbx2f(m,k+1,j,i) - cbx2f(m,k  ,j,i);
    dvar3 = 0.125*(SIGN(el) + SIGN(er))*fmin(fabs(el), fabs(er));
  }
  const Real b0 = cbx2f(m,k,j,i);
  if (three_d) {
    const Real a00 = ax2f(m,fk  ,fj,fi  ), a10 = ax2f(m,fk  ,fj,fi+1);
    const Real a01 = ax2f(m,fk+1,fj,fi  ), a11 = ax2f(m,fk+1,fj,fi+1);
    const Real d00 = -dvar1 - dvar3, d10 =  dvar1 - dvar3;
    const Real d01 = -dvar1 + dvar3, d11 =  dvar1 + dvar3;
    const Real w = (a00*d00 + a10*d10 + a01*d01 + a11*d11)/(a00 + a10 + a01 + a11);
    bx2f(m,fk  ,fj,fi  ) = b0 + d00 - w;
    bx2f(m,fk  ,fj,fi+1) = b0 + d10 - w;
    bx2f(m,fk+1,fj,fi  ) = b0 + d01 - w;
    bx2f(m,fk+1,fj,fi+1) = b0 + d11 - w;
  } else {
    const Real a0 = ax2f(m,fk,fj,fi), a1 = ax2f(m,fk,fj,fi+1);
    const Real w = (a1 - a0)*dvar1/(a0 + a1);
    bx2f(m,fk,fj,fi  ) = b0 - dvar1 - w;
    bx2f(m,fk,fj,fi+1) = b0 + dvar1 - w;
  }
  return;
}

KOKKOS_INLINE_FUNCTION
void ProlongFCSharedX3FaceCurvi(const int m, const int k, const int j, const int i,
                   const int fk, const int fj, const int fi, const bool multi_d,
                   const DvceArray4D<Real> &cbx3f, const DvceArray4D<Real> &bx3f,
                   const DvceArray4D<Real> &ax3f) {
  Real dl = cbx3f(m,k,j,i  ) - cbx3f(m,k,j,i-1);
  Real dr = cbx3f(m,k,j,i+1) - cbx3f(m,k,j,i  );
  Real dvar1 = 0.125*(SIGN(dl) + SIGN(dr))*fmin(fabs(dl), fabs(dr));
  Real dvar2 = 0.0;
  if (multi_d) {
    Real el = cbx3f(m,k,j  ,i) - cbx3f(m,k,j-1,i);
    Real er = cbx3f(m,k,j+1,i) - cbx3f(m,k,j  ,i);
    dvar2 = 0.125*(SIGN(el) + SIGN(er))*fmin(fabs(el), fabs(er));
  }
  const Real b0 = cbx3f(m,k,j,i);
  if (multi_d) {
    const Real a00 = ax3f(m,fk,fj  ,fi  ), a10 = ax3f(m,fk,fj  ,fi+1);
    const Real a01 = ax3f(m,fk,fj+1,fi  ), a11 = ax3f(m,fk,fj+1,fi+1);
    const Real d00 = -dvar1 - dvar2, d10 =  dvar1 - dvar2;
    const Real d01 = -dvar1 + dvar2, d11 =  dvar1 + dvar2;
    const Real w = (a00*d00 + a10*d10 + a01*d01 + a11*d11)/(a00 + a10 + a01 + a11);
    bx3f(m,fk,fj  ,fi  ) = b0 + d00 - w;
    bx3f(m,fk,fj  ,fi+1) = b0 + d10 - w;
    bx3f(m,fk,fj+1,fi  ) = b0 + d01 - w;
    bx3f(m,fk,fj+1,fi+1) = b0 + d11 - w;
  } else {
    const Real a0 = ax3f(m,fk,fj,fi), a1 = ax3f(m,fk,fj,fi+1);
    const Real w = (a1 - a0)*dvar1/(a0 + a1);
    bx3f(m,fk,fj,fi  ) = b0 - dvar1 - w;
    bx3f(m,fk,fj,fi+1) = b0 + dvar1 - w;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn ProlongInternalFC()
//! \brief 2nd-order prolongation operator for face-centered variables on internal edges
//! of new fine cells within one coarse cell using divergence-preserving interpolation
//! scheme of Toth & Roe, JCP 180, 736 (2002).

//----------------------------------------------------------------------------------------
//! \fn ProlongInternalFC()
//! \brief 2nd-order prolongation operator for face-centered variables on internal edges
//! of new fine cells within one coarse cell using divergence-preserving interpolation
//! scheme of Toth & Roe, JCP 180, 736 (2002).

KOKKOS_INLINE_FUNCTION
void ProlongFCInternal(const int m, const int fk, const int fj, const int fi,
                       const bool three_d, const DvceFaceFld4D<Real> &b) {
  // Prolongate internal fields in 3D
  if (three_d) {
    Real Uxx  = 0.0, Vyy  = 0.0, Wzz  = 0.0;
    Real Uxyz = 0.0, Vxyz = 0.0, Wxyz = 0.0;
    for (int jj=0; jj<2; jj++) {
      int jsgn = 2*jj - 1;
      int fjj  = fj + jj, fjp = fj + 2*jj;
      for (int ii=0; ii<2; ii++) {
        int isgn = 2*ii - 1;
        int fii = fi + ii, fip = fi + 2*ii;
        Uxx += isgn*(jsgn*(b.x2f(m,fk  ,fjp,fii) + b.x2f(m,fk+1,fjp,fii)) +
                          (b.x3f(m,fk+2,fjj,fii) - b.x3f(m,fk  ,fjj,fii)));

        Vyy += jsgn*(     (b.x3f(m,fk+2,fjj,fii) - b.x3f(m,fk  ,fjj,fii)) +
                     isgn*(b.x1f(m,fk  ,fjj,fip) + b.x1f(m,fk+1,fjj,fip)));

        Wzz +=       isgn*(b.x1f(m,fk+1,fjj,fip) - b.x1f(m,fk  ,fjj,fip)) +
                     jsgn*(b.x2f(m,fk+1,fjp,fii) - b.x2f(m,fk  ,fjp,fii));

        Uxyz += isgn*jsgn*(b.x1f(m,fk+1,fjj,fip) - b.x1f(m,fk  ,fjj,fip));
        Vxyz += isgn*jsgn*(b.x2f(m,fk+1,fjp,fii) - b.x2f(m,fk  ,fjp,fii));
        Wxyz += isgn*jsgn*(b.x3f(m,fk+2,fjj,fii) - b.x3f(m,fk  ,fjj,fii));
      }
    }
    Uxx *= 0.125;  Vyy *= 0.125;  Wzz *= 0.125;
    Uxyz *= 0.0625; Vxyz *= 0.0625; Wxyz *= 0.0625;

    b.x1f(m,fk  ,fj  ,fi+1) = 0.5*(b.x1f(m,fk  ,fj  ,fi  ) + b.x1f(m,fk  ,fj  ,fi+2))
                            + Uxx - Vxyz - Wxyz;
    b.x1f(m,fk  ,fj+1,fi+1) = 0.5*(b.x1f(m,fk  ,fj+1,fi  ) + b.x1f(m,fk  ,fj+1,fi+2))
                            + Uxx - Vxyz + Wxyz;
    b.x1f(m,fk+1,fj  ,fi+1) = 0.5*(b.x1f(m,fk+1,fj  ,fi  ) + b.x1f(m,fk+1,fj  ,fi+2))
                            + Uxx + Vxyz - Wxyz;
    b.x1f(m,fk+1,fj+1,fi+1) = 0.5*(b.x1f(m,fk+1,fj+1,fi  ) + b.x1f(m,fk+1,fj+1,fi+2))
                            + Uxx + Vxyz + Wxyz;
    b.x2f(m,fk  ,fj+1,fi  ) = 0.5*(b.x2f(m,fk  ,fj  ,fi  ) + b.x2f(m,fk  ,fj+2,fi  ))
                            + Vyy - Uxyz - Wxyz;
    b.x2f(m,fk  ,fj+1,fi+1) = 0.5*(b.x2f(m,fk  ,fj  ,fi+1) + b.x2f(m,fk  ,fj+2,fi+1))
                            + Vyy - Uxyz + Wxyz;
    b.x2f(m,fk+1,fj+1,fi  ) = 0.5*(b.x2f(m,fk+1,fj  ,fi  ) + b.x2f(m,fk+1,fj+2,fi  ))
                            + Vyy + Uxyz - Wxyz;
    b.x2f(m,fk+1,fj+1,fi+1) = 0.5*(b.x2f(m,fk+1,fj  ,fi+1) + b.x2f(m,fk+1,fj+2,fi+1))
                            + Vyy + Uxyz + Wxyz;
    b.x3f(m,fk+1,fj  ,fi  ) = 0.5*(b.x3f(m,fk+2,fj  ,fi  ) + b.x3f(m,fk  ,fj  ,fi  ))
                            + Wzz - Uxyz - Vxyz;
    b.x3f(m,fk+1,fj  ,fi+1) = 0.5*(b.x3f(m,fk+2,fj  ,fi+1) + b.x3f(m,fk  ,fj  ,fi+1))
                            + Wzz - Uxyz + Vxyz;
    b.x3f(m,fk+1,fj+1,fi  ) = 0.5*(b.x3f(m,fk+2,fj+1,fi  ) + b.x3f(m,fk  ,fj+1,fi  ))
                            + Wzz + Uxyz - Vxyz;
    b.x3f(m,fk+1,fj+1,fi+1) = 0.5*(b.x3f(m,fk+2,fj+1,fi+1) + b.x3f(m,fk  ,fj+1,fi+1))
                            + Wzz + Uxyz + Vxyz;

  // Prolongate internal fields in 2D
  } else {
    Real tmp1 = 0.25*(b.x2f(m,fk,fj+2,fi+1) - b.x2f(m,fk,fj,  fi+1)
                    - b.x2f(m,fk,fj+2,fi  ) + b.x2f(m,fk,fj,  fi  ));
    Real tmp2 = 0.25*(b.x1f(m,fk,fj,  fi  ) - b.x1f(m,fk,fj,  fi+2)
                    - b.x1f(m,fk,fj+1,fi  ) + b.x1f(m,fk,fj+1,fi+2));
    b.x1f(m,fk,fj  ,fi+1) = 0.5*(b.x1f(m,fk,fj,  fi  ) + b.x1f(m,fk,fj,  fi+2)) + tmp1;
    b.x1f(m,fk,fj+1,fi+1) = 0.5*(b.x1f(m,fk,fj+1,fi  ) + b.x1f(m,fk,fj+1,fi+2)) + tmp1;
    b.x2f(m,fk,fj+1,fi  ) = 0.5*(b.x2f(m,fk,fj,  fi  ) + b.x2f(m,fk,fj+2,fi  )) + tmp2;
    b.x2f(m,fk,fj+1,fi+1) = 0.5*(b.x2f(m,fk,fj,  fi+1) + b.x2f(m,fk,fj+2,fi+1)) + tmp2;
  }
  return;
}
//----------------------------------------------------------------------------------------
//! \fn ProlongFCInternalCurvi()
//! \brief The Toth & Roe internal-face operator on a CURVILINEAR grid.
//!
//! Term for term the same algebra as ProlongFCInternal, but carried out on the face FLUX
//! A*B instead of on B. That is the whole change, and it is what makes the operator do
//! its job here: what Toth & Roe enforce is the discrete solenoidality of the
//! interpolated field, which on a uniform Cartesian mesh reads sum(+/-B) = 0 because
//! every face has the same area, and on a gnomonic grid reads sum(+/-A*B) = 0. Feeding
//! it B directly therefore enforces the wrong constraint by an O(1) factor -- measured,
//! that lost five orders of magnitude of magnetic energy in a single cycle and reached
//! NaN on the next.
//!
//! Each internal face is then converted back with B = Phi / A at its own location.

KOKKOS_INLINE_FUNCTION
void ProlongFCInternalCurvi(const int m, const int fk, const int fj, const int fi,
                       const bool three_d, const DvceFaceFld4D<Real> &b,
                       const DvceArray4D<Real> &a1,
                       const DvceArray4D<Real> &a2,
                       const DvceArray4D<Real> &a3) {
  // Flux accessors: everything below is the Cartesian algebra with A*B in place of B.
  auto P1 = [&](const int kk, const int jj, const int ii) {
    return b.x1f(m,kk,jj,ii)*a1(m,kk,jj,ii);
  };
  auto P2 = [&](const int kk, const int jj, const int ii) {
    return b.x2f(m,kk,jj,ii)*a2(m,kk,jj,ii);
  };
  auto P3 = [&](const int kk, const int jj, const int ii) {
    return b.x3f(m,kk,jj,ii)*a3(m,kk,jj,ii);
  };
  // Prolongate internal fields in 3D
  if (three_d) {
    Real Uxx  = 0.0, Vyy  = 0.0, Wzz  = 0.0;
    Real Uxyz = 0.0, Vxyz = 0.0, Wxyz = 0.0;
    for (int jj=0; jj<2; jj++) {
      int jsgn = 2*jj - 1;
      int fjj  = fj + jj, fjp = fj + 2*jj;
      for (int ii=0; ii<2; ii++) {
        int isgn = 2*ii - 1;
        int fii = fi + ii, fip = fi + 2*ii;
        Uxx += isgn*(jsgn*(P2(fk  ,fjp,fii) + P2(fk+1,fjp,fii)) +
                          (P3(fk+2,fjj,fii) - P3(fk  ,fjj,fii)));

        Vyy += jsgn*(     (P3(fk+2,fjj,fii) - P3(fk  ,fjj,fii)) +
                     isgn*(P1(fk  ,fjj,fip) + P1(fk+1,fjj,fip)));

        Wzz +=       isgn*(P1(fk+1,fjj,fip) - P1(fk  ,fjj,fip)) +
                     jsgn*(P2(fk+1,fjp,fii) - P2(fk  ,fjp,fii));

        Uxyz += isgn*jsgn*(P1(fk+1,fjj,fip) - P1(fk  ,fjj,fip));
        Vxyz += isgn*jsgn*(P2(fk+1,fjp,fii) - P2(fk  ,fjp,fii));
        Wxyz += isgn*jsgn*(P3(fk+2,fjj,fii) - P3(fk  ,fjj,fii));
      }
    }
    Uxx *= 0.125;  Vyy *= 0.125;  Wzz *= 0.125;
    Uxyz *= 0.0625; Vxyz *= 0.0625; Wxyz *= 0.0625;

    b.x1f(m,fk  ,fj  ,fi+1) = (0.5*(P1(fk  ,fj  ,fi  ) + P1(fk  ,fj  ,fi+2))
                            + Uxx - Vxyz - Wxyz)
        / a1(m,fk  ,fj  ,fi+1);
    b.x1f(m,fk  ,fj+1,fi+1) = (0.5*(P1(fk  ,fj+1,fi  ) + P1(fk  ,fj+1,fi+2))
                            + Uxx - Vxyz + Wxyz)
        / a1(m,fk  ,fj+1,fi+1);
    b.x1f(m,fk+1,fj  ,fi+1) = (0.5*(P1(fk+1,fj  ,fi  ) + P1(fk+1,fj  ,fi+2))
                            + Uxx + Vxyz - Wxyz)
        / a1(m,fk+1,fj  ,fi+1);
    b.x1f(m,fk+1,fj+1,fi+1) = (0.5*(P1(fk+1,fj+1,fi  ) + P1(fk+1,fj+1,fi+2))
                            + Uxx + Vxyz + Wxyz)
        / a1(m,fk+1,fj+1,fi+1);
    b.x2f(m,fk  ,fj+1,fi  ) = (0.5*(P2(fk  ,fj  ,fi  ) + P2(fk  ,fj+2,fi  ))
                            + Vyy - Uxyz - Wxyz)
        / a2(m,fk  ,fj+1,fi  );
    b.x2f(m,fk  ,fj+1,fi+1) = (0.5*(P2(fk  ,fj  ,fi+1) + P2(fk  ,fj+2,fi+1))
                            + Vyy - Uxyz + Wxyz)
        / a2(m,fk  ,fj+1,fi+1);
    b.x2f(m,fk+1,fj+1,fi  ) = (0.5*(P2(fk+1,fj  ,fi  ) + P2(fk+1,fj+2,fi  ))
                            + Vyy + Uxyz - Wxyz)
        / a2(m,fk+1,fj+1,fi  );
    b.x2f(m,fk+1,fj+1,fi+1) = (0.5*(P2(fk+1,fj  ,fi+1) + P2(fk+1,fj+2,fi+1))
                            + Vyy + Uxyz + Wxyz)
        / a2(m,fk+1,fj+1,fi+1);
    b.x3f(m,fk+1,fj  ,fi  ) = (0.5*(P3(fk+2,fj  ,fi  ) + P3(fk  ,fj  ,fi  ))
                            + Wzz - Uxyz - Vxyz)
        / a3(m,fk+1,fj  ,fi  );
    b.x3f(m,fk+1,fj  ,fi+1) = (0.5*(P3(fk+2,fj  ,fi+1) + P3(fk  ,fj  ,fi+1))
                            + Wzz - Uxyz + Vxyz)
        / a3(m,fk+1,fj  ,fi+1);
    b.x3f(m,fk+1,fj+1,fi  ) = (0.5*(P3(fk+2,fj+1,fi  ) + P3(fk  ,fj+1,fi  ))
                            + Wzz + Uxyz - Vxyz)
        / a3(m,fk+1,fj+1,fi  );
    b.x3f(m,fk+1,fj+1,fi+1) = (0.5*(P3(fk+2,fj+1,fi+1) + P3(fk  ,fj+1,fi+1))
                            + Wzz + Uxyz + Vxyz)
        / a3(m,fk+1,fj+1,fi+1);

  // Prolongate internal fields in 2D
  } else {
    Real tmp1 = 0.25*(P2(fk,fj+2,fi+1) - P2(fk,fj,  fi+1)
                    - P2(fk,fj+2,fi  ) + P2(fk,fj,  fi  ));
    Real tmp2 = 0.25*(P1(fk,fj,  fi  ) - P1(fk,fj,  fi+2)
                    - P1(fk,fj+1,fi  ) + P1(fk,fj+1,fi+2));
    b.x1f(m,fk,fj  ,fi+1) = (0.5*(P1(fk,fj,  fi  ) + P1(fk,fj,  fi+2)) + tmp1)
        / a1(m,fk,fj  ,fi+1);
    b.x1f(m,fk,fj+1,fi+1) = (0.5*(P1(fk,fj+1,fi  ) + P1(fk,fj+1,fi+2)) + tmp1)
        / a1(m,fk,fj+1,fi+1);
    b.x2f(m,fk,fj+1,fi  ) = (0.5*(P2(fk,fj,  fi  ) + P2(fk,fj+2,fi  )) + tmp2)
        / a2(m,fk,fj+1,fi  );
    b.x2f(m,fk,fj+1,fi+1) = (0.5*(P2(fk,fj,  fi+1) + P2(fk,fj+2,fi+1)) + tmp2)
        / a2(m,fk,fj+1,fi+1);
  }
  return;
}

template <int NGHOST>
KOKKOS_INLINE_FUNCTION
Real ProlongInterpolation(const int m, const int v, int k, int j, int i,
                            const int nx1, const int nx2, const int nx3,
                            const bool offsetk, const bool offsetj, const bool offseti,
                        const DvceArray5D<Real> &ca, const DualArray3D<Real> &weights) {
  // interpolated value at new grid point
  Real ivals = 0;

  for (int kk=0; kk<NGHOST+1; kk++) {
    for (int jj=0; jj<NGHOST+1; jj++) {
      for (int ii=0; ii<NGHOST+1; ii++) {
        int wghti = (offseti) ? NGHOST-ii : ii;
        int wghtj = (offsetj) ? NGHOST-jj : jj;
        int wghtk = (offsetk) ? NGHOST-kk : kk;
        ivals += weights.d_view(wghtk,wghtj,wghti)*ca(m,v,
                    k-NGHOST/2+kk,j-NGHOST/2+jj,i-NGHOST/2+ii);
      }
    }
  }

  return ivals;
}

//----------------------------------------------------------------------------------------
//! \fn HighOrderProlongCC()
//! \brief high-order prolongation operator for cell-centered variables

template <int NGHOST>
KOKKOS_INLINE_FUNCTION
void HighOrderProlongCC(const int m, const int v, const int k, const int j, const int i,
               const int fk, const int fj, const int fi, const int nx1, const int nx2,
               const int nx3, const DvceArray5D<Real> &ca, const DvceArray5D<Real> &a,
               const DualArray3D<Real> &weights) {
  // stencil size for interpolator
  a(m,v,fk  ,fj  ,fi  ) = ProlongInterpolation<NGHOST>(m,v,k,j,i, nx1, nx2, nx3,
                                                        false,false,false, ca, weights);
  a(m,v,fk  ,fj  ,fi+1) = ProlongInterpolation<NGHOST>(m,v,k,j,i, nx1, nx2, nx3,
                                                        false,false, true, ca, weights);
  a(m,v,fk  ,fj+1,fi  ) = ProlongInterpolation<NGHOST>(m,v,k,j,i, nx1, nx2, nx3,
                                                        false, true,false, ca, weights);
  a(m,v,fk  ,fj+1,fi+1) = ProlongInterpolation<NGHOST>(m,v,k,j,i, nx1, nx2, nx3,
                                                        false, true, true, ca, weights);
  a(m,v,fk+1,fj  ,fi  ) = ProlongInterpolation<NGHOST>(m,v,k,j,i, nx1, nx2, nx3,
                                                         true,false,false, ca, weights);
  a(m,v,fk+1,fj  ,fi+1) = ProlongInterpolation<NGHOST>(m,v,k,j,i, nx1, nx2, nx3,
                                                         true,false, true, ca, weights);
  a(m,v,fk+1,fj+1,fi  ) = ProlongInterpolation<NGHOST>(m,v,k,j,i, nx1, nx2, nx3,
                                                         true, true,false, ca, weights);
  a(m,v,fk+1,fj+1,fi+1) = ProlongInterpolation<NGHOST>(m,v,k,j,i, nx1, nx2, nx3,
                                                         true, true, true, ca, weights);
  return;
}

#endif // MESH_PROLONGATION_HPP_

