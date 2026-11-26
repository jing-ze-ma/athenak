#ifndef RAD_VET_RAD_VET_ROWS_HPP_
#define RAD_VET_RAD_VET_ROWS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_rows.hpp
//! \brief device functions of the operator's rows and of the sources' divergence.

#include "rad_vet/rad_vet_faces.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void CleanStAdd
//! \brief adds f D_ab E of the cell at offset (OK, OJ, OI) from (k,j,i) to its slot of
//! row c

template <int OK, int OJ, int OI>
KOKKOS_FORCEINLINE_FUNCTION
void CleanStAdd(const CleanArgs &x, const int ab, const int m, const int k, const int j,
                const int i, const Real f, Real (&c)[CLN_ST]) {
  c[CleanStSlot(OK, OJ, OI)] += f*x.cs(m,ab,k+OK,j+OJ,i+OI);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStCross
//! \brief row terms (times cr) of the transverse derivative along B of the off-diagonal
//! stress at face SD (0 low, 1 high) normal to D: CleanFaceStress with constant slots

template <int D, int SD, int B>
KOKKOS_FORCEINLINE_FUNCTION
void CleanStCross(const CleanArgs &x, const int m, const int k, const int j, const int i,
                  const Real cr, Real (&c)[CLN_ST]) {
  constexpr int ek = (D == 2), ej = (D == 1), ei = (D == 0);
  constexpr int bk = (B == 2), bj = (B == 1), bi = (B == 0);
  constexpr int ab = (D + B == 1) ? 3 : ((D + B == 2) ? 4 : 5);   // D12, D13, D23
  // (the cells L, R of the face)
  constexpr int lk = (SD - 1)*ek, lj = (SD - 1)*ej, li = (SD - 1)*ei;
  constexpr int rk = SD*ek, rj = SD*ej, ri = SD*ei;
  const int kf = k + rk, jf = j + rj, iff = i + ri;
  // (one sided where a side's ghosts do not carry the field, CleanTransSide)
  const int up = CleanTransSide(x, B, m, kf, jf, iff, 1);
  const int dn = CleanTransSide(x, B, m, kf, jf, iff, 0);
  if constexpr (D == 2) {
    if (x.sph) {   // (torque form, CleanFaceStress)
      const Real vr = x.vol(m,kf,jf,iff)*x.x1v(m,iff)*x.tg(m,2,jf);
      if constexpr (B == 0) {
        const Real s = x.tg(m,2,jf)/vr*cr*0.25;
        const Real fp = x.area[0](m,kf,jf,iff+1)*x.x1f(m,iff+1)*s;
        const Real fm = -x.area[0](m,kf,jf,iff)*x.x1f(m,iff)*s;
        const Real f0 = (up ? fp : 2.0*fp) + (dn ? fm : 2.0*fm);
        CleanStAdd<rk, 0, 1>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        CleanStAdd<lk, 0, 1>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        CleanStAdd<rk, 0, 0>(x, ab, m, k, j, i, f0, c);
        CleanStAdd<lk, 0, 0>(x, ab, m, k, j, i, f0, c);
        CleanStAdd<rk, 0, -1>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
        CleanStAdd<lk, 0, -1>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
      } else {
        const Real s = x.x1v(m,iff)/vr*cr*0.25;
        const Real fp = x.area[1](m,kf,jf+1,iff)*x.tg(m,3,jf+1)*s;
        const Real fm = -x.area[1](m,kf,jf,iff)*x.tg(m,3,jf)*s;
        const Real f0 = (up ? fp : 2.0*fp) + (dn ? fm : 2.0*fm);
        CleanStAdd<rk, 1, 0>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        CleanStAdd<lk, 1, 0>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        CleanStAdd<rk, 0, 0>(x, ab, m, k, j, i, f0, c);
        CleanStAdd<lk, 0, 0>(x, ab, m, k, j, i, f0, c);
        CleanStAdd<rk, -1, 0>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
        CleanStAdd<lk, -1, 0>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
      }
      return;
    }
  }
  if (up + dn == 0) return;
  const Real rb = !x.sph ? 1.0 : ((D == 0) ? x.x1f(m,iff) : x.x1v(m,iff));
  const Real w = (up + dn == 2) ? 1.0 : 2.0;
  const Real ib = !x.sph ? w*x.idx[B] :
                  ((B == 0) ? 2.0/(x.x1v(m,iff+up) - x.x1v(m,iff-dn)) :
                  ((B == 1) ? 2.0/(rb*(x.x2v(m,jf+up) - x.x2v(m,jf-dn))) :
                  w*x.idx[2]/(rb*sin((D == 1) ? x.x2f(m,jf) : x.x2v(m,jf)))));
  const Real q = 0.25*ib*cr;
  if (up) {
    CleanStAdd<rk + bk, rj + bj, ri + bi>(x, ab, m, k, j, i, q, c);
    CleanStAdd<lk + bk, lj + bj, li + bi>(x, ab, m, k, j, i, q, c);
  } else {
    CleanStAdd<rk, rj, ri>(x, ab, m, k, j, i, q, c);
    CleanStAdd<lk, lj, li>(x, ab, m, k, j, i, q, c);
  }
  if (dn) {
    CleanStAdd<rk - bk, rj - bj, ri - bi>(x, ab, m, k, j, i, -q, c);
    CleanStAdd<lk - bk, lj - bj, li - bi>(x, ab, m, k, j, i, -q, c);
  } else {
    CleanStAdd<rk, rj, ri>(x, ab, m, k, j, i, -q, c);
    CleanStAdd<lk, lj, li>(x, ab, m, k, j, i, -q, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStFace
//! \brief row terms of face SD normal to D in h div F: al, ar and the cross factor
//! (CLF_QL) with the off-diagonal stress

template <int D, int SD>
KOKKOS_FORCEINLINE_FUNCTION
void CleanStFace(const CleanArgs &x, const int m, const int k, const int j, const int i,
                 Real (&c)[CLN_ST]) {
  constexpr int ek = (D == 2), ej = (D == 1), ei = (D == 0);
  const int kf = k + SD*ek, jf = j + SD*ej, iff = i + SD*ei;
  const Real wd = (SD ? x.h : -x.h)*CleanDivW(x, D, m, k, j, i, SD);
  c[CleanStSlot((SD - 1)*ek, (SD - 1)*ej, (SD - 1)*ei)] += wd*x.fc[D](m,CLF_QL,kf,jf,iff);
  c[CleanStSlot(SD*ek, SD*ej, SD*ei)] += wd*x.fc[D](m,CLF_QL+1,kf,jf,iff);
  // (the cubed sphere's lateral stress comes through its sums Y, CleanApply)
  if (!x.multi_d || x.cub || CleanBoundaryFace(x, D, m, kf, jf, iff)) return;
  const Real cr = wd*x.fc[D](m,CLF_QL+2,kf,jf,iff);
  if constexpr (D != 0) CleanStCross<D, SD, 0>(x, m, k, j, i, cr, c);
  if constexpr (D != 1) CleanStCross<D, SD, 1>(x, m, k, j, i, cr, c);
  if constexpr (D != 2) {
    if (x.three_d) CleanStCross<D, SD, 2>(x, m, k, j, i, cr, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStencilRow
//! \brief c = the row of cell (k,j,i) of CleanApply's operator (x.fwc)

KOKKOS_INLINE_FUNCTION
void CleanStencilRow(const CleanArgs &x, const int m, const int k, const int j,
                     const int i, Real (&c)[CLN_ST]) {
  for (int s = 0; s < CLN_ST; ++s) c[s] = 0.0;
  const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
  c[0] = 1.0 + x.h*x.c*x.op(m,VET_OP_E,k,j,i)*qr + x.w(m,CLW_PD,k,j,i);
  CleanStFace<0, 0>(x, m, k, j, i, c);
  CleanStFace<0, 1>(x, m, k, j, i, c);
  if (x.multi_d) {
    CleanStFace<1, 0>(x, m, k, j, i, c);
    CleanStFace<1, 1>(x, m, k, j, i, c);
  }
  if (x.three_d) {
    CleanStFace<2, 0>(x, m, k, j, i, c);
    CleanStFace<2, 1>(x, m, k, j, i, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanSourceDivergence
//! \brief h div of the stored face sources (CLF_SR) of cell (k,j,i)

KOKKOS_INLINE_FUNCTION
Real CleanSourceDivergence(const CleanArgs &x, const int m, const int k, const int j,
                           const int i) {
  Real d = 0.0;
  for (int n = 0; n < 3; ++n) {
    if (n >= x.nf) break;
    const int dk = (n == 2), dj = (n == 1), di = (n == 0);
    d += CleanDivF(x, n, m, k, j, i, x.fc[n](m,CLF_SR,k,j,i),
                   x.fc[n](m,CLF_SR,k+dk,j+dj,i+di));
  }
  return x.h*d;
}

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_ROWS_HPP_
