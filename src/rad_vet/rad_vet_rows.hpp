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
//! \fn void StAdd
//! \brief adds f D_ab E of the cell at offset (OK, OJ, OI) from (k,j,i) to its slot of
//! row c

template <int OK, int OJ, int OI>
KOKKOS_FORCEINLINE_FUNCTION
void StAdd(const SolveArgs &x, const int ab, const int m, const int k, const int j,
           const int i, const Real f, Real (&c)[CLN_ST]) {
  c[StSlot(OK, OJ, OI)] += f*x.cs(m,ab,k+OK,j+OJ,i+OI);
}

//----------------------------------------------------------------------------------------
//! \fn void StCross
//! \brief row terms (times cr) of the transverse derivative along B of the off-diagonal
//! stress at face SD (0 low, 1 high) normal to D: FaceStress with constant slots

template <int D, int SD, int B>
KOKKOS_FORCEINLINE_FUNCTION
void StCross(const SolveArgs &x, const int m, const int k, const int j, const int i,
             const Real cr, Real (&c)[CLN_ST]) {
  constexpr int ek = (D == 2), ej = (D == 1), ei = (D == 0);
  constexpr int bk = (B == 2), bj = (B == 1), bi = (B == 0);
  constexpr int ab = (D + B == 1) ? 3 : ((D + B == 2) ? 4 : 5);   // D12, D13, D23
  // (the cells L, R of the face)
  constexpr int lk = (SD - 1)*ek, lj = (SD - 1)*ej, li = (SD - 1)*ei;
  constexpr int rk = SD*ek, rj = SD*ej, ri = SD*ei;
  const int kf = k + rk, jf = j + rj, iff = i + ri;
  // (one sided where a side's ghosts do not carry the field, TransSide)
  const int up = TransSide(x, B, m, kf, jf, iff, 1);
  const int dn = TransSide(x, B, m, kf, jf, iff, 0);
  if constexpr (D == 2) {
    if (x.sph) {   // (torque form, FaceStress)
      const Real vr = x.vol(m,kf,jf,iff)*x.x1v(m,iff)*x.tg(m,2,jf);
      if constexpr (B == 0) {
        const Real s = x.tg(m,2,jf)/vr*cr*0.25;
        const Real fp = x.area[0](m,kf,jf,iff+1)*x.x1f(m,iff+1)*s;
        const Real fm = -x.area[0](m,kf,jf,iff)*x.x1f(m,iff)*s;
        const Real f0 = (up ? fp : 2.0*fp) + (dn ? fm : 2.0*fm);
        StAdd<rk, 0, 1>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        StAdd<lk, 0, 1>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        StAdd<rk, 0, 0>(x, ab, m, k, j, i, f0, c);
        StAdd<lk, 0, 0>(x, ab, m, k, j, i, f0, c);
        StAdd<rk, 0, -1>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
        StAdd<lk, 0, -1>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
      } else {
        const Real s = x.x1v(m,iff)/vr*cr*0.25;
        const Real fp = x.area[1](m,kf,jf+1,iff)*x.tg(m,3,jf+1)*s;
        const Real fm = -x.area[1](m,kf,jf,iff)*x.tg(m,3,jf)*s;
        const Real f0 = (up ? fp : 2.0*fp) + (dn ? fm : 2.0*fm);
        StAdd<rk, 1, 0>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        StAdd<lk, 1, 0>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        StAdd<rk, 0, 0>(x, ab, m, k, j, i, f0, c);
        StAdd<lk, 0, 0>(x, ab, m, k, j, i, f0, c);
        StAdd<rk, -1, 0>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
        StAdd<lk, -1, 0>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
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
    StAdd<rk + bk, rj + bj, ri + bi>(x, ab, m, k, j, i, q, c);
    StAdd<lk + bk, lj + bj, li + bi>(x, ab, m, k, j, i, q, c);
  } else {
    StAdd<rk, rj, ri>(x, ab, m, k, j, i, q, c);
    StAdd<lk, lj, li>(x, ab, m, k, j, i, q, c);
  }
  if (dn) {
    StAdd<rk - bk, rj - bj, ri - bi>(x, ab, m, k, j, i, -q, c);
    StAdd<lk - bk, lj - bj, li - bi>(x, ab, m, k, j, i, -q, c);
  } else {
    StAdd<rk, rj, ri>(x, ab, m, k, j, i, -q, c);
    StAdd<lk, lj, li>(x, ab, m, k, j, i, -q, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void StFace
//! \brief row terms of face SD normal to D in h div F: al, ar and the cross factor
//! (CLF_QL) with the off-diagonal stress

template <int D, int SD>
KOKKOS_FORCEINLINE_FUNCTION
void StFace(const SolveArgs &x, const int m, const int k, const int j, const int i,
            Real (&c)[CLN_ST]) {
  constexpr int ek = (D == 2), ej = (D == 1), ei = (D == 0);
  const int kf = k + SD*ek, jf = j + SD*ej, iff = i + SD*ei;
  const Real wd = (SD ? x.h : -x.h)*DivW(x, D, m, k, j, i, SD);
  c[StSlot((SD - 1)*ek, (SD - 1)*ej, (SD - 1)*ei)] += wd*x.fc[D](m,CLF_QL,kf,jf,iff);
  c[StSlot(SD*ek, SD*ej, SD*ei)] += wd*x.fc[D](m,CLF_QL+1,kf,jf,iff);
  // (the cubed sphere's lateral stress comes through its sums Y, ApplyOperator)
  if (!x.multi_d || x.cub || IsBoundaryFace(x, D, m, kf, jf, iff)) return;
  const Real cr = wd*x.fc[D](m,CLF_QL+2,kf,jf,iff);
  if constexpr (D != 0) StCross<D, SD, 0>(x, m, k, j, i, cr, c);
  if constexpr (D != 1) StCross<D, SD, 1>(x, m, k, j, i, cr, c);
  if constexpr (D != 2) {
    if (x.three_d) StCross<D, SD, 2>(x, m, k, j, i, cr, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void StencilRow
//! \brief c = the row of cell (k,j,i) of ApplyOperator's operator (x.fwc)

KOKKOS_INLINE_FUNCTION
void StencilRow(const SolveArgs &x, const int m, const int k, const int j,
                const int i, Real (&c)[CLN_ST]) {
  for (int s = 0; s < CLN_ST; ++s) c[s] = 0.0;
  const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
  c[0] = 1.0 + x.h*x.c*x.op(m,VET_OP_E,k,j,i)*qr + x.w(m,CLW_PD,k,j,i);
  StFace<0, 0>(x, m, k, j, i, c);
  StFace<0, 1>(x, m, k, j, i, c);
  if (x.multi_d) {
    StFace<1, 0>(x, m, k, j, i, c);
    StFace<1, 1>(x, m, k, j, i, c);
  }
  if (x.three_d) {
    StFace<2, 0>(x, m, k, j, i, c);
    StFace<2, 1>(x, m, k, j, i, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real SourceDivergence
//! \brief h div of the stored face sources (CLF_SR) of cell (k,j,i)

KOKKOS_INLINE_FUNCTION
Real SourceDivergence(const SolveArgs &x, const int m, const int k, const int j,
                      const int i) {
  Real d = 0.0;
  for (int n = 0; n < 3; ++n) {
    if (n >= x.nf) break;
    const int dk = (n == 2), dj = (n == 1), di = (n == 0);
    d += DivF(x, n, m, k, j, i, x.fc[n](m,CLF_SR,k,j,i),
              x.fc[n](m,CLF_SR,k+dk,j+dj,i+di));
  }
  return x.h*d;
}

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_ROWS_HPP_
