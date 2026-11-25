#ifndef RAD_VET_RAD_VET_CLOSURE_HPP_
#define RAD_VET_RAD_VET_CLOSURE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_closure.hpp
//! \brief Levermore's closure factor chi, the physical-boundary ghost fill (device inline).
//! E, F_i lab-frame moments, f_i = F_i/(c E), |f| <= 1, n_i = f_i/|f|:
//!   P_ij = E [ (1-chi)/2 delta_ij + (3 chi - 1)/2 n_i n_j ]

#include <math.h>
#include "athena.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn VetChi
//! \brief Levermore's (1984) chi(|f|) = (3 + 4 f^2)/(5 + 2 sqrt(4 - 3 f^2)) (problem
//! generators: initial states and exact solutions)

KOKKOS_INLINE_FUNCTION
Real VetChi(const Real fnorm) {
  Real f2 = fnorm*fnorm;
  if (f2 > 1.0) {f2 = 1.0;}
  Real s = sqrt(fmax(4.0 - 3.0*f2, 0.0));
  return (3.0 + 4.0*f2)/(5.0 + 2.0*s);
}

//----------------------------------------------------------------------------------------
//! \fn VetFillGhost
//! \brief fill one physical-boundary ghost cell from the conserved moments of active
//! cell (ka,ja,ia); shared by RadVetBCs and user BCs.
//! mode 0 = copy (outflow), 1 = reflect (flip normal flux), 2 = vacuum (dark ghost).
//! inorm = 1,2,3 is the face-normal flux component; sgn = +1 if outgoing is +x.

KOKKOS_INLINE_FUNCTION
void VetFillGhost(const DvceArray5D<Real> &u, const int m,
                 const int kg, const int jg, const int ig,
                 const int ka, const int ja, const int ia,
                 const int inorm, const int mode, const Real sgn,
                 const Real cl, const Real efl) {
  Real e  = u(m,VET_E, ka,ja,ia);
  Real f1 = u(m,VET_F1,ka,ja,ia);
  Real f2 = u(m,VET_F2,ka,ja,ia);
  Real f3 = u(m,VET_F3,ka,ja,ia);
  if (mode == 1) {
    if (inorm == 1) {
      f1 = -f1;
    } else if (inorm == 2) {
      f2 = -f2;
    } else {
      f3 = -f3;
    }
  } else if (mode == 2) {
    // vacuum: dark ghost (E = e_floor, F = 0), no incoming radiation. A beam at f = 1
    // still leaves (HLL b_L = 0 for a dark state); the HLL flux then gives a
    // Marshak-like F-E relation at the face. Copying the interior would let E float.
    (void) sgn;
    e = efl;
    f1 = 0.0;
    f2 = 0.0;
    f3 = 0.0;
  }
  (void) cl;
  u(m,VET_E, kg,jg,ig) = e;
  u(m,VET_F1,kg,jg,ig) = f1;
  u(m,VET_F2,kg,jg,ig) = f2;
  u(m,VET_F3,kg,jg,ig) = f3;
}

} // namespace radvet
#endif // RAD_VET_RAD_VET_CLOSURE_HPP_
