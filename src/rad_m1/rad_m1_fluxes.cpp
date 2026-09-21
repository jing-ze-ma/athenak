//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_fluxes.cpp
//! \brief closure/limit pass and the transport fluxes of the M1 moments.
//!
//! Design sect. 3: reconstruct (E, f_i = F_i/(c E)) with PLM and rebuild F = c E f, so
//! that |f| <= 1 survives reconstruction (QUOKKA).  The face flux is then plain HLL with
//! the Skinner & Ostriker (2013) closed-form wave speeds (thick_flux = none).

#include <math.h>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "reconstruct/plm.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn M1Prim
//! \brief load the reconstruction variables (E, f_1, f_2, f_3) of one cell

KOKKOS_INLINE_FUNCTION
void M1Prim(const DvceArray5D<Real> &u, const int m, const int k, const int j,
            const int i, const Real cl, const Real efl, Real *q) {
  Real e = fmax(u(m,M1_E,k,j,i), efl);
  Real r1, r2, r3, rn;
  M1ReducedFlux(cl, e, u(m,M1_F1,k,j,i), u(m,M1_F2,k,j,i), u(m,M1_F3,k,j,i),
                r1, r2, r3, rn);
  q[0] = e;
  q[1] = r1;
  q[2] = r2;
  q[3] = r3;
}

//----------------------------------------------------------------------------------------
//! \fn M1Rebuild
//! \brief turn a reconstructed (E, f_i) state back into (E, F_i), enforcing the floor
//! and |f| <= 1 on the FACE state

KOKKOS_INLINE_FUNCTION
void M1Rebuild(const Real *q, const Real cl, const Real efl,
               Real &e, Real &ff1, Real &ff2, Real &ff3) {
  e = fmax(q[0], efl);
  Real r1 = q[1], r2 = q[2], r3 = q[3];
  Real rn = sqrt(r1*r1 + r2*r2 + r3*r3);
  if (rn > 1.0) {
    Real s = 1.0/rn;
    r1 *= s;
    r2 *= s;
    r3 *= s;
  }
  ff1 = cl*e*r1;
  ff2 = cl*e*r2;
  ff3 = cl*e*r3;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::ApplyClosureLimits
//! \brief enforce E >= e_floor and |F| <= c E on the evolved state, over the whole
//! array (active cells and ghosts alike) so that the reconstruction below only ever
//! sees admissible states.

TaskStatus RadiationM1::ApplyClosureLimits(Driver *pdrive, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  auto u0_ = u0;
  Real cl = c_light;
  Real efl = e_floor;

  par_for("m1_limits", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = u0_(m,M1_E,k,j,i);
    Real f1 = u0_(m,M1_F1,k,j,i);
    Real f2 = u0_(m,M1_F2,k,j,i);
    Real f3 = u0_(m,M1_F3,k,j,i);
    if (M1ApplyLimits(cl, efl, e, f1, f2, f3)) {
      u0_(m,M1_E,k,j,i) = e;
      u0_(m,M1_F1,k,j,i) = f1;
      u0_(m,M1_F2,k,j,i) = f2;
      u0_(m,M1_F3,k,j,i) = f3;
    }
  });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::CalculateFluxes
//! \brief PLM + HLL fluxes of (E, F_i) on all cell faces

TaskStatus RadiationM1::CalculateFluxes(Driver *pdrive, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;

  auto u0_ = u0;
  auto flx1 = uflx.x1f;
  auto flx2 = uflx.x2f;
  auto flx3 = uflx.x3f;
  Real cl = c_light;
  Real ch = chat;
  Real efl = e_floor;
  bool edd = eddington;

  //--------------------------------------------------------------------------------- x1
  par_for("m1_flx1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real qa[4], qb[4], qc[4], qd[4], ql[4], qr[4], dum;
    M1Prim(u0_, m, k, j, i-2, cl, efl, qa);
    M1Prim(u0_, m, k, j, i-1, cl, efl, qb);
    M1Prim(u0_, m, k, j, i  , cl, efl, qc);
    M1Prim(u0_, m, k, j, i+1, cl, efl, qd);
    for (int n=0; n<M1_NVAR; ++n) {
      PLM(qa[n], qb[n], qc[n], ql[n], dum);
      PLM(qb[n], qc[n], qd[n], dum, qr[n]);
    }
    Real el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx[4];
    M1Rebuild(ql, cl, efl, el, fl1, fl2, fl3);
    M1Rebuild(qr, cl, efl, er, fr1, fr2, fr3);
    M1HLLFlux(1, cl, ch, edd, el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx);
    for (int n=0; n<M1_NVAR; ++n) {
      flx1(m,n,k,j,i) = flx[n];
    }
  });

  //--------------------------------------------------------------------------------- x2
  if (multi_d) {
    par_for("m1_flx2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real qa[4], qb[4], qc[4], qd[4], ql[4], qr[4], dum;
      M1Prim(u0_, m, k, j-2, i, cl, efl, qa);
      M1Prim(u0_, m, k, j-1, i, cl, efl, qb);
      M1Prim(u0_, m, k, j  , i, cl, efl, qc);
      M1Prim(u0_, m, k, j+1, i, cl, efl, qd);
      for (int n=0; n<M1_NVAR; ++n) {
        PLM(qa[n], qb[n], qc[n], ql[n], dum);
        PLM(qb[n], qc[n], qd[n], dum, qr[n]);
      }
      Real el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx[4];
      M1Rebuild(ql, cl, efl, el, fl1, fl2, fl3);
      M1Rebuild(qr, cl, efl, er, fr1, fr2, fr3);
      M1HLLFlux(2, cl, ch, edd, el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx);
      for (int n=0; n<M1_NVAR; ++n) {
        flx2(m,n,k,j,i) = flx[n];
      }
    });
  }

  //--------------------------------------------------------------------------------- x3
  if (three_d) {
    par_for("m1_flx3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real qa[4], qb[4], qc[4], qd[4], ql[4], qr[4], dum;
      M1Prim(u0_, m, k-2, j, i, cl, efl, qa);
      M1Prim(u0_, m, k-1, j, i, cl, efl, qb);
      M1Prim(u0_, m, k  , j, i, cl, efl, qc);
      M1Prim(u0_, m, k+1, j, i, cl, efl, qd);
      for (int n=0; n<M1_NVAR; ++n) {
        PLM(qa[n], qb[n], qc[n], ql[n], dum);
        PLM(qb[n], qc[n], qd[n], dum, qr[n]);
      }
      Real el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx[4];
      M1Rebuild(ql, cl, efl, el, fl1, fl2, fl3);
      M1Rebuild(qr, cl, efl, er, fr1, fr2, fr3);
      M1HLLFlux(3, cl, ch, edd, el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx);
      for (int n=0; n<M1_NVAR; ++n) {
        flx3(m,n,k,j,i) = flx[n];
      }
    });
  }

  return TaskStatus::complete;
}

} // namespace radm1
