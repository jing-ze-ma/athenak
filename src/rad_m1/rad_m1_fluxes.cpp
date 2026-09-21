//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_fluxes.cpp
//! \brief closure/limit pass and the transport fluxes of the M1 moments.
//!
//! Design sect. 3: reconstruct (E, f_i = F_i/(c E)) and rebuild F = c E f, so that
//! |f| <= 1 survives reconstruction (QUOKKA).  The face flux is then plain HLL with
//! the Skinner & Ostriker (2013) closed-form wave speeds (thick_flux = none).
//!
//! MILESTONE 1c-B: <rad_m1>/reconstruct = dc | plm | ppm4 | ppmx | wenoz, using the
//! scalar kernels of src/reconstruct/.  The three high-order ones read a 5-cell stencil
//! per face state, so the face at i needs cells i-3..i+2 and <mesh>/nghost >= 3; dc and
//! plm keep the 4-cell load of 1c-A and are bit-identical to it.  In the ap_hll flux
//! only `dc` selects Berthon's form: every reconstruction that puts a second-order (or
//! better) polynomial on the face makes the face jump O(dx^2) and needs the alpha2
//! blend, exactly as plm does.
//!
//! MILESTONE 1c-B, the advective split: the velocity that builds the enthalpy flux
//! A = v E + v.P is reconstructed to the face with the SAME method as (E, f_i) when
//! <rad_m1>/split_vel = recon (the default), instead of each side's cell velocity
//! (= split_vel is cell, the 1c-A behaviour).  The upwind direction is the mean of the
//! two reconstructed face-normal velocities.

#include <math.h>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "reconstruct/plm.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn M1ReconFace
//! \brief One scalar reconstructed to the two sides of the face between cell c (= i-1)
//! and cell d (= i), from the six cell values a..f = q(i-3) .. q(i+2).  ql is the state
//! on the LEFT of the face (the right edge of cell i-1), qr the state on the RIGHT (the
//! left edge of cell i).  a and f are read only by the high-order methods.

KOKKOS_INLINE_FUNCTION
void M1ReconFace(const int meth, const Real a, const Real b, const Real c,
                 const Real d, const Real e, const Real f, Real &ql, Real &qr) {
  Real dum;
  if (meth == M1_RECON_DC) {
    ql = c;
    qr = d;
  } else if (meth == M1_RECON_PLM) {
    PLM(b, c, d, ql, dum);
    PLM(c, d, e, dum, qr);
  } else if (meth == M1_RECON_PPM4) {
    PPM4(a, b, c, d, e, ql, dum);
    PPM4(b, c, d, e, f, dum, qr);
  } else if (meth == M1_RECON_PPMX) {
    PPMX(a, b, c, d, e, ql, dum);
    PPMX(b, c, d, e, f, dum, qr);
  } else {
    WENOZ(a, b, c, d, e, ql, dum);
    WENOZ(b, c, d, e, f, dum, qr);
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1Vel
//! \brief the three velocity components of one cell, from hydro's CONSERVED u0 (the
//! coupling writes u0 inside the substep, so w0 is one ConToPrim behind -- the same rule
//! the opacity kernel follows).

KOKKOS_INLINE_FUNCTION
void M1Vel(const DvceArray5D<Real> &uh, const int m, const int k, const int j,
           const int i, Real *v) {
  Real id = 1.0/fmax(uh(m,IDN,k,j,i), 1.0e-300);
  v[0] = uh(m,IM1,k,j,i)*id;
  v[1] = uh(m,IM2,k,j,i)*id;
  v[2] = uh(m,IM3,k,j,i)*id;
}

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
//! \brief reconstructed HLL fluxes of (E, F_i) on all cell faces

TaskStatus RadiationM1::CalculateFluxes(Driver *pdrive, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;

  auto u0_ = u0;
  auto opac_ = opac;
  auto flx1 = uflx.x1f;
  auto flx2 = uflx.x2f;
  auto flx3 = uflx.x3f;
  auto &mbsize = pmy_pack->pmb->mb_size;
  Real cl = c_light;
  Real ch = chat;
  Real efl = e_floor;
  bool edd = eddington;
  int rmeth = recon_code;
  bool dc = (rmeth == M1_RECON_DC);
  // stencil slots actually loaded, of the six (i-3 .. i+2).  dc and plm read four, so
  // with <mesh>/nghost = 2 the outer two are never touched (they stay zero).
  int s0 = (rmeth >= M1_RECON_PPM4) ? 0 : 1;
  int s1 = (rmeth >= M1_RECON_PPM4) ? 5 : 4;
  // the thick-limit correction needs a face opacity; with none stored there is nothing
  // to correct and the flux is plain HLL (and bit-identical to milestone 1a)
  int thick = (opac_zero) ? M1_THICK_NONE : thick_flux;
  Real spref = scaled_pref;
  int apform = ap_form;
  // the advective enthalpy-flux split needs the hydro velocity field; with no <hydro>
  // there is no medium to move and the switch is forced off (u0 is captured as a dummy
  // in that case, and never read)
  bool split = advect_split && (pmy_pack->phydro != nullptr);
  bool vrec = split && split_vel_recon;
  auto uh = (pmy_pack->phydro != nullptr) ? pmy_pack->phydro->u0 : u0;

  //--------------------------------------------------------------------------------- x1
  par_for("m1_flx1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real qs[6][4], vs[6][3];
    for (int s=0; s<6; ++s) {
      for (int n=0; n<4; ++n) {
        qs[s][n] = 0.0;
      }
      for (int n=0; n<3; ++n) {
        vs[s][n] = 0.0;
      }
    }
    for (int s=s0; s<=s1; ++s) {
      M1Prim(u0_, m, k, j, i-3+s, cl, efl, qs[s]);
      if (vrec) {
        M1Vel(uh, m, k, j, i-3+s, vs[s]);
      }
    }
    if (split && !vrec) {
      M1Vel(uh, m, k, j, i-1, vs[2]);
      M1Vel(uh, m, k, j, i, vs[3]);
    }
    Real ql[4], qr[4];
    for (int n=0; n<M1_NVAR; ++n) {
      M1ReconFace(rmeth, qs[0][n], qs[1][n], qs[2][n], qs[3][n], qs[4][n], qs[5][n],
                  ql[n], qr[n]);
    }
    Real vl[3] = {0.0, 0.0, 0.0}, vr[3] = {0.0, 0.0, 0.0};
    if (split) {
      for (int n=0; n<3; ++n) {
        if (vrec) {
          M1ReconFace(rmeth, vs[0][n], vs[1][n], vs[2][n], vs[3][n], vs[4][n], vs[5][n],
                      vl[n], vr[n]);
        } else {
          vl[n] = vs[2][n];
          vr[n] = vs[3][n];
        }
      }
    }
    Real tauf = 0.0;
    if (thick != M1_THICK_NONE) {
      tauf = 0.5*(opac_(m,M1_OP_T,k,j,i-1) + opac_(m,M1_OP_T,k,j,i))
             *mbsize.d_view(m).dx1;
    }
    Real el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx[4];
    M1Rebuild(ql, cl, efl, el, fl1, fl2, fl3);
    M1Rebuild(qr, cl, efl, er, fr1, fr2, fr3);
    M1HLLFlux(1, cl, ch, edd, el, fl1, fl2, fl3, er, fr1, fr2, fr3,
              thick, tauf, spref, dc, qs[2][0], qs[3][0], apform, split, vl, vr, flx);
    for (int n=0; n<M1_NVAR; ++n) {
      flx1(m,n,k,j,i) = flx[n];
    }
  });

  //--------------------------------------------------------------------------------- x2
  if (multi_d) {
    par_for("m1_flx2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real qs[6][4], vs[6][3];
      for (int s=0; s<6; ++s) {
        for (int n=0; n<4; ++n) {
          qs[s][n] = 0.0;
        }
        for (int n=0; n<3; ++n) {
          vs[s][n] = 0.0;
        }
      }
      for (int s=s0; s<=s1; ++s) {
        M1Prim(u0_, m, k, j-3+s, i, cl, efl, qs[s]);
        if (vrec) {
          M1Vel(uh, m, k, j-3+s, i, vs[s]);
        }
      }
      if (split && !vrec) {
        M1Vel(uh, m, k, j-1, i, vs[2]);
        M1Vel(uh, m, k, j, i, vs[3]);
      }
      Real ql[4], qr[4];
      for (int n=0; n<M1_NVAR; ++n) {
        M1ReconFace(rmeth, qs[0][n], qs[1][n], qs[2][n], qs[3][n], qs[4][n], qs[5][n],
                    ql[n], qr[n]);
      }
      Real vl[3] = {0.0, 0.0, 0.0}, vr[3] = {0.0, 0.0, 0.0};
      if (split) {
        for (int n=0; n<3; ++n) {
          if (vrec) {
            M1ReconFace(rmeth, vs[0][n], vs[1][n], vs[2][n], vs[3][n], vs[4][n],
                        vs[5][n], vl[n], vr[n]);
          } else {
            vl[n] = vs[2][n];
            vr[n] = vs[3][n];
          }
        }
      }
      Real tauf = 0.0;
      if (thick != M1_THICK_NONE) {
        tauf = 0.5*(opac_(m,M1_OP_T,k,j-1,i) + opac_(m,M1_OP_T,k,j,i))
               *mbsize.d_view(m).dx2;
      }
      Real el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx[4];
      M1Rebuild(ql, cl, efl, el, fl1, fl2, fl3);
      M1Rebuild(qr, cl, efl, er, fr1, fr2, fr3);
      M1HLLFlux(2, cl, ch, edd, el, fl1, fl2, fl3, er, fr1, fr2, fr3,
                thick, tauf, spref, dc, qs[2][0], qs[3][0], apform, split, vl, vr, flx);
      for (int n=0; n<M1_NVAR; ++n) {
        flx2(m,n,k,j,i) = flx[n];
      }
    });
  }

  //--------------------------------------------------------------------------------- x3
  if (three_d) {
    par_for("m1_flx3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real qs[6][4], vs[6][3];
      for (int s=0; s<6; ++s) {
        for (int n=0; n<4; ++n) {
          qs[s][n] = 0.0;
        }
        for (int n=0; n<3; ++n) {
          vs[s][n] = 0.0;
        }
      }
      for (int s=s0; s<=s1; ++s) {
        M1Prim(u0_, m, k-3+s, j, i, cl, efl, qs[s]);
        if (vrec) {
          M1Vel(uh, m, k-3+s, j, i, vs[s]);
        }
      }
      if (split && !vrec) {
        M1Vel(uh, m, k-1, j, i, vs[2]);
        M1Vel(uh, m, k, j, i, vs[3]);
      }
      Real ql[4], qr[4];
      for (int n=0; n<M1_NVAR; ++n) {
        M1ReconFace(rmeth, qs[0][n], qs[1][n], qs[2][n], qs[3][n], qs[4][n], qs[5][n],
                    ql[n], qr[n]);
      }
      Real vl[3] = {0.0, 0.0, 0.0}, vr[3] = {0.0, 0.0, 0.0};
      if (split) {
        for (int n=0; n<3; ++n) {
          if (vrec) {
            M1ReconFace(rmeth, vs[0][n], vs[1][n], vs[2][n], vs[3][n], vs[4][n],
                        vs[5][n], vl[n], vr[n]);
          } else {
            vl[n] = vs[2][n];
            vr[n] = vs[3][n];
          }
        }
      }
      Real tauf = 0.0;
      if (thick != M1_THICK_NONE) {
        tauf = 0.5*(opac_(m,M1_OP_T,k-1,j,i) + opac_(m,M1_OP_T,k,j,i))
               *mbsize.d_view(m).dx3;
      }
      Real el, fl1, fl2, fl3, er, fr1, fr2, fr3, flx[4];
      M1Rebuild(ql, cl, efl, el, fl1, fl2, fl3);
      M1Rebuild(qr, cl, efl, er, fr1, fr2, fr3);
      M1HLLFlux(3, cl, ch, edd, el, fl1, fl2, fl3, er, fr1, fr2, fr3,
                thick, tauf, spref, dc, qs[2][0], qs[3][0], apform, split, vl, vr, flx);
      for (int n=0; n<M1_NVAR; ++n) {
        flx3(m,n,k,j,i) = flx[n];
      }
    });
  }

  return TaskStatus::complete;
}

} // namespace radm1
