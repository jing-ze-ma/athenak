#ifndef BVALS_CS_SEAM_STATES_HPP_
#define BVALS_CS_SEAM_STATES_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cs_seam_states.hpp
//! \brief the cubed-sphere SHARED seam face states (csedge-1019 hydro, ucthshare-1102
//! MHD), one copy for both modules.
//!
//! Each panel's OWN interior-side face state at its panel-seam faces -- the PLM of its
//! edge cell over (ghost, edge, inner), exactly the reconstruction the sweep makes there
//! -- is written into a face record (DvceFaceFld5D, x2f / x3f at the seam faces) and
//! exchanged with a seam_raw MeshBoundaryValuesCC; the sweep then solves each seam face
//! with (neighbour's own, this panel's own) states, i.e. ONE Riemann problem per face.
//!   CsSeamPLMFill        PLM of nv variables of q into the record at offset out_off
//!                        (rows i in [ia,ib]); floor_n >= 0 floors that variable at pfl.
//!                        mlat: the hydro merged-latitude group PLM (hydro_mergelat.hpp).
//!   CsSeamCovariant      the tangential pair at slot p -> COVARIANT components on the
//!                        face's unit basis (mode 0: contravariant velocity, mode 1: the
//!                        orthonormal (B.e_xi, B.f2) of bcc), which the seam exchange
//!                        re-expresses on the neighbour's basis.
//!   CsSeamGhostPrim      in a sweep (after the gnomonic face rotation): the ghost-side
//!                        primitive / derived state of a seam face <- the record.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "reconstruct/plm.hpp"
#include "reconstruct/plm_c1.hpp"
#include "hydro/hydro_mergelat.hpp"

namespace cs_seam {

inline void CsSeamPLMFill(MeshBlockPack *pp, const DvceArray5D<Real> &q, const int nv,
                          const int out_off, const int floor_n, const Real pfl,
                          const bool plm_c1, DvceFaceFld5D<Real> &fs, const bool mlat_,
                          const DvceArray2D<int> &mlfi_, const int ia, const int ib) {
  if (nv <= 0) return;
  auto &indcs = pp->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pp->nmb_thispack - 1;
  auto &size = pp->pmb->mb_size;
  auto f2 = fs.x2f;
  auto f3 = fs.x3f;
  auto q_ = q;
  const int mlng_ = indcs.ng;
  const int mln2_ = indcs.nx2 + 2*indcs.ng;
  const int mln3_ = indcs.nx3 + 2*indcs.ng;
  // x2 (xi) seams: side 0 = the face js (x2min = -1), side 1 = the face je+1 (x2max = 1)
  par_for("cs_seamfs_x2", DevExeSpace(), 0, nmb1, ks, ke, 0, 1, ia, ib,
  KOKKOS_LAMBDA(const int m, const int k, const int sd, const int i) {
    const bool lo = (sd == 0);
    if (lo && !(size.d_view(m).x2min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x2max >= 1.0 - 1.0e-12)) return;
    const int jc = lo ? js : je, jf = lo ? js : je+1;
    const int f = mlat_ ? mlfi_(m,i) : 1;
    if (f > 1) {
      // the x2 sweep's merged branch at row jc (hydro_fluxes.cpp, mrglat-1019)
      const int jlo = js + ((jc - js)/f)*f;
      const int jm = jlo - 1, jp = jlo + f;
      const bool narrow = (f*((f + 1)/f) > mlng_ - 1);
      const bool flat = (jm < 0 || jp > mln2_ - 1) ||
                        (narrow && (jlo <= js || jlo + f - 1 >= je));
      const int sc = jc - jlo;
      for (int n=0; n<nv; ++n) {
        Real qlo, qhi;
        hydro::MlatCoarsePLM(flat ? q_(m,n,k,jc,i) : q_(m,n,k,jm,i), q_(m,n,k,jc,i),
                             flat ? q_(m,n,k,jc,i) : q_(m,n,k,jp,i), sc, f, qlo, qhi);
        Real qq = lo ? qlo : qhi;
        if (n == floor_n) qq = fmax(qq, pfl);
        f2(m,out_off+n,k,jf,i) = qq;
      }
    } else {
      for (int n=0; n<nv; ++n) {
        Real ql, qr;
        if (plm_c1) {
          PLMC1(q_(m,n,k,jc-1,i), q_(m,n,k,jc,i), q_(m,n,k,jc+1,i), ql, qr);
        } else {
          PLM(q_(m,n,k,jc-1,i), q_(m,n,k,jc,i), q_(m,n,k,jc+1,i), ql, qr);
        }
        Real qq = lo ? qr : ql;
        if (n == floor_n) qq = fmax(qq, pfl);
        f2(m,out_off+n,k,jf,i) = qq;
      }
    }
  });
  par_for("cs_seamfs_x3", DevExeSpace(), 0, nmb1, 0, 1, js, je, ia, ib,
  KOKKOS_LAMBDA(const int m, const int sd, const int j, const int i) {
    const bool lo = (sd == 0);
    if (lo && !(size.d_view(m).x3min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x3max >= 1.0 - 1.0e-12)) return;
    const int kc = lo ? ks : ke, kf = lo ? ks : ke+1;
    const int f = mlat_ ? mlfi_(m,i) : 1;
    if (f > 1) {
      // the x3 sweep's merged branch at row kc (hydro_fluxes.cpp, mrglat-1019)
      const int klo = ks + ((kc - ks)/f)*f;
      const int km = klo - 1, kp = klo + f;
      const bool narrow = (f*((f + 1)/f) > mlng_ - 1);
      const bool flat = (km < 0 || kp > mln3_ - 1) ||
                        (narrow && (klo <= ks || klo + f - 1 >= ke));
      const int sc = kc - klo;
      for (int n=0; n<nv; ++n) {
        Real qlo, qhi;
        hydro::MlatCoarsePLM(flat ? q_(m,n,kc,j,i) : q_(m,n,km,j,i), q_(m,n,kc,j,i),
                             flat ? q_(m,n,kc,j,i) : q_(m,n,kp,j,i), sc, f, qlo, qhi);
        Real qq = lo ? qlo : qhi;
        if (n == floor_n) qq = fmax(qq, pfl);
        f3(m,out_off+n,kf,j,i) = qq;
      }
    } else {
      for (int n=0; n<nv; ++n) {
        Real ql, qr;
        if (plm_c1) {
          PLMC1(q_(m,n,kc-1,j,i), q_(m,n,kc,j,i), q_(m,n,kc+1,j,i), ql, qr);
        } else {
          PLM(q_(m,n,kc-1,j,i), q_(m,n,kc,j,i), q_(m,n,kc+1,j,i), ql, qr);
        }
        Real qq = lo ? qr : ql;
        if (n == floor_n) qq = fmax(qq, pfl);
        f3(m,out_off+n,kf,j,i) = qq;
      }
    }
  });
}

// the tangential pair (p, p+1) of the record -> covariant on the face's unit basis.
// mode 0: contravariant (v^xi, v^eta) -> (v.e_xi, v.e_eta) = (v2 + c v3, v3 + c v2);
// mode 1: orthonormal (B.e_xi, B.f2) -> (B.e_xi, B.e_eta) = (b1, c b1 + s b2).
inline void CsSeamCovariant(MeshBlockPack *pp, DvceFaceFld5D<Real> &fs, const int p,
                            const int mode, const int ia, const int ib) {
  auto &indcs = pp->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pp->nmb_thispack - 1;
  auto &size = pp->pmb->mb_size;
  auto gtx = pp->pcoord->GnomonicTrigFaceXi();
  auto gte = pp->pcoord->GnomonicTrigFaceEta();
  auto f2 = fs.x2f;
  auto f3 = fs.x3f;
  par_for("cs_seamfs_cov2", DevExeSpace(), 0, nmb1, ks, ke, 0, 1, ia, ib,
  KOKKOS_LAMBDA(const int m, const int k, const int sd, const int i) {
    const bool lo = (sd == 0);
    if (lo && !(size.d_view(m).x2min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x2max >= 1.0 - 1.0e-12)) return;
    const int jf = lo ? js : je+1;
    const Real c = gtx.cs(m,k,jf);
    const Real v2 = f2(m,p,k,jf,i), v3 = f2(m,p+1,k,jf,i);
    if (mode == 0) {
      f2(m,p,k,jf,i) = v2 + c*v3;     // covariant on the face's (unit) basis
      f2(m,p+1,k,jf,i) = v3 + c*v2;
    } else {
      f2(m,p+1,k,jf,i) = c*v2 + gtx.sn(m,k,jf)*v3;
    }
  });
  par_for("cs_seamfs_cov3", DevExeSpace(), 0, nmb1, 0, 1, js, je, ia, ib,
  KOKKOS_LAMBDA(const int m, const int sd, const int j, const int i) {
    const bool lo = (sd == 0);
    if (lo && !(size.d_view(m).x3min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x3max >= 1.0 - 1.0e-12)) return;
    const int kf = lo ? ks : ke+1;
    const Real c = gte.cs(m,kf,j);
    const Real v2 = f3(m,p,kf,j,i), v3 = f3(m,p+1,kf,j,i);
    if (mode == 0) {
      f3(m,p,kf,j,i) = v2 + c*v3;
      f3(m,p+1,kf,j,i) = v3 + c*v2;
    } else {
      f3(m,p+1,kf,j,i) = c*v2 + gte.sn(m,kf,j)*v3;
    }
  });
}

// In a sweep, after the gnomonic face rotation: the ghost-side state of a seam face
// (wl at the lower face lo, wr at the upper) <- the neighbour's own record.  ivn is the
// face-normal velocity slot (x2: IVY, x3: IVZ) and ivt the along-seam one: the record
// holds covariant components, v.n = (v.e_n - c v.e_t)/s; every other slot is in the
// sweep's frame already (IVX = v.r; the along-seam slot IS the covariant one).
// soff: record offset of this row's set (0, or the radial-ghost-row set).
template <class TM, class SA, class SB>
KOKKOS_INLINE_FUNCTION
void CsSeamGhostPrim(const TM &member, const DvceArray5D<Real> &sfs, const int m,
                     const int k, const int j, const int ia, const int ib, const bool lo,
                     const Real sf, const Real cf, const int nvars, const int nder,
                     const int ivn, const int ivt, SA &wl, SA &wr, SB &dl, SB &dr) {
  par_for_inner(member, ia, ib, [&](const int i) {
    for (int n=0; n<nvars; ++n) {
      Real q = sfs(m,n,k,j,i);
      if (n == ivn) q = (sfs(m,ivn,k,j,i) - cf*sfs(m,ivt,k,j,i))/sf;
      if (lo) {wl(n,i) = q;} else {wr(n,i) = q;}
    }
    for (int n=0; n<nder; ++n) {
      if (lo) {dl(n,i) = sfs(m,nvars+n,k,j,i);} else {dr(n,i) = sfs(m,nvars+n,k,j,i);}
    }
  });
}

}  // namespace cs_seam
#endif  // BVALS_CS_SEAM_STATES_HPP_
