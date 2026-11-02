//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_cs_seam.cpp
//! \brief ucthshare-1102: the SHARED cubed-sphere seam treatment for MHD,
//! <mhd>/cs_seam_mhd = faces | shared (read only when named).
//!
//! (1) faces: every panel-seam face is ONE Riemann problem.  The ghost-side state the
//!     sweep would build from the resampled seam ghosts is replaced by the neighbour
//!     panel's OWN PLM state of its edge cell (primitives, derived thermodynamics, bcc),
//!     the same machinery as the hydro <mesh>/cs_seam_flux = shared
//!     (bvals/cs_seam_states.hpp), extended to the radial ghost rows is-1 / ie+1 (they
//!     feed the seam EMFs on the radial boundary edges).  Both panels then solve mirror
//!     images; the existing mean of the two outward fluxes stays (conservation).
//! (2)+(3) shared (emf_scheme = uct | uct_hlld): the seam-edge and cube-vertex UCT EMFs
//!     are evaluated by every panel from the SAME inputs: its own faces on its side, the
//!     neighbour panel's OWN faces on the other (never the resampled ghosts), the shared
//!     seam faces between.  One general rule for an edge along t with faces f (unit
//!     normals n_f, L -> R), exact for any set of face directions:
//!       (v x B).t = sum_f (v.n_f) B.(t x g_f),   g_f = S^-1 n_f,  S = sum_f n_f n_f^T
//!     (S restricted to the plane normal to t), upwinded per face with its own fan,
//!       E.t = -sum_f [ a_f^L vn^L X^L + a_f^R vn^R X^R - (d_f^R X^R - d_f^L X^L) ],
//!     X^s = B^s.(t x g_f), B^s the in-plane field of the cell on side s solved from the
//!     staggered B.n of its two faces at the edge, vn^s = V.n_f with V the face velocity
//!     of the cell's OTHER face carried to the edge.  For the four faces of a regular
//!     edge (two normal pairs) this is MD21 eq. 35 with the 1/sin(p,q) of the cubed-
//!     sphere port exactly (MHD::CornerE, UctEdge); for three faces at 120 deg (a cube
//!     vertex: the gnomonic grid has the exact 3-fold symmetry there) g_f = (2/3) n_f
//!     and it is the 3-face rule E.t = -(2/3) sum_f [...] (B_w only).  At a seam radial
//!     edge the two own faces are kinked (mirror images about the seam plane): the
//!     general S^-1 takes the actual normals.
//!     Neighbour data (UctSeamPackets): per seam face position, the radial faces and
//!     along-seam faces of its first row with a, d^L, d^R, B.n carried to the seam edge,
//!     the face velocity carried to the edge (Cartesian), its normal (Cartesian), and its
//!     cells' half slopes of the Cartesian velocity (radial and along the seam).  Moved
//!     verbatim by the seam_raw face exchange (pbval_us); the along-seam reversal of the
//!     index map is undone here (lower <-> upper, the along-seam EMF sign).  A stencil
//!     that would reach past a cube vertex (the three panels' corner ghosts differ) takes
//!     the ONE-SIDED slope of the two in-panel values instead: B.n of the seam faces
//!     along a seam (values shared by both panels), each panel's own corner-cell
//!     velocity; the face velocities of the three vertex faces are used uncarried.  (A
//!     zero slope there was measured first order at the vertex: cs_test iprob=9 per-step
//!     vertex/interior L1(B) 3.0 / 6.7 / 12.8 at nx2 = 16/32/64; one-sided 2.1/2.7/3.1.)
//!     Four-block seam edges (block corners on a seam that are not a cube vertex) keep
//!     the per-block value; the edge-EMF sum/average reconciles them as before.
//!     The existing seam/vertex EMF exchange + averaging stays as a guard; with shared
//!     inputs it is a no-op to round-off (<mhd>/cs_seam_emf_avg = false skips it).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "bvals/bvals.hpp"
#include "bvals/cs_seam_states.hpp"
#include "reconstruct/plm.hpp"
#include "mhd.hpp"

namespace mhd {

namespace cs_uct {
// packet layout (per seam face position: along-seam cell c, radial cell i)
constexpr int kNU = 54;
constexpr int kA = 0;    // radial face i: a, dL, dR, Bn, V(3), Hr_up(i-1)(3), Hr_lo(i)(3), Es
constexpr int kB = 14;   // the same for radial face ie+1 (at i = ie only)
constexpr int kCL = 28;  // lower along-seam face: a, dL, dR, Bn, n(3), V(3)
constexpr int kCU = 38;  // upper along-seam face
constexpr int kHL = 48;  // along-seam half slope of the cell velocity, lower side (3)
constexpr int kHU = 51;  // upper side (3)

struct SeamCtx {
  DvceArray5D<Real> u1, u2, u3, w;
  DvceArray4D<Real> b1, b2, b3, e21, e31;
  DvceArray5D<Real> o2, o3, n2, n3;            // own / neighbour packets (x2f, x3f)
  DualArray1D<int> panel;
  DualArray1D<RegionSize> size;
  DualArray2D<NeighborBlock> nghbr;
  int is, ie, js, je, ks, ke, nx2, nx3;
};

struct UF {
  Real n[3], a, dl, dr, bn, v[3];
};

KOKKOS_INLINE_FUNCTION Real Dot3(const Real *a, const Real *b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
KOKKOS_INLINE_FUNCTION Real Half(const Real qm, const Real q0, const Real qp,
                                 const bool upper) {
  Real ql, qr;
  PLM(qm, q0, qp, ql, qr);
  return upper ? (ql - q0) : (qr - q0);
}
KOKKOS_INLINE_FUNCTION Real Edge(const Real qm, const Real q0, const Real qp,
                                 const bool upper) {
  Real ql, qr;
  PLM(qm, q0, qp, ql, qr);
  return upper ? ql : qr;
}

// local frame at (xi, eta) of panel p: rhat, unit tangents e1 (xi), e2 (eta), and the
// unit face normals nxi = (e1 - c e2)/s, neta = (e2 - c e1)/s
struct Frame {
  Real r[3], e1[3], e2[3], nxi[3], neta[3];
};
KOKKOS_INLINE_FUNCTION void MakeFrame(const int p, const Real xi, const Real eta,
                                      Frame &f) {
  cubed_sphere::PanelTangents(p, xi, eta, f.e1, f.e2);
  cubed_sphere::PanelToCart(p, xi, eta, f.r);
  const Real c = Dot3(f.e1, f.e2);
  const Real s = sqrt(1.0 - c*c);
  for (int q=0; q<3; ++q) {
    f.nxi[q] = (f.e1[q] - c*f.e2[q])/s;
    f.neta[q] = (f.e2[q] - c*f.e1[q])/s;
  }
}
// angles of cell / face positions of block m (x2 = xi, x3 = eta)
KOKKOS_INLINE_FUNCTION Real XiC(const SeamCtx &d, const int m, const int j) {
  const RegionSize &z = d.size.d_view(m);
  return 0.25*M_PI*CellCenterX(j - d.js, d.nx2, z.x2min, z.x2max);
}
KOKKOS_INLINE_FUNCTION Real XiF(const SeamCtx &d, const int m, const int j) {
  const RegionSize &z = d.size.d_view(m);
  return 0.25*M_PI*LeftEdgeX(j - d.js, d.nx2, z.x2min, z.x2max);
}
KOKKOS_INLINE_FUNCTION Real EtaC(const SeamCtx &d, const int m, const int k) {
  const RegionSize &z = d.size.d_view(m);
  return 0.25*M_PI*CellCenterX(k - d.ks, d.nx3, z.x3min, z.x3max);
}
KOKKOS_INLINE_FUNCTION Real EtaF(const SeamCtx &d, const int m, const int k) {
  const RegionSize &z = d.size.d_view(m);
  return 0.25*M_PI*LeftEdgeX(k - d.ks, d.nx3, z.x3min, z.x3max);
}
// Cartesian velocity of cell (k,j,i): w0 holds v.rhat and the CONTRAVARIANT (v^xi,
// v^eta) on the unit tangents
KOKKOS_INLINE_FUNCTION void CellV(const SeamCtx &d, const Frame &f, const int m,
                                  const int k, const int j, const int i, Real v[3]) {
  const Real a = d.w(m,IVX,k,j,i), b = d.w(m,IVY,k,j,i), c = d.w(m,IVZ,k,j,i);
  for (int q=0; q<3; ++q) v[q] = a*f.r[q] + b*f.e1[q] + c*f.e2[q];
}
// Cartesian a-weighted face velocity of a UctFaceStore slot triple (3,4,5)
//   dir 1: x1 face {rhat, e_xi, f2 = neta}; 2: x2 face {nxi, e_eta, rhat};
//   3: x3 face {neta, rhat, e_xi}
KOKKOS_INLINE_FUNCTION void FaceV(const DvceArray5D<Real> &u, const Frame &f,
                                  const int dir, const int m, const int k, const int j,
                                  const int i, Real v[3]) {
  const Real s3 = u(m,3,k,j,i), s4 = u(m,4,k,j,i), s5 = u(m,5,k,j,i);
  for (int q=0; q<3; ++q) {
    if (dir == 1) {
      v[q] = s3*f.r[q] + s4*f.e1[q] + s5*f.neta[q];
    } else if (dir == 2) {
      v[q] = s3*f.nxi[q] + s4*f.e2[q] + s5*f.r[q];
    } else {
      v[q] = s3*f.neta[q] + s4*f.r[q] + s5*f.e1[q];
    }
  }
}

// the general per-face UCT edge EMF (see the file header).  oL[f], oR[f]: the OTHER face
// of the cell on the L / R side of face f.
KOKKOS_INLINE_FUNCTION Real EdgeEmfN(const Real t[3], const int nf, const UF *f,
                                     const int *oL, const int *oR) {
  Real u1[3], u2[3];
  const Real tn = Dot3(f[0].n, t);
  for (int q=0; q<3; ++q) u1[q] = f[0].n[q] - tn*t[q];
  const Real un = sqrt(Dot3(u1, u1));
  for (int q=0; q<3; ++q) u1[q] /= un;
  u2[0] = t[1]*u1[2] - t[2]*u1[1];
  u2[1] = t[2]*u1[0] - t[0]*u1[2];
  u2[2] = t[0]*u1[1] - t[1]*u1[0];
  Real p0[4], p1[4];
  Real s00 = 0.0, s01 = 0.0, s11 = 0.0;
  for (int a=0; a<nf; ++a) {
    p0[a] = Dot3(f[a].n, u1);
    p1[a] = Dot3(f[a].n, u2);
    s00 += p0[a]*p0[a]; s01 += p0[a]*p1[a]; s11 += p1[a]*p1[a];
  }
  const Real sdet = s00*s11 - s01*s01;
  Real e = 0.0;
  for (int a=0; a<nf; ++a) {
    const Real g0 = ( s11*p0[a] - s01*p1[a])/sdet;
    const Real g1 = (-s01*p0[a] + s00*p1[a])/sdet;
    const Real tg0 = -g1, tg1 = g0;          // t x g in the (u1, u2) basis
    Real x[2], vn[2];
    for (int sd=0; sd<2; ++sd) {
      const int o = (sd == 0) ? oL[a] : oR[a];
      const Real dd = p0[a]*p1[o] - p1[a]*p0[o];
      const Real bb0 = (f[a].bn*p1[o] - p1[a]*f[o].bn)/dd;
      const Real bb1 = (p0[a]*f[o].bn - p0[o]*f[a].bn)/dd;
      x[sd] = bb0*tg0 + bb1*tg1;
      vn[sd] = Dot3(f[o].v, f[a].n);
    }
    e -= f[a].a*vn[0]*x[0] + (1.0 - f[a].a)*vn[1]*x[1] - (f[a].dr*x[1] - f[a].dl*x[0]);
  }
  return e;
}

// whether the along-seam index runs reversed between block m's seam sd and its
// neighbour panel (the seam_raw index map reverses it; the packets' lower/upper and the
// along-seam EMF sign must follow)
KOKKOS_INLINE_FUNCTION bool SeamRev(const SeamCtx &d, const int m, const int sd) {
  const int nid = (sd == 0) ? 8 : (sd == 1) ? 12 : (sd == 2) ? 24 : 28;
  const int pp = d.panel.d_view(m);
  const int qp = d.nghbr.d_view(m,nid).panel;
  const RegionSize &z = d.size.d_view(m);
  Real xi, eta;
  if (sd < 2) {
    xi = (sd == 0) ? -0.25*M_PI : 0.25*M_PI;
    eta = 0.25*M_PI*0.5*(z.x3min + z.x3max);
  } else {
    eta = (sd == 2) ? -0.25*M_PI : 0.25*M_PI;
    xi = 0.25*M_PI*0.5*(z.x2min + z.x2max);
  }
  Real e1[3], e2[3], r[3];
  cubed_sphere::PanelTangents(pp, xi, eta, e1, e2);
  cubed_sphere::PanelToCart(pp, xi, eta, r);
  const Real *sp = (sd < 2) ? e2 : e1;
  Real sv[3] = {sp[0], sp[1], sp[2]};
  Real xq, eq;
  cubed_sphere::CartToPanel(qp, r, xq, eq);
  Real f1[3], f2[3];
  cubed_sphere::PanelTangents(qp, xq, eq, f1, f2);
  const Real *sq = (fabs(xq) > fabs(eq)) ? f2 : f1;
  return (sv[0]*sq[0] + sv[1]*sq[1] + sv[2]*sq[2]) < 0.0;
}
KOKKOS_INLINE_FUNCTION bool OnPanelEdge(const SeamCtx &d, const int m, const int sd) {
  const RegionSize &z = d.size.d_view(m);
  if (sd == 0) return z.x2min <= -1.0 + 1.0e-12;
  if (sd == 1) return z.x2max >= 1.0 - 1.0e-12;
  if (sd == 2) return z.x3min <= -1.0 + 1.0e-12;
  return z.x3max >= 1.0 - 1.0e-12;
}

}  // namespace cs_uct
using namespace cs_uct;  // NOLINT

// the radial ghost rows is-1 / ie+1 of the record <-> its second set at is / ie (the seam
// exchange covers the active radial rows only)
void CsSeamMoveRadialGhosts(MeshBlockPack *pp, DvceFaceFld5D<Real> &fs, const int nfs1,
                            const bool fwd) {
  auto &indcs = pp->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pp->nmb_thispack - 1;
  auto f2 = fs.x2f;
  auto f3 = fs.x3f;
  auto &size = pp->pmb->mb_size;
  par_for("cs_seamfs_rg2", DevExeSpace(), 0, nmb1, ks, ke, 0, 3, 0, nfs1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int c, const int n) {
    const bool lo = (c & 1) == 0;
    const bool up = (c & 2) != 0;    // 0: row is-1 <-> set B at is; 1: ie+1 <-> ie
    if (lo && !(size.d_view(m).x2min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x2max >= 1.0 - 1.0e-12)) return;
    const int jf = lo ? js : je+1;
    const int ig = up ? ie+1 : is-1, ib = up ? ie : is;
    if (fwd) {f2(m,nfs1+n,k,jf,ib) = f2(m,n,k,jf,ig);}
    else {f2(m,n,k,jf,ig) = f2(m,nfs1+n,k,jf,ib);}
  });
  par_for("cs_seamfs_rg3", DevExeSpace(), 0, nmb1, js, je, 0, 3, 0, nfs1-1,
  KOKKOS_LAMBDA(const int m, const int j, const int c, const int n) {
    const bool lo = (c & 1) == 0;
    const bool up = (c & 2) != 0;
    if (lo && !(size.d_view(m).x3min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x3max >= 1.0 - 1.0e-12)) return;
    const int kf = lo ? ks : ke+1;
    const int ig = up ? ie+1 : is-1, ib = up ? ie : is;
    if (fwd) {f3(m,nfs1+n,kf,j,ib) = f3(m,n,kf,j,ig);}
    else {f3(m,n,kf,j,ig) = f3(m,nfs1+n,kf,j,ib);}
  });
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::CsSeamInit
//! \brief <mhd>/cs_seam_mhd = none | faces | shared, <mhd>/cs_seam_emf_avg (diagnostic)

void MHD::CsSeamInit(ParameterInput *pin) {
  auto *pp = pmy_pack;
  if (pin->DoesParameterExist("mhd","cs_seam_mhd")) {
    const std::string v = pin->GetString("mhd","cs_seam_mhd");
    if (v == "none") {
      cs_seam_mhd = 0;
    } else if (v == "faces") {
      cs_seam_mhd = 1;
    } else if (v == "shared") {
      cs_seam_mhd = 2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/cs_seam_mhd = '" << v << "' (none | faces | shared)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (pin->DoesParameterExist("mhd","cs_seam_emf_avg")) {
    cs_seam_noavg = !pin->GetBoolean("mhd","cs_seam_emf_avg");
    pbval_b->skip_emf_reconcile = cs_seam_noavg;
  }
  if (cs_seam_mhd == 0) return;
  if (!pp->pmesh->use_cubed_sphere || !pp->pmesh->three_d ||
      recon_method != ReconstructionMethod::plm ||
      (use_wellbalance_dynamic && use_wb_x2) || use_wellbalance_static_reconst_perturb ||
      pp->pmesh->multilevel || pp->pcoord->is_special_relativistic ||
      pp->pcoord->is_general_relativistic || use_fofc) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/cs_seam_mhd needs a 3-D cubed sphere, reconstruct = plm, no x2/x3 "
              << "well-balanced or static-perturbation reconstruction, a uniform mesh, "
              << "Newtonian MHD, no fofc" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (cs_seam_mhd == 2 && !(emf_uct || emf_uct_hlld)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/cs_seam_mhd = shared needs <mhd>/emf_scheme = uct or uct_hlld "
              << "(use faces with GS07)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int nmb = std::max((pp->nmb_thispack), (pp->pmesh->nmb_maxperrank));
  auto &ind = pp->pmesh->mb_indcs;
  const int n1 = ind.nx1 + 2*ind.ng, n2 = ind.nx2 + 2*ind.ng, n3 = ind.nx3 + 2*ind.ng;
  const int nvars = nmhd + nscalars;
  const int nfs1 = nvars + NDERIVED + 3;
  cs_seam_nfs = 2*nfs1;
  pbval_fs = new MeshBoundaryValuesCC(pp, pin, false);
  pbval_fs->seam_raw = true;
  pbval_fs->seam_xform_b = nvars + NDERIVED + 1;
  pbval_fs->seam_xform_stride = nfs1;
  pbval_fs->InitializeBuffers(cs_seam_nfs);
  pseamfs = new DvceFaceFld5D<Real>("cs_seamfs_mhd", nmb, cs_seam_nfs, n3, n2, n1);
  Kokkos::realloc(pseamfs->x1f, nmb, cs_seam_nfs, 1, 1, 1);
  if (cs_seam_mhd == 2) {
    pbval_us = new MeshBoundaryValuesCC(pp, pin, false);
    pbval_us->seam_raw = true;
    pbval_us->seam_xform_off = true;
    pbval_us->InitializeBuffers(kNU);
    pusown = new DvceFaceFld5D<Real>("cs_uctseam_own", nmb, kNU, n3, n2, n1);
    pusnbr = new DvceFaceFld5D<Real>("cs_uctseam_nbr", nmb, kNU, n3, n2, n1);
    Kokkos::realloc(pusown->x1f, nmb, kNU, 1, 1, 1);
    Kokkos::realloc(pusnbr->x1f, nmb, kNU, 1, 1, 1);
  }
  if (global_variable::my_rank == 0) {
    std::cout << "MHD: cubed-sphere seam = "
              << ((cs_seam_mhd == 2) ? "shared (faces + seam-edge/vertex UCT EMF)"
                                     : "faces (shared seam Riemann states)")
              << (cs_seam_noavg ? ", edge-EMF reconciliation OFF (diagnostic)" : "")
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::CsSeamFaceStates
//! \brief (1): fill pseamfs with this panel's own seam face states (rows is-1..ie+1), move
//! the radial ghost rows into the second record set, exchange, move them back.  After it
//! pseamfs holds, at every seam face of this rank, the NEIGHBOUR's own state.

void MHD::CsSeamFaceStates() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nvars = nmhd + nscalars;
  auto &eos_ = peos->eos_data;
  const int nder = eos_.IsGeneral() ? NDERIVED : 0;
  const int ob = nvars + NDERIVED;
  const int nfs1 = cs_seam_nfs/2;
  DvceArray2D<int> nomlat;
  cs_seam::CsSeamPLMFill(pmy_pack, w0, nvars, 0, -1, eos_.pfloor, eos_.plm_c1, *pseamfs,
                         false, nomlat, is-1, ie+1);
  cs_seam::CsSeamPLMFill(pmy_pack, wder, nder, nvars, IDPR, eos_.pfloor, eos_.plm_c1,
                         *pseamfs, false, nomlat, is-1, ie+1);
  cs_seam::CsSeamPLMFill(pmy_pack, bcc0, 3, ob, -1, eos_.pfloor, eos_.plm_c1, *pseamfs,
                         false, nomlat, is-1, ie+1);
  cs_seam::CsSeamCovariant(pmy_pack, *pseamfs, IVY, 0, is-1, ie+1);
  cs_seam::CsSeamCovariant(pmy_pack, *pseamfs, ob+1, 1, is-1, ie+1);
  CsSeamMoveRadialGhosts(pmy_pack, *pseamfs, nfs1, true);
  pbval_fs->InitFluxSeamRecv(cs_seam_nfs);
  pbval_fs->PackAndSendFluxSeamCC(*pseamfs);
  while (pbval_fs->RecvAndUnpackFluxSeamCC(*pseamfs) != TaskStatus::complete) {}
  pbval_fs->ClearFluxSend();
  pbval_fs->ClearFluxRecv();
  CsSeamMoveRadialGhosts(pmy_pack, *pseamfs, nfs1, false);
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::UctSeamPackets
//! \brief (2)+(3): this panel's own first-row data at every seam face position (pusown,
//! layout kA..kHU above), then the neighbour's (pusnbr) through the seam_raw exchange.

void MHD::UctSeamPackets() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  SeamCtx d;
  d.u1 = uct1; d.u2 = uct2; d.u3 = uct3; d.w = w0;
  d.b1 = b0.x1f; d.b2 = b0.x2f; d.b3 = b0.x3f; d.e21 = e2x1; d.e31 = e3x1;
  d.o2 = pusown->x2f; d.o3 = pusown->x3f; d.n2 = pusnbr->x2f; d.n3 = pusnbr->x3f;
  d.panel = pmy_pack->pmb->mb_panel; d.size = pmy_pack->pmb->mb_size;
  d.nghbr = pmy_pack->pmb->nghbr;
  d.is = indcs.is; d.ie = indcs.ie; d.js = indcs.js; d.je = indcs.je;
  d.ks = indcs.ks; d.ke = indcs.ke; d.nx2 = indcs.nx2; d.nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nc = indcs.nx2;   // == nx3 on the cubed sphere
  par_for("uct_seam_pack", DevExeSpace(), 0, nmb1, 0, 3, 0, nc-1, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int c, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const int p = d.panel.d_view(m);
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const int r0 = x2t ? (lo ? d.js : d.je) : (lo ? d.ks : d.ke);
    const int cb = x2t ? d.ks : d.js;          // first along-seam cell index
    // (k, j) of the cell at along-seam index a (block-relative) in cross-seam row r
    auto K = [&](const int a, const int r) { return x2t ? (cb + a) : r; };
    auto J = [&](const int a, const int r) { return x2t ? r : (cb + a); };
    auto FrC = [&](const int a, const int r, Frame &f) {
      MakeFrame(p, XiC(d, m, J(a,r)), EtaC(d, m, K(a,r)), f);
    };
    auto out = [&](const int v) -> Real & {
      return x2t ? d.o2(m, v, cb + c, lo ? d.js : d.je+1, i)
                 : d.o3(m, v, lo ? d.ks : d.ke+1, cb + c, i);
    };
    // cross-seam half slope (toward the seam) of the Cartesian velocity of cell (a, r0, ii)
    auto Hx = [&](const int a, const int ii, Real h[3]) {
      Frame fm, f0, fp;
      FrC(a, r0-1, fm); FrC(a, r0, f0); FrC(a, r0+1, fp);
      Real vm[3], v0[3], vp[3];
      CellV(d, fm, m, K(a,r0-1), J(a,r0-1), ii, vm);
      CellV(d, f0, m, K(a,r0), J(a,r0), ii, v0);
      CellV(d, fp, m, K(a,r0+1), J(a,r0+1), ii, vp);
      for (int q=0; q<3; ++q) h[q] = Half(vm[q], v0[q], vp[q], !lo);
    };
    // field value of a 4-D face array carried across the rows to the seam
    auto Bx = [&](const DvceArray4D<Real> &b, const int kf, const int jf, const int ii) {
      if (x2t) return Edge(b(m,kf,jf-1,ii), b(m,kf,jf,ii), b(m,kf,jf+1,ii), !lo);
      return Edge(b(m,kf-1,jf,ii), b(m,kf,jf,ii), b(m,kf+1,jf,ii), !lo);
    };
    Frame f0;
    FrC(c, r0, f0);
    const int k0 = K(c, r0), j0 = J(c, r0);
    // A / B: radial face(s)
    for (int s=0; s<2; ++s) {
      if (s == 1 && i != d.ie) break;
      const int ifc = (s == 0) ? i : d.ie+1, o = (s == 0) ? kA : kB;
      out(o+0) = d.u1(m,0,k0,j0,ifc);
      out(o+1) = d.u1(m,1,k0,j0,ifc);
      out(o+2) = d.u1(m,2,k0,j0,ifc);
      out(o+3) = Bx(d.b1, k0, j0, ifc);
      Real v[3], h0[3], h1[3];
      FaceV(d.u1, f0, 1, m, k0, j0, ifc, v);
      Hx(c, ifc-1, h0); Hx(c, ifc, h1);
      for (int q=0; q<3; ++q) out(o+4+q) = v[q] + 0.5*(h0[q] + h1[q]);
      Real va[3], vb[3], vc[3], vd[3];
      CellV(d, f0, m, k0, j0, ifc-2, va); CellV(d, f0, m, k0, j0, ifc-1, vb);
      CellV(d, f0, m, k0, j0, ifc, vc);   CellV(d, f0, m, k0, j0, ifc+1, vd);
      for (int q=0; q<3; ++q) {
        out(o+7+q) = Half(va[q], vb[q], vc[q], true);
        out(o+10+q) = Half(vb[q], vc[q], vd[q], false);
      }
      out(o+13) = x2t ? d.e31(m,k0,j0,ifc) : d.e21(m,k0,j0,ifc);
    }
    // C: the two along-seam faces of cell c (x2 seam: x3 faces; x3 seam: x2 faces)
    const bool vlo = OnPanelEdge(d, m, x2t ? 2 : 0) && (c == 0);
    const bool vhi = OnPanelEdge(d, m, x2t ? 3 : 1) && (c == nc-1);
    for (int s=0; s<2; ++s) {
      const int a = c + s;                 // face index along the seam (lower face of a)
      const bool vtx = (s == 0) ? vlo : vhi;
      const int o = (s == 0) ? kCL : kCU;
      const int kf = x2t ? (cb + a) : r0, jf = x2t ? r0 : (cb + a);
      const DvceArray5D<Real> &u = x2t ? d.u3 : d.u2;
      const DvceArray4D<Real> &b = x2t ? d.b3 : d.b2;
      Frame ff;
      if (x2t) {
        MakeFrame(p, XiC(d, m, jf), EtaF(d, m, kf), ff);
      } else {
        MakeFrame(p, XiF(d, m, jf), EtaC(d, m, kf), ff);
      }
      out(o+0) = u(m,0,kf,jf,i);
      out(o+1) = u(m,1,kf,jf,i);
      out(o+2) = u(m,2,kf,jf,i);
      if (!vtx) {
        out(o+3) = Bx(b, kf, jf, i);
      } else {   // one-sided along the shared seam faces of the other seam
        const Real b0 = b(m,kf,jf,i);
        const Real bi = x2t ? b(m,kf,jf + (lo ? 1 : -1),i) : b(m,kf + (lo ? 1 : -1),jf,i);
        out(o+3) = b0 + 0.5*(b0 - bi);
      }
      const Real *nn = x2t ? ff.neta : ff.nxi;
      for (int q=0; q<3; ++q) out(o+4+q) = nn[q];
      Real v[3];
      FaceV(u, ff, x2t ? 3 : 2, m, kf, jf, i, v);
      if (!vtx) {
        Real h0[3], h1[3];
        Hx(a-1, i, h0); Hx(a, i, h1);
        for (int q=0; q<3; ++q) v[q] += 0.5*(h0[q] + h1[q]);
      }
      for (int q=0; q<3; ++q) out(o+7+q) = v[q];
    }
    // H: along-seam half slopes of cell c (zero for a cube-vertex corner cell)
    if ((vlo || vhi) && nc > 1) {   // one-sided from the in-panel neighbour
      const int cn = vlo ? c+1 : c-1;
      Frame fn;
      FrC(cn, r0, fn);
      Real v0[3], vn[3];
      CellV(d, f0, m, k0, j0, i, v0);
      CellV(d, fn, m, K(cn,r0), J(cn,r0), i, vn);
      for (int q=0; q<3; ++q) {
        const Real sl = vlo ? (vn[q] - v0[q]) : (v0[q] - vn[q]);   // per cell, increasing
        out(kHL+q) = -0.5*sl;
        out(kHU+q) = 0.5*sl;
      }
    } else if (vlo || vhi) {
      for (int q=0; q<3; ++q) {out(kHL+q) = 0.0; out(kHU+q) = 0.0;}
    } else {
      Frame fm, fp;
      FrC(c-1, r0, fm); FrC(c+1, r0, fp);
      Real vm[3], v0[3], vp[3];
      CellV(d, fm, m, K(c-1,r0), J(c-1,r0), i, vm);
      CellV(d, f0, m, k0, j0, i, v0);
      CellV(d, fp, m, K(c+1,r0), J(c+1,r0), i, vp);
      for (int q=0; q<3; ++q) {
        out(kHL+q) = Half(vm[q], v0[q], vp[q], false);
        out(kHU+q) = Half(vm[q], v0[q], vp[q], true);
      }
    }
  });
  pbval_us->InitFluxSeamRecv(kNU);
  pbval_us->PackAndSendFluxSeamCC(*pusown);
  while (pbval_us->RecvAndUnpackFluxSeamCC(*pusnbr) != TaskStatus::complete) {}
  pbval_us->ClearFluxSend();
  pbval_us->ClearFluxRecv();
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::UctSeamEdges
//! \brief (2)+(3): overwrite the seam-edge and cube-vertex EMFs of efld with the shared
//! evaluation (file header).  Called at the end of the 3-D UCT branch of CornerE.

void MHD::UctSeamEdges() {
  UctSeamPackets();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  SeamCtx d;
  d.u1 = uct1; d.u2 = uct2; d.u3 = uct3; d.w = w0;
  d.b1 = b0.x1f; d.b2 = b0.x2f; d.b3 = b0.x3f; d.e21 = e2x1; d.e31 = e3x1;
  d.o2 = pusown->x2f; d.o3 = pusown->x3f; d.n2 = pusnbr->x2f; d.n3 = pusnbr->x3f;
  d.panel = pmy_pack->pmb->mb_panel; d.size = pmy_pack->pmb->mb_size;
  d.nghbr = pmy_pack->pmb->nghbr;
  d.is = indcs.is; d.ie = indcs.ie; d.js = indcs.js; d.je = indcs.je;
  d.ks = indcs.ks; d.ke = indcs.ke; d.nx2 = indcs.nx2; d.nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nc = indcs.nx2;
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  // the radial walls of CornerE's UCT branch (mean of the two wall-face EMFs)
  const bool wsf_ = wall_closed_ix1 && wall_stress_free;
  const BoundaryFlag rflag_ = BoundaryFlag::reflect;
  const BoundaryFlag uflag_ = BoundaryFlag::user;

  // (2a) along-seam edges (x2 seam: e3, x3 seam: e2), radial faces is..ie+1
  par_for("uct_seam_e_t", DevExeSpace(), 0, nmb1, 0, 3, 0, nc-1, d.is, d.ie+1,
  KOKKOS_LAMBDA(const int m, const int sd, const int c, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const int p = d.panel.d_view(m);
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const bool rev = SeamRev(d, m, sd);
    const int cb = x2t ? d.ks : d.js;
    const int kf = x2t ? (cb + c) : (lo ? d.ks : d.ke+1);   // seam face / edge indices
    const int jf = x2t ? (lo ? d.js : d.je+1) : (cb + c);
    const int ci = (i <= d.ie) ? i : d.ie;                  // packet cell
    const int ob = (i <= d.ie) ? kA : kB;                   // packet set
    auto P = [&](const DvceArray5D<Real> &a2, const DvceArray5D<Real> &a3, const int v) {
      return x2t ? a2(m, v, kf, jf, ci) : a3(m, v, kf, jf, ci);
    };
    const bool wall_in = i == d.is &&
        ((mb_bcs.d_view(m,BoundaryFace::inner_x1) == rflag_) ||
         (wsf_ && mb_bcs.d_view(m,BoundaryFace::inner_x1) == uflag_));
    const bool wall_out = i == d.ie+1 && mb_bcs.d_view(m,BoundaryFace::outer_x1) == rflag_;
    Real e;
    if (wall_in || wall_out) {
      e = 0.5*(P(d.o2, d.o3, ob+13) + (rev ? -1.0 : 1.0)*P(d.n2, d.n3, ob+13));
    } else {
      Frame fe;   // at the edge: (seam xi or eta, along-seam cell centre)
      if (x2t) {
        MakeFrame(p, XiF(d, m, jf), EtaC(d, m, kf), fe);
      } else {
        MakeFrame(p, XiC(d, m, jf), EtaF(d, m, kf), fe);
      }
      const Real *t = x2t ? fe.e2 : fe.e1;
      const Real *ns = x2t ? fe.nxi : fe.neta;
      const DvceArray5D<Real> &u = x2t ? d.u2 : d.u3;
      const DvceArray4D<Real> &b = x2t ? d.b2 : d.b3;
      UF f[4];
      // seam faces at radial cells i-1 (0) and i (2)
      for (int s=0; s<2; ++s) {
        const int ii = (s == 0) ? i-1 : i;
        UF &g = f[2*s];
        for (int q=0; q<3; ++q) g.n[q] = ns[q];
        g.a = u(m,0,kf,jf,ii); g.dl = u(m,1,kf,jf,ii); g.dr = u(m,2,kf,jf,ii);
        g.bn = (s == 0) ? Edge(b(m,kf,jf,ii-1), b(m,kf,jf,ii), b(m,kf,jf,ii+1), true)
                        : Edge(b(m,kf,jf,ii-1), b(m,kf,jf,ii), b(m,kf,jf,ii+1), false);
        FaceV(u, fe, x2t ? 2 : 3, m, kf, jf, ii, g.v);
        const int ho = (s == 0) ? ob+7 : ob+10;
        for (int q=0; q<3; ++q) {
          g.v[q] += 0.5*(P(d.o2, d.o3, ho+q) + P(d.n2, d.n3, ho+q));
        }
      }
      // radial faces: neighbour (1), own (3)
      for (int s=0; s<2; ++s) {
        UF &g = f[1 + 2*s];
        const DvceArray5D<Real> &a2 = (s == 0) ? d.n2 : d.o2;
        const DvceArray5D<Real> &a3 = (s == 0) ? d.n3 : d.o3;
        for (int q=0; q<3; ++q) g.n[q] = fe.r[q];
        g.a = P(a2, a3, ob); g.dl = P(a2, a3, ob+1); g.dr = P(a2, a3, ob+2);
        g.bn = P(a2, a3, ob+3);
        for (int q=0; q<3; ++q) g.v[q] = P(a2, a3, ob+4+q);
      }
      // cells: (0,1) nbr i-1, (1,2) nbr i, (2,3) own i, (3,0) own i-1
      int oL[4], oR[4];
      oL[0] = lo ? 1 : 3; oR[0] = lo ? 3 : 1;
      oL[2] = lo ? 1 : 3; oR[2] = lo ? 3 : 1;
      oL[1] = 0; oR[1] = 2;
      oL[3] = 0; oR[3] = 2;
      e = EdgeEmfN(t, 4, f, oL, oR);
    }
    if (x2t) {e3(m,kf,jf,i) = e;} else {e2(m,kf,jf,i) = e;}
  });

  // (2b) radial seam edges between two along-seam cells of this block (e1)
  par_for("uct_seam_e_r", DevExeSpace(), 0, nmb1, 0, 3, 1, nc-1, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int a, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const int p = d.panel.d_view(m);
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const bool rev = SeamRev(d, m, sd);
    const int cb = x2t ? d.ks : d.js;
    const int kf = x2t ? (cb + a) : (lo ? d.ks : d.ke+1);
    const int jf = x2t ? (lo ? d.js : d.je+1) : (cb + a);
    Frame fe;   // the edge point (seam, along-seam face a)
    if (x2t) {
      MakeFrame(p, XiF(d, m, jf), EtaF(d, m, kf), fe);
    } else {
      MakeFrame(p, XiF(d, m, jf), EtaF(d, m, kf), fe);
    }
    const Real *t = fe.r;
    const Real *sp = x2t ? fe.e2 : fe.e1;
    const DvceArray5D<Real> &u = x2t ? d.u2 : d.u3;
    const DvceArray4D<Real> &b = x2t ? d.b2 : d.b3;
    const bool vlo = OnPanelEdge(d, m, x2t ? 2 : 0);
    const bool vhi = OnPanelEdge(d, m, x2t ? 3 : 1);
    // packet of along-seam cell cc at this seam
    auto P = [&](const DvceArray5D<Real> &a2, const DvceArray5D<Real> &a3, const int cc,
                 const int v) {
      return x2t ? a2(m, v, cb + cc, jf, i) : a3(m, v, kf, cb + cc, i);
    };
    // seam face of cell cc (along-seam index), as (k, j)
    auto SK = [&](const int cc) { return x2t ? (cb + cc) : kf; };
    auto SJ = [&](const int cc) { return x2t ? jf : (cb + cc); };
    auto Bs = [&](const int cc) {   // seam-face field of cell cc
      return x2t ? b(m,cb+cc,jf,i) : b(m,kf,cb+cc,i);
    };
    UF f[4];
    // S1 = seam face of cell a-1 (0), S2 = of cell a (2)
    for (int s=0; s<2; ++s) {
      const int cc = (s == 0) ? a-1 : a;
      UF &g = f[2*s];
      Frame ff;
      if (x2t) {
        MakeFrame(p, XiF(d, m, jf), EtaC(d, m, SK(cc)), ff);
      } else {
        MakeFrame(p, XiC(d, m, SJ(cc)), EtaF(d, m, kf), ff);
      }
      const Real *ns = x2t ? ff.nxi : ff.neta;
      for (int q=0; q<3; ++q) g.n[q] = ns[q];
      g.a = u(m,0,SK(cc),SJ(cc),i); g.dl = u(m,1,SK(cc),SJ(cc),i);
      g.dr = u(m,2,SK(cc),SJ(cc),i);
      const bool corner = (cc == 0 && vlo) || (cc == nc-1 && vhi);
      if (!corner) {
        g.bn = Edge(Bs(cc-1), Bs(cc), Bs(cc+1), s == 0);
      } else if (nc > 1) {
        const Real sl = (cc == 0 && vlo) ? (Bs(cc+1) - Bs(cc)) : (Bs(cc) - Bs(cc-1));
        g.bn = Bs(cc) + ((s == 0) ? 0.5 : -0.5)*sl;
      } else {
        g.bn = Bs(cc);
      }
      FaceV(u, ff, x2t ? 2 : 3, m, SK(cc), SJ(cc), i, g.v);
      // toward the edge: the upper half of cell a-1, the lower half of cell a
      const int ho = (s == 0) ? kHU : kHL, hn = ((s == 0) != rev) ? kHU : kHL;
      for (int q=0; q<3; ++q) {
        g.v[q] += 0.5*(P(d.o2, d.o3, cc, ho+q) + P(d.n2, d.n3, cc, hn+q));
      }
    }
    // E = own along-seam face (3): the lower face of cell a; W = neighbour's (1)
    for (int s=0; s<2; ++s) {
      UF &g = f[(s == 0) ? 3 : 1];
      const DvceArray5D<Real> &a2 = (s == 0) ? d.o2 : d.n2;
      const DvceArray5D<Real> &a3 = (s == 0) ? d.o3 : d.n3;
      const int o = (s == 0 || !rev) ? kCL : kCU;
      g.a = P(a2, a3, a, o); g.dl = P(a2, a3, a, o+1); g.dr = P(a2, a3, a, o+2);
      g.bn = P(a2, a3, a, o+3);
      for (int q=0; q<3; ++q) {
        g.n[q] = P(a2, a3, a, o+4+q);
        g.v[q] = P(a2, a3, a, o+7+q);
      }
    }
    // cells: (0,1) nbr a-1, (1,2) nbr a, (2,3) own a, (3,0) own a-1
    int oL[4], oR[4];
    oL[0] = lo ? 1 : 3; oR[0] = lo ? 3 : 1;
    oL[2] = lo ? 1 : 3; oR[2] = lo ? 3 : 1;
    oL[3] = 0; oR[3] = 2;
    if (Dot3(f[1].n, sp) > 0.0) {oL[1] = 0; oR[1] = 2;} else {oL[1] = 2; oR[1] = 0;}
    e1(m,kf,jf,i) = EdgeEmfN(t, 4, f, oL, oR);
  });

  // (3) cube-vertex radial edges: the three seam faces at 120 deg
  par_for("uct_seam_vtx", DevExeSpace(), 0, nmb1, 0, 3, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int cn, const int i) {
    const bool jlo = (cn & 1) == 0, klo = (cn & 2) == 0;
    const int sj = jlo ? 0 : 1, sk = klo ? 2 : 3;
    if (!OnPanelEdge(d, m, sj) || !OnPanelEdge(d, m, sk)) return;
    const int p = d.panel.d_view(m);
    const int jc = jlo ? d.js : d.je+1, kc = klo ? d.ks : d.ke+1;
    const int jv = jlo ? d.js : d.je, kv = klo ? d.ks : d.ke;
    Frame fv, fa, fb;
    MakeFrame(p, XiF(d, m, jc), EtaF(d, m, kc), fv);
    MakeFrame(p, XiF(d, m, jc), EtaC(d, m, kv), fa);
    MakeFrame(p, XiC(d, m, jv), EtaF(d, m, kc), fb);
    UF f[3];
    // a = the x2 seam face (0), b = the x3 seam face (2), c = the third (1)
    for (int q=0; q<3; ++q) {f[0].n[q] = fa.nxi[q]; f[2].n[q] = fb.neta[q];}
    f[0].a = d.u2(m,0,kv,jc,i); f[0].dl = d.u2(m,1,kv,jc,i); f[0].dr = d.u2(m,2,kv,jc,i);
    f[0].bn = d.b2(m,kv,jc,i);
    if (d.nx2 > 1) {f[0].bn += 0.5*(d.b2(m,kv,jc,i) - d.b2(m,klo ? kv+1 : kv-1,jc,i));}
    FaceV(d.u2, fa, 2, m, kv, jc, i, f[0].v);
    f[2].a = d.u3(m,0,kc,jv,i); f[2].dl = d.u3(m,1,kc,jv,i); f[2].dr = d.u3(m,2,kc,jv,i);
    f[2].bn = d.b3(m,kc,jv,i);
    if (d.nx2 > 1) {f[2].bn += 0.5*(d.b3(m,kc,jv,i) - d.b3(m,kc,jlo ? jv+1 : jv-1,i));}
    FaceV(d.u3, fb, 3, m, kc, jv, i, f[2].v);
    // e_a: along face a away from the vertex
    Real ea[3];
    for (int q=0; q<3; ++q) ea[q] = klo ? fa.e2[q] : -fa.e2[q];
    Real e = 0.0;
    for (int s=0; s<2; ++s) {
      // the third face from the x2-seam neighbour (s = 0) or the x3-seam one (s = 1)
      const int sd = (s == 0) ? sj : sk;
      const bool rev = SeamRev(d, m, sd);
      const bool want_lo = (s == 0) ? klo : jlo;
      const int o = (want_lo != rev) ? kCL : kCU;
      for (int v=0; v<10; ++v) {
        const Real x = (s == 0) ? d.n2(m, o+v, kv, jc, i) : d.n3(m, o+v, kc, jv, i);
        if (v == 0) f[1].a = x;
        else if (v == 1) f[1].dl = x;
        else if (v == 2) f[1].dr = x;
        else if (v == 3) f[1].bn = x;
        else if (v < 7) f[1].n[v-4] = x;
        else f[1].v[v-7] = x;
      }
      int oL[3], oR[3];
      // cells: P = (a,b), Q = (a,c), R = (b,c)
      oL[0] = jlo ? 1 : 2; oR[0] = jlo ? 2 : 1;
      oL[2] = klo ? 1 : 0; oR[2] = klo ? 0 : 1;
      if (Dot3(ea, f[1].n) < 0.0) {oL[1] = 0; oR[1] = 2;} else {oL[1] = 2; oR[1] = 0;}
      e += 0.5*EdgeEmfN(fv.r, 3, f, oL, oR);
    }
    e1(m,kc,jc,i) = e;
  });
}

}  // namespace mhd
