//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file flux_seam_cc.cpp
//! \brief Makes the flux of cell-centred variables through a cubed-sphere panel seam
//! single-valued, so the seam is exactly conservative. The seam ghosts are interpolated,
//! so the two panels compute slightly different fluxes through the shared face. Each
//! panel sends its outward flux sgn*flx through its copy of the face, and both replace
//! their own by F_out <- 0.5*(F_out - F_out_neighbour). The tangential momentum fluxes
//! are transformed between the panel bases (cubed_sphere::TransformMomentum); scalar and
//! radial fluxes are copied. Seam faces pair one-to-one, so no resample is needed. The
//! flux buffers are shared with the fine/coarse flux correction, which never serves the
//! same (MeshBlock, neighbour).

#include <cstdlib>
#include <map>
#include <vector>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn DvceArray1D<int> CrossPanelPairs
//! \brief (block, buffer) pairs whose neighbour is on another panel, as m*nnghbr + n,
//! cached per boundary object (the cubed-sphere neighbour table is static). The seam
//! kernels launch one team per (pair, variable) over this list.

namespace {
DvceArray1D<int> CrossPanelPairs(MeshBlockPack *pp, const void *key) {
  static std::map<const void*, DvceArray1D<int>> cache;
  // the Views must be freed before Kokkos::finalize
  static const bool hook_ = [] {
    Kokkos::push_finalize_hook([] {cache.clear();});
    return true;
  }();
  (void) hook_;
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;
  const int nmb = pp->nmb_thispack, nnghbr = pp->pmb->nnghbr;
  auto &nghbr = pp->pmb->nghbr;
  auto &mbpanel = pp->pmb->mb_panel;
  std::vector<int> l;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid >= 0 && nghbr.h_view(m,n).panel != mbpanel.h_view(m)) {
        l.push_back(m*nnghbr + n);
      }
    }
  }
  const int np = static_cast<int>(l.size());
  DvceArray1D<int> d("cs_xpanel_pairs", (np > 0) ? np : 1);
  auto h = Kokkos::create_mirror_view(d);
  for (int q=0; q<np; ++q) {
    h(q) = l[q];
  }
  if (np == 0) {
    h(0) = -1;
  }
  Kokkos::deep_copy(d, h);
  cache[key] = d;
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn DvceArray1D<int> SamePanelFacePairs
//! \brief (block, buffer) pairs of a seam_same_panel exchange: every same-panel,
//! same-level, non-polar neighbour. Cached per boundary object.
DvceArray1D<int> SamePanelFacePairs(MeshBlockPack *pp, const void *key) {
  static std::map<const void*, DvceArray1D<int>> cache;
  static const bool hook_ = [] {
    Kokkos::push_finalize_hook([] {cache.clear();});
    return true;
  }();
  (void) hook_;
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;
  const int nmb = pp->nmb_thispack, nnghbr = pp->pmb->nnghbr;
  auto &nghbr = pp->pmb->nghbr;
  auto &mbpanel = pp->pmb->mb_panel;
  auto &mblev = pp->pmb->mb_lev;
  std::vector<int> l;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      const auto nb = nghbr.h_view(m,n);
      if (nb.gid >= 0 && nb.panel == mbpanel.h_view(m) && nb.lev == mblev.h_view(m) &&
          nb.polar <= 0) {
        l.push_back(m*nnghbr + n);
      }
    }
  }
  const int np = static_cast<int>(l.size());
  DvceArray1D<int> d("cs_spanel_pairs", (np > 0) ? np : 1);
  auto h = Kokkos::create_mirror_view(d);
  for (int q=0; q<np; ++q) {
    h(q) = l[q];
  }
  if (np == 0) {
    h(0) = -1;
  }
  Kokkos::deep_copy(d, h);
  cache[key] = d;
  return d;
}

//! \brief Host check: is (m,n) a neighbour served by this exchange (cross-panel seam, or
//! same-panel same-level non-polar neighbour with seam_same_panel)?
inline bool SeamPairWantedH(MeshBlockPack *pp, const int m, const int n, const bool ssp) {
  const auto nb = pp->pmb->nghbr.h_view(m,n);
  const int myp = pp->pmb->mb_panel.h_view(m);
  if (!ssp) {
    return nb.panel != myp;
  }
  return nb.panel == myp && nb.lev == pp->pmb->mb_lev.h_view(m) && nb.polar <= 0;
}
}  // namespace

namespace {
//! \fn bool SeamFaceGeom
//! \brief Classifies a boundary buffer from its index range: true for an x2 or x3 face
//! buffer (x2face tells which), with sgn the sign turning the stored +x2/+x3 flux into
//! the outward flux. x1 faces, edges and corners return false.
KOKKOS_INLINE_FUNCTION
bool SeamFaceGeom(const MeshBufferIndcs &ix, const int is, const int ie,
                  const int js, const int ks, bool &x2face, Real &sgn) {
  // a face buffer spans the full active range in x1 (an x1x2 edge buffer is also thin in
  // j and full in k)
  if (ix.bis != is || ix.bie != ie) {
    return false;
  }
  const bool jthin = (ix.bjs == ix.bje);
  const bool kthin = (ix.bks == ix.bke);
  // x1 face, or an x2x3 edge
  if (jthin == kthin) {
    return false;
  }
  x2face = jthin;
  sgn = (x2face ? (ix.bjs == js) : (ix.bks == ks)) ? -1.0 : 1.0;
  return true;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn static int SeamFluxVl
//! \brief Vector length of the seam-flux team kernels (env BVALS_SEAM_VL; default 32 on
//! CUDA, 1 otherwise).
static int SeamFluxVl() {
  static int v = -1;
  if (v < 0) {
    const char *e = std::getenv("BVALS_SEAM_VL");
#if defined(KOKKOS_ENABLE_CUDA)
    v = (e != nullptr) ? std::atoi(e) : 32;
#else
    v = (e != nullptr) ? std::atoi(e) : 1;
#endif
    if (v < 1) {
      v = 1;
    }
  }
  return v;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::InitFluxSeamRecv
//! \brief Post non-blocking receives for the seam flux exchange.

TaskStatus MeshBoundaryValuesCC::InitFluxSeamRecv(const int nvars) {
#if MPI_PARALLEL_ENABLED
  // wait for queued kernels that still read the receive buffers before posting receives
  // (GPU-aware MPI writes the device buffers directly)
  if (global_variable::nranks > 1) {
    Kokkos::fence();
  }
  int &nmb = pmy_pack->nmb_thispack;
  int &nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const int is_ = pmy_pack->pmesh->mb_indcs.is;
  const int ie_ = pmy_pack->pmesh->mb_indcs.ie;
  const int js_ = pmy_pack->pmesh->mb_indcs.js;
  const int ks_ = pmy_pack->pmesh->mb_indcs.ks;

  bool no_errors=true;
  bool x2f; Real sg;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if ((nghbr.h_view(m,n).gid >= 0) &&
          SeamFaceGeom(recvbuf[n].iflux_same[0], is_, ie_, js_, ks_, x2f, sg) &&
          SeamPairWantedH(pmy_pack, m, n, seam_same_panel)) {
        int drank = nghbr.h_view(m,n).rank;
        if (drank != global_variable::my_rank) {
          // tag uses the local ID and buffer index of the *receiving* MeshBlock
          int tag = CreateBvals_MPI_Tag(m, n);
          int data_size = nvars*(recvbuf[n].iflxs_ndat);
          auto recv_ptr = Kokkos::subview(recvbuf[n].flux, m, Kokkos::ALL);
          int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                               comm_flux, &(recvbuf[n].flux_req[m]));
          if (ierr != MPI_SUCCESS) {
            no_errors=false;
          }
        }
      }
    }
  }
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in posting non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::PackAndSendFluxSeamCC
//! \brief Packs the outward flux through each panel-seam face and sends it to the panel
//! on the other side, using the signed index map of the halo exchange (buffer strides
//! transposed across a swap seam).

TaskStatus MeshBoundaryValuesCC::PackAndSendFluxSeamCC(DvceFaceFld5D<Real> &flx) {
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  int nvar = flx.x1f.extent_int(1);

  int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const bool sraw = seam_raw;   // (bvals.hpp)
  const int sxs_ = seam_xform_stride, sxb_ = seam_xform_b;
  const bool sxoff_ = seam_xform_off;
  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;
  auto rbuf = RecvBufDv();
  // only variable sfv_ (see bvals.hpp); averaging modes only
  const int sfm_ = pmy_pack->pmesh->cs_seam_flux;
  const int sfv_ = (seam_flux_var >= 0 && !sraw) ? seam_flux_var : -1;

  const bool ssp_ = seam_same_panel;
  auto xpp_ = ssp_ ? SamePanelFacePairs(pmy_pack, this)
                   : CrossPanelPairs(pmy_pack, this);
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (xpp_.extent_int(0)*nvar), Kokkos::AUTO,
                              SeamFluxVl());
  BvalsTeamFor("SeamFluxSendLocal", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int mn_ = xpp_(tmember.league_rank()/nvar);
    if (mn_ < 0) {
      return;
    }
    const int m = mn_/nnghbr;
    const int n = mn_ - m*nnghbr;
    const int v = tmember.league_rank() - (tmember.league_rank()/nvar)*nvar;
    if (sfv_ >= 0 && v != sfv_) {
      return;
    }

    if (nghbr.d_view(m,n).gid < 0) {
      return;
    }
    bool x2face; Real sgn;
    if (!SeamFaceGeom(rbuf[n].iflux_same[0], cs_indcs.is, cs_indcs.ie,
                      cs_indcs.js, cs_indcs.ks, x2face, sgn)) {
      return;
    }
    const int my_panel = mbpanel.d_view(m);
    const int dst_panel = nghbr.d_view(m,n).panel;
    // seam_same_panel
    if ((dst_panel == my_panel) != ssp_) {
      return;
    }
    // off-rank neighbours: second kernel
    if (nghbr.d_view(m,n).rank != my_rank) {
      return;
    }

    const int il = rbuf[n].iflux_same[0].bis;
    const int iu = rbuf[n].iflux_same[0].bie;
    const int jl = rbuf[n].iflux_same[0].bjs;
    const int ju = rbuf[n].iflux_same[0].bje;
    const int kl = rbuf[n].iflux_same[0].bks;
    const int ku = rbuf[n].iflux_same[0].bke;
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;
    const int nkj = nk*nj;

    const int dm = nghbr.d_view(m,n).gid - mbgid.d_view(0);
    const int dn = nghbr.d_view(m,n).dest;

    // the seam index map, exactly as in MeshBoundaryValuesCC::PackAndSendCC
    int aj = 1, bj = 0, ak = 1, bk = 0, sj = 1, sk = nj;
    if (!ssp_) {   // one chart (seam_same_panel) = the identity map
      const PanelBoundaries pb = GetPanelBoundary(my_panel, dst_panel);
      const int rev_a_preswap = (pb.swap_ax == 1) ? pb.rev_b : pb.rev_a;
      const int rev_b_preswap = (pb.swap_ax == 1) ? pb.rev_a : pb.rev_b;
      if (pb.swap_ax == 1) {
        sj = nk;
        sk = 1;
      }
      if (rev_a_preswap) {
        aj = -1;
        bj = jl + ju;
      }
      if (rev_b_preswap) {
        ak = -1;
        bk = kl + ku;
      }
    }

      // mass, energy and radial momentum fluxes are scalars across the seam; only the
      // tangential momentum pair is transformed (pairs given by seam_xform_*)
    const int vv = (sxs_ > 0) ? (v % sxs_) : v;
    const int vb = v - vv;
    const bool vpair = (vv == IM2) || (vv == IM3);
    const bool bpair = (sxb_ >= 0) && ((vv == sxb_) || (vv == sxb_ + 1));
    const bool xform = !sxoff_ && !ssp_ && (vpair || bpair);
    const int p0 = vb + (bpair ? sxb_ : IM2);
    const bool first_ = bpair ? (vv == sxb_) : (vv == IM2);

    Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
      int k = idx / nj;
      int j = (idx - k * nj) + jl;
      k += kl;
      const int jj = aj*j + bj;
      const int kk = ak*k + bk;
      // seam geometry depends on (xi,eta) only: computed once per (k,j), outside the
      // radial loop
      cubed_sphere::SeamXform xf;
      if (xform) {
        // (xi,eta) of the FACE CENTRE: the seam-normal coordinate is a cell EDGE.
        const int js_ = cs_indcs.js, ks_ = cs_indcs.ks;
        const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
        const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;
        const Real xi = 0.25*M_PI*(x2face
            ? LeftEdgeX(jj-js_, cs_indcs.nx2, x2mn, x2mx)
            : CellCenterX(jj-js_, cs_indcs.nx2, x2mn, x2mx));
        const Real eta = 0.25*M_PI*(x2face
            ? CellCenterX(kk-ks_, cs_indcs.nx3, x3mn, x3mx)
            : LeftEdgeX(kk-ks_, cs_indcs.nx3, x3mn, x3mx));
        cubed_sphere::SeamXformAt(my_panel, dst_panel, xi, eta, xf);
      }
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1), [&](const int i) {
        Real val;
        if (!xform) {
          val = x2face ? flx.x2f(m,v,kk,jj,i) : flx.x3f(m,v,kk,jj,i);
        } else {
          const Real f2 = x2face ? flx.x2f(m,p0,kk,jj,i) : flx.x3f(m,p0,kk,jj,i);
          const Real f3 = x2face ? flx.x2f(m,p0+1,kk,jj,i) : flx.x3f(m,p0+1,kk,jj,i);
          Real f2o, f3o;
          cubed_sphere::ApplyMomentumXform(xf, f2, f3, f2o, f3o);
          val = first_ ? f2o : f3o;
        }
        const int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*v);
        rbuf[dn].flux(dm, index) = sraw ? val : sgn*val;
      });
    });
  });

#if MPI_PARALLEL_ENABLED
  // off-rank part as a separate kernel, so no closure captures both buffer arrays (too
  // large for the kernel argument space)
  auto sbuf = SendBufDv();
  BvalsTeamFor("SeamFluxSendMPI", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int mn_ = xpp_(tmember.league_rank()/nvar);
    if (mn_ < 0) {
      return;
    }
    const int m = mn_/nnghbr;
    const int n = mn_ - m*nnghbr;
    const int v = tmember.league_rank() - (tmember.league_rank()/nvar)*nvar;
    if (sfv_ >= 0 && v != sfv_) {
      return;
    }

    if (nghbr.d_view(m,n).gid < 0) {
      return;
    }
    bool x2face; Real sgn;
    if (!SeamFaceGeom(sbuf[n].iflux_same[0], cs_indcs.is, cs_indcs.ie,
                      cs_indcs.js, cs_indcs.ks, x2face, sgn)) {
      return;
    }
    const int my_panel = mbpanel.d_view(m);
    const int dst_panel = nghbr.d_view(m,n).panel;
    // seam_same_panel
    if ((dst_panel == my_panel) != ssp_) {
      return;
    }
    // on-rank neighbours: first kernel
    if (nghbr.d_view(m,n).rank == my_rank) {
      return;
    }

    const int il = sbuf[n].iflux_same[0].bis;
    const int iu = sbuf[n].iflux_same[0].bie;
    const int jl = sbuf[n].iflux_same[0].bjs;
    const int ju = sbuf[n].iflux_same[0].bje;
    const int kl = sbuf[n].iflux_same[0].bks;
    const int ku = sbuf[n].iflux_same[0].bke;
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;
    const int nkj = nk*nj;


    // the seam index map, exactly as in MeshBoundaryValuesCC::PackAndSendCC
    int aj = 1, bj = 0, ak = 1, bk = 0, sj = 1, sk = nj;
    if (!ssp_) {   // one chart (seam_same_panel) = the identity map
      const PanelBoundaries pb = GetPanelBoundary(my_panel, dst_panel);
      const int rev_a_preswap = (pb.swap_ax == 1) ? pb.rev_b : pb.rev_a;
      const int rev_b_preswap = (pb.swap_ax == 1) ? pb.rev_a : pb.rev_b;
      if (pb.swap_ax == 1) {
        sj = nk;
        sk = 1;
      }
      if (rev_a_preswap) {
        aj = -1;
        bj = jl + ju;
      }
      if (rev_b_preswap) {
        ak = -1;
        bk = kl + ku;
      }
    }

      // mass, energy and radial momentum fluxes are scalars across the seam; only the
      // tangential momentum pair is transformed (pairs given by seam_xform_*)
    const int vv = (sxs_ > 0) ? (v % sxs_) : v;
    const int vb = v - vv;
    const bool vpair = (vv == IM2) || (vv == IM3);
    const bool bpair = (sxb_ >= 0) && ((vv == sxb_) || (vv == sxb_ + 1));
    const bool xform = !sxoff_ && !ssp_ && (vpair || bpair);
    const int p0 = vb + (bpair ? sxb_ : IM2);
    const bool first_ = bpair ? (vv == sxb_) : (vv == IM2);

    Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
      int k = idx / nj;
      int j = (idx - k * nj) + jl;
      k += kl;
      const int jj = aj*j + bj;
      const int kk = ak*k + bk;
      // seam geometry depends on (xi,eta) only: computed once per (k,j), outside the
      // radial loop
      cubed_sphere::SeamXform xf;
      if (xform) {
        // (xi,eta) of the FACE CENTRE: the seam-normal coordinate is a cell EDGE.
        const int js_ = cs_indcs.js, ks_ = cs_indcs.ks;
        const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
        const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;
        const Real xi = 0.25*M_PI*(x2face
            ? LeftEdgeX(jj-js_, cs_indcs.nx2, x2mn, x2mx)
            : CellCenterX(jj-js_, cs_indcs.nx2, x2mn, x2mx));
        const Real eta = 0.25*M_PI*(x2face
            ? CellCenterX(kk-ks_, cs_indcs.nx3, x3mn, x3mx)
            : LeftEdgeX(kk-ks_, cs_indcs.nx3, x3mn, x3mx));
        cubed_sphere::SeamXformAt(my_panel, dst_panel, xi, eta, xf);
      }
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1), [&](const int i) {
        Real val;
        if (!xform) {
          val = x2face ? flx.x2f(m,v,kk,jj,i) : flx.x3f(m,v,kk,jj,i);
        } else {
          const Real f2 = x2face ? flx.x2f(m,p0,kk,jj,i) : flx.x3f(m,p0,kk,jj,i);
          const Real f3 = x2face ? flx.x2f(m,p0+1,kk,jj,i) : flx.x3f(m,p0+1,kk,jj,i);
          Real f2o, f3o;
          cubed_sphere::ApplyMomentumXform(xf, f2, f3, f2o, f3o);
          val = first_ ? f2o : f3o;
        }
        const int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*v);
        sbuf[n].flux(m, index) = sraw ? val : sgn*val;
      });
    });
  });

#endif

#if MPI_PARALLEL_ENABLED
  // fence before the first MPI_Isend only (on-rank exchanges share the device stream)
  bool fenced_ = false;
  bool no_errors=true;
  auto &nghbr_h = pmy_pack->pmb->nghbr;
  auto &mbpanel_h = pmy_pack->pmb->mb_panel;
  const int is_h = pmy_pack->pmesh->mb_indcs.is;
  const int ie_h = pmy_pack->pmesh->mb_indcs.ie;
  const int js_h = pmy_pack->pmesh->mb_indcs.js;
  const int ks_h = pmy_pack->pmesh->mb_indcs.ks;
  bool x2f_h; Real sg_h;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if ((nghbr_h.h_view(m,n).gid >= 0) &&
          SeamFaceGeom(sendbuf[n].iflux_same[0], is_h, ie_h, js_h, ks_h, x2f_h, sg_h) &&
          SeamPairWantedH(pmy_pack, m, n, seam_same_panel)) {
        int dn = nghbr_h.h_view(m,n).dest;
        int drank = nghbr_h.h_view(m,n).rank;
        if (drank != my_rank) {
          if (!fenced_) {
            Kokkos::fence();
            fenced_ = true;
          }
          int lid = nghbr_h.h_view(m,n).gid - pmy_pack->pmesh->gids_eachrank[drank];
          int tag = CreateBvals_MPI_Tag(lid, dn);
          int data_size = nvar*(sendbuf[n].iflxs_ndat);
          auto send_ptr = Kokkos::subview(sendbuf[n].flux, m, Kokkos::ALL);
          int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                               comm_flux, &(sendbuf[n].flux_req[m]));
          if (ierr != MPI_SUCCESS) {
            no_errors=false;
          }
        }
      }
    }
  }
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in posting sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::RecvAndUnpackFluxSeamCC
//! \brief Replaces the flux on each seam face by F_out <- 0.5*(F_out - F_out_neighbour),
//! so the outward fluxes of the two panels are exact negatives (with seam_raw, copies the
//! neighbour's face record instead).

TaskStatus MeshBoundaryValuesCC::RecvAndUnpackFluxSeamCC(DvceFaceFld5D<Real> &flx) {
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto rbuf = RecvBufDv();
  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;

#if MPI_PARALLEL_ENABLED
  bool bflag = false;
  bool no_errors=true;
  bool x2f_h; Real sg_h;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if ((nghbr.h_view(m,n).gid >= 0) &&
          SeamFaceGeom(recvbuf[n].iflux_same[0], cs_indcs.is, cs_indcs.ie,
                       cs_indcs.js, cs_indcs.ks, x2f_h, sg_h) &&
          SeamPairWantedH(pmy_pack, m, n, seam_same_panel) &&
          (nghbr.h_view(m,n).rank != global_variable::my_rank)) {
        int test;
        int ierr = MPI_Test(&(recvbuf[n].flux_req[m]), &test, MPI_STATUS_IGNORE);
        if (ierr != MPI_SUCCESS) {
          no_errors=false;
        }
        if (!(static_cast<bool>(test))) {
          bflag = true;
        }
      }
    }
  }
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in testing non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (bflag) {
    return TaskStatus::incomplete;
  }
#endif

  int nvar = flx.x1f.extent_int(1);
  const int sfmode = seam_raw ? -1 : pmy_pack->pmesh->cs_seam_flux;
  if (sfmode == -1) {
    // seam_raw: every variable of the face takes the neighbour's value (the pack already
    // transformed the tangential pair); one team per (MeshBlock, buffer, variable)
    const bool ssp_ = seam_same_panel;
    auto xpp_ = ssp_ ? SamePanelFacePairs(pmy_pack, this)
                     : CrossPanelPairs(pmy_pack, this);
    Kokkos::TeamPolicy<> pol_r(DevExeSpace(), (xpp_.extent_int(0)*nvar), Kokkos::AUTO);
    BvalsTeamFor("SeamStateRecvRaw", pol_r, KOKKOS_LAMBDA(TeamMember_t tmember) {
      const int mn_ = xpp_(tmember.league_rank()/nvar);
      if (mn_ < 0) {
        return;
      }
      const int m = mn_/nnghbr;
      const int n = mn_ - m*nnghbr;
      const int v = tmember.league_rank() - (tmember.league_rank()/nvar)*nvar;
      if (nghbr.d_view(m,n).gid < 0) {
        return;
      }
      bool x2face; Real sgn;
      if (!SeamFaceGeom(rbuf[n].iflux_same[0], cs_indcs.is, cs_indcs.ie,
                        cs_indcs.js, cs_indcs.ks, x2face, sgn)) {
        return;
      }
      if ((nghbr.d_view(m,n).panel == mbpanel.d_view(m)) != ssp_) {
        return;
      }
      const int il = rbuf[n].iflux_same[0].bis;
      const int iu = rbuf[n].iflux_same[0].bie;
      const int jl = rbuf[n].iflux_same[0].bjs;
      const int ju = rbuf[n].iflux_same[0].bje;
      const int kl = rbuf[n].iflux_same[0].bks;
      const int ku = rbuf[n].iflux_same[0].bke;
      const int ni = iu - il + 1;
      const int nj = ju - jl + 1;
      const int nk = ku - kl + 1;
      const int nvs = ni*nj*nk;
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nvs), [&](const int idx) {
        const int k = idx/(ni*nj);
        const int j = (idx - k*ni*nj)/ni;
        const int i = idx - k*ni*nj - j*ni;
        const Real recv = rbuf[n].flux(m, idx + nvs*v);
        if (x2face) {
          flx.x2f(m,v,k+kl,j+jl,i+il) = recv;
        } else {
          flx.x3f(m,v,k+kl,j+jl,i+il) = recv;
        }
      });
    });
    return TaskStatus::complete;
  }
  // only variable sfv_ (see bvals.hpp; this is an averaging mode)
  const int sfv_ = seam_flux_var;
  auto xpp_ = CrossPanelPairs(pmy_pack, this);
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (xpp_.extent_int(0)*nvar), Kokkos::AUTO,
                              SeamFluxVl());
  BvalsTeamFor("SeamFluxRecv", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int mn_ = xpp_(tmember.league_rank()/nvar);
    if (mn_ < 0) {
      return;
    }
    const int m = mn_/nnghbr;
    const int n = mn_ - m*nnghbr;
    const int v = tmember.league_rank() - (tmember.league_rank()/nvar)*nvar;
    if (sfv_ >= 0 && v != sfv_) {
      return;
    }

    if (nghbr.d_view(m,n).gid < 0) {
      return;
    }
    bool x2face; Real sgn;
    if (!SeamFaceGeom(rbuf[n].iflux_same[0], cs_indcs.is, cs_indcs.ie,
                      cs_indcs.js, cs_indcs.ks, x2face, sgn)) {
      return;
    }
    if (nghbr.d_view(m,n).panel == mbpanel.d_view(m)) {
      return;
    }

    const int il = rbuf[n].iflux_same[0].bis;
    const int iu = rbuf[n].iflux_same[0].bie;
    const int jl = rbuf[n].iflux_same[0].bjs;
    const int ju = rbuf[n].iflux_same[0].bje;
    const int kl = rbuf[n].iflux_same[0].bks;
    const int ku = rbuf[n].iflux_same[0].bke;
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;
    const int nkj = nk*nj;

    Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
      int k = idx / nj;
      int j = (idx - k * nj) + jl;
      k += kl;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1), [&](const int i) {
        const int index = i-il + ni*((j-jl) + nj*(k-kl) + nk*nj*v);
        // rbuf holds the neighbour's OUTWARD flux; ours is sgn*flx. Averaging the two
        // outward values and writing back gives 0.5*flx - 0.5*sgn*recv.
        const Real recv = rbuf[n].flux(m, index);
        if (x2face) {
          flx.x2f(m,v,k,j,i) = 0.5*flx.x2f(m,v,k,j,i) - 0.5*sgn*recv;
        } else {
          flx.x3f(m,v,k,j,i) = 0.5*flx.x3f(m,v,k,j,i) - 0.5*sgn*recv;
        }
      });
    });
  });
  return TaskStatus::complete;
}

