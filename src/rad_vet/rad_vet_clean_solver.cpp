//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_clean_solver.cpp
//! \brief Linear algebra of the clean implicit radiation solve: the operator
//! A E = E (1 + h c kE (1-r)) + h div F_lin(E), x1-line tridiagonal / multigrid
//! preconditioners, right-preconditioned BiCGStab, ghost exchange, layout-invariant sums.

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "rad_vet/rad_vet.hpp"
#include "utils/repro_sum.hpp"
#include "rad_vet/rad_vet_clean.hpp"
#include "utils/scratch_view.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

// layout-invariant reductions of the clean solve (repro sums)
namespace radvet {

//! up to three sums and one max (VetLayReduce)
struct VetBcgVal {
  Real s0, s1, s2, mx;
};

constexpr int VET_DR_NT = 256;    // slots per team (a power of 2)
// block mode: level-1 slots per MeshBlock (fills the GPU)
constexpr int VET_LB_NS = 32768;

//----------------------------------------------------------------------------------------
//! \fn VetGlobalMax
inline void VetGlobalMax(Real &a) {
#if MPI_PARALLEL_ENABLED
  Real g;
  MPI_Allreduce(&a, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  a = g;
#else
  (void)a;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn VetLayReduce
//! \brief <mesh>/layout_invariant: the VetBcgVal reduction of fn(idx, v) over idx in
//! [0, n) and all ranks with exact order-independent sums (repro::ReproSum).
//! Cell mode: every cell term deposited exactly (independent of layout and backend).
//! Block mode (default): each MeshBlock's (or tile's) partial in a fixed order
//! (VET_LB_NS strided slots, then a fixed pairwise tree), partials summed exactly on the
//! host and merged over ranks.  bpart / hpart: device / host scratch.

template <class F>
void VetLayReduce(const char *name, const int n, const int nkji, DvceArray1D<Real> &bpart,
                 Kokkos::View<Real*, Kokkos::HostSpace> &hpart, const F &fn,
                 VetBcgVal &out, const int *nxb = nullptr) {
  if (repro::block_mode) {
    // level 1: slot s of VET_LB_NS sums q = s, s + VET_LB_NS, ...; level 2: one team per
    // unit sums slots t, t + VET_DR_NT, ..., then the pairwise tree.  Tile mode (nxb):
    // the units are N x N (x2, x3) tiles (a MeshBlock-sized tile is bitwise the block)
    const bool tl = repro::tile_n > 0 && nxb != nullptr;
    const int t2 = (tl && nxb[1] > 1) ? repro::tile_n : 1;
    const int t3 = (tl && nxb[2] > 1) ? repro::tile_n : 1;
    const int n1 = tl ? nxb[0] : 1, n2 = tl ? nxb[1] : 1, nt2 = n2/t2;
    const int ntb = tl ? nt2*(nxb[2]/t3) : 1;
    const int nmb = (n/nkji)*ntb;   // (the units)
    const int ntc = n1*t2*t3;   // (tile mode: the tile's cells)
    const int ns = (tl && ntc < VET_LB_NS) ? ntc : VET_LB_NS;
    const int off2 = 4*nmb*ns;
    if (bpart.extent_int(0) < off2 + 4*nmb) {
      bpart = DvceArray1D<Real>(
          Kokkos::view_alloc("vet_lay_part", Kokkos::WithoutInitializing), off2 + 4*nmb);
    }
    if (hpart.extent_int(0) < 4*nmb) {
      hpart = Kokkos::View<Real*, Kokkos::HostSpace>("vet_lay_part_h", 4*nmb);
    }
    auto part = bpart;
    Kokkos::parallel_for(name, Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nmb*ns),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/ns, sl = idx - m*ns;
      VetBcgVal v;
      v.s0 = v.s1 = v.s2 = v.mx = 0.0;
      if (tl) {
        const int b = m/ntb, u = m - b*ntb, k0 = (u/nt2)*t3, j0 = (u%nt2)*t2;
#pragma unroll 4
        for (int q = sl; q < ntc; q += ns) {
          const int r = q/n1, k = k0 + r/t2, j = j0 + r%t2;
          fn(b*nkji + (k*n2 + j)*n1 + q - r*n1, v);
        }
      } else {
#pragma unroll 4
        for (int q = sl; q < nkji; q += ns) {fn(m*nkji + q, v);}
      }
      part(4*idx) = v.s0;
      part(4*idx + 1) = v.s1;
      part(4*idx + 2) = v.s2;
      part(4*idx + 3) = v.mx;
    });
    using TP = Kokkos::TeamPolicy<DevExeSpace>;
    using Mem = Kokkos::View<Real*, ScratchMemSpace,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    const size_t sb = Mem::shmem_size(4*VET_DR_NT);
    Kokkos::parallel_for("vet_lay_red2", TP(DevExeSpace(), nmb, Kokkos::AUTO)
                         .set_scratch_size(0, Kokkos::PerTeam(sb)),
    KOKKOS_LAMBDA(const TP::member_type &tm) {
      Mem sh(tm.team_scratch(0), 4*VET_DR_NT);
      const int m = tm.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, VET_DR_NT), [&](const int t) {
        Real a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        // (kU slots' values loaded at once, then added in the slot order: bitwise)
        constexpr int kU = 8;
        int q = t;
        for (; q + (kU - 1)*VET_DR_NT < ns; q += kU*VET_DR_NT) {
          Real p[kU][4];
          for (int u = 0; u < kU; ++u) {
            for (int c = 0; c < 4; ++c) p[u][c] = part(4*(m*ns + q + u*VET_DR_NT) + c);
          }
          for (int u = 0; u < kU; ++u) {
            a0 += p[u][0];
            a1 += p[u][1];
            a2 += p[u][2];
            a3 = (p[u][3] > a3) ? p[u][3] : a3;
          }
        }
        for (; q < ns; q += VET_DR_NT) {
          const int b = 4*(m*ns + q);
          a0 += part(b);
          a1 += part(b + 1);
          a2 += part(b + 2);
          a3 = (part(b + 3) > a3) ? part(b + 3) : a3;
        }
        sh(4*t) = a0;
        sh(4*t + 1) = a1;
        sh(4*t + 2) = a2;
        sh(4*t + 3) = a3;
      });
      tm.team_barrier();
      for (int w = VET_DR_NT/2; w > 0; w /= 2) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, w), [&](const int t) {
          sh(4*t) += sh(4*(t + w));
          sh(4*t + 1) += sh(4*(t + w) + 1);
          sh(4*t + 2) += sh(4*(t + w) + 2);
          sh(4*t + 3) = (sh(4*(t + w) + 3) > sh(4*t + 3)) ? sh(4*(t + w) + 3)
                                                          : sh(4*t + 3);
        });
        tm.team_barrier();
      }
      Kokkos::single(Kokkos::PerTeam(tm), [&]() {
        for (int c = 0; c < 4; ++c) {part(off2 + 4*m + c) = sh(c);}
      });
    });
    auto hp = Kokkos::subview(hpart, Kokkos::make_pair(0, 4*nmb));
    // (stream-ordered copy + one fence)
    Kokkos::deep_copy(DevExeSpace(), hp, Kokkos::subview(part, Kokkos::make_pair(off2,
                      off2 + 4*nmb)));
    DevExeSpace().fence();
    repro::ReproSum rs[3];
    for (int c = 0; c < 3; ++c) {rs[c].Init();}
    Real mx = 0.0;
    for (int m = 0; m < nmb; ++m) {
      for (int c = 0; c < 3; ++c) {rs[c].Add(hp(4*m + c));}
      mx = (hp(4*m + 3) > mx) ? hp(4*m + 3) : mx;
    }
    repro::AllMergeMax(rs, 3, mx);
    out.s0 = rs[0].Value();
    out.s1 = rs[1].Value();
    out.s2 = rs[2].Value();
    out.mx = mx;
    return;
  }
  repro::ReproVal<3> rv;
  Kokkos::parallel_reduce(name, Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, n),
  KOKKOS_LAMBDA(const int idx, repro::ReproVal<3> &a) {
    VetBcgVal t;
    t.s0 = t.s1 = t.s2 = t.mx = 0.0;
    fn(idx, t);
    a.s[0].Add(t.s0);
    a.s[1].Add(t.s1);
    a.s[2].Add(t.s2);
    a.mx = (t.mx > a.mx) ? t.mx : a.mx;
  }, repro::ReproRed<3, Kokkos::HostSpace>(rv));
  repro::AllMerge(rv.s, 3);
  VetGlobalMax(rv.mx);
  out.s0 = rv.s[0].Value();
  out.s1 = rv.s[1].Value();
  out.s2 = rv.s[2].Value();
  out.mx = rv.mx;
}

}  // namespace radvet

namespace radvet {
//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::MeshBlock0Dx
//! \brief cell width along x(d+1) of the mesh's MeshBlock 0: the same on every rank
//! (a rank's own first MeshBlock differs by round-off)

Real RadiationVet::MeshBlock0Dx(const int d) const {
  auto *pm = pmy_pack->pmesh;
  auto &ms = pm->mesh_size;
  auto &id = pm->mb_indcs;
  const int lv = pm->lloc_eachmb[0].level - pm->root_level;
  const int nr = ((d == 0) ? pm->nmb_rootx1 :
                  ((d == 1) ? pm->nmb_rootx2 : pm->nmb_rootx3)) << lv;
  const Real a = (d == 0) ? ms.x1min : ((d == 1) ? ms.x2min : ms.x3min);
  const Real b = (d == 0) ? ms.x1max : ((d == 1) ? ms.x2max : ms.x3max);
  const int n = (d == 0) ? id.nx1 : ((d == 1) ? id.nx2 : id.nx3);
  return (((nr > 1) ? LeftEdgeX(1, nr, a, b) : b) - a)/static_cast<Real>(n);
}

//----------------------------------------------------------------------------------------
//! \fn CleanArgs RadiationVet::CleanMakeArgs
//! \brief collects the arrays and constants the device functions need

CleanArgs RadiationVet::CleanMakeArgs() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  CleanArgs x;
  x.cs = cln_s;
  x.w = cln_w;
  x.fc[0] = cln_fc1;
  x.fc[1] = cln_fc2;
  x.fc[2] = cln_fc3;
  x.beam = ibm_f;
  x.bstress = ibm_s;
  x.h = dt_sub;
  x.c = c_light;
  for (int d = 0; d < 3; ++d) x.idx[d] = 1.0/MeshBlock0Dx(d);
  const int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  // an outflow face is the open one without an entering bath
  for (int s = 0; s < 6; ++s) {
    x.bc[s] = (bc[s] == VET_IBC_X2FLUX) ? VET_IBC_FLUX :
              ((bc[s] == VET_IBC_OUTFLOW) ? VET_IBC_MARSHAK : bc[s]);
  }
  x.iflux[0] = iflux_x1min;
  x.iflux[1] = iflux_x1max;
  x.ifx2[0] = ifx2_lo;
  x.ifx2[1] = ifx2_hi;
  x.mbc = pmy_pack->pmb->mb_bcs.d_view;
  for (int n = 0; n < 6; ++n) x.ebath[n] = 0.0;
  x.ebath[0] = iebath_x1min;
  x.ebath[1] = iebath_x1max;
  for (int s = 0; s < 2; ++s) x.ebath[s] *= (bc[s] != VET_IBC_OUTFLOW);
  x.marshak_q = marshak_q;
  x.marshak = (cln_bdry_form == 1);
  x.use_beam = (ibm_n > 0);
  x.feedback = cln_feedback;
  x.vzero = cln_vzero;
  x.hi = false;
  x.mcl = false;
  x.fa = false;
  x.fwc = false;
  x.fs = false;
  x.tf = false;
  x.fd = cln_fsolve && cln_fkin;
  x.gf = cln_gflux;
  x.fre = cln_freal;
  x.fdg = cln_fdg;
  x.xq = cln_xd_on ? cln_xqs : -1;
  x.pd = cln_ucp;
  x.xlag = cln_xlag;
  x.xsrc = cln_xsrc;
  x.drag = cln_vimp && cln_feedback;
  x.upw = cln_upw;
  x.ucen = cln_ucen;
  x.ucp = cln_ucp;
  x.sar = !cln_gflux || cln_farea;
  x.bar = cln_bcarea;
  x.op = opac;
  x.nf = pmy_pack->pmesh->three_d ? 3 : (pmy_pack->pmesh->multi_d ? 2 : 1);
  x.is = indcs.is;
  x.ie = indcs.ie;
  x.js = indcs.js;
  x.je = indcs.je;
  x.ks = indcs.ks;
  x.ke = indcs.ke;
  x.multi_d = pmy_pack->pmesh->multi_d;
  x.three_d = pmy_pack->pmesh->three_d;
  x.sph = pmy_pack->pmesh->use_spherical_polar;
  auto &ms = pmy_pack->pmesh->mesh_size;
  x.cseam = (pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic)
            ? 0.5/tan(ms.x2min) + 0.5/tan(ms.x2max) : 0.0;   // (read on that wrap only)
  x.cub = pmy_pack->pmesh->use_cubed_sphere;
  x.ys = !x.cub && cln_ysum;
  x.lwb = x.cub && cln_latwb;
  x.yl = cln_yl;
  if (x.sph) {
    auto *pc = pmy_pack->pcoord;
    x.vol = pc->volume;
    x.area[0] = pc->area.x1f;
    x.area[1] = pc->area.x2f;
    x.area[2] = pc->area.x3f;
    x.dxf[0] = pc->dxface.x1f;
    x.dxf[1] = pc->dxface.x2f;
    x.dxf[2] = pc->dxface.x3f;
    x.x1v = pc->x1v;
    x.x1f = pc->xx1f;
    x.x2v = pc->x2v;
    x.x2f = pc->xx2f;
    const int nmb = pmy_pack->nmb_thispack, n2 = pc->x2v.extent_int(1);
    if (x.sph && (cln_tg.extent_int(0) != nmb || cln_tg.extent_int(2) != n2 + 1)) {
      Kokkos::realloc(cln_tg, nmb, 4, n2 + 1);   // (once: the trig of theta, per j)
      auto tg = cln_tg;
      auto x2v = pc->x2v, x2f = pc->xx2f;
      par_for("cln_tg", DevExeSpace(), 0, nmb - 1, 0, n2,
      KOKKOS_LAMBDA(const int m, const int j) {
        if (j < n2) {
          tg(m,0,j) = tan(x2v(m,j));
          tg(m,2,j) = sin(x2v(m,j));
        }
        tg(m,1,j) = tan(x2f(m,j));
        tg(m,3,j) = sin(x2f(m,j));
      });
    }
    x.tg = cln_tg;
    if (vgd_wst) {
      if (cln_sgv.extent_int(0) != nmb) {CleanSeamImage();}
      x.vol = cln_sgv;
      x.area[0] = cln_sga0;
      x.area[1] = cln_sga1;
      x.dxf[1] = cln_sgd1;
      x.dxf[2] = cln_sgd2;
      x.x2v = cln_sgx;
      x.tg = cln_sgtg;
    }
  }
  if (x.cub) {
    x.vol = pmy_pack->pcoord->volume;
    x.ae0 = pmy_pack->pcoord->area.x1f;
    for (int d = 0; d < 3; ++d) {
      x.area[d] = cln_ga[d];
      x.dxf[d] = cln_gd[d];
    }
    x.g = cln_gn;
    x.x1v = pmy_pack->pcoord->x1v;
    x.x1f = pmy_pack->pcoord->xx1f;
  }
  if (!x.sph && !x.cub && cln_lzsp) {   // (clean_lz_spin, Cartesian: z x r of cells)
    const int nmb = pmy_pack->nmb_thispack, ng = indcs.ng;
    const int n1 = indcs.nx1 + 2*ng, n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
    if (cln_cx[0].extent_int(0) != nmb) {
      Kokkos::realloc(cln_cx[0], nmb, n1);
      Kokkos::realloc(cln_cx[1], nmb, n2);
      auto cx = cln_cx[0], cy = cln_cx[1];
      auto size = pmy_pack->pmb->mb_size.d_view;
      const int is = indcs.is, js = indcs.js, nx1 = indcs.nx1, nx2 = indcs.nx2;
      par_for("cln_cx", DevExeSpace(), 0, nmb - 1, 0, std::max(n1, n2) - 1,
      KOKKOS_LAMBDA(const int m, const int i) {
        if (i < n1) cx(m,i) = CellCenterX(i - is, nx1, size(m).x1min, size(m).x1max);
        if (i < n2) cy(m,i) = CellCenterX(i - js, nx2, size(m).x2min, size(m).x2max);
      });
    }
    x.x1v = cln_cx[0];
    x.x2v = cln_cx[1];
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  x.etg = fl.on && fl.use_etotgrav;
  x.mhd = fl.on && fl_mhd;
  x.phi = fl.phicc0;
  x.emag = emag0;
  return x;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanSeamImage
//! \brief vet_gd_theta_seam_sym: the geometry the clean scheme reads, with the ghost rows
//! of a theta-periodic seam (MeshBlocks at x2min / x2max) taking that of the cells they
//! hold, their periodic images.  On a theta wedge centred on the equator the image of the
//! ghost at depth g beyond a wall is the equatorial mirror of the row at depth g inside
//! it: vol, the x1 / x3 face lengths and areas of that row, the x2 faces mirrored, the
//! centroid reflected about the wall (as BuildX2vRecon).  The seam face's two sides then
//! have mirror-paired stencils, as they have across any interior face; without it the
//! one copy of the face computed on the x2min side sees a ghost nearer the pole than the
//! cell it holds (a net F_theta that grows with the wedge)

void RadiationVet::CleanSeamImage() {
  auto *pc = pmy_pack->pcoord;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je, ng = indcs.ng, nmb = pmy_pack->nmb_thispack;
  if (indcs.nx2 < ng) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "vet_gd_theta_seam_sym needs nx2 >= nghost per MeshBlock" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto cp = [](const DvceArray4D<Real> &a, const char *nm) {
    DvceArray4D<Real> b(std::string(nm), a.extent(0), a.extent(1), a.extent(2),
                        a.extent(3));
    Kokkos::deep_copy(b, a);
    return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);
  };
  auto v = cp(pc->volume, "cln_sgv"), a0 = cp(pc->area.x1f, "cln_sga0");
  auto a1 = cp(pc->area.x2f, "cln_sga1"), d1 = cp(pc->dxface.x2f, "cln_sgd1");
  auto d2 = cp(pc->dxface.x3f, "cln_sgd2");
  cln_sgx = DvceArray2D<Real>("cln_sgx", pc->x2v.extent(0), pc->x2v.extent(1));
  cln_sgtg = DvceArray3D<Real>("cln_sgtg", cln_tg.extent(0), cln_tg.extent(1),
                               cln_tg.extent(2));
  Kokkos::deep_copy(cln_sgx, pc->x2v);
  Kokkos::deep_copy(cln_sgtg, cln_tg);
  auto x2v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cln_sgx);
  auto tg = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cln_sgtg);
  auto x2f = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->xx2f);
  auto x1v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->x1v);
  auto bc = pmy_pack->pmb->mb_bcs.h_view;
  const BoundaryFlag pf = BoundaryFlag::periodic;
  for (int m = 0; m < nmb; ++m) {
    for (int sd = 0; sd < 2; ++sd) {
      if (bc(m, sd ? BoundaryFace::outer_x2 : BoundaryFace::inner_x2) != pf) {continue;}
      const int w = sd ? je + 1 : js;   // the seam face
      for (int g = 0; g < ng; ++g) {
        const int jg = sd ? je + 1 + g : js - 1 - g, jm = sd ? je - g : js + g;
        // (cells; the x2 face beyond the ghost, depth g + 1 from the wall)
        const int fg = sd ? je + 2 + g : js - 1 - g, fm = sd ? je - g : js + 1 + g;
        x2v(m,jg) = 2.0*x2f(m,w) - x2v(m,jm);
        tg(m,0,jg) = tg(m,0,jm);   // (tan, sin of the cell and of its outer x2 face)
        tg(m,2,jg) = tg(m,2,jm);
        tg(m,1,fg) = tg(m,1,fm);
        tg(m,3,fg) = tg(m,3,fm);
        for (int k = 0; k < static_cast<int>(v.extent(1)); ++k) {
          for (int i = 0; i < static_cast<int>(a0.extent(3)); ++i) {
            if (i < static_cast<int>(v.extent(3))) {
              v(m,k,jg,i) = v(m,k,jm,i);
              a1(m,k,fg,i) = a1(m,k,fm,i);
              d2(m,k,jg,i) = d2(m,k,jm,i);
            }
            a0(m,k,jg,i) = a0(m,k,jm,i);
          }
        }
        for (int i = 0; i < static_cast<int>(v.extent(3)); ++i) {   // (d2: k = nk too)
          d2(m,d2.extent(1)-1,jg,i) = d2(m,d2.extent(1)-1,jm,i);
        }
      }
      // the centre distances across the seam face and the ghosts' own x2 faces
      for (int g = 0; g <= ng - 1; ++g) {
        const int fg = sd ? w + g : w - g;
        for (int k = 0; k < static_cast<int>(d1.extent(1)); ++k) {
          for (int i = 0; i < static_cast<int>(d1.extent(3)); ++i) {
            d1(m,k,fg,i) = x1v(m,i)*(x2v(m,fg) - x2v(m,fg-1));
          }
        }
      }
    }
  }
  auto dv = [](auto &h, DvceArray4D<Real> &d, const char *nm) {
    d = DvceArray4D<Real>(std::string(nm), h.extent(0), h.extent(1), h.extent(2),
                          h.extent(3));
    Kokkos::deep_copy(d, h);
  };
  dv(v, cln_sgv, "cln_sgv");
  dv(a0, cln_sga0, "cln_sga0");
  dv(a1, cln_sga1, "cln_sga1");
  dv(d1, cln_sgd1, "cln_sgd1");
  dv(d2, cln_sgd2, "cln_sgd2");
  Kokkos::deep_copy(cln_sgx, x2v);
  Kokkos::deep_copy(cln_sgtg, tg);
}

//----------------------------------------------------------------------------------------
//! \fn const void *RadiationVet::CleanDevCopy
//! \brief device copy (slot s < 8, n <= 8 KB) of an argument struct for non-inlined
//! device calls: read in place, not copied to every thread's stack as a capture is

const void *RadiationVet::CleanDevCopy(const void *p, const size_t n, const int s) {
  // (host builds too: the copy outlives the struct)
  constexpr size_t slot = 8192;
  if (n > slot || s < 0 || s >= 8) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "CleanDevCopy: " << n << " bytes in slot " << s
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (cln_xdev.extent(0) == 0) Kokkos::realloc(cln_xdev, 8*slot);
  Kokkos::View<const char *, HostMemSpace, Kokkos::MemoryUnmanaged>
      h(static_cast<const char *>(p), n);
  auto d = Kokkos::subview(cln_xdev, std::make_pair(s*slot, s*slot + n));
#if defined(KOKKOS_ENABLE_CUDA)
  // (CUDA: in stream order, no fence; a small pageable copy is staged at the call)
  Kokkos::deep_copy(DevExeSpace(), d, h);
#else
  Kokkos::deep_copy(d, h);
#endif
  return cln_xdev.data() + s*slot;
}

const CleanArgs *RadiationVet::CleanDevArgs(const CleanArgs &x, const int s) {
  return CleanDev(x, s);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanGhostFill
//! \brief fill the ghosts of components [c0, c0+nc) by the cell-centred exchange
//! (cubed-sphere seams without resample), then (phys) copy the edge cell on physical
//! boundaries.  part 1 starts, part 2 finishes (inner cells may run between).  Deferred
//! fills (CleanGhostDefer) join the next exchange where they fit.  D12, D13, D23 are odd
//! across a reflecting wall normal to one of their axes and D_r,theta, D_r,phi across
//! the pole, so d_b(D E) puts no shear traction on the wall

void RadiationVet::CleanGhostFill(DvceArray5D<Real> a, const int c0, const int nc,
                                  const bool phys, const int part, const int dep,
                                  const std::array<uint32_t, 4> fm) {
  CleanGhostFillT(a, c0, nc, phys, part, dep, fm);
}

void RadiationVet::CleanGhostDefer(DvceArray5D<Real> a, const int c0, const int nc,
                                   const int dep) {
  int np = nc;
  for (auto &s : cln_hpend) np += s.nc;
  if (np > cln_hmax) CleanGhostFlush();
  cln_hpend.push_back({a, c0, nc, true, dep});
}

// (a range lies in one of the two arrays: CLW_NS)
void RadiationVet::CleanGhostFill(const CleanW &a, const int c0, const int nc,
                                  const int part, const int dep,
                                  const std::array<uint32_t, 4> fm) {
  CleanGhostFillT((c0 < a.ns) ? a.s : a.o, (c0 < a.ns) ? c0 : c0 - a.ns, nc, true, part,
                  dep, fm);
}

void RadiationVet::CleanGhostDefer(const CleanW &a, const int c0, const int nc,
                                   const int dep) {
  CleanGhostDefer((c0 < a.ns) ? a.s : a.o, (c0 < a.ns) ? c0 : c0 - a.ns, nc, dep);
}

void RadiationVet::CleanGhostFlush() {
  if (cln_hpend.empty()) return;
  std::vector<CleanGSeg<Real>> sg;
  sg.swap(cln_hpend);
  CleanGhostExchange(sg, 0);
}

template <typename T>
void RadiationVet::CleanGhostFillT(DvceArray5D<T> a, const int c0, const int nc,
                                   const bool phys, const int part, const int dep,
                                   const std::array<uint32_t, 4> fm) {
  if (nc > cln_hmax) {   // (the exchange moves at most cln_hmax variables at once)
    const int h = cln_hmax;
    // (a full mask stays full for every chunk)
    auto sh = [&](const uint32_t m) {return (m == ~0u || h >= 32) ? m : m >> h;};
    CleanGhostFillT(a, c0, h, phys, 0, dep, fm);
    CleanGhostFillT(a, c0 + h, nc - h, phys, 0, dep, {sh(fm[0]), sh(fm[1]), sh(fm[2]),
                                                       sh(fm[3])});
    return;
  }
  std::vector<CleanGSeg<T>> sg;
  if constexpr (std::is_same<T, Real>::value) {
    int np = nc;
    for (auto &s : cln_hpend) np += s.nc;
    // (part 2: the segments part 1 left in cln_hpend, nothing deferred in between)
    if (part == 2 || np <= cln_hmax) sg.swap(cln_hpend);
    if (part == 2) return CleanGhostExchange(sg, 2);
  }
  CleanGhostFlush();   // (those that do not fit: on their own first)
  sg.push_back({a, c0, nc, phys, dep, fm});
  CleanGhostExchange(sg, part);
  if constexpr (std::is_same<T, Real>::value) {
    if (part == 1) cln_hpend.swap(sg);
  }
}

template <typename T>
void RadiationVet::CleanGhostExchange(const std::vector<CleanGSeg<T>> &sg,
                                      const int part) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ng = indcs.ng;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = sg[0].a.extent_int(4), n2 = sg[0].a.extent_int(3);
  const int n3 = sg[0].a.extent_int(2);
  int nc = 0;
  for (auto &s : sg) nc += s.nc;
  if (!cln_hpb) {
    cln_hpb.reset(new MeshBoundaryValuesCC(pmy_pack, nullptr, false));
    cln_hpb->InitializeBuffers(cln_hmax);
    cln_hpb->SetVectorPairs(cln_hmax, {});
    cln_hpb->cs_noresample = true;
    // (16 teams per buffer instead of one per MeshBlock face; bitwise)
    cln_hpb->cc_split = 16;
  }
  if (cln_hdl.state == 0) CleanHaloInit(n3, n2, n1);
  const bool dl = (cln_hdl.state == 1);
  if (dl && (n1 != cln_hdl.n1 || n2 != cln_hdl.n2 || n3 != cln_hdl.n3)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "CleanGhostExchange: an array of another shape than the "
              << "lists of clean_halo_direct" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (dl && part != 2) {   // (the depth: the deepest its segments need)
    int d = 0;
    for (auto &s : sg) d = (s.dep <= 0) ? ng : std::max(d, s.dep);
    cln_hdl.dcur = std::min(std::max(d, 1), ng);
    CleanHaloPost(sg);
  }
  if (dl && part == 1) return;
  if (dl) CleanHaloFinish(sg);
  const bool cs = pmy_pack->pmesh->use_cubed_sphere;
  // one MeshBlock and one Real array: exchange on the array itself, no copies
  bool own = false;
  if constexpr (std::is_same<T, Real>::value) own = (nmb1 == 0 && sg.size() == 1);
  // (clean_halo_direct, doubles: the cube-vertex corners on the arrays themselves)
  const bool dc = dl && cs && std::is_same<T, Real>::value;
  // (exchanges run one at a time: their scratch shares one allocation, grown on demand)
  if (!own && (dl ? cs && !dc : true) && cln_hbuf.extent_int(1) < nc) {
    Kokkos::realloc(cln_hbuf, nmb1 + 1, nc, n3, n2, n1);
    cln_halo.clear();
  }
  CleanHalo &hl = cln_halo[nc];
  if (hl.ca.extent_int(0) == 0) {
    hl.a = DvceArray5D<Real>(cln_hbuf.data(), nmb1 + 1, nc, n3, n2, n1);
    Kokkos::realloc(hl.ca, nmb1 + 1, nc, 1, 1, 1);
  }
  auto h = hl.a;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  MeshBoundaryValuesCC *pb = cln_hpb.get();
  DvceArray5D<Real> ha = hl.a;
  if constexpr (std::is_same<T, Real>::value) {
    const size_t off = static_cast<size_t>(sg[0].c0)*n3*n2*n1;
    if (own) ha = DvceArray5D<Real>(sg[0].a.data() + off, 1, nc, n3, n2, n1);
  }
  // (only the shell outside the inner box is copied; clean_halo_direct: only on the
  // cubed sphere, for the cube-vertex corners)
  if (!own && !dc && (dl ? cs : part != 2)) {
    const CleanShell sc(n3, n2, n1, ng, is, ie, js, je, ks, ke, md, td);
    int o = 0;
    for (auto &s : sg) {
      auto a = s.a;
      const int c0 = s.c0;
      if (sc.Empty()) {
        par_for("cln_hcopy", DevExeSpace(), 0, nmb1, 0, s.nc-1, 0, n3-1, 0, n2-1, 0, n1-1,
        KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
          h(m,o+n,k,j,i) = a(m,c0+n,k,j,i);
        });
      } else {
        par_for("cln_hcopy", DevExeSpace(), 0, nmb1, 0, s.nc-1, 0, sc.ncell-1,
        KOKKOS_LAMBDA(const int m, const int n, const int q) {
          int k, j, i;
          CleanShellCell(sc, q, k, j, i);
          h(m,o+n,k,j,i) = a(m,c0+n,k,j,i);
        });
      }
      o += s.nc;
    }
  }
  if (dc) {
    if constexpr (std::is_same<T, Real>::value) {
      for (auto &s : sg) {
        DvceArray5D<Real> a = s.a;
        pb->FillPanelCornersCC(a, false, 0, s.c0, s.nc);
      }
    }
  } else if (dl && cs) {
    pb->FillPanelCornersCC(ha);
  } else if (!dl) {
    if (part != 2) {
      while (pb->InitRecv(nc) == TaskStatus::incomplete) {}
      while (pb->PackAndSendCC(ha, hl.ca) == TaskStatus::incomplete) {}
      if (part == 1) return;
    }
    while (pb->RecvAndUnpackCC(ha, hl.ca) == TaskStatus::incomplete) {}
    while (pb->ClearSend() == TaskStatus::incomplete) {}
    while (pb->ClearRecv() == TaskStatus::incomplete) {}
  }
  const CleanShell sb(n3, n2, n1, 0, is, ie, js, je, ks, ke, md, td);
  int o = 0;
  for (auto &s : sg) {
    if (own || (dl && !cs) || dc) break;
    auto a = s.a;
    const int c0 = s.c0;
    par_for("cln_hback", DevExeSpace(), 0, nmb1, 0, s.nc-1, 0, sb.ncell-1,
    KOKKOS_LAMBDA(const int m, const int n, const int q) {
      int k, j, i;
      CleanShellCell(sb, q, k, j, i);
      a(m,c0+n,k,j,i) = h(m,o+n,k,j,i);
    });
    o += s.nc;
  }
  for (auto &s : sg) {
    if (s.phys) CleanGhostPhys(s.a, s.c0, s.nc);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanHaloInit
//! \brief clean_halo_direct: build the ghost exchange's cell lists by sending cell ids
//! through cln_hpb: on-rank pairs make the copy kernel, off-rank ghosts are grouped by
//! rank (pack lists in the receiver's order, no headers).  Self-check: a test field must
//! come out bitwise equal through both paths, else the run stops

void RadiationVet::CleanHaloInit(const int n3, const int n2, const int n1) {
  CleanHaloList &q = cln_hdl;
  auto *pm = pmy_pack->pmesh;
  const int me = global_variable::my_rank;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const bool md = pm->multi_d, td = pm->three_d;
  const int nmb = pmy_pack->nmb_thispack, nrk = global_variable::nranks;
  int mx = 0;
  for (int r = 0; r < nrk; ++r) mx = std::max(mx, pm->nmb_eachrank[r]);
  const int n21 = n2*n1, n321 = n3*n21;
  const int64_t nid = static_cast<int64_t>(mx)*n321;
  q.n1 = n1; q.n2 = n2; q.n3 = n3;
  const int ng = indcs.ng;
  q.ng = ng;
  DvceArray5D<Real> x("cln_hdl_x", nmb, 1, n3, n2, n1),
                    y("cln_hdl_y", nmb, 1, n3, n2, n1);
  const std::vector<CleanGSeg<Real>> sx = {{x, 0, 1, false}}, sy = {{y, 0, 1, false}};
  // (pass 0: the ids; pass 1: other values, the check)
  for (int pass = 0; pass < 2; ++pass) {
    par_for("cln_hd_id", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real v = static_cast<Real>((me*nid + m*n321) + (k*n2 + j)*n1 + i);
      x(m,0,k,j,i) = (pass == 0) ? v : 2.0 + Kokkos::sin(1.0e-3*v);
      y(m,0,k,j,i) = x(m,0,k,j,i);
    });
    q.state = -1;
    CleanGhostExchange(sx, 0);
    if (pass == 1) {
      q.state = 1;
      CleanGhostExchange(sy, 0);
      break;
    }
    auto hx = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x);
    std::vector<int> ld, ls, lfm;
    std::vector<std::vector<int>> rg(nrk), rs(nrk), rf(nrk);   // off-rank ghosts, their
                                                                // sources and families
    // (family of a ghost copying id w: the axis along which its neighbour one layer in
    // differs by one cell; 3 at edges and corners)
    auto fam = [&](const int m, const int k, const int j, const int i, const int64_t w) {
      const int ox = (i < is || i > ie), oy = md && (j < js || j > je);
      const int oz = td && (k < ks || k > ke);
      if (ox + oy + oz != 1 || ng < 2) return 3;
      const int a = oy ? 1 : (oz ? 2 : 0), c = (a == 0) ? i : ((a == 1) ? j : k);
      const int lo = (a == 0) ? is : ((a == 1) ? js : ks);
      const int d = (c < lo) ? lo - c : c - ((a == 0) ? ie : ((a == 1) ? je : ke));
      const int t = ((c < lo) == (d > 1)) ? 1 : -1;
      const int64_t v = static_cast<int64_t>(hx(m,0,k+t*(a == 2),j+t*(a == 1),
                                                i+t*(a == 0)));
      const int64_t dw = (v > w) ? v - w : w - v;
      return (dw == 1) ? 0 : ((dw == n1) ? 1 : ((dw == n21) ? 2 : 3));
    };
    for (int m = 0; m < nmb; ++m) {
      for (int k = 0; k < n3; ++k) {
        for (int j = 0; j < n2; ++j) {
          for (int i = 0; i < n1; ++i) {
            const int g = (k*n2 + j)*n1 + i;
            const Real v = hx(m,0,k,j,i);
            const int64_t w = static_cast<int64_t>(v);
            if (v == static_cast<Real>(me*nid + m*n321 + g) || v < 0.0 ||
                static_cast<Real>(w) != v) continue;   // (not filled, or not a copy)
            const int r = static_cast<int>(w/nid);
            const int ms = static_cast<int>((w - r*nid)/n321);
            const int f = static_cast<int>(w - r*nid - ms*n321);
            const int kk = f/n21, jj = (f - kk*n21)/n1, ii = f - kk*n21 - jj*n1;
            if (r >= nrk || ms >= pm->nmb_eachrank[r] || ii < is || ii > ie ||
                (md && (jj < js || jj > je)) || (td && (kk < ks || kk > ke))) continue;
            if (r == me) {
              ld.push_back(m*n321 + g);
              ls.push_back(ms*n321 + f);
              lfm.push_back(fam(m, k, j, i, w));
            } else {
              rg[r].push_back(m*n321 + g);
              rs[r].push_back(ms*n321 + f);
              rf[r].push_back(fam(m, k, j, i, w));
            }
          }
        }
      }
    }
    // (each list in order of depth, the ghost layer of the receiving cell)
    auto dep = [&](const int f) {
      const int r = f%n321, k = r/n21, j = (r - k*n21)/n1, i = r - k*n21 - j*n1;
      const int d = std::max({is - i, i - ie, js - j, j - je, ks - k, k - ke});
      return std::min(std::max(d, 1), ng);
    };
    std::vector<int> nrd(static_cast<size_t>(nrk)*ng, 0), lrd, ldd;
    auto srt = [&](std::vector<int> &g, std::vector<int> &c, std::vector<int> &fv, int *n,
                   std::vector<int> &dv) {
      std::vector<int> o(g.size()), g2, c2, f2;
      for (size_t a = 0; a < o.size(); ++a) o[a] = static_cast<int>(a);
      std::stable_sort(o.begin(), o.end(),
                       [&](int a, int b) {return dep(g[a]) < dep(g[b]);});
      for (int a : o) {
        g2.push_back(g[a]);
        c2.push_back(c[a]);
        f2.push_back(fv[a]);
        dv.push_back(dep(g[a]));
        for (int d = dv.back(); d <= ng; ++d) ++n[d-1];
      }
      g.swap(g2);
      c.swap(c2);
      fv.swap(f2);
    };
    q.nld.assign(ng, 0);
    srt(ld, ls, lfm, q.nld.data(), ldd);
    for (int r = 0; r < nrk; ++r) srt(rg[r], rs[r], rf[r], nrd.data() + r*ng, lrd);
    std::vector<int> sc, rc;
#if MPI_PARALLEL_ENABLED
    MPI_Comm_dup(MPI_COMM_WORLD, &q.comm);
    std::vector<int> nr(nrk), ns(nrk);
    for (int r = 0; r < nrk; ++r) nr[r] = static_cast<int>(rs[r].size());
    MPI_Alltoall(nr.data(), 1, MPI_INT, ns.data(), 1, MPI_INT, q.comm);
    std::vector<std::vector<int>> ss(nrk), sf(nrk);
    std::vector<MPI_Request> rq;
    for (int r = 0; r < nrk; ++r) {
      if (nr[r] == 0 && ns[r] == 0) continue;
      q.rank.push_back(r);
      ss[r].resize(ns[r]);
      sf[r].resize(ns[r]);
      rq.resize(rq.size() + 4);
      MPI_Irecv(ss[r].data(), ns[r], MPI_INT, r, 4244, q.comm, &rq[rq.size() - 4]);
      MPI_Isend(rs[r].data(), nr[r], MPI_INT, r, 4244, q.comm, &rq[rq.size() - 3]);
      MPI_Irecv(sf[r].data(), ns[r], MPI_INT, r, 4243, q.comm, &rq[rq.size() - 2]);
      MPI_Isend(rf[r].data(), nr[r], MPI_INT, r, 4243, q.comm, &rq[rq.size() - 1]);
    }
    MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    std::vector<int> nsd(nrd.size());
    MPI_Alltoall(nrd.data(), ng, MPI_INT, nsd.data(), ng, MPI_INT, q.comm);
    for (int r : q.rank) {
      q.rlend.insert(q.rlend.end(), nrd.begin() + r*ng, nrd.begin() + (r + 1)*ng);
      q.slend.insert(q.slend.end(), nsd.begin() + r*ng, nsd.begin() + (r + 1)*ng);
      q.soff.push_back(static_cast<int>(sc.size()));
      q.slen.push_back(ns[r]);
      sc.insert(sc.end(), ss[r].begin(), ss[r].end());
      q.hsf.insert(q.hsf.end(), sf[r].begin(), sf[r].end());
      q.roff.push_back(static_cast<int>(rc.size()));
      q.rlen.push_back(nr[r]);
      rc.insert(rc.end(), rg[r].begin(), rg[r].end());
      q.hrf.insert(q.hrf.end(), rf[r].begin(), rf[r].end());
    }
#endif
    auto mk = [](const char *nm, const std::vector<int> &v) {
      DvceArray1D<int> a(std::string(nm), std::max(v.size(), static_cast<size_t>(1)));
      Kokkos::deep_copy(Kokkos::subview(a, std::make_pair(static_cast<size_t>(0),
                                                          v.size())),
                        HostArray1D<int>(const_cast<int *>(v.data()), v.size()));
      return a;
    };
    q.nl = static_cast<int>(ld.size());
    q.ns = static_cast<int>(sc.size());
    q.nr = static_cast<int>(rc.size());
    q.ld = mk("cln_hdl_ld", ld);
    q.ls = mk("cln_hdl_ls", ls);
    q.sc = mk("cln_hdl_sc", sc);
    q.rc = mk("cln_hdl_rc", rc);
    q.rd = mk("cln_hdl_rd", lrd);   // (rg in rank order: rc's order)
    q.ldp = mk("cln_hdl_ldp", ldd);
    q.lf = mk("cln_hdl_lf", lfm);
    q.sf = mk("cln_hdl_sf", q.hsf);
    std::vector<int> sdv(sc.size(), ng);   // (each message's entries in order of depth)
    for (size_t s = 0; s < q.soff.size(); ++s) {
      for (int e = 0, d = 1; e < q.slen[s]; ++e) {
        while (d < ng && e >= q.slend[s*ng + d - 1]) ++d;
        sdv[q.soff[s] + e] = d;
      }
    }
    q.sdp = mk("cln_hdl_sdp", sdv);
    q.rf = mk("cln_hdl_rf", q.hrf);
    q.sbuf = decltype(q.sbuf)("cln_hdl_sbuf", std::max(q.ns, 1)*cln_hmax);
    q.rbuf = decltype(q.rbuf)("cln_hdl_rbuf", 2*std::max(q.nr, 1)*cln_hmax);
  }
  int bad = 0;
  Kokkos::parallel_reduce("cln_hd_chk", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          static_cast<int>(x.size())),
  KOKKOS_LAMBDA(const int p, int &b) {
    b += (x.data()[p] != y.data()[p]) ? 1 : 0;
  }, bad);
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (bad != 0) {
    if (me == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "CleanHaloInit: the self-check of clean_halo_direct "
                << "failed: its cell lists give " << bad << " cells other than the "
                << "MeshBoundaryValuesCC exchange (array " << n3 << " x " << n2 << " x "
                << n1 << ", " << global_variable::nranks << " ranks)" << std::endl;
    }
    std::exit(EXIT_FAILURE);
  }
  q.state = 1;
  // (the lists replace cln_hpb's buffers from here on)
  Kokkos::fence();
  for (int b = 0; b < 56; ++b) {
    cln_hpb->sendbuf[b].vars = cln_hpb->recvbuf[b].vars = DvceArray2D<Real>();
    cln_hpb->sendbuf[b].flux = cln_hpb->recvbuf[b].flux = DvceArray2D<Real>();
  }
  if (me == 0) {
    std::cout << "<rad_vet> clean_halo_direct: ON, rank 0 copies " << q.nl
              << " ghosts, sends " << q.ns << " cells, receives " << q.nr << " in "
              << q.rank.size() << " messages; on-rank ghosts by depth";
    for (int d = 0; d < ng; ++d) std::cout << " " << q.nld[d];
    std::cout << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::CleanHaloKey, CleanHaloMask
//! \brief clean_halo_direct: the masks of an exchange (CleanGSeg::fm per family, then
//! nc) and their packed buffer layout

template <typename T>
std::array<uint32_t, 5> RadiationVet::CleanHaloKey(const std::vector<CleanGSeg<T>> &sg) {
  std::array<uint32_t, 5> k = {0, 0, 0, 0, 0};
  for (auto &s : sg) {
    const uint32_t a = (s.nc < 32) ? (1u << s.nc) - 1 : ~0u;
    for (int f = 0; f < 4; ++f) k[f] |= (s.fm[f] & a) << k[4];
    if (((s.fm[0] | s.fm[1] | s.fm[2] | s.fm[3]) & a) != a) {   // (no silent skip)
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "CleanHaloKey: a ghost fill whose masks leave a "
                << "component out of every face family" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    k[4] += s.nc;
  }
  return k;
}

RadiationVet::CleanHaloList::Mask &RadiationVet::CleanHaloMask(
    const std::array<uint32_t, 5> &key) {
  CleanHaloList &q = cln_hdl;
  auto &h = q.mask[key];
  if (h.so.extent(0) > 0) return h;
  const int nc = static_cast<int>(key[4]), nseg = static_cast<int>(q.rank.size());
  auto lay = [&](const std::vector<int> &fv, const std::vector<int> &off,
                 const std::vector<int> &len, const std::vector<int> &lend,
                 std::vector<int> &lm, const char *nm) {
    std::vector<int> o(std::max(fv.size(), static_cast<size_t>(1)), 0), c;
    lm.assign(static_cast<size_t>(nseg)*q.ng, 0);
    for (int s = 0; s < nseg; ++s) {
      c.assign(len[s] + 1, 0);
      for (int e = 0; e < len[s]; ++e) {
        o[off[s] + e] = off[s]*nc + c[e];
        c[e+1] = c[e] + static_cast<int>(std::bitset<32>(key[fv[off[s] + e]]).count());
      }
      for (int d = 0; d < q.ng; ++d) lm[s*q.ng + d] = c[lend[s*q.ng + d]];
    }
    DvceArray1D<int> a(std::string(nm), o.size());
    Kokkos::deep_copy(a, HostArray1D<int>(o.data(), o.size()));
    return a;
  };
  h.so = lay(q.hsf, q.soff, q.slen, q.slend, h.slm, "cln_hdl_so");
  h.ro = lay(q.hrf, q.roff, q.rlen, q.rlend, h.rlm, "cln_hdl_ro");
  return h;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanHaloPost
//! \brief clean_halo_direct, first half: post receives (alternate buffer halves), pack,
//! fence, on-rank copy, one send per neighbour rank.  Entry p, component n at p nc + n

template <typename T>
void RadiationVet::CleanHaloPost(const std::vector<CleanGSeg<T>> &sg) {
  CleanHaloList &q = cln_hdl;
  int nc = 0;
  for (auto &s : sg) nc += s.nc;
  const int n1 = q.n1, n21 = q.n2*q.n1, n321 = q.n3*n21;
  q.par = 1 - q.par;
  // (fm: a ghost takes component b if bit b of its family's mask is set, msk)
  const std::array<uint32_t, 5> ky = CleanHaloKey(sg);
  const uint32_t k0 = ky[0], k1 = ky[1], k2 = ky[2], k3 = ky[3];
#if MPI_PARALLEL_ENABLED
  const uint32_t all = (nc < 32) ? (1u << nc) - 1 : ~0u;
  const bool msk = (k0 & k1 & k2 & k3) != all;
  CleanHaloList::Mask *mh = msk ? &CleanHaloMask(ky) : nullptr;
  const int nseg = static_cast<int>(q.rank.size());
  q.rq.assign(nseg, MPI_REQUEST_NULL);
  q.sq.assign(nseg, MPI_REQUEST_NULL);
  Real *rb = q.rbuf.data() + static_cast<size_t>(q.par)*std::max(q.nr, 1)*cln_hmax;
  for (int s = 0; s < nseg; ++s) {
    // (the entries of depth <= dcur)
    const int nrs = msk ? mh->rlm[s*q.ng + q.dcur - 1] : q.rlend[s*q.ng + q.dcur - 1]*nc;
    if (nrs > 0) {
      MPI_Irecv(rb + static_cast<size_t>(q.roff[s])*nc, nrs, MPI_ATHENA_REAL,
                q.rank[s], 4245, q.comm, &q.rq[s]);
    }
  }
  auto sb = q.sbuf;
  auto so = msk ? mh->so : q.sc, sf = q.sf, sdp = q.sdp;
  const int dsd = q.dcur;   // (entries deeper than the exchange are not sent: not packed)
  int o = 0;
  for (auto &s : sg) {
    auto a = s.a;
    auto e = q.sc;
    const int c0 = s.c0;
    if (q.ns > 0) {
      par_for("cln_hd_pack", DevExeSpace(), 0, q.ns-1, 0, s.nc-1,
      KOKKOS_LAMBDA(const int p, const int n) {
        if (sdp(p) > dsd) return;
        const int f = e(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
        int b = p*nc + o + n;
        if (msk) {
          const int fp = sf(p);
          const uint32_t w = (fp == 0) ? k0 : ((fp == 1) ? k1 : ((fp == 2) ? k2 : k3));
          if (!((w >> (o + n)) & 1u)) return;
          b = so(p);
          for (int t = 0; t < o + n; ++t) b += (w >> t) & 1u;
        }
        sb(b) = a(m,c0+n,k,j,r - k*n21 - j*n1);
      });
    }
    o += s.nc;
  }
  if (q.ns > 0 || q.nr > 0) DevExeSpace().fence();
#endif
  // (on-rank copy runs while the sends are posted; ghosts deeper than dcur kept)
  const int dc = q.dcur;
  int oc = 0;
  for (auto &s : sg) {
    auto a = s.a;
    auto d = q.ld, c = q.ls, lf = q.lf;
    const int c0 = s.c0, o = oc;
    oc += s.nc;
    if (q.nl == 0) break;
    par_for("cln_hd_copy", DevExeSpace(), 0, q.nld[dc-1]-1, 0, s.nc-1,
    KOKKOS_LAMBDA(const int p, const int n) {
      const int f = d(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
      const int fp = lf(p);
      const uint32_t w = (fp == 0) ? k0 : ((fp == 1) ? k1 : ((fp == 2) ? k2 : k3));
      if (!((w >> (o + n)) & 1u)) return;
      const int g = c(p), mg = g/n321, t = g - mg*n321, kg = t/n21, jg = (t - kg*n21)/n1;
      a(m,c0+n,k,j,r - k*n21 - j*n1) = a(mg,c0+n,kg,jg,t - kg*n21 - jg*n1);
    });
  }
#if MPI_PARALLEL_ENABLED
  for (int s = 0; s < nseg; ++s) {
    const int nss = msk ? mh->slm[s*q.ng + q.dcur - 1] : q.slend[s*q.ng + q.dcur - 1]*nc;
    if (nss > 0) {
      MPI_Isend(sb.data() + static_cast<size_t>(q.soff[s])*nc, nss,
                MPI_ATHENA_REAL, q.rank[s], 4245, q.comm, &q.sq[s]);
    }
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanHaloFinish
//! \brief clean_halo_direct, second half: wait, unpack, wait for the sends

template <typename T>
void RadiationVet::CleanHaloFinish(const std::vector<CleanGSeg<T>> &sg) {
#if MPI_PARALLEL_ENABLED
  CleanHaloList &q = cln_hdl;
  int nc = 0;
  for (auto &s : sg) nc += s.nc;
  const int n1 = q.n1, n21 = q.n2*q.n1, n321 = q.n3*n21;
  const int nseg = static_cast<int>(q.rank.size());
  MPI_Waitall(nseg, q.rq.data(), MPI_STATUSES_IGNORE);
  const std::array<uint32_t, 5> ky = CleanHaloKey(sg);
  const uint32_t all = (nc < 32) ? (1u << nc) - 1 : ~0u;
  const bool msk = (ky[0] & ky[1] & ky[2] & ky[3]) != all;
  const uint32_t k0 = ky[0], k1 = ky[1], k2 = ky[2], k3 = ky[3];
  auto ro = msk ? CleanHaloMask(ky).ro : q.rc, rf = q.rf;
  auto rb = q.rbuf;
  const size_t b = static_cast<size_t>(q.par)*std::max(q.nr, 1)*cln_hmax;
  const int dc = q.dcur;
  int o = 0;
  for (auto &s : sg) {
    auto a = s.a;
    auto e = q.rc, dp = q.rd;
    const int c0 = s.c0;
    if (q.nr > 0) {
      par_for("cln_hd_unpack", DevExeSpace(), 0, q.nr-1, 0, s.nc-1,
      KOKKOS_LAMBDA(const int p, const int n) {
        const int f = e(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
        int u = p*nc + o + n;
        const int fp = rf(p);
        const uint32_t w = (fp == 0) ? k0 : ((fp == 1) ? k1 : ((fp == 2) ? k2 : k3));
        if (dp(p) > dc || !((w >> (o + n)) & 1u)) return;   // (not sent)
        if (msk) {
          u = ro(p);
          for (int t = 0; t < o + n; ++t) u += (w >> t) & 1u;
        }
        a(m,c0+n,k,j,r - k*n21 - j*n1) = static_cast<T>(rb(b + u));
      });
    }
    o += s.nc;
  }
  MPI_Waitall(nseg, q.sq.data(), MPI_STATUSES_IGNORE);
#else
  (void)sg;
#endif
}

template <typename T>
void RadiationVet::CleanGhostPhys(DvceArray5D<T> a, const int c0, const int nc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ng = indcs.ng;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = a.extent_int(4), n2 = a.extent_int(3), n3 = a.extent_int(2);
  auto bc = pmy_pack->pmb->mb_bcs.d_view;
  auto pf = [] KOKKOS_FUNCTION (const BoundaryFlag f) {
    return f != BoundaryFlag::block && f != BoundaryFlag::periodic &&
           f != BoundaryFlag::panel;
  };
  // (the cubed sphere's D is Cartesian: a radial wall mirrors it about n, not done here)
  const bool ts = (static_cast<const void *>(a.data()) ==
                   static_cast<const void *>(cln_s.data())) &&
                  !pmy_pack->pmesh->use_cubed_sphere;
  const int ib[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  int rf = 0;   // (bit s: side s reflects)
  for (int s = 0; s < 6; ++s) rf |= (ib[s] == VET_IBC_REFLECT) ? (1 << s) : 0;
  // (x1min without vet_bc_reflect = mirror: the rays do not mirror there, so the D_12
  // ghost stays even (no spurious shear); D_13 stays odd (no wall torque: L_z exact))
  const bool rm = vet_bc_mirror;
  auto sg = [=] KOKKOS_FUNCTION (const int c, const int s, const int m) {
    const int p = c - CLN_D11 - 3;   // (D12, D13, D23: axes 01, 02, 12; p != 2 - d)
    if (!ts || p < 0 || p > 2) return 1.0;
    if (bc(m,s) == BoundaryFlag::polar) return (p < 2) ? -1.0 : 1.0;
    if (s == 0 && !rm && p == 0) return 1.0;
    return (((rf >> s) & 1) && p != 2 - s/2) ? -1.0 : 1.0;
  };
  // (no launch for a direction without a physical side on this rank)
  auto bh = pmy_pack->pmb->mb_bcs.h_view;
  auto phys = [&](const int d) {
    for (int m = 0; m <= nmb1; ++m)
      if (pf(bh(m,2*d)) || pf(bh(m,2*d+1))) return true;
    return false;
  };
  // (x1, x2, x3 in turn: a corner takes the product of its sides' parities)
  if (phys(0)) {
    par_for("cln_ghost1", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, n2-1, 0, ng-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int g) {
      if (pf(bc(m,BoundaryFace::inner_x1))) {
        a(m,c0+n,k,j,is-1-g) = sg(c0+n, 0, m)*a(m,c0+n,k,j,is);
      }
      if (pf(bc(m,BoundaryFace::outer_x1))) {
        a(m,c0+n,k,j,ie+1+g) = sg(c0+n, 1, m)*a(m,c0+n,k,j,ie);
      }
    });
  }
  if (pmy_pack->pmesh->multi_d && phys(1)) {
    par_for("cln_ghost2", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, ng-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int g, const int i) {
      if (pf(bc(m,BoundaryFace::inner_x2))) {
        a(m,c0+n,k,js-1-g,i) = sg(c0+n, 2, m)*a(m,c0+n,k,js,i);
      }
      if (pf(bc(m,BoundaryFace::outer_x2))) {
        a(m,c0+n,k,je+1+g,i) = sg(c0+n, 3, m)*a(m,c0+n,k,je,i);
      }
    });
  }
  if (pmy_pack->pmesh->three_d && phys(2)) {
    par_for("cln_ghost3", DevExeSpace(), 0, nmb1, 0, nc-1, 0, ng-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int g, const int j, const int i) {
      if (pf(bc(m,BoundaryFace::inner_x3))) {
        a(m,c0+n,ks-1-g,j,i) = sg(c0+n, 4, m)*a(m,c0+n,ks,j,i);
      }
      if (pf(bc(m,BoundaryFace::outer_x3))) {
        a(m,c0+n,ke+1+g,j,i) = sg(c0+n, 5, m)*a(m,c0+n,ke,j,i);
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::CleanReduce
//! \brief maximum (max) or sum of v over all ranks

Real RadiationVet::CleanReduce(Real v, const bool max) {
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_ATHENA_REAL, max ? MPI_MAX : MPI_SUM,
                MPI_COMM_WORLD);
#endif
  return v;
}

//----------------------------------------------------------------------------------------
//! \fn std::array<Real, 2> CleanSums
//! \brief the n <= 2 sums f(m, k, j, i, s0, s1) adds over the active cells and ranks
//! (f may also update vectors); layout_invariant: via VetLayReduce

template <typename F>
std::array<Real, 2> CleanSums(MeshBlockPack *pp, const char *name, const int n,
                              const F &f) {
  auto &indcs = pp->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pp->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  std::array<Real, 2> d = {0.0, 0.0};
  if (pp->pmesh->layout_invariant) {
    VetBcgVal v;
    auto *r = pp->pradvet;
    const int nxb[3] = {nx1, nx2, nx3};   // (layout_invariant_mode = tile)
    VetLayReduce(name, nmb*nkji, nkji, r->lay_part, r->lay_part_h,
    KOKKOS_LAMBDA(const int idx, VetBcgVal &a) {
      int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/nx1;
      f(m, k+ks, j+js, idx - m*nkji - k*nji - j*nx1 + is, a.s0, a.s1);
    }, v, nxb);
    return {v.s0, v.s1};
  }
  Kokkos::parallel_reduce(name, Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &s0, Real &s1) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    f(m, k+ks, j+js, i+is, s0, s1);
  }, Kokkos::Sum<Real>(d[0]), Kokkos::Sum<Real>(d[1]));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, d.data(), n, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanSumExact
//! \brief CleanSums' s0 by the exact accumulator (utils/repro_sum.hpp): no round-off of
//! the sum, the same for every layout and rank count

template <typename F>
Real CleanSumExact(MeshBlockPack *pp, const char *name, const F &f) {
  auto &indcs = pp->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nji = indcs.nx2*nx1, nkji = indcs.nx3*nji;
  repro::ReproVal<1> v;
  Kokkos::parallel_reduce(name, Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int idx, repro::ReproVal<1> &a) {
    const int m = idx/nkji, k = (idx%nkji)/nji, j = (idx%nji)/nx1;
    Real s0 = 0.0, s1 = 0.0;
    f(m, k+ks, j+js, idx%nx1 + is, s0, s1);
    a.s[0].Add(s0);
  }, repro::ReproRed<1, Kokkos::HostSpace>(v));
  repro::AllMerge(v.s, 1);
  return v.s[0].Value();
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::CleanSumCells
//! \brief the sum of a over the active cells (CleanSums)

Real RadiationVet::CleanSumCells(const DvceArray4D<Real> &a) {
  auto a_ = a;
  return CleanSums(pmy_pack, "cln_sumc", 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0, Real &s1) {
    s0 += a_(m,k,j,i);
  })[0];
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::CleanLzTotal
//! \brief history L_z: sum V R (p_phi + F_phi/c^2) of gas and radiation plus the L_z
//! ledger (CleanLzBook); 0 before the first solve; clean_cs_lz_tel: face momenta

Real RadiationVet::CleanLzTotal() {
  if (!cln_init || !(sph_geom || cs_geom || cln_lzsp)) return 0.0;
  CleanArgs x = CleanMakeArgs();
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool fo = fl.on;
  auto uh = fl.u0;
  auto ur = u0;
  auto f2 = f0x2, f3 = f0x3;
  const bool tel = cln_lzt;
  const Real ic2 = 1.0/(x.c*x.c);
  const Real s = CleanSumExact(pmy_pack, "cln_lztot",
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0, Real &s1) {
    Real p[3] = {0.0, 0.0, 0.0};
    if (fo) CleanGasMomentum(x, uh, m, k, j, i, p);
    const Real r = CleanLzArm(x, m, k, j, i);
    const Real v = CleanVol(x, m, k, j, i);
    if (!tel) {   // (Cartesian: x p_y - y p_x)
      s0 += (x.sph || x.cub) ? v*r*(p[2] + ur(m,VET_F3,k,j,i)*ic2) :
            v*(x.x1v(m,i)*(p[1] + ur(m,VET_F2,k,j,i)*ic2) -
               x.x2v(m,j)*(p[0] + ur(m,VET_F1,k,j,i)*ic2));
      return;
    }
    Real lr = 0.0;   // (the cell's halves of its lateral faces' V F lambda/c^2)
    for (int qd = 1; qd < 3; ++qd) {
      for (int sd = 0; sd < 2; ++sd) {
        const int kf = k + sd*(qd == 2), jf = j + sd*(qd == 1);
        lr += 0.5*v*((qd == 1) ? f2(m,kf,jf,i) : f3(m,kf,jf,i))*ic2*
              CleanCsArm(x, qd, m, kf, jf, i);
      }
    }
    s0 += v*r*p[2] + lr;
  });
  return s + (cln_lzbook ? cln_lzl + cln_lzp : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn std::array<Real, 2> RadiationVet::CleanDots
//! \brief the n <= 2 dot products of w(ab[2q]) and w(ab[2q+1]) over nc components, one
//! pass and one reduction

std::array<Real, 2> RadiationVet::CleanDots(CleanW w, const int nc, const int n,
                                            const std::array<int, 4> ab) {
  const int a0 = ab[0], a1 = ab[1], a2 = ab[2], a3 = ab[3];   // (std::array is host-only)
  return CleanSums(pmy_pack, "cln_dot", n,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0, Real &s1) {
    for (int c = 0; c < nc; ++c) {
      s0 += w(m,a0+c,k,j,i)*w(m,a1+c,k,j,i);
      if (n > 1) s1 += w(m,a2+c,k,j,i)*w(m,a3+c,k,j,i);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanFaceCoefficients
//! \brief per face: blend weight w(C1), theta, old-vector flux G and h+_f = H+_f/J_L,
//! h-_f = H-_f/J_R (H+- minmod-reconstructed from the upwind side; cubed sphere: about
//! the face normal).  |h+-_f| <= 1 keeps the operator an M-matrix

void RadiationVet::CleanFaceCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real h = dt_sub, c = c_light;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  const bool stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  auto cs = cln_s;
  auto op = opac;
  auto inc = t2inc;
  const bool fsk = cln_fsolve;
  // (cubed sphere: seam ghost opacities from the neighbour's cells: single-valued faces)
  if (pmy_pack->pmesh->use_cubed_sphere) {
    CleanGhostFill(opac, 0, opac.extent_int(1), false);
  }
  auto gn = cln_gn;
  auto fi = cln_fvi;
  auto dv = cln_fvd;
  auto fh = cln_fh;
  auto x1f = pmy_pack->pcoord->xx1f;
  auto mbc = pmy_pack->pmb->mb_bcs.d_view;
  const bool hrm = cln_hface;
  const bool fsc = cln_cssc;   // (the sc rays: their folded sums, CleanFaceHalf)
  const bool wfs = cln_wfse;   // (clean_w_fse)
  const bool wft = cln_wfst;   // (clean_w_fse_tau)
  const bool fjs = cln_jse;    // (clean_face_jse)
  const bool fjc = cln_jsec;   // (clean_face_jse_cell)
  const bool fho = cln_rfh;    // (vet_rotate_avg_ema_wb: the cells' own)
  const CleanArgs *px = (fho || fjs || wft) ? CleanDevArgs(CleanMakeArgs(), 0) : nullptr;
  // (clean_closure = eddington | sphere sweeps no rays: the cells' half ranges)
  const int nd = (pmy_pack->pmesh->use_cubed_sphere && !cln_edd && !cln_sphere)
                 ? cln_fvn : 0;

  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = VET_T2_F1 + dir;

    par_for("cln_facecoef", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      const Real chi = CleanFaceChi(op, x1f, mbc, hrm, dir, m, k, j, i, is, ie);
      // (clean_w_fse: g = max(|H|/J, eps C1/1.5) of the cell, else its C1)
      auto gw = [&](const int kk, const int jj, const int ii) {
        if (!wfs) return cs(m,CLN_C1,kk,jj,ii);
        Real h2 = 0.0;
        for (int a = 0; a < 3; ++a) {
          h2 += SQR(cs(m,CLN_HP+a,kk,jj,ii) + cs(m,CLN_HM+a,kk,jj,ii));
        }
        const Real ct = op(m,VET_OP_T,kk,jj,ii);
        Real ep = (ct > 0.0) ? fmin(op(m,VET_OP_P,kk,jj,ii)/ct, 1.0) : 1.0;
        if (wft && ep < 1.0) {   // (clean_w_fse_tau: tau_d = chi V/A_f, the faces' mean)
          Real td = 2.0*ct/(CleanDivW(*px, dir, m, kk, jj, ii, 0) +
                            CleanDivW(*px, dir, m, kk, jj, ii, 1));
          if (px->cub && dir > 0) {   // (the cell's own family: one w at a seam)
            Real s;
            const int q = CleanCsSlot(*px, dir, m, k, j, i, kk, jj, ii, false, s)/2;
            td = ct*px->g(m,CLG_F+3,kk,jj,ii)/px->g.At<CLG_A>(m,q,kk,jj,ii);
          }
          ep = 1.0 - (1.0 - ep)*(1.0 - exp(-td));
        }
        return fmax(sqrt(h2), ep*cs(m,CLN_C1,kk,jj,ii)/1.5);
      };
      const Real c1 = fmax(gw(k, j, i), gw(kl, jl, il));
      fc(m,CLF_W,k,j,i) = CleanSmoothStep((c1 - c1lo)/(c1hi - c1lo));
      fc(m,CLF_TH,k,j,i) = 1.0/(1.0 + h*c*chi);
      fc(m,CLF_G,k,j,i) = f0(m,k,j,i) + (stage ? inc(m,tf,k,j,i) : 0.0);
      if (fsk) fc(m,CLF_N+2,k,j,i) = fc(m,CLF_G,k,j,i);   // (fsolve: start of its F)

      // rays' half-range fluxes H+ = h+ J (from the left) and H- = h- J (from the right)
      auto hp = [&](const int kk, const int jj, const int ii) {
        return cs(m,CLN_HP+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);
      };
      auto hm = [&](const int kk, const int jj, const int ii) {
        return cs(m,CLN_HM+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);
      };
      auto minmod = [](const Real a, const Real b) {
        return (a*b <= 0.0) ? 0.0 : ((fabs(a) < fabs(b)) ? a : b);
      };
      const Real hpl = hp(kl,jl,il);
      const Real hpf = hpl + 0.5*minmod(hpl - hp(kl-dk,jl-dj,il-di), hp(k,j,i) - hpl);
      const Real hmr = hm(k,j,i);
      const Real hmf = hmr - 0.5*minmod(hmr - hm(kl,jl,il), hm(k+dk,j+dj,i+di) - hmr);
      const Real jl_ = cs(m,CLN_J,kl,jl,il), jr_ = cs(m,CLN_J,k,j,i);
      fc(m,CLF_HP,k,j,i) = (jl_ > 0.0) ? fmin(fmax(hpf/jl_, 0.0), 1.0)
                                       : cs(m,CLN_HP+dir,kl,jl,il);
      fc(m,CLF_HM,k,j,i) = (jr_ > 0.0) ? fmax(fmin(hmf/jr_, 0.0), -1.0)
                                       : cs(m,CLN_HM+dir,k,j,i);
      // (clean_face_jse: scale h+ by H+_f/H+_L, H+_f = wl H+_L + wr H+_R (CleanHcTauW),
      // factor in [1/2, 2], blended by eps = min(kappa_P/chi, 1) of the face (h- alike);
      // radial faces on curved meshes, all on Cartesian.  clean_face_jse_cell: scale the
      // donor cell's h+- instead of the face value)
      auto jse = [&]() {
        if (!fjs || !(dir == 0 || !(px->sph || px->cub))) return;
        auto ep = [&](const int kk, const int jj, const int ii) {
          const Real ct = op(m,VET_OP_T,kk,jj,ii);
          return (ct > 0.0) ? fmin(op(m,VET_OP_P,kk,jj,ii)/ct, 1.0) : 1.0;
        };
        const Real wt = fmin(ep(kl,jl,il), ep(k,j,i));
        Real wl, wr;
        CleanHcTauW(*px, dir, m, k, j, i, wl, wr);
        const Real ppl = hp(kl,jl,il), ppr = fmax(hp(k,j,i), 0.0);
        const Real mml = fmin(hm(kl,jl,il), 0.0), mmr = hm(k,j,i);
        auto &fp = fc(m,CLF_HP,k,j,i), &fm = fc(m,CLF_HM,k,j,i);
        if (ppl > 0.0 && fp > 0.0) {
          const Real ph = fmin(fmax(wl + wr*ppr/ppl, 0.5), 2.0);
          fp += wt*(fmin((fjc ? cs(m,CLN_HP+dir,kl,jl,il) : fp)*ph, 1.0) - fp);
        }
        if (mmr < 0.0 && fm < 0.0) {
          const Real ph = fmin(fmax(wr + wl*mml/mmr, 0.5), 2.0);
          fm += wt*(fmax((fjc ? cs(m,CLN_HM+dir,k,j,i) : fm)*ph, -1.0) - fm);
        }
      };
      if (nd == 0) {
        jse();
        return;
      }
      // (cubed sphere: about the face normal from the FV intensities of L and R, so a
      // seam face takes the same h+- on both sides)
      // (J_L, J_R, H+_L, H-_R as scalars: nvcc -O3 miscompiles a runtime-indexed local)
      if (fho) {   // (each cell's value out through this face, in its own axes)
        Real s;
        const int nl = CleanCsSlot(*px, dir, m, k, j, i, kl, jl, il, false, s);
        const int nr = CleanCsSlot(*px, dir, m, k, j, i, k, j, i, false, s);
        fc(m,CLF_HP,k,j,i) = (nl & 1) ? cs(m,CLN_HO+nl/2,kl,jl,il)
                                      : -cs(m,CLN_HO+3+nl/2,kl,jl,il);
        fc(m,CLF_HM,k,j,i) = (nr & 1) ? -cs(m,CLN_HO+nr/2,k,j,i)
                                      : cs(m,CLN_HO+3+nr/2,k,j,i);
        jse();
        return;
      }
      Real jl0 = 0.0, jr0 = 0.0, hp0 = 0.0, hm0 = 0.0;
      if (fsc) {
        jl0 = fh(m,4*dir,k,j,i);
        hp0 = fh(m,4*dir+2,k,j,i);
        hm0 = fh(m,4*dir+3,k,j,i);
        jr0 = fh(m,4*dir+1,k,j,i);
      }
      for (int d = 0; d < nd && !fsc; ++d) {
        Real mu = 0.0;
        for (int b = 0; b < 3; ++b) mu += dv(d,b)*gn.At<CLG_N>(m,3*dir+b,k,j,i);
        const Real wl = dv(d,3)*fi(m,d,kl,jl,il), wr = dv(d,3)*fi(m,d,k,j,i);
        jl0 += wl;
        jr0 += wr;
        if (mu > 0.0) hp0 += wl*mu;
        if (!(mu > 0.0)) hm0 += wr*mu;
      }
      fc(m,CLF_HP,k,j,i) = (jl0 > 0.0) ? hp0/jl0 : 0.25;
      fc(m,CLF_HM,k,j,i) = (jr0 > 0.0) ? hm0/jr0 : -0.25;
      jse();
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanApply
//! \brief w(dst) = A w(src) over the active cells; src's exchange overlaps the inner
//! cells (cubed sphere: only with the xdonor, else the lateral sums need ghosts first)

void RadiationVet::CleanApply(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const bool ov = !pmy_pack->pmesh->use_cubed_sphere || cln_xd_on;
  auto wa = (src < cln_w.ns) ? cln_w.s : cln_w.o;   // (src's array, CLW_NS)
  const int ws = (src < cln_w.ns) ? src : src - cln_w.ns;
  // (cubed sphere: the stencil and lateral sums read the first ghost layer)
  const bool cub0 = pmy_pack->pmesh->use_cubed_sphere;
  CleanGhostFill(wa, ws, 1, true, ov ? 1 : 0, cub0 ? 1 : 0);
  const std::function<void()> fin = [&]() {if (ov) CleanGhostFill(wa, ws, 1, true, 2);};
  // (cubed sphere: the exchange of src's lateral sums Y overlaps the inner cross terms)
  const bool yov = pmy_pack->pmesh->use_cubed_sphere && !cln_xd_on;
  if (!cln_xd_on) CleanStagCsLateral(src, true, yov ? 1 : 0);   // (xdonor: none)
  CleanArgs x = CleanMakeArgs();
  x.fwc = true;
  // (the operator from its rows, CleanPreconditionerSetup)
  auto w = cln_w;
  auto st = cln_st;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  // (cubed sphere: lateral stress from the sums Y; none with the xdonor)
  const bool cub = x.cub, cxf = cub && !(x.xq >= 0 && !x.hi);
  const CleanArgs *px = cxf ? CleanDevArgs(x, 0) : nullptr;
  // (cubed sphere: each face's cross term once, not per adjacent cell)
  DvceArray5D<Real> xf;
  if (cxf) {
    const int ie = indcs.ie, je = indcs.je, ke = indcs.ke;
    xf = ScratchView<DvceArray5D<Real>>("cln_xf", pmy_pack->nmb_thispack, 3,
                                        w.s.extent_int(2), w.s.extent_int(3),
                                        w.s.extent_int(4));
    // (pass 0: faces between two active cells; pass 1, after the exchange: end planes)
    for (int ps = 0; ps < 2; ++ps) {
      if (ps == 1 && yov) CleanGhostFill(cln_w, CLW_L, 9, 2, 1);
      const bool xon0 = x.multi_d && !(x.xq >= 0 && !x.hi);
      if (ps == 1 || !xon0) {
        // all directions' faces in one launch: pass 1 the end planes, pass 0 (no inner
        // kernel) the inner faces
        const int ne[3] = {ie - is + 1, je - js + 1, ke - ks + 1};
        Kokkos::Array<int, 7> c0;   // (segment s: faces c0[s] .. c0[s+1]-1 per block)
        Kokkos::Array<int, 6> sd, sk0, sj0, si0, snj, sni;
        int ns = 0;
        c0[0] = 0;
        for (int d = 0; d < x.nf; ++d) {
          const int dk = (d == 2), dj = (d == 1), di = (d == 0);
          for (int pl = 0; pl < (ps ? 2 : 1); ++pl) {
            const int o = pl*ne[d];
            const int k0 = ps ? ks + o*dk : ks + dk, k1 = ps ? (dk ? ks + o : ke) : ke;
            const int j0 = ps ? js + o*dj : js + dj, j1 = ps ? (dj ? js + o : je) : je;
            const int i0 = ps ? is + o*di : is + di, i1 = ps ? (di ? is + o : ie) : ie;
            sd[ns] = d;
            sk0[ns] = k0;
            sj0[ns] = j0;
            si0[ns] = i0;
            snj[ns] = j1 - j0 + 1;
            sni[ns] = i1 - i0 + 1;
            c0[ns+1] = c0[ns] + (k1 - k0 + 1)*snj[ns]*sni[ns];
            ++ns;
          }
        }
        const int nper = c0[ns];
        Kokkos::parallel_for("cln_xfe", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                             pmy_pack->nmb_thispack*nper),
        KOKKOS_LAMBDA(const int t) {
          const int m = t/nper;
          int c = t - m*nper, q = 0;
          while (c >= c0[q+1]) {++q;}
          c -= c0[q];
          const int i = si0[q] + c % sni[q];
          c /= sni[q];
          const int j = sj0[q] + c % snj[q];
          const int k = sk0[q] + c/snj[q];
          xf(m,sd[q],k,j,i) = CleanFaceCross(*px, src, sd[q], m, k, j, i);
        });
        continue;
      }
      for (int d = 0; d < x.nf; ++d) {
        const int dk = (d == 2), dj = (d == 1), di = (d == 0);
        auto gn = x.g;
        auto fcd = x.fc[d];
        // (inside: CleanFaceCross with no boundary or seam, its own kernel)
        CleanFor("cln_xf", pmy_pack->nmb_thispack - 1, ks + dk, ke, js + dj, je, is + di,
                 ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          const int kl = k - dk, jl = j - dj, il = i - di;
          auto lat = [&](const int kc, const int jc, const int ic) {
            Real s;   // (CleanStagCsLat)
            const int q = d ? CleanCsSlot(*px, d, m, k, j, i, kc, jc, ic, false, s)/2
                            : 0;
            const int p = (q == 0) ? 1 : 0, pp = (q == 2) ? 1 : 2;
            Real r = 0.0;
            for (int c = 0; c < 3; ++c) {
              r += gn.At<CLG_N>(m,3*d+c,k,j,i)*(w(m,CLW_L+3*p+c,kc,jc,ic) +
                                           w(m,CLW_L+3*pp+c,kc,jc,ic));
            }
            return r;
          };
          Real dv = 0.0;
          dv += lat(kl, jl, il) + lat(k, j, i);
          xf(m,d,k,j,i) = fcd(m,CLF_QL+2,k,j,i)*
              (dv/(0.5*(gn(m,CLG_F+3,k,j,i) + gn(m,CLG_F+3,kl,jl,il))));
        });
      }
    }
  }
  auto f = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real a = 0.0;
#pragma unroll
    for (int s = 0; s < CLN_ST; ++s) {
      int dk, dj, di;
      CleanStOff(s, dk, dj, di);
      if ((dj != 0 && !md) || (dk != 0 && !td) || (cub && dk*dk + dj*dj + di*di > 1)) {
        continue;   // (no such neighbour)
      }
      a += st(m,s,k-ks,j-js,i-is)*w(m,src,k+dk,j+dj,i+di);
    }
    if (cxf) {
      const CleanArgs &x = *px;   // (in device memory: no per-thread copy)
      for (int d = 0; d < 3; ++d) {
        for (int sd = 0; sd < 2 && d < x.nf; ++sd) {
          const int kf = k + sd*(d == 2), jf = j + sd*(d == 1), iff = i + sd*(d == 0);
          a += (sd ? x.h : -x.h)*CleanDivW(x, d, m, k, j, i, sd)*xf(m,d,kf,jf,iff);
        }
      }
    }
    w(m,dst,k,j,i) = a;
  };
  CleanOverlapFor("cln_apply_st", pmy_pack, ov, f, fin);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanPreconditionerSetup
//! \brief preconditioner coefficients: the x1-line tridiagonal part of A plus the
//! diagonal of its x2, x3 parts (cross terms left to Krylov); for the mg also the x2, x3
//! couplings and the coarse levels.  Also the face weights and the operator rows

void RadiationVet::CleanPreconditionerSetup() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto op = opac;
  // the face weights of the pass, once per face for CleanApply (x.fwc)
  CLEAN_HIPX(const CleanArgs *pxw = CleanDevArgs(x, 0);)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    CleanFor<0, 3>("cln_fwc", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *pxw;)   // (in device memory: no per-thread copy)
      CleanFaceWeights(x, dir, m, k, j, i, fc(m,CLF_QL,k,j,i), fc(m,CLF_QL+1,k,j,i));
      const Real t0 = fc(m,CLF_TH,k,j,i);   // (clean_stag_vimp: theta_d, CleanFaceCross)
      const Real th = t0 + (x.drag ? fc(m,CLF_DB,k,j,i) : 0.0)*(1.0 - t0);
      fc(m,CLF_QL+2,k,j,i) = (x.xq >= 0) ? 0.0
                                         : -(1.0 - fc(m,CLF_W,k,j,i))*th*x.h*x.c*x.c;
    });
  }
  x.fwc = true;
  const CleanArgs *px = CleanDevArgs(x, 0);
  // the operator's rows (CleanApply)
  auto st = cln_st;
  CleanFor<0, 3>("cln_stencil", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px;   // (in device memory: no per-thread copy)
    Real c[CLN_ST];
    CleanStencilRow(x, m, k, j, i, c);
    for (int s = 0; s < CLN_ST; ++s) st(m,s,k-ks,j-js,i-is) = c[s];
  });
  // (lateral stress sums' coefficients: geometry and closure only, built once per
  // closure, cln_csyok)
  if (x.cub && !cln_csyok) {
    cln_csyok = true;
    auto cy = cln_csy;
    par_for("cln_csyb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;   // (in device memory: no per-thread copy)
      for (int q = 0; q < 3; ++q) {
        const int dk = (q == 2), dj = (q == 1), di = (q == 0);
        Real y0[3] = {0.0, 0.0, 0.0}, yn[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        const int vl = CleanTransSide(x, q, m, k, j, i, 0);
        const int vh = CleanTransSide(x, q, m, k, j, i, 1);
        // (seam side: ghost extrapolated in panel, CleanCsSeamCoef)
        const bool o2 = (vl + vh == 1) && x.mbc(m,2*q+vl) == BoundaryFlag::panel;
        const Real wq = (vl + vh == 2 || o2) ? 0.25 : 0.5;
        for (int s = -1; s <= 1; s += 2) {
          const int kf = k + (s > 0)*dk, jf = j + (s > 0)*dj, iff = i + (s > 0)*di;
          const bool ex = o2 && !CleanTransSide(x, q, m, k, j, i, s > 0);
          if (!ex && !CleanTransSide(x, q, m, k, j, i, s > 0)) continue;
          const Real a = wq*s*x.area[q](m,kf,jf,iff);
          for (int b = 0; b < 3; ++b) {
            for (int c = 0; c < 3; ++c) {
              const int n = CLN_D11 + ((b == c) ? b : 2 + b + c);
              const Real ag = a*x.g.At<CLG_N>(m,3*q+c,kf,jf,iff);
              if (ex) {   // (c0 P0 + c1 P1 + c2 P2 of the other side's cells)
                Real cf[3];
                CleanCsSeamCoef(x, n, q, m, k, j, i, -s, cf);
                y0[b] += cf[0]*ag*x.cs(m,n,k,j,i);
                yn[s < 0][b] += cf[1]*ag*x.cs(m,n,k-s*dk,j-s*dj,i-s*di);
                yn[s > 0][b] += cf[2]*ag*x.cs(m,n,k-2*s*dk,j-2*s*dj,i-2*s*di);
                continue;
              }
              yn[s > 0][b] += ag*x.cs(m,n,k+s*dk,j+s*dj,i+s*di);
              y0[b] -= ag*x.cs(m,n,k,j,i);
            }
          }
        }
        for (int b = 0; b < 3; ++b) {
          cy(m,6*q+b,k-ks,j-js,i-is) = yn[0][b];
          cy(m,6*q+3+b,k-ks,j-js,i-is) = yn[1][b];
        }
        for (int b = 0; b < 3; ++b) cy(m,18+3*q+b,k-ks,j-js,i-is) = y0[b];
      }
    });
  }
  par_for("cln_pcsetup", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    CLEAN_HIPX(const CleanArgs &x = *px;)   // (in device memory: no per-thread copy)
    Real al_p, ar_p, al_m, ar_m;
    CleanFaceWeights(x, 0, m, k, j, i+1, al_p, ar_p);
    CleanFaceWeights(x, 0, m, k, j, i, al_m, ar_m);
    const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
    const Real pd = x.pd ? x.w(m,CLW_PD,k,j,i) : 0.0;   // (clean_stag_central)
    Real a0 = 1.0 + x.h*x.c*op(m,VET_OP_E,k,j,i)*qr + pd + x.h*x.idx[0]*(al_p - ar_m);
    x.w(m,CLW_AM,k,j,i) = -x.h*x.idx[0]*al_m;
    x.w(m,CLW_AP,k,j,i) = x.h*x.idx[0]*ar_p;
    if (x.sph || x.cub) {   // (A_f/V weights, CleanDivW)
      const Real wp = CleanDivW(x, 0, m, k, j, i, 1), wm = CleanDivW(x, 0, m, k, j, i, 0);
      a0 = 1.0 + x.h*x.c*op(m,VET_OP_E,k,j,i)*qr + pd + x.h*(wp*al_p - wm*ar_m);
      x.w(m,CLW_AM,k,j,i) = -x.h*wm*al_m;
      x.w(m,CLW_AP,k,j,i) = x.h*wp*ar_p;
    }
    if (x.multi_d) {
      CleanFaceWeights(x, 1, m, k, j+1, i, al_p, ar_p);
      CleanFaceWeights(x, 1, m, k, j, i, al_m, ar_m);
      a0 += (x.sph || x.cub) ? x.h*(CleanDivW(x, 1, m, k, j, i, 1)*al_p -
                                    CleanDivW(x, 1, m, k, j, i, 0)*ar_m)
                             : x.h*x.idx[1]*(al_p - ar_m);
      const bool cv = x.sph || x.cub;
      x.w(m,CLW_TM,k,j,i) = -x.h*al_m*(cv ? CleanDivW(x, 1, m, k, j, i, 0) : x.idx[1]);
      x.w(m,CLW_TM+1,k,j,i) = x.h*ar_p*(cv ? CleanDivW(x, 1, m, k, j, i, 1) : x.idx[1]);
    }
    if (x.three_d) {
      CleanFaceWeights(x, 2, m, k+1, j, i, al_p, ar_p);
      CleanFaceWeights(x, 2, m, k, j, i, al_m, ar_m);
      const bool cv = x.sph || x.cub;
      a0 += cv ? x.h*(CleanDivW(x, 2, m, k, j, i, 1)*al_p -
                      CleanDivW(x, 2, m, k, j, i, 0)*ar_m)
               : x.h*x.idx[2]*(al_p - ar_m);
      x.w(m,CLW_TM+2,k,j,i) = -x.h*al_m*(cv ? CleanDivW(x, 2, m, k, j, i, 0) : x.idx[2]);
      x.w(m,CLW_TM+3,k,j,i) = x.h*ar_p*(cv ? CleanDivW(x, 2, m, k, j, i, 1) : x.idx[2]);
    }
    x.w(m,CLW_A0,k,j,i) = a0;
  });
  CleanStagCoarseSetup();
}

//----------------------------------------------------------------------------------------
//! \fn CleanLevel RadiationVet::CleanStagLevel
//! \brief level l of the staggered mg: 0 the mesh (cln_w), l > 0 cln_mg[l-1]; scalar
//! rows in cln_w's slot order (AM, A0, AP, scratch) on every level

CleanLevel RadiationVet::CleanStagLevel(const int l) {
  CleanLevel lv;
  if (l == 0) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    lv.u = cln_w;
    lv.o[0] = indcs.is;
    lv.o[1] = indcs.js;
    lv.o[2] = indcs.ks;
    lv.n[0] = indcs.nx1;
    lv.n[1] = indcs.nx2;
    lv.n[2] = indcs.nx3;
    lv.a0 = CLW_A0;
  } else {
    lv.u = cln_mg[l-1];
    for (int d = 0; d < 3; ++d) {
      lv.o[d] = 0;
      lv.n[d] = cln_mgn[l-1][d];
    }
    lv.a0 = CLS_A0;
  }
  const int nb[3] = {pmy_pack->pmesh->mb_indcs.nx1, pmy_pack->pmesh->mb_indcs.nx2,
                     pmy_pack->pmesh->mb_indcs.nx3};
  for (int d = 0; d < 3; ++d) {   // (layout_tile; default the MeshBlock)
    lv.t[d] = cln_mgt[0] ? cln_mgt[d]*lv.n[d]/nb[d] : lv.n[d];
  }
  const int tm = (l == 0) ? CLW_TM : CLS_TM;
  lv.sm[0] = lv.a0 - 1;
  lv.sp[0] = lv.a0 + 1;
  for (int d = 1; d < 3; ++d) {
    lv.sm[d] = tm + 2*(d-1);
    lv.sp[d] = tm + 2*(d-1) + 1;
  }
  return lv;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagCoarseSetup
//! \brief coarse rows of the staggered mg, Galerkin R A P: P constant over 2 x 2 (x2,x3)
//! aggregates at fixed x1 within a MeshBlock, R = P^T; in-aggregate couplings go onto
//! the diagonal.  Curved meshes: the mesh rows weighted by cell volume

void RadiationVet::CleanStagCoarseSetup() {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nlev = static_cast<int>(cln_mg.size()) + 1;
  const bool curv = pmy_pack->pmesh->use_spherical_polar ||
                    pmy_pack->pmesh->use_cubed_sphere;
  auto vol = pmy_pack->pcoord->volume;
  for (int lc = 1; lc < nlev; ++lc) {
    const CleanLevel f = CleanStagLevel(lc-1), c = CleanStagLevel(lc);
    const int r[3] = {1, f.n[1]/c.n[1], f.n[2]/c.n[2]};
    const bool vw = (lc == 1) && curv;

    par_for("cls_galerkin", DevExeSpace(), 0, nmb1, 0, c.n[2]-1, 0, c.n[1]-1,
            0, c.n[0]-1,
    KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int ic) {
      Real a0 = 0.0, am[3] = {0.0, 0.0, 0.0}, ap[3] = {0.0, 0.0, 0.0};
      for (int a2 = 0; a2 < r[2]; ++a2) {
        for (int a1 = 0; a1 < r[1]; ++a1) {
          const int k = f.o[2] + r[2]*kc + a2, j = f.o[1] + r[1]*jc + a1;
          const int i = f.o[0] + ic;
          const int a[3] = {0, a1, a2};
          const Real v = vw ? vol(m,k,j,i) : 1.0;
          a0 += v*f.u(m,f.a0,k,j,i);
          for (int d = 0; d < 3; ++d) {
            if (f.n[d] == 1) continue;
            const Real lo = v*f.u(m,f.sm[d],k,j,i), hi = v*f.u(m,f.sp[d],k,j,i);
            if (a[d] > 0) {
              a0 += lo;
            } else {
              am[d] += lo;
            }
            if (a[d] < r[d] - 1) {
              a0 += hi;
            } else {
              ap[d] += hi;
            }
          }
        }
      }
      c.u(m,c.a0,kc,jc,ic) = a0;
      for (int d = 0; d < 3; ++d) {
        c.u(m,c.sm[d],kc,jc,ic) = am[d];
        c.u(m,c.sp[d],kc,jc,ic) = ap[d];
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagLineSolve
//! \brief u(dst) = L^{-1} u(src), L the x1-line tridiagonal part, on the lines of colour
//! col = (j + k) mod 2 (-1: all), with the neighbours of colour != col and <= lim on the
//! RHS.  clean_stag_mg_color = mirror: j, k counted from the nearer edge, so the sweep is
//! mirror-symmetric about the x2, x3 midplanes.  GPU: parallel cyclic reduction, one team
//! per line (cln_pcr); host: Thomas

void RadiationVet::CleanStagLineSolve(const CleanLevel &lv, const int src, const int dst,
                                      const int col, const int lim) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const CleanLevel l = lv;
  const int is = l.o[0], ie = l.o[0] + l.n[0] - 1, cp = l.a0 + 2;
  auto w = l.u;
  const int n1 = l.t[1], n2 = l.t[2];   // (positions in the line's tile)
  const bool mir = (cln_stagmg == 2);
  const int c0 = col;
  auto cl = KOKKOS_LAMBDA(int pj, int pk) {   // (the line's colour)
    pj %= n1;
    pk %= n2;
    return (((mir && 2*pj >= n1) ? n1 - 1 - pj : pj) +
            ((mir && 2*pk >= n2) ? n2 - 1 - pk : pk)) % 2;
  };
  auto on = KOKKOS_LAMBDA(const int pj, const int pk) {   // (a neighbour on the rhs)
    const int c = cl(pj, pk);
    return c != c0 && c <= lim;
  };

  if (!std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace> && cln_pcr) {
    // (n2 the lines along x2; t2 the tile's along x3, as the host's Thomas)
    const int nx = l.n[0], n2 = l.n[1], n32 = l.n[2]*l.n[1], t2 = l.t[2];
    int nround = 0;
    while ((1 << nround) < nx) ++nround;
    // (one row per lane: a team of nx rounded to a warp, double-buffered in scratch, one
    // barrier per round; the rows' arithmetic is that of the host's PCR, so bitwise)
    const int ts = std::min(1024, (nx + 31)/32*32);
    // (a colour: one team per line of that colour, from a list built once per level)
    int nt = (nmb1 + 1)*n32;
    DvceArray1D<int> ll;
    if (col >= 0) {
      const std::array<int, 8> key = {nt, n2, n1, l.t[2], mir, 2, c0, l.o[1]};
      auto it = cln_pcrl.find(key);
      if (it == cln_pcrl.end()) {
        std::vector<int> h;
        for (int q = 0; q < nt; ++q) {
          const int pk = (q % n32)/n2;
          if (cl(q % n32 - pk*n2, pk) == c0) h.push_back(q);
        }
        DvceArray1D<int> d("cln_pcrl", h.size());
        auto hm = Kokkos::create_mirror_view(d);
        for (size_t q = 0; q < h.size(); ++q) hm(q) = h[q];
        Kokkos::deep_copy(d, hm);
        it = cln_pcrl.emplace(key, d).first;
      }
      ll = it->second;
      nt = ll.extent_int(0);
    }
    const bool lst = (col >= 0);
    const size_t scr = ScrArray1D<Real>::shmem_size(8*nx);
    Kokkos::TeamPolicy<DevExeSpace> pol(DevExeSpace(), nt, ts);
    const int lv0 = TeamScratchLevel(scr, 0, pol.team_size());
    pol.set_scratch_size(lv0, Kokkos::PerTeam(scr));
    // (the line of team r: m, j, k, and its rows loaded into the team's scratch sw)
    auto ld = KOKKOS_LAMBDA(const TeamMember_t &tm, const ScrArray1D<Real> &sw, int &m,
                            int &k, int &j) {
      const int q = lst ? ll(tm.league_rank()) : tm.league_rank();
      m = q/n32;
      const int pk = (q - m*n32)/n2, pj = q - m*n32 - pk*n2;
      k = l.o[2] + pk;
      j = l.o[1] + pj;
      if (col >= 0 && cl(pj, pk) != c0) return false;   // (team-uniform)
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
        const int i = is + q;
        Real r = w(m,src,k,j,i);
        if (lim >= 0) {
          if (pj%n1 > 0 && on(pj-1, pk)) r -= w(m,l.sm[1],k,j,i)*w(m,dst,k,j-1,i);
          if (pj%n1 < n1 - 1 && on(pj+1, pk)) r -= w(m,l.sp[1],k,j,i)*w(m,dst,k,j+1,i);
          if (pk%t2 > 0 && on(pj, pk-1)) r -= w(m,l.sm[2],k,j,i)*w(m,dst,k-1,j,i);
          if (pk%t2 < t2 - 1 && on(pj, pk+1)) r -= w(m,l.sp[2],k,j,i)*w(m,dst,k+1,j,i);
        }
        sw(q) = (q > 0) ? w(m,l.sm[0],k,j,i) : 0.0;
        sw(nx + q) = w(m,l.a0,k,j,i);
        sw(2*nx + q) = (q < nx - 1) ? w(m,l.sp[0],k,j,i) : 0.0;
        sw(3*nx + q) = r;
      });
      tm.team_barrier();
      return true;
    };
    Kokkos::parallel_for("cln_pcr2", pol, KOKKOS_LAMBDA(TeamMember_t tm) {
      ScrArray1D<Real> sw(tm.team_scratch(lv0), 8*nx);
      int m, k, j;
      if (!ld(tm, sw, m, k, j)) return;
      int o = 0;
      for (int rd = 0, s = 1; rd < nround; ++rd, s *= 2) {
        const int d = 4*nx - o;
        Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
          Real a = 0.0, b = sw(o + nx + q), c = 0.0, r = sw(o + 3*nx + q);
          if (q >= s) {
            const Real f = -sw(o + q)/sw(o + nx + q - s);
            a = f*sw(o + q - s);
            b += f*sw(o + 2*nx + q - s);
            r += f*sw(o + 3*nx + q - s);
          }
          if (q + s < nx) {
            const Real g = -sw(o + 2*nx + q)/sw(o + nx + q + s);
            c = g*sw(o + 2*nx + q + s);
            b += g*sw(o + q + s);
            r += g*sw(o + 3*nx + q + s);
          }
          sw(d + q) = a;
          sw(d + nx + q) = b;
          sw(d + 2*nx + q) = c;
          sw(d + 3*nx + q) = r;
        });
        tm.team_barrier();
        o = d;
      }
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
        w(m,dst,k,j,is+q) = sw(o + 3*nx + q)/sw(o + nx + q);
      });
    });
    return;
  }
  par_for("cln_thomas", DevExeSpace(), 0, nmb1, l.o[2], l.o[2]+l.n[2]-1,
          l.o[1], l.o[1]+l.n[1]-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int pj = j - l.o[1], pk = k - l.o[2];
    if (col >= 0 && cl(pj, pk) != c0) return;
    Real cprev = 0.0, dprev = 0.0;
    for (int i = is; i <= ie; ++i) {
      const Real am = (i > is) ? w(m,l.sm[0],k,j,i) : 0.0;
      const Real ap = (i < ie) ? w(m,l.sp[0],k,j,i) : 0.0;
      Real r = w(m,src,k,j,i);
      if (lim >= 0) {
        if (pj%n1 > 0 && on(pj-1, pk)) r -= w(m,l.sm[1],k,j,i)*w(m,dst,k,j-1,i);
        if (pj%n1 < n1 - 1 && on(pj+1, pk)) r -= w(m,l.sp[1],k,j,i)*w(m,dst,k,j+1,i);
        if (pk%n2 > 0 && on(pj, pk-1)) r -= w(m,l.sm[2],k,j,i)*w(m,dst,k-1,j,i);
        if (pk%n2 < n2 - 1 && on(pj, pk+1)) r -= w(m,l.sp[2],k,j,i)*w(m,dst,k+1,j,i);
      }
      const Real den = w(m,l.a0,k,j,i) - am*cprev;
      cprev = ap/den;
      dprev = (r - am*dprev)/den;
      w(m,cp,k,j,i) = cprev;
      w(m,dst,k,j,i) = dprev;
    }
    for (int i = ie-1; i >= is; --i) {
      w(m,dst,k,j,i) -= w(m,cp,k,j,i)*w(m,dst,k,j,i+1);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagVCycle
//! \brief u(dst) = M^{-1} u(src) on level l of the staggered mg: z0 = red-black x1-line
//! sweep, r1 = R (r - A z0) (volume-weighted on curved meshes), z = z0 + P z1 with z1
//! from level l+1; the coarsest level only sweeps.  clean_pc_sweep: &1 a second even
//! pass, &2 a post-sweep (odd, then even)

void RadiationVet::CleanStagVCycle(const int l, const int src, const int dst) {
  const CleanLevel f = CleanStagLevel(l);
  CleanStagLineSolve(f, src, dst, 0, -1);
  CleanStagLineSolve(f, src, dst, 1, 1);
  if (cln_pcsweep & 1) CleanStagLineSolve(f, src, dst, 0, 1);
  if (l == static_cast<int>(cln_mg.size())) return;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const CleanLevel c = CleanStagLevel(l+1);
  const int r[3] = {1, f.n[1]/c.n[1], f.n[2]/c.n[2]};
  const bool vw = (l == 0) && (pmy_pack->pmesh->use_spherical_polar ||
                                pmy_pack->pmesh->use_cubed_sphere);
  auto vol = pmy_pack->pcoord->volume;

  par_for("cls_restrict", DevExeSpace(), 0, nmb1, 0, c.n[2]-1, 0, c.n[1]-1, 0, c.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int ic) {
    Real s = 0.0;
    for (int a2 = 0; a2 < r[2]; ++a2) {
      for (int a1 = 0; a1 < r[1]; ++a1) {
        const int k = f.o[2] + r[2]*kc + a2, j = f.o[1] + r[1]*jc + a1;
        const int i = f.o[0] + ic;
        const int pos[3] = {ic, (j - f.o[1])%f.t[1], (k - f.o[2])%f.t[2]};
        Real y = f.u(m,f.a0,k,j,i)*f.u(m,dst,k,j,i);
        for (int d = 0; d < 3; ++d) {
          const int dk = (d == 2), dj = (d == 1), di = (d == 0);
          if (pos[d] > 0) y += f.u(m,f.sm[d],k,j,i)*f.u(m,dst,k-dk,j-dj,i-di);
          if (pos[d] < f.t[d] - 1) y += f.u(m,f.sp[d],k,j,i)*f.u(m,dst,k+dk,j+dj,i+di);
        }
        s += (vw ? vol(m,k,j,i) : 1.0)*(f.u(m,src,k,j,i) - y);
      }
    }
    c.u(m,CLS_B,kc,jc,ic) = s;
  });
  CleanStagVCycle(l+1, CLS_B, CLS_X);
  par_for("cls_prolong", DevExeSpace(), 0, nmb1, f.o[2], f.o[2]+f.n[2]-1,
          f.o[1], f.o[1]+f.n[1]-1, f.o[0], f.o[0]+f.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    f.u(m,dst,k,j,i) += c.u(m,CLS_X,(k - f.o[2])/r[2],(j - f.o[1])/r[1],i - f.o[0]);
  });
  if (cln_pcsweep & 2) {
    CleanStagLineSolve(f, src, dst, 1, 1);
    CleanStagLineSolve(f, src, dst, 0, 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanPrecondition
//! \brief w(dst) = M^{-1} w(src) (CleanStagVCycle)

void RadiationVet::CleanPrecondition(const int src, const int dst) {
  CleanStagVCycle(0, src, dst);
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::CleanBiCGStab
//! \brief right-preconditioned BiCGStab for A X = B from the given X (cln_w); returns
//! the iteration count, negative if not converged

int RadiationVet::CleanBiCGStab() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanW w = cln_w;
  const int nc = 1;
  // clean_conv = budget, in the Newton passes (a single solve keeps clean_tol in its own
  // norm): the rows over their gross terms S (CleanRightHandSide), W A X = W B with
  // W = 1/S, right-preconditioned by M^-1 W^-1, so that every norm below is the budget's
  // (cln_ewr -3: clean_clos1_budget's pass, in that norm to 0.1 clean_conv_max)
  const bool bud = cln_bon && (cln_ewr > -1.5 || cln_ewr < -2.5);
  const int sx = CLW_X, sb = bud ? CLW_BW : CLW_B;
  const int sr = CLW_R, sr0 = CLW_R0;
  const int sp = CLW_P, sv = CLW_V;
  const int ss = CLW_S, sq = CLW_Q;
  const int sy = CLW_Y, sz = CLW_Z;

  // w(dst) = sa w(a) + sb w(b)
  auto axpby = [&](const int dst, const int a, const Real sa, const int b,
                   const Real sb_) {
    par_for("cln_axpby", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,dst+c,k,j,i) = sa*w(m,a+c,k,j,i) + sb_*w(m,b+c,k,j,i);
      }
    });
  };
  // w(dst) = w(a) * w(s)^p
  auto scale = [&](const int dst, const int a, const int s, const Real p) {
    par_for("cln_scale", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,dst,k,j,i) = (p > 0.0) ? w(m,a,k,j,i)*w(m,s,k,j,i) : w(m,a,k,j,i)/w(m,s,k,j,i);
    });
  };
  auto apply = [&](const int src, const int dst) {
    CleanApply(src, dst);
    if (bud) scale(dst, dst, CLW_SE, -1.0);
  };
  auto precondition = [&](const int src, const int dst) {
    if (bud) {
      scale(CLW_PW, src, CLW_SE, 1.0);
      CleanPrecondition(CLW_PW, dst);
    } else {
      CleanPrecondition(src, dst);
    }
  };
  if (bud) scale(CLW_BW, CLW_B, CLW_SE, -1.0);
  auto dots = [&](const int n, const std::array<int, 4> ab) {
    return CleanDots(w, nc, n, ab);
  };

  // (re)start from the true residual with r0 = r.  A breakdown (r0.r at rounding level,
  // r0.v or omega zero) restarts; a second one right after a restart fails.  Updates and
  // their dots share a pass (CleanSums); four reductions per iteration.  Convergence of
  // the recursive residual is confirmed on the true one (<= 10 tol), else restart
  Real rho = 1.0, alpha = 1.0, omega = 1.0, rr = 0.0, r0r0 = 0.0, rho_new = 0.0, bb = 0.0;
  int it0 = 0;
  auto restart = [&](const int it, const bool force = false) {
    if (it == it0 + 1 && !force) return false;
    apply(sx, sr);
    const int n = (it < 0) ? 2 : 1;
    const auto d = CleanSums(pmy_pack, "cln_dot_rst", n,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0,
                  Real &s1) {
      for (int c = 0; c < nc; ++c) {   // (r = b - A x, r0 = r, p = v = 0)
        w(m,sr+c,k,j,i) = 1.0*w(m,sb+c,k,j,i) + -1.0*w(m,sr+c,k,j,i);
        w(m,sr0+c,k,j,i) = 1.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        w(m,sp+c,k,j,i) = 0.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        w(m,sv+c,k,j,i) = 0.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        s0 += w(m,sr+c,k,j,i)*w(m,sr+c,k,j,i);
        if (n > 1) s1 += w(m,sb+c,k,j,i)*w(m,sb+c,k,j,i);
      }
    });
    rho = alpha = omega = 1.0;
    rr = r0r0 = rho_new = d[0];
    if (it < 0) bb = d[1];
    it0 = it;
    return true;
  };
  restart(-1);
  Real tol = cln_tol*fmax(std::sqrt(bb), 1.0e-300);
  // clean_bud_maxnorm: the budget's tests on max_i |r_i|, the norm of the pass's own test
  // (an L2 norm above sqrt(N) t has its max above t: not reduced)
  const bool mn = bud && cln_bmn;
  const Real sqn = std::sqrt(static_cast<Real>(pmy_pack->pmesh->nmb_total)*
                             (ie - is + 1)*(je - js + 1)*(ke - ks + 1));
  auto vmax = [&](const int s) {
    Real a = 0.0;
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Kokkos::parallel_reduce("cln_vmax", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &mx) {
      const int m = idx/nkji, k = ks + (idx - m*nkji)/nji;
      const int j = js + (idx - m*nkji - (k - ks)*nji)/ni, i = is + idx%ni;
      mx = fmax(mx, fabs(w(m,s,k,j,i)));
    }, Kokkos::Max<Real>(a));
    return CleanReduce(a, true);
  };
  // Eisenstat-Walker forcing (choice 2, safeguarded) on |r0|, the Newton residual of the
  // pass in the budget's norm
  if (bud && cln_ewr < -2.5) {
    tol = 0.1*cln_bmax;
  } else if (bud) {
    const Real r0 = mn ? vmax(sr) : std::sqrt(rr);
    Real eta = 0.1;
    if (cln_ewr > 0.0) {
      eta = 0.9*SQR(r0/cln_ewr);
      if (0.9*SQR(cln_eta) > 0.1) eta = fmax(eta, 0.9*SQR(cln_eta));
      eta = fmin(eta, 0.1);
    }
    cln_ewr = r0;
    cln_eta = eta;
    tol = fmax(eta*r0, 0.1*cln_tolb);
  }
  // not with the budget: also
  // require max_i |r_i|/(a_i E_i) <= cnorm, a_i the diagonal, E_i the max of E^k and |X|
  // over the stencil (an M-matrix row's round-off bound, Oettli-Prager), since the |b|
  // test alone leaves faint regions unconverged.  If this stalls with the global test
  // met, it is at round-off: converged
  constexpr Real cnt = 1.0e-10;
  const Real efl = e_floor;
  auto conv = [&](const int s, const Real r2, const Real t) {
    if (!(std::sqrt(r2) <= t) && !(mn && std::sqrt(r2) <= sqn*t && vmax(s) <= t)) {
      return false;
    }
    if (bud) return true;
    CleanGhostFill(w, sx, 1, 0, 1);   // (X's ghosts are those of the last apply of X)
    Real a = 0.0;
    const int dj = pmy_pack->pmesh->multi_d, dk = pmy_pack->pmesh->three_d;
    // (a flat range; a max is order-independent)
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Kokkos::parallel_reduce("cln_cnorm", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &mx) {
      const int m = idx/nkji, k = ks + (idx - m*nkji)/nji;
      const int j = js + (idx - m*nkji - (k - ks)*nji)/ni, i = is + idx%ni;
      Real e = fmax(w(m,CLW_EK,k,j,i), efl);
      for (int o = -1; o <= 1; o += 2) {
        e = fmax(e, fmax(fabs(w(m,sx,k,j,i+o)), fabs(w(m,sx,k,j+o*dj,i))));
        e = fmax(e, fmax(fabs(w(m,sx,k+o*dk,j,i)), fabs(w(m,sx,k,j,i))));
      }
      mx = fmax(mx, fabs(w(m,s,k,j,i))/fmax(cnt*fabs(w(m,CLW_A0,k,j,i))*e, 1.0e-300));
    }, Kokkos::Max<Real>(a));
    return CleanReduce(a, true) <= 1.0;
  };
  if (conv(sr, rr, tol)) return 0;
  auto done = [&](const int it) {return restart(it, true) && conv(sr, rr, 10.0*tol);};

  const Real eps = std::numeric_limits<Real>::epsilon();
  for (int it = 1; it <= cln_maxit; ++it) {
    if (!(fabs(rho_new) > eps*std::sqrt(r0r0*rr))) {
      if (restart(it)) continue;
      return -it;
    }
    const Real beta = (rho_new/rho)*(alpha/omega);
    rho = rho_new;
    par_for("cln_pupdate", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,sp+c,k,j,i) = w(m,sr+c,k,j,i)
                          + beta*(w(m,sp+c,k,j,i) - omega*w(m,sv+c,k,j,i));
      }
    });
    precondition(sp, sy);
    apply(sy, sv);
    const Real r0v = dots(1, {sr0, sv, 0, 0})[0];
    if (!std::isfinite(r0v) || r0v == 0.0) {
      if (restart(it)) continue;
      return -it;
    }
    alpha = rho/r0v;
    const Real ma = -alpha;
    const Real s2 = CleanSums(pmy_pack, "cln_dot_s", 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0, Real &) {
      for (int c = 0; c < nc; ++c) {   // (s = r - alpha v)
        w(m,ss+c,k,j,i) = 1.0*w(m,sr+c,k,j,i) + ma*w(m,sv+c,k,j,i);
        s0 += w(m,ss+c,k,j,i)*w(m,ss+c,k,j,i);
      }
    })[0];
    if (conv(ss, s2, tol)) {
      axpby(sx, sx, 1.0, sy, alpha);
      if (done(it)) return it;
      continue;
    }
    precondition(ss, sz);
    apply(sz, sq);
    const auto d = dots(2, {sq, sq, sq, ss});
    omega = (d[0] > 0.0) ? d[1]/d[0] : 0.0;
    const Real al_ = alpha, om = omega;
    const auto e = CleanSums(pmy_pack, "cln_dot_x", 2,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0,
                  Real &s1) {
      for (int c = 0; c < nc; ++c) {
        w(m,sx+c,k,j,i) += al_*w(m,sy+c,k,j,i) + om*w(m,sz+c,k,j,i);
        w(m,sr+c,k,j,i) = w(m,ss+c,k,j,i) - om*w(m,sq+c,k,j,i);
        s0 += w(m,sr+c,k,j,i)*w(m,sr+c,k,j,i);
        s1 += w(m,sr0+c,k,j,i)*w(m,sr+c,k,j,i);
      }
    });
    rr = e[0];
    rho_new = e[1];
    if (conv(sr, rr, tol) && done(it)) return it;
    if (omega == 0.0 && !restart(it)) return -it;
  }
  return -cln_maxit;
}

}  // namespace radvet
