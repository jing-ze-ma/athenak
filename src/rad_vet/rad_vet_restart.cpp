//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_restart.cpp
//! \brief restart records of the module and the dt_retry hooks.
//!
//!  called from: driver/dt_retry.cpp, outputs/restart.cpp, pgen/pgen.cpp
//!  calls into: rad_vet_closure_ring.cpp, rad_vet_fct.cpp, rays/sc_tables.cpp

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_args.hpp"
#include "utils/deep_copy_across.hpp"

namespace radvet {

extern Real g_vgd_tab_alpha;   // rays/sc_tables.cpp
extern int64_t g_vgd_csgen[16];       // rays/sc_cs_setup.cpp; VET_VGD_CSCL_MAX
extern int64_t g_vgd_csgen_save[16];

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::RingRst*
//! \brief ring in restart files (version >= 4): offsets cln_rr then Delta
//! cln_rpd on interior cells in the file layout (RingFields; 0 where cln_rr holds
//! none; ema: 3 extra entries), flat per MeshBlock in RingRstN() ghosted slabs.
//! RingRstCopy packs (pk; zeros without have) or unpacks slabs c .. of array a
//! (c < 0: channel 0 is slab -c) for MeshBlocks m0 .. m0 + nmb - 1

int RadiationVet::RingRstN() const {
  auto &id = pmy_pack->pmesh->mb_indcs;
  int fm[CLN_N];
  const int nff = RingFields(fm, true);
  const int64_t nv = static_cast<int64_t>(RingNr(nff) + RingNd(nff))*id.nx1*
                     id.nx2*id.nx3;
  const int64_t ns = static_cast<int64_t>(id.nx1 + 2*id.ng)*
                     ((id.nx2 > 1) ? id.nx2 + 2*id.ng : 1)*
                     ((id.nx3 > 1) ? id.nx3 + 2*id.ng : 1);
  return static_cast<int>((nv + ns - 1)/ns);
}

void RadiationVet::RingRstCopy(DvceArray5D<Real> &a, const int c, const int nmb,
                               const int m0, const bool pk, const bool have) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  int fm[CLN_N], ff[CLN_N];
  const int nf = RingFields(fm), nff = RingFields(ff, true);
  Kokkos::Array<int, CLN_N> mp;   // (the file's field slot -> cln_rr's, -1: not held)
  for (int n = 0; n < CLN_N; ++n) {
    if (ff[n] >= 0) mp[ff[n]] = fm[n];
  }
  if (!pk) {   // (the writers of cln_rr bump vra_gen, DtRetry)
    vra_gen += 1;
    if (cln_rr.extent_int(0) != nmb) {
      Kokkos::realloc(cln_rr, nmb, RingNr(nf), id.nx3, id.nx2, id.nx1);
      Kokkos::realloc(cln_rpd, nmb, RingNd(nf), id.nx3, id.nx2, id.nx1);
    }
  }
  // (cln_rr's field groups: K - 1 offsets or the EMA's 1; later entries copied as is)
  const int gr = vra_ema ? 1 : vra_k - 1;
  const int64_t ni = static_cast<int64_t>(id.nx1)*id.nx2*id.nx3;
  const int64_t nr = RingNr(nf)*ni, ndm = RingNd(nf)*ni;
  const int64_t nro = RingNr(nff)*ni, nd = nro + RingNd(nff)*ni;
  const int64_t ns = static_cast<int64_t>(a.extent(2))*a.extent(3)*a.extent(4);
  const int64_t q0 = std::max(c, 0) - c, q1 = std::min<int64_t>(c + RingRstN(),
                                                                 a.extent(1)) - c;
  if (q1 <= q0) return;
  const int64_t nq = (q1 - q0)*ns, sa = a.stride(0);
  Real *pa = a.data();
  Real *pr = (have && cln_rr.extent_int(0) >= m0 + nmb) ? cln_rr.data() : nullptr;
  Real *pd = (pr != nullptr) ? cln_rpd.data() : nullptr;
  Kokkos::parallel_for("vet_ring_rst",
  Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<int64_t>>(DevExeSpace(), 0, nmb*nq),
  KOKKOS_LAMBDA(const int64_t p) {
    const int64_t m = p/nq, q = q0*ns + p - m*nq;
    Real *x = pa + m*sa + c*ns + q;
    Real *y = nullptr;
    if (pr != nullptr && q < nd) {
      const bool o = (q < nro);
      const int64_t u = o ? q : q - nro, ic = u%ni;
      const int g = static_cast<int>(u/ni), ng = o ? gr : 1;
      const int f = (g >= ng*nff) ? g - ng*(nff - nf) :
                    ((mp[g%nff] < 0) ? -1 : (g/nff)*nf + mp[g%nff]);
      if (f >= 0) y = (o ? pr + (m + m0)*nr : pd + (m + m0)*ndm) + f*ni + ic;
    }
    if (pk) {
      *x = (y != nullptr) ? *y : 0.0;
    } else if (y != nullptr) {
      *y = *x;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::FaxRstPack / FaxRstUnpack
//! \brief time2_face_fct: t2pf1 to / from cell-sized restart slabs (with ghosts), one
//! per direction; files with two slabs per direction: the first one is read

void RadiationVet::FaxRstPack(DvceArray5D<Real> out, int nmb) {
  const int n3 = out.extent_int(2), n2 = out.extent_int(3), n1 = out.extent_int(4);
  const int nd = FaxRstNslab();
  auto pf = t2pf1;
  par_for("pf_rst_pack", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int d = 0; d < nd; ++d) out(m,d,k,j,i) = pf(m,d,k,j,i);
  });
}

bool RadiationVet::FaxRstUnpack(DvceArray5D<Real> in, int nmb, int nsl) {
  const int n3 = in.extent_int(2), n2 = in.extent_int(3), n1 = in.extent_int(4);
  // file without them: the first stage takes the stage-start face F
  const int nd = trans_on ? (trans_x3 ? 3 : 2) : 1;
  if (t2_pfst && t2pf1.extent_int(0) >= nmb && (nsl == nd || nsl == 2*nd)) {
    auto pf = t2pf1;
    const int st = (nsl == nd) ? 1 : 2;
    par_for("pf_rst_unpack", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int d = 0; d < nd; ++d) pf(m,d,k,j,i) = in(m,st*d,k,j,i);
    });
    t2pf1_ok = true;
    return true;
  }
  if (t2_pfst && global_variable::my_rank == 0) {
    std::cout << "NOTE <rad_vet> restart: no time2_face_fct face fluxes in the file; the "
              << "first stage takes the stage-start face F" << std::endl;
  }
  return false;
}

using deep_copy_across::DeepCopyAcross;

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::RstNq
//! \brief closure components per record of the clean extrap history and ring
//! (RstSlot order): all CLN_N

int RadiationVet::RstNq() const {
  return CLN_N;
}

// closure_form = iterate: the closure in place is restart state (none in a file of
// another closure, Time2RstAccept)
int RadiationVet::RstNcl() const {
  if (t2_rcl >= 0) return t2_rcl;  // a file's (Time2RstAccept)
  return cln_closit ? NumClosureSlots(pmy_pack->pmesh->use_cubed_sphere, cln_rfh) : 0;
}

int RadiationVet::RstNfv() const {
  // cubed sphere sc rays: their cln_fvi slots (before the first sweep: the file's);
  // version >= 6: the face half ranges cln_fh instead (none where every solve rebuilds
  // them)
  if (cln_cssc && t2_rstv >= 6) return FhBuilt() ? 0 : 12;
  if (cln_cssc && cln_fvi.extent_int(0) > 0) return cln_fvi.extent_int(1);
  if (cln_cssc) return (cln_fvn > 0) ? cln_fvn : t2_nfv;
  return cln_scmem ? cln_top.extent_int(1) : 0;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::Time2Rst*
//! \brief the restart channels: K1 (VET_T2_NK), vet_prev (VET_T2_NVET, vet_sc), then the
//! clean scheme's state; files of version < 4 also hold ipred2 and extrap2 channels
//! (read past)

int RadiationVet::Time2RstNchWant() {
  // version >= 4: no ipred2 / extrap2 channels, ring offsets dense (RingRstCopy)
  const bool sk = t2_rstv >= 4;
  const int np = (impl_pord == 2) ? 5 : 3;
  const int nt = (time_scheme != VET_TIME_SDIRK2) ? 0 :
         VET_T2_NK + (sk ? 0 : np) + (vet_sc ? VET_T2_NVET : 0)
         + ((t2_x2old && !sk) ? VET_X2_NCH : 0);   // extrap2
  return nt + ((force_ref == VET_FREF_CLEAN) ? 4 + cln_lzt : 0)   // the reference, set
         + ((force_ref == VET_FREF_WB_MD) ? 10 : 0)           // wb_arad_md: 9, set
         + ((cln_extrap > 0) ? 2*RstNq() + 1 : 0)    // extrap history
         + RstNcl()                                     // iterate: the closure
         + (cln_ring ? (t2_rstv < 4 ? vra_k*RstNq() : RingRstN()) + 1
                                : 0)                          // ring
         + ((cln_fsolve && cln_fkin) ? 6 : 0)       // stag fkin: d and its K1
         + RstNfv();                                    // the rays' intensities
}

bool RadiationVet::Time2RstAccept(const int nfile, const int ver) {
  // version 1 holds d without its slope K1_d: K1_d = 0, backward-Euler first step
  t2_rstd3 = (ver < 2) && cln_fsolve && cln_fkin;
  // version 5 on the cubed sphere (sc ray intensities) is refused
  if (cln_cssc && ver == 5) {
    std::cout << "### FATAL ERROR <rad_vet> restart: a version-5 clean file on the cubed "
              << "sphere (the sc rays' intensities) is no longer read" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  t2_rstv = ver;
  t2_rcl = -1;
  int nb = Time2RstNchWant() - RstNfv() - (t2_rstd3 ? 3 : 0);
  t2_nfv = (cln_cssc || cln_scmem) ? nfile - nb : 0;
  // a file of another closure_form: the iterate closure is skipped, or missing (the
  // first solve takes the predicted closure)
  const int ns = NumClosureSlots(pmy_pack->pmesh->use_cubed_sphere, cln_rfh);
  auto fits = [&](const int n) {
    return cln_cssc ? (nfile - n == 0 || nfile - n == 12) : (nfile == n);
  };
  if (!cln_scmem && !fits(nb) && fits(cln_closit ? nb - ns : nb + ns)) {
    t2_rcl = cln_closit ? 0 : ns;
    nb += cln_closit ? -ns : ns;
    t2_nfv = cln_cssc ? nfile - nb : 0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> NOTE: restart from a file of another closure_form: "
                << (cln_closit ? "no closure in it, the first solve takes the predicted "
                    "one" : "its closure (closure_form = iterate) is not read")
                << std::endl;
    }
  }
  // (the ring's half ranges, cln_rfh: no cln_fh channels; a file of a run with them is
  // refused without them, and one without them with them: 1 .. 11 channels over)
  const bool fho = cln_rfh ? t2_nfv != 0 : (t2_nfv > 0 && t2_nfv < 12);
  if (cln_cssc && ver >= 6 && fho) {
    std::cout << "### FATAL ERROR <rad_vet> restart: a file " << (cln_rfh ? "not " : "")
              << "written with the ring's half ranges (vet_rotate_avg_ema_wb on the "
              << "cubed sphere; " << t2_nfv << " channels over)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // cubed sphere, version >= 6: the face half ranges are the last 0 or 12 channels; a
  // surplus is a ring file: with the ring off its channels are skipped, else
  // refused
  t2_rskip = 0;
  if (cln_cssc && ver >= 6 && t2_nfv > 12) {
    t2_rskip = cln_ring ? -1 : RingRstN() + 1;
    const int nh = t2_nfv - t2_rskip;
    if (t2_rskip < 0 || (nh != 0 && nh != 12)) {
      std::cout << "### FATAL ERROR <rad_vet> restart: the clean block has " << nfile
                << " channels, not this run's layout (" << nb << " + 0 or 12) nor it "
                << "with the ring's offsets" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    t2_nfv = nh;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> NOTE: restart with ring off: the file's ring "
                << "offsets (" << t2_rskip << " channels) are not read" << std::endl;
    }
  }
  return (nfile == nb + t2_nfv + t2_rskip) && (t2_nfv >= 0);
}

int RadiationVet::Time2RstNch() {
  t2_rstv = 7;   // (the writer: version 7, after a restart from an older file too)
  t2_rcl = -1;
  // the clean state (and K1 with cln_k1ok) is written in backward-Euler steps too;
  // restart.cpp stores t2_ok
  return Time2RstNchWant();
}

// m0 > 0 packs MeshBlocks m0 .. m0+nmb-1 into a(0 .. nmb-1) (per-MB staging)
void RadiationVet::Time2RstPack(DvceArray5D<Real> &a, int nmb, int m0, int w0) {
  auto mb = std::make_pair(0, nmb);
  auto ms = std::make_pair(m0, m0 + nmb);
  auto AL = Kokkos::ALL;
  const bool h2 = (time_scheme == VET_TIME_SDIRK2);
  // a holds the record's channels w0 .. w1-1 (the writer streams the record in windows);
  // put copies the record channels c0 .. c0+n-1 = src's s0 .. s0+n-1 that fall inside,
  // one sets channel cc to v
  const int w1 = w0 + a.extent_int(1);
  auto put = [&](const int c0, const int n, const auto &src, const int s0) {
    const int lo = std::max(c0, w0), hi = std::min(c0 + n, w1);
    if (lo < hi) {
      Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(lo - w0, hi - w0), AL, AL,
                                        AL),
                        Kokkos::subview(src, ms,
                                        std::make_pair(s0 + lo - c0, s0 + hi - c0),
                                        AL, AL, AL));
    }
  };
  auto one = [&](const int cc, const Real v) {
    if (cc >= w0 && cc < w1) {
      Kokkos::deep_copy(Kokkos::subview(a, mb, cc - w0, AL, AL, AL), v);
    }
  };
  const int na3 = static_cast<int>(a.extent(2)), na2 = static_cast<int>(a.extent(3));
  const int na1 = static_cast<int>(a.extent(4)), nw = a.extent_int(1);
  int c = 0;
  if (h2) {
    put(c, VET_T2_NK, t2k1, 0);
    c += VET_T2_NK;
  }
  // version >= 4 layout: no predictor / extrap2 channels
  if (h2 && vet_sc) {
    put(c, VET_T2_NVET, vet_prev, 0);
    c += VET_T2_NVET;
  }
  // force_reference = clean (wb_arad_md): the reference force (and potentials) and
  // whether it is set
  if (force_ref == VET_FREF_CLEAN || force_ref == VET_FREF_WB_MD) {
    const int nr = (force_ref == VET_FREF_WB_MD) ? 9 : 3 + cln_lzt;
    const bool have = (cln_fref.extent_int(0) >= m0 + nmb);
    if (have) put(c, nr, cln_fref, 0);
    one(c + nr, (have && cln_fref_set) ? 1.0 : 0.0);
    c += nr + 1;
  }
  // closure_form = extrap*: the ray moments at t^n and at the end of the step, then
  // t^n, that end time and the flags at i = 0..4
  const int nq = RstNq();
  if (cln_extrap > 0) {
    const bool have = (cln_mn.extent_int(0) >= m0 + nmb);
    auto a_ = a;
    auto &id = pmy_pack->pmesh->mb_indcs;
    const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
    for (int n = 0; n < 2*nq && have; ++n) {   // (ghosts and slots this mesh lacks: 0)
      auto mx = (n < nq) ? cln_mn : cln_mc;
      const int s = RstSlot(n%nq), cc = c + n - w0;
      if (cc < 0 || cc >= nw) continue;
      if (s >= mx.extent_int(1)) {
        one(c + n, 0.0);
        continue;
      }
      par_for("vet_t2_clnrsm", DevExeSpace(), 0, nmb - 1, 0, na3 - 1, 0, na2 - 1, 0,
              na1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const bool in = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
        a_(m,cc,k,j,i) = in ? mx(m+m0,s,k-ks,j-js,i-is) : 0.0;
      });
    }
    const int cq = c + 2*nq - w0;
    const Real v0 = cln_tn, v1 = cln_tm, v2 = (have && cln_have_n) ? 1.0 : 0.0;
    const Real v3 = (have && cln_have_c) ? 1.0 : 0.0, v4 = cln_cand_be ? 1.0 : 0.0;
    if (cq >= 0 && cq < nw) {
      par_for("vet_t2_clnrst", DevExeSpace(), 0, nmb - 1, 0, na3 - 1, 0, na2 - 1, 0,
              na1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        a_(m,cq,k,j,i) = (i == 0) ? v0 : ((i == 1) ? v1 : ((i == 2) ? v2 :
                         ((i == 3) ? v3 : ((i == 4) ? v4 : 0.0))));
      });
    }
    c += 2*nq + 1;
  }
  // closure_form = iterate: the closure in place, the last iterate's (RstNcl)
  if (RstNcl() > 0) {
    if (cln_s.extent_int(0) >= m0 + nmb) put(c, RstNcl(), cln_s, 0);
    c += RstNcl();
  }
  // ring: the offsets and the pending Delta of the fields the ring carries
  // (interior cells, dense: RingRstCopy), then the ring's flag, block, the pending
  // block and its time at i = 0..3
  if (cln_ring) {
    const int nd = RingRstN();
    const bool have = cln_rok && (cln_rr.extent_int(0) >= m0 + nmb);
    RingRstCopy(a, c - w0, nmb, m0, true, have);
    auto a_ = a;
    const int cq = c + nd - w0;
    const Real v0 = have ? 1.0 : 0.0, v1 = static_cast<Real>(cln_rb);
    const Real v2 = static_cast<Real>(cln_rpb), v3 = cln_rpt;
    if (cq >= 0 && cq < nw) {
      par_for("vet_t2_ringrst", DevExeSpace(), 0, nmb - 1, 0, na3 - 1, 0, na2 - 1, 0,
              na1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        a_(m,cq,k,j,i) = (i == 0) ? v0 : ((i == 1) ? v1 : ((i == 2) ? v2 :
                         ((i == 3) ? v3 : 0.0)));
      });
    }
    c += nd + 1;
  }
  // clean_stag_fkin: the energy flux memory d of the three face families (CLF_N + 3),
  // then its slope K1 (CLF_KD); before the first solve, those read from the file
  if (cln_fsolve && cln_fkin) {
    for (int d = 0; d < 3; ++d) {
      auto fc = (d == 0) ? cln_fc1 : ((d == 1) ? cln_fc2 : cln_fc3);
      if (!cln_init && cln_rstd.extent_int(0) >= m0 + nmb) {
        put(c, 6, cln_rstd, 0);
        break;
      }
      if (fc.extent_int(0) < m0 + nmb) continue;
      put(c + d, 1, fc.s, CLF_N + 3);
      put(c + 3 + d, 1, fc.s, CLF_KD);
    }
    c += 6;
  }
  // (a window that does not reach the rays' intensities: restart.cpp windows only a
  // record without them)
  if (w0 > 0 || w1 < c + RstNfv()) return;
  // the sc rays' intensities the next build reads, last (unless a has no room for them:
  // restart.cpp copies them itself)
  if (a.extent_int(1) < c + RstNfv()) return;
  if (cln_cssc && RstNfv() > 0) {   // (version 6: cln_fh)
    auto af = Kokkos::subview(a, mb, std::make_pair(c, c + 12), AL, AL, AL);
    if (cln_fh.extent_int(0) >= m0 + nmb) {
      Kokkos::deep_copy(af, Kokkos::subview(cln_fh, ms, AL, AL, AL, AL));
    } else {
      Kokkos::deep_copy(af, 0.0);
    }
  } else if (RstNfv() > 0 && RstFv().extent_int(0) >= m0 + nmb) {
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c + RstNfv()), AL, AL,
                                      AL),
                      Kokkos::subview(RstFv(), ms, AL, AL, AL, AL));
  }
}

void RadiationVet::Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb) {
  auto mb = std::make_pair(0, nmb);
  auto AL = Kokkos::ALL;
  const bool h2 = (time_scheme == VET_TIME_SDIRK2);
  // a file's predictor / extrap2 channels are read past (a NOTE once if nonzero)
  auto skip = [&w]() {
    static bool said = false;
    for (size_t q = 0; q < w.span() && !said && global_variable::my_rank == 0; ++q) {
      if (w.data()[q] == 0.0) continue;
      std::cout << "<rad_vet> NOTE: the restart file's implicit predictor / extrap2 "
                << "channels are skipped" << std::endl;
      said = true;
    }
  };
  // the last t2_nfv channels: the sc rays' intensities (their number is not
  // known before the first sweep, Time2RstAccept)
  const int nfv0 = Time2RstNchWant() - RstNfv() - (t2_rstd3 ? 3 : 0) + t2_rskip;
  if (cln_cssc && t2_rstv >= 6 && ch >= nfv0) {   // (cln_fh)
    if (cln_fh.extent_int(0) != nmb) {
      Kokkos::realloc(cln_fh, nmb, 12, w.extent_int(1), w.extent_int(2), w.extent_int(3));
    }
    DeepCopyAcross(Kokkos::subview(cln_fh, mb, ch - nfv0, AL, AL, AL), w);
    return;
  }
  if ((cln_cssc || cln_scmem) && ch >= nfv0) {
    auto &fv = RstFv();
    if (fv.extent_int(0) != nmb || fv.extent_int(1) != t2_nfv) {
      Kokkos::realloc(fv, nmb, t2_nfv, w.extent_int(1), w.extent_int(2), w.extent_int(3));
    }
    DeepCopyAcross(Kokkos::subview(fv, mb, ch - nfv0, AL, AL, AL), w);
    return;
  }
  if (h2 && ch < VET_T2_NK) {
    DeepCopyAcross(Kokkos::subview(t2k1, mb, ch, AL, AL, AL), w);
    return;
  }
  if (h2) ch -= VET_T2_NK;
  if (h2 && t2_rstv < 4) {
    const int np = (impl_pord == 2) ? 5 : 3;
    if (ch < np) {
      skip();
      return;
    }
    ch -= np;
  }
  if (h2 && vet_sc) {
    if (ch < VET_T2_NVET) {
      DeepCopyAcross(Kokkos::subview(vet_prev, mb, ch, AL, AL, AL), w);
      return;
    }
    ch -= VET_T2_NVET;
  }
  if (force_ref == VET_FREF_CLEAN || force_ref == VET_FREF_WB_MD) {
    const int nr = (force_ref == VET_FREF_WB_MD) ? 9 : 3 + cln_lzt;
    if (ch < nr) {
      if (cln_fref.extent_int(0) != nmb) {
        Kokkos::realloc(cln_fref, nmb, nr, w.extent_int(1), w.extent_int(2),
                        w.extent_int(3));
      }
      DeepCopyAcross(Kokkos::subview(cln_fref, mb, ch, AL, AL, AL), w);
      return;
    }
    if (ch == nr) {
      cln_fref_set = (w(0,0,0,0) > 0.5);
      return;
    }
    ch -= nr + 1;
  }
  const int nq = RstNq();
  if (cln_extrap > 0) {
    if (ch < 2*nq) {
      auto &mx = (ch < nq) ? cln_mn : cln_mc;
      const int ns = NumClosureSlots(pmy_pack->pmesh->use_cubed_sphere, cln_rfh);
      const int s = RstSlot(ch%nq);
      auto &id = pmy_pack->pmesh->mb_indcs;
      const int is = id.is, js = id.js, ks = id.ks;
      if (mx.extent_int(0) != nmb) {
        Kokkos::realloc(mx, nmb, ns, id.nx3, id.nx2, id.nx1);
      }
      if (s >= ns) return;
      DvceArray4D<Real> t("vet_t2_clnrd", nmb, w.extent_int(1), w.extent_int(2),
                          w.extent_int(3));
      DeepCopyAcross(t, w);
      auto mx_ = mx;
      par_for("vet_t2_clnset", DevExeSpace(), 0, nmb - 1, 0, id.nx3 - 1, 0, id.nx2 - 1, 0,
              id.nx1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        mx_(m,s,k,j,i) = t(m,k+ks,j+js,i+is);
      });
      return;
    }
    if (ch == 2*nq) {
      cln_tn = w(0,0,0,0);
      cln_tm = w(0,0,0,1);
      cln_have_n = (w(0,0,0,2) > 0.5);
      cln_have_c = (w(0,0,0,3) > 0.5);
      cln_cand_be = (w(0,0,0,4) > 0.5);
      return;
    }
    ch -= 2*nq + 1;
  }
  if (RstNcl() > 0) {   // (closure_form = iterate: cln_s, with its ghosts)
    if (ch < RstNcl()) {
      if (!cln_closit) return;   // (an iterate run's file, read by another closure)
      if (cln_s.extent_int(0) != nmb) {
        Kokkos::realloc(cln_s, nmb, RstNcl(), w.extent_int(1), w.extent_int(2),
                        w.extent_int(3));
      }
      DeepCopyAcross(Kokkos::subview(cln_s, mb, ch, AL, AL, AL), w);
      cln_cip = true;
      return;
    }
    ch -= RstNcl();
  }
  if (t2_rskip > 0) {  // the file's ring with the ring off: not read
    if (ch < t2_rskip) return;
    ch -= t2_rskip;
  }
  if (cln_ring && t2_rstv >= 4) {   // (version 4, RingRstCopy)
    if (ch < RingRstN()) {
      DvceArray5D<Real> t("vet_t2_ringrd", nmb, 1, w.extent_int(1), w.extent_int(2),
                          w.extent_int(3));
      DeepCopyAcross(Kokkos::subview(t, AL, 0, AL, AL, AL), w);
      RingRstCopy(t, -ch, nmb, 0, false, true);
      return;
    }
    ch -= RingRstN() - vra_k*nq;  // the flag slab
  }
  if (cln_ring) {
    const int nr = (vra_k - 1)*nq;
    if (ch < nr + nq) {
      // (the interior cells of the carried fields, RingMoments)
      auto &rx = (ch < nr) ? cln_rr : cln_rpd;
      vra_gen += 1;
      auto &id = pmy_pack->pmesh->mb_indcs;
      const int is = id.is, js = id.js, ks = id.ks;
      Kokkos::Array<int, CLN_N> fm;
      const int nf = RingFields(fm.data());
      const int f = fm[RstSlot((ch < nr) ? ch%nq : ch - nr)];
      const int nx = (ch < nr) ? (vra_k - 1)*nf : nf;
      if (rx.extent_int(0) != nmb) {
        Kokkos::realloc(rx, nmb, nx, id.nx3, id.nx2, id.nx1);
      }
      if (f < 0) return;
      const int q = (ch < nr) ? (ch/nq)*nf + f : f;
      DvceArray4D<Real> t("vet_t2_ringrd", nmb, w.extent_int(1), w.extent_int(2),
                          w.extent_int(3));
      DeepCopyAcross(t, w);
      auto rx_ = rx;
      par_for("vet_t2_ringset", DevExeSpace(), 0, nmb - 1, 0, id.nx3 - 1, 0, id.nx2 - 1,
              0, id.nx1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        rx_(m,q,k,j,i) = t(m,k+ks,j+js,i+is);
      });
      return;
    }
    if (ch == nr + nq) {
      cln_rok = (w(0,0,0,0) > 0.5);
      cln_rb = static_cast<int64_t>(std::llround(w(0,0,0,1)));
      cln_rpb = static_cast<int64_t>(std::llround(w(0,0,0,2)));
      cln_rpt = w(0,0,0,3);
      return;
    }
    ch -= nr + nq + 1;
  }
  if (cln_fsolve && cln_fkin) {   // (InitArrays copies it to the face arrays)
    const int nd = t2_rstd3 ? 3 : 6;
    if (ch < nd) {
      if (cln_rstd.extent_int(0) != nmb) {   // (zeroed: K1_d = 0 for a version-1 file)
        Kokkos::realloc(cln_rstd, nmb, 6, w.extent_int(1), w.extent_int(2),
                        w.extent_int(3));
      }
      DeepCopyAcross(Kokkos::subview(cln_rstd, mb, ch, AL, AL, AL), w);
      return;
    }
    ch -= nd;
  }
  // extrap2 (version < 4): read past
  if (h2 && t2_x2old) skip();
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::X1eRst*
//! \brief vet_x1_periodic (adaptive): the exit plane the next build starts from, to /
//! from cell-sized restart slabs (kVetX1eRstMagic): ray r of row (k, j) at slab r / n1,
//! cell i = r % n1

void RadiationVet::X1eRstPack(DvceArray5D<Real> &dst, const int nmb) const {
  auto e_ = vet_x1ent;
  auto d_ = dst;
  const int nray = e_.extent_int(1), n1 = dst.extent_int(4);
  Kokkos::deep_copy(dst, 0.0);
  par_for("vet_x1e_pack", DevExeSpace(), 0, nmb - 1, 0, nray - 1,
          0, dst.extent_int(2) - 1, 0, dst.extent_int(3) - 1,
  KOKKOS_LAMBDA(const int m, const int r, const int k, const int j) {
    d_(m,r/n1,k,j,r%n1) = e_(m,r,k,j);
  });
}

bool RadiationVet::X1eRstUnpack(const DvceArray5D<Real> &src, const int nmb,
                                const int nsl) {
  auto e_ = vet_x1ent;
  auto s_ = src;
  const int nray = e_.extent_int(1), n1 = src.extent_int(4);
  if (!vet_x1per || vet_x1ptol <= 0.0 || e_.extent_int(0) < nmb ||
      nsl != (nray + n1 - 1)/n1) {return false;}
  par_for("vet_x1e_unpack", DevExeSpace(), 0, nmb - 1, 0, nray - 1, 0,
          src.extent_int(2) - 1, 0, src.extent_int(3) - 1,
  KOKKOS_LAMBDA(const int m, const int r, const int k, const int j) {
    e_(m,r,k,j) = s_(m,r/n1,k,j,r%n1);
  });
  vet_x1_warm = true;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::DtRetryQuiesce
//! \brief before the dt-retry checkpoint is taken (restoring = false) or put back
//! (true): join a pending symmetric-reference solve and, before a restore, the table
//! prefetch. Changes no number (vet_gd_async is refused with dt_retry).

void RadiationVet::DtRetryQuiesce(bool restoring) {
  if (!restoring) {
    // compact-list build generations at the checkpoint
    for (int c = 0; c < vgd_csncl; ++c) {g_vgd_csgen_save[c] = g_vgd_csgen[c];}
    RingPreShift();
    FluxMemoryRetry(true);
    return;
  }
  if (vgd_tab_state == 1 && vgd_tabf.valid()) {vgd_tabf.wait();}
  if (vgd_on) {vgd_tex.fence();}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::DtRetryAfterRestore
//! \brief after the checkpoint is put back: reset the table prefetch state (not
//! checkpointed) to "nothing prefetched"; the next rotation computes its tables
//! synchronously (ScTabCompute: same numbers).

void RadiationVet::DtRetryAfterRestore() {
  FluxMemoryRetry(false);
  if (vgd_tabf.valid()) {vgd_tabf.get();}
  if (vgd_on) {
    vgd_tex.fence();
    vgd_tab_state = 0;
    vgd_tabh.alpha = -9.0;
    // the sweep's halo caches (vgd_hle, vgd_hxe, mask slots vgd_htag, overlap groups)
    // have checkpointed keys but arrays remade in place: a failed attempt may leave a
    // restored key naming another cut's lists. Drop every record; the next sweep
    // remakes the ones it needs (same masks: same numbers).
    for (auto &x : vgd_hle) {x.gen = -1; x.alpha = -9.0; x.nofit = false;}
    for (auto &x : vgd_hxe) {x.gen = -1; x.alpha = -9.0; x.nofit = false;}
    for (auto &x : vgd_hxe_alt) {x.gen = -1; x.alpha = -9.0; x.nofit = false;}
    vgd_shl_alpha = vgd_shl_alpha_xb = -3.0;   // the shell lists (vgd_shl[_xb]) too
    for (int s = 0; s < 2; ++s) {vgd_htag[s][0] = -1;}
    // the cubed-sphere compact lists (vgd_cscl) are the same kind of record: drop an
    // entry the failed attempt rebuilt (LRU eviction), else the restored key names the
    // evicted set's arrays; the next sweep remakes it (same numbers)
    for (int c = 0; c < vgd_csncl; ++c) {
      if (g_vgd_csgen[c] != g_vgd_csgen_save[c]) {
        vgd_cscl[c].alpha = -9.0;
        vgd_cscl[c].use = 0;
      }
    }
    // the two read masks may hold other sets than the restored handles say: every
    // cached set remakes its mask at its next use (ScCsCompact)
    for (auto &e : vgd_cscl) {e.msk = decltype(e.msk)();}
    // a failed attempt on the last step of a rotation block left the tables on the
    // next set while the restored vgd_alpha names the current one: remake the tables
    // of the restored angle (ScTabCompute: same numbers)
    if (vgd_alpha >= 0.0 && g_vgd_tab_alpha != vgd_alpha) {
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> dt_retry: vet_gd tables back from angle "
                  << g_vgd_tab_alpha << " to the restored " << vgd_alpha << std::endl;
      }
      ScTables(vgd_alpha);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::RotAvgRstLoad
//! \brief the vet_rotate_avg history of an older restart file: not read (retired with the
//! old Cartesian sweep); under the ring an EMA history marks its closure as converged
//! (cln_rhist), as before

bool RadiationVet::RotAvgRstLoad(const int64_t *h, const DvceArray5D<Real> &,
                                 const int) {
  cln_rhist = cln_ring && vra_ema && h[0] == 1;
  if (!cln_ring && h[0] == 1 && global_variable::my_rank == 0) {
    std::cout << "### NOTE: <rad_vet> the vet_rotate_avg history of this restart file "
              << "(the retired moment averaging of the old Cartesian sweep) is not read"
              << std::endl;
  }
  return false;
}

}  // namespace radvet
