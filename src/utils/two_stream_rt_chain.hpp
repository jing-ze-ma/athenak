//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
#ifndef UTILS_TWO_STREAM_RT_CHAIN_HPP_
#define UTILS_TWO_STREAM_RT_CHAIN_HPP_
//! \file two_stream_rt_chain.hpp
//! \brief the BODIES of TsrtCkChain and TsrtGreyRt (declared in two_stream_rt.hpp).
//! tusplit-1102: moved out of two_stream_rt.hpp verbatim, so that their instantiations
//! (the bulk of the device compile of the old single translation unit) are compiled in
//! their own translation units, utils/two_stream_rt_chain_*.cpp.  Include it ONLY there.

#include <algorithm>
#include <iostream>
#include <tuple>
#include <type_traits>

#include "utils/two_stream_rt.hpp"

namespace two_stream_rt {

//----------------------------------------------------------------------------------------
//! \fn TsrtCkChain
//! \brief the correlated-k chain launch of picket_fence_two_stream_RT_pass
//! (its launch_ck_chain), one instantiation per compile-time tag set.  A function
//! template rather than a generic lambda only because nvcc forbids extended
//! lambdas inside generic lambdas; ctx_ holds references to the locals of
//! picket_fence_two_stream_RT_pass the body reads (see ck_chain_ctx there).
template <typename Ctx, typename NnTag, typename SphTag, typename BspTag, typename CchTag,
          typename FrmTag, typename FopTag, typename JacTag>
void TsrtCkChain(Ctx &ctx_, NnTag nn_tag, SphTag sph_tag, BspTag bsp_tag,
                 CchTag cch_tag, FrmTag frm_tag, FopTag fop_tag, JacTag jac_tag) {
  auto &ACC = std::get<0>(ctx_);
  auto &AFC = std::get<1>(ctx_);
  auto &albedo = std::get<2>(ctx_);
  auto &ap = std::get<3>(ctx_);
  auto &Bb_g = std::get<4>(ctx_);
  auto &bface_on = std::get<5>(ctx_);
  auto &boltz_sigma = std::get<6>(ctx_);
  auto &BTF = std::get<7>(ctx_);
  auto &cbt_ = std::get<8>(ctx_);
  auto &cf_g = std::get<9>(ctx_);
  auto &ckc0_g = std::get<10>(ctx_);
  auto &ckci_g = std::get<11>(ctx_);
  auto &ckco_g = std::get<12>(ctx_);
  auto &ckdb_g = std::get<13>(ctx_);
  auto &ckdif_ = std::get<14>(ctx_);
  auto &ckdone_g = std::get<15>(ctx_);
  auto &ckfcf_ = std::get<16>(ctx_);
  auto &ckfst_ = std::get<17>(ctx_);
  auto &ckfus_ = std::get<18>(ctx_);
  auto &ckgw = std::get<19>(ctx_);
  auto &ckjac_g = std::get<20>(ctx_);
  auto &ckjacp_ = std::get<21>(ctx_);
  auto &ckkro_g = std::get<22>(ctx_);
  auto &cklk = std::get<23>(ctx_);
  auto &ck_nq_ = std::get<24>(ctx_);
  auto &ckpf = std::get<25>(ctx_);
  auto &ckskip_ = std::get<26>(ctx_);
  auto &ckswf = std::get<27>(ctx_);
  auto &cktpf_g = std::get<28>(ctx_);
  auto &difg_g = std::get<29>(ctx_);
  auto &dx1 = std::get<30>(ctx_);
  auto &Em_g = std::get<31>(ctx_);
  auto &Fb_g = std::get<32>(ctx_);
  auto &Fstar = std::get<33>(ctx_);
  auto &grav = std::get<34>(ctx_);
  auto &grav_pmass = std::get<35>(ctx_);
  auto &icc_g = std::get<36>(ctx_);
  auto &icut_g = std::get<37>(ctx_);
  auto &ie = std::get<38>(ctx_);
  auto &int_at_cut = std::get<39>(ctx_);
  auto &is = std::get<40>(ctx_);
  auto &je = std::get<41>(ctx_);
  auto &js = std::get<42>(ctx_);
  auto &kc_g = std::get<43>(ctx_);
  auto &ke = std::get<44>(ctx_);
  auto &ks = std::get<45>(ctx_);
  auto &layer_legacy = std::get<46>(ctx_);
  auto &MIXF = std::get<47>(ctx_);
  auto &mug = std::get<48>(ctx_);
  auto &n1 = std::get<49>(ctx_);
  auto &nblk = std::get<50>(ctx_);
  auto &nmb1 = std::get<51>(ctx_);
  auto &omega = std::get<52>(ctx_);
  auto &pb_g = std::get<53>(ctx_);
  auto &pfid = std::get<54>(ctx_);
  auto &pfl0 = std::get<55>(ctx_);
  auto &Qb_g = std::get<56>(ctx_);
  auto &rhoN = std::get<57>(ctx_);
  auto &Src_g = std::get<58>(ctx_);
  auto &T_g = std::get<59>(ctx_);
  auto &tide = std::get<60>(ctx_);
  auto &Tint = std::get<61>(ctx_);
  auto &Tint4 = std::get<62>(ctx_);
  auto &wg = std::get<63>(ctx_);
  auto &X1F = std::get<64>(ctx_);
  auto &x1v_ = std::get<65>(ctx_);
  auto &xP_g = std::get<66>(ctx_);
  auto &xT_g = std::get<67>(ctx_);
  const int fce = std::get<68>(ctx_);    // problem/ck_sph_face
  const bool dil = (fce == 1 || fce == 3);   // the dilute face (ck_sph_dilute)
  const bool fbl = (fce == 2 || fce == 4);   // the thick/thin blend
  const bool pth = (fce == 3 || fce == 4);  // + the up ray's cone path factor
  const int tpm = std::get<69>(ctx_);        // problem/ck_sph_top
  auto &vfs_g = std::get<70>(ctx_);          // ck_sph_face = 5: s = sqrt(3 f)
  auto &vfg_g = std::get<71>(ctx_);          // ... gamma per face, zeta at ie+2
  const bool vef = (fce == 5);
  const bool wex = std::get<72>(ctx_);       // problem/ck_wall_flux_exact
  constexpr bool SPH = decltype(sph_tag)::value;
  constexpr bool BSP = decltype(bsp_tag)::value;
  constexpr int CCH = decltype(cch_tag)::value;    // see ck_sweep_cache
  // FRM, JAC and gJ are static: nvcc counts a local constexpr that the chain kernel
  // first names inside `if constexpr` as a capture, which it does not allow
  static constexpr int FRM = decltype(frm_tag)::value;   // see ck_sweep_form
  // problem/ck_impl_frozen_op: a COMPILE-TIME tag, exactly like CCH.  The frozen
  // operator is reachable only under ck_implicit, so with FOP = 0 the store and
  // re-apply branches below (about eight per (cell, chain)) and the loads from
  // the five stored-operator Views fold away and this kernel is the one the
  // explicit sweep compiled to before the flag existed.  A default-off switch
  // that stays a runtime bool inside this kernel costs ~28 % of the whole RT
  // loop; that is measured in bench/bisect_cost/README.md.
  constexpr bool FOP = decltype(fop_tag)::value;   // see ck_impl_frozen_op
  // problem/ck_implicit on the tm sweep (ck_sweep_form = 1): the tridiagonal
  // assembly of the probe-free recurrence, a COMPILE-TIME tag for the same reason
  // as FOP.  It is true only in the FRM = 1 instantiation launched on a pass that
  // builds the Jacobian (ckjacp_), so the production kernel (implicit off) and a
  // pass that reuses the Jacobian run the JAC = 0 code, which is the tm sweep
  // exactly as it was.  See the JAC blocks in the tm body.
  static constexpr bool JAC = decltype(jac_tag)::value;
  // ck-scratch: the whole-column arrays' groups in rt_ckscr_ptr (see CkScrColI).
  // gI..gJ are this instantiation's group numbers; a switched-off array has none.
  constexpr int gI = 0;
  constexpr int gC = 1;
  constexpr int gK = 1 + (SPH ? 1 : 0);
  constexpr int gP = gK + (BSP ? 1 : 0);
  static constexpr int gJ = gP + ((CCH >= 2) ? 4 : 0);
  constexpr int nGrp = gJ + (JAC ? 3 : 0);
  static constexpr int nScS = gJ + (JAC ? 3 : 0);   // = nGrp: see CkScrColI
  const int scn1 = n1;
  const int scn2 = je - js + 1, scn3 = ke - ks + 1;
  const size_t grp1 = static_cast<size_t>(RT_NB)*scn1*(nmb1 + 1)*scn3*scn2;
  // chain blocks per launch: all of them unless the buffer would exceed the cap
  int nbc = nblk;
  const double gb1 = 8.0e-9*static_cast<double>(grp1)*nGrp;   // per chain block
  if (gb1*nblk > rt_ckscr_gb) {
    nbc = std::max(1, static_cast<int>(rt_ckscr_gb/gb1));
  }
  const size_t sctsz = static_cast<size_t>(nGrp)*RT_NB*scn1*CKS_W;
  const size_t scntl = (static_cast<size_t>(nmb1 + 1)*nbc*scn3*scn2 + CKS_W - 1)
                       /CKS_W;
  CkScrEnsure(scntl*sctsz);
  Real *scbuf = rt_ckscr_ptr->data();
  auto chain_body =
  KOKKOS_LAMBDA(const int m, const int blk, const int k, const int j) {
    constexpr int NC = RT_NB;
#if defined(KOKKOS_ENABLE_CUDA)
    // nvcc: an extended lambda may not capture a variable for the first time inside an
    // `if constexpr` branch.  Name up front the ones only such branches read (a
    // discarded expression: no code; other backends capture them as before).
    (void)ACC; (void)AFC; (void)BTF; (void)X1F; (void)cbt_; (void)ckdif_;
    (void)difg_g; (void)ckjac_g; (void)ckdb_g; (void)icut_g;
#endif
    // this thread's number in its launch (the par_for flattening of m, blk, k, j)
    const int scb0 = blk - blk%nbc;
    const int scnb = (nblk - scb0 < nbc) ? (nblk - scb0) : nbc;
    const int sct = ((m*scnb + (blk - scb0))*scn3 + (k - ks))*scn2 + (j - js);
    // the frozen-operator pass selectors.  With the tag off these are
    // compile-time false and every frozen branch below is dead code; see FOP.
    const bool ckfst = FOP && ckfst_;   // this pass STORES the operator
    const bool ckfus = FOP && ckfus_;   // this pass RE-APPLIES it
    const bool ckfcf = FOP && ckfcf_;   // ... including the coefficient triple
    // ck_impl_colskip: a converged column keeps the fluxes and the source its
    // last pass left, so the zeroing below must not run for it either
    if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
    for (int i=is; i<ie+2; ++i) {
      Fb_g(m,blk,i,k,j) = 0.0;
      // ck_impl_frozen_op: the direct beam is exactly temperature-independent at
      // frozen opacity, so a frozen pass keeps the deposit pass 0 left and never
      // re-runs the ray integration
      if (!ckfus) Qb_g(m,blk,i,k,j) = 0.0;
      Em_g(m,blk,i,k,j) = 0.0;
      Src_g(m,blk,i,k,j) = 0.0;
      // ck-determinism: this block's own slots of the tridiagonal (ck_jacq_ptr)
      if constexpr (JAC) {
        ckjac_g(m,3*blk,k,j,i) = 0.0;
        ckjac_g(m,3*blk+1,k,j,i) = 0.0;
        ckjac_g(m,3*blk+2,k,j,i) = 0.0;
      }
    }
    const int icut = icc_g(m,k,j);
    if (icut > ie) return;                  // whole column deeper than the cut
    // Shortwave. This is the one part of the scheme that genuinely restructures:
    // the longwave only ever needs a LAYER optical depth, which is local, but the
    // direct stellar beam needs the CUMULATIVE depth from the top, so each chain
    // carries its own downward recurrence and it cannot live in the per-column
    // rt_pre. It rides along in the longwave down-sweep because the two share the
    // same kappa lookup at every cell -- doing it in a second kernel would pay
    // for
    // that lookup twice.
    const Real mu0 = cf_g(m,k,j,3);
    const Real facsw = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
    const bool lit = (mu0 > 0.0);
    // problem/ck_beam_sph: no secant clamp and no dayside test -- the twilight
    // columns are lit too, and whether a given ray reaches the target is decided
    // per target by its own tangent radius.  sinz is sin(theta0).
    const Real sinz = sqrt((mu0*mu0 < 1.0) ? (1.0 - mu0*mu0) : 0.0);
    const bool lit_sph = true;
    Real tausw[NC];
    Real transw[NC];      // beam transmission at the face above the current cell

    int bandc[NC], gc[NC];
    Real muc[NC], wfc[NC], wgc[NC];
    for (int cc=0; cc<NC; ++cc) {
      const int c = blk*NC + cc;
      if (ck_nq_ == 1) {
        gc[cc] = c % CK_NG;
        bandc[cc] = c/CK_NG;
        muc[cc] = 1.0/CK_DIFFUSIVITY;
        wfc[cc] = M_PI*ckgw(gc[cc]);          // F = pi I
      } else {
        const int nq = c % 2;
        gc[cc] = (c/2) % CK_NG;
        bandc[cc] = c/(2*CK_NG);
        muc[cc] = mug[nq];
        wfc[cc] = 2.0*M_PI*wg[nq]*mug[nq]*ckgw(gc[cc]);
      }
      // the shortwave weights by the g-point alone: it is a direct beam, not an
      // angular quadrature, and with nquad = 2 each angular point would otherwise
      // double-count the incident flux
      // ck-nq2: with nquad = 2 the two angles of a (band, g) pair are
      // adjacent chains of the same block (c = 2 (band CK_NG + g) + q, RT_NB
      // even), and the beam depends on neither, so the q = 0 chain carries the
      // pair's whole weight and the q = 1 chain skips the beam entirely (tm
      // sweep; the other forms still run it, at weight 0).  nquad = 1: gw.
      wgc[cc] = (ck_nq_ == 1) ? ckgw(gc[cc])/static_cast<Real>(ck_nq_)
                              : ((cc & 1) ? 0.0 : ckgw(gc[cc]));
    }
    const bool bpair = (ck_nq_ == 2);   // skip the beam on odd chains
    const int bst = bpair ? 2 : 1;       // the beam's chain stride
    // A WARNING ABOUT MEASURING THIS KERNEL. Seven optimisations were measured
    // against it while its per-cell arrays were laid out (m,slot,k,j,i), which
    // put
    // adjacent lanes 544 bytes apart. All seven failed, and several of those
    // verdicts
    // were artefacts of that: with the wave starved on scattered loads, nothing
    // done
    // to the arithmetic could show up. Re-measured on the (m,slot,i,k,j) layout:
    //
    //                                     starved layout    coalesced layout
    //   FP32 recurrence                        +2.8 %            +1.42x
    //   RT_NB = 2 instead of 4                  -16 %      faster on rt_chain,
    //                                                      slower on the total
    //   dropping this private column            -22 %            neutral
    //   RT_NB = 8                              slower            slower
    //
    // So: do not trust a null result on this kernel without checking that memory
    // is
    // not the thing in the way. Still genuinely useless, both layouts: the
    // k-table
    // layout (0.5 %), blocking the lookup by band (0.4 %), and precomputing kappa
    // per (cell, chain) (-1 %, because those four table loads are cache hits).
    //
    // Storing the intensity column costs 2304 bytes of scratch per thread, and
    // accumulating each sweep's flux contribution separately instead would remove
    // it entirely. On the starved layout that cost 22 % (1585 -> 1927 ms),
    // because
    // the extra flux traffic it adds was uncoalesced. On the fixed layout it is
    // neutral (454 vs 457 ms), which confirms the diagnosis. Neutral is not a
    // reason to change it, so the column stays.  (ck-scratch: it, and the other
    // whole-column arrays below, now live in rt_ckscr_ptr, not on the stack.)
    auto I_down = CkScrGetI<RtF, true, nScS>(scbuf, gI, sctsz, sct);
    // THE FACE MIXING, STORED.  Conservation needs the two rays to use the SAME
    // c at a face -- that, and only that, is what makes A_below (u_b - d_b) =
    // A_above (u_a - d_a) hold and the deposit telescope.  The down-sweep forms
    // c from its exact d_above and the probe's u_below and parks it here; the
    // up-sweep REUSES it instead of re-forming it from its own u_below, which
    // would be a different number and would break the telescoping at O(beta) per
    // face (measured: 0.7 % per cell, 19 % over the production column).  The
    // accuracy of c is then the probe's; its CONSISTENCY is exact.
    auto Cmx = CkScrGetI<RtF, SPH, nScS>(scbuf, gC, sctsz, sct);
    // problem/ck_beam_sph: the (kappa rho) column, filled by the down-sweep and
    // read by the ray integration.  tau_ray is NOT a running sum -- every target
    // radius has its own chord set -- so the profile has to be kept.  One element
    // when the switch is off.
    auto Krs = CkScrGetI<RtF, BSP, nScS>(scbuf, gK, sctsz, sct);

    // Top: the column above the domain, using the top cell's opacity over the
    // hydrostatic column p/g -- the same construction the grey scheme uses.
    {
      const Real Ttop = T_g(m,k,j,ie+1);
      const Real ptop = pb_g(m,k,j,ie+1);
      int iT, iP;
      Real fT, fP;
      const Real xTv = xT_g(m,k,j,ie+1);
      const Real xPv = xP_g(m,k,j,ie+1);
      iT = static_cast<int>(xTv); fT = xTv - static_cast<Real>(iT);
      iP = static_cast<int>(xPv); fP = xPv - static_cast<Real>(iP);
      for (int cc=0; cc<NC; ++cc) {
        // ck_impl_frozen_op: the top layer contributes (1 - e^-dtau) B_ghost, and
        // the factor depends on the opacity alone
        // problem/ck_sph_top: 1 = nothing, 2 = the top active cell's kappa and B
        const int itb = (tpm == 2) ? ie : ie+1;
        if (ckfus) {
          I_down[cc][ie+1] = static_cast<RtF>(cktpf_g(m,blk*NC+cc,k,j))
                           * static_cast<RtF>(Bb_g(m,bandc[cc],itb,k,j));
          if (tpm == 1) I_down[cc][ie+1] = static_cast<RtF>(0.0);
          tausw[cc] = 0.0;
          transw[cc] = static_cast<RtF>(1.0);
          continue;
        }
        const Real kap = ck_kappa(cklk, iT, fT, iP, fP, bandc[cc], gc[cc])
                       + kc_g(m,bandc[cc],ie+1,k,j);
        const Real dtau = RTTopDtau(kap, ptop*1.0e6,
                                    EffGravAt(grav, ap, x1v_(m,ie+1),
                                              grav_pmass, omega, mu0, tide)
                                    *cf_g(m,k,j,4));
        RtF trans = RT_EXP(-static_cast<RtF>(dtau/muc[cc]));
        if (tpm == 2) {
          const Real xTa = xT_g(m,k,j,ie);
          const Real xPa = xP_g(m,k,j,ie);
          const int iTa = static_cast<int>(xTa);
          const int iPa = static_cast<int>(xPa);
          const Real kpa = ck_kappa(cklk, iTa, xTa - static_cast<Real>(iTa), iPa,
                                    xPa - static_cast<Real>(iPa), bandc[cc], gc[cc])
                         + kc_g(m,bandc[cc],ie,k,j);
          const Real dta = RTTopDtau(kpa, ptop*1.0e6,
                                     EffGravAt(grav, ap, x1v_(m,ie+1),
                                               grav_pmass, omega, mu0, tide)
                                     *cf_g(m,k,j,4));
          trans = RT_EXP(-static_cast<RtF>(dta/muc[cc]));
        }
        if (ckfst) {
          cktpf_g(m,blk*NC+cc,k,j) =
              static_cast<Real>(static_cast<RtF>(1.0)-trans);
          if (tpm == 1) cktpf_g(m,blk*NC+cc,k,j) = 0.0;   // ck_sph_top = 1: no datum
        }
        I_down[cc][ie+1] = (static_cast<RtF>(1.0)-trans)
                         * static_cast<RtF>(Bb_g(m,bandc[cc],itb,k,j));
        if (tpm == 1) I_down[cc][ie+1] = static_cast<RtF>(0.0);
        // BSP: nothing above the domain top, so the ray enters unattenuated;
        // the plane-parallel path instead charges it the ghost column's depth
        tausw[cc] = BSP ? 0.0 : dtau;
        transw[cc] = BSP ? static_cast<RtF>(1.0)
                         : RT_EXP(-static_cast<RtF>(dtau*facsw));
      }
    }
    if (layer_legacy) {
      // down-sweep
      for (int i=ie; i>icut-1; --i) {
        const Real rho = rhoN(m,k,j,i);
        const Real drho = rho*dx1(m,k,j,i);
        int iT, iP;
        Real fT, fP;
        const Real xTv = xT_g(m,k,j,i);
        const Real xPv = xP_g(m,k,j,i);
        iT = static_cast<int>(xTv); fT = xTv - static_cast<Real>(iT);
        iP = static_cast<int>(xPv); fP = xPv - static_cast<Real>(iP);
        // the far endpoint of this layer's source function, weighted by emitting
        // matter (see BFace).  The neighbour's kappa costs a second table lookup,
        // so it is paid only where the two Planck functions differ by more than
        // 4x
        // -- a smooth column takes the old expression and is bit-identical.
        const int ir1 = (i+1 > ie) ? ie : i+1;
        const Real rhf = rhoN(m,k,j,ir1);
        const Real xTf = xT_g(m,k,j,i+1);
        const Real xPf = xP_g(m,k,j,i+1);
        const int iTf = static_cast<int>(xTf);
        const int iPf = static_cast<int>(xPf);
        const Real fTf = xTf - static_cast<Real>(iTf);
        const Real fPf = xPf - static_cast<Real>(iPf);
        for (int cc=0; cc<NC; ++cc) {
          const int b = bandc[cc];
          const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                           + kc_g(m,b,i,k,j);
          const Real bown = Bb_g(m,b,i,k,j);
          Real bfar = Bb_g(m,b,i+1,k,j);
          if (bfar > 4.0*bown || bown > 4.0*bfar) {
            const Real kapf = ck_kappa(cklk, iTf, fTf, iPf, fPf, b, gc[cc])
                            + kc_g(m,b,i+1,k,j);
            bfar = BFace(kap*rho, kapf*rhf, bown, bfar, bface_on);
          }
          const RtF x = static_cast<RtF>(kap*drho/muc[cc]);
          const RtF e0 = -RT_EXPM1(-x);
          const RtF one = static_cast<RtF>(1.0);
          const RtF alp = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                                        : (x/2 - x*x/3);
          const RtF bet = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                        : (x/2 - x*x/6);
          Src_g(m,blk,i,k,j) += wfc[cc]/dx1(m,k,j,i)
              *static_cast<Real>(e0*I_down[cc][i+1]
                - (alp*static_cast<RtF>(bfar)
                   + bet*static_cast<RtF>(bown)));
          I_down[cc][i] = (one-e0)*I_down[cc][i+1]
                        + alp*static_cast<RtF>(bfar)
                        + bet*static_cast<RtF>(bown);
          // Direct beam. Deposit the flux DIFFERENCE across the cell, not
          // kappa rho F exp(-tau) evaluated at one face. The latter is what the
          // grey
          // picket fence does, and it under-deposits badly once a layer is not
          // thin:
          // the ratio of deposited to absorbed is u e^-u/(1 - e^-u) with u =
          // dtau/mu,
          // which is 0.95 at u = 0.1 but 0.58 at u = 1 and 0.31 at u = 2. Summed
          // down
          // a column with u ~ 0.5 it loses a quarter of the incident flux --
          // measured
          // against Exo-FMS on an identical column, which is how this was found.
          // Written this way the column integral is F mu (1 - e^-tau_total) by
          // construction, and it still reduces to the old form as dtau -> 0.
          tausw[cc] += kap*drho;
          if (lit) {
            const Real tnew = RT_EXP(-static_cast<RtF>(tausw[cc]*facsw));
            Qb_g(m,blk,i,k,j) += (1.0-albedo)*Fstar*ckswf(b)*wgc[cc]/facsw
                      * (transw[cc] - tnew)/dx1(m,k,j,i);
            transw[cc] = tnew;
          }
        }
      }
      // Bottom of the CORRELATED-K DOMAIN, not of the column. At the cut the grey
      // optical depth is of order 1e4, so the layer is thermalised to e^-tau and
      // the
      // upward intensity is its own Planck function; the planet's internal flux
      // is
      // delivered here as an extra band-weighted source. Below the cut nothing
      // radiative is applied -- that region is optically thick and convective,
      // and
      // the flux simply passes through it.
      RtF I_up[NC];
      for (int cc=0; cc<NC; ++cc) {
        const int b = bandc[cc];
        const Real Iint_b = (int_at_cut ? boltz_sigma/M_PI*Tint4
                          * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
        I_up[cc] = static_cast<RtF>(Bb_g(m,b,icut,k,j) + Iint_b);
        Fb_g(m,blk,icut,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][icut]);
      }
      // up-sweep
      for (int i=icut+1; i<ie+2; ++i) {
        const Real rho = rhoN(m,k,j,i-1);
        const Real drho = rho*dx1(m,k,j,i-1);
        int iT, iP;
        Real fT, fP;
        const Real xTv = xT_g(m,k,j,i-1);
        const Real xPv = xP_g(m,k,j,i-1);
        iT = static_cast<int>(xTv); fT = xTv - static_cast<Real>(iT);
        iP = static_cast<int>(xPv); fP = xPv - static_cast<Real>(iP);
        // same emissivity-weighted far endpoint, and the same 4x guard: see BFace
        const int ir1 = (i > ie) ? ie : i;
        const Real rhf = rhoN(m,k,j,ir1);
        const Real xTf = xT_g(m,k,j,i);
        const Real xPf = xP_g(m,k,j,i);
        const int iTf = static_cast<int>(xTf);
        const int iPf = static_cast<int>(xPf);
        const Real fTf = xTf - static_cast<Real>(iTf);
        const Real fPf = xPf - static_cast<Real>(iPf);
        for (int cc=0; cc<NC; ++cc) {
          const int b = bandc[cc];
          const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                         + kc_g(m,b,i-1,k,j);
          const Real bown = Bb_g(m,b,i-1,k,j);
          Real bfar = Bb_g(m,b,i,k,j);
          if (bfar > 4.0*bown || bown > 4.0*bfar) {
            const Real kapf = ck_kappa(cklk, iTf, fTf, iPf, fPf, b, gc[cc])
                            + kc_g(m,b,i,k,j);
            bfar = BFace(kap*rho, kapf*rhf, bown, bfar, bface_on);
          }
          const RtF x = static_cast<RtF>(kap*drho/muc[cc]);
          const RtF e0 = -RT_EXPM1(-x);
          const RtF one = static_cast<RtF>(1.0);
          const RtF bet = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                        : (x/2 - x*x/6);
          const RtF gm = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                                       : (x/2 - x*x/3);
          const RtF Iup_in = I_up[cc];
          Src_g(m,blk,i-1,k,j) += wfc[cc]/dx1(m,k,j,i-1)
              *static_cast<Real>(e0*Iup_in
                - (bet*static_cast<RtF>(bfar)
                   + gm*static_cast<RtF>(bown)));
          I_up[cc] = (one-e0)*Iup_in
                   + bet*static_cast<RtF>(bfar)
                   + gm*static_cast<RtF>(bown);
          Fb_g(m,blk,i,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][i]);
          // the cell's OWN emission, both hemispheres: in the thin limit each
          // stream adds wfc*(kap*rho*dr/mu)*B over the layer, so per unit volume
          // the two together give 2*(wfc/mu)*kap*rho*B.  Summed over bands and g
          // this is 2*CK_DIFFUSIVITY*sigma*kappa_P*rho*T^4, i.e. the exact
          // 4 sigma kappa_P rho T^4 with 1.66 in place of 2 for the hemispheric
          // mean.  17 % low, which is well inside what a rate estimate needs.
          Em_g(m,blk,i-1,k,j) += 2.0*(wfc[cc]/muc[cc])*kap*rho
                               * 0.5*(bfar + bown);
        }
      }
    } else {
      // ---- the centre-to-centre layers (see rt_layer_legacy) ----
      // Every layer now runs from one cell CENTRE to the next, split at the face
      // into its two halves, so the interval the optical depth measures and the
      // interval the source runs over are the same one.  The column is closed by
      // two half layers, the upper half of cell ie and the lower half of cell
      // icut, each with its own cell's Planck function at both ends.
      //
      // NO EXTRA TABLE LOOKUP.  The far cell's kappa rho is CARRIED in kfar from
      // the previous iteration of the sweep, so this pays exactly one ck_kappa
      // per cell and chain -- one FEWER than the staggered layers, which looked
      // the neighbour up again whenever the two Planck functions differed by 4x.
      // kfar and the running intensity Idn are the only new per-chain state.
      Real kfar[NC];
      RtF Idn[NC];
      // ---- problem/ck_sweep_cache: the per-(cell, chain) operator, kept in the
      // thread instead of rebuilt.  See the flag's note for why this is the
      // whole of the redundancy and for what each mode costs.
      //
      // mode >= 1: cry0/cryi/cryo carry the coefficient triple of the half layer
      // the NEXT iteration of this pass is about to cross again -- REGISTERS,
      // NC deep, exactly like kfar, so the private segment does not grow.
      // mode >= 2: Kpc additionally holds kappa (the k-table value plus kc_g,
      // BEFORE the multiplication by rho, so that the beam's kappa and every
      // kappa rho are still formed by the expression that formed them) for the
      // whole column, which is what lets the three later passes skip the table.
      // That one costs NC x n1 of column storage (in rt_ckscr_ptr since
      // ck-scratch) and is a one-element private array otherwise.
      RtF cry0[NC], cryi[NC], cryo[NC];
      constexpr bool CC2 = (CCH >= 2);
      auto Kpc = CkScrGetI<Real, CC2, nScS>(scbuf, gP, sctsz, sct);
      auto Cc0 = CkScrGetI<RtF, CC2, nScS>(scbuf, gP+1, sctsz, sct);
      auto Cci = CkScrGetI<RtF, CC2, nScS>(scbuf, gP+2, sctsz, sct);
      auto Cco = CkScrGetI<RtF, CC2, nScS>(scbuf, gP+3, sctsz, sct);
      // one half-layer step, in the chain's own precision: the same exponential
      // coefficients the staggered layers used, handed the half interval.  dsrc
      // comes back as absorbed minus emitted, which is what Src_g wants and is
      // exactly the change in the intensity.
      //
      // ---- problem/ck_impl_frozen_op ------------------------------------------
      // EVERY half layer in this kernel is the half of ONE cell taken at ONE
      // chain's mu -- the two probe passes, the down-sweep and the up-sweep all
      // step across the same six half layers per cell pair -- so one coefficient
      // triple per (cell, chain) describes the whole operator.  A frozen pass
      // loads it instead of taking the exponential; a storing pass computes it
      // exactly as it always did and writes it down.  The expression is otherwise
      // untouched, so with the switch off nothing here changes.
      // ck_sweep_cache: two flags, because the two modes reuse over different
      // spans.  `first` marks the ONE call in the whole COLUMN that first visits
      // this half layer (mode 2 stores there and loads everywhere else);
      // `fill` marks the one call per PASS that computes it fresh, the other
      // call in that pass crossing the same half layer one iteration later and
      // reading the carry (mode 1).  With the switch off both are dead and every
      // call computes, as it always did.
      auto step = [&](const int cc, const int ic, const Real dtau,
                      const Real s_in, const Real s_out, RtF &I, Real &dsrc,
                      const bool fill, const bool first) {
        RtF e0, cin, cout;
        if (ckfus && ckfcf) {
          e0 = static_cast<RtF>(ckc0_g(m,blk*NC+cc,ic,k,j));
          cin = static_cast<RtF>(ckci_g(m,blk*NC+cc,ic,k,j));
          cout = static_cast<RtF>(ckco_g(m,blk*NC+cc,ic,k,j));
        } else if (CCH >= 2 && !first) {
          e0 = Cc0[cc][ic];
          cin = Cci[cc][ic];
          cout = Cco[cc][ic];
        } else if (CCH == 1 && !fill) {
          e0 = cry0[cc];
          cin = cryi[cc];
          cout = cryo[cc];
        } else {
          const RtF x = static_cast<RtF>(dtau/muc[cc]);
          e0 = -RT_EXPM1(-x);
          const RtF one = static_cast<RtF>(1.0);
          cin = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                               : (x/2 - x*x/3);
          cout = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                : (x/2 - x*x/6);
          if (ckfst && ckfcf) {
            ckc0_g(m,blk*NC+cc,ic,k,j) = static_cast<Real>(e0);
            ckci_g(m,blk*NC+cc,ic,k,j) = static_cast<Real>(cin);
            ckco_g(m,blk*NC+cc,ic,k,j) = static_cast<Real>(cout);
          }
          if (CCH >= 2) {
            Cc0[cc][ic] = e0;
            Cci[cc][ic] = cin;
            Cco[cc][ic] = cout;
          } else if (CCH == 1) {
            cry0[cc] = e0;
            cryi[cc] = cin;
            cryo[cc] = cout;
          }
        }
        const RtF em = cin*static_cast<RtF>(s_in)
                     + cout*static_cast<RtF>(s_out);
        dsrc = static_cast<Real>(e0*I - em);
        I = (static_cast<RtF>(1.0) - e0)*I + em;
      };
      // ---- problem/ck_sweep_form: THE SAME TRIPLE, HANDED BACK INSTEAD OF
      // APPLIED.  The probe-free forms do not push an intensity through a half
      // layer in pass 1 -- they compose the layer's 2x2 map -- so they need
      // (e0, cin, cout) as numbers.  The branch order below is `step`'s, slot
      // for slot, so that the two passes and the four passes read the very same
      // cache entries and produce the very same coefficients; a change to the
      // cache logic in `step` has to be repeated here.
      auto cofs = [&](const int cc, const int ic, const Real dtau,
                      RtF &e0, RtF &cin, RtF &cout,
                      const bool fill, const bool first) {
        if (ckfus && ckfcf) {
          e0 = static_cast<RtF>(ckc0_g(m,blk*NC+cc,ic,k,j));
          cin = static_cast<RtF>(ckci_g(m,blk*NC+cc,ic,k,j));
          cout = static_cast<RtF>(ckco_g(m,blk*NC+cc,ic,k,j));
        } else if (CCH >= 2 && !first) {
          e0 = Cc0[cc][ic];
          cin = Cci[cc][ic];
          cout = Cco[cc][ic];
        } else if (CCH == 1 && !fill) {
          e0 = cry0[cc];
          cin = cryi[cc];
          cout = cryo[cc];
        } else {
          const RtF x = static_cast<RtF>(dtau/muc[cc]);
          e0 = -RT_EXPM1(-x);
          const RtF one = static_cast<RtF>(1.0);
          cin = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                               : (x/2 - x*x/3);
          cout = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                : (x/2 - x*x/6);
          if (ckfst && ckfcf) {
            ckc0_g(m,blk*NC+cc,ic,k,j) = static_cast<Real>(e0);
            ckci_g(m,blk*NC+cc,ic,k,j) = static_cast<Real>(cin);
            ckco_g(m,blk*NC+cc,ic,k,j) = static_cast<Real>(cout);
          }
          if (CCH >= 2) {
            Cc0[cc][ic] = e0;
            Cci[cc][ic] = cin;
            Cco[cc][ic] = cout;
          } else if (CCH == 1) {
            cry0[cc] = e0;
            cryi[cc] = cin;
            cryo[cc] = cout;
          }
        }
      };
      // ck_impl_frozen_op: kappa rho of cell ic, chain cc -- loaded on a frozen
      // pass (which is what removes the correlated-k table look-up, the dominant
      // cost of a pass), computed and stored otherwise.
      // ck_sweep_cache: `first` as above -- in mode 2 the column's first pass
      // pays the k-table look-up and the other three read Kpc.  kappa is cached,
      // not kappa rho, so the product is formed where it always was.  Mode 1
      // keeps every pass's own look-up (one per cell, carried in kdp/kpr/kfar)
      // because the column-wide store is the part that costs private memory.
      auto kapof = [&](const int cc, const int ic, const int b, const int iT,
                       const Real fT, const int iP, const Real fP,
                       const bool first) {
        if (CCH >= 2 && !first) return Kpc[cc][ic];
        // ck-nq2: the q = 1 chain of a (band, g) pair has the SAME kappa as the
        // q = 0 chain the same loop filled just before it (every first-pass
        // call site walks cc upward at a fixed cell), so it is copied, not
        // looked up again: bitwise the same value, half the table traffic
        if (CCH >= 2 && bpair && (cc & 1)) {
          const Real kq = Kpc[cc-1][ic];
          Kpc[cc][ic] = kq;
          return kq;
        }
        const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                       + kc_g(m,b,ic,k,j);
        if (CCH >= 2) Kpc[cc][ic] = kap;
        return kap;
      };
      auto krof = [&](const int cc, const int ic, const int b, const int iT,
                      const Real fT, const int iP, const Real fP,
                      const Real rho, const bool fill) {
        if (ckfus) return ckkro_g(m,blk*NC+cc,ic,k,j);
        const Real kr = kapof(cc, ic, b, iT, fT, iP, fP, fill)*rho;
        if (ckfst) ckkro_g(m,blk*NC+cc,ic,k,j) = kr;
        return kr;
      };
      // ---- problem/ck_implicit: THE SAME HALF LAYER, DIFFERENTIATED ---------
      // The three coefficients of `step`, recomputed in Real (the Jacobian does
      // not need the chain's working precision) and OUTSIDE the sweep's own
      // expression DAG, so that with the switch off not one instruction of the
      // sweep changes.  A half layer costs one extra expm1 when the switch is
      // on; that, and the two per-face BFaceW calls below, are the whole cost of
      // the assembly.
      //
      // The layer source is a convex combination of the two centre Planck
      // functions (BFace) and of the two half-layer endpoints (the dt_l/dtc face
      // interpolation), so six weights describe both halves completely:
      //   ds_l/dB_l = pl, ds_l/dB_u = 1-pl, ds_u/dB_u = pu, ds_u/dB_l = 1-pu,
      //   ds_f/dB_l = df_l = (1-f) pl + f (1-pu),  df_u = 1 - df_l.
      // dsrc/dB = e0 dI_in/dB - dem/dB and dI_out/dB = (1-e0) dI_in/dB + dem/dB,
      // and the carried dI/dB is truncated to the LAST cell the ray crossed --
      // everything further is attenuated by (1-e0) and belongs outside a
      // tridiagonal.  See utils/two_stream_column_ck.hpp.
      // problem/ck_impl_reuse_jac: a reusing pass assembles NOTHING -- it only
      // re-forms the residual, which is the sweep it was going to run anyway.
      const bool jck = ckjacp_;
      // ck_impl_frozen_op: in double precision (the default) the stored triple IS
      // this triple, so a frozen pass loads it here too and the assembly costs no
      // exponential either.  Under RT_FP32 the sweep's coefficients are floats
      // and the Jacobian's are not, so the store is not usable and the triple is
      // recomputed -- which is what it always did.
      const bool jfrz = ckfus && ckfcf && (RT_FP32 == 0);
      auto jcof = [&](const int cc, const int ic, const Real dtau, const Real mu,
                      Real &e0, Real &cin, Real &cout) {
        if (jfrz) {
          e0 = ckc0_g(m,blk*NC+cc,ic,k,j);
          cin = ckci_g(m,blk*NC+cc,ic,k,j);
          cout = ckco_g(m,blk*NC+cc,ic,k,j);
          return;
        }
        const Real x = dtau/mu;
        e0 = -expm1(-x);
        cin = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - x*x/3.0);
        cout = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - x*x/6.0);
      };
      // the running dI/dB of the cell the ray has just left, down ray and up ray
      Real jdn[NC], jup[NC];
      if (jck) {
        for (int cc=0; cc<NC; ++cc) {
          jdn[cc] = 0.0;
          jup[cc] = 0.0;
        }
      }
      // ================ problem/ck_sweep_form: THE PROBE-FREE FORMS =========
      // Two passes, exact.  PASS 1 climbs from the cut face carrying the linear
      // relation between the two rays -- (R, Sc) in the transfer-matrix form,
      // (G, H) in the S/D one -- and parks it at every face; PASS 2 comes back
      // down, steps the down ray, closes the up ray on the stored relation, and
      // writes the face fluxes and the deposit.  See the flag's note for the
      // recurrences, the boundary data and why the Riccati variable is the
      // stable way to say "sweep in S and D".
      //
      // NO NEW PRIVATE MEMORY.  The probe form needs I_down (the down column)
      // and Cmx (the face mixing); this one needs neither -- pass 2 has the down
      // intensity and the up intensity at the same face in the same iteration --
      // so the two columns are BORROWED to hold the Riccati pair instead.  The
      // thread's private segment is the same size it was.
      // nvcc: an extended lambda may not FIRST capture a variable inside an
      // if-constexpr branch, so name every one this branch alone uses here.
      (void)fbl; (void)vef; (void)vfg_g; (void)vfs_g; (void)wex; (void)dil; (void)pth;
      if constexpr (FRM != 0) {
        const RtF one = static_cast<RtF>(1.0);
        const RtF two = static_cast<RtF>(2.0);
        // the top face's incoming datum, taken out before the column is borrowed
        RtF Itp[NC];
        for (int cc=0; cc<NC; ++cc) {
          Itp[cc] = I_down[cc][ie+1];
        }
        // per-chain carries: the cell's kappa rho, the source at the face below
        // it and the source at its centre on the lower half's side
        Real kcu[NC], sfc[NC], suc[NC];
        // ck_sph_face = 2: the radial optical depth of the two half cells that meet
        // at the NEXT face up (the face the pass-1 loop crosses first), and the
        // blend map of a face.  e is the Eddington weight; see ck_sph_face.
        Real dtf[NC];
        auto fblm = [&](const int cc, const Real rj, const Real bt, const Real dtc,
                        Real &al, Real &be, Real &c1, Real &r1) {
          const Real we = -expm1(-dtc/muc[cc]);
          const Real rh = (1.0 - bt)/(1.0 + bt);
          const Real dc = 1.0 + we*rj + we*rh*(1.0 - rj);
          al = (1.0 + we)/dc;
          be = we*(1.0 - rh)/dc;
          c1 = rh*al;
          r1 = 1.0 - rh*(1.0 - rj)*al;
        };
        // ck_sph_face = 3, 4: R just above a face (the value map, bitwise the one
        // pass 1 applies), and the up ray's half-cell transmission with the cone
        // path factor f = 1/(1 + cos th_c), cos^2 th_c = T^4 R (see ck_sph_face)
        auto rfm = [&](const int cc, const RtF rr, const Real bt, const Real dtc) {
          if (fbl) {
            Real alx, bex, c1x, r1x;
            fblm(cc, static_cast<Real>(rr), bt, dtc, alx, bex, c1x, r1x);
            return static_cast<RtF>(r1x);
          }
          const RtF bb = static_cast<RtF>(bt);
          const RtF rh = (one - bb)/(one + bb);
          return one - rh*(one - rr);
        };
        auto upt = [&](const int cc, const Real rf, const Real tr, const Real dth) {
          const Real t2 = tr*tr;
          const Real cq = t2*t2*rf;
          const Real cth = sqrt((cq > 0.0) ? cq : 0.0);
          return exp(-dth/(muc[cc]*(1.0 + cth)));
        };
        RtF Rc[NC], Hc[NC];
        // ---- problem/ck_implicit (JAC): THE TRIDIAGONAL OF THE tm SWEEP -------
        // At frozen opacity the tm sweep is affine in the band Planck functions
        // and R does not depend on B at all, so the Jacobian only has to carry
        // the Sc relation and the rays.  The deposit of cell i is exactly
        //   Src_i = (wfc/dz_i) [rat_i D_i - D_{i+1}],
        //   D_f = u_f - d_f below face f = (1+b)(Sc_f - (1-R_f) d_f^+)/(1+R_f b),
        // d_f^+ the down ray arriving at face f from above (the up rays and the
        // closure telescope out of the four half-layer steps).  So the three
        // entries of row i need dD_f/dB_j for j = f-2 .. f+1 at f = i, i+1:
        //   PASS 1 carries dSc_f/dB_{f-2,f-1,f} (the window jS) up the column
        //     through the same Moebius / half-layer maps as Sc itself, and parks
        //     it at every face (Js0..2);
        //   PASS 2 carries dd^+/dB_{i-1,i,i+1} (v) down the column, and ONE
        //     scalar for everything further below: any B_j with j <= f-2 reaches
        //     d_f^+ only through reflection of Sc_f, i.e. dd_f^+/dB_j =
        //     -W_f dSc_f/dB_j with W_f = T_f^2 c_f (alpha_{f+1} W_{f+1} +
        //     beta_{f+1}), c_f = T_f^2 (1-b_f)/(1+R_f b_f), alpha = (1+b)/(1+Rb),
        //     beta = b/(1+Rb), W = 0 at the top (the ghost datum is fixed).
        // The nearest-neighbour entries are therefore EXACT derivatives of this
        // sweep (reflection included); entries two cells away or more stay in the
        // residual, as in the four-pass assembly.  The O(beta) reflection terms
        // can make an off-diagonal contribution negative; those are dropped (the
        // M-matrix the Thomas sweep relies on), which makes it quasi-Newton with
        // the same root.  All of it is compiled only into the JAC instantiation.
        constexpr int NJ = JAC ? NC : 1;
        auto Js0 = CkScrGetI<Real, JAC, nScS>(scbuf, gJ, sctsz, sct);
        auto Js1 = CkScrGetI<Real, JAC, nScS>(scbuf, gJ+1, sctsz, sct);
        auto Js2 = CkScrGetI<Real, JAC, nScS>(scbuf, gJ+2, sctsz, sct);
        Real jS[NJ][3], jfc[NJ][2], juc[NJ][2];
        // ---- PASS 1: upward, accumulating the relation at every face --------
        {
          const Real rho = rhoN(m,k,j,icut);
          const Real xTv = xT_g(m,k,j,icut);
          const Real xPv = xP_g(m,k,j,icut);
          const int iT = static_cast<int>(xTv);
          const int iP = static_cast<int>(xPv);
          const Real fT = xTv - static_cast<Real>(iT);
          const Real fP = xPv - static_cast<Real>(iP);
          for (int cc=0; cc<NC; ++cc) {
            const int b = bandc[cc];
            kcu[cc] = krof(cc, icut, b, iT, fT, iP, fP, rho, true);
            const Real bcut = Bb_g(m,b,icut,k,j);
            // the bottom datum: u_below is given, i.e. no reflection at all
            const Real Iint_b = (int_at_cut ? boltz_sigma/M_PI*Tint4
                            * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
            sfc[cc] = bcut;
            suc[cc] = bcut;
            dtf[cc] = 0.5*kcu[cc]*dx1(m,k,j,icut);
            Rc[cc] = (FRM == 2) ? one : static_cast<RtF>(0.0);
            Hc[cc] = static_cast<RtF>(bcut + Iint_b);
            // ck_dif_dtau: a handover column starts on the deep diffusion
            // flux instead, u - d = g_c (B_{cut-1} - B_cut) (R = 1)
            const Real gdf = (FRM == 1 && ckdif_)
                           ? difg_g(m,blk*NC+cc,k,j) : -1.0;
            if (gdf >= 0.0) {
              Rc[cc] = one;
              Hc[cc] = static_cast<RtF>(gdf*(Bb_g(m,b,icut-1,k,j) - bcut));
            }
            // ck_sph_face = 5: the black wall as a Marshak-type condition on the
            // moments, J + beta H = B (1/2 + beta/4), beta = J^-/H^- of the formal
            // solution's incoming field (2 for an isotropic one): the outgoing
            // half-range flux is B/4 whatever comes in.  With the pair at mu_eff:
            // u = R0 d + S0, R0 = -(1 - beta mu)/(1 + beta mu),
            // S0 = (B + Iint)(1 + beta/2)/(1 + beta mu)
            Real vs0 = 1.0;
            if (vef) {
              const int cv = b*CK_NG + gc[cc];
              const Real bw = vfg_g(m,cv,icut-1,k,j);
              const Real mw = vfs_g(m,cv,icut,k,j)/SQRT3_VEF;
              Rc[cc] = static_cast<RtF>(-(1.0 - bw*mw)/(1.0 + bw*mw));
              vs0 = (1.0 + 0.5*bw)/(1.0 + bw*mw);
              // ck_wall_flux_exact: I_int scaled so the wall's net flux is pi I_int
              const Real iwb = wex ? Iint_b*(1.0 + bw*mw)/(2.0*mw*(1.0 + 0.5*bw))
                                   : Iint_b;
              Hc[cc] = static_cast<RtF>((bcut + iwb)*vs0);
            }
            if constexpr (JAC) {
              // Sc_cut = B_cut + Iint: window (B_cut-2, B_cut-1, B_cut)
              jS[cc][0] = 0.0;
              jS[cc][1] = (gdf >= 0.0) ? gdf : 0.0;
              jS[cc][2] = (gdf >= 0.0) ? -gdf : (vef ? vs0 : 1.0);
              // sfc = suc = B_cut: slots (B_{i-1}, B_i) of the first layer
              jfc[cc][0] = 0.0;
              jfc[cc][1] = 1.0;
              juc[cc][0] = 0.0;
              juc[cc][1] = 1.0;
            }
          }
        }
        for (int i=icut; i<ie+1; ++i) {
          const Real dz = dx1(m,k,j,i);
          const int iu = (i < ie) ? i+1 : i;
          const Real dzu = dx1(m,k,j,iu);
          const Real rhou = rhoN(m,k,j,iu);
          const Real xTu = xT_g(m,k,j,iu);
          const Real xPu = xP_g(m,k,j,iu);
          const int iTu = static_cast<int>(xTu);
          const int iPu = static_cast<int>(xPu);
          const Real fTu = xTu - static_cast<Real>(iTu);
          const Real fPu = xPu - static_cast<Real>(iPu);
          const Real bt = BTF(m,k,j,i,icut);
          const Real rat = (i == icut)
              ? AFC(m,k,j,icut)/ACC(m,k,j,icut)
              : ACC(m,k,j,i-1)/ACC(m,k,j,i);
          for (int cc=0; cc<NC; ++cc) {
            const int b = bandc[cc];
            const Real kro = kcu[cc];
            const Real bown = Bb_g(m,b,i,k,j);
            // the layer joining cells i and i+1: its endpoint at i's centre
            // (slv), at (i+1)'s centre (suu) and at the face between (sfv).  The
            // top half layer of cell ie has no neighbour above and keeps its own
            // Planck function at both ends, exactly as the four-pass form does.
            Real slv = bown, sfv = bown, suu = bown, kru = kro;
            // the face depth this cell's face map uses; the next face's, for the
            // top face the top half cell alone (as pass 2 forms them)
            const Real dtfc = dtf[cc];
            Real dtn = 0.5*kro*dz;
            if (i < ie) {
              kru = krof(cc, i+1, b, iTu, fTu, iPu, fPu, rhou, true);
              const Real bfar = Bb_g(m,b,i+1,k,j);
              slv = BFace(kru, kro, bfar, bown, bface_on);
              suu = BFace(kro, kru, bown, bfar, bface_on);
              const Real dt_l = 0.5*kro*dz;
              const Real dt_u = 0.5*kru*dzu;
              const Real dtc = dt_l + dt_u;
              sfv = (dtc > 0.0) ? (slv + (suu - slv)*(dt_l/dtc))
                                : (0.5*(slv + suu));
              dtn = dtc;
            }
            // both halves of a cell have the same optical thickness, so ONE
            // triple describes the pair (this is the ck_impl_frozen_op note)
            // ck_sph_face = 5: s below face i (vsb) and in cell i (vsa), gamma
            Real vsb = 1.0, vsa = 1.0, vgm = 1.0;
            if (vef) {
              const int cv = b*CK_NG + gc[cc];
              vsb = vfs_g(m,cv,i,k,j);
              vsa = vfs_g(m,cv,i+1,k,j);
              vgm = vfg_g(m,cv,i,k,j);
            }
            RtF e0, cin, cout;
            // (ck_sph_face = 5: BOTH chains of the pair at mu_eff = sqrt(f) = s/sqrt 3,
            // i.e. dtau/mu_eff = (dtau/mu) sqrt(3) mu/s: the pair is then exactly
            // the moment system with K = f J inside the cell)
            cofs(cc, i, vef ? 0.5*kro*dz*(SQRT3_VEF*muc[cc])/vsa : 0.5*kro*dz, e0, cin,
                 cout, true, true);
            const RtF tr = one - e0;
            // the emission along an up crossing (p) and a down crossing (q) of
            // the lower half (l) and the upper half (u) of this cell
            const RtF pl = cin*static_cast<RtF>(sfc[cc])
                         + cout*static_cast<RtF>(suc[cc]);
            const RtF ql = cin*static_cast<RtF>(suc[cc])
                         + cout*static_cast<RtF>(sfc[cc]);
            const RtF pu = cin*static_cast<RtF>(slv)
                         + cout*static_cast<RtF>(sfv);
            const RtF qu = cin*static_cast<RtF>(sfv)
                         + cout*static_cast<RtF>(slv);
            // park the relation on the BELOW side of face i, where pass 2 wants
            // it, and only then cross the face
            I_down[cc][i] = Rc[cc];
            Cmx[cc][i] = Hc[cc];
            if constexpr (JAC) {
              // park dSc_i/dB_{i-2,i-1,i}, then carry it across face i and the
              // two halves of cell i: slots below are (B_{i-1}, B_i, B_{i+1})
              Js0[cc][i] = jS[cc][0];
              Js1[cc][i] = jS[cc][1];
              Js2[cc][i] = jS[cc][2];
              Real dslv[3] = {0.0, 1.0, 0.0};
              Real dsuu[3] = {0.0, 1.0, 0.0};
              Real dsfw[3] = {0.0, 1.0, 0.0};
              if (i < ie) {
                const Real wl = BFaceW(kru, kro, bface_on);   // dslv/dB_i
                const Real wu = BFaceW(kro, kru, bface_on);   // dsuu/dB_{i+1}
                const Real dtl = 0.5*kro*dz;
                const Real dtcj = dtl + 0.5*kru*dzu;
                const Real ffj = (dtcj > 0.0) ? (dtl/dtcj) : 0.5;
                const Real dfl = (1.0 - ffj)*wl + ffj*(1.0 - wu);
                dslv[1] = wl;
                dslv[2] = 1.0 - wl;
                dsuu[1] = 1.0 - wu;
                dsuu[2] = wu;
                dsfw[1] = dfl;
                dsfw[2] = 1.0 - dfl;
              }
              const Real ci = static_cast<Real>(cin);
              const Real co = static_cast<Real>(cout);
              const Real tj = 1.0 - static_cast<Real>(e0);
              const Real rj = static_cast<Real>(Rc[cc]);
              const Real dnj = 1.0 + rj*bt;
              const Real rhd = (1.0 - bt)/(1.0 + bt);   // ck_sph_dilute: A_b/A_a
              Real r1 = dil ? (1.0 - rhd*(1.0 - rj)) : (rj + bt)/dnj;
              Real c1 = dil ? rhd : (1.0 - bt)/dnj;
              if (fbl) {
                Real alx, bex;
                fblm(cc, rj, bt, dtfc, alx, bex, c1, r1);
              }
              if (vef) {
                const Real rpj = rhd*vsb/vsa;
                const Real evj = (1.0 + rj) + vgm*rpj*(1.0 - rj);
                r1 = 1.0 - 2.0*vgm*rpj*(1.0 - rj)/evj;
                c1 = 2.0*rpj/evj;
              }
              const Real tuj = pth ? upt(cc, r1, tj, 0.5*kro*dz) : tj;
              const Real r2 = tuj*tj*r1;
              const Real sold[3] = {jS[cc][1], jS[cc][2], 0.0};
              for (int s=0; s<3; ++s) {
                const Real dsf = (s < 2) ? jfc[cc][s] : 0.0;
                const Real dsu = (s < 2) ? juc[cc][s] : 0.0;
                const Real dpl = ci*dsf + co*dsu;
                const Real dql = ci*dsu + co*dsf;
                const Real dpu = ci*dslv[s] + co*dsfw[s];
                const Real dqu = ci*dsfw[s] + co*dslv[s];
                const Real s2 = tuj*(r1*dql + c1*sold[s]) + dpl;
                jS[cc][s] = tuj*(r2*dqu + s2) + dpu;
              }
              jfc[cc][0] = dsfw[1];
              jfc[cc][1] = dsfw[2];
              juc[cc][0] = dsuu[1];
              juc[cc][1] = dsuu[2];
            }
            RtF rr = Rc[cc];
            RtF ss = Hc[cc];
            if (FRM == 1) {
              // tm: reflection addition at the face, T^2 through each half
              const RtF bb = static_cast<RtF>(bt);
              if (dil) {
                // ck_sph_dilute: d passes unchanged, the up ray is diluted by
                // A_b/A_a and gains the down rays that miss the sphere below
                const RtF rh = (one - bb)/(one + bb);
                rr = one - rh*(one - rr);
                ss = rh*ss;
              } else if (fbl) {
                Real alx, bex, c1x, r1x;
                fblm(cc, static_cast<Real>(rr), bt, dtfc, alx, bex, c1x, r1x);
                ss = static_cast<RtF>(c1x)*ss;
                rr = static_cast<RtF>(r1x);
              } else if (vef) {
                // ck_sph_face = 5: Phi f J and A s (u - d) continuous
                const RtF rpv = static_cast<RtF>((1.0 - bt)/(1.0 + bt)*vsb/vsa);
                const RtF gmv = static_cast<RtF>(vgm);
                const RtF ev = (one + rr) + gmv*rpv*(one - rr);
                ss = two*rpv*ss/ev;
                rr = one - two*gmv*rpv*(one - rr)/ev;
              } else {
                const RtF dn = one + rr*bb;
                const RtF rn = (rr + bb)/dn;
                ss = (one - bb)*ss/dn;
                rr = rn;
              }
              // ck_sph_face = 3, 4: the up crossing absorbs along the cone path
              // (tu), the down crossing and both emissions are unchanged
              const RtF tu = pth ? static_cast<RtF>(upt(cc, static_cast<Real>(rr),
                                                        static_cast<Real>(tr),
                                                        0.5*kro*dz))
                                 : tr;
              ss = tu*(rr*ql + ss) + pl;
              rr = tu*tr*rr;
              ss = tu*(rr*qu + ss) + pu;
              rr = tu*tr*rr;
            } else {
              // sd: the face is a rescaling by A_below/A_above, the layer is the
              // Moebius map with the shared denominator
              const RtF rh = static_cast<RtF>(rat);
              rr *= rh;
              ss *= rh;
              const RtF t2 = tr*tr;
              RtF dn = (one + rr) + t2*(one - rr);
              RtF rn = ((one - t2) + rr*(one + t2))/dn;
              ss = (tr*((one - rr)*ql + two*ss) + pl*(one + rr))/dn;
              rr = rn;
              dn = (one + rr) + t2*(one - rr);
              rn = ((one - t2) + rr*(one + t2))/dn;
              ss = (tr*((one - rr)*qu + two*ss) + pu*(one + rr))/dn;
              rr = rn;
            }
            Rc[cc] = rr;
            Hc[cc] = ss;
            kcu[cc] = kru;
            sfc[cc] = sfv;
            suc[cc] = suu;
            dtf[cc] = dtn;
          }
        }
        for (int cc=0; cc<NC; ++cc) {
          I_down[cc][ie+1] = Rc[cc];
          Cmx[cc][ie+1] = Hc[cc];
          if constexpr (JAC) {
            Js0[cc][ie+1] = jS[cc][0];
            Js1[cc][ie+1] = jS[cc][1];
            Js2[cc][ie+1] = jS[cc][2];
          }
        }
        // ---- the face solve: the stored relation meets the down ray ---------
        // Returns the pair (u, d) on the BELOW side of face f, given the down
        // intensity ABOVE it.  tm inverts u = R d_b + Sc against
        // d_b = d_a + beta (d_a - u); sd combines D = -G S + H, which the face
        // merely rescales, with S - D = d_a.
        auto fsolve = [&](const int cc, const int f, const Real bt,
                          const Real rat, const RtF da, RtF &ub, RtF &db,
                          const Real dtc) {
          if (FRM == 1) {
            const RtF rr = I_down[cc][f];
            const RtF ss = Cmx[cc][f];
            const RtF bb = static_cast<RtF>(bt);
            if (dil) {
              ub = rr*da + ss;
              db = da;
            } else if (fbl) {
              Real alx, bex, c1x, r1x;
              fblm(cc, static_cast<Real>(rr), bt, dtc, alx, bex, c1x, r1x);
              db = static_cast<RtF>(alx)*da - static_cast<RtF>(bex)*ss;
              ub = rr*db + ss;
            } else {
              ub = (rr*(one + bb)*da + ss)/(one + rr*bb);
              db = da + bb*(da - ub);
            }
          } else {
            const RtF gg = I_down[cc][f];
            const RtF hh = Cmx[cc][f];
            const RtF rh = static_cast<RtF>(rat);
            const RtF sv = (da + rh*hh)/(one + rh*gg);
            const RtF dd = hh - gg*sv;
            ub = sv + dd;
            db = sv - dd;
          }
        };
        // ck_sph_face = 5: the face solve of the VEF jump conditions (see
        // ck_sph_face): gamma (u_a + d_a) = u_b + d_b, u_a - d_a = rho' (u_b - d_b)
        auto vsolve = [&](const RtF rr, const RtF ss, const Real rp, const Real gm,
                          const RtF da, RtF &ub, RtF &db) {
          const RtF rpv = static_cast<RtF>(rp);
          const RtF gmv = static_cast<RtF>(gm);
          const RtF ev = (one + rr) + gmv*rpv*(one - rr);
          db = (two*gmv*da - (one - gmv*rpv)*ss)/ev;
          ub = rr*db + ss;
        };
        // ---- PASS 2: downward ----------------------------------------------
        RtF dcu[NC], ubf[NC];
        Real slc[NC], sfu[NC];
        // JAC: W's successor Q = alpha W + beta at the face above, dD_{i+1}
        // (slots B_{i-1}, B_i, B_{i+1}), dd below face i+1 (same slots), and the
        // two carried endpoint sources sfu, slc (slots B_i, B_{i+1})
        Real jQ[NJ], jDu[NJ][3], jdm[NJ][3], jfu[NJ][2], jlc[NJ][2];
        {
          const Real bt = BTF(m,k,j,ie+1,icut);
          // at the TOP face A_above is the face area itself, so the frame change
          // and the per-unit-area rescaling coincide
          const Real fsc = ACC(m,k,j,ie)/AFC(m,k,j,ie+1);
          for (int cc=0; cc<NC; ++cc) {
            RtF ub, db;
            // ck_sph_face = 5: s and gamma of the top face, and the reflection
            // d_a = zeta u_a + datum of Auer's H = h J
            Real vsb = 1.0, vrp = 1.0, vgm = 1.0, vzt = 0.0, vkp = 0.0;
            if (vef) {
              const int cv = bandc[cc]*CK_NG + gc[cc];
              vsb = vfs_g(m,cv,ie+1,k,j);
              vrp = (1.0 - bt)/(1.0 + bt)*vsb/vfs_g(m,cv,ie+2,k,j);
              vgm = vfg_g(m,cv,ie+1,k,j);
              vzt = vfg_g(m,cv,ie+2,k,j);
              const RtF rr = I_down[cc][ie+1];
              const RtF rpv = static_cast<RtF>(vrp);
              const RtF gmv = static_cast<RtF>(vgm);
              const RtF ev = (one + rr) + gmv*rpv*(one - rr);
              const RtF rup = one - two*gmv*rpv*(one - rr)/ev;
              const RtF sup = two*rpv*Cmx[cc][ie+1]/ev;
              const RtF zt = static_cast<RtF>(vzt);
              const RtF dat = (zt*sup + Itp[cc])/(one - zt*rup);
              vkp = static_cast<Real>(zt*(two*rpv/ev)/(one - zt*rup));   // dd_a/dSc
              vsolve(rr, Cmx[cc][ie+1], vrp, vgm, dat, ub, db);
            } else {
              fsolve(cc, ie+1, bt, fsc, Itp[cc], ub, db, dtf[cc]);
            }
            Fb_g(m,blk,ie+1,k,j) += (vef ? wfc[cc]*vsb/(SQRT3_VEF*muc[cc]) : wfc[cc])
                                  *static_cast<Real>(ub - db)*fsc;
            dcu[cc] = db;
            ubf[cc] = ub;
            const Real btop = Bb_g(m,bandc[cc],ie,k,j);
            slc[cc] = btop;
            sfu[cc] = btop;
            if constexpr (JAC) {
              // the top datum is fixed: d^+ carries nothing, W = 0
              const Real rj = static_cast<Real>(I_down[cc][ie+1]);
              const Real dnj = 1.0 + rj*bt;
              Real al = dil ? 1.0 : (1.0 + bt)/dnj;
              Real be = dil ? 0.0 : bt/dnj;
              if (fbl) {
                Real c1x, r1x;
                fblm(cc, rj, bt, dtf[cc], al, be, c1x, r1x);
              }
              // ck_sph_top = 2: the datum is (1 - e^-dtau) B_ie, so d^+ at the top
              // face depends on B_ie (slot 1)
              Real tfd = 0.0;
              if (tpm == 2) {
                const Real bie = Bb_g(m,bandc[cc],ie,k,j);
                tfd = (bie > 0.0) ? static_cast<Real>(Itp[cc])/bie : 0.0;
              }
              const Real s0 = Js0[cc][ie+1];
              const Real s1 = Js1[cc][ie+1];
              const Real s2 = Js2[cc][ie+1];
              Real aD = al;
              if (vef) {
                // the reflected datum depends on Sc_top through u_a:
                // dd_b/dSc = al kp - be, dD/dSc = aD (1 - gamma (1 - R) kp)
                const Real ev = (1.0 + rj) + vgm*vrp*(1.0 - rj);
                al = 2.0*vgm/ev;
                be = (1.0 - vgm*vrp)/ev - al*vkp;
                aD = 2.0/ev*(1.0 - vgm*(1.0 - rj)*vkp);
              }
              jDu[cc][0] = aD*s0;
              jDu[cc][1] = aD*s1;
              jDu[cc][2] = aD*s2;
              jdm[cc][0] = -be*s0;
              jdm[cc][1] = -be*s1;
              jdm[cc][2] = -be*s2;
              if (tpm == 2) {
                jDu[cc][1] -= al*(1.0 - rj)*tfd;
                jdm[cc][1] += al*tfd;
              }
              jQ[cc] = be;
              jfu[cc][0] = 1.0;
              jfu[cc][1] = 0.0;
              jlc[cc][0] = 1.0;
              jlc[cc][1] = 0.0;
            }
          }
        }
        Real xTc = xT_g(m,k,j,ie);
        Real xPc = xP_g(m,k,j,ie);
        for (int i=ie; i>icut-1; --i) {
          const Real dz = dx1(m,k,j,i);
          const Real rho = rhoN(m,k,j,i);
          const Real drho = rho*dz;
          const int iTc = static_cast<int>(xTc);
          const int iPc = static_cast<int>(xPc);
          const Real fTc = xTc - static_cast<Real>(iTc);
          const Real fPc = xPc - static_cast<Real>(iPc);
          const int il = (i > icut) ? i-1 : i;
          const Real dzl = dx1(m,k,j,il);
          const Real rhol = rhoN(m,k,j,il);
          const Real xTl = xT_g(m,k,j,il);
          const Real xPl = xP_g(m,k,j,il);
          const int iTl = static_cast<int>(xTl);
          const int iPl = static_cast<int>(xPl);
          const Real fTl = xTl - static_cast<Real>(iTl);
          const Real fPl = xPl - static_cast<Real>(iPl);
          const Real bt = BTF(m,k,j,i,icut);
          const Real rat = (i == icut)
              ? AFC(m,k,j,icut)/ACC(m,k,j,icut)
              : ACC(m,k,j,i-1)/ACC(m,k,j,i);
          const Real fsc = (i == icut)
              ? 1.0 : ACC(m,k,j,i-1)/AFC(m,k,j,i);
          for (int cc=0; cc<NC; ++cc) {
            const int b = bandc[cc];
            const Real kro = kcu[cc];
            const Real bown = Bb_g(m,b,i,k,j);
            const Real dth = 0.5*kro*dz;
            // the layer joining cells i-1 and i, from below this time
            Real suv = bown, sfv = bown, snl = bown, krl = kro;
            Real dtfc = dth;   // ck_sph_face = 2: face i's depth, as pass 1 formed it
            if (i > icut) {
              // ck-store: krof without its store -- pass 1 stored this very
              // product (the same cached kappa times the same rho) at cell il
              krl = ckfus ? ckkro_g(m,blk*NC+cc,il,k,j)
                          : kapof(cc, il, b, iTl, fTl, iPl, fPl, false)*rhol;
              const Real bfar = Bb_g(m,b,il,k,j);
              suv = BFace(krl, kro, bfar, bown, bface_on);
              snl = BFace(kro, krl, bown, bfar, bface_on);
              const Real dt_l = 0.5*krl*dzl;
              const Real dtc = dt_l + dth;
              sfv = (dtc > 0.0) ? (snl + (suv - snl)*(dt_l/dtc))
                                : (0.5*(snl + suv));
              dtfc = dtc;
            }
            // ck_sph_face = 5: s of the frame below face i and of cell i, gamma of
            // face i; the cell's half steps along dtau/s, its weight wfc s
            Real vsb = 1.0, vsa = 1.0, vgm = 1.0;
            if (vef) {
              const int cv = b*CK_NG + gc[cc];
              vsb = vfs_g(m,cv,i,k,j);
              vsa = vfs_g(m,cv,i+1,k,j);
              vgm = vfg_g(m,cv,i,k,j);
            }
            const Real dthe = vef ? dth*(SQRT3_VEF*muc[cc])/vsa : dth;
            const Real wfi = vef ? wfc[cc]*vsa/(SQRT3_VEF*muc[cc]) : wfc[cc];
            const Real ratv = vef ? rat*vsb/vsa : rat;
            Real dsrc;
            // the down ray: the upper half of cell i, then the lower half
            step(cc, i, dthe, sfu[cc], slc[cc], dcu[cc], dsrc, true, false);
            Src_g(m,blk,i,k,j) += wfi/dz*dsrc;
            step(cc, i, dthe, suv, sfv, dcu[cc], dsrc, false, false);
            Src_g(m,blk,i,k,j) += wfi/dz*dsrc;
            // face i: the up ray is recovered here, and A D is made continuous
            // by construction -- the flux is A_below (u_b - d_b)/A_face and the
            // up ray is handed d_above + (u_b - d_b) A_below/A_above
            RtF ub, db;
            if (vef) {
              vsolve(I_down[cc][i], Cmx[cc][i], ratv, vgm, dcu[cc], ub, db);
            } else {
              fsolve(cc, i, bt, rat, dcu[cc], ub, db, dtfc);
            }
            const RtF dm = ub - db;
            Fb_g(m,blk,i,k,j) += (vef ? wfc[cc]*vsb/(SQRT3_VEF*muc[cc]) : wfc[cc])
                                 *static_cast<Real>(dm)*fsc;
            RtF ua = dcu[cc] + dm*static_cast<RtF>(ratv);
            dcu[cc] = db;
            // the up ray through cell i: the lower half, then the upper half
            if (pth) {
              // ck_sph_face = 3, 4: pass 1's up crossing, cone-path absorption
              RtF e0x, cix, cox;
              cofs(cc, i, dth, e0x, cix, cox, false, false);
              const RtF rfx = rfm(cc, I_down[cc][i], bt, dtfc);
              const RtF tux = static_cast<RtF>(upt(cc, static_cast<Real>(rfx),
                                                   static_cast<Real>(one - e0x), dth));
              RtF emx = cix*static_cast<RtF>(sfv) + cox*static_cast<RtF>(suv);
              dsrc = static_cast<Real>((one - tux)*ua - emx);
              ua = tux*ua + emx;
              Src_g(m,blk,i,k,j) += wfi/dz*dsrc;
              emx = cix*static_cast<RtF>(slc[cc]) + cox*static_cast<RtF>(sfu[cc]);
              dsrc = static_cast<Real>((one - tux)*ua - emx);
              ua = tux*ua + emx;
              Src_g(m,blk,i,k,j) += wfi/dz*dsrc;
            } else {
              step(cc, i, dthe, sfv, suv, ua, dsrc, false, false);
              Src_g(m,blk,i,k,j) += wfi/dz*dsrc;
              step(cc, i, dthe, slc[cc], sfu[cc], ua, dsrc, false, false);
              Src_g(m,blk,i,k,j) += wfi/dz*dsrc;
            }
            // CLOSE THE CELL ON THE FLUX IT REPORTED.  The up ray leaving cell i
            // is the u_below the face above it already used; the two agree to
            // round-off, and adding the difference here is what makes the
            // per-cell budget telescope with nothing left over.
            Src_g(m,blk,i,k,j) += wfi/dz
                *static_cast<Real>(ua - ubf[cc]);
            ubf[cc] = ub;
            Em_g(m,blk,i,k,j) += 2.0*(wfc[cc]/muc[cc])*kro*bown;
            // the direct beam, UNCHANGED (see the note in the four-pass sweep):
            // it crosses whole cells, carries no source and is not touched by
            // the face mixing, so it rides this pass exactly as it rode the
            // four-pass down-sweep.
            if (!ckfus) {
              if (BSP) Krs[cc][i] = static_cast<RtF>(kro);
              if (!(bpair && (cc & 1))) {
                const Real kap = kapof(cc, i, b, iTc, fTc, iPc, fPc, false);
                tausw[cc] += kap*drho;
                if (lit && !BSP) {
                  const Real tnew = RT_EXP(-static_cast<RtF>(tausw[cc]*facsw));
                  Qb_g(m,blk,i,k,j) += (1.0-albedo)*Fstar*ckswf(b)*wgc[cc]
                            /facsw*(transw[cc] - tnew)/dz;
                  transw[cc] = tnew;
                }
              }
            }
            if constexpr (JAC) {
              // row i of the tridiagonal; slots (B_{i-1}, B_i, B_{i+1})
              RtF e0r, cir, cor;
              cofs(cc, i, dthe, e0r, cir, cor, false, false);
              const Real ci = static_cast<Real>(cir);
              const Real co = static_cast<Real>(cor);
              const Real tj = 1.0 - static_cast<Real>(e0r);
              Real dsuv[3] = {0.0, 1.0, 0.0};
              Real dsnl[3] = {0.0, 1.0, 0.0};
              Real dsfw[3] = {0.0, 1.0, 0.0};
              if (i > icut) {
                const Real wa = BFaceW(krl, kro, bface_on);   // dsuv/dB_i
                const Real wc = BFaceW(kro, krl, bface_on);   // dsnl/dB_{i-1}
                const Real dtl = 0.5*krl*dzl;
                const Real dtcj = dtl + dth;
                const Real ffj = (dtcj > 0.0) ? (dtl/dtcj) : 0.5;
                const Real s0 = (1.0 - ffj)*wc + ffj*(1.0 - wa);
                dsuv[0] = 1.0 - wa;
                dsuv[1] = wa;
                dsnl[0] = wc;
                dsnl[1] = 1.0 - wc;
                dsfw[0] = s0;
                dsfw[1] = 1.0 - s0;
              }
              // dd_i^+: the down ray's two half steps through cell i
              Real v[3];
              for (int s=0; s<3; ++s) {
                const Real dfu = (s > 0) ? jfu[cc][s-1] : 0.0;
                const Real dlc = (s > 0) ? jlc[cc][s-1] : 0.0;
                v[s] = tj*jdm[cc][s] + ci*dfu + co*dlc;
                v[s] = tj*v[s] + ci*dsuv[s] + co*dsfw[s];
              }
              const Real rj = static_cast<Real>(I_down[cc][i]);
              const Real dnj = 1.0 + rj*bt;
              Real al = dil ? 1.0 : (1.0 + bt)/dnj;
              Real be = dil ? 0.0 : bt/dnj;
              Real wr = dil ? tj*tj*tj*tj*(1.0 - bt)/(1.0 + bt)*jQ[cc]
                            : tj*tj*tj*tj*(1.0 - bt)/dnj*jQ[cc];
              if (fbl) {
                Real c1x, r1x;
                fblm(cc, rj, bt, dtfc, al, be, c1x, r1x);
                wr = tj*tj*tj*tj*c1x*jQ[cc];
                if (pth) {
                  const Real tuj = upt(cc, r1x, tj, dth);
                  wr = tuj*tuj*tj*tj*c1x*jQ[cc];
                }
              } else if (pth) {
                const Real rhd = (1.0 - bt)/(1.0 + bt);
                const Real tuj = upt(cc, 1.0 - rhd*(1.0 - rj), tj, dth);
                wr = tuj*tuj*tj*tj*rhd*jQ[cc];
              }
              const Real s0 = Js0[cc][i];
              const Real s1 = Js1[cc][i];
              const Real s2 = Js2[cc][i];
              Real omr = 1.0 - rj;
              // D = aD (Sc - omr d^+): aD = al and omr = 1 - R except under
              // ck_sph_face = 5, where D = (2/E) (Sc - gamma (1 - R) d^+) and d
              // crosses the face with al = 2 gamma/E, be = (1 - gamma rho')/E
              Real aD = al;
              if (vef) {
                const Real ev = (1.0 + rj) + vgm*ratv*(1.0 - rj);
                al = 2.0*vgm/ev;
                be = (1.0 - vgm*ratv)/ev;
                aD = 2.0/ev;
                wr = tj*tj*tj*tj*(2.0*ratv/ev)*jQ[cc];
                omr = vgm*(1.0 - rj);
              }
              // dD_i/dB_{i-2 .. i+1}
              const Real dd0 = aD*(1.0 + omr*wr)*s0;
              const Real dd1 = aD*(s1 - omr*v[0]);
              const Real dd2 = aD*(s2 - omr*v[1]);
              const Real dd3 = -aD*omr*v[2];
              const Real wj = wfi/dz;
              const Real jm = wj*(ratv*dd1 - jDu[cc][0]);
              const Real j0 = wj*(ratv*dd2 - jDu[cc][1]);
              const Real jp = wj*(ratv*dd3 - jDu[cc][2]);
              CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,i), j0*ckdb_g(m,b,i,k,j));
              if ((i > icut || (ckdif_ && icut > icut_g(m,k,j)))
                  && jm > 0.0) {
                CkJacQAdd(&ckjac_g(m,3*blk+0,k,j,i),
                                   jm*ckdb_g(m,b,i-1,k,j));
              }
              if (i < ie && jp > 0.0) {
                CkJacQAdd(&ckjac_g(m,3*blk+2,k,j,i),
                                   jp*ckdb_g(m,b,i+1,k,j));
              }
              // shift everything down one cell: slots become (B_{i-2}, B_{i-1},
              // B_i), which is (B_{i'-1}, B_i', B_{i'+1}) of the next cell
              jDu[cc][0] = dd0;
              jDu[cc][1] = dd1;
              jDu[cc][2] = dd2;
              const Real qn = al*wr + be;
              jdm[cc][0] = -qn*s0;
              jdm[cc][1] = al*v[0] - be*s1;
              jdm[cc][2] = al*v[1] - be*s2;
              jQ[cc] = qn;
              jfu[cc][0] = dsfw[0];
              jfu[cc][1] = dsfw[1];
              jlc[cc][0] = dsnl[0];
              jlc[cc][1] = dsnl[1];
            }
            kcu[cc] = krl;
            slc[cc] = snl;
            sfu[cc] = sfv;
          }
          xTc = xTl;
          xPc = xPl;
        }
        // ---- problem/ck_beam_sph: the pseudo-spherical beam.  Identical to the
        // block in the four-pass sweep below, which carries the derivation and
        // the sketch; it is repeated rather than shared because the two forms
        // reach it from different places.
        if (BSP && lit_sph && !ckfus && !(ckfst && cbt_)) {
          const Real rcut = X1F(m,icut);
          RtF thi[NC];
          Real tauh[NC];
          for (int cc=0; cc<NC; ++cc) {
            tauh[cc] = 0.0;
            thi[cc] = static_cast<RtF>(1.0);
          }
          // mu0 < 0: the ray reaches the TOP face itself through its tangent
          // point, so tau there is not 0 (see BeamTopTau)
          if (mu0 < 0.0) BeamTopTau(rcut, sinz, icut, ie, m, X1F, Krs, tauh, thi);
          for (int i=ie; i>icut-1; --i) {
            const Real rf = X1F(m,i);
            const Real bb = rf*sinz;
            const Real b2 = bb*bb;
            Real taul[NC];
            for (int cc=0; cc<NC; cc+=bst) taul[cc] = 0.0;
            bool dark = false;
            int jlo = i;
            if (mu0 < 0.0) {
              if (bb <= rcut) {
                dark = true;
              } else {
                jlo = icut;
                for (int jj=i-1; jj>=icut; --jj) {
                  if (X1F(m,jj) <= bb) {
                    jlo = jj;
                    break;
                  }
                }
              }
            }
            if (!dark) {
              const Real rl = X1F(m,jlo);
              Real prev = (rl*rl > b2) ? sqrt(rl*rl - b2) : 0.0;
              for (int jj=jlo; jj<ie+1; ++jj) {
                const Real ru = X1F(m,jj+1);
                const Real cur = sqrt(ru*ru - b2);
                const Real ds = (jj < i) ? 2.0*(cur - prev) : (cur - prev);
                prev = cur;
                for (int cc=0; cc<NC; cc+=bst) {
                  taul[cc] += ds*static_cast<Real>(Krs[cc][jj]);
                }
              }
            }
            for (int cc=0; cc<NC; cc+=bst) {
              const Real dtl = dark ? 1.0e30 : (taul[cc] - tauh[cc]);
              const RtF tlo = dark ? static_cast<RtF>(0.0)
                                   : RT_EXP(-static_cast<RtF>(taul[cc]));
              const Real dif = static_cast<Real>(thi[cc])
                             - static_cast<Real>(tlo);
              const Real fac = (fabs(dtl) > 1.0e-3)
                  ? (dif/dtl)
                  : (static_cast<Real>(thi[cc])*(1.0 - 0.5*dtl));
              Qb_g(m,blk,i,k,j) += (1.0-albedo)*Fstar*ckswf(bandc[cc])*wgc[cc]
                                 * fac*static_cast<Real>(Krs[cc][i]);
              tauh[cc] = dark ? 1.0e30 : taul[cc];
              thi[cc] = tlo;
            }
            Real tmin = tauh[0];
            for (int cc=bst; cc<NC; cc+=bst) {
              tmin = (tauh[cc] < tmin) ? tauh[cc] : tmin;
            }
            if (tmin > 60.0) break;
          }
        }
        return;
      }
      // ---- problem/ck_spherical: THE TWO PROBE PASSES.  See the flag's note.
      //
      // The real down-sweep's face mixing c = beta (d_above - u_below) needs the
      // UPWARD intensity at every face, which it does not have.  CONSERVATION,
      // though, does not care whether c is accurate -- it cares only that the
      // two rays use the SAME c, which is what Cmx enforces.  So the probe's
      // only job is to make c ACCURATE, and for that it must reproduce the one
      // thing that dominates u on a deep grid: the r^-2 dilution.  A bare
      // plane-parallel up-sweep does not (it is the A_top/A_cut = 3.8 error this
      // switch exists to remove), and feeding its u into c over-dilutes the real
      // up ray by the same factor -- MEASURED, A F then varies by 26 % over the
      // top twelve faces.  So the probe is two passes:
      //   P1  a plane-parallel DOWN pass, d0, parked in I_down;
      //   P2  an UP pass carrying the face mixing against d0, whose u_below at
      //       each face is parked in Cmx.
      // P2 therefore dilutes correctly and u_probe - u_true is O(beta) in d
      // alone, i.e. c is right to O(beta^2), which is the grey mode-0 bound.
      // Neither pass accumulates anything.  Itop holds the top face's incoming
      // datum across the borrowing of I_down.
      RtF Itop[NC];
      for (int cc=0; cc<NC; ++cc) {
        Itop[cc] = I_down[cc][ie+1];
      }
      const Real bt_top = SPH ? BTF(m,k,j,ie+1,icut) : 0.0;
      if (SPH) {
        // ---- P1: the plane-parallel down probe, d0 -> I_down
        RtF Idp[NC];
        Real kdp[NC];
        for (int cc=0; cc<NC; ++cc) {
          Idp[cc] = Itop[cc];
          I_down[cc][ie+1] = Itop[cc];
        }
        for (int i=ie; i>icut-1; --i) {
          const Real rho = rhoN(m,k,j,i);
          const Real dz = dx1(m,k,j,i);
          const Real dzf = (i < ie) ? dx1(m,k,j,i+1) : dz;
          const Real xTv = xT_g(m,k,j,i);
          const Real xPv = xP_g(m,k,j,i);
          const int iT = static_cast<int>(xTv);
          const int iP = static_cast<int>(xPv);
          const Real fT = xTv - static_cast<Real>(iT);
          const Real fP = xPv - static_cast<Real>(iP);
          for (int cc=0; cc<NC; ++cc) {
            const int b = bandc[cc];
            const Real kro = krof(cc, i, b, iT, fT, iP, fP, rho, true);
            const Real bown = Bb_g(m,b,i,k,j);
            Real dsrc;
            if (i == ie) {
              step(cc, i, 0.5*kro*dz, bown, bown, Idp[cc], dsrc, true, true);
            } else {
              const Real kru = kdp[cc];
              const Real bfar = Bb_g(m,b,i+1,k,j);
              const Real dt_u = 0.5*kru*dzf;
              const Real dt_l = 0.5*kro*dz;
              const Real s_l = BFace(kru, kro, bfar, bown, bface_on);
              const Real s_u = BFace(kro, kru, bown, bfar, bface_on);
              const Real dtc = dt_l + dt_u;
              const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dt_l/dtc))
                                           : (0.5*(s_l + s_u));
              step(cc, i+1, dt_u, s_u, s_f, Idp[cc], dsrc, false, false);
              I_down[cc][i+1] = Idp[cc];
              step(cc, i, dt_l, s_f, s_l, Idp[cc], dsrc, true, true);
            }
            kdp[cc] = kro;
          }
        }
        {
          const Real dz = dx1(m,k,j,icut);
          for (int cc=0; cc<NC; ++cc) {
            const Real bcut = Bb_g(m,bandc[cc],icut,k,j);
            Real dsrc;
            step(cc, icut, 0.5*kdp[cc]*dz, bcut, bcut, Idp[cc],
                 dsrc, false, false);
            I_down[cc][icut] = Idp[cc];
          }
        }
        // ---- P2: the up probe, WITH the face mixing against d0.  u_below at
        // each face goes to Cmx, which the real down-sweep reads and then
        // overwrites with the c it forms from it.
        RtF Iup_p[NC];
        Real kpr[NC];
        {
          const Real dz = dx1(m,k,j,icut);
          const Real rho = rhoN(m,k,j,icut);
          const Real xTv = xT_g(m,k,j,icut);
          const Real xPv = xP_g(m,k,j,icut);
          const int iT = static_cast<int>(xTv);
          const int iP = static_cast<int>(xPv);
          const Real fT = xTv - static_cast<Real>(iT);
          const Real fP = xPv - static_cast<Real>(iP);
          const Real bt_c = BTF(m,k,j,icut,icut);
          for (int cc=0; cc<NC; ++cc) {
            const int b = bandc[cc];
            const Real Iint_b = (int_at_cut ? boltz_sigma/M_PI*Tint4
                            * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
            const Real bcut = Bb_g(m,b,icut,k,j);
            Iup_p[cc] = static_cast<RtF>(bcut + Iint_b);
            Iup_p[cc] += static_cast<RtF>(MIXF(bt_c,
                             static_cast<Real>(Iup_p[cc]),
                             static_cast<Real>(I_down[cc][icut])));
            kpr[cc] = krof(cc, icut, b, iT, fT, iP, fP, rho, false);
            Real dsrc;
            step(cc, icut, 0.5*kpr[cc]*dz, bcut, bcut, Iup_p[cc],
                 dsrc, true, false);
          }
        }
        for (int i=icut; i<ie; ++i) {
          const Real rhou = rhoN(m,k,j,i+1);
          const Real dzu = dx1(m,k,j,i+1);
          const Real dzl = dx1(m,k,j,i);
          const Real xTv = xT_g(m,k,j,i+1);
          const Real xPv = xP_g(m,k,j,i+1);
          const int iT = static_cast<int>(xTv);
          const int iP = static_cast<int>(xPv);
          const Real fT = xTv - static_cast<Real>(iT);
          const Real fP = xPv - static_cast<Real>(iP);
          const Real bt_f = BTF(m,k,j,i+1,icut);
          for (int cc=0; cc<NC; ++cc) {
            const int b = bandc[cc];
            const Real kru = krof(cc, i+1, b, iT, fT, iP, fP, rhou,
                                  false);
            const Real krl = kpr[cc];
            const Real bl = Bb_g(m,b,i,k,j);
            const Real bu = Bb_g(m,b,i+1,k,j);
            const Real dt_l = 0.5*krl*dzl;
            const Real dt_u = 0.5*kru*dzu;
            const Real s_l = BFace(kru, krl, bu, bl, bface_on);
            const Real s_u = BFace(krl, kru, bl, bu, bface_on);
            const Real dtc = dt_l + dt_u;
            const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dt_l/dtc))
                                         : (0.5*(s_l + s_u));
            Real dsrc;
            step(cc, i, dt_l, s_l, s_f, Iup_p[cc], dsrc, false, false);
            Cmx[cc][i+1] = Iup_p[cc];        // u_below, for the real down-sweep
            Iup_p[cc] += static_cast<RtF>(MIXF(bt_f,
                             static_cast<Real>(Iup_p[cc]),
                             static_cast<Real>(I_down[cc][i+1])));
            step(cc, i+1, dt_u, s_f, s_u, Iup_p[cc], dsrc, true, false);
            kpr[cc] = kru;
          }
        }
        {
          const Real dz = dx1(m,k,j,ie);
          for (int cc=0; cc<NC; ++cc) {
            const Real btop = Bb_g(m,bandc[cc],ie,k,j);
            Real dsrc;
            step(cc, ie, 0.5*kpr[cc]*dz, btop, btop, Iup_p[cc],
                 dsrc, false, false);
            Cmx[cc][ie+1] = Iup_p[cc];       // u_below at the top face
          }
        }
      }
      for (int cc=0; cc<NC; ++cc) {
        Idn[cc] = Itop[cc];
        if (SPH) {
          // c from the probe's u_below at the top face and the boundary datum
          const RtF cmx = static_cast<RtF>(MIXF(bt_top,
                              static_cast<Real>(Cmx[cc][ie+1]),
                              static_cast<Real>(Itop[cc])));
          Cmx[cc][ie+1] = cmx;              // the up-sweep reuses this exact c
          Idn[cc] += cmx;
          I_down[cc][ie+1] = Itop[cc];      // give the slot back
        }
      }
      // down-sweep, recording the FACE intensity in between the two halves
      for (int i=ie; i>icut-1; --i) {
        const Real rho = rhoN(m,k,j,i);
        const Real dz = dx1(m,k,j,i);
        const Real drho = rho*dz;
        const Real dzf = (i < ie) ? dx1(m,k,j,i+1) : dz;
        // the face between cells i+1 and i
        const Real bt_d = SPH ? BTF(m,k,j,i+1,icut) : 0.0;
        const Real xTv = xT_g(m,k,j,i);
        const Real xPv = xP_g(m,k,j,i);
        const int iT = static_cast<int>(xTv);
        const int iP = static_cast<int>(xPv);
        const Real fT = xTv - static_cast<Real>(iT);
        const Real fP = xPv - static_cast<Real>(iP);
        for (int cc=0; cc<NC; ++cc) {
          const int b = bandc[cc];
          // ck_impl_frozen_op: kap itself is needed only by the direct beam,
          // which a frozen pass does not re-run.
          // ck_sweep_cache: under the spherical form P1 has already been here,
          // so this reads Kpc; the plane-parallel path makes this the filling
          // pass.  Same expression, same value, either way.
          const Real kap = ckfus ? 0.0
              : kapof(cc, i, b, iT, fT, iP, fP, !SPH);
          Real kro;
          if (ckfus) {
            kro = ckkro_g(m,blk*NC+cc,i,k,j);
          } else {
            kro = kap*rho;
            if (ckfst) ckkro_g(m,blk*NC+cc,i,k,j) = kro;
          }
          const Real bown = Bb_g(m,b,i,k,j);
          Real dsrc;
          if (i == ie) {
            // the top half layer: the upper half of cell ie, entered at the top
            // face, its source held at the cell's own value -- there is nothing
            // above it to interpolate towards
            step(cc, i, 0.5*kro*dz, bown, bown, Idn[cc], dsrc, true, !SPH);
            Src_g(m,blk,i,k,j) += wfc[cc]/dz*dsrc;
            if (jck) {
              // both endpoints are B_ie and the entering datum is a boundary
              // condition, so the whole row entry is the layer's own emission
              Real e0j, cij, coj;
              jcof(cc, i, 0.5*kro*dz, muc[cc], e0j, cij, coj);
              CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,i),
                  -(wfc[cc]/dz)*(cij + coj)*ckdb_g(m,b,i,k,j));
              jdn[cc] = cij + coj;
            }
          } else {
            // the layer between the centres of cells i+1 and i.  BFace is
            // applied to both endpoints and symmetrically, so a radiatively
            // inert neighbour leaves a layer emitting with its own Planck
            // function; with equal opacities this is the straight line.
            const Real kru = kfar[cc];
            const Real bfar = Bb_g(m,b,i+1,k,j);
            const Real dt_u = 0.5*kru*dzf;
            const Real dt_l = 0.5*kro*dz;
            const Real s_l = BFace(kru, kro, bfar, bown, bface_on);
            const Real s_u = BFace(kro, kru, bown, bfar, bface_on);
            const Real dtc = dt_l + dt_u;
            const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dt_l/dtc))
                                         : (0.5*(s_l + s_u));
            // the probe's u_below at this face, read before the real
            // down-sweep reclaims the slot for c (see Cmx)
            const Real ubf = SPH ? static_cast<Real>(Cmx[cc][i+1]) : 0.0;
            step(cc, i+1, dt_u, s_u, s_f, Idn[cc], dsrc, false, false);
            Src_g(m,blk,i+1,k,j) += wfc[cc]/dzf*dsrc;
            I_down[cc][i+1] = Idn[cc];
            if (SPH) {
              // d_below = d_above + beta (d_above - u_below): into cell i's frame
              const RtF cmx = static_cast<RtF>(MIXF(bt_d, ubf,
                                  static_cast<Real>(Idn[cc])));
              Cmx[cc][i+1] = cmx;           // the up-sweep reuses this exact c
              Idn[cc] += cmx;
            }
            step(cc, i, dt_l, s_f, s_l, Idn[cc], dsrc, true, !SPH);
            Src_g(m,blk,i,k,j) += wfc[cc]/dz*dsrc;
            if (jck) {
              Real e0u, ciu, cou, e0l, cil, col;
              jcof(cc, i+1, dt_u, muc[cc], e0u, ciu, cou);
              jcof(cc, i, dt_l, muc[cc], e0l, cil, col);
              const Real pl = BFaceW(kru, kro, bface_on);   // ds_l/dB_i
              const Real pu = BFaceW(kro, kru, bface_on);   // ds_u/dB_{i+1}
              const Real ff = (dtc > 0.0) ? (dt_l/dtc) : 0.5;
              const Real dfl = (1.0 - ff)*pl + ff*(1.0 - pu);
              const Real dfu = 1.0 - dfl;
              const Real dbl = ckdb_g(m,b,i,k,j);
              const Real dbu = ckdb_g(m,b,i+1,k,j);
              const Real wu = wfc[cc]/dzf;
              const Real wl = wfc[cc]/dz;
              // the upper half lies inside cell i+1 and is that cell's row
              CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,i+1),
                  wu*(e0u*jdn[cc] - (ciu*pu + cou*dfu))*dbu);
              CkJacQAdd(&ckjac_g(m,3*blk+0,k,j,i+1),
                  -wu*(ciu*(1.0 - pu) + cou*dfl)*dbl);
              Real au = (1.0 - e0u)*jdn[cc] + ciu*pu + cou*dfu;
              Real al = ciu*(1.0 - pu) + cou*dfl;
              // the face mixing multiplies the ray by (1 + beta): the probe's
              // u_below is frozen, so this factor is the whole of its derivative
              if (SPH) {
                au *= (1.0 + bt_d);
                al *= (1.0 + bt_d);
              }
              CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,i),
                  wl*(e0l*al - (cil*dfl + col*pl))*dbl);
              CkJacQAdd(&ckjac_g(m,3*blk+2,k,j,i),
                  wl*(e0l*au - (cil*dfu + col*(1.0 - pl)))*dbu);
              jdn[cc] = (1.0 - e0l)*al + cil*dfl + col*pl;
            }
          }
          if (BSP && !ckfus) Krs[cc][i] = static_cast<RtF>(kro);
          kfar[cc] = kro;
          // Direct beam, UNCHANGED: it crosses whole cells and carries no
          // source, so the layer construction does not touch it.  Deposit the
          // flux DIFFERENCE across the cell, not kappa rho F exp(-tau) at one
          // face: the latter is what the grey picket fence does, and it
          // under-deposits badly once a layer is not thin (the ratio of
          // deposited to absorbed is u e^-u/(1 - e^-u), 0.95 at u = 0.1 but 0.58
          // at u = 1), losing a quarter of the incident flux down a column with
          // u ~ 0.5 -- measured against Exo-FMS on an identical column.
          if (ckfus) continue;    // the beam is frozen: see ck_impl_frozen_op
          tausw[cc] += kap*drho;
          if (lit && !BSP) {
            const Real tnew = RT_EXP(-static_cast<RtF>(tausw[cc]*facsw));
            // THE STELLAR BEAM IS NOT TOUCHED BY ck_spherical.  It is a parallel
            // pencil, and a parallel pencil is NOT confined to a spherically
            // diverging column: the column widens with r, the beam does not, so
            // rays that enter the top face leave through the column's sides.  The
            // power a column intercepts at its TOP face is therefore NOT the
            // power it absorbs, and the only exact local statement is the volume
            // absorption rate kappa rho F* exp(-tau_ray) -- which is what the
            // expression below already is, written as a flux difference across
            // the cell so that the column integral is F* mu0 (1 - e^-tau_tot) by
            // construction (see the note above on u e^-u/(1 - e^-u)).
            //
            // RETRACTED: an earlier revision of this switch carried the beam as a
            // conserved POWER A(ie+1) F* mu0 T(f), on the argument that a column
            // must deposit every watt it intercepts.  It must not.  That form
            // made every dayside column absorb the power crossing the DOMAIN TOP,
            // the planet absorbed pi r_top^2 F* instead of pi r_abs^2 F* with
            // r_abs where the slant optical depth reaches 1 -- an over-heating by
            // (r_top/r_abs)^2 ~ 3 on the production grid.  See
            // tests_ck_sph/README.md, which keeps the retraction.
            //
            // So this is bit-for-bit the plane-parallel expression under BOTH
            // instantiations, and ck_spherical touches the THERMAL two-stream
            // alone.  The slant path stays the plane-parallel dz/mu0 it was.
            Qb_g(m,blk,i,k,j) += (1.0-albedo)*Fstar*ckswf(b)*wgc[cc]/facsw
                      * (transw[cc] - tnew)/dz;
            transw[cc] = tnew;
          }
        }
      }
      // the lower half of cell icut, from its centre down to the cut face
      {
        const Real dz = dx1(m,k,j,icut);
        for (int cc=0; cc<NC; ++cc) {
          const Real bcut = Bb_g(m,bandc[cc],icut,k,j);
          Real dsrc;
          step(cc, icut, 0.5*kfar[cc]*dz, bcut, bcut, Idn[cc],
               dsrc, false, false);
          Src_g(m,blk,icut,k,j) += wfc[cc]/dz*dsrc;
          I_down[cc][icut] = Idn[cc];
          if (jck) {
            Real e0j, cij, coj;
            jcof(cc, icut, 0.5*kfar[cc]*dz, muc[cc], e0j, cij, coj);
            CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,icut),
                (wfc[cc]/dz)*(e0j*jdn[cc] - (cij + coj))
                *ckdb_g(m,bandc[cc],icut,k,j));
          }
        }
      }
      // ================= problem/ck_beam_sph: THE PSEUDO-SPHERICAL BEAM =====
      // Run here, after the down-sweep has left (kappa rho)_j in Krs for every
      // shell of the correlated-k domain.  See the flag's note for the geometry;
      // the sketch is in tests_ck_sph/README_beam.md.
      //
      //        star
      //          \   theta0                     mu0 >= 0: the ray only climbs,
      //           \                             shells j >= f, weight 1.
      //   r_top ---\-------------------
      //             \   * target r_f          mu0 < 0: NEAR leg r_f -> b and FAR
      //   ...  -------\-*-----------------    leg b -> r_top, so the shells
      //                \* b (tangent)         between b and r_f are crossed
      //   r_cut --------o----------------     TWICE (weight 2) and those above
      //         (opaque below; dark if        r_f once.  Dark if b <= r_cut.
      //          b <= r_cut)
      //
      // ONE sqrt per (face, shell) pair, hoisted out of the g-point loop because
      // the impact parameter does not depend on the chain.  Faces are walked from
      // the top down so the transmission of the face above is already in hand and
      // the deposit is the flux difference across the cell.
      // ck_impl_frozen_op: at frozen opacity this whole ray integration is
      // temperature-independent, so a frozen pass keeps the deposit of pass 0
      if (BSP && lit_sph && !ckfus) {
        const Real rcut = X1F(m,icut);
        RtF thi[NC];
        Real tauh[NC];
        // tau at the top face is 0 for mu0 >= 0: nothing above the domain (see
        // ck_beam_sph).  For mu0 < 0 the ray reaches the top face through its
        // tangent point, near + far leg, which is BeamTopTau.
        for (int cc=0; cc<NC; ++cc) {
          tauh[cc] = 0.0;
          thi[cc] = static_cast<RtF>(1.0);
        }
        if (mu0 < 0.0) BeamTopTau(rcut, sinz, icut, ie, m, X1F, Krs, tauh, thi);
        for (int i=ie; i>icut-1; --i) {
          // ---- tau_ray at the LOWER face of cell i
          const Real rf = X1F(m,i);
          const Real bb = rf*sinz;               // impact parameter
          const Real b2 = bb*bb;
          Real taul[NC];
          for (int cc=0; cc<NC; ++cc) taul[cc] = 0.0;
          bool dark = false;
          int jlo = i;
          if (mu0 < 0.0) {
            if (bb <= rcut) {
              dark = true;            // the ray passes through the body
            } else {
              jlo = icut;                        // the shell holding the tangent
              for (int jj=i-1; jj>=icut; --jj) {
                if (X1F(m,jj) <= bb) {
                  jlo = jj;
                  break;
                }
              }
            }
          }
          if (!dark) {
            const Real rl = X1F(m,jlo);
            Real prev = (rl*rl > b2) ? sqrt(rl*rl - b2) : 0.0;
            for (int jj=jlo; jj<ie+1; ++jj) {
              const Real ru = X1F(m,jj+1);
              const Real cur = sqrt(ru*ru - b2);
              const Real ds = (jj < i) ? 2.0*(cur - prev) : (cur - prev);
              prev = cur;
              for (int cc=0; cc<NC; ++cc) {
                taul[cc] += ds*static_cast<Real>(Krs[cc][jj]);
              }
            }
          }
          // ---- deposit in cell i: the flux difference along the ray, divided
          // by the path length the two optical depths themselves imply.
          // Identical to the plane-parallel expression when tau = tau_vert/mu0.
          for (int cc=0; cc<NC; ++cc) {
            const Real dtl = dark ? 1.0e30 : (taul[cc] - tauh[cc]);
            const RtF tlo = dark ? static_cast<RtF>(0.0)
                                 : RT_EXP(-static_cast<RtF>(taul[cc]));
            const Real dif = static_cast<Real>(thi[cc]) - static_cast<Real>(tlo);
            // dif/dtau -> e^-tau as dtau -> 0; the guard is the same 1e-3 the
            // layer coefficients use
            const Real fac = (fabs(dtl) > 1.0e-3)
                ? (dif/dtl)
                : (static_cast<Real>(thi[cc])*(1.0 - 0.5*dtl));
            Qb_g(m,blk,i,k,j) += (1.0-albedo)*Fstar*ckswf(bandc[cc])*wgc[cc]
                               * fac*static_cast<Real>(Krs[cc][i]);
            tauh[cc] = dark ? 1.0e30 : taul[cc];
            thi[cc] = tlo;
          }
          // the beam is dead: everything below it gets nothing
          Real tmin = tauh[0];
          for (int cc=1; cc<NC; ++cc) tmin = (tauh[cc] < tmin) ? tauh[cc] : tmin;
          if (tmin > 60.0) break;
        }
      }
      // Bottom of the CORRELATED-K DOMAIN, not of the column, and now AT the cut
      // face.  There the grey optical depth is of order 1e4, so the layer is
      // thermalised to e^-tau and the upward intensity is its own Planck
      // function; the planet's internal flux is delivered here as an extra
      // band-weighted source.  Below the cut nothing radiative is applied.
      //
      // ck_spherical NOTE.  The internal flux is delivered in the CUT face's own
      // frame, so under the spherical form the luminosity it represents is
      // A(icut) sigma Tint^4 and what emerges at the top is diluted to
      // (r_cut/r_top)^2 of it -- Tint is a boundary datum at the cut, not at the
      // top.  That is the physically consistent reading of "the internal flux
      // enters at the bottom of the radiative domain"; it is also a real change
      // in the emergent internal flux on a deep grid, and it is measured in
      // tests_ck_sph/README.md rather than quietly rescaled here.
      const Real bt_cut = SPH ? BTF(m,k,j,icut,icut) : 0.0;
      RtF I_up[NC];
      for (int cc=0; cc<NC; ++cc) {
        const int b = bandc[cc];
        const Real Iint_b = (int_at_cut ? boltz_sigma/M_PI*Tint4
                          * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
        I_up[cc] = static_cast<RtF>(Bb_g(m,b,icut,k,j) + Iint_b);
        if (SPH) {
          // Phi = A(icut) (u_below - d_below), reported per unit area of that
          // same face, so the scale factor is A(icut)/A(icut) = 1
          const Real cmx = MIXF(bt_cut, static_cast<Real>(I_up[cc]),
                                static_cast<Real>(I_down[cc][icut]));
          Fb_g(m,blk,icut,k,j) += wfc[cc]*(static_cast<Real>(I_up[cc])
              - (static_cast<Real>(I_down[cc][icut]) + cmx));
          I_up[cc] += static_cast<RtF>(cmx);   // -> cell icut's own frame
        } else {
          Fb_g(m,blk,icut,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][icut]);
        }
        // the internal flux is a boundary datum, so the only B in the starting
        // intensity is the cut cell's own, scaled by the face mixing
        if (jck) jup[cc] = SPH ? (1.0 - bt_cut) : 1.0;
      }
      // up-sweep: the cut half layer first, on the opacity the down-sweep left
      // in kfar, then centre to centre.  Em_g is now the cell's OWN emission,
      // 2 (wfc/mu) kappa rho B summed over the chains, instead of the
      // centre-to-centre average of two Planck functions.
      {
        const Real dz = dx1(m,k,j,icut);
        for (int cc=0; cc<NC; ++cc) {
          const Real bcut = Bb_g(m,bandc[cc],icut,k,j);
          Real dsrc;
          step(cc, icut, 0.5*kfar[cc]*dz, bcut, bcut, I_up[cc],
               dsrc, true, false);
          Src_g(m,blk,icut,k,j) += wfc[cc]/dz*dsrc;
          Em_g(m,blk,icut,k,j) += 2.0*(wfc[cc]/muc[cc])*kfar[cc]*bcut;
          if (jck) {
            Real e0j, cij, coj;
            jcof(cc, icut, 0.5*kfar[cc]*dz, muc[cc], e0j, cij, coj);
            CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,icut),
                (wfc[cc]/dz)*(e0j*jup[cc] - (cij + coj))
                *ckdb_g(m,bandc[cc],icut,k,j));
            jup[cc] = (1.0 - e0j)*jup[cc] + cij + coj;
          }
        }
      }
      for (int i=icut; i<ie; ++i) {
        const Real rhou = rhoN(m,k,j,i+1);
        const Real dzu = dx1(m,k,j,i+1);
        const Real dzl = dx1(m,k,j,i);
        const Real xTv = xT_g(m,k,j,i+1);
        const Real xPv = xP_g(m,k,j,i+1);
        const int iT = static_cast<int>(xTv);
        const int iP = static_cast<int>(xPv);
        const Real fT = xTv - static_cast<Real>(iT);
        const Real fP = xPv - static_cast<Real>(iP);
        // the factor that turns Phi built on the LOWER side, A_cell(i) F, into
        // Phi/A(i+1), and the frame change A_below/A_above for the up ray
        const Real fsc_u = SPH ? ACC(m,k,j,i)/AFC(m,k,j,i+1) : 1.0;
        const Real rat_u = SPH ? ACC(m,k,j,i)/ACC(m,k,j,i+1) : 1.0;
        for (int cc=0; cc<NC; ++cc) {
          const int b = bandc[cc];
          const Real kru = krof(cc, i+1, b, iT, fT, iP, fP, rhou, false);
          const Real krl = kfar[cc];
          const Real bl = Bb_g(m,b,i,k,j);
          const Real bu = Bb_g(m,b,i+1,k,j);
          const Real dt_l = 0.5*krl*dzl;
          const Real dt_u = 0.5*kru*dzu;
          const Real s_l = BFace(kru, krl, bu, bl, bface_on);
          const Real s_u = BFace(krl, kru, bl, bu, bface_on);
          const Real dtc = dt_l + dt_u;
          const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dt_l/dtc))
                                       : (0.5*(s_l + s_u));
          Real dsrc;
          step(cc, i, dt_l, s_l, s_f, I_up[cc], dsrc, false, false);
          Src_g(m,blk,i,k,j) += wfc[cc]/dzl*dsrc;
          if (SPH) {
            // A D CONTINUOUS, ENFORCED, NOT INFERRED.  c is the one the down ray
            // used (Cmx), formed from the probe's u_below; forming a second c
            // here from the up ray's own u_below would make the two sides of the
            // face disagree by beta times the probe's lag, and the deposit would
            // then telescope only to that accuracy (MEASURED: 8.8 % of the
            // emergent luminosity over the production column).  So the up ray's
            // frame change is taken straight from the conservation law,
            //     A_above (u_a - d_a) = A_below (u_b - d_b),
            // which is what the mixing is FOR.  Exact to round-off whatever the
            // probe did; the probe's lag survives only as an O(beta) shift in
            // u_a itself, i.e. in S, which nothing conserves.
            const Real db = static_cast<Real>(I_down[cc][i+1] + Cmx[cc][i+1]);
            const Real dif = static_cast<Real>(I_up[cc]) - db;
            Fb_g(m,blk,i+1,k,j) += wfc[cc]*dif*fsc_u;
            I_up[cc] = static_cast<RtF>(static_cast<Real>(I_down[cc][i+1])
                                        + dif*rat_u);
          } else {
            Fb_g(m,blk,i+1,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][i+1]);
          }
          step(cc, i+1, dt_u, s_f, s_u, I_up[cc], dsrc, true, false);
          Src_g(m,blk,i+1,k,j) += wfc[cc]/dzu*dsrc;
          Em_g(m,blk,i+1,k,j) += 2.0*(wfc[cc]/muc[cc])*kru*bu;
          if (jck) {
            Real e0l, cil, col, e0u, ciu, cou;
            jcof(cc, i, dt_l, muc[cc], e0l, cil, col);
            jcof(cc, i+1, dt_u, muc[cc], e0u, ciu, cou);
            const Real pl = BFaceW(kru, krl, bface_on);   // ds_l/dB_i
            const Real pu = BFaceW(krl, kru, bface_on);   // ds_u/dB_{i+1}
            const Real ff = (dtc > 0.0) ? (dt_l/dtc) : 0.5;
            const Real dfl = (1.0 - ff)*pl + ff*(1.0 - pu);
            const Real dfu = 1.0 - dfl;
            const Real dbl = ckdb_g(m,b,i,k,j);
            const Real dbu = ckdb_g(m,b,i+1,k,j);
            const Real wl = wfc[cc]/dzl;
            const Real wu = wfc[cc]/dzu;
            CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,i),
                wl*(e0l*jup[cc] - (cil*pl + col*dfl))*dbl);
            CkJacQAdd(&ckjac_g(m,3*blk+2,k,j,i),
                -wl*(cil*(1.0 - pl) + col*dfu)*dbu);
            Real bd = (1.0 - e0l)*jup[cc] + cil*pl + col*dfl;
            Real bpu = cil*(1.0 - pl) + col*dfu;
            // the up ray's frame change A_below/A_above, taken from the
            // conservation law: I_down and Cmx are frozen, so it is a factor
            if (SPH) {
              bd *= rat_u;
              bpu *= rat_u;
            }
            CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,i+1),
                wu*(e0u*bpu - (ciu*dfu + cou*pu))*dbu);
            CkJacQAdd(&ckjac_g(m,3*blk+0,k,j,i+1),
                wu*(e0u*bd - (ciu*dfl + cou*(1.0 - pu)))*dbl);
            jup[cc] = (1.0 - e0u)*bpu + ciu*dfu + cou*pu;
          }
          kfar[cc] = kru;
        }
      }
      // the top half of cell ie, back out through the top face
      {
        const Real dz = dx1(m,k,j,ie);
        // at the TOP face A_above is the face area itself: the boundary datum
        // lives in the face's own frame, so fsc and the frame change coincide
        const Real fsc_t = SPH ? ACC(m,k,j,ie)/AFC(m,k,j,ie+1) : 1.0;
        for (int cc=0; cc<NC; ++cc) {
          const Real btop = Bb_g(m,bandc[cc],ie,k,j);
          Real dsrc;
          step(cc, ie, 0.5*kfar[cc]*dz, btop, btop, I_up[cc], dsrc, false, false);
          Src_g(m,blk,ie,k,j) += wfc[cc]/dz*dsrc;
          if (jck) {
            Real e0j, cij, coj;
            jcof(cc, ie, 0.5*kfar[cc]*dz, muc[cc], e0j, cij, coj);
            CkJacQAdd(&ckjac_g(m,3*blk+1,k,j,ie),
                (wfc[cc]/dz)*(e0j*jup[cc] - (cij + coj))
                *ckdb_g(m,bandc[cc],ie,k,j));
          }
          if (SPH) {
            const Real db = static_cast<Real>(I_down[cc][ie+1] + Cmx[cc][ie+1]);
            const Real dif = static_cast<Real>(I_up[cc]) - db;
            Fb_g(m,blk,ie+1,k,j) += wfc[cc]*dif*fsc_t;
            I_up[cc] = static_cast<RtF>(static_cast<Real>(I_down[cc][ie+1])
                                        + dif*fsc_t);
          } else {
            Fb_g(m,blk,ie+1,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][ie+1]);
          }
        }
      }
    }
  };
  for (int b0=0; b0<nblk; b0+=nbc) {
    par_for("rt_chain_ck", DevExeSpace(), 0, nmb1, b0, std::min(nblk, b0+nbc)-1,
            ks, ke, js, je, chain_body);
  }
}

//----------------------------------------------------------------------------------------
//! \fn TsrtGreyRt
//! \brief the grey two-stream column launch of picket_fence_two_stream_RT_pass
//! (its launch_grey_rt), one instantiation per radial tier.  A function template
//! rather than a generic lambda only because nvcc forbids extended lambdas inside
//! generic lambdas; ctx_ holds references to the locals of
//! picket_fence_two_stream_RT_pass the body reads (see grey_rt_ctx there).
template <typename Ctx, typename NnTag>
void TsrtGreyRt(Ctx &ctx_, NnTag nn_tag) {
  auto &albedo = std::get<0>(ctx_);
  auto &ap = std::get<1>(ctx_);
  auto &area1 = std::get<2>(ctx_);
  auto &bdt = std::get<3>(ctx_);
  auto &boltz_sigma = std::get<4>(ctx_);
  auto &demax_grey = std::get<5>(ctx_);
  auto &dx1 = std::get<6>(ctx_);
  auto &eiN = std::get<7>(ctx_);
  auto &eos = std::get<8>(ctx_);
  auto &Fstar = std::get<9>(ctx_);
  auto &gm1 = std::get<10>(ctx_);
  auto &grav = std::get<11>(ctx_);
  auto &grav_pmass = std::get<12>(ctx_);
  auto &iap = std::get<13>(ctx_);
  auto &ie = std::get<14>(ctx_);
  auto &Iint = std::get<15>(ctx_);
  auto &is = std::get<16>(ctx_);
  auto &je = std::get<17>(ctx_);
  auto &js = std::get<18>(ctx_);
  auto &ke = std::get<19>(ctx_);
  auto &ks = std::get<20>(ctx_);
  auto &mbpanel_ = std::get<21>(ctx_);
  auto &met = std::get<22>(ctx_);
  auto &mug = std::get<23>(ctx_);
  auto &nanrep_g = std::get<24>(ctx_);
  auto &nchain_rt = std::get<25>(ctx_);
  auto &nclip_grey = std::get<26>(ctx_);
  auto &nmb1 = std::get<27>(ctx_);
  auto &nrcnt = std::get<28>(ctx_);
  auto &nrrec = std::get<29>(ctx_);
  auto &omega = std::get<30>(ctx_);
  auto &Rgas = std::get<31>(ctx_);
  auto &rhoN = std::get<32>(ctx_);
  auto &stcl_g = std::get<33>(ctx_);
  auto &Teq = std::get<34>(ctx_);
  auto &tide = std::get<35>(ctx_);
  auto &Tint4 = std::get<36>(ctx_);
  auto &Tirr4 = std::get<37>(ctx_);
  auto &u0 = std::get<38>(ctx_);
  auto &use_cubed_sphere_ = std::get<39>(ctx_);
  auto &use_spherical_polar = std::get<40>(ctx_);
  auto &volume = std::get<41>(ctx_);
  auto &wg = std::get<42>(ctx_);
  auto &wtemp_ = std::get<43>(ctx_);
  auto &x1f_ = std::get<44>(ctx_);
  auto &x1v_ = std::get<45>(ctx_);
  auto &x2v_ = std::get<46>(ctx_);
  auto &x3v_ = std::get<47>(ctx_);
  // the compile-time switches picket_fence_two_stream_RT_pass declares next to the
  // launcher (constants, so they are not part of grey_rt_ctx): keep the two in step
  const bool correct_spherical = false;
  const bool test_oned = false;
  constexpr int NN = decltype(nn_tag)::value;
  par_reduce_clip3("2stream_rt", 0, nmb1, ks, ke, js, je, nclip_grey,
  KOKKOS_LAMBDA(const int m, const int k, const int j, int &nc) {
//        ScrArray1D<Real> tau_down_r_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> F_v_down_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> B(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> I_ir_down_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> I_ir_up_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> F_ir_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> Q_v(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> kapJ_ir(member.team_scratch(scr_level), n1);
      Real tau_down_r_f[NN];
//        Real F_v_down_f[NN];
      Real B[NN];
      Real F_ir_f[NN];
      Real Q_v[NN];
//        Real kapJ_ir[NN];

      Real x2v = x2v_(m,j);
      Real x3v = x3v_(m,k);

      Real rtop = x1v_(m,ie+1);
      Real rbot = x1v_(m,is);

      Real lam, phi, theta;
      if (use_spherical_polar) {
        theta = x2v;
        lam = -theta+M_PI/2.0;
        phi = x3v-M_PI;
      } else if (use_cubed_sphere_) {
        CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
      } else {
        lam = x3v*iap;
        theta = -lam+M_PI/2.0;
        phi = x2v*iap;
      }
      Real ex = sin(theta)*cos(phi);
      Real ex0 = 1.0;
      Real mu0 = ex*ex0;
      if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);

      Real mus = (mu0 > 0.0) ? mu0 : 0.0;
      Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
      Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
      get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);

      // 3 V Bands
      // top
      Real rho = rhoN(m,k,j,ie+1);
      Real p = PresFromEint(eos,gm1,rho,eiN(m,k,j,ie+1));
      Real T = TempKelvin(eos,Rgas,rho,eiN(m,k,j,ie+1),p);
      bool badtop = RTBadState(p, T);    // see RTBadState: one NaN kills the column
      if (badtop) {
        Kokkos::atomic_fetch_add(&stcl_g(0), 1);
        T = 0.0;
        p = 0.0;
      }
      B[ie+1] = badtop ? 0.0 : boltz_sigma/M_PI*SQR(SQR(T));
      Real kapr = 0.0;
      if (!badtop) get_kapr(T, p, met, kapr);
      Real tau_r_f = RTTopDtau(kapr, p,
                               EffGravAt(grav, ap, rtop, grav_pmass, omega, mu0,
                                         tide));
      tau_down_r_f[ie+1] = tau_r_f;
      Real drtop = (kapr*rho != 0.0) ? tau_r_f/(kapr*rho) : 0.0;
      Real delta = drtop/rtop;
      Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
      fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
      Real tausl = tau_r_f*fac;
      Real trans1 = exp(-gamv1*tau_down_r_f[ie+1]*fac);
      Real trans2 = exp(-gamv2*tau_down_r_f[ie+1]*fac);
      Real trans3 = exp(-gamv3*tau_down_r_f[ie+1]*fac);
      // beam transmission at the face ABOVE the cell being filled, carried down the
      // sweep so that differencing the flux across a cell costs no extra exp
      Real trp1 = trans1;
      Real trp2 = trans2;
      Real trp3 = trans3;
//        F_v_down_f[ie+1] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
//        F_v_down_f(ie+1) = (mu0 > 0.0)? F_v_down_f(ie+1) : 0.0;
      // down-sweep
      for (int i=ie; i>is-1; --i) {
        Real rho = rhoN(m,k,j,i);
        Real p, T;
        PresTempFromEint(eos,gm1,Rgas,rho,eiN(m,k,j,i),
                         TGuess(wtemp_, m, k, j, i),p,T);
        bool badcell = RTBadState(p, T); // see RTBadState: one NaN kills the column
        if (badcell) {
          Kokkos::atomic_fetch_add(&stcl_g(0), 1);
          T = 0.0;
          p = 0.0;
        }
        B[i] = badcell ? 0.0 : boltz_sigma/M_PI*SQR(SQR(T));
        Real kapr = 0.0;
        if (!badcell) get_kapr(T, p, met, kapr);
        Real dr = dx1(m,k,j,i);
        tau_down_r_f[i] = tau_down_r_f[i+1] + kapr*rho*dr;
        Real r = x1f_(m,i);
////          Real delta = (drtop+(rtop-r))/r;
//          Real delta = dr/r;
//          Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
//          tausl += kapr*rho*r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
        Real fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
        Real trans1 = exp(-gamv1*tau_down_r_f[i]*fac);
        Real trans2 = exp(-gamv2*tau_down_r_f[i]*fac);
        Real trans3 = exp(-gamv3*tau_down_r_f[i]*fac);
//          Real trans1 = exp(-gamv1*tausl);
//          Real trans2 = exp(-gamv2*tausl);
//          Real trans3 = exp(-gamv3*tausl);
//          F_v_down_f[i] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
        Real mucr = 0.0; //sqrt(1.0-SQR(r0/r));
        // Deposit the flux DIFFERENCE across the cell. The old form,
        // kappa rho F exp(-tau) with tau at the lower face, is right only for a thin
        // layer: it returns u e^-u / (1 - e^-u) of what the cell actually absorbs,
        // with
        // u = dtau/mu, which is 0.95 at u = 0.1 but 0.58 at u = 1. On this grid,
        // about
        // 0.46 scale heights per cell, that put dtau/mu near one wherever it mattered
        // and lost about 24 % of the incident stellar flux. Found by comparing the
        // correlated-k version of the same expression against Exo-FMS on an identical
        // column. Written this way the column integral telescopes to
        // mu F (1 - e^-tau_total) exactly, and it still reduces to the old expression
        // as dtau -> 0.
        Real Qv = (1.0-albedo)*Fstar*(1.0/3.0)
                * ((trp1-trans1)+(trp2-trans2)+(trp3-trans3))/(fac*dr);
        Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
        trp1 = trans1;
        trp2 = trans2;
        trp3 = trans3;
      }

      // 2 IR Bands x two quadrature points, interleaved, in blocks of NC.
      //
      // Each (band, quadrature) combination is an independent pair of linear
      // recurrences in radius, and running them one after another leaves the
      // wavefront
      // stalled on a single dependency chain: this kernel has only nmb*nx3*nx2/64
      // wavefronts for 912 SIMDs, so there is no other wave to hide that latency and
      // VALUBusy sits near 3 %. Stepping NC combinations inside one radial loop gives
      // the chain NC independent strands to overlap, and blocking keeps the private
      // I_ir_down_c footprint at NC columns however many chains are requested.
      //
      // This is the grey path, so nchain_rt is 4 and there is exactly one block; the
      // blocking survives because the correlated-k kernel shares the structure.
      constexpr int NC = RT_NB;
      const int nblk_rt = (nchain_rt + NC - 1)/NC;
      for (int i=is; i<ie+2; ++i) {
        F_ir_f[i] = 0.0;
      }

      for (int blk=0; blk<nblk_rt; ++blk) {
        Real gamirc[NC], fbc[NC], muggc[NC], wggc[NC];
        for (int cc=0; cc<NC; ++cc) {
          const int c = blk*NC + cc;
          const int n = (c/2) % 2;
          const int vir = c % 2;
          muggc[cc] = mug[n];
          wggc[cc] = wg[n];
          gamirc[cc] = (vir == 0) ? gamir1 : gamir2;
          fbc[cc] = (vir == 0) ? beta : (1.0-beta);
        }
        Real I_ir_down_c[NC][NN];

        // top
        for (int cc=0; cc<NC; ++cc) {
          Real dtauir = gamirc[cc]*tau_down_r_f[ie+1];
          Real trans = exp(-dtauir/muggc[cc]);
          I_ir_down_c[cc][ie+1] = (1.0-trans)*(fbc[cc]*B[ie+1]);
        }
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real dtau_i = tau_down_r_f[i]-tau_down_r_f[i+1];
          for (int cc=0; cc<NC; ++cc) {
            Real dtauir = gamirc[cc]*dtau_i;
            Real x = dtauir/muggc[cc];
            Real e0 = -expm1(-x);
            Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
            Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
            I_ir_down_c[cc][i] = (1.0-e0)*I_ir_down_c[cc][i+1]
                               + alp*fbc[cc]*B[i+1] + bet*fbc[cc]*B[i];
          }
        }

        // bottom
        Real I_ir_up_c[NC];
        for (int cc=0; cc<NC; ++cc) {
          I_ir_up_c[cc] = Iint + I_ir_down_c[cc][is];
          Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][is];
          Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
          F_ir_f[is] += (F_ir_up_f - F_ir_down_f);
        }
        // up-sweep, accumulating the band flux as it goes
        for (int i=is+1; i<ie+2; ++i) {
          Real dtau_i = tau_down_r_f[i-1]-tau_down_r_f[i];
          for (int cc=0; cc<NC; ++cc) {
            Real dtauir = gamirc[cc]*dtau_i;
            Real x = dtauir/muggc[cc];
            Real e0 = -expm1(-x);
            Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
            Real gm = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
            I_ir_up_c[cc] = (1.0-e0)*I_ir_up_c[cc]
                          + bet*fbc[cc]*B[i] + gm*fbc[cc]*B[i-1];
            Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][i];
            Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
            F_ir_f[i] += (F_ir_up_f - F_ir_down_f);
          }
        }
      }

//        // Sync all threads in the team so that scratch memory is consistent
//        member.team_barrier();

//        par_for_inner(member, is, ie, [&](const int i) {
      for (int i=is; i<ie+1; ++i) {
        // source term as flux divergence
        Real area_t = area1(m,k,j,i+1);
        Real area_b = area1(m,k,j,i);
        Real vol = volume(m,k,j,i);
          Real Ft = F_ir_f[i+1];//-F_v_down_f(i+1);
          Real Fb = F_ir_f[i];//-F_v_down_f(i);
        Real src = -(Ft-Fb)/dx1(m,k,j,i);
        if (correct_spherical) {
          src = -(Ft*area_t-Fb*area_b)/vol;
        }
          src += Q_v[i];
        Real du_flux = src*bdt;

//          // source term semi-implicit
//          Real p = PresFromEint(eos,gm1,rhoN(m,k,j,i),eiN(m,k,j,i));
//          Real rho = rhoN(m,k,j,i);
//          Real T = TempKelvin(eos,Rgas,rho,eiN(m,k,j,i),p);
//          Real kapr;
//          get_kapr(T, p, met, kapr);
//          Real cv = Rgas*rho*igm1;
//          Real e0 = eos.IsGeneral() ? eiN(m,k,j,i) : cv*T;
//          Real kk = 0.0;
//          Real bb = du_flux + e0;
////          Real bb = Q_v(i)*bdt + e0;
//          for (int vir=0; vir<2; ++vir) {
//            Real gamir, fb;
//            if (vir == 0) {
//              gamir = gamir1;
//              fb = beta;
//            } else {
//              gamir = gamir2;
//              fb = 1.0-beta;
//            }
//            kk += -4.0*M_PI*gamir*kapr*rho*fb*boltz_sigma/M_PI*bdt;
//            bb += 4.0*M_PI*gamir*kapr*rho*fb*B[i]*bdt;
//          }
////          bb += 4.0*M_PI*rho*kapJ_ir(i)*bdt;
//          int ierr=0;
//          Real e = e0;
//          // Newton-Raphson. A general EOS has no e = c_v T with constant c_v, so the
//          // iteration runs on the internal energy directly rather than on T:
//          // F(e) = e - kk T(e)^4 - bb, with dT/de = temp_cgs/(d c_v) since T is in K.
//          for (int n=0; n<100; ++n) {
//            Real de;
//            if (eos.IsGeneral()) {
//              Real dTde = eos.temp_cgs/(rho*eos.SpecificHeatCv(rho,e));
//              de = e - kk*SQR(SQR(T)) - bb;
//              e -= de / (1.0 - 4.0*kk*T*T*T*dTde);
//              T = eos.Temperature(rho,e)*eos.temp_cgs;
//            } else {
//              e = cv*T;
//              de = e - kk*SQR(SQR(T)) - bb;
//              T -= de / (cv - 4.0*kk*T*T*T);
//            }
//            if (T < 0.0) {
//              e = e0;
//              ierr = 1;
//              break;
//            }
//            if (fabs(de) <= 1.0e-10*e)
//              break;
//          }
//          Real du_src = e-e0;
//
//          Real du = (fabs(du_flux) < e0 && ierr == 1) ? du_flux : du_src;
        Real du = du_flux;
        if (demax_grey > 0.0) {
          const Real dl = LimitRTSource(du, eiN(m,k,j,i), demax_grey);
          if (dl != du) { ++nc; du = dl; }
        }
        u0(m,IEN,k,j,i) += du;
        // --- the nan_report catcher.  Record the cell the grey apply just made
        // non-finite or non-positive, with the state that produced it.
        if (nanrep_g) {
          const Real enew = u0(m,IEN,k,j,i);
          if (!(enew > 0.0) || !isfinite(enew) || !isfinite(du)) {
            if (Kokkos::atomic_fetch_add(&nrcnt(0), 1) == 0) {
              Real pc, tc, pm1, tm1, pp1, tp1, kc;
              PresTempFromEint(eos, gm1, Rgas, rhoN(m,k,j,i), eiN(m,k,j,i),
                               TGuess(wtemp_, m, k, j, i), pc, tc);
              PresTempFromEint(eos, gm1, Rgas, rhoN(m,k,j,i-1), eiN(m,k,j,i-1),
                               TGuess(wtemp_, m, k, j, i-1), pm1, tm1);
              PresTempFromEint(eos, gm1, Rgas, rhoN(m,k,j,i+1), eiN(m,k,j,i+1),
                               TGuess(wtemp_, m, k, j, i+1), pp1, tp1);
              get_kapr(tc, pc, met, kc);
              nrrec(0) = static_cast<Real>(m);
              nrrec(1) = static_cast<Real>(k);
              nrrec(2) = static_cast<Real>(j);
              nrrec(3) = static_cast<Real>(i);
              nrrec(4) = rhoN(m,k,j,i);
              nrrec(5) = eiN(m,k,j,i);
              nrrec(6) = tc;
              nrrec(7) = kc*rhoN(m,k,j,i);
              nrrec(8) = du;
              nrrec(9) = du_flux;
              nrrec(10) = enew;
              nrrec(11) = eiN(m,k,j,i-1);
              nrrec(12) = tm1;
              nrrec(13) = eiN(m,k,j,i+1);
              nrrec(14) = tp1;
              nrrec(15) = tau_down_r_f[i];
            }
          }
        }
      }
//        });

//        // Sync all threads in the team so that scratch memory is consistent
//        member.team_barrier();

  });
}

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_RT_CHAIN_HPP_
