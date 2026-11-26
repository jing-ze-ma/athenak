//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_ledger.cpp
//! \brief ledgers of L_z and boundary energy, timers, end-of-run reports.
//!
//!  called from: rad_vet.cpp, rad_vet_gas.cpp, rad_vet_imex.cpp, rad_vet_oldvector.cpp,
//!    rad_vet_tasks.cpp, rad_vet_update.cpp, driver/driver.cpp, outputs/eventlog.cpp,
//!    pgen/star_vet.cpp, pgen/tests/rad_vet_tests2.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_fct.cpp, rad_vet_init.cpp

#include <cstdio>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "rad_vet/rad_vet_reduce.hpp"
#include "utils/imex_sdirk2.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ReportCounters
//! \brief one line at the end of the run: cost and health of the implicit energy
//! exchange solve (failures: cells that hit the iteration cap above the tolerance)

void RadiationVet::ReportCounters() {
  if (!coupling) return;
  cnt.modify_device();
  cnt.sync_host();
  Real nsolve = cnt.h_view(VET_CNT_NSOLVE);
#if MPI_PARALLEL_ENABLED
  {
    Real loc[VET_NCNT], glob[VET_NCNT];
    for (int n=0; n<VET_NCNT; ++n) {loc[n] = cnt.h_view(n);}
    MPI_Reduce(loc, glob, VET_NCNT, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    Real lmax = cnt.h_view(VET_CNT_ITMAX), gmax;
    MPI_Reduce(&lmax, &gmax, 1, MPI_ATHENA_REAL, MPI_MAX, 0, MPI_COMM_WORLD);
    for (int n=0; n<VET_NCNT; ++n) {cnt.h_view(n) = glob[n];}
    cnt.h_view(VET_CNT_ITMAX) = gmax;
    nsolve = cnt.h_view(VET_CNT_NSOLVE);
  }
#endif
  if (global_variable::my_rank != 0) return;
  Real mean = (nsolve > 0.0) ? (cnt.h_view(VET_CNT_ITSUM)/nsolve) : 0.0;
  std::cout << "<rad_vet> implicit energy solve: cells=" << nsolve
            << " iterations mean=" << mean << " max=" << cnt.h_view(VET_CNT_ITMAX)
            << " FAILURES=" << cnt.h_view(VET_CNT_NFAIL)
            << " (no bracket: " << cnt.h_view(VET_CNT_NBRAK) << ")" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanLzBook
//! \brief clean_cs_lz_abs = 2: books the cells' sum s of the L_z the gas did not take (a
//! write-back, or with start the stage start's increment) as its gas momentum enters the
//! step: stage-1 part weighted as in the stage-2 start (Heun 1/2), stage-2 parts pending
//! until stage 2 is accepted (s = 0; a backward-Euler redo drops them, Time2Restore).
//! Total and K1's sum in <rad_vet> (restart).  3-D spherical polar: s is the L_z that
//! crossed the boundaries (CleanStagMomentum's dlz), <rad_vet>/clean_lz_ledger.

void RadiationVet::CleanLzBook(const Real s, const bool start) {
  if (t2_inprobe) return;   // a helper solve, put back (Time2RotCorrect)
  if (t2_solve == VET_T2S_STAGE1) {
    cln_lzp = (start ? 0.0 : cln_lzp) + ((t2_tab == 0 || t2_tab == 2) ? 0.5 : t2_w0)*s;
    return;
  }
  if (t2_solve == VET_T2S_STAGE2 && start) {
    cln_lzp += s;
    return;
  }
  cln_lzl += (t2_solve == VET_T2S_STAGE2) ? cln_lzp : cln_lzp + s;
  cln_lzp = 0.0;
  const bool cs = (cln_lza == 2);
  cln_lzpin->SetReal("rad_vet", VetLzKey(cln_lzpin, cs, false), cln_lzl);
  cln_lzpin->SetReal("rad_vet", VetLzKey(cln_lzpin, cs, true), cln_lzs1);
  if (global_variable::my_rank == 0 && cln_verbose > 0 && !cln_insteady) {
    std::printf("CLEAN LZ%s cycle %d t2s %d ledger %.16e\n", cs ? "" : "B",
                pmy_pack->pmesh->ncycle, t2_solve, cln_lzl);
  }
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
//! \fn void RadiationVet::TmrMark
//! \brief implicit_timers: fence, charge the time since the last mark to category c
//! (c < 0: only set the mark)

void RadiationVet::TmrMark(int c) {
  if (!tmr_on) return;
  Kokkos::fence();
  const double t = tmr_t.seconds();
  if (c >= 0) {tmr_acc[c] += t - tmr_last;}
  tmr_last = t;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::EventTotals
//! \brief limiter totals for the event log, summed over ranks (collective).
//! Slots: 0..6 unused (0); 7 opac-Newton face terms dropped; 8 E floor raises;
//! 9 floor energy created (x volume); 10 of which solve undershoot; 11-12 closure clips

void RadiationVet::EventTotals(Real *v) {
  for (int q = 0; q < 13; ++q) {v[q] = 0.0;}
  if (opn_nskip_d.extent_int(0) > 0) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), opn_nskip_d);
    v[7] = h(0);
  }
  if (pos_cnt_d.extent_int(0) >= VET_POS_N) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    v[8] = h(VET_POS_FLR);
    v[9] = h(VET_POS_FLR_UN);
    v[10] = h(VET_POS_FLRS);
  }
  if (acl_cnt_d.extent_int(0) == 2) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), acl_cnt_d);
    v[11] = h(0);
    v[12] = h(1);
  }
#if MPI_PARALLEL_ENABLED
  Real g[6];
  MPI_Allreduce(&v[7], g, 6, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int q = 0; q < 6; ++q) {v[7+q] = g[q];}
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplicitReport
//! \brief end-of-run report: timers and positivity bookings

void RadiationVet::ImplicitReport() {
  // positivity counters and hydro floor energy (collective: same pos keys on all ranks)
  Real hfl = pmy_pack->pmesh->ecounter.flr_de_v;
  {
    auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    for (int q = 0; q < VET_POS_N; ++q) {pos_cnt[q] = hpc(q);}
#if MPI_PARALLEL_ENABLED
    Real g[VET_POS_N];
    MPI_Allreduce(pos_cnt, g, VET_POS_N, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < VET_POS_N; ++q) {pos_cnt[q] = g[q];}
    Real gd = 0.0;
    MPI_Allreduce(&pos_un_drop, &gd, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    pos_un_drop = gd;
    MPI_Allreduce(MPI_IN_PLACE, &hfl, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  }
  if (global_variable::my_rank != 0) return;
  if (tmr_c0 > 0) {
    // seconds per category (see TmrMark)
    std::cout << "<rad_vet> timers from cycle " << tmr_c0 << ": stages=" << tmr_cnt[7]
              << std::endl;
    const char *nm[11] = {"closure", "opacity", "solve_pre", "tensor", "pred",
                          "pass_setup", "krylov", "pass_post", "solve_end", "hyd_c2p",
                          "vet_bvals"};
    std::cout << "<rad_vet> timers(s):";
    for (int q = 0; q < 11; ++q) {std::cout << " " << nm[q] << "=" << tmr_acc[q];}
    std::cout << std::endl;
  }
  std::cout << "<rad_vet> vet-positivity: gas cell-solves=" << pos_cnt[VET_POS_GAS]
            << " energy moved rad->gas=" << pos_cnt[VET_POS_GAS_DE]
            << " | E floor raises cell-solves=" << pos_cnt[VET_POS_FLR]
            << " energy from gas=" << pos_cnt[VET_POS_FLR_DE]
            << " energy created=" << pos_cnt[VET_POS_FLR_UN] << " (code units x volume)"
            << " [in discarded attempts, rolled back: " << pos_un_drop << "]"
            << " taken from the column's radiation (implicit_pos_floor_cons)="
            << pos_cnt[VET_POS_FLR_MV]
            << " | hydro EOS floors energy created=" << hfl << std::endl;
}

namespace {

const Real kT2G = imex_sdirk2::Gamma();   // shared tableau

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanBndAcc
//! \brief the history's x1 boundary energy cln_lacc, booking the flux that moved E: stage
//! 0 starts a backward-Euler step (its solves add h A Phi); stage 1 / 2 the old vector of
//! a hesdirk2 stage, w x (energy of Y2) + A (c1 t2pf1 + c2 t2pf2) at the donors' limiter
//! R; the stage solve adds g dt A Phi

void RadiationVet::CleanBndAcc(int stage, Real c1, Real c2, Real w, Real dt) {
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  if (!(sph_geom || cs_geom)) return;
  if (cln_lacc.extent_int(0) <= nmb1) {
    Kokkos::realloc(cln_lacc, nmb1 + 1, 4, u0.extent_int(2), u0.extent_int(3));
  }
  if (stage == 0 || stage == 1) {cln_ldt = dt;}
  if (stage == 0) {Kokkos::deep_copy(DevExeSpace(), cln_lacc, 0.0); return;}
  const bool two = (stage != 1), ok1 = t2pf1_ok;
  if (t2pf1.extent_int(0) <= nmb1 || (two && !t2pf2_ok)) {cln_ldt = 0.0; return;}
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie;
  auto la = cln_lacc;
  auto p1 = t2pf1, p2 = t2pf2;
  auto g1 = t2f1;
  auto iw_ = iw;
  auto ar = pmy_pack->pcoord->area.x1f;   // (the x1 face area of the clean fluxes)
  const bool lim = ok1;   // (the limiter's R in iw slot Time2FctSlot(0))
  const int sr = Time2FctSlot(0);
  const Real crp = chat/c_light, cc2 = two ? c2 : 0.0;
  par_for("cln_bndacc", DevExeSpace(), 0, nmb1, 0, 1, id.ks, id.ke, id.js, id.je,
  KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
    const int i = s ? ie + 1 : is;
    const Real g = ar(m,k,j,i)*(c1*(ok1 ? p1(m,0,k,j,i) : crp*g1(m,k,j,i)) +
                                cc2*p2(m,0,k,j,i));
    const Real rd = lim ? iw_(m,sr,k,j,(g > 0.0) ? i - 1 : i) : 1.0;
    la(m,2+s,k,j) = w*la(m,s,k,j) + rd*g;
  });
}

//! \fn Real RadiationVet::EbndEsrc
//! \brief ebnd_split: (chat/c) sum esrc V over the active cells, all ranks (0 without
//! esrc)

Real RadiationVet::EbndEsrc() {
  if (!esrc_on) return 0.0;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, js = id.js, ks = id.ks;
  const int ni = id.nx1, nj = id.nx2, nk = id.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  auto es_ = esrc;
  auto &vol = pmy_pack->pcoord->volume;
  Real l = 0.0;
  Kokkos::parallel_reduce("vet_ebnd_src", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nk*nj*ni),
  KOKKOS_LAMBDA(const int n, Real &acc) {
    const int m = n/(nk*nj*ni);
    const int k = (n/(nj*ni))%nk + ks, j = (n/ni)%nj + js, i = n%ni + is;
    acc += es_(m,k,j,i)*vol(m,k,j,i);
  }, Kokkos::Sum<Real>(l));
#if MPI_PARALLEL_ENABLED
  Real g;
  MPI_Allreduce(&l, &g, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  l = g;
#endif
  return (chat/c_light)*l;
}

//! \fn void RadiationVet::EbndStart / EbndStage1 / EbndAccept
//! \brief Ebnd: L at the step start (face flux of K1), of the accepted stage-1 solve
//! (K2), and the step's applied energy at acceptance (after the dt_retry check)

void RadiationVet::EbndStart() {
  if (!ebnd_on) return;
  if (ebnd_split) ebnd_s1 = ebnd_s3_ok ? ebnd_s3 : EbndEsrc();
}

void RadiationVet::EbndStage1() {
  if (!ebnd_on) return;
  if (ebnd_split) ebnd_s2 = EbndEsrc();
}

void RadiationVet::EbndAccept(Real dt, bool t2step) {
  if (!ebnd_on) return;
  // the step's applied energy of the slopes l1 (K1), l2 (K2), l3 (K3)
  auto applied = [&](const Real l1, const Real l2, const Real l3) {
    if (t2step) {
      const Real g = kT2G;
      const Real wy = (t2_tab != 1) ? 0.5 : t2_w0;   // weight of Y1 in the stage-2 start
      return dt*((wy*t2_a21 + t2_c31)*l1 + (wy*g + t2_c32)*l2 + g*l3);
    }
    return dt*l3;
  };
  // clean: f0x1 is the face momentum; the energy through the physical x1 faces is
  // CleanBndAcc's
  Real lc[2] = {0.0, 0.0};
  if (!(cln_ldt > 0.0) && global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d: Ebnd without the step's boundary energy "
                "(CleanBndAcc), its faces take 0\n", pmy_pack->pmesh->ncycle);
  }
  if (cln_ldt > 0.0) {
    auto la = cln_lacc;
    auto &bcs = pmy_pack->pmb->mb_bcs;
    auto &id = pmy_pack->pmesh->mb_indcs;
    Kokkos::parallel_reduce("vet_ebnd_cln", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        DevExeSpace(), {0, id.ks, id.js}, {pmy_pack->nmb_thispack, id.ke + 1, id.je + 1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, Real &si, Real &so) {
      for (int s = 0; s < 2; ++s) {
        const BoundaryFlag f = bcs.d_view(m, s ? BoundaryFace::outer_x1 :
                                                 BoundaryFace::inner_x1);
        if (f == BoundaryFlag::block || f == BoundaryFlag::periodic) continue;
        (s ? so : si) += la(m,s,k,j);
      }
    }, Kokkos::Sum<Real>(lc[0]), Kokkos::Sum<Real>(lc[1]));
    for (int s = 0; s < 2; ++s) lc[s] = CleanReduce(lc[s], false);
  }
  if (ebnd_split) {
    ebnd_acc_in += lc[0];
    ebnd_acc_out += lc[1];
    const Real s3 = EbndEsrc();
    ebnd_acc_src += applied(ebnd_s1, ebnd_s2, s3);
    ebnd_s3 = s3;
    ebnd_s3_ok = true;
    ebnd_acc = ebnd_acc_in - ebnd_acc_out;
    return;
  }
  ebnd_acc += lc[0] - lc[1];
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Report
//! \brief end-of-run counters of the hesdirk2 steps (rank 0)

void RadiationVet::Time2Report() {
  if (time_scheme != VET_TIME_HESDIRK2) return;
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_vet> time_scheme=" << (t2_beonly ? "hesdirk2(be_only)" : "hesdirk2")
            << ": stage steps=" << t2_nstep
            << " backward-Euler steps=" << t2_nbe << " stage fallbacks=" << t2_nfall
            << " vet clips=" << t2_nclip;
  // global sums
  std::cout << " face_fct calls=" << t2_ffct_calls << " skipped(no face flux)="
            << t2_ffct_nop << " K1 flux from U^n=" << t2_ffct_fill << " limited cells="
            << t2_ffct_cells << " max rounds="
            << t2_ffct_rmax << " withheld E x vol=" << t2_ffct_de << " left (Q<=0)="
            << t2_ffct_left << " emptied (E < " << t2_ffct_emin << " E^n) cells="
            << t2_ffct_nemp << " min E/E^n="
            << t2_ffct_qmin << " (rank 0: rounds=" << zl_rounds << " R exchanges="
            << zl_halos;
  std::cout << ")";
  if (vsc_rot && !cln_ring) {
    std::cout << " set-change slope solves=" << vgd_rot_ncorr << " failed="
              << vgd_rot_ncfail << " Picard passes=" << vgd_rot_ncpass;
  }
  // VetGdAlignDt calls that shortened dt (NewTimeStep calls, incl. the start-up ones)
  if (vgd_rot_dt > 0.0) {std::cout << " rotate_align shortened dt=" << vgd_rot_nalign;}
  if (cln_ring) {
    std::cout << " clean_rot_ring guard cells=" << cln_rguard_sum
              << " off-boundary shifts=" << cln_roff_sum;
  }
  std::cout << std::endl;
  if (t2_abe > 0) {
    std::cout << "<rad_vet> time2_auto_be=" << t2_abe << ": skipped-stage BE steps="
              << t2_abe_n << " stage retries=" << t2_abe_ntry << " (succeeded "
              << t2_abe_nok << ")" << std::endl;
  }
  // the wall-time split (rank 0, fenced; accepted cycles by kind)
  std::cout << "<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps " << t2w_clean
            << " | failed stage attempts " << t2w_fail << " | their BE redo " << t2w_redo
            << " | auto_be steps " << t2w_abe << " | other BE steps " << t2w_be
            << " | discarded dt_retry attempts " << t2w_retry << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetReport

void RadiationVet::VetReport() {
  if (global_variable::my_rank != 0 || vet_ncall <= 0.0) return;
  std::cout << "<rad_vet> vet_sc: formal solutions=" << vet_ncall
            << " SC seconds=" << vet_time << " (" << vet_time/vet_ncall << " per call)"
            << " implicit-solve seconds incl. SC=" << vet_itime
            << " SC/(implicit solve without SC)="
            << vet_time/std::fmax(vet_itime - vet_time, 1.0e-300) << std::endl;
  if (vet_x1_nbuild > 0.0) {   // adaptive periodic x1 sweep
    std::cout << "<rad_vet> vet_sc periodic x1 sweep (vet_x1_ptol=" << vet_x1ptol
              << ", cap vet_x1_npass=" << vet_x1npass << "): builds=" << vet_x1_nbuild
              << " passes mean=" << vet_x1_npsum/vet_x1_nbuild << " max=" << vet_x1_npmax
              << " builds that hit the cap=" << vet_x1_ncap << std::endl;
  }
  if (vra_on && vsc_rot) {
    std::cout << "<rad_vet> vet_rotate_avg = " << vra_k << ": snapshot sweeps "
              << vra_nsweep << " (rotations " << vra_nrot << ", inits " << vra_ninit
              << "), " << vra_time << " s incl. offsets (rank 0, fenced)" << std::endl;
  }
  if (impl_hr_sf) {
    std::cout << "<rad_vet> implicit_hr_face=sf: cells with Jh > J (sf face ratios) "
              << hsf_nsf << " of " << hsf_ncell << " cell-solves (rank 0)" << std::endl;
  }
  if (vet_full) {
    std::cout << "<rad_vet> vet_tensor=full: guard triggered in " << vet_nguard << " of "
              << vet_ncell << " cell-calls (eigenvalue floor " << vet_eig_min
              << "); max |D_guarded - K/J|=" << vet_guard_max
              << " (incl. trace normalisation); min D_11=" << vet_d11_min << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRayReport
//! \brief one cost line at the end of the run (rank 0)

void RadiationVet::VetRayReport() {
  if (global_variable::my_rank != 0 || !vgd_on) return;
  std::cout << "<rad_vet> vet_gd: " << vray_ncall << " builds checked; lateral "
            << "reads clamped to the ghost band (rank 0, all sweeps): " << vray_nclamp
            << std::endl;
  std::cout << "<rad_vet> vet_gd: rank 0 fenced seconds: source " << vgd_tsrc
            << ", sweeps " << vgd_tswp << " (of which per-shell halo+walls "
            << vgd_thalo << ", of which MPI wait " << vgd_tmpi
            << "), tables (rotation) " << vgd_ttab
            << "; in-solve rebuilds " << vgd_nrb << " ("
            << vgd_trb << " s)" << std::endl;
  std::cout << "<rad_vet> vet_gd halo (rank 0): exchanges " << vgd_nexch
            << ", bytes sent "
            << vgd_nbyte << " (" << ((vgd_nexch > 0.0) ? vgd_nbyte/vgd_nexch : 0.0)
            << " per exchange), post " << vgd_tpost << " s, partners "
            << vgd_prk.size() << "; band-clamped lateral reads (ALL ranks, all "
            << "sweeps) "
            << vgd_nclamp_all << std::endl;
  if (vgd_csg || vgd_gbc) {   // cs or sp global band
    std::cout << "<rad_vet> vet_gd " << (vgd_csg ? "cs" : "global")
              << " band (rank 0): " << ((vgd_cscomp || vgd_gbc) ? "compact" : "dense")
              << ", intensity values sent "
              << vgd_csvsent << " of dense " << vgd_csvdense << " ("
              << ((vgd_csvdense > 0.0) ? vgd_csvsent/vgd_csvdense : 0.0)
              << "); list builds " << vgd_csnbuild << " (" << vgd_cstbuild << " s)"
              << std::endl;
  }
  if (vgd_csg && vgd_cscomp) {   // where the list builds go
    std::cout << "<rad_vet> cs list builds (rank 0 s): mask kernel "
              << vgd_cstb_k << ", download+check " << vgd_cstb_d << ", mark exchange "
              << vgd_cstb_m << ", host lists " << vgd_cstb_h << ", upload "
              << vgd_cstb_u << "; mask stride " << vgd_csper << " (padded)"
              << std::endl;
  }
  {
    // ragged band: accesses outside a shell's slab (must be 0)
    Real no = 0.0;
    if (vgd_i.oos.extent(0) > 0) {
      auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_i.oos);
      no = h(0);
    }
    std::cout << "<rad_vet> vet_gd ragged band (rank 0): out-of-slab accesses "
              << no << std::endl;
    if (vgd_ring) {
      std::cout << "<rad_vet> vet_gd_band_ring (rank 0): sweeps " << vgd_rg_nk
                << ", dense band saves " << vgd_rg_ndense << ", save store "
                << vgd_rg_store.extent(0)*sizeof(Real)/1.0e9 << " GB" << std::endl;
    }
    if (vgd_xff) {
      double oos2 = 0.0;
      if (vgd_i_xb.oos.extent(0) > 0) {
        auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_i_xb.oos);
        oos2 = h(0);
      }
      std::cout << "<rad_vet> vet_gd_rotate_xfade_fuse (rank 0): fused builds "
                << vgd_xff_n << ", two-sweep fallbacks " << vgd_xff_nfb
                << ", merged exchanges " << vgd_xff_nmsg << ", out-of-slab accesses "
                << "(other set) " << oos2 << std::endl;
    }
  }
  if (vgd_hcomp > 0) {
    std::cout << "<rad_vet> (rank 0) vet_gd halo lists "
              << vgd_hl_nmade << " made, " << vgd_hl_bytes/1048576.0 << " MB of "
              << vgd_hl_cap/1048576.0 << " (key " << vgd_hl_mb
              << ", GPU: <= free/2); exchanges by list " << vgd_hl_nuse << ", dense "
              << vgd_hl_ndense << std::endl;
  }
  if (vgd_hx_on) {
    std::cout << "<rad_vet> (rank 0) vet_gd_halo_exact lists made "
              << vgd_hx_nmade << ", fallbacks " << vgd_hx_nfall
              << " (marked reads outside the masks' region " << vgd_hx_nout
              << "), read ghosts with an unsent source " << vgd_hx_nun
              << "; received values per list: exact/mask "
              << ((vgd_hx_vmask > 0.0) ? vgd_hx_vexact/vgd_hx_vmask : 0.0)
              << std::endl;
  }
  std::cout << "<rad_vet> (rank 0) exact-list batches " << vgd_hxb_ngroup
            << " (group " << vgd_hxb_g << ", " << vgd_hxb_t << " s, list reallocs "
            << vgd_hxb_nrealloc << "); table prefetch 1"
            << " swapped in " << vgd_tab_nswap << " of the rotations" << std::endl;
}

}  // namespace radvet
