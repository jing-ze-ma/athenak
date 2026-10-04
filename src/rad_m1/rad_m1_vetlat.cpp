//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetlat.cpp
//! \brief <rad_m1>/vet_col_lat = true (m1-vetcol-lat; design rt_design_1003/SP_SC_G0.md
//! sect. 4, plan C1 = variant V5): a LATERAL CORRECTION of closure = vet_col from a
//! few-angle short-characteristics (SC) formal solution on the mesh itself.
//!
//! WHAT IS SOLVED.  At every vet_col build (rad_m1_vetcol.cpp: the same state, source S
//! and extinction chi as the column solution), the grey transfer equation
//!     n . grad I = chi (S - I)
//! is swept shell by shell along x1 for vet_col_lat_nmu Gauss nodes in mu = |n . r_hat|
//! on [0, 1] per hemisphere times vet_col_lat_npsi uniform azimuths psi about r_hat
//! (the LOCAL frame (r_hat, theta_hat, phi_hat) of each cell): inward rays from the top
//! (vacuum) down to the shell vlat_icut, then outward rays back up.  A ray arriving at a
//! cell centre of shell i is traced back to the sphere of the upwind shell (i+1 inward,
//! i-1 outward, or its own shell for a ray that turned between the two); there I, ln chi
//! and ln S are interpolated BILINEARLY in (theta, phi) index space, I also linearly in
//! mu between the bracketing nodes and in psi, the azimuth being rebuilt from the
//! Cartesian direction in each lateral neighbour's own frame; the source is linear in tau
//! along the segment (the positive first-order SC weights of vet_col, VcolW).  This is
//! "naive A" of scsp_0925/sc.py, adapted to the periodic wedge in
//! rt_design_1003/g0_he/sc_he.py, against which the C++ is gated.
//!
//! THE TWIN AND THE CORRECTION.  Every direction is swept a second time through the
//! column's OWN laterally homogeneous field: the same kernel, the same geometry, the same
//! interpolation weights and azimuths, with every lateral read taken at the own column.
//! The moments of the two sweeps give, per cell,
//!     dD_rr = K_rr/J (3-D) - K_rr/J (twin),   D_ra = K_ra/J (3-D) - K_ra/J (twin),
//! i.e. what the lateral structure changes, with the angular and sphericity
//! discretisation errors of the few-angle sweep cancelling to first order (vet_col keeps
//! the exact spherical column part).  All interpolations are written as nested lerps
//! (a + u (b - a)), which return a common value exactly, so on a LATERALLY UNIFORM state
//! the two sweeps agree bit for bit and the correction is exactly 0 (gate b).
//!
//! HAND-OVER (the interface shared with the cubed sphere, CS2): tau_ten slots
//! M1_TT_LAT0..+5 (rr, r-theta, r-phi, then the tangential block, 0 here) in the mesh
//! basis, ghosts filled by the cell-centred exchange.  dD_rr is FOLDED into slot 0
//! (f_K -> clamp(f_K + dD_rr, vet_col_fk_min, 1)): the implicit diagonal, no new term.
//! D_r,lat enters the face equations as the LAGGED Picard term M1SphLat
//! (rad_m1_implicit.hpp), evaluated from the iterate E (vet_col_lat_offdiag, default
//! true).
//!
//! SEVERAL MESHBLOCKS.  Each build sweeps the TWIN first (it reads only its own
//! column) and exchanges its ghost band; then the 3-D sweep runs on every MeshBlock of
//! the pack at once.  An upwind point inside the block reads the current 3-D intensity;
//! one in the lateral ghost band reads the neighbour column's CURRENT twin plus the
//! LAGGED lateral difference (3-D minus twin) of the previous sweep, exchanged after
//! each 3-D sweep.  Only the lateral part of the inflow across block boundaries is
//! lagged (by one build, i.e. one step in science mode), and a laterally uniform state
//! stays exact.  The first build iterates the 3-D sweep vet_col_lat_init_iter times.  A
//! ray whose upwind point lies beyond the ghost band (very oblique rays through
//! laterally thin cells) reads the last ghost (counted, reported at the end).
//!
//! CADENCE.  Built with vet_col (every step; under hesdirk2 + time2_vet_col = predict at
//! U^n + dt K1).  vet_col_lat_every = N > 1 keeps the correction of the last sweep for N
//! steps (the f_rr fold is re-applied after every vet_col build) and is allowed only with
//! time2_vet_col = lag, i.e. in relaxation.
//!
//! LIMITS (fatal at start-up): the spherical-polar wedge (sph_geom), vet_col_axis =
//! radial; the cubed sphere waits for CS2, Cartesian boxes have vet_sc.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

namespace {
void VlatFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_m1>/vet_col_lat: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

// ray types of vlat_geo(l, 2 a + dir, 0)
constexpr int VL_TOP = 0;     // inward, top shell: from the vacuum top face
constexpr int VL_NORM = 1;    // from the upwind shell (i+1 inward, i-1 outward)
constexpr int VL_START = 2;   // outward, first shell: from the inner face of shell icut
constexpr int VL_TURN = 3;    // outward, turned between the shells: own shell, inward I

KOKKOS_INLINE_FUNCTION
Real VlatLerp(const Real a, const Real b, const Real u) {return a + u*(b - a);}

// the positive first-order SC weights (VcolW of rad_m1_vetcol.cpp, term for term)
KOKKOS_INLINE_FUNCTION
void VlatW(const Real dtau, Real &ex, Real &w0, Real &wu) {
  ex = exp(-dtau);
  if (dtau < 1.0e-3) {
    w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
    wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
  } else {
    const Real g = (1.0 - ex)/dtau;
    w0 = 1.0 - g;
    wu = g - ex;
  }
}

// vet_col_source = relaxed: VcolRelaxedSource of rad_m1_vetcol.cpp, term for term
template <class EosT>
KOKKOS_INLINE_FUNCTION
Real VlatRelaxedSource(const EosT &eos, const Real d, const Real tb, const Real eb,
                       const Real rkp, const Real rke, const Real eth, const Real ar,
                       const Real cl, const Real ch, const Real dt) {
  const Real ix = 1.0/(1.0 + ch*dt*rke);
  M1EosDirect<EosT> th{eos};
  Real ee, cv;
  th(d, tb, ee, cv);
  Real ts = tb;
  bool ok = true;
  (void) M1ImplTemperatureT(th, d, tb, ee, cl*dt*rkp*ar*ix, cl*dt*rke*eb*ix, ts, ok);
  if (!(ok && ts > 0.0)) {ts = tb;}
  const Real t2 = ts*ts;
  const Real t4 = t2*t2;
  const Real es = (eb + ch*dt*rkp*ar*t4)*ix;
  return eth*ar*t4 + (1.0 - eth)*es;
}

// mu label of a ray at the upwind shell: position between the descending nodes
void VlatLabel(const std::vector<double> &nd, const double mu, double &a0, double &fa) {
  const int n = static_cast<int>(nd.size());
  double f;
  if (mu >= nd[0]) {
    f = 0.0;
  } else if (mu <= nd[n-1]) {
    f = n - 1;
  } else {
    int q = 0;
    while (q < n-2 && mu < nd[q+1]) {++q;}
    f = q + (nd[q] - mu)/(nd[q] - nd[q+1]);
  }
  a0 = std::min(std::floor(f), static_cast<double>(n - 1));
  fa = f - a0;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetLatInit
//! \brief checks, buffers, exchange objects (called at the end of VetColInit)

void RadiationM1::VetLatInit() {
  Mesh *pm = pmy_pack->pmesh;
  if (pm->use_cubed_sphere) {
    VlatFatal("the cubed sphere waits for stage CS2 (m1-cs-implicit); the tau_ten slots "
              "M1_TT_LAT0.. are the shared interface");
  }
  if (!sph_geom) {
    VlatFatal("needs the spherical-polar wedge (Cartesian: closure = vet_sc)");
  }
  if (vcol_axis_flux) {VlatFatal("needs vet_col_axis = radial");}
  if (vlat_nmu < 1 || vlat_npsi < 2 || vlat_every < 1 || vlat_iinit < 1) {
    VlatFatal("vet_col_lat_nmu >= 1, _npsi >= 2, _every >= 1, _init_iter >= 1");
  }
  if (vlat_every > 1 && t2_vcmode != 0) {
    VlatFatal("vet_col_lat_every > 1 is a RELAXATION option: it needs "
              "time2_vet_col = lag (science mode rebuilds every step)");
  }
  auto &indcs = pm->mb_indcs;
  if (indcs.nx1 < 3) {VlatFatal("needs at least 3 cells along x1");}
  const int nmb = pmy_pack->nmb_thispack;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int c2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int c3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  const int nd = 2*vlat_nmu*vlat_npsi;
  Kokkos::realloc(vlat_i, nmb, nd, c3, c2, c1);
  Kokkos::deep_copy(vlat_i, 0.0);
  Kokkos::realloc(vlat_t, nmb, nd, c3, c2, c1);
  Kokkos::deep_copy(vlat_t, 0.0);
  Kokkos::realloc(vlat_t_c, nmb, nd, 1, 1, 1);
  Kokkos::realloc(vlat_d, nmb, nd, c3, c2, c1);
  Kokkos::deep_copy(vlat_d, 0.0);
  Kokkos::realloc(vlat_d_c, nmb, nd, 1, 1, 1);
  Kokkos::realloc(vlat_cs, nmb, 2, c3, c2, c1);
  Kokkos::deep_copy(vlat_cs, 0.0);
  Kokkos::realloc(vlat_cs_c, nmb, 2, 1, 1, 1);
  const int ntt = static_cast<int>(tau_ten.extent(1));
  Kokkos::realloc(vlat_tt_c, nmb, ntt, 1, 1, 1);
  Kokkos::realloc(vlat_geo, indcs.nx1, 2*vlat_nmu, 4);
  Kokkos::realloc(vlat_mu, vlat_nmu);
  Kokkos::realloc(vlat_w, vlat_nmu);
  Kokkos::realloc(vlat_cnt, 1);
  Kokkos::deep_copy(vlat_cnt, 0.0);
  // Gauss-Legendre nodes on [0, 1] (descending) and weights summing to 1/2 per hemisphere
  {
    std::vector<double> x(vlat_nmu), w(vlat_nmu);
    const int n = vlat_nmu;
    for (int q = 0; q < n; ++q) {   // Newton on P_n
      double z = std::cos(M_PI*(q + 0.75)/(n + 0.5));
      double pp = 1.0;
      for (int it = 0; it < 100; ++it) {
        double p1 = 1.0, p2 = 0.0;
        for (int l = 1; l <= n; ++l) {
          const double p3 = p2;
          p2 = p1;
          p1 = ((2.0*l - 1.0)*z*p2 - (l - 1.0)*p3)/l;
        }
        pp = n*(z*p1 - p2)/(z*z - 1.0);
        const double dz = p1/pp;
        z -= dz;
        if (std::fabs(dz) < 1.0e-15) {break;}
      }
      x[q] = 0.5*(z + 1.0);                 // z descending in q -> x descending
      w[q] = 0.5*2.0/((1.0 - z*z)*pp*pp);   // 0.5 x (the [-1,1] weight)
    }
    vlat_nodes = x;
    auto mu_h = Kokkos::create_mirror_view(vlat_mu);
    auto w_h = Kokkos::create_mirror_view(vlat_w);
    for (int q = 0; q < n; ++q) {
      mu_h(q) = x[q];
      w_h(q) = 0.5*w[q];   // sum over the nodes 1/2: J = sum_hemispheres sum_q w I
    }
    Kokkos::deep_copy(vlat_mu, mu_h);
    Kokkos::deep_copy(vlat_w, w_h);
  }
  pbval_vl = new MeshBoundaryValuesCC(pmy_pack, nullptr, false);
  pbval_vl->InitializeBuffers(nd);
  pbval_vl->SetVectorPairs(nd, {});
  pbval_vs = new MeshBoundaryValuesCC(pmy_pack, nullptr, false);
  pbval_vs->InitializeBuffers(2);
  pbval_vs->SetVectorPairs(2, {});
  pbval_vt = new MeshBoundaryValuesCC(pmy_pack, nullptr, false);
  pbval_vt->InitializeBuffers(ntt);
  pbval_vt->SetVectorPairs(ntt, {});
  vlat_geo_icut = -1;
  vlat_ready = true;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> vet_col_lat: lateral SC correction of vet_col, " << vlat_nmu
              << " mu x " << vlat_npsi << " psi per hemisphere (" << nd
              << " directions), "
              << "3-D minus twin; f_rr folded into the diagonal"
              << (vlat_od ? ", D_r,lat as a lagged Picard term" : ", D_r,lat NOT used")
              << "; sweep from tau_top >= " << vlat_taucut << "; rebuilt every "
              << vlat_every << " step(s); first build " << vlat_iinit
              << " inflow iterations" << std::endl;
    if (vcol_rtop) {
      std::cout << "<rad_m1> vet_col_lat: NOTE the sweep's top is vacuum (vet_col's "
                << "reflecting top is not mirrored; both sweeps alike)" << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetLatExchange
//! \brief the ordinary cell-centred ghost exchange of a (m, n, k, j, i) array

void RadiationM1::VetLatExchange(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                                 MeshBoundaryValuesCC *pb) {
  const int nv = static_cast<int>(a.extent(1));
  while (pb->InitRecv(nv) == TaskStatus::incomplete) {}
  while (pb->PackAndSendCC(a, ac) == TaskStatus::incomplete) {}
  while (pb->RecvAndUnpackCC(a, ac) == TaskStatus::incomplete) {}
  while (pb->ClearSend() == TaskStatus::incomplete) {}
  while (pb->ClearRecv() == TaskStatus::incomplete) {}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetLatBuild
//! \brief after every vet_col build: (re)build the correction when due, fold dD_rr into
//! tau_ten slot 0 and fill the tau_ten ghosts

void RadiationM1::VetLatBuild() {
  if (!vlat_ready) {VetLatInit();}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int cyc = pmy_pack->pmesh->ncycle;
  const bool due = (vlat_nbuild == 0) || (vlat_every <= 1) || ((cyc % vlat_every) == 0);
  if (due) {
    Kokkos::fence();
    Kokkos::Timer timer;
    // the twin first (its ghost band is the current part of the 3-D inflow), then the
    // 3-D sweep(s), each followed by the exchange of the lagged lateral difference
    VetLatSweep(0);
    VetLatExchange(vlat_t, vlat_t_c, pbval_vl);
    const int nit = (vlat_nbuild == 0) ? vlat_iinit : 1;
    for (int it = 0; it < nit; ++it) {
      VetLatSweep(1);
      VetLatExchange(vlat_d, vlat_d_c, pbval_vl);
      vlat_ncall += 1.0;
    }
    VetLatSweep(2);
    vlat_nbuild += 1;
    Kokkos::fence();
    vlat_time += timer.seconds();
  }
  // fold dD_rr into f_K (slot 0, rebuilt by vet_col every step), with vet_col's clamp
  auto tt_ = tau_ten;
  const Real fkm = vcol_fkmin;
  const int ilo = is + vlat_icut;
  par_for("m1_vlat_fold", DevExeSpace(), 0, nmb1, ks, ke, js, je, ilo, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    tt_(m,0,k,j,i) = fmin(fmax(tt_(m,0,k,j,i) + tt_(m,M1_TT_LAT0,k,j,i), fkm), 1.0);
  });
  VetLatExchange(tau_ten, vlat_tt_c, pbval_vt);
  if (due && !vlat_dump.empty() && (vlat_nbuild == 1 || (vlat_dump_every > 0 &&
      (vlat_nbuild % vlat_dump_every) == 0))) {VetLatDump();}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetLatSweep
//! \brief stage 0: source, first shell, geometry and the TWIN sweep; stage 1: the
//! 3-D sweep and its lateral difference; stage 2: the moments -> tau_ten LAT

void RadiationM1::VetLatSweep(const int stage) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, ng = indcs.ng;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool thrd = trans_x3;
  const int nmu = vlat_nmu, npsi = vlat_npsi, nh = vlat_nmu*vlat_npsi;
  auto iw_ = iw;
  auto opac_ = opac;
  auto cs_ = vlat_cs;
  auto vi_ = vlat_i;
  auto vt_ = vlat_t;
  auto tt_ = tau_ten;
  auto geo_ = vlat_geo;
  auto mu_ = vlat_mu;
  auto w_ = vlat_w;
  auto cnt_ = vlat_cnt;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx2v = pmy_pack->pcoord->x2v;
  auto cx3v = pmy_pack->pcoord->x3v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  auto &mbsize = pmy_pack->pmb->mb_size;

  if (stage == 0) {
  // (1) ln chi, ln S of every active cell (vet_col's extinction and source), ghosts
  const bool thermal = fl_on && coupling && !opac_zero;
  const bool srx = thermal && vcol_srelax;
  FluidRef flv = FluidRef::Get(pmy_pack);
  auto eos = flv.eos;
  auto uh = flv.u0;
  const Real chs = chat;
  const Real dts = mr_on ? mr_dt : pm->dt;
  const Real cl = c_light, ar = arad, efl = e_floor;
  par_for("m1_vlat_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real chi = fmax(iw_(m,M1_IW_KT,k,j,i), 1.0e-300);
    Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
    Real s = e;
    if (thermal) {
      Real tg = iw_(m,M1_IW_TP,k,j,i);
      Real t2 = tg*tg;
      Real eth = fmin(opac_(m,M1_OP_P,k,j,i)/chi, 1.0);
      if (srx) {
        s = VlatRelaxedSource(eos, uh(m,IDN,k,j,i), tg, e, opac_(m,M1_OP_P,k,j,i),
                              opac_(m,M1_OP_E,k,j,i), eth, ar, cl, chs, dts);
      } else {
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
    }
    cs_(m,0,k,j,i) = log(chi);
    cs_(m,1,k,j,i) = log(fmax(s, 1.0e-300));
  });
  VetLatExchange(vlat_cs, vlat_cs_c, pbval_vs);

  // (2) the first shell of the sweep: every column has tau_top >= vet_col_lat_taucut
  // there (one global min); the shells below keep dD = 0
  {
    const Real tcut = vlat_taucut;
    const int nk = ke - ks + 1, nj = je - js + 1;
    int lmin = n1 - 2;
    Kokkos::parallel_reduce("m1_vlat_icut",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nk*nj),
      KOKKOS_LAMBDA(const int idx, int &lm) {
        const int m = idx/(nk*nj);
        const int k = ks + (idx/nj) % nk;
        const int j = js + idx % nj;
        Real tau = 0.0;
        int lc = 0;
        for (int i = ie; i >= is; --i) {
          const Real dt_ = exp(cs_(m,0,k,j,i))*(cx1f(m,i+1) - cx1f(m,i));
          if (tau + 0.5*dt_ >= tcut) {lc = i - is; break;}
          tau += dt_;
        }
        lm = (lc < lm) ? lc : lm;
      }, Kokkos::Min<int>(lmin));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &lmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
    vlat_icut = std::max(0, std::min(lmin, n1 - 2));
  }
  }
  const int lcut = vlat_icut;

  // (3) the per-shell ray geometry (host; rebuilt when the first shell moves)
  if (stage == 0 && lcut != vlat_geo_icut) {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cx1v);
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cx1f);
    std::vector<double> rc(n1);
    for (int l = 0; l < n1; ++l) {rc[l] = x1v_h(0, is + l);}
    const double rtop = x1f_h(0, ie + 1);
    const double rin = x1f_h(0, is + lcut);
    auto g_h = Kokkos::create_mirror_view(vlat_geo);
    Kokkos::deep_copy(g_h, 0.0);
    const std::vector<double> &nd = vlat_nodes;
    for (int l = lcut; l < n1; ++l) {
      const double r = rc[l];
      for (int a = 0; a < nmu; ++a) {
        const double mu = nd[a];
        const double p = r*std::sqrt(std::max(1.0 - mu*mu, 0.0));
        const double zr = std::sqrt(std::max(r*r - p*p, 0.0));
        double a0 = 0.0, fa = 0.0;
        // inward (dir 0)
        if (l == n1 - 1) {
          g_h(l, 2*a, 0) = VL_TOP;
          g_h(l, 2*a, 1) = std::sqrt(std::max(rtop*rtop - p*p, 0.0)) - zr;
        } else {
          const double ru = rc[l+1];
          g_h(l, 2*a, 0) = VL_NORM;
          g_h(l, 2*a, 1) = std::sqrt(ru*ru - p*p) - zr;
          VlatLabel(nd, std::sqrt(std::max(1.0 - (p/ru)*(p/ru), 0.0)), a0, fa);
          g_h(l, 2*a, 2) = a0;
          g_h(l, 2*a, 3) = fa;
        }
        // outward (dir 1)
        const double rd = (l > lcut) ? rc[l-1] : rin;
        if (l > lcut && p <= rd*(1.0 + 1.0e-13)) {
          g_h(l, 2*a+1, 0) = VL_NORM;
          g_h(l, 2*a+1, 1) = zr - std::sqrt(std::max(rd*rd - p*p, 0.0));
          VlatLabel(nd, std::sqrt(std::max(1.0 - (p/rd)*(p/rd), 0.0)), a0, fa);
          g_h(l, 2*a+1, 2) = a0;
          g_h(l, 2*a+1, 3) = fa;
        } else if (l == lcut && p <= rin) {
          g_h(l, 2*a+1, 0) = VL_START;
          g_h(l, 2*a+1, 1) = zr - std::sqrt(std::max(rin*rin - p*p, 0.0));
          g_h(l, 2*a+1, 3) = std::sqrt(std::max(1.0 - (p/rin)*(p/rin), 0.0));
        } else {
          g_h(l, 2*a+1, 0) = VL_TURN;
          g_h(l, 2*a+1, 1) = 2.0*zr;
          g_h(l, 2*a+1, 2) = a;
          g_h(l, 2*a+1, 3) = 0.0;
        }
      }
    }
    Kokkos::deep_copy(vlat_geo, g_h);
    vlat_geo_icut = lcut;
  }

  // (4) the sweeps: inward top-down, then outward bottom-up; one launch per shell, one
  // thread per (block, cell, direction).  stage 0: the TWIN (every lateral read at the
  // own column); stage 1: the 3-D sweep, whose lateral reads take the current 3-D
  // intensity inside the block and, in the ghost band, the CURRENT twin of the
  // neighbour column plus the LAGGED lateral difference (3-D minus twin) of the last
  // sweep: a lagged inflow that is exact on a laterally uniform state.
  if (stage != 2) {
  const int jlo = js - ng, jhi = je + ng - 1;   // j0 range with j0+1 in the array
  const int klo = thrd ? (ks - ng) : ks, khi = thrd ? (ke + ng - 1) : ks;
  const Real twopi = 2.0*M_PI;
  const bool tw = (stage == 0);
  auto vd_ = vlat_d;
  for (int pass = 0; pass < 2; ++pass) {
    const int dir = pass;                     // 0 inward, 1 outward
    for (int q = 0; q < n1 - lcut; ++q) {
      const int l = (dir == 0) ? (n1 - 1 - q) : (lcut + q);
      const int i = is + l;
      par_for("m1_vlat_shell", DevExeSpace(), 0, nmb1, ks, ke, js, je, 0, nh - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
        const int a = d/npsi, p = d - a*npsi;
        const int typ = static_cast<int>(geo_(l, 2*a+dir, 0));
        const Real ds = geo_(l, 2*a+dir, 1);
        const Real mu = mu_(a);
        const Real sg = (dir == 0) ? -1.0 : 1.0;
        const Real r = cx1v(m,i);
        const Real th = cx2v(m,j), ph = cx3v(m,k);
        const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
        const Real sm = sqrt(fmax(1.0 - mu*mu, 0.0));
        const Real psi = (p + 0.5)*twopi/npsi;
        const Real c1 = sm*cos(psi), c2 = sm*sin(psi);
        // direction n (Cartesian) = sg mu e_r + c1 e_th + c2 e_ph
        const Real nx = sg*mu*st*cp + c1*ct*cp - c2*sp;
        const Real ny = sg*mu*st*sp + c1*ct*sp + c2*cp;
        const Real nz = sg*mu*ct - c1*st;
        const Real ch0 = exp(cs_(m,0,k,j,i));
        const Real s0 = exp(cs_(m,1,k,j,i));
        const int slot = dir*nh + d;
        Real iv;
        if (typ == VL_TOP) {
          const Real stp = fmax(1.5*exp(cs_(m,1,k,j,i)) - 0.5*exp(cs_(m,1,k,j,i-1)),
                                1.0e-300);
          Real ex, w0, wu;
          VlatW(ch0*ds, ex, w0, wu);
          iv = wu*stp + w0*s0;
        } else {
          // the upwind point, its lateral cell and bilinear fractions (index space)
          const Real xu = r*st*cp - ds*nx;
          const Real yu = r*st*sp - ds*ny;
          const Real zu = r*ct - ds*nz;
          const Real ru = sqrt(xu*xu + yu*yu + zu*zu);
          const Real thu = acos(fmin(fmax(zu/ru, -1.0), 1.0));
          Real dph = atan2(yu, xu) - ph;
          dph -= twopi*floor((dph + M_PI)/twopi);
          const Real fj = j + (thu - th)/mbsize.d_view(m).dx2;
          const Real fk = thrd ? (k + dph/mbsize.d_view(m).dx3) : static_cast<Real>(k);
          int j0 = static_cast<int>(floor(fj));
          int k0 = static_cast<int>(floor(fk));
          Real uj = fj - j0, uk = fk - k0;
          bool clp = false;
          if (j0 < jlo) {j0 = jlo; uj = 0.0; clp = true;}
          if (j0 > jhi) {j0 = jhi; uj = 1.0; clp = true;}
          if (thrd) {
            if (k0 < klo) {k0 = klo; uk = 0.0; clp = true;}
            if (k0 > khi) {k0 = khi; uk = 1.0; clp = true;}
          } else {
            k0 = k;
            uk = 0.0;
          }
          if (clp && !tw) {Kokkos::atomic_add(&cnt_(0), 1.0);}
          const int k1 = thrd ? (k0 + 1) : k0;
          // the lateral cells read: the four neighbours (3-D) or the own column (twin)
          const int ka = tw ? k : k0, kb = tw ? k : k1;
          const int ja = tw ? j : j0, jb = tw ? j : (j0 + 1);
          if (typ == VL_START) {
            // the diffusion intensity at the inner face of shell icut (sc_he.py)
            const Real muf = geo_(l, 2*a+dir, 3);
            const int i1 = i + 1;
            auto ib = [&](const int kk, const int jj) {
              const Real sa = exp(cs_(m,1,kk,jj,i)), sb = exp(cs_(m,1,kk,jj,i1));
              const Real g0 = -(sb - sa)/(cx1v(m,i1) - cx1v(m,i))/exp(cs_(m,0,kk,jj,i));
              return sa + g0*muf;
            };
            auto sbt = [&](const int kk, const int jj) {
              return fmax(1.5*exp(cs_(m,1,kk,jj,i)) - 0.5*exp(cs_(m,1,kk,jj,i1)),
                          1.0e-300);
            };
            Real ex, w0, wu;
            VlatW(ch0*ds, ex, w0, wu);
            const Real bv = fmax(VlatLerp(VlatLerp(ib(ka,ja), ib(ka,jb), uj),
                                          VlatLerp(ib(kb,ja), ib(kb,jb), uj), uk), 0.0);
            const Real qv = VlatLerp(VlatLerp(sbt(ka,ja), sbt(ka,jb), uj),
                                     VlatLerp(sbt(kb,ja), sbt(kb,jb), uj), uk);
            iv = bv*ex + wu*qv + w0*s0;
          } else {
            // upwind shell and the intensity family read there
            const int iu = (typ == VL_TURN) ? i : ((dir == 0) ? (i + 1) : (i - 1));
            const int dsrc = (typ == VL_TURN) ? 0 : dir;
            const int a0 = static_cast<int>(geo_(l, 2*a+dir, 2));
            const Real fa = geo_(l, 2*a+dir, 3);
            const int a1 = (a0 + 1 < nmu) ? (a0 + 1) : a0;
            const int b0 = dsrc*nh;
            // azimuth of n in the frame of lateral cell (kk, jj) (always the 3-D
            // neighbour's frame: the twin shares the arithmetic), lerp in mu and psi
            auto azi = [&](const int kk, const int jj, int &p0, Real &up) {
              const Real tq = cx2v(m,jj), pq = cx3v(m,kk);
              const Real ctq = cos(tq), stq = sin(tq), cpq = cos(pq), spq = sin(pq);
              const Real nt = nx*ctq*cpq + ny*ctq*spq - nz*stq;
              const Real np = -nx*spq + ny*cpq;
              Real ps = atan2(np, nt);
              ps -= twopi*floor(ps/twopi);
              const Real f = ps/(twopi/npsi) - 0.5;
              p0 = static_cast<int>(floor(f));
              up = f - p0;
            };
            // the intensity of slot c at lateral cell (kk, jj) of the upwind shell
            auto rd = [&](const int c, const int kk, const int jj) {
              if (tw) {return vt_(m,c,kk,jj,iu);}
              if (kk >= ks && kk <= ke && jj >= js && jj <= je) {
                return vi_(m,c,kk,jj,iu);
              }
              return vt_(m,c,kk,jj,iu) + vd_(m,c,kk,jj,iu);
            };
            auto val = [&](const int kk, const int jj, const int p0, const Real up) {
              const int pa = ((p0 % npsi) + npsi) % npsi;
              const int pb = (pa + 1) % npsi;
              const Real v0 = VlatLerp(rd(b0+a0*npsi+pa,kk,jj), rd(b0+a0*npsi+pb,kk,jj),
                                       up);
              const Real v1 = VlatLerp(rd(b0+a1*npsi+pa,kk,jj), rd(b0+a1*npsi+pb,kk,jj),
                                       up);
              return VlatLerp(v0, v1, fa);
            };
            int p00, p01, p10, p11;
            Real u00, u01, u10, u11;
            azi(k0, j0, p00, u00);
            azi(k0, j0+1, p01, u01);
            azi(k1, j0, p10, u10);
            azi(k1, j0+1, p11, u11);
            const Real v0 = VlatLerp(val(ka,ja,p00,u00), val(ka,jb,p01,u01), uj);
            const Real v1 = VlatLerp(val(kb,ja,p10,u10), val(kb,jb,p11,u11), uj);
            const Real ivu = VlatLerp(v0, v1, uk);
            auto lin = [&](const int c) {
              return VlatLerp(VlatLerp(cs_(m,c,ka,ja,iu), cs_(m,c,ka,jb,iu), uj),
                              VlatLerp(cs_(m,c,kb,ja,iu), cs_(m,c,kb,jb,iu), uj), uk);
            };
            const Real chu = exp(lin(0));
            const Real su = exp(lin(1));
            Real ex, w0, wu;
            VlatW(0.5*(chu + ch0)*ds, ex, w0, wu);
            iv = fmax(ivu*ex + wu*su + w0*s0, 0.0);
          }
        }
        if (tw) {
          vt_(m,slot,k,j,i) = iv;
        } else {
          vi_(m,slot,k,j,i) = iv;
        }
      });
    }
  }
  if (!tw) {
    // the lateral difference of this sweep: the lagged part of the next sweep's inflow
    const int nd = 2*nh;
    const int ilo = is + lcut;
    par_for("m1_vlat_diff", DevExeSpace(), 0, nmb1, ks, ke, js, je, ilo, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nd; ++c) {vd_(m,c,k,j,i) = vi_(m,c,k,j,i) - vt_(m,c,k,j,i);}
    });
  }
  return;
  }
  // (5) the moments: the correction dD = D(3-D) - D(twin) -> tau_ten LAT slots
  const int ilo = is + lcut;
  const Real twopi = 2.0*M_PI;
  par_for("m1_vlat_mom", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int c = 0; c < M1_TT_NLAT; ++c) {tt_(m,M1_TT_LAT0+c,k,j,i) = 0.0;}
    if (i < ilo) {return;}
    Real j3 = 0.0, krr3 = 0.0, krt3 = 0.0, krp3 = 0.0;
    Real jt = 0.0, krrt = 0.0, krtt = 0.0, krpt = 0.0;
    for (int dir = 0; dir < 2; ++dir) {
      const Real sg = (dir == 0) ? -1.0 : 1.0;
      for (int a = 0; a < nmu; ++a) {
        const Real mu = mu_(a);
        const Real wv = w_(a)/npsi;
        const Real sm = sqrt(fmax(1.0 - mu*mu, 0.0));
        for (int p = 0; p < npsi; ++p) {
          const Real psi = (p + 0.5)*twopi/npsi;
          const Real nr = sg*mu, nt = sm*cos(psi), np = sm*sin(psi);
          const int slot = dir*nh + a*npsi + p;
          const Real v3 = vi_(m,slot,k,j,i), vt = vt_(m,slot,k,j,i);
          j3 += wv*v3;
          krr3 += wv*nr*nr*v3;
          krt3 += wv*nr*nt*v3;
          krp3 += wv*nr*np*v3;
          jt += wv*vt;
          krrt += wv*nr*nr*vt;
          krtt += wv*nr*nt*vt;
          krpt += wv*nr*np*vt;
        }
      }
    }
    if (j3 > 0.0 && jt > 0.0) {
      tt_(m,M1_TT_LAT0,k,j,i) = krr3/j3 - krrt/jt;
      tt_(m,M1_TT_LAT0+1,k,j,i) = krt3/j3 - krtt/jt;
      tt_(m,M1_TT_LAT0+2,k,j,i) = krp3/j3 - krpt/jt;
    }
  });
  auto cnt_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vlat_cnt);
  vlat_nclamp = cnt_h(0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetLatDump
//! \brief vet_col_lat_dump = prefix: one binary file per rank and dumped build,
//! prefix.<build>.<rank>.bin (the first build, and every vet_col_lat_dump_every-th): int64 count, then per active cell 9 doubles
//! (x1v, x2v, x3v, ln chi, ln S, f_K after the fold, dD_rr, D_rt, D_rp)

void RadiationM1::VetLatDump() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  auto tt_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), tau_ten);
  auto cs_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vlat_cs);
  auto x1_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x1v);
  auto x2_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x2v);
  auto x3_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x3v);
  std::vector<double> buf;
  for (int m = 0; m < nmb; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          buf.push_back(x1_h(m,i));
          buf.push_back(x2_h(m,j));
          buf.push_back(x3_h(m,k));
          buf.push_back(cs_h(m,0,k,j,i));
          buf.push_back(cs_h(m,1,k,j,i));
          buf.push_back(tt_h(m,0,k,j,i));
          for (int c = 0; c < 3; ++c) {buf.push_back(tt_h(m,M1_TT_LAT0+c,k,j,i));}
        }
      }
    }
  }
  const std::string fn = vlat_dump + "." + std::to_string(vlat_nbuild) + "."
                         + std::to_string(global_variable::my_rank) + ".bin";
  FILE *fp = std::fopen(fn.c_str(), "wb");
  if (fp == nullptr) {VlatFatal("cannot open " + fn);}
  int64_t n = static_cast<int64_t>(buf.size()/9);
  std::fwrite(&n, sizeof(int64_t), 1, fp);
  std::fwrite(buf.data(), sizeof(double), buf.size(), fp);
  std::fclose(fp);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> vet_col_lat: dumped build " << vlat_nbuild << " to " << vlat_dump
              << "." << vlat_nbuild << ".<rank>.bin (first shell of the sweep i = is + "
              << vlat_icut << ")" << std::endl;
  }
}

} // namespace radm1
