//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ry_per_accretor_raydonor.hpp
//! \brief accretor-raydonor-1017 (USER-approved fix (A), 2026-10-17): RAY-TRACED EXTERNAL
//! DONOR IRRADIATION for ry_per_accretor.cpp, problem/donor_irr_mode = ray (default m1 =
//! the lobe holds E = s a T_loc^4 for the M1 transport, bitwise unchanged).
//! Included INSIDE the anonymous namespace of ry_per_accretor.cpp (it uses RocheParams,
//! RochePotXYZ, RocheGeffXYZ, DonorSurfT and the file's globals).
//!
//! WHY.  vet_col's tensor is per radial column, diag(f_K, (1-f_K)/2, (1-f_K)/2) about r_hat,
//! and the sp curvature carries no r-phi term (M1SphCurv), so the donor's lateral light
//! wrapped round the star as an azimuthal beam (night/facing dE ~1, |dF_phi|/(c dE) up to
//! 1.9; accbig_1016 NOTE "Fix 3").  Here the donor's light never enters the M1 field: it is
//! a beam from the donor's surface, attenuated along straight rays, deposited where it is
//! removed; only the light the gas re-emits (and the scattered part) is M1 radiation.
//!
//! GEOMETRY (2-D band, the production geometry).  The donor is the 3-D Roche L1 lobe
//! (exact surface Phi = Phi_L1, x > d_L1), each surface point a grey isotropic emitter with
//! I = (c/4 pi) a T_loc^4, T_loc^4 = max(C g^(4 beta), floor) (the donor_gdark law, same C
//! and floor as donor_irr = m1) x problem/donor_ray_scale (default 1: NO dilution or
//! calibration factor; the 3-D solid angle is integrated directly).  For every (r, phi)
//! column P of the band the incident field is
//!     q = int dOmega I(n) exp(-tau(P, n)) x (absorption in the cell),
//! over the 3-D directions n toward the donor.  Directions: azimuth psi about P (npsi
//! midpoints over the lobe's in-plane angular span seen from P, computed from a 720-vertex
//! polygon of the in-plane lobe cut) x elevation e (nel Gauss-Legendre nodes on
//! [0, e_max(psi)], doubled for -e by symmetry), where e_max is where the 3-D ray grazes the
//! lobe (from the lobe's surface height z_s(x, y) sampled at 32 points along the in-plane
//! chord); the hit point of each (psi, e) ray (first s with z_s(s) >= s tan e) gives g and
//! so I.  This geometry is static: it is tabulated once at setup per column (float W =
//! I cos e w_e dpsi, cos e).
//! The band is z-uniform (reflecting theta walls), so a ray at elevation e crosses the same
//! gas as its in-plane projection, with 1/cos e the path: tau_3D = tau_proj / cos e.  Near
//! the photosphere this equals the plane-parallel slab of a sphere (tangent plane of the
//! band's cylinder = tangent plane of the sphere at the equator), so the horizon and the
//! limb are those of the 3-D star.  Beyond the in-plane lobe entry the projected ray is in
//! the masked lobe (floor gas, tau ~ 0), beyond r_out vacuum, below r_in opaque.
//!
//! ATTENUATION.  Long characteristics: one march per (column, psi) from the cell's
//! donor-side face along the in-plane direction through a GLOBAL extinction map chi(k, i)
//! (the theta-mean of max(rho kappa_T, rho kappa_E) of <rad_m1>'s opac, i.e. the module's
//! own Rosseland total incl. electron scattering, gathered over ranks every update), steps
//! h = donor_ray_step x min(dr, r dphi) (default 0.5), stopped at tau > donor_ray_taumax
//! (40), at the in-plane lobe entry, at r_out (vacuum) or r_in (opaque).
//! DEPOSITION (photon conserving): per ray, I dOmega exp(-tau_up/cos e) (1 - exp(-chi L/cos e))
//! / (L/cos e), L = the in-plane chord of the cell along psi (1/max(|n_r|/dr, |n_phi|/(r dphi))).
//! problem/donor_ray_deposit = rad (DEFAULT): the whole removed power q is added to the M1 E
//! at the last hydro stage (isotropic, F unchanged); the implicit gas-radiation coupling then
//! hands the absorbed part to the gas at the iterate's kappa(T).  = gas: eps = rho kappa_E/chi
//! of q heats u(IEN) explicitly every stage (RK weighted), the rest goes to E.  The explicit
//! gas split is UNSTABLE in the thin atmosphere (10-17 smoke: kappa_P ~ 2000 cm^2/g at
//! rho 1e-13 g/cc took e/rho 2e3 -> 8e5 in one step, the hot gas then cannot cool through
//! its now tiny kappa_P; dt 7e-6 -> 6e-7).  With rad the beam's local energy density
//! (|F_inc|/c ~ 0.1-0.3) is < 1 % of the gainer's own E near the photosphere, so the thin gas
//! T barely depends on where the deposit goes; in thick cells the deposit is absorbed in place.
//! Momentum: the removed beam momentum q/c along the propagation direction x rho kappa_T/chi
//! goes to the gas (problem/donor_ray_force, default true), with its work v.f in u(IEN): the
//! flux-mean (transport) opacity, as the M1 module couples the gainer's own flux.  With the
//! Planck mean (smoke 10-17, total q/c) the 1e-7 thin atmosphere got 28x gravity inward
//! (v_r -290 km/s in 60 cycles).
//! EXTINCTION: problem/donor_ray_opac = total (default; max(rho kappa_T, rho kappa_E), the
//! M1 module's own absorption + Rosseland) | rosseland (rho kappa_T: Rosseland incl. electron
//! scattering, a continuum beam).  NOTE: with total, the 1e-8 floor gas (20-30 kK, Planck mean
//! 300-1600 cm^2/g of the table's low-rho extension) has tau ~ 1 between the stars.
//! The march step grows x2 per step up to donor_ray_grow (4) cell sizes while chi h < 1e-3.
//! The incident flux at the cell centre (exp(-chi L/2 cos e) inside the own cell) is kept as
//! a diagnostic (dump).
//! LEDGER: everything the beam gives the domain (gas heat + force work + scattered E) is
//! booked as radiation IN: Eirr (hst column 24, as donor_irr = m1) and Ebnd (register 10).
//! The donor lobe (and the stream_sector_light = surface sector) no longer hold light: their
//! radiation is left to the transport (as donor_irr = false); the lobe is then a transparent
//! sink for the gainer's light (no reflection effect).
//! UPDATE: every problem/donor_ray_every cycles (default 1), at the first hydro stage, from
//! the opacities of the last M1 solve (first update of a run: Opacity() of the current
//! state).  Between updates q, f are held per volume.
//! 3-D (wide theta, g3d_): NOT implemented (FATAL): rays must then be marched in 3-D per cell
//! and direction (no z-uniformity); see accretor_rhd_design.md 17.
//! DUMP (problem/donor_ray_dump_dt > 0): at the first update after each multiple, rank 0
//! writes ryper.raydonor.NNNNN.bin in the run directory: int32 n3g, n1g, nf (= 6), float64
//! time, float64 x1v[n1g], float64 x3v[n3g], float32 data[nf][n3g][n1g] with fields
//! q (removed power per volume), f_r, f_phi (force per volume), F_r, F_phi (incident flux
//! at the centre, vector, in the direction of propagation), chi (extinction map).

constexpr int kRdNelMax = 8;
constexpr int kRdNs = 32;          // lobe-height samples along the in-plane chord
constexpr int kRdNpoly = 720;      // in-plane lobe cut polygon
constexpr int kRdNf = 6;           // dump fields

struct RdPar {
  RocheParams p;
  Real phl1, xl1, rb;              // Phi_L1, d_L1, bounding radius of the lobe about the donor
  Real ec, ef, eun, gdb;           // E = max(ec g^gdb, ef) (gdark) or eun; x c/(4 pi) = I
  Real cfac;                       // c/(4 pi)
  int gdk;
  int nel;
  Real glx[kRdNelMax], glw[kRdNelMax];   // Gauss-Legendre on [0, 1]
};

bool rd_on_ = false;
RdPar rd_par_;
int rd_npsi_ = 24, rd_every_ = 1;
Real rd_fh_ = 0.5, rd_taumax_ = 40.0, rd_scale_ = 1.0;
bool rd_force_ = true;
bool rd_dgas_ = false;             // donor_ray_deposit = gas (eps q into u(IEN)); default rad
bool rd_oros_ = false;             // donor_ray_opac = rosseland (rho kappa_T only); default total
int rd_grow_ = 4;                  // donor_ray_grow: max step growth in thin gas
int rd_nupd_ = 0;
bool rd_have_ = false;
Real rd_dump_dt_ = 0.0, rd_dump_next_ = 0.0;
int rd_dump_n_ = 0;
int rd_n1g_ = 0, rd_n3g_ = 0, rd_n2g_ = 1;
Real rd_x3min_ = 0.0, rd_x3max_ = 0.0, rd_dph_ = 0.0, rd_rin_ = 0.0, rd_rout_ = 0.0;
DvceArray1D<Real> rd_rf_;          // global x1 faces (n1g + 1)
DvceArray2D<float> rd_chi_;        // global extinction map (n3g, n1g), theta mean
DvceArray1D<int> rd_koff_, rd_ioff_;   // per local MeshBlock: global k, i offsets
DvceArray2D<Real> rd_span_;        // (ncol, 2): psi_lo, dpsi (dpsi = 0: no rays)
DvceArray4D<float> rd_geo_;        // (ncol, npsi, nel, 2): W = I cos e w dpsi, cos e
DvceArray2D<Real> rd_q_;           // (ncol, 5): q, f_r, f_phi, F_r, F_phi
DvceArray2D<Real> rd_poly_;        // (kRdNpoly + 1, 2): in-plane lobe cut (x, y)
std::vector<Real> rd_x1v_h_, rd_x3v_h_;
Real rd_qtot_ = 0.0;               // last update: total removed power (band, code), diag

KOKKOS_INLINE_FUNCTION
bool RdInLobe2(const RdPar &q, const Real x, const Real y) {
  return (x > q.xl1) && (RochePotXYZ(q.p, x, y, 0.0) < q.phl1);
}

//! the lobe surface height above (x, y) inside the in-plane cut (Phi grows with |z|)
KOKKOS_INLINE_FUNCTION
Real RdZsurf(const RdPar &q, const Real x, const Real y) {
  if (!RdInLobe2(q, x, y)) return 0.0;
  Real lo = 0.0, hi = q.rb;
  for (int n=0; n<34; ++n) {
    const Real zm = 0.5*(lo + hi);
    if (RochePotXYZ(q.p, x, y, zm) < q.phl1) lo = zm; else hi = zm;
  }
  return 0.5*(lo + hi);
}

KOKKOS_INLINE_FUNCTION
Real RdIntens(const RdPar &q, const Real x, const Real y, const Real z) {
  Real e = q.eun;
  if (q.gdk) {
    const Real g = RocheGeffXYZ(q.p, x, y, z);
    e = fmax(q.ec*pow(g, q.gdb), q.ef);
  }
  return q.cfac*e;
}

//! the in-plane azimuth span of the lobe seen from (x0, y0), relative to the direction to
//! the donor centre; returns false (no rays) for a point inside the lobe cut
KOKKOS_INLINE_FUNCTION
bool RdSpan(const RdPar &q, const Real *poly, const int npoly, const Real x0, const Real y0,
            Real &plo, Real &phi) {
  if (RdInLobe2(q, x0, y0)) return false;
  const Real dx = q.p.asep - x0, dy = -y0, dd = sqrt(dx*dx + dy*dy);
  const Real cx = dx/dd, cy = dy/dd;
  plo = 1.0e30; phi = -1.0e30;
  for (int n=0; n<npoly; ++n) {
    const Real vx = poly[2*n] - x0, vy = poly[2*n+1] - y0;
    const Real ps = atan2(cx*vy - cy*vx, cx*vx + cy*vy);
    plo = fmin(plo, ps); phi = fmax(phi, ps);
  }
  return (phi > plo) && (phi - plo < M_PI);
}

//! the elevation quadrature of one in-plane direction (nx, ny) from (x0, y0): W[b] = I(hit)
//! cos e_b w_b dpsi and cos e_b for the nel GL nodes on [0, e_max]; W = 0 if the line misses
//! the lobe
KOKKOS_INLINE_FUNCTION
void RdElev(const RdPar &q, const Real x0, const Real y0, const Real nx, const Real ny,
            const Real dpsi, float *W, float *ce) {
  for (int b=0; b<q.nel; ++b) {W[b] = 0.0f; ce[b] = 1.0f;}
  const Real sc = (q.p.asep - x0)*nx - y0*ny;
  const Real slo = fmax(sc - q.rb, 0.0), shi = sc + q.rb;
  if (shi <= slo) return;
  const int nm = 400;
  const Real hs = (shi - slo)/nm;
  Real sen = -1.0, sex = -1.0;
  bool prev = RdInLobe2(q, x0 + slo*nx, y0 + slo*ny);
  if (prev) return;
  for (int n=1; n<=nm; ++n) {
    const Real s = slo + n*hs;
    const bool in = RdInLobe2(q, x0 + s*nx, y0 + s*ny);
    if (in != prev) {
      Real a = s - hs, c = s;
      for (int it=0; it<40; ++it) {
        const Real m = 0.5*(a + c);
        if (RdInLobe2(q, x0 + m*nx, y0 + m*ny) == prev) a = m; else c = m;
      }
      if (in) {
        sen = 0.5*(a + c);
      } else {
        sex = 0.5*(a + c);
        break;
      }
    }
    prev = in;
  }
  if (sen < 0.0) return;
  if (sex < 0.0) sex = shi;
  Real ss[kRdNs + 2], zz[kRdNs + 2];
  ss[0] = sen; zz[0] = 0.0; ss[kRdNs + 1] = sex; zz[kRdNs + 1] = 0.0;
  Real tmx = 0.0;
  for (int j=1; j<=kRdNs; ++j) {
    const Real s = sen + (j - 0.5)*(sex - sen)/kRdNs;
    ss[j] = s;
    zz[j] = RdZsurf(q, x0 + s*nx, y0 + s*ny);
    tmx = fmax(tmx, zz[j]/s);
  }
  if (!(tmx > 0.0)) return;
  const Real emax = atan(tmx);
  for (int b=0; b<q.nel; ++b) {
    const Real e = emax*q.glx[b], w = emax*q.glw[b], t = tan(e);
    int jh = -1;
    for (int j=1; j<=kRdNs + 1; ++j) {
      if (zz[j] - t*ss[j] >= 0.0) {jh = j; break;}
    }
    Real sh;
    if (jh < 0) {
      // e within roundoff of e_max: the sample of the maximum
      jh = 1;
      for (int j=2; j<=kRdNs; ++j) if (zz[j]/ss[j] > zz[jh]/ss[jh]) jh = j;
      sh = ss[jh];
    } else {
      const Real g0 = zz[jh-1] - t*ss[jh-1], g1 = zz[jh] - t*ss[jh];
      const Real f = (g0 < 0.0 && g1 > g0) ? g0/(g0 - g1) : 0.0;
      sh = ss[jh-1] + f*(ss[jh] - ss[jh-1]);
    }
    const Real ib = RdIntens(q, x0 + sh*nx, y0 + sh*ny, t*sh);
    W[b] = static_cast<float>(ib*cos(e)*w*dpsi);
    ce[b] = static_cast<float>(cos(e));
  }
}

//! unattenuated incident flux (vector, propagation direction) at (x0, y0), host or device,
//! with npsi azimuths (no tabulation): the setup self-check
KOKKOS_INLINE_FUNCTION
void RdFluxFree(const RdPar &q, const Real *poly, const int npoly, const Real x0,
                const Real y0, const int npsi, Real &fx, Real &fy) {
  fx = 0.0; fy = 0.0;
  Real plo, phi;
  if (!RdSpan(q, poly, npoly, x0, y0, plo, phi)) return;
  const Real dpsi = (phi - plo)/npsi;
  const Real b0 = atan2(-y0, q.p.asep - x0);
  float W[kRdNelMax], ce[kRdNelMax];
  for (int a=0; a<npsi; ++a) {
    const Real be = b0 + plo + (a + 0.5)*dpsi, nx = cos(be), ny = sin(be);
    RdElev(q, x0, y0, nx, ny, dpsi, W, ce);
    for (int b=0; b<q.nel; ++b) {
      fx -= 2.0*W[b]*ce[b]*nx;
      fy -= 2.0*W[b]*ce[b]*ny;
    }
  }
}

void RdGaussLegendre01(const int n, Real *x, Real *w) {
  for (int i=0; i<n; ++i) {
    Real z = std::cos(M_PI*(i + 0.75)/(n + 0.5)), pp = 0.0;
    for (int it=0; it<100; ++it) {
      Real p1 = 1.0, p2 = 0.0;
      for (int j=0; j<n; ++j) {
        const Real p3 = p2; p2 = p1;
        p1 = ((2.0*j + 1.0)*z*p2 - j*p3)/(j + 1.0);
      }
      pp = n*(z*p1 - p2)/(z*z - 1.0);
      const Real dz = p1/pp;
      z -= dz;
      if (std::fabs(dz) < 1.0e-15) break;
    }
    x[i] = 0.5*(1.0 - z);
    w[i] = 1.0/((1.0 - z*z)*pp*pp);   // = 0.5 x the [-1, 1] weight 2/((1-z^2) p'^2)
  }
}

//----------------------------------------------------------------------------------------
//! setup (after the donor_irr constants): the global grid, the lobe polygon, the static
//! per-column geometry table, the self-check
void RayDonorSetup(MeshBlockPack *pmbp, Mesh *pm) {
  const bool root = (global_variable::my_rank == 0);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ks = indcs.ks;
  const int ni = indcs.nx1, nk = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  rd_n1g_ = pm->mesh_indcs.nx1;
  rd_n2g_ = pm->mesh_indcs.nx2;
  rd_n3g_ = pm->mesh_indcs.nx3;
  rd_x3min_ = pm->mesh_size.x3min;
  rd_x3max_ = pm->mesh_size.x3max;
  rd_dph_ = (rd_x3max_ - rd_x3min_)/rd_n3g_;
  rd_rin_ = pm->mesh_size.x1min;
  rd_rout_ = pm->mesh_size.x1max;
  // global offsets of the local MeshBlocks (uniform root grid)
  Kokkos::realloc(rd_koff_, nmb);
  Kokkos::realloc(rd_ioff_, nmb);
  auto koh = Kokkos::create_mirror_view(rd_koff_);
  auto ioh = Kokkos::create_mirror_view(rd_ioff_);
  auto &gid = pmbp->pmb->mb_gid;
  for (int m=0; m<nmb; ++m) {
    const LogicalLocation &ll = pm->lloc_eachmb[gid.h_view(m)];
    koh(m) = static_cast<int>(ll.lx3)*nk;
    ioh(m) = static_cast<int>(ll.lx1)*ni;
    if (ll.level != pm->root_level) {
      std::cout << "### FATAL ERROR ry_per_accretor: donor_irr_mode = ray needs a uniform "
                << "root grid (no mesh refinement)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  Kokkos::deep_copy(rd_koff_, koh);
  Kokkos::deep_copy(rd_ioff_, ioh);
  // global x1 faces and centres, x3 centres (MAX over ranks: each entry from its owner)
  std::vector<Real> rf(rd_n1g_ + 1, -1.0e300), xv(rd_n1g_, -1.0e300), pv(rd_n3g_, -1.0e300);
  {
    auto x1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
    auto x1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
    auto x3v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);
    for (int m=0; m<nmb; ++m) {
      for (int i=0; i<=ni; ++i) rf[ioh(m) + i] = x1f(m, is + i);
      for (int i=0; i<ni; ++i) xv[ioh(m) + i] = x1v(m, is + i);
      for (int k=0; k<nk; ++k) pv[koh(m) + k] = x3v(m, ks + k);
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, rf.data(), rd_n1g_ + 1, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, xv.data(), rd_n1g_, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, pv.data(), rd_n3g_, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  rd_x1v_h_ = xv;
  rd_x3v_h_ = pv;
  Kokkos::realloc(rd_rf_, rd_n1g_ + 1);
  {
    auto h = Kokkos::create_mirror_view(rd_rf_);
    for (int i=0; i<=rd_n1g_; ++i) h(i) = rf[i];
    Kokkos::deep_copy(rd_rf_, h);
  }
  Kokkos::realloc(rd_chi_, rd_n3g_, rd_n1g_);
  // the in-plane lobe cut: rays from the donor centre (DonorSurfT), plus the L1 point
  const RdPar q = rd_par_;
  std::vector<Real> poly(2*(kRdNpoly + 1));
  int np = 0;
  for (int n=0; n<kRdNpoly; ++n) {
    const Real al = (n + 0.5)*2.0*M_PI/kRdNpoly;
    const Real nx = -std::cos(al), ny = std::sin(al);
    const Real t = DonorSurfT(q.p, q.phl1, q.xl1, nx, ny, 0.0, 0.0);
    if (t <= 0.0) continue;
    poly[2*np] = q.p.asep + t*nx; poly[2*np+1] = t*ny; ++np;
  }
  poly[2*np] = q.xl1; poly[2*np+1] = 0.0; ++np;
  const int npoly = np;
  Kokkos::realloc(rd_poly_, npoly, 2);
  {
    auto h = Kokkos::create_mirror_view(rd_poly_);
    for (int n=0; n<npoly; ++n) {h(n, 0) = poly[2*n]; h(n, 1) = poly[2*n+1];}
    Kokkos::deep_copy(rd_poly_, h);
  }
  // per-column geometry table (static)
  const int ncol = nmb*nk*ni, npsi = rd_npsi_, nel = q.nel;
  Kokkos::realloc(rd_span_, ncol, 2);
  Kokkos::realloc(rd_geo_, ncol, npsi, nel, 2);
  Kokkos::realloc(rd_q_, ncol, 5);
  Kokkos::deep_copy(rd_q_, 0.0);
  {
    auto span = rd_span_;
    auto geo = rd_geo_;
    auto pl = rd_poly_;
    auto &x1v = pmbp->pcoord->x1v;
    auto &x3v = pmbp->pcoord->x3v;
    const int nki = nk*ni;
    Kokkos::parallel_for("ryper_rd_geo", Kokkos::RangePolicy<>(DevExeSpace(), 0, ncol),
    KOKKOS_LAMBDA(const int c) {
      const int m = c/nki, k = (c - m*nki)/ni + ks, i = c%ni + is;
      const Real r = x1v(m,i), ph = x3v(m,k);
      const Real x0 = r*cos(ph), y0 = r*sin(ph);
      Real plo, phi;
      span(c, 0) = 0.0; span(c, 1) = 0.0;
      for (int a=0; a<npsi; ++a) for (int b=0; b<nel; ++b) {
        geo(c, a, b, 0) = 0.0f; geo(c, a, b, 1) = 1.0f;
      }
      if (!RdSpan(q, pl.data(), npoly, x0, y0, plo, phi)) return;
      const Real dpsi = (phi - plo)/npsi;
      span(c, 0) = plo; span(c, 1) = dpsi;
      const Real b0 = atan2(-y0, q.p.asep - x0);
      float W[kRdNelMax], ce[kRdNelMax];
      for (int a=0; a<npsi; ++a) {
        const Real be = b0 + plo + (a + 0.5)*dpsi;
        RdElev(q, x0, y0, cos(be), sin(be), dpsi, W, ce);
        for (int b=0; b<nel; ++b) {geo(c, a, b, 0) = W[b]; geo(c, a, b, 1) = ce[b];}
      }
    });
    Kokkos::fence();
  }
  if (root) {
    std::printf("ry_per_accretor: donor_irr_mode = ray: npsi %d, nel %d (GL on [0, e_max]), "
                "every %d cycle(s), step %.3g x min(dr, r dphi), tau_max %.3g, force %s, "
                "scale %.4g; lobe polygon %d vertices, bounding radius %.4f; table %.1f MB "
                "per rank\n", npsi, nel, rd_every_, rd_fh_, rd_taumax_,
                rd_force_ ? "on" : "off", rd_scale_, npoly, q.rb,
                (4.0*ncol*npsi*nel*2 + 8.0*ncol*7)/1.0e6);
    // self-check: the unattenuated incident flux on the line of centres and at the
    // quadrature, vs a 16x finer one and the sphere estimate L_d/(4 pi D^2)
    RdPar qf = q;
    qf.nel = kRdNelMax;
    RdGaussLegendre01(kRdNelMax, qf.glx, qf.glw);
    const Real rr[6] = {9.4, 9.6, 10.0, 10.4, 12.0, 14.0};
    for (int n=0; n<6; ++n) {
      Real fx, fy, gx, gy;
      RdFluxFree(q, poly.data(), npoly, rr[n], 0.0, npsi, fx, fy);
      RdFluxFree(qf, poly.data(), npoly, rr[n], 0.0, 16*npsi, gx, gy);
      std::printf("ry_per_accretor: raydonor self-check (no gas): r %.2f phi 0: F_inc,r = "
                  "%.6g (table quadrature) / %.6g (16 npsi x %d nel) code\n", rr[n], -fx,
                  -gx, kRdNelMax);
    }
    const Real ph2[3] = {0.5*M_PI, 0.75*M_PI, M_PI};
    for (int n=0; n<3; ++n) {
      Real fx, fy;
      RdFluxFree(q, poly.data(), npoly, 10.0*std::cos(ph2[n]), 10.0*std::sin(ph2[n]), npsi,
                 fx, fy);
      std::printf("ry_per_accretor: raydonor self-check (no gas): r 10 phi %.0f deg: |F_inc| "
                  "%.6g code (geometric only: the star's shadow comes from tau)\n",
                  ph2[n]*180.0/M_PI, std::sqrt(fx*fx + fy*fy));
    }
  }
}

//----------------------------------------------------------------------------------------
//! the update: global extinction map, then one march per (column, psi)
void RayDonorUpdate(MeshBlockPack *pmbp, Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int ni = indcs.nx1, nj = indcs.nx2, nk = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto *pm1 = pmbp->pradm1;
  if (!rd_have_ && !pm1->opac_freeze) (void) pm1->Opacity(nullptr, 1);
  auto opac = pm1->opac;
  auto chi = rd_chi_;
  auto koff = rd_koff_, ioff = rd_ioff_;
  const int n1g = rd_n1g_, n3g = rd_n3g_;
  const Real rn2 = 1.0/nj;
  const bool oros = rd_oros_;
  const auto tw0 = std::chrono::steady_clock::now();
  Kokkos::deep_copy(chi, 0.0f);
  const int nki = nk*ni, ncol = nmb*nki;
  Kokkos::parallel_for("ryper_rd_map", Kokkos::RangePolicy<>(DevExeSpace(), 0, ncol),
  KOKKOS_LAMBDA(const int c) {
    const int m = c/nki, k = (c - m*nki)/ni + ks, i = c%ni + is;
    Real s = 0.0;
    for (int j=js; j<js+nj; ++j) {
      s += oros ? opac(m,radm1::M1_OP_T,k,j,i)
                : fmax(opac(m,radm1::M1_OP_T,k,j,i), opac(m,radm1::M1_OP_E,k,j,i));
    }
    chi(koff(m) + k - ks, ioff(m) + i - is) = static_cast<float>(s*rn2);
  });
#if MPI_PARALLEL_ENABLED
  {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), chi);
    MPI_Allreduce(MPI_IN_PLACE, h.data(), n3g*n1g, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    Kokkos::deep_copy(chi, h);
  }
#endif
  const RdPar q = rd_par_;
  auto span = rd_span_;
  auto geo = rd_geo_;
  auto qq = rd_q_;
  auto rf = rd_rf_;
  auto &x1v = pmbp->pcoord->x1v;
  auto &x3v = pmbp->pcoord->x3v;
  const int npsi = rd_npsi_, nel = q.nel;
  const Real fh = rd_fh_, tmx = rd_taumax_, rin = rd_rin_, rout = rd_rout_;
  const Real x3a = rd_x3min_, x3b = rd_x3max_, dph = rd_dph_;
  const bool per = std::fabs(rd_x3max_ - rd_x3min_ - 2.0*M_PI) < 1.0e-10;
  const Real cl = pm1->c_light;
  const bool frc = rd_force_;
  const Real gmx = static_cast<Real>(rd_grow_);
  Real qsum = 0.0;
  auto &vol = pmbp->pcoord->volume;
  Kokkos::parallel_reduce("ryper_rd_march", Kokkos::RangePolicy<>(DevExeSpace(), 0, ncol),
  KOKKOS_LAMBDA(const int c, Real &qs) {
    const int m = c/nki, k = (c - m*nki)/ni + ks, i = c%ni + is;
    for (int n=0; n<5; ++n) qq(c, n) = 0.0;
    const Real dpsi = span(c, 1);
    if (!(dpsi > 0.0)) return;
    const Real r0 = x1v(m,i), ph0 = x3v(m,k);
    const Real cp0 = cos(ph0), sp0 = sin(ph0);
    const Real x0 = r0*cp0, y0 = r0*sp0;
    const int ig0 = ioff(m) + i - is, kg0 = koff(m) + k - ks;
    const Real dr0 = rf(ig0 + 1) - rf(ig0), chic = chi(kg0, ig0);
    const Real b0 = atan2(-y0, q.p.asep - x0), plo = span(c, 0);
    Real qt = 0.0, fx = 0.0, fy = 0.0, gx = 0.0, gy = 0.0;
    for (int a=0; a<npsi; ++a) {
      const Real be = b0 + plo + (a + 0.5)*dpsi, nx = cos(be), ny = sin(be);
      const Real nr = nx*cp0 + ny*sp0, nph = -nx*sp0 + ny*cp0;
      const Real lp = 1.0/fmax(fabs(nr)/dr0, fabs(nph)/(r0*dph));
      // march from the donor-side face of the own cell
      Real s = 0.5*lp, tau = 0.0, gf = 1.0;
      int ig = ig0;
      for (int it=0; it<100000; ++it) {
        const Real x = x0 + s*nx, y = y0 + s*ny;
        const Real r = sqrt(x*x + y*y);
        if (r <= rin) {tau = 1.0e30; break;}
        if (r >= rout) break;
        if (x > q.xl1 && RochePotXYZ(q.p, x, y, 0.0) < q.phl1) break;
        while (ig > 0 && r < rf(ig)) --ig;
        while (ig < n1g - 1 && r >= rf(ig + 1)) ++ig;
        Real ph = atan2(y, x);
        if (ph < x3a) ph += 2.0*M_PI;
        // the step grows (x2 per step, up to donor_ray_grow) while it stays optically thin
        // (chi h < 1e-3) and falls back to one cell size as soon as it is not
        const Real h0 = fh*fmin(rf(ig + 1) - rf(ig), r*dph);
        Real h = gf*h0, dt = 0.0;
        if (ph < x3b) {
          int kg = static_cast<int>((ph - x3a)/dph);
          kg = (kg < 0) ? 0 : ((kg > n3g - 1) ? n3g - 1 : kg);
          dt = chi(kg, ig)*h;
          if (dt > 1.0e-3 && gf > 1.0) {h = h0; dt = chi(kg, ig)*h;}
        } else if (per) {
          dt = 0.0;
        }
        tau += dt;
        gf = (dt < 1.0e-3) ? fmin(2.0*gf, gmx) : 1.0;
        s += h;
        if (tau > tmx) break;
      }
      if (tau > tmx) continue;
      for (int b=0; b<nel; ++b) {
        const Real w = geo(c, a, b, 0), ce = geo(c, a, b, 1);
        if (!(w > 0.0)) continue;
        const Real t3 = tau/ce, d3 = chic*lp/ce, l3 = lp/ce;
        const Real att = exp(-t3);
        const Real dep = (d3 > 1.0e-6) ? (1.0 - exp(-d3))/l3 : chic*(1.0 - 0.5*d3);
        const Real qa = 2.0*w*att*dep;
        qt += qa;
        fx -= qa*ce*nx;
        fy -= qa*ce*ny;
        const Real fa = 2.0*w*att*exp(-0.5*d3)*ce;
        gx -= fa*nx;
        gy -= fa*ny;
      }
    }
    qq(c, 0) = qt;
    if (frc) {
      qq(c, 1) = (fx*cp0 + fy*sp0)/cl;
      qq(c, 2) = (-fx*sp0 + fy*cp0)/cl;
    }
    qq(c, 3) = gx*cp0 + gy*sp0;
    qq(c, 4) = -gx*sp0 + gy*cp0;
    for (int j=js; j<js+nj; ++j) qs += qt*vol(m,k,j,i);
  }, Kokkos::Sum<Real>(qsum));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &qsum, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  rd_qtot_ = qsum;
  rd_have_ = true;
  rd_nupd_++;
  if (rd_nupd_ <= 3 && global_variable::my_rank == 0) {
    const double tw = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                                    - tw0).count();
    std::printf("ry_per_accretor: raydonor update %d (cycle %d): %.3f s wall (rank 0), "
                "removed power %.6g code\n", rd_nupd_, pm->ncycle, tw, qsum);
  }
  // the dump
  if (rd_dump_dt_ > 0.0 && pm->time >= rd_dump_next_) {
    const size_t nn = static_cast<size_t>(n3g)*n1g;
    DvceArray3D<float> gd("ryper_rd_dump", kRdNf, n3g, n1g);
    Kokkos::parallel_for("ryper_rd_dumpf", Kokkos::RangePolicy<>(DevExeSpace(), 0, ncol),
    KOKKOS_LAMBDA(const int c) {
      const int m = c/nki, k = (c - m*nki)/ni + ks, i = c%ni + is;
      const int kg = koff(m) + k - ks, ig = ioff(m) + i - is;
      for (int n=0; n<5; ++n) gd(n, kg, ig) = static_cast<float>(qq(c, n));
    });
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), gd);
#if MPI_PARALLEL_ENABLED
    if (global_variable::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, h.data(), 5*nn, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(h.data(), nullptr, 5*nn, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
    if (global_variable::my_rank == 0) {
      auto hc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), chi);
      for (size_t n=0; n<nn; ++n) h.data()[5*nn + n] = hc.data()[n];
      char fn[64];
      std::snprintf(fn, sizeof(fn), "ryper.raydonor.%05d.bin",
                    static_cast<int>(std::floor(pm->time/rd_dump_dt_ + 1.0e-6)));
      FILE *fp = std::fopen(fn, "wb");
      if (fp != nullptr) {
        const int hd[3] = {n3g, n1g, kRdNf};
        const double tt = pm->time;
        std::fwrite(hd, sizeof(int), 3, fp);
        std::fwrite(&tt, sizeof(double), 1, fp);
        std::vector<double> xv(rd_x1v_h_.begin(), rd_x1v_h_.end());
        std::vector<double> pv(rd_x3v_h_.begin(), rd_x3v_h_.end());
        std::fwrite(xv.data(), sizeof(double), n1g, fp);
        std::fwrite(pv.data(), sizeof(double), n3g, fp);
        std::fwrite(h.data(), sizeof(float), kRdNf*nn, fp);
        std::fclose(fp);
        std::printf("ry_per_accretor: raydonor dump %s at t %.6e (total removed power %.6g "
                    "code, band)\n", fn, pm->time, rd_qtot_);
      }
    }
    rd_dump_n_++;
    rd_dump_next_ += rd_dump_dt_*std::floor((pm->time - rd_dump_next_)/rd_dump_dt_ + 1.0);
  }
}

//! the per-stage application (before the nozzle/mask resets): eps q + v.f into u(IEN),
//! f into u(IM1), u(IM3), and at the last stage (1 - eps) q dt into the M1 E.  Returns the
//! energy given to the domain this stage (per rank).
Real RayDonorApply(MeshBlockPack *pmbp, Mesh *pm, const Real bdt, const bool last) {
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ni = indcs.nx1, nk = indcs.nx3, nki = nk*ni;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto u0 = pmbp->phydro->u0;
  auto w0 = pmbp->phydro->w0;
  auto *pm1 = pmbp->pradm1;
  auto opac = pm1->opac;
  auto urd = pm1->u0;
  auto qq = rd_q_;
  auto &vol = pmbp->pcoord->volume;
  const Real dtl = pm->dt;
  const bool dgas = rd_dgas_;
  Real esum = 0.0;
  const int ntot = (nmb1 + 1)*(ke - ks + 1)*(je - js + 1)*ni;
  const int nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("ryper_rd_apply", Kokkos::RangePolicy<>(DevExeSpace(), 0, ntot),
  KOKKOS_LAMBDA(const int idx, Real &sum) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
    const int i = idx%ni + is;
    const int c = m*nki + (k - ks)*ni + (i - is);
    const Real qt = qq(c, 0);
    if (!(qt > 0.0)) return;
    const Real ka = opac(m,radm1::M1_OP_E,k,j,i);
    const Real ch = fmax(opac(m,radm1::M1_OP_T,k,j,i), ka);
    const Real eps = (ch > 0.0) ? fmin(ka/ch, 1.0) : 1.0;
    // momentum through the TRANSPORT opacity only (rho kappa_T, as the M1 flux coupling):
    // the removed beam momentum x kappa_T/chi_ext (= 1 for donor_ray_opac = rosseland)
    const Real ft = (ch > 0.0) ? fmin(opac(m,radm1::M1_OP_T,k,j,i)/ch, 1.0) : 1.0;
    const Real fr = ft*qq(c, 1), fp = ft*qq(c, 2);
    const Real wk = fr*w0(m,IVX,k,j,i) + fp*w0(m,IVZ,k,j,i);
    const Real epg = dgas ? eps : 0.0;
    const Real de = bdt*(epg*qt + wk);
    u0(m,IEN,k,j,i) += de;
    u0(m,IM1,k,j,i) += bdt*fr;
    u0(m,IM3,k,j,i) += bdt*fp;
    Real dr = 0.0;
    if (last) {
      dr = dtl*(1.0 - epg)*qt;
      urd(m,radm1::M1_E,k,j,i) += dr;
    }
    sum += (de + dr)*vol(m,k,j,i);
  }, Kokkos::Sum<Real>(esum));
  (void) ie;
  return esum;
}

void RayDonorFinal() {
  rd_rf_ = DvceArray1D<Real>();
  rd_chi_ = DvceArray2D<float>();
  rd_koff_ = DvceArray1D<int>();
  rd_ioff_ = DvceArray1D<int>();
  rd_span_ = DvceArray2D<Real>();
  rd_geo_ = DvceArray4D<float>();
  rd_q_ = DvceArray2D<Real>();
  rd_poly_ = DvceArray2D<Real>();
}
