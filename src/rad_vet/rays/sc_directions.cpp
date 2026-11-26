//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_directions.cpp
//! \brief curved-mesh rays: direction sets and their rotation.
//!
//!  called from: rad_vet_closure.cpp, rad_vet_closure_ring.cpp,
//!    rad_vet_params_closure.cpp, rays/sc_build.cpp, rays/sc_cart_setup.cpp,
//!    rays/sc_cart_sweep.cpp, rays/sc_setup.cpp, rays/sc_tables.cpp, driver/driver.cpp,
//!    driver/dt_retry.cpp, mesh/mesh.cpp
//!  calls into: rad_vet_imex.cpp

#include <math.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_halo.hpp"
#include "utils/recipe.hpp"

namespace radvet {

// <rad_vet>/vet_gd_twin_srw: weight the twin's pattern subtraction per cell by the
// symref weight w = min(J_twin/J, 1) (ScSrW); w = 1 cells are bitwise. Parsed in
// rad_vet.cpp.
// <rad_vet>/vet_rotate_advance: with vet_rotate_avg = K the angle set of block b does
// not repeat with period K, so the averaged window advances through SO(3) / the azimuth
// (parsed in VetRotInit)
bool g_vra_adv = false;

namespace {

// HEALPix ring-scheme centres (12 ns^2) with the moment-fixed weights
void ScDirections(const int ns, std::vector<double> &d) {
  d.clear();
  for (int i = 1; i < 4*ns; ++i) {
    double z, ph0, dph;
    int cnt;
    if (i < ns) {
      z = 1.0 - static_cast<double>(i*i)/(3.0*ns*ns);
      cnt = 4*i; dph = M_PI/(2.0*i); ph0 = -0.5*dph;
    } else if (i <= 3*ns) {
      z = 4.0/3.0 - 2.0*i/(3.0*ns);
      const int s = (i - ns + 1) % 2;
      cnt = 4*ns; dph = M_PI/(2.0*ns); ph0 = -0.5*s*dph;
    } else {
      const int ii = 4*ns - i;
      z = -(1.0 - static_cast<double>(ii*ii)/(3.0*ns*ns));
      cnt = 4*ii; dph = M_PI/(2.0*ii); ph0 = -0.5*dph;
    }
    const double st = std::sqrt(std::max(1.0 - z*z, 0.0));
    for (int j = 1; j <= cnt; ++j) {
      const double p = ph0 + j*dph;
      d.push_back(st*std::cos(p));
      d.push_back(st*std::sin(p));
      d.push_back(z);
      d.push_back(0.0);
    }
  }
  const int n = static_cast<int>(d.size()/4);
  double m2 = 0.0, m4 = 0.0;
  for (int q = 0; q < n; ++q) {
    const double z2 = d[4*q+2]*d[4*q+2];
    m2 += z2/n;
    m4 += z2*z2/n;
  }
  const double c = (1.0/3.0 - m2)/(m4 - m2*m2);
  for (int q = 0; q < n; ++q) {
    d[4*q+3] = (1.0 + c*(d[4*q+2]*d[4*q+2] - m2))/n;
  }
}

// product set: Gauss-Legendre in n_z (nmu nodes) x nphi uniform azimuths per ring,
// alternate rings shifted by half a step (stagger); positive weights, sum w = 1,
// sum w n = 0, sum w nn = I/3 exactly
void ScProduct(const int nmu, const int nphi, const bool stagger,
               std::vector<double> &d) {
  d.clear();
  for (int q = 0; q < nmu; ++q) {       // Newton on P_nmu
    double z = std::cos(M_PI*(q + 0.75)/(nmu + 0.5));
    double pp = 1.0;
    for (int it = 0; it < 100; ++it) {
      double p1 = 1.0, p2 = 0.0;
      for (int l = 1; l <= nmu; ++l) {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0*l - 1.0)*z*p2 - (l - 1.0)*p3)/l;
      }
      pp = nmu*(z*p1 - p2)/(z*z - 1.0);
      const double dz = p1/pp;
      z -= dz;
      if (std::fabs(dz) < 1.0e-15) {break;}
    }
    const double wq = 2.0/((1.0 - z*z)*pp*pp);      // [-1, 1] weight, sum 2
    const double st = std::sqrt(std::max(1.0 - z*z, 0.0));
    const double off = (stagger && (q % 2 == 1)) ? 0.5 : 0.0;
    for (int k = 0; k < nphi; ++k) {
      const double ph = (k + off)*2.0*M_PI/nphi;
      d.push_back(st*std::cos(ph));
      d.push_back(st*std::sin(ph));
      d.push_back(z);
      d.push_back(0.5*wq/nphi);
    }
  }
}

// vet_gd_geo = n >= 1 (stored as vgd_ls = -n): the geodesic set, the vertices of the
// class-I subdivided icosahedron (10 n^2 + 2 directions); weight = solid angle of the
// vertex's spherical Voronoi cell / 4 pi. Full icosahedral symmetry: sum w = 1,
// sum w n = 0, sum w nn = I/3, 4th moments isotropic to round-off (ScMomentCheck)
double ScSphTri(const double *a, const double *b, const double *c) {
  // the solid angle of the spherical triangle of unit vectors a, b, c (Van Oosterom)
  const double tp = a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) +
                    a[2]*(b[0]*c[1] - b[1]*c[0]);
  const double ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  const double bc = b[0]*c[0] + b[1]*c[1] + b[2]*c[2];
  const double ca = c[0]*a[0] + c[1]*a[1] + c[2]*a[2];
  return 2.0*std::fabs(std::atan2(tp, 1.0 + ab + bc + ca));
}

void ScUnit(double *v) {
  const double s = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] /= s; v[1] /= s; v[2] /= s;
}

bool ScGeodesic(const int nl, std::vector<double> &d) {
  if (nl < 1 || nl > 12) {return false;}
  const double g = 0.5*(1.0 + std::sqrt(5.0));
  const double iv[12][3] = {{-1, g, 0}, {1, g, 0}, {-1, -g, 0}, {1, -g, 0},
                            {0, -1, g}, {0, 1, g}, {0, -1, -g}, {0, 1, -g},
                            {g, 0, -1}, {g, 0, 1}, {-g, 0, -1}, {-g, 0, 1}};
  const int fc[20][3] = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
                         {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                         {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
                         {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
  // the vertices: barycentric points (i, j, k)/n of every face, merged by position
  std::vector<std::array<double, 3>> vx;
  std::vector<std::array<int, 3>> tr;
  auto find = [&](double *p) {
    ScUnit(p);
    for (size_t q = 0; q < vx.size(); ++q) {
      const double dx = vx[q][0] - p[0], dy = vx[q][1] - p[1], dz = vx[q][2] - p[2];
      if (dx*dx + dy*dy + dz*dz < 1.0e-20) {return static_cast<int>(q);}
    }
    vx.push_back({p[0], p[1], p[2]});
    return static_cast<int>(vx.size()) - 1;
  };
  for (int f = 0; f < 20; ++f) {
    const double *A = iv[fc[f][0]], *B = iv[fc[f][1]], *C = iv[fc[f][2]];
    std::vector<int> id((nl + 1)*(nl + 2)/2);
    auto at = [&](const int i, const int j) {return i*(nl + 1) - i*(i - 1)/2 + j;};
    for (int i = 0; i <= nl; ++i) {
      for (int j = 0; j <= nl - i; ++j) {
        const int k = nl - i - j;
        double p[3];
        for (int c = 0; c < 3; ++c) {p[c] = (i*A[c] + j*B[c] + k*C[c])/nl;}
        id[at(i, j)] = find(p);
      }
    }
    for (int i = 0; i < nl; ++i) {
      for (int j = 0; j < nl - i; ++j) {
        tr.push_back({id[at(i, j)], id[at(i + 1, j)], id[at(i, j + 1)]});
        if (j + 1 < nl - i) {
          tr.push_back({id[at(i + 1, j)], id[at(i + 1, j + 1)], id[at(i, j + 1)]});
        }
      }
    }
  }
  const int nv = static_cast<int>(vx.size());
  if (nv != 10*nl*nl + 2) {return false;}
  std::vector<double> w(nv, 0.0);
  for (const auto &t : tr) {
    double cc[3];
    const double *a = vx[t[0]].data(), *b = vx[t[1]].data(), *c = vx[t[2]].data();
    const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    cc[0] = u[1]*v[2] - u[2]*v[1];
    cc[1] = u[2]*v[0] - u[0]*v[2];
    cc[2] = u[0]*v[1] - u[1]*v[0];
    if (cc[0]*a[0] + cc[1]*a[1] + cc[2]*a[2] < 0.0) {
      cc[0] = -cc[0]; cc[1] = -cc[1]; cc[2] = -cc[2];
    }
    ScUnit(cc);
    for (int e = 0; e < 3; ++e) {
      const double *p = vx[t[e]].data();
      const double *q1 = vx[t[(e + 1) % 3]].data(), *q2 = vx[t[(e + 2) % 3]].data();
      double m1[3] = {p[0] + q1[0], p[1] + q1[1], p[2] + q1[2]};
      double m2[3] = {p[0] + q2[0], p[1] + q2[1], p[2] + q2[2]};
      ScUnit(m1);
      ScUnit(m2);
      w[t[e]] += ScSphTri(p, m1, cc) + ScSphTri(p, cc, m2);
    }
  }
  double ws = 0.0;
  for (int q = 0; q < nv; ++q) {ws += w[q];}
  d.clear();
  for (int q = 0; q < nv; ++q) {
    d.push_back(vx[q][0]);
    d.push_back(vx[q][1]);
    d.push_back(vx[q][2]);
    d.push_back(w[q]/ws);
  }
  return true;
}

// vet_gd_ls: level-symmetric LQ_N (Carlson; mu_1 of Lewis & Miller), positive weights
// per permutation class from the even moments (sum w nn = I/3 by octahedral symmetry)
bool ScLevelSym(const int nl, std::vector<double> &d) {
  double m1;
  switch (nl) {
    case 6: m1 = 0.2666355; break;
    case 8: m1 = 0.2182179; break;
    case 10: m1 = 0.1893213; break;
    case 12: m1 = 0.1672126; break;
    default: return false;
  }
  const int h = nl/2;
  const double dl = 2.0*(1.0 - 3.0*m1*m1)/(nl - 2);
  std::vector<double> mu(h);
  for (int q = 0; q < h; ++q) {mu[q] = std::sqrt(m1*m1 + dl*q);}
  std::vector<std::vector<int>> cls;     // sorted index triples
  std::vector<int> pc;                   // class of each (unsigned) point
  std::vector<double> pt;
  for (int a = 0; a < h; ++a) {
    for (int b = 0; b < h; ++b) {
      const int c = h - 1 - a - b;
      if (c < 0) {continue;}
      std::vector<int> t = {a, b, c};
      std::sort(t.begin(), t.end());
      int id = -1;
      for (int q = 0; q < static_cast<int>(cls.size()); ++q) {
        if (cls[q] == t) {id = q;}
      }
      if (id < 0) {cls.push_back(t); id = static_cast<int>(cls.size()) - 1;}
      pc.push_back(id);
      pt.push_back(mu[a]); pt.push_back(mu[b]); pt.push_back(mu[c]);
    }
  }
  const int nc = static_cast<int>(cls.size());
  const int np = static_cast<int>(pc.size());
  // A w = rhs over one octant (weights per point sum to 1/8 there)
  std::vector<double> A(nc*nc, 0.0), r(nc);
  for (int row = 0; row < nc; ++row) {
    const int p2 = (row == 0) ? 0 : (2*(row + 1));
    for (int q = 0; q < np; ++q) {A[row*nc + pc[q]] += std::pow(pt[3*q], p2)*8.0;}
    r[row] = 1.0/(p2 + 1.0);
  }
  for (int c = 0; c < nc; ++c) {       // Gauss elimination with partial pivoting
    int pv = c;
    for (int q = c + 1; q < nc; ++q) {
      if (std::fabs(A[q*nc + c]) > std::fabs(A[pv*nc + c])) {pv = q;}
    }
    for (int q = 0; q < nc; ++q) {std::swap(A[c*nc + q], A[pv*nc + q]);}
    std::swap(r[c], r[pv]);
    for (int q = c + 1; q < nc; ++q) {
      const double f = A[q*nc + c]/A[c*nc + c];
      for (int t = c; t < nc; ++t) {A[q*nc + t] -= f*A[c*nc + t];}
      r[q] -= f*r[c];
    }
  }
  std::vector<double> w(nc);
  for (int c = nc - 1; c >= 0; --c) {
    double v = r[c];
    for (int t = c + 1; t < nc; ++t) {v -= A[c*nc + t]*w[t];}
    w[c] = v/A[c*nc + c];
  }
  d.clear();
  for (int q = 0; q < np; ++q) {
    for (int sg = 0; sg < 8; ++sg) {
      d.push_back(((sg & 1) ? -1.0 : 1.0)*pt[3*q]);
      d.push_back(((sg & 2) ? -1.0 : 1.0)*pt[3*q+1]);
      d.push_back(((sg & 4) ? -1.0 : 1.0)*pt[3*q+2]);
      d.push_back(w[pc[q]]);
    }
  }
  return true;
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScBuildSet
//! \brief the vet_gd direction set of the keys: product set, vet_gd_ls, vet_gd_geo or
//! HEALPix (vet_gd_nside); shared by vet_gd and vet_sc_set = gd

void RadiationVet::ScBuildSet(std::vector<double> &d) const {
  if (vgd_gl_nmu > 0) {
    ScProduct(vgd_gl_nmu, vgd_gl_nphi, vgd_gl_stag, d);
  } else if (vgd_ls > 0) {
    if (!ScLevelSym(vgd_ls, d)) {ScFatal("vet_gd_ls must be 6, 8, 10 or 12");}
  } else if (vgd_ls < 0) {   // vet_gd_geo
    if (!ScGeodesic(-vgd_ls, d)) {ScFatal("vet_gd_geo must be 1 .. 12");}
  } else {
    ScDirections(vgd_nside, d);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRotInit
//! \brief rotation keys of the angle set, shared by vet_gd and vet_sc (one clock, one
//! key set); called from the closure parse, curved and Cartesian

void RadiationVet::VetRotInit(ParameterInput *pin) {
  const bool fresh = !global_variable::restart_run;
  // the rotated sets do not repeat (the averaged last-K window advances) on the curved
  // meshes and on the ring
  g_vra_adv = vgd_on || cln_ring > 0;
  {
    const bool rdt_named = pin->DoesParameterExist("rad_vet","vet_gd_rotate_dt");
    const bool rev_named = pin->DoesParameterExist("rad_vet","vet_gd_rotate_every");
    // vet_sc rotates on the same clock with the same keys, unless its set is pinned
    // by boundary beams (implicit_beam_n)
    const bool sc_pin = vet_sc &&
        (pin->DoesParameterExist("rad_vet","implicit_beam_n") &&
         pin->GetInteger("rad_vet","implicit_beam_n") > 0);
    // a pinned vet_sc set gets no default rotation and a named one is fatal (the beams
    // are injected along fixed rays of the set)
    const bool rot_cl = vgd_on || vet_sc;
    bool rot_phys = false;
    if (rot_cl && rdt_named) {
      const std::string rv = pin->GetString("rad_vet","vet_gd_rotate_dt");
      if (rv.compare("auto") == 0) {
        vgd_rot_auto = true;
      } else {
        vgd_rot_dt = pin->GetReal("rad_vet","vet_gd_rotate_dt");
      }
      rot_phys = vgd_rot_auto || (vgd_rot_dt > 0.0);
    } else if (vgd_on && !rev_named && fresh) {
      vgd_rot_auto = vgd_rot_unc = true;
      rot_phys = true;
    } else if (vet_sc && !sc_pin && !rev_named && fresh) {
      vgd_rot_auto = vgd_rot_unc = true;
      rot_phys = true;
    }
    if (rot_phys && vet_sc && sc_pin) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_gd_rotate_dt: the vet_sc set is pinned by "
                << "boundary beams (implicit_beam_n); no rotation"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (rot_phys) {
      vsc_rot = vet_sc;
      vgd_rot = 0;
      pin_rot_ = pin;
      if (vet_sc && !cln_ring) {
        recipe::Rule(pin, "rad_vet", "vet_gd_rotate_slope", "correct");
        recipe::Rule(pin, "rad_vet", "vet_gd_rotate_slope_frac",
                     ParameterInput::RealToString(1.0 - 1.0/std::sqrt(2.0)));
      }
      if (vgd_rot_auto) {
        vgd_rot_auto_fac = 20.0;
        // recorded once resolved; until then the key reads "auto"
        pin->SetString("rad_vet","vet_gd_rotate_dt","auto");
      }
    }
  }
  if (ScRotOn()) {
    recipe::Rule(pin, "rad_vet", "vet_rotate_advance", g_vra_adv ? "true" : "false");
  }
  // the axis: sp turns the set about the polar axis (exact wall maps for HEALPix), the
  // cubed sphere about a quasi-random axis per rotation block, so no direction stays
  // aligned with panel seams
  if (vgd_on) {
    vgd_rot3 = pmy_pack->pmesh->use_cubed_sphere;
    recipe::Rule(pin, "rad_vet", "vet_gd_rotate_axis", vgd_rot3 ? "3d" : "z");
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::RotAngle
//! \brief vet_gd_rotate_every = N > 0: the z-rotation of the set for cycle c, constant
//! per block of N cycles: alpha = P frac((c/N) g), g = (sqrt 5 - 1)/2, P the set's
//! period (deterministic: reruns and restarts bitwise)

Real RadiationVet::RotAngle(const int cyc) const {
  if (!ScRotOn()) {return 0.0;}
  return RotAngleBlk(RotBlock(cyc));
}

Real RadiationVet::RotAngleBlk(const int64_t blk) const {
  // vet_gd_rotate_axis = 3d: the key of block b (decoded by ScRot3) is zigzag(b), or
  // zigzag(b mod K) with vet_rotate_avg = K. vet_rotate_advance: always zigzag(b), so
  // any K consecutive blocks are consecutive R_3 points and the window moves by one set
  // per block. vet_rotate_avg_ema: a never-repeating sequence (3d: key of b; z: golden
  // phase of b)
  if (vra_ema) {
    if (vgd_rot3) {
      const int64_t z = (blk >= 0) ? 2*blk : (-2*blk - 1);
      return static_cast<Real>(z);
    }
    const double per = vgd_zper;
    const double b = static_cast<double>(blk);
    const double g = 0.5*(std::sqrt(5.0) - 1.0);
    const double a = per*(b*g - std::floor(b*g));
    return static_cast<Real>(a);
  }
  if (vgd_rot3) {
    int64_t r = blk;
    if (vra_on && vra_k >= 2 && !g_vra_adv) {r = ((blk % vra_k) + vra_k) % vra_k;}
    const int64_t z = (r >= 0) ? 2*r : (-2*r - 1);
    return static_cast<Real>(z);
  }
  // vet_rotate_avg = K: angle = period x RotPhase, (b mod K)/K, so any K consecutive
  // blocks are K equispaced rotations (the K-fold refined azimuth set)
  if (vra_on) {
    const double per = vgd_zper;
    return static_cast<Real>(per*RotPhase(blk));
  }
  const double b = static_cast<double>(blk);
  const double g = 0.5*(std::sqrt(5.0) - 1.0);
  // a product set keeps pi/2 here
  return ((vgd_gl_nmu > 0) ? 0.5*M_PI : vgd_zper)*(b*g - std::floor(b*g));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScRot3
//! \brief vet_gd_rotate_axis = 3d: the rotation matrix (row major) of block key `key`.
//! Shoemake (1992) uniform rotation from the R_3 point u_b = frac(b (1/g, 1/g^2,
//! 1/g^3)) (Roberts 2018), times the inverse of b = 0's (block 0 unrotated). A pure
//! function of the key: same bits on every rank and restart.

void RadiationVet::ScRot3(const Real key, double r[9]) const {
  const int64_t z = static_cast<int64_t>(std::llround(static_cast<double>(key)));
  const int64_t b = (z % 2 == 0) ? (z/2) : (-(z + 1)/2);
  auto quat = [](const int64_t bb, double q[4]) {
    const double g = 1.2207440846057594754;   // the real root of x^4 = x + 1
    const double a[3] = {1.0/g, 1.0/(g*g), 1.0/(g*g*g)};
    double u[3];
    for (int c = 0; c < 3; ++c) {
      const double x = static_cast<double>(bb)*a[c];
      u[c] = x - std::floor(x);
    }
    const double s1 = std::sqrt(1.0 - u[0]), s2 = std::sqrt(u[0]);
    q[0] = s1*std::sin(2.0*M_PI*u[1]);
    q[1] = s1*std::cos(2.0*M_PI*u[1]);
    q[2] = s2*std::sin(2.0*M_PI*u[2]);
    q[3] = s2*std::cos(2.0*M_PI*u[2]);
  };
  auto mat = [](const double q[4], double m[9]) {   // q = (x, y, z, w), unit
    const double x = q[0], y = q[1], zz = q[2], w = q[3];
    m[0] = 1.0 - 2.0*(y*y + zz*zz); m[1] = 2.0*(x*y - zz*w);  m[2] = 2.0*(x*zz + y*w);
    m[3] = 2.0*(x*y + zz*w);  m[4] = 1.0 - 2.0*(x*x + zz*zz); m[5] = 2.0*(y*zz - x*w);
    m[6] = 2.0*(x*zz - y*w);  m[7] = 2.0*(y*zz + x*w);  m[8] = 1.0 - 2.0*(x*x + y*y);
  };
  if (b == 0) {
    for (int c = 0; c < 9; ++c) {r[c] = (c % 4 == 0) ? 1.0 : 0.0;}
    return;
  }
  double qb[4], q0[4], mb[9], m0[9];
  quat(b, qb);
  quat(0, q0);
  mat(qb, mb);
  mat(q0, m0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double v = 0.0;
      for (int k = 0; k < 3; ++k) {v += mb[3*i + k]*m0[3*j + k];}   // mb m0^T
      r[3*i + j] = v;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::RotPhase
//! \brief rotation of block blk as a fraction of the set's period in [0, 1): golden
//! frac(b g), or with vet_rotate_avg = K the stratified (b mod K)/K

double RadiationVet::RotPhase(const int64_t blk) const {
  // vet_rotate_advance: frac(b (1 + gam)/K), gam = g/(4 (K-1)): any K consecutive phases
  // are equispaced to within 0.16/K of the period, and the window drifts by gam/K per
  // block without repeating
  if (vra_on && vra_k >= 2 && g_vra_adv) {
    const double g = 0.5*(std::sqrt(5.0) - 1.0);
    const double kk = static_cast<double>(vra_k);
    const double x = static_cast<double>(blk)*(1.0 + g/(4.0*(kk - 1.0)))/kk;
    return x - std::floor(x);
  }
  if (vra_on && vra_k >= 2) {
    const int64_t r = ((blk % vra_k) + vra_k) % vra_k;
    return static_cast<double>(r)/static_cast<double>(vra_k);
  }
  const double b = static_cast<double>(blk);
  const double g = 0.5*(std::sqrt(5.0) - 1.0);
  return b*g - std::floor(b*g);
}

//----------------------------------------------------------------------------------------
//! \fn int64_t RadiationVet::RotBlock
//! \brief rotation block of the step of cycle cyc. Cycle clock (vet_gd_rotate_every =
//! N): cyc/N. Physical clock (vet_gd_rotate_dt > 0): floor((t + dt/2 - t0)/dt_rot),
//! the block of the step midpoint, shared by all stages and builds of the step. With
//! vet_gd_rotate_align every jump sits on a step boundary. Restart-safe.

int64_t RadiationVet::RotBlock(const int cyc) const {
  if (vgd_rot_dt > 0.0) {
    const Mesh *pm = pmy_pack->pmesh;
    const Real dt = (pm->dt > 0.0 && pm->dt < vgd_rot_dt) ? pm->dt : 0.0;
    return static_cast<int64_t>(std::floor((pm->time + 0.5*dt)/vgd_rot_dt));
  }
  return (vgd_rot > 0) ? static_cast<int64_t>(cyc/vgd_rot) : 0;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::RotLastStep
//! \brief is the current step the last of its rotation block? Physical clock: t + dt
//! reaches t0 + (b+1) dt_rot; cycle clock: (cyc+1) % N == 0. Shared with vet_sc.

bool RadiationVet::RotLastStep() const {
  const Mesh *pm = pmy_pack->pmesh;
  if (vgd_rot_dt > 0.0) {
    const int64_t b = RotBlock(pm->ncycle);
    const double tk = (static_cast<double>(b) + 1.0)*vgd_rot_dt;
    return pm->time + pm->dt >= tk - 1.0e-9*vgd_rot_dt;
  }
  if (vgd_rot > 0) {return ((pm->ncycle + 1) % vgd_rot) == 0;}
  return false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotNewBlock
//! \brief (Driver, before Time2RotCorrect) on the first step of a new rotation block
//! the stored hesdirk2 slope K1 belongs to the old set; be resets it (t2_ok = false).
//! The first call takes the block before t as the last one (restart on a boundary).

void RadiationVet::RotNewBlock() {
  if (!(vgd_rot_dt > 0.0) || !(vgd_on || vsc_rot)) {return;}
  const int64_t b = RotBlock(pmy_pack->pmesh->ncycle);
  if (vgd_rot_lastblk == INT64_MIN) {
    const Mesh *pm = pmy_pack->pmesh;
    vgd_rot_lastblk = static_cast<int64_t>(
        std::floor(pm->time/vgd_rot_dt - 1.0e-9));
  }
  // a Cartesian set without the ring keeps the slope and asks Time2RotCorrect for its
  // jump; the cross-fade has no jump; clean_rot_ring = jump: a backward-Euler step
  if (vgd_xf && vgd_on) {
    vgd_rot_lastblk = b;
    return;
  }
  if (b != vgd_rot_lastblk && time_scheme == VET_TIME_HESDIRK2) {
    if (vet_sc && !cln_ring && Time2Active()) t2_rotcorr = true;
    if (cln_ring == 2) t2_ok = false;
  }
  vgd_rot_lastblk = b;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::RotAlignDt
//! \brief vet_gd_rotate_align (Mesh::NewTimeStep, after the global min): reach the next
//! boundary k dt_rot in n = ceil(rem/dt) equal steps, dt' = rem/n >= dt/2.
//! vet_gd_rotate_stretch = s > 0: n - 1 steps when dt' <= (1 + s) dt.

Real RadiationVet::RotAlignDt(const Real time, const Real dt) {
  // vet_gd_rotate_dt = auto: vgd_rot_auto_fac x the first global dt, recorded as a
  // number (restarts keep it)
  if (vgd_rot_auto && !(vgd_rot_dt > 0.0) && dt > 0.0) {
    vgd_rot_dt = vgd_rot_auto_fac*(vgd_rot_unc ? pmy_pack->pmesh->dt_uncapped : dt);
    vgd_rot_auto = false;
    if (pin_rot_ != nullptr) {pin_rot_->SetReal("rad_vet","vet_gd_rotate_dt",vgd_rot_dt);}
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_rotate_dt = auto: " << vgd_rot_auto_fac
                << " x the first "
                << (vgd_rot_unc ? "CFL dt (no dt_fixed / dt_retry cap) " : "dt ")
                << (vgd_rot_unc ? pmy_pack->pmesh->dt_uncapped : dt) << " = "
                << ParameterInput::RealToString(vgd_rot_dt)
                << " (recorded)" << std::endl;
    }
  }
  if (!(vgd_rot_dt > 0.0) || !(dt > 0.0)) {return dt;}
  const double x = time/vgd_rot_dt;
  double k = std::floor(x) + 1.0;
  // already on a boundary (to round-off): the next one
  if ((k - x)*vgd_rot_dt < 1.0e-9*dt) {k += 1.0;}
  const double rem = k*vgd_rot_dt - time;
  double n = std::ceil(rem/dt*(1.0 - 1.0e-9));
  if (n > 1.0e6) {return dt;}
  if (n > 1.0 && rem/(n - 1.0) <= 1.02*dt) n -= 1.0;   // (dt may grow by 2%)
  const Real dtn = static_cast<Real>(rem/std::max(n, 1.0));
  if (dtn < dt) {vgd_rot_nalign += 1.0;}
  return dtn;
}

//----------------------------------------------------------------------------------------
//! \fn std::vector<double> RadiationVet::ScRotate
//! \brief the set at angle (or block key) alpha, packed as vgd_base (x, y, z, w):
//! rotated about z, or with vet_gd_rotate_axis = 3d by R_b (ScRot3)

std::vector<double> RadiationVet::ScRotate(const Real alpha) const {
  const int n = vgd_n;
  std::vector<double> d(vgd_base);
  if (vgd_rot3) {
    double rm[9];
    ScRot3(alpha, rm);
    for (int q = 0; q < n; ++q) {
      const double x = vgd_base[4*q], y = vgd_base[4*q+1], z = vgd_base[4*q+2];
      double v[3];
      for (int c = 0; c < 3; ++c) {v[c] = rm[3*c]*x + rm[3*c+1]*y + rm[3*c+2]*z;}
      const double nv = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      for (int c = 0; c < 3; ++c) {d[4*q+c] = v[c]/nv;}
    }
  } else {
    const double ca = std::cos(alpha), sa = std::sin(alpha);
    for (int q = 0; q < n; ++q) {
      const double x = vgd_base[4*q], y = vgd_base[4*q+1];
      d[4*q] = ca*x - sa*y;
      d[4*q+1] = sa*x + ca*y;
    }
  }
  return d;
}

}  // namespace radvet
