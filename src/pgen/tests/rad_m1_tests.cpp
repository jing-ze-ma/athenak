//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_tests.cpp
//! \brief test problems for the grey photon M1 module, selected by <problem>/m1_test.
//!
//!   beam    (design sect. 8, T1): 2-D unit box, no opacity, a beam entering at 45
//!           degrees through a patch of the inner-x1 face.  The rest of the domain
//!           starts (and the other faces stay) at the vacuum floor.
//!   pulse1d (convergence check): 1-D periodic Gaussian in E with f = |F|/(cE) = 1
//!           EXACTLY.  At f = 1 the M1 closure gives chi = 1, P = E, the two
//!           eigenvalues merge at +c, and the system reduces to dE/dt + c dE/dx = 0:
//!           the pulse translates rigidly at c, so the exact solution after one box
//!           crossing is the initial condition.  That is the only free-streaming case
//!           that is exactly solvable under M1 (with f < 1 the profile spreads), which
//!           is why f = 1 is used here rather than 0.9.

#include <math.h>

#include <algorithm>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "pgen/pgen.hpp"

namespace {
// beam parameters, needed by the user boundary function
Real m1_beam_e = 1.0;
Real m1_beam_n1 = 1.0/M_SQRT2;
Real m1_beam_n2 = 1.0/M_SQRT2;
Real m1_beam_f = 1.0 - 1.0e-6;
Real m1_beam_y0 = 0.0;
Real m1_beam_y1 = 0.125;
int  m1_test_id = 0;   // 0 = beam, 1 = pulse1d
} // namespace

// prototype for the user BC
void RadM1BeamBC(Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::RadiationM1Tests()

void ProblemGenerator::RadiationM1Tests(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pradm1 == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "pgen_name = rad_m1_beam requires a <rad_m1> block" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string test = pin->GetOrAddString("problem","m1_test","beam");
  Real cl = pmbp->pradm1->c_light;
  Real efl = pmbp->pradm1->e_floor;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int nx1 = indcs.nx1;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto u0 = pmbp->pradm1->u0;

  if (test.compare("beam") == 0) {
    m1_test_id = 0;
    // beam direction and patch, all optional
    Real theta = pin->GetOrAddReal("problem","beam_angle",0.25*M_PI);
    m1_beam_n1 = cos(theta);
    m1_beam_n2 = sin(theta);
    m1_beam_e = pin->GetOrAddReal("problem","beam_e",1.0);
    m1_beam_f = pin->GetOrAddReal("problem","beam_f",1.0-1.0e-6);
    m1_beam_y0 = pin->GetOrAddReal("problem","beam_y0",0.0);
    m1_beam_y1 = pin->GetOrAddReal("problem","beam_y1",0.125);
    user_bcs_func = RadM1BeamBC;

    // vacuum everywhere
    par_for("m1_beam_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,radm1::M1_E,k,j,i) = efl;
      u0(m,radm1::M1_F1,k,j,i) = 0.0;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("pulse1d") == 0) {
    m1_test_id = 1;
    Real amp = pin->GetOrAddReal("problem","pulse_amp",1.0);
    Real wid = pin->GetOrAddReal("problem","pulse_width",0.1);
    Real x0 = pin->GetOrAddReal("problem","pulse_x0",0.5);
    Real bg = pin->GetOrAddReal("problem","pulse_bg",1.0e-4);

    par_for("m1_pulse_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real e = bg + amp*exp(-SQR((x1v - x0)/wid));
      u0(m,radm1::M1_E,k,j,i) = e;
      u0(m,radm1::M1_F1,k,j,i) = cl*e;   // f = 1 exactly: rigid translation at +c
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/m1_test = '" << test << "' not implemented "
      << "(beam | pulse1d)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1BeamBC()
//! \brief inner-x1 inflow PATCH for the beam test: inside [beam_y0, beam_y1] the ghost
//! cells hold E = beam_e and F = c E f n (f just below 1, so the closure is at chi ~ 1
//! without ever dividing by a vanishing sqrt(4-3f^2)); outside the patch, and on the
//! inner-x2 face, they hold the vacuum floor so nothing enters.  Outer faces are left to
//! the <rad_m1> physical BCs (vacuum) set in the input file.

void RadM1BeamBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pradm1 == nullptr) return;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &js = indcs.js;
  int nx2 = indcs.nx2;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto u0 = pmbp->pradm1->u0;

  Real cl = pmbp->pradm1->c_light;
  Real efl = pmbp->pradm1->e_floor;
  Real be = m1_beam_e, bf = m1_beam_f, bn1 = m1_beam_n1, bn2 = m1_beam_n2;
  Real by0 = m1_beam_y0, by1 = m1_beam_y1;

  par_for("m1_beam_bc", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),
  KOKKOS_LAMBDA(int m, int k, int j) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      bool inbeam = (x2v >= by0 && x2v <= by1);
      for (int i=0; i<ng; ++i) {
        if (inbeam) {
          u0(m,radm1::M1_E,k,j,is-i-1) = be;
          u0(m,radm1::M1_F1,k,j,is-i-1) = cl*be*bf*bn1;
          u0(m,radm1::M1_F2,k,j,is-i-1) = cl*be*bf*bn2;
          u0(m,radm1::M1_F3,k,j,is-i-1) = 0.0;
        } else {
          u0(m,radm1::M1_E,k,j,is-i-1) = efl;
          u0(m,radm1::M1_F1,k,j,is-i-1) = 0.0;
          u0(m,radm1::M1_F2,k,j,is-i-1) = 0.0;
          u0(m,radm1::M1_F3,k,j,is-i-1) = 0.0;
        }
      }
    }
  });
  return;
}
