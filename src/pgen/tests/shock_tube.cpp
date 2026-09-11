//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file shock_tube.cpp
//! \brief Problem generator for shock tube (1-D Riemann) problems in both hydro and MHD.
//! Works for both non-relativistic and relativistic dynamics in flat (Minkowski)
//! spacetimes.  Can be used to test GR, but metric must be Minkowski.
//!
//! Works by initializing plane-parallel shock along x1 (in 1D, 2D, 3D), along x2
//! (in 2D, 3D), and along x3 (in 3D).  Shock must be along a coordinate directions,
//! i.e. shocks propagating along an angle inclined to grid are not implemented.

#include <iostream>
#include <sstream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "coordinates/adm.hpp"
#include "pgen/pgen.hpp"

namespace {

int shk_dir;
void SetADMVariablesToSchwarzschild(MeshBlockPack *pmbp);

}

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::ShockTube_()
//! \brief Problem Generator for the shock tube (Riemann problem) tests

void ProblemGenerator::ShockTube(ParameterInput *pin, const bool restart) {
  if (restart) return;
  // parse shock direction: {1,2,3} -> {x1,x2,x3}
  shk_dir = pin->GetInteger("problem","shock_dir");
  if (shk_dir < 1 || shk_dir > 3) {
    // Invaild input value for shk_dir
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "shock_dir=" <<shk_dir<< " must be either 1,2, or 3" << std::endl;
    exit(EXIT_FAILURE);
  }
  // set indices of parallel and perpendicular velocities
  int ivx = shk_dir;
  int ivy = IVX + ((ivx - IVX) + 1)%3;
  int ivz = IVX + ((ivx - IVX) + 2)%3;

  // parse shock location (must be inside grid; set L/R states equal to each other to
  // initialize uniform initial conditions))
  Real xshock = pin->GetReal("problem","xshock");
  if (shk_dir == 1 && (xshock < pmy_mesh_->mesh_size.x1min ||
                       xshock > pmy_mesh_->mesh_size.x1max)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "xshock=" << xshock << " lies outside x1 domain" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (shk_dir == 2 && (xshock < pmy_mesh_->mesh_size.x2min ||
                       xshock > pmy_mesh_->mesh_size.x2max)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "xshock=" << xshock << " lies outside x2 domain" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (shk_dir == 3 && (xshock < pmy_mesh_->mesh_size.x3min ||
                       xshock > pmy_mesh_->mesh_size.x3max)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "xshock=" << xshock << " lies outside x3 domain" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for the kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;

  // Initialize Hydro variables -------------------------------
  if (pmbp->phydro != nullptr) {
    auto &eos = pmbp->phydro->peos->eos_data;
    // Parse left state read from input file: d,vx,vy,vz,[P]
    HydPrim1D wl,wr;
    wl.d  = pin->GetReal("problem","dl");
    wl.vx = pin->GetReal("problem","ul");
    wl.vy = pin->GetReal("problem","vl");
    wl.vz = pin->GetReal("problem","wl");
    // <problem>/pl is a PRESSURE. Under a general EOS the internal energy that goes with
    // it is not p/(gamma-1) and must come from the EOS; that conversion is done inside
    // the kernel below, because the tabulated EOS lives in device memory.
    Real pl_ = pin->GetReal("problem","pl");
    wl.e  = pl_/(eos.gamma - 1.0);
    // compute Lorentz factor (needed for SR/GR)
    Real u0l = 1.0;
    if (pmbp->pcoord->is_special_relativistic || pmbp->pcoord->is_general_relativistic ||
        pmbp->pcoord->is_dynamical_relativistic) {
      u0l = 1.0/sqrt( 1.0 - (SQR(wl.vx) + SQR(wl.vy) + SQR(wl.vz)) );
    }
    if (pmbp->pcoord->is_dynamical_relativistic) {
      wl.e = pin->GetReal("problem","pl");
    }
    Real yl = pin->GetOrAddReal("problem","yl",1.0);

    // Parse right state read from input file: d,vx,vy,vz,[P]
    wr.d  = pin->GetReal("problem","dr");
    wr.vx = pin->GetReal("problem","ur");
    wr.vy = pin->GetReal("problem","vr");
    wr.vz = pin->GetReal("problem","wr");
    Real pr_ = pin->GetReal("problem","pr");
    wr.e  = pr_/(eos.gamma - 1.0);
    // compute Lorentz factor (needed for SR/GR)
    Real u0r = 1.0;
    if (pmbp->pcoord->is_special_relativistic || pmbp->pcoord->is_general_relativistic ||
        pmbp->pcoord->is_dynamical_relativistic) {
      u0r = 1.0/sqrt( 1.0 - (SQR(wr.vx) + SQR(wr.vy) + SQR(wr.vz)) );
    }
    if (pmbp->pcoord->is_dynamical_relativistic) {
      wr.e = pin->GetReal("problem","pr");
    }
    Real yr = pin->GetOrAddReal("problem","yr",0.0);

    auto &w0 = pmbp->phydro->w0;
    auto &nscal = pmbp->phydro->nscalars;
    auto shk_dir_ = shk_dir;
    par_for("pgen_shock1", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m,int k, int j, int i) {
      Real x;
      if (shk_dir_ == 1) {
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        int nx1 = indcs.nx1;
        x = CellCenterX(i-is, nx1, x1min, x1max);
      } else if (shk_dir_ == 2) {
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        int nx2 = indcs.nx2;
        x = CellCenterX(j-js, nx2, x2min, x2max);
      } else {
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        x = CellCenterX(k-ks, nx3, x3min, x3max);
      }

      // in SR/GR, primitive variables use spatial components of 4-vel u^i = gamma * v^i
      if (x < xshock) {
        w0(m,IDN,k,j,i) = wl.d;
        w0(m,ivx,k,j,i) = wl.vx*u0l;
        w0(m,ivy,k,j,i) = wl.vy*u0l;
        w0(m,ivz,k,j,i) = wl.vz*u0l;
        w0(m,IEN,k,j,i) = eos.IsGeneral() ? eos.EnergyFromPressure(wl.d, pl_) : wl.e;
        for (int r=0; r<nscal; ++r) {
          w0(m,IYF+r,k,j,i) = yl;
        }
      } else {
        w0(m,IDN,k,j,i) = wr.d;
        w0(m,ivx,k,j,i) = wr.vx*u0r;
        w0(m,ivy,k,j,i) = wr.vy*u0r;
        w0(m,ivz,k,j,i) = wr.vz*u0r;
        w0(m,IEN,k,j,i) = eos.IsGeneral() ? eos.EnergyFromPressure(wr.d, pr_) : wr.e;
        for (int r=0; r<nscal; ++r) {
          w0(m,IYF+r,k,j,i) = yr;
        }
      }
    });

    // Convert primitives to conserved
    auto &u0 = pmbp->phydro->u0;
    if (pmbp->padm == nullptr) {
      pmbp->phydro->peos->PrimToCons(w0, u0, is, ie, js, je, ks, ke);
    }
  } // End initialization of Hydro variables

  // Initialize MHD variables -------------------------------
  if (pmbp->pmhd != nullptr) {
    auto &eos = pmbp->pmhd->peos->eos_data;
    // Parse left state read from input file: d,vx,vy,vz,[P]
    MHDPrim1D wl,wr;
    wl.d  = pin->GetReal("problem","dl");
    wl.vx = pin->GetReal("problem","ul");
    wl.vy = pin->GetReal("problem","vl");
    wl.vz = pin->GetReal("problem","wl");
    // see the hydro block above: <problem>/pl is a pressure, converted to an internal
    // energy in the kernel when the EOS is general
    Real pl_ = pin->GetReal("problem","pl");
    wl.e  = pl_/(eos.gamma - 1.0);
    wl.by = pin->GetReal("problem","byl");
    wl.bz = pin->GetReal("problem","bzl");
    Real bx_l = pin->GetReal("problem","bxl");
    // compute Lorentz factor (needed for SR/GR)
    Real u0l = 1.0;
    if (pmbp->pcoord->is_special_relativistic || pmbp->pcoord->is_general_relativistic ||
        pmbp->pcoord->is_dynamical_relativistic) {
      u0l = 1.0/sqrt( 1.0 - (SQR(wl.vx) + SQR(wl.vy) + SQR(wl.vz)) );
    }
    if (pmbp->pcoord->is_dynamical_relativistic) {
      wl.e = pin->GetReal("problem", "pl");
    }
    Real yl = pin->GetOrAddReal("problem","yl",1.0);

    // Parse right state read from input file: d,vx,vy,vz,[P]
    wr.d  = pin->GetReal("problem","dr");
    wr.vx = pin->GetReal("problem","ur");
    wr.vy = pin->GetReal("problem","vr");
    wr.vz = pin->GetReal("problem","wr");
    Real pr_ = pin->GetReal("problem","pr");
    wr.e  = pr_/(eos.gamma - 1.0);
    wr.by = pin->GetReal("problem","byr");
    wr.bz = pin->GetReal("problem","bzr");
    Real bx_r = pin->GetReal("problem","bxr");
    // compute Lorentz factor (needed for SR/GR)
    Real u0r = 1.0;
    if (pmbp->pcoord->is_special_relativistic || pmbp->pcoord->is_general_relativistic ||
        pmbp->pcoord->is_dynamical_relativistic) {
      u0r = 1.0/sqrt( 1.0 - (SQR(wr.vx) + SQR(wr.vy) + SQR(wr.vz)) );
    }
    if (pmbp->pcoord->is_dynamical_relativistic) {
      wr.e = pin->GetReal("problem", "pr");
    }
    Real yr = pin->GetOrAddReal("problem","yr",0.0);

    auto &w0 = pmbp->pmhd->w0;
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc0 = pmbp->pmhd->bcc0;
    auto &nscal = pmbp->pmhd->nscalars;
    auto shk_dir_ = shk_dir;
    par_for("pgen_shock1", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m,int k, int j, int i) {
      Real x,bxl,byl,bzl,bxr,byr,bzr;
      if (shk_dir_ == 1) {
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        int nx1 = indcs.nx1;
        x = CellCenterX(i-is, nx1, x1min, x1max);
        bxl = bx_l; byl = wl.by; bzl = wl.bz;
        bxr = bx_r; byr = wr.by; bzr = wr.bz;
      } else if (shk_dir_ == 2) {
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        int nx2 = indcs.nx2;
        x = CellCenterX(j-js, nx2, x2min, x2max);
        bxl = wl.bz; byl = bx_l; bzl = wl.by;
        bxr = wr.bz; byr = bx_r; bzr = wr.by;
      } else {
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        x = CellCenterX(k-ks, nx3, x3min, x3max);
        bxl = wl.by; byl = wl.bz; bzl = bx_l;
        bxr = wr.by; byr = wr.bz; bzr = bx_r;
      }

      // in SR/GR, primitive variables use spatial components of 4-vel u^i = gamma * v^i
      if (x < xshock) {
        w0(m,IDN,k,j,i) = wl.d;
        w0(m,ivx,k,j,i) = wl.vx*u0l;
        w0(m,ivy,k,j,i) = wl.vy*u0l;
        w0(m,ivz,k,j,i) = wl.vz*u0l;
        w0(m,IEN,k,j,i) = eos.IsGeneral() ? eos.EnergyFromPressure(wl.d, pl_) : wl.e;
        for (int r=0; r<nscal; ++r) {
          w0(m,IYF+r,k,j,i) = yl;
        }
        b0.x1f(m,k,j,i) = bxl;
        b0.x2f(m,k,j,i) = byl;
        b0.x3f(m,k,j,i) = bzl;
        if (i==ie) {b0.x1f(m,k,j,i+1) = bxl;}
        if (j==je) {b0.x2f(m,k,j+1,i) = byl;}
        if (k==ke) {b0.x3f(m,k+1,j,i) = bzl;}
        bcc0(m,IBX,k,j,i) = bxl;
        bcc0(m,IBY,k,j,i) = byl;
        bcc0(m,IBZ,k,j,i) = bzl;
      } else {
        w0(m,IDN,k,j,i) = wr.d;
        w0(m,ivx,k,j,i) = wr.vx*u0r;
        w0(m,ivy,k,j,i) = wr.vy*u0r;
        w0(m,ivz,k,j,i) = wr.vz*u0r;
        w0(m,IEN,k,j,i) = eos.IsGeneral() ? eos.EnergyFromPressure(wr.d, pr_) : wr.e;
        for (int r=0; r<nscal; ++r) {
          w0(m,IYF+r,k,j,i) = yr;
        }
        b0.x1f(m,k,j,i) = bxr;
        b0.x2f(m,k,j,i) = byr;
        b0.x3f(m,k,j,i) = bzr;
        if (i==ie) {b0.x1f(m,k,j,i+1) = bxr;}
        if (j==je) {b0.x2f(m,k,j+1,i) = byr;}
        if (k==ke) {b0.x3f(m,k+1,j,i) = bzr;}
        bcc0(m,IBX,k,j,i) = bxr;
        bcc0(m,IBY,k,j,i) = byr;
        bcc0(m,IBZ,k,j,i) = bzr;
      }
    });
    // CUBED SPHERE: replace the field with a divergence-free AZIMUTHAL one -------------
    //
    // The uniform face values written above are NOT divergence-free on this grid (every
    // face area varies from face to face), and a field with B.rhat != 0 makes a
    // reflecting radial wall an ill-posed Riemann problem that leaks mass at O(B^2).  So
    // on the cubed sphere <problem>/bazi replaces them with B = bazi*(-y, x, 0), purely
    // azimuthal about the z axis, hence B.rhat = 0 exactly.  It is laid down as the
    // DISCRETE CURL of its vector potential on the very edges the CT update integrates
    // over, so div B is round-off BY CONSTRUCTION; the gauge is
    // A = (1/3)*bazi*r^2*sin(theta) thetahat, which has no radial component, so each
    // tangential face's circulation is exact over its two tangential edges alone.  This
    // mirrors CSTestBlastFaces in the cs_test problem generator, where the construction
    // is gated on four independent numbers.
    // Strictly ADDITIVE: with <problem>/bazi absent or zero nothing here runs and the
    // uniform bxl/byl/bzl above stand exactly as they did (they are not
    // divergence-free on this grid -- that is then the caller's problem, as it always
    // was).
    const Real bazi = pin->GetOrAddReal("problem","bazi",0.0);
    if (pmy_mesh_->use_cubed_sphere && bazi != 0.0) {
      auto &mbp = pmbp->pmb->mb_panel;
      auto &x1f_ = pmbp->pcoord->xx1f;
      auto &x1v_ = pmbp->pcoord->x1v;
      auto &ar1 = pmbp->pcoord->area.x1f;
      auto &ar2 = pmbp->pcoord->area.x2f;
      auto &ar3 = pmbp->pcoord->area.x3f;
      auto &dxe2 = pmbp->pcoord->dxedge.x2e;
      auto &dxe3 = pmbp->pcoord->dxedge.x3e;
      auto &ccell = pmbp->pcoord->cos_cell;
      auto &scell = pmbp->pcoord->sin_cell;
      const int ng = indcs.ng;
      const int n1m1 = indcs.nx1 + 2*ng - 1;
      const int n2m1 = indcs.nx2 + 2*ng - 1;
      const int n3m1 = indcs.nx3 + 2*ng - 1;
      par_for("pgen_shock_csb", DevExeSpace(), 0, (pmbp->nmb_thispack-1),
              0, n3m1, 0, n2m1, 0, n1m1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        const int p = mbp.d_view(m);
        const Real x2mn = size.d_view(m).x2min, x2mx = size.d_view(m).x2max;
        const Real x3mn = size.d_view(m).x3min, x3mx = size.d_view(m).x3max;
        // A on one edge, projected on that edge's own UNIT tangent -- the component
        // dxedge*A consumes.  along_xi selects the xi edge (r face, xi CENTRE, eta face)
        // or the eta edge (r face, xi face, eta CENTRE).
        auto Aedge = [&](const int ii, const int jj, const int kk, const bool along_xi) {
          const Real rf = x1f_(m,ii);
          const Real xi = 0.25*M_PI*(along_xi
              ? CellCenterX(jj-js, indcs.nx2, x2mn, x2mx)
              : LeftEdgeX(jj-js, indcs.nx2, x2mn, x2mx));
          const Real et = 0.25*M_PI*(along_xi
              ? LeftEdgeX(kk-ks, indcs.nx3, x3mn, x3mx)
              : CellCenterX(kk-ks, indcs.nx3, x3mn, x3mx));
          Real qh[3], e1[3], e2[3];
          cubed_sphere::PanelToCart(p, xi, et, qh);
          cubed_sphere::PanelTangents(p, xi, et, e1, e2);
          const Real aa = bazi*rf*rf/3.0;
          const Real ax =  aa*qh[0]*qh[2];
          const Real ay =  aa*qh[1]*qh[2];
          const Real az = -aa*(qh[0]*qh[0] + qh[1]*qh[1]);
          const Real *t = along_xi ? e1 : e2;
          const Real tn = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
          return (ax*t[0] + ay*t[1] + az*t[2])/tn;
        };
        // B.n = (1/area) * circulation of A around the face, exactly as mhd_ct.cpp
        auto setb1 = [&](const int ifc) {
          b0.x1f(m,k,j,ifc) = (dxe3(m,k,j+1,ifc)*Aedge(ifc,j+1,k,false)
                             - dxe3(m,k,j  ,ifc)*Aedge(ifc,j  ,k,false)
                             - dxe2(m,k+1,j,ifc)*Aedge(ifc,j,k+1,true)
                             + dxe2(m,k  ,j,ifc)*Aedge(ifc,j,k  ,true))/ar1(m,k,j,ifc);
        };
        auto setb2 = [&](const int jfc) {
          b0.x2f(m,k,jfc,i) = -(dxe3(m,k,jfc,i+1)*Aedge(i+1,jfc,k,false)
                              - dxe3(m,k,jfc,i  )*Aedge(i  ,jfc,k,false))/ar2(m,k,jfc,i);
        };
        auto setb3 = [&](const int kfc) {
          b0.x3f(m,kfc,j,i) = (dxe2(m,kfc,j,i+1)*Aedge(i+1,j,kfc,true)
                             - dxe2(m,kfc,j,i  )*Aedge(i  ,j,kfc,true))/ar3(m,kfc,j,i);
        };
        setb1(i);
        setb2(j);
        setb3(k);
        if (i == n1m1) { setb1(i+1); }
        if (j == n2m1) { setb2(j+1); }
        if (k == n3m1) { setb3(k+1); }
        // the cell-centred field, in the ORTHONORMAL frame {rhat, e_xi, (e_eta -
        // c e_xi)/s} that bcc0 is stored in on this grid (see
        // Coordinates::GnomonicEquiangleRaiseVelMHD), so that the PrimToCons below adds
        // the right magnetic energy.  It is rebuilt from the faces at the first
        // ConsToPrim anyway; only |B|^2 has to be right here.
        const Real rc = x1v_(m,i);
        const Real xc = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, x2mn, x2mx);
        const Real ec = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, x3mn, x3mx);
        Real qh[3], e1[3], e2[3];
        cubed_sphere::PanelToCart(p, xc, ec, qh);
        cubed_sphere::PanelTangents(p, xc, ec, e1, e2);
        const Real bx = -bazi*rc*qh[1], by = bazi*rc*qh[0], bz = 0.0;
        const Real br  = bx*qh[0] + by*qh[1] + bz*qh[2];
        const Real bxi = bx*e1[0] + by*e1[1] + bz*e1[2];
        const Real bet = bx*e2[0] + by*e2[1] + bz*e2[2];
        bcc0(m,IBX,k,j,i) = br;
        bcc0(m,IBY,k,j,i) = bxi;
        bcc0(m,IBZ,k,j,i) = (bet - ccell(m,k,j)*bxi)/scell(m,k,j);
      });
    }

    // Convert primitives to conserved
    auto &u0 = pmbp->pmhd->u0;
    if (!pmbp->pcoord->is_dynamical_relativistic) {
      pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
    }
  } // End initialization of MHD variables

  // Initialize ADM variables -------------------------------
  if (pmbp->padm != nullptr) {
    // Assume Minkowski space for now
    bool schwarzschild = pin->GetOrAddBoolean("problem", "schwarzschild", false);
    if (schwarzschild) {
      pmbp->padm->SetADMVariables = &SetADMVariablesToSchwarzschild;
    }
    pmbp->padm->SetADMVariables(pmbp);

    // If we're using the ADM variables, then we've got dynamic GR enabled.
    // Because we need the metric, we can't initialize the conserved variables
    // until we've filled out the ADM variables.
    pmbp->pdyngr->PrimToConInit(is, ie, js, je, ks, ke);
  }


  return;
}

namespace {

void SetADMVariablesToSchwarzschild(MeshBlockPack *pmbp) {
  auto &adm = pmbp->padm->adm;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pmbp->pmesh->mb_indcs;
  int ivx = shk_dir;
  int ivy = IVX + ((ivx - IVX) + 1)%3;
  int ivz = IVX + ((ivx - IVX) + 2)%3;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nmb1 = pmbp->nmb_thispack - 1;
  int ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  auto shk_dir_ = shk_dir;
  par_for("pgen_adm_vars", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    // Set ADM to Schwarzschild coordinates
    Real r = 0;
    if (shk_dir_ == 1) {
      r = x1v;
    } else if (shk_dir_ == 2) {
      r = x2v;
    } else if (shk_dir_ == 3) {
      r = x3v;
    }
    Real fac = 1. - 2./r;

    adm.alpha(m, k, j, i) = sqrt(fac);
    adm.beta_u(m, 0, k, j, i) = 0.0;
    adm.beta_u(m, 1, k, j, i) = 0.0;
    adm.beta_u(m, 2, k, j, i) = 0.0;

    adm.psi4(m, k, j, i) = 1.0/sqrt(adm.alpha(m, k, j, i));

    adm.g_dd(m, ivx-IVX, ivx-IVX, k, j, i) = 1.0/fac;
    adm.g_dd(m, ivy-IVX, ivy-IVX, k, j, i) = r*r;
    adm.g_dd(m, ivz-IVX, ivz-IVX, k, j, i) = r*r;
    adm.g_dd(m, 0, 1, k, j, i) = 0.0;
    adm.g_dd(m, 0, 2, k, j, i) = 0.0;
    adm.g_dd(m, 1, 2, k, j, i) = 0.0;

    adm.vK_dd(m, 0, 0, k, j, i) = 0.0;
    adm.vK_dd(m, 0, 1, k, j, i) = 0.0;
    adm.vK_dd(m, 0, 2, k, j, i) = 0.0;
    adm.vK_dd(m, 1, 1, k, j, i) = 0.0;
    adm.vK_dd(m, 1, 2, k, j, i) = 0.0;
    adm.vK_dd(m, 2, 2, k, j, i) = 0.0;
  });
}

} // namespace
