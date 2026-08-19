//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_corner_e.cpp
//  \brief
//  Also includes contributions to electric field from "source terms" such as the
//  shearing box.

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "diffusion/resistivity.hpp"
#include "mhd.hpp"

#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn  void MHD::CornerE
//  \brief calculate the corner electric fields.

TaskStatus MHD::CornerE(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &size = pmy_pack->pmb->mb_size;
  auto &flat = pmy_pack->pcoord->coord_data.is_minkowski;
  auto &spin = pmy_pack->pcoord->coord_data.bh_spin;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;

  //---- 1-D problem:
  //  copy face-centered E-fields to edges and return.
  //  Note e2[is:ie+1,js:je,  ks:ke+1]
  //       e3[is:ie+1,js:je+1,ks:ke  ]

  if (pmy_pack->pmesh->one_d) {
    // capture class variables for the kernels
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto e2x1_ = e2x1;
    auto e3x1_ = e3x1;
    par_for("emf1", DevExeSpace(), 0, nmb1, is, ie+1,
    KOKKOS_LAMBDA(int m, int i) {
      e2(m,ks  ,js  ,i) = e2x1_(m,ks,js,i);
      e2(m,ke+1,js  ,i) = e2x1_(m,ks,js,i);
      e3(m,ks  ,js  ,i) = e3x1_(m,ks,js,i);
      e3(m,ks  ,je+1,i) = e3x1_(m,ks,js,i);
    });
  }

  //---- 2-D problem:
  // Copy face-centered E1 and E2 to edges, use GS07 algorithm to compute E3

  if (pmy_pack->pmesh->two_d) {
    // Compute cell-centered E3 = -(v X B) = VyBx-VxBy
    auto w0_ = w0;
    auto bcc_ = bcc0;
    auto e3cc_ = e3_cc;

    // compute cell-centered EMF in dynamical GRMHD
    if (pmy_pack->padm != nullptr) {
      auto &adm = pmy_pack->padm->adm;
      par_for("e_cc_2d", DevExeSpace(), 0, nmb1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int j, int i) {
        // Calculate the spatial components of the three-velocity
        const Real &ux = w0_(m,IVX,ks,j,i);
        const Real &uy = w0_(m,IVY,ks,j,i);
        const Real &uz = w0_(m,IVZ,ks,j,i);
        Real iW = 1.0/sqrt(1.0
                    + adm.g_dd(m,0,0,ks,j,i)*ux*ux + 2.0*adm.g_dd(m,0,1,ks,j,i)*ux*uy
                    + 2.0*adm.g_dd(m,0,2,ks,j,i)*ux*uz + adm.g_dd(m,1,1,ks,j,i)*uy*uy
                    + 2.0*adm.g_dd(m,1,2,ks,j,i)*uy*uz + adm.g_dd(m,2,2,ks,j,i)*uz*uz);
        Real v1 = ux*iW;
        Real v2 = uy*iW;
        //Real v3 = uz*iW;

        const Real &alpha = adm.alpha(m,ks,j,i);
        e3cc_(m,ks,j,i) = bcc_(m,IBX,ks,j,i)*(alpha*v2 - adm.beta_u(m, 1, ks, j, i))
                        - bcc_(m,IBY,ks,j,i)*(alpha*v1 - adm.beta_u(m, 0, ks, j, i));
      });
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      // compute cell-centered EMF in GR MHD
      par_for("e_cc_2d", DevExeSpace(), 0, nmb1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int j, int i) {
        // Extract components of metric
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(0, indcs.nx3, x3min, x3max);

        Real glower[4][4], gupper[4][4];
        ComputeMetricAndInverse(x1v, x2v, x3v, flat, spin, glower, gupper);

        const Real &ux = w0_(m,IVX,ks,j,i);
        const Real &uy = w0_(m,IVY,ks,j,i);
        const Real &uz = w0_(m,IVZ,ks,j,i);
        const Real &bx = bcc_(m,IBX,ks,j,i);
        const Real &by = bcc_(m,IBY,ks,j,i);
        const Real &bz = bcc_(m,IBZ,ks,j,i);
        // Calculate 4-velocity
        Real tmp = glower[1][1]*ux*ux + 2.0*glower[1][2]*ux*uy + 2.0*glower[1][3]*ux*uz
                 + glower[2][2]*uy*uy + 2.0*glower[2][3]*uy*uz
                 + glower[3][3]*uz*uz;
        Real alpha = sqrt(-1.0/gupper[0][0]);
        Real gamma = sqrt(1.0 + tmp);
        Real u0 = gamma / alpha;
        Real u1 = ux - alpha * gamma * gupper[0][1];
        Real u2 = uy - alpha * gamma * gupper[0][2];
        Real u3 = uz - alpha * gamma * gupper[0][3];
        // lower vector indices
        Real u_1 = glower[1][0]*u0 + glower[1][1]*u1 + glower[1][2]*u2 + glower[1][3]*u3;
        Real u_2 = glower[2][0]*u0 + glower[2][1]*u1 + glower[2][2]*u2 + glower[2][3]*u3;
        Real u_3 = glower[3][0]*u0 + glower[3][1]*u1 + glower[3][2]*u2 + glower[3][3]*u3;
        // calculate 4-magnetic field
        Real b0 = u_1*bx + u_2*by + u_3*bz;
        Real b1 = (bx + b0 * u1) / u0;
        Real b2 = (by + b0 * u2) / u0;
        //Real b3 = (bz + b0 * u3) / u0;

        e3cc_(m,ks,j,i) = b1 * u2 - b2 * u1;
      });

    // compute cell-centered EMF in SR MHD
    } else if (pmy_pack->pcoord->is_special_relativistic) {
      par_for("e_cc_2d", DevExeSpace(), 0, nmb1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int j, int i) {
        const Real &u1 = w0_(m,IVX,ks,j,i);
        const Real &u2 = w0_(m,IVY,ks,j,i);
        const Real &u3 = w0_(m,IVZ,ks,j,i);
        Real u0 = sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
        e3cc_(m,ks,j,i) = (u2 * bcc_(m,IBX,ks,j,i) - u1 * bcc_(m,IBY,ks,j,i)) / u0;
      });

    // compute cell-centered EMF in Newtonian MHD
    } else {
      par_for("e_cc_2d", DevExeSpace(), 0, nmb1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int j, int i) {
        e3cc_(m,ks,j,i) = w0_(m,IVY,ks,j,i)*bcc_(m,IBX,ks,j,i) -
                          w0_(m,IVX,ks,j,i)*bcc_(m,IBY,ks,j,i);
      });
    }

    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto e2x1_ = e2x1;
    auto e3x1_ = e3x1;
    auto e1x2_ = e1x2;
    auto e3x2_ = e3x2;
    auto flx1 = uflx.x1f;
    auto flx2 = uflx.x2f;

    // integrate E3 to corner using SG07
    //  Note e1[is:ie,  js:je+1,ks:ke+1]
    //       e2[is:ie+1,js:je,  ks:ke+1]
    //       e3[is:ie+1,js:je+1,ks:ke  ]
    par_for("emf2", DevExeSpace(), 0, nmb1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      e2(m,ks  ,j,i) = e2x1_(m,ks,j,i);
      e2(m,ke+1,j,i) = e2x1_(m,ks,j,i);
      e1(m,ks  ,j,i) = e1x2_(m,ks,j,i);
      e1(m,ke+1,j,i) = e1x2_(m,ks,j,i);

      Real e3_l2, e3_r2, e3_l1, e3_r1;
      if (flx1(m,IDN,ks,j-1,i) >= 0.0) {
        e3_l2 = e3x2_(m,ks,j,i-1) - e3cc_(m,ks,j-1,i-1);
      } else {
        e3_l2 = e3x2_(m,ks,j,i  ) - e3cc_(m,ks,j-1,i  );
      }
      if (flx1(m,IDN,ks,j,i) >= 0.0) {
        e3_r2 = e3x2_(m,ks,j,i-1) - e3cc_(m,ks,j  ,i-1);
      } else {
        e3_r2 = e3x2_(m,ks,j,i  ) - e3cc_(m,ks,j  ,i  );
      }
      if (flx2(m,IDN,ks,j,i-1) >= 0.0) {
        e3_l1 = e3x1_(m,ks,j-1,i) - e3cc_(m,ks,j-1,i-1);
      } else {
        e3_l1 = e3x1_(m,ks,j  ,i) - e3cc_(m,ks,j  ,i-1);
      }
      if (flx2(m,IDN,ks,j,i) >= 0.0) {
        e3_r1 = e3x1_(m,ks,j-1,i) - e3cc_(m,ks,j-1,i  );
      } else {
        e3_r1 = e3x1_(m,ks,j  ,i) - e3cc_(m,ks,j  ,i  );
      }
      e3(m,ks,j,i) = 0.25*(e3_l1 + e3_r1 + e3_l2 + e3_r2 +
             e3x2_(m,ks,j,i-1) + e3x2_(m,ks,j,i) + e3x1_(m,ks,j-1,i) + e3x1_(m,ks,j,i));
    });
  }

  //---- 3-D problem:
  // Use GS07 algorithm to compute all three of E1, E2, and E3

  if (pmy_pack->pmesh->three_d) {
    // Compute cell-centered electric fields
    // E1=-(v X B)=VzBy-VyBz
    // E2=-(v X B)=VxBz-VzBx
    // E3=-(v X B)=VyBx-VxBy
    auto w0_ = w0;
    auto bcc_ = bcc0;
    auto e1cc_ = e1_cc;
    auto e2cc_ = e2_cc;
    auto e3cc_ = e3_cc;

    // compute cell-centered EMFs in dynamical GRMHD
    if (pmy_pack->padm != nullptr) {
      auto &adm = pmy_pack->padm->adm;
      par_for("e_cc_3d", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        // Calculate something that resembles the spatial components of the four-velocity
        // normalized by W.
        const Real &ux = w0_(m,IVX,k,j,i);
        const Real &uy = w0_(m,IVY,k,j,i);
        const Real &uz = w0_(m,IVZ,k,j,i);
        const Real &bx = bcc_(m,IBX,k,j,i);
        const Real &by = bcc_(m,IBY,k,j,i);
        const Real &bz = bcc_(m,IBZ,k,j,i);
        Real iW = 1.0/sqrt(1.0
                         + adm.g_dd(m,0,0,k,j,i)*ux*ux + 2.0*adm.g_dd(m,0,1,k,j,i)*ux*uy
                         + 2.0*adm.g_dd(m,0,2,k,j,i)*ux*uz + adm.g_dd(m,1,1,k,j,i)*uy*uy
                         + 2.0*adm.g_dd(m,1,2,k,j,i)*uy*uz + adm.g_dd(m,2,2,k,j,i)*uz*uz);
        const Real &alpha = adm.alpha(m, k, j, i);
        Real v1c = alpha*ux*iW - adm.beta_u(m, 0, k, j, i);
        Real v2c = alpha*uy*iW - adm.beta_u(m, 1, k, j, i);
        Real v3c = alpha*uz*iW - adm.beta_u(m, 2, k, j, i);

        e1cc_(m,k,j,i) = by * v3c - bz * v2c;
        e2cc_(m,k,j,i) = bz * v1c - bx * v3c;
        e3cc_(m,k,j,i) = bx * v2c - by * v1c;
      });
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      // compute cell-centered EMFs in GR MHD
      par_for("e_cc_3d", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        // Extract components of metric
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

        Real glower[4][4], gupper[4][4];
        ComputeMetricAndInverse(x1v, x2v, x3v, flat, spin, glower, gupper);

        const Real &ux = w0_(m,IVX,k,j,i);
        const Real &uy = w0_(m,IVY,k,j,i);
        const Real &uz = w0_(m,IVZ,k,j,i);
        const Real &bx = bcc_(m,IBX,k,j,i);
        const Real &by = bcc_(m,IBY,k,j,i);
        const Real &bz = bcc_(m,IBZ,k,j,i);
        // Calculate 4-velocity
        Real tmp = glower[1][1]*ux*ux + 2.0*glower[1][2]*ux*uy + 2.0*glower[1][3]*ux*uz
                 + glower[2][2]*uy*uy + 2.0*glower[2][3]*uy*uz
                 + glower[3][3]*uz*uz;
        Real alpha = sqrt(-1.0/gupper[0][0]);
        Real gamma = sqrt(1.0 + tmp);
        Real u0 = gamma / alpha;
        Real u1 = ux - alpha * gamma * gupper[0][1];
        Real u2 = uy - alpha * gamma * gupper[0][2];
        Real u3 = uz - alpha * gamma * gupper[0][3];
        // lower vector indices
        Real u_1 = glower[1][0]*u0 + glower[1][1]*u1 + glower[1][2]*u2 + glower[1][3]*u3;
        Real u_2 = glower[2][0]*u0 + glower[2][1]*u1 + glower[2][2]*u2 + glower[2][3]*u3;
        Real u_3 = glower[3][0]*u0 + glower[3][1]*u1 + glower[3][2]*u2 + glower[3][3]*u3;
        // calculate 4-magnetic field
        Real b0 = u_1*bx + u_2*by + u_3*bz;
        Real b1 = (bx + b0 * u1) / u0;
        Real b2 = (by + b0 * u2) / u0;
        Real b3 = (bz + b0 * u3) / u0;

        e1cc_(m,k,j,i) = b2 * u3 - b3 * u2;
        e2cc_(m,k,j,i) = b3 * u1 - b1 * u3;
        e3cc_(m,k,j,i) = b1 * u2 - b2 * u1;
      });

    // compute cell-centered EMFs in SR MHD
    } else if (pmy_pack->pcoord->is_special_relativistic) {
      par_for("e_cc_3d", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        const Real &u1 = w0_(m,IVX,k,j,i);
        const Real &u2 = w0_(m,IVY,k,j,i);
        const Real &u3 = w0_(m,IVZ,k,j,i);
        Real u0 = sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
        e1cc_(m,k,j,i) = (u3 * bcc_(m,IBY,k,j,i) - u2 * bcc_(m,IBZ,k,j,i)) / u0;
        e2cc_(m,k,j,i) = (u1 * bcc_(m,IBZ,k,j,i) - u3 * bcc_(m,IBX,k,j,i)) / u0;
        e3cc_(m,k,j,i) = (u2 * bcc_(m,IBX,k,j,i) - u1 * bcc_(m,IBY,k,j,i)) / u0;
      });

    // compute cell-centered EMFs in Newtonian MHD
    } else {
      par_for("e_cc_3d", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        e1cc_(m,k,j,i) = w0_(m,IVZ,k,j,i)*bcc_(m,IBY,k,j,i) -
                         w0_(m,IVY,k,j,i)*bcc_(m,IBZ,k,j,i);
        e2cc_(m,k,j,i) = w0_(m,IVX,k,j,i)*bcc_(m,IBZ,k,j,i) -
                         w0_(m,IVZ,k,j,i)*bcc_(m,IBX,k,j,i);
        e3cc_(m,k,j,i) = w0_(m,IVY,k,j,i)*bcc_(m,IBX,k,j,i) -
                         w0_(m,IVX,k,j,i)*bcc_(m,IBY,k,j,i);
      });
    }

    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto e2x1_ = e2x1;
    auto e3x1_ = e3x1;
    auto e1x2_ = e1x2;
    auto e3x2_ = e3x2;
    auto e1x3_ = e1x3;
    auto e2x3_ = e2x3;
    auto flx1 = uflx.x1f;
    auto flx2 = uflx.x2f;
    auto flx3 = uflx.x3f;

    // Integrate E1, E2, E3 to corners
    //  Note e1[is:ie,  js:je+1,ks:ke+1]
    //       e2[is:ie+1,js:je,  ks:ke+1]
    //       e3[is:ie+1,js:je+1,ks:ke  ]
    par_for("emf3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // integrate E1 to corner using SG07
      Real e1_l3, e1_r3, e1_l2, e1_r2;
      if (flx2(m,IDN,k-1,j,i) >= 0.0) {
        e1_l3 = e1x3_(m,k,j-1,i) - e1cc_(m,k-1,j-1,i);
      } else {
        e1_l3 = e1x3_(m,k,j  ,i) - e1cc_(m,k-1,j  ,i);
      }
      if (flx2(m,IDN,k,j,i) >= 0.0) {
        e1_r3 = e1x3_(m,k,j-1,i) - e1cc_(m,k  ,j-1,i);
      } else {
        e1_r3 = e1x3_(m,k,j  ,i) - e1cc_(m,k  ,j  ,i);
      }
      if (flx3(m,IDN,k,j-1,i) >= 0.0) {
        e1_l2 = e1x2_(m,k-1,j,i) - e1cc_(m,k-1,j-1,i);
      } else {
        e1_l2 = e1x2_(m,k  ,j,i) - e1cc_(m,k  ,j-1,i);
      }
      if (flx3(m,IDN,k,j,i) >= 0.0) {
        e1_r2 = e1x2_(m,k-1,j,i) - e1cc_(m,k-1,j  ,i);
      } else {
        e1_r2 = e1x2_(m,k  ,j,i) - e1cc_(m,k  ,j  ,i);
      }
      e1(m,k,j,i) = 0.25*(e1_l3 + e1_r3 + e1_l2 + e1_r2 +
                e1x2_(m,k-1,j,i) + e1x2_(m,k,j,i) + e1x3_(m,k,j-1,i) + e1x3_(m,k,j,i));

      // integrate E2 to corner using SG07
      Real e2_l3, e2_r3, e2_l1, e2_r1;
      if (flx1(m,IDN,k-1,j,i) >= 0.0) {
        e2_l3 = e2x3_(m,k,j,i-1) - e2cc_(m,k-1,j,i-1);
      } else {
        e2_l3 = e2x3_(m,k,j,i  ) - e2cc_(m,k-1,j,i  );
      }
      if (flx1(m,IDN,k,j,i) >= 0.0) {
        e2_r3 = e2x3_(m,k,j,i-1) - e2cc_(m,k  ,j,i-1);
      } else {
        e2_r3 = e2x3_(m,k,j,i  ) - e2cc_(m,k  ,j,i  );
      }
      if (flx3(m,IDN,k,j,i-1) >= 0.0) {
        e2_l1 = e2x1_(m,k-1,j,i) - e2cc_(m,k-1,j,i-1);
      } else {
        e2_l1 = e2x1_(m,k  ,j,i) - e2cc_(m,k  ,j,i-1);
      }
      if (flx3(m,IDN,k,j,i) >= 0.0) {
        e2_r1 = e2x1_(m,k-1,j,i) - e2cc_(m,k-1,j,i  );
      } else {
        e2_r1 = e2x1_(m,k  ,j,i) - e2cc_(m,k  ,j,i  );
      }
      e2(m,k,j,i) = 0.25*(e2_l3 + e2_r3 + e2_l1 + e2_r1 +
                e2x3_(m,k,j,i-1) + e2x3_(m,k,j,i) + e2x1_(m,k-1,j,i) + e2x1_(m,k,j,i));

      // integrate E3 to corner using SG07
      Real e3_l2, e3_r2, e3_l1, e3_r1;
      if (flx1(m,IDN,k,j-1,i) >= 0.0) {
        e3_l2 = e3x2_(m,k,j,i-1) - e3cc_(m,k,j-1,i-1);
      } else {
        e3_l2 = e3x2_(m,k,j,i  ) - e3cc_(m,k,j-1,i  );
      }
      if (flx1(m,IDN,k,j,i) >= 0.0) {
        e3_r2 = e3x2_(m,k,j,i-1) - e3cc_(m,k,j  ,i-1);
      } else {
        e3_r2 = e3x2_(m,k,j,i  ) - e3cc_(m,k,j  ,i  );
      }
      if (flx2(m,IDN,k,j,i-1) >= 0.0) {
        e3_l1 = e3x1_(m,k,j-1,i) - e3cc_(m,k,j-1,i-1);
      } else {
        e3_l1 = e3x1_(m,k,j  ,i) - e3cc_(m,k,j  ,i-1);
      }
      if (flx2(m,IDN,k,j,i) >= 0.0) {
        e3_r1 = e3x1_(m,k,j-1,i) - e3cc_(m,k,j-1,i  );
      } else {
        e3_r1 = e3x1_(m,k,j  ,i) - e3cc_(m,k,j  ,i  );
      }
      e3(m,k,j,i) = 0.25*(e3_l1 + e3_r1 + e3_l2 + e3_r2 +
                e3x2_(m,k,j,i-1) + e3x2_(m,k,j,i) + e3x1_(m,k,j-1,i) + e3x1_(m,k,j,i));
    });
  }
  if (pmy_pack->pmesh->use_polar_boundary) PolarAzimuthalAverageEr();
  return TaskStatus::complete;
}

void MHD::PolarAzimuthalAverageEr() {
  auto e1 = efld.x1e;

  auto &indcs = pmy_pack->pmesh->mb_indcs;

  int nx1 = indcs.nx1;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &ng = indcs.ng;
  const int ncells1 = nx1 + 2*ng;

  const int mesh_nx3 = pmy_pack->pmesh->mesh_indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;

  auto &mb_bcs = pmy_pack->pmb->mb_bcs;

  // Work arrays are members, sized once.  The previous version allocated four Views, a
  // pair of host mirrors and two std::vectors on EVERY call -- twice per cycle.
  if (static_cast<int>(inner_local.extent(0)) != ncells1) {
    Kokkos::realloc(inner_local, ncells1);
    Kokkos::realloc(outer_local, ncells1);
    Kokkos::realloc(polar_inner_h, ncells1);
    Kokkos::realloc(polar_outer_h, ncells1);
  }
  if (static_cast<int>(polar_part_in.extent(0)) != nmb ||
      static_cast<int>(polar_part_in.extent(1)) != ncells1) {
    Kokkos::realloc(polar_part_in,  nmb, ncells1);
    Kokkos::realloc(polar_part_out, nmb, ncells1);
  }
  auto part_in  = polar_part_in;
  auto part_out = polar_part_out;
  auto inner_   = inner_local;
  auto outer_   = outer_local;

  // =====================================================
  // stage 1: per-(MeshBlock, radius) partial sums over phi
  // =====================================================
  // One thread per (m,i) rather than one per i.  The old kernel launched ncells1 threads
  // -- 68 on this problem, well under a single MI300A workgroup -- and made each of them
  // walk every MeshBlock and every phi cell serially.  The sum over phi within one block
  // is left serial and in ascending k, exactly as before.
  par_for("polar_partial", DevExeSpace(), 0, nmb-1, 0, ncells1-1,
  KOKKOS_LAMBDA(const int m, const int i) {
    bool do_inner = (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::polar);
    bool do_outer = (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::polar);

    Real inner_sum = 0.0;
    Real outer_sum = 0.0;
    if ((do_inner || do_outer) && (i >= is) && (i <= ie)) {
      for (int k=ks+1; k<=ke+1; ++k) {
        if (do_inner) inner_sum += e1(m,k,js,i);
        if (do_outer) outer_sum += e1(m,k,je+1,i);
      }
    }
    part_in(m,i)  = inner_sum;
    part_out(m,i) = outer_sum;
  });

  // =====================================================
  // stage 2: combine the per-block partials
  // =====================================================
  // Ascending block order, so the result does not depend on how the runtime scheduled
  // stage 1.  (It is NOT bit-identical to the old single-loop accumulation, which carried
  // one running total across block boundaries -- see the branch notes.)
  par_for("polar_combine", DevExeSpace(), 0, ncells1-1,
  KOKKOS_LAMBDA(const int i) {
    Real inner_sum = 0.0;
    Real outer_sum = 0.0;
    for (int m=0; m<nmb; ++m) {
      inner_sum += part_in(m,i);
      outer_sum += part_out(m,i);
    }
    inner_(i) = inner_sum;
    outer_(i) = outer_sum;
  });

#if MPI_PARALLEL_ENABLED
  // Only a genuinely distributed run needs the device->host->device round trip.  On one
  // rank the old code paid it anyway, twice per cycle, for an Allreduce over itself.
  if (global_variable::nranks > 1) {
    Kokkos::deep_copy(polar_inner_h, inner_local);
    Kokkos::deep_copy(polar_outer_h, outer_local);
    MPI_Allreduce(MPI_IN_PLACE, polar_inner_h.data(), ncells1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, polar_outer_h.data(), ncells1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    Kokkos::deep_copy(inner_local, polar_inner_h);
    Kokkos::deep_copy(outer_local, polar_outer_h);
  }
#endif

  // =====================================================
  // refill polar EMFs with the azimuthal average
  // =====================================================
  // Parallel over (m,k,i); the old kernel ran one thread per MeshBlock and looped k and i
  // inside it.  Pure writes, so this is bit-identical to before.
  par_for("polar_refill", DevExeSpace(), 0, nmb-1, ks, ke+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    bool do_inner = (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::polar);
    bool do_outer = (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::polar);

    if (do_inner) e1(m,k,js  ,i) = inner_(i)/static_cast<Real>(mesh_nx3);
    if (do_outer) e1(m,k,je+1,i) = outer_(i)/static_cast<Real>(mesh_nx3);
  });
}

} // namespace mhd
