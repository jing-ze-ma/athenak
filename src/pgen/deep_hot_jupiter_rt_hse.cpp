//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter_rt_hse.cpp
//! \brief deep_hot_jupiter_rt problem generator (tusplit-1102 piece):
//! HydrostaticEquilibrium, the boundary/HSE hook.  Moved verbatim from deep_hot_jupiter_rt.cpp.

#include "deep_hot_jupiter_rt.hpp"

//----------------------------------------------------------------------------------------
//! \fn HseMaxwellBlockend
//! \brief seamexact-1103: top_bc = hse with bc_outer_maxwell and cs_seam_blockend =
//! extend.  The stress divergence of a top ghost cell reads the cell-centred field one
//! column further out (j+1, k+1); in the OUTERMOST tangential ghost row that column is
//! not in the block and the difference is one-sided -- next to an interior MeshBlock end
//! the same physical column is an active one, differenced centrally, when the panel side
//! is one block.  The row's PLM state enters the x2/x3 flux at the block-boundary face
//! of the first radial ghost and through it the top edge EMF (~1e-7 in B between
//! layouts).  Each block evaluates the divergence on its ACTIVE columns (one value per
//! top ghost layer, stored at i = ie), a same-panel ghost exchange (plain scalars) hands
//! the owning neighbour's values to the outermost row, and the boundary kernel uses them
//! there (only at an interior high-side block end, not in the corner of two).

namespace {
MeshBoundaryValuesCC *hse_mx_pbv_ = nullptr;
}  // namespace

template <typename F>
void HseMaxwellBlockend(MeshBlockPack *pmbp, const F &mxdiv, DvceArray5D<Real> &q,
                        DvceArray1D<int> &hj, DvceArray1D<int> &hk) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
  const int n1 = indcs.nx1 + 2*ng, n2 = indcs.nx2 + 2*ng, n3 = indcs.nx3 + 2*ng;
  const int nmb = pmbp->nmb_thispack;
  q = DvceArray5D<Real>("hse_mxq", nmb, ng, n3, n2, n1);
  hj = DvceArray1D<int>("hse_mxhj", nmb);
  hk = DvceArray1D<int>("hse_mxhk", nmb);
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto qq = q;
  auto hjj = hj;
  auto hkk = hk;
  par_for("hse_mx_end", DevExeSpace(), 0, nmb-1, KOKKOS_LAMBDA(const int m) {
    hjj(m) = (size.d_view(m).x2max < 1.0 - 1.0e-12) ? 1 : 0;   // interior high j end
    hkk(m) = (size.d_view(m).x3max < 1.0 - 1.0e-12) ? 1 : 0;
  });
  par_for("hse_mx_own", DevExeSpace(), 0, nmb-1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::user) return;
    for (int i=0; i<ng; ++i) qq(m,i,k,j,ie) = mxdiv(m,k,j,i);
  });
  if (hse_mx_pbv_ == nullptr) {
    hse_mx_pbv_ = new MeshBoundaryValuesCC(pmbp, nullptr, false);
    hse_mx_pbv_->InitializeBuffers(ng);
    hse_mx_pbv_->SetVectorPairs(ng, {});   // plain scalars; only same-panel ghosts are read
    hse_mx_pbv_->cs_noresample = true;
  }
  while (hse_mx_pbv_->InitRecv(ng) == TaskStatus::incomplete) {}
  while (hse_mx_pbv_->PackAndSendCC(q, q) == TaskStatus::incomplete) {}
  while (hse_mx_pbv_->RecvAndUnpackCC(q, q) == TaskStatus::incomplete) {}
  while (hse_mx_pbv_->ClearSend() == TaskStatus::incomplete) {}
}

void HydrostaticEquilibrium(Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;  int &ie  = indcs.ie;
  int &js = indcs.js;  int &je  = indcs.je;
  int &ks = indcs.ks;  int &ke  = indcs.ke;
  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
  int nmb = pm->pmb_pack->nmb_thispack;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &size = pmbp->pmb->mb_size;

    // THE BOUNDARY READS THE CONSERVED STATE ONLY -- u0 and the face fields b0, never
    // w0, bcc0 or wtemp.  Those three are NOT restart state: restart.cpp writes u0 and
    // b0 (the full nx+2*ng arrays, ghosts included) and nothing else, and
    // Driver::InitBoundaryValuesAndPrimitives runs RecvU -> ApplyPhysicalBCs ->
    // Prolongate -> ConToPrim, i.e. the user boundary is called BEFORE the first
    // ConToPrim of a restarted run, when w0/bcc0/wtemp are still exactly zero.  Anything
    // read from them there is not the value the running boundary sees, so the restart is
    // not a continuation of the run that wrote the file.  Reading u0/b0 also drops the
    // one-stage lag: ApplyPhysicalBCs runs after RKUpdate, the source terms and the RT,
    // and w0 predates all of them.
    DvceArray5D<Real> u0_;
    Real gamma;
    bool use_etotgrav = false;
    bool use_wellbalance_dynamic = false;
    
    DvceArray4D<Real> phi0_x1f;
    DvceArray4D<Real> phicc0;
    DvceArray5D<Real> bcc0;
    DvceArray4D<Real> b0_x1f;
    DvceArray4D<Real> b0_x2f;
    DvceArray4D<Real> b0_x3f;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x2f_ = pmbp->pcoord->xx2f;
    auto &x3v_ = pmbp->pcoord->x3v;
    auto &x3f_ = pmbp->pcoord->xx3f;
    auto &area1 = pmbp->pcoord->area.x1f;
    auto &area2 = pmbp->pcoord->area.x2f;
    auto &area3 = pmbp->pcoord->area.x3f;
    auto &volume = pmbp->pcoord->volume;
    auto &z_ov_rE = pmbp->pcoord->z_ov_rE;
    Real grav_acc = -pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;

    EOS_Data eos;
    if (pmbp->phydro != nullptr) {
      u0_ = pmbp->phydro->u0;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
      use_etotgrav = pmbp->phydro->use_etotgrav;
      use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
      phi0_x1f = pmbp->phydro->phi0.x1f;
      phicc0 = pmbp->phydro->phicc0;
    } else if (pmbp->pmhd != nullptr) {
      u0_ = pmbp->pmhd->u0;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
      use_etotgrav = pmbp->pmhd->use_etotgrav;
      use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
      phi0_x1f = pmbp->pmhd->phi0.x1f;
      phicc0 = pmbp->pmhd->phicc0;
      bcc0 = pmbp->pmhd->bcc0;
      b0_x1f = pmbp->pmhd->b0.x1f;
      b0_x2f = pmbp->pmhd->b0.x2f;
      b0_x3f = pmbp->pmhd->b0.x3f;
    }
    Real bbot = pm->pgen->hot_jupiter_param.bbot;
    
    int nvar = u0_.extent_int(1);
    
    Real igm1 = 1.0/(gamma-1.0);
    Real gigm1 = gamma*igm1;
    Real gm1ig = (gamma-1.0)/gamma;
    Real ig = 1.0/gamma;
    
    // Local copy of the file-scope flag. Reading the global directly from the kernel is
    // a reference to a __host__ variable in device code, which hipcc rejects outright --
    // the switch has to be captured by value like any other host state.
    const bool bc_outer_maxwell_ = bc_outer_maxwell;
    const bool top_out_ = (dhj_top_bc != 0);   // problem/top_bc = outflow1 | outflow0
    const bool top_ss_ = top_out_ && (dhj_top_bt == 1);   // mhd/top_bc_bt sourcesurface
    // for the metric-correct ghost kinetic energy below (cubed sphere only)
    const bool cs_bc_ = pm->use_cubed_sphere;
    auto &ccell_bc_ = pmbp->pcoord->cos_cell;
    auto &scell_bc_ = pmbp->pcoord->sin_cell;
    // whether MHD is on, as a VALUE: `pmbp->pmhd != nullptr` inside a device lambda is a
    // host-pointer dereference (see the note in CLAUDE.md on capturing `this`).
    const bool is_mhd_ = (pmbp->pmhd != nullptr);

    // THE INNER-x1 GHOST MAGNETIC ENERGY IS REMOVED HERE, BEFORE THE FIELD IS REFLECTED.
    // The inner ghosts keep the internal energy and the momenta the problem generator
    // gave them; only their MAGNETIC energy follows the reflected field, so every call
    // swaps one for the other:  E -= e_mag(old ghost field);  reflect;  E += e_mag(new).
    // The subtraction used to read bcc0, which is NOT restart state -- it is rebuilt by
    // ConToPrim, and the user boundary is called BEFORE the first ConToPrim of a
    // restarted run (Driver::InitBoundaryValuesAndPrimitives), where bcc0 is still
    // exactly zero.  A restart therefore subtracted NOTHING and added the full field
    // energy, so every inner radial ghost carried a spurious e_mag that the run which
    // wrote the file never had.
    // The face field b0 IS in the restart file, ghosts included (restart.cpp writes the
    // full nx+2*ng arrays), so taking e_mag from the faces -- with the same averaging and
    // the same cubed-sphere frame rotation the rebuild below uses, which is what
    // ConToPrim built bcc0 from in the first place -- is the same number in a running
    // stage and the right one after a restart.  It is also idempotent: a second call in
    // a row subtracts exactly what the first added.
    if (pmbp->pmhd != nullptr) {
      par_for("usrboundaryx1_demag_inner", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
      KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          for (int i=0; i<ng; ++i) {
            const int ig = is-i-1;
            Real lw, rw;
            lw = (x1f_(m,ig+1)-x1v_(m,ig))/(x1f_(m,ig+1)-x1f_(m,ig));
            rw = (x1v_(m,ig)-x1f_(m,ig))/(x1f_(m,ig+1)-x1f_(m,ig));
            Real bxg = lw*b0_x1f(m,k,j,ig) + rw*b0_x1f(m,k,j,ig+1);
            lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            Real byg = lw*b0_x2f(m,k,j,ig) + rw*b0_x2f(m,k,j+1,ig);
            lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            Real bzg = lw*b0_x3f(m,k,j,ig) + rw*b0_x3f(m,k+1,j,ig);
            if (cs_bc_) {
              const Real cg = ccell_bc_(m,k,j), sg = scell_bc_(m,k,j);
              byg = (byg + cg*bzg)/sg;
            }
            u0_(m,IEN,k,j,ig) -= 0.5*(SQR(bxg)+SQR(byg)+SQR(bzg));
          }
        }
      });
    }

    if (pmbp->pmhd != nullptr) {
      // stress-free wall on a Cartesian mesh (wallsf-1028): inner ghost B_t odd with unit
      // ratios (the r^2 / area weights below are spherical, and a Cartesian mesh has no
      // area arrays; legacy stays bitwise)
      const bool cart_sf_ = (pmbp->pmhd->wall_stress_free ||
                             pmbp->pmhd->wall_conducting_ix1) &&
                            !pm->use_spherical_polar && !pm->use_cubed_sphere;
      par_for("usrboundaryx1_bfield", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
      KOKKOS_LAMBDA(int m, int k, int j) {
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              int iin = is+i;
              int iex = is-i-1;
//              Real fac2 = -area2(m,k,j,iin)/area2(m,k,j,iex);
//              Real fac2p = -area2(m,k,j+1,iin)/area2(m,k,j+1,iex);
              Real fac2 = -(SQR(x1f_(m,iin+1))-SQR(x1f_(m,iin)))/(SQR(x1f_(m,iex+1))-SQR(x1f_(m,iex)));
              Real fac2p = fac2;
              // (a Cartesian mesh allocates no area arrays: unit ratios there)
              Real fac3 = cart_sf_ ? -1.0 : -area3(m,k,j,iin)/area3(m,k,j,iex);
              Real fac3p = cart_sf_ ? -1.0 : -area3(m,k+1,j,iin)/area3(m,k+1,j,iex);
              if (cart_sf_) {  // r^2 weights mean nothing on a Cartesian x1: B_t odd
                fac2 = -1.0;
                fac2p = -1.0;
              }
              b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,iin)*fac2;
              if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,iin)*fac2p;}
              b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,iin)*fac3;
              if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,iin)*fac3p;}
              Real fac1 = cart_sf_ ? 1.0 : area1(m,k,j,iin+1)/area1(m,k,j,iex);
              b0_x1f(m,k,j,iex) = b0_x1f(m,k,j,iin+1)*fac1;
            }
//            for (int i=0; i<ng; ++i) {
//              int iex = is-i-1;
//              Real div_rest = b0_x1f(m,k,j,iex+1)*area1(m,k,j,iex+1) + (b0_x2f(m,k,j+1,iex)*area2(m,k,j+1,iex)-b0_x2f(m,k,j,iex)*area2(m,k,j,iex)) + (b0_x3f(m,k+1,j,iex)*area3(m,k+1,j,iex)-b0_x3f(m,k,j,iex)*area3(m,k,j,iex));
//              b0_x1f(m,k,j,iex) = div_rest/area1(m,k,j,iex);
//            }
          }
          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user && top_out_) {
            // problem/top_bc = outflow1 | outflow0: TANGENTIAL face field zero-gradient
            // (the last active cell's faces copied into every ghost); the ghost NORMAL
            // faces follow from div B = 0 in usrboundaryx1_bfield_divb below
            // ssbc-1111, mhd/top_bc_bt = sourcesurface: the tangential faces are ODD
            // about the top face instead (ghost ie+1+i = -cell ie-i), so B_t = 0 on the
            // face to second order and the field leaving the top is radial (potential-
            // field source surface).  The face EMF the Riemann solver builds from the
            // interior B_t and the opposite ghost B_t carries the wound flux out; B_r
            // follows from div B = 0 as for zerograd.  Unit ratios: the face-normal
            // projections live on the same angular faces at every r, and the
            // r-weighting is O(dr/r) at the face.
            if (top_ss_) {
              for (int i=0; i<ng; ++i) {
                const int iex = ie+i+1;
                const int iin = ie-i;
                b0_x2f(m,k,j,iex) = -b0_x2f(m,k,j,iin);
                if (j == n2-1) {b0_x2f(m,k,j+1,iex) = -b0_x2f(m,k,j+1,iin);}
                b0_x3f(m,k,j,iex) = -b0_x3f(m,k,j,iin);
                if (k == n3-1) {b0_x3f(m,k+1,j,iex) = -b0_x3f(m,k+1,j,iin);}
              }
            } else {
              for (int i=0; i<ng; ++i) {
                int iex = ie+i+1;
                b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,ie);
                if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,ie);}
                b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,ie);
                if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,ie);}
              }
            }
          } else if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              int iin = ie-i;
              int iex = ie+i+1;
//              Real fac2 = -area2(m,k,j,iin)/area2(m,k,j,iex);
//              Real fac2p = -area2(m,k,j+1,iin)/area2(m,k,j+1,iex);
              Real fac2 = -(SQR(x1f_(m,iin+1))-SQR(x1f_(m,iin)))/(SQR(x1f_(m,iex+1))-SQR(x1f_(m,iex)));
              Real fac2p = fac2;
              Real fac3 = -area3(m,k,j,iin)/area3(m,k,j,iex);
              Real fac3p = -area3(m,k+1,j,iin)/area3(m,k+1,j,iex);
              b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,iin)*fac2;
              if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,iin)*fac2p;}
              b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,iin)*fac3;
              if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,iin)*fac3p;}
              Real fac1 = area1(m,k,j,iin)/area1(m,k,j,iex+1);
              b0_x1f(m,k,j,iex+1) = b0_x1f(m,k,j,iin)*fac1;
            }
//            for (int i=0; i<ng; ++i) {
//              int iex = ie+i+1;
//              Real div_rest = b0_x1f(m,k,j,iex)*area1(m,k,j,iex) - (b0_x2f(m,k,j+1,iex)*area2(m,k,j+1,iex)-b0_x2f(m,k,j,iex)*area2(m,k,j,iex)) - (b0_x3f(m,k+1,j,iex)*area3(m,k+1,j,iex)-b0_x3f(m,k,j,iex)*area3(m,k,j,iex));
//              b0_x1f(m,k,j,iex+1) = div_rest/area1(m,k,j,iex+1);
//            }
//              for (int i=0; i<ng; ++i) {
//                b0_x1f(m,k,j,ie+i+2) = -b0_x1f(m,k,j,ie-i);
//                b0_x2f(m,k,j,ie+i+1) =  b0_x2f(m,k,j,ie-i);
//                if (j == n2-1) {b0_x2f(m,k,j+1,ie+i+1) = b0_x2f(m,k,j+1,ie-i);}
//                b0_x3f(m,k,j,ie+i+1) =  b0_x3f(m,k,j,ie-i);
//                if (k == n3-1) {b0_x3f(m,k+1,j,ie+i+1) = b0_x3f(m,k+1,j,ie-i);}
//              }
        }
      });
      // problem/top_bc = outflow1 | outflow0: the ghost NORMAL faces from div B = 0,
      // ghost cell by ghost cell outward (each needs the face below it), with the
      // tangential faces written by the kernel above -- a separate launch, because the
      // j+1 / k+1 faces of a cell belong to other threads of that kernel
      if (top_out_) {
        par_for("usrboundaryx1_bfield_divb", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
        KOKKOS_LAMBDA(int m, int k, int j) {
          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              const int iex = ie+i+1;
              const Real flx = b0_x1f(m,k,j,iex)*area1(m,k,j,iex)
                  - (b0_x2f(m,k,j+1,iex)*area2(m,k,j+1,iex)
                     - b0_x2f(m,k,j,iex)*area2(m,k,j,iex))
                  - (b0_x3f(m,k+1,j,iex)*area3(m,k+1,j,iex)
                     - b0_x3f(m,k,j,iex)*area3(m,k,j,iex));
              b0_x1f(m,k,j,iex+1) = flx/area1(m,k,j,iex+1);
            }
          }
        });
      }
    }
    
//    if (pmbp->phydro != nullptr) {
//      if (use_etotgrav) {
//        pmbp->phydro->RemoveGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//      pmbp->phydro->peos->ConsToPrim(u0_, w0_, false, 0, is-1, 0, (n2-1), 0, (n3-1));
//      if (use_etotgrav) {
//        pmbp->phydro->AddGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//    }
//    else if (pmbp->pmhd != nullptr) {
//      auto b0 = pmbp->pmhd->b0;
//      if (use_etotgrav) {
//        pmbp->pmhd->RemoveGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//      pmbp->pmhd->peos->ConsToPrim(u0_, b0, w0_, bcc0, false, 0, is-1, 0, (n2-1), 0, (n3-1));
//      if (use_etotgrav) {
//        pmbp->pmhd->AddGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//    }
    // (The ConsToPrim on the ie column that used to stand here is gone.  It existed only
    // to refresh w0/wtemp/bcc0 for the outer ghost extrapolation below, which now takes
    // rho, e and T from u0 and b0 itself.  It could not stay: it is called on an ACTIVE
    // column, so its floors and -- under etotgrav -- its RemoveGravEtot/AddGravEtot round
    // trip edited u0(IEN,ie) by a rounding error on EVERY boundary call, and a restarted
    // run makes one extra call (Driver::InitBoundaryValuesAndPrimitives) that the run
    // which wrote the file did not, so the two could not agree bit for bit.  Flooring an
    // active cell is the job of the ConToPrim that follows this boundary in the same task
    // list, and it still does it.)

//    par_for("usrboundaryx1", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n2-1),
//    KOKKOS_LAMBDA(int m, int n, int k, int j) {
//        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
//          for (int i=0; i<ng; ++i) {
//            u0_(m,n,k,j,ie+i+1) = u0_(m,n,k,j,ie);
//          }
//        }
//    });
    
    // The cell-centred field in the OUTER-x1 ghost zones is built here, in its own kernel,
    // and not in the extrapolation kernel below.  It used to be computed inline there, in
    // the same launch that reads bcc0 at (k,j+1) and (k+1,j) for the Maxwell stress -- cells
    // owned by OTHER threads of that launch.  Whether a thread saw its neighbour's new value
    // or the previous cycle's stale one depended on how the wavefronts happened to be
    // scheduled, which made the whole run non-reproducible on a GPU: two runs of the same
    // binary on one rank diverged within a single cycle, worst in the outer half of the
    // domain, where the two-stream RT then spread it along each radial column.  On a CPU the
    // ascending j loop always lost the race the same way, so this was deterministic (and
    // deterministically wrong) until the problem first ran on an accelerator.
    //
    // Splitting the write from the read is the whole fix: every bcc0 the stress term reads
    // is now a completed write from a previous kernel.
    //
    // The INNER-x1 ghosts are deliberately NOT hoisted: that branch subtracts the magnetic
    // energy using the OLD bcc0 and adds it back with the new one, so moving the update
    // ahead of it would change what it subtracts.  It has no cross-thread read and was
    // measured deterministic on its own.
    par_for("usrboundaryx1_bcc_outer", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user &&
            is_mhd_) {
          for (int i=0; i<ng; ++i) {
            Real lw, rw;
            lw = (x1f_(m,(ie+i+1)+1)-x1v_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
            rw = (x1v_(m,(ie+i+1))-x1f_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
            bcc0(m,IBX,k,j,(ie+i+1)) = lw*b0_x1f(m,k,j,(ie+i+1))
                                     + rw*b0_x1f(m,k,j,(ie+i+1)+1);
            lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            bcc0(m,IBY,k,j,(ie+i+1)) = lw*b0_x2f(m,k,j,(ie+i+1))
                                     + rw*b0_x2f(m,k,j+1,(ie+i+1));
            // (the cubed-sphere frame rotation is applied after IBZ below)
            lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            bcc0(m,IBZ,k,j,(ie+i+1)) = lw*b0_x3f(m,k,j,(ie+i+1))
                                     + rw*b0_x3f(m,k+1,j,(ie+i+1));
            if (cs_bc_) {
              // CUBED SPHERE: bcc is stored in the ORTHONORMAL frame that
              // GnomonicEquiangleRaiseVelMHD builds (B.e_xi = (b_xi + c b_eta)/s, eta slot
              // unchanged), and ConsToPrim subtracts 0.5*|bcc|^2 in that frame.  The
              // plain face averages are components on the non-orthogonal tangent pair,
              // whose sum of squares is NOT |B|^2: the two differ by the cross term
              // (c^2 (b2^2 + b3^2) + 2 c b2 b3)/s^2, of order B_tan^2 at a panel corner.
              const Real cg = ccell_bc_(m,k,j), sg = scell_bc_(m,k,j);
              bcc0(m,IBY,k,j,(ie+i+1)) = (bcc0(m,IBY,k,j,(ie+i+1))
                                          + cg*bcc0(m,IBZ,k,j,(ie+i+1)))/sg;
            }
          }
        }
    });

    // seamexact-1103: the Maxwell-stress divergence of the hse top ghost (i = 0..ng-1),
    // shared by the boundary kernel and the block-end exchange below
    auto mxdiv = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i) {
      Real dM1mag;
      // (k,j+1) and (k+1,j) are read here. They belong to other threads, which is
      // why the write had to move to its own kernel -- see the note on that kernel.
      // At the outermost ghost row j+1 and k+1 leave the cell-centred array (the
      // FACE arrays have the extra slot, the cell-centred one does not), so the
      // index is held at the edge: the correction degenerates to a one-sided
      // difference in that row rather than reading past the end, which is what the
      // previous version did.
      int jp1 = (j+1 < n2) ? (j+1) : j;
      int kp1 = (k+1 < n3) ? (k+1) : k;
      // Same one past the end in x1: on the LAST ghost cell (i = ng-1) the cell
      // centred index ie+i+2 is n1, one beyond the array. The face array b0_x1f has
      // the extra slot and is indexed as before; only bcc0 is held back.
      int ip2 = (ie+i+2 < n1) ? (ie+i+2) : (n1-1);
      Real pb = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
      Real pbp1 = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+2)))+SQR(bcc0(m,IBY,k,j,ip2))+SQR(bcc0(m,IBZ,k,j,ip2)));
      Real M11 = pb - SQR(b0_x1f(m,k,j,(ie+i+1)));
      Real M11p1 = pbp1 - SQR(b0_x1f(m,k,j,(ie+i+2)));
      Real M12 = - b0_x2f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
      Real M12p1 = - b0_x2f(m,k,j+1,(ie+i+1)) * bcc0(m,IBX,k,jp1,(ie+i+1));
      Real M13 = - b0_x3f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
      Real M13p1 = - b0_x3f(m,k+1,j,(ie+i+1)) * bcc0(m,IBX,kp1,j,(ie+i+1));
      dM1mag = -( (M11p1*area1(m,k,j,(ie+i+2))-M11*area1(m,k,j,(ie+i+1))) + (M12p1*area2(m,k,j+1,(ie+i+1))-M12*area2(m,k,j,(ie+i+1))) + (M13p1*area3(m,k+1,j,(ie+i+1))-M13*area3(m,k,j,(ie+i+1))) )/volume(m,k,j,(ie+i+1));
      return dM1mag;
    };
    // seamexact-1103 (cs_seam_blockend = extend): the outermost tangential ghost row's
    // stress divergence from the same-panel neighbour that owns those columns
    const bool mxe_ = is_mhd_ && bc_outer_maxwell_ && !top_out_ &&
                      pmbp->pmhd->pbval_u->cs_blockend_ext_;
    DvceArray5D<Real> mxq;
    DvceArray1D<int> mxhj, mxhk;
    if (mxe_) {
      HseMaxwellBlockend(pmbp, mxdiv, mxq, mxhj, mxhk);
    }
    par_for("usrboundaryx1_bfieldc", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          Real rho_i = u0_(m,IDN,k,j,is);
//          Real e_i = w0_(m,IEN,k,j,is);
//          Real phi_i = phicc0(m,k,j,is);
//          Real q0_i = log(e_i);
//          Real factor_i = rho_i/e_i*igm1;
          for (int i=0; i<ng; ++i) {
            if (is_mhd_) {
              // the OLD ghost magnetic energy was removed in usrboundaryx1_demag_inner,
              // before the reflection overwrote the faces it was built from; here only
              // the NEW one is added back.  See the note on that kernel.
              Real lw, rw;
              lw = (x1f_(m,(is-i-1)+1)-x1v_(m,(is-i-1)))/(x1f_(m,(is-i-1)+1)-x1f_(m,(is-i-1)));
              rw = (x1v_(m,(is-i-1))-x1f_(m,(is-i-1)))/(x1f_(m,(is-i-1)+1)-x1f_(m,(is-i-1)));
              bcc0(m,IBX,k,j,(is-i-1)) = lw*b0_x1f(m,k,j,(is-i-1)) + rw*b0_x1f(m,k,j,(is-i-1)+1);
              lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              bcc0(m,IBY,k,j,(is-i-1)) = lw*b0_x2f(m,k,j,(is-i-1)) + rw*b0_x2f(m,k,j+1,(is-i-1));
              lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              bcc0(m,IBZ,k,j,(is-i-1)) = lw*b0_x3f(m,k,j,(is-i-1)) + rw*b0_x3f(m,k+1,j,(is-i-1));
              if (cs_bc_) {
                // CUBED SPHERE: same orthonormal-frame rotation as the outer ghosts.  THIS
                // ONE MATTERED.  The energy update above subtracts 0.5*|bcc0|^2 with the
                // bcc0 that ConsToPrim left here (orthonormal) and adds it back with the
                // one built here (plain averages): on cs the two differ by the cross term,
                // so the ghost TOTAL energy drifted by that mismatch on every call.  The
                // ghost internal energy, uniform over (j,k) at each depth by construction,
                // had a spread of 1e-3 at rot 6, 8 % at rot 16 and ranged 3e8-2e9 around
                // 1.5e9 by rot 26 of cs_prod_mhd_rot -- bottom cells under a low-pressure
                // ghost drained at 2 km/s into the boundary, the domain lost 1 % of its
                // mass and the run went NaN at rot 32.  Spherical polar has c = 0.
                const Real cg = ccell_bc_(m,k,j), sg = scell_bc_(m,k,j);
                bcc0(m,IBY,k,j,(is-i-1)) = (bcc0(m,IBY,k,j,(is-i-1))
                                            + cg*bcc0(m,IBZ,k,j,(is-i-1)))/sg;
              }
              u0_(m,IEN,k,j,is-i-1) += 0.5*(SQR(bcc0(m,IBX,k,j,is-i-1))+SQR(bcc0(m,IBY,k,j,is-i-1))+SQR(bcc0(m,IBZ,k,j,is-i-1)));
            }
//              Real rho0_ip = u0_(m,IDN,k,j,(is-i-1));
//              u0_(m,IM2,k,j,(is-i-1)) = u0_(m,IM2,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IM3,k,j,(is-i-1)) = u0_(m,IM3,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IM1,k,j,(is-i-1)) = 0.0;//-u0_(m,IM1,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IEN,k,j,(is-i-1)) = w0_(m,IEN,k,j,(is-i-1)) + 0.5*(SQR(u0_(m,IM1,k,j,(is-i-1)))+SQR(u0_(m,IM2,k,j,(is-i-1)))+SQR(u0_(m,IM3,k,j,(is-i-1))))/rho0_ip;
////            Real dphi_i = phicc0(m,k,j,(is-i-1))-phi_i;
////            Real q0_ip = q0_i - factor_i * dphi_i;
////            Real e0_ip = exp(q0_ip);
////            if (e0_ip < 0.0) e0_ip = e_i;
////            Real rho0_ip = e0_ip/e_i*rho_i;
////            u0_(m,IDN,k,j,(is-i-1)) = rho0_ip;
////            u0_(m,IM2,k,j,(is-i-1)) = u0_(m,IM2,k,j,is)/rho_i*rho0_ip;
////            u0_(m,IM3,k,j,(is-i-1)) = u0_(m,IM3,k,j,is)/rho_i*rho0_ip;
////            Real mom = u0_(m,IM1,k,j,is)/rho_i*rho0_ip; // fmax(0.0,u0_(m,IM1,k,j,is)/rho_i*rho0_ip); //
////            u0_(m,IM1,k,j,(is-i-1)) = mom;
////            u0_(m,IEN,k,j,(is-i-1)) = e0_ip + 0.5*(SQR(u0_(m,IM1,k,j,(is-i-1)))+SQR(u0_(m,IM2,k,j,(is-i-1)))+SQR(u0_(m,IM3,k,j,(is-i-1))))/rho0_ip;
//            if (use_etotgrav) u0_(m,IEN,k,j,(is-i-1)) += rho0_ip*phicc0(m,k,j,(is-i-1));
//            if (pmbp->pmhd != nullptr) u0_(m,IEN,k,j,(is-i-1)) +=  0.5*(SQR(bcc0(m,IBX,k,j,(is-i-1)))+SQR(bcc0(m,IBY,k,j,(is-i-1)))+SQR(bcc0(m,IBZ,k,j,(is-i-1))));
          }
        }
        // AN OPEN GHOST NEEDS AN INTERIOR TO CONTINUE.  This boundary runs twice before
        // the problem generator has filled anything, with u0 exactly zero; the
        // hydrostatic continuation of THAT is a division by zero.  It used to produce
        // inf/NaN ghosts that the generator then overwrote; the gate makes it leave them
        // alone instead.  Once u0 has a state the continuation is used, and at t = 0 it
        // reproduces the initial column anyway.
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user &&
            u0_(m,IDN,k,j,ie) > 0.0) {
          // rho, e AND T OF THE LAST ACTIVE CELL, FROM THE CONSERVED STATE.  These three
          // used to come from w0(IDN,ie), w0(IEN,ie) and wtemp(ie), refreshed by a
          // ConsToPrim on this column a few lines above.  None of the three is restart
          // state, and the warm start wtemp gave that inversion is not either: the root
          // find stops at |dz| < logtol, so its answer depends on where it started, and
          // at the first boundary call of a restarted run wtemp is zero while in a
          // running stage it is the previous stage's temperature.  Extracting the state
          // here from u0 and b0 -- both restored bitwise by the restart file -- makes the
          // whole ghost fill a pure function of the restart state, hence idempotent, and
          // drops the one-stage lag w0 carried (ApplyPhysicalBCs runs after RKUpdate, the
          // source terms and the RT).  EintFromCons performs exactly ConsToPrim's
          // extraction, in the same order, so in a smooth cell nothing moves.
          Real rho_i = u0_(m,IDN,k,j,ie);
          if (rho_i < eos.dfloor) rho_i = eos.dfloor;
          Real emag_i = 0.0;
          if (is_mhd_) {
            // the cell-centred field ConsToPrim would have built from these faces: the
            // same positional average and the same cubed-sphere orthonormal rotation the
            // ghost kernels above use (see usrboundaryx1_bcc_outer).
            Real lw, rw;
            lw = (x1f_(m,ie+1)-x1v_(m,ie))/(x1f_(m,ie+1)-x1f_(m,ie));
            rw = (x1v_(m,ie)-x1f_(m,ie))/(x1f_(m,ie+1)-x1f_(m,ie));
            Real bxi = lw*b0_x1f(m,k,j,ie) + rw*b0_x1f(m,k,j,ie+1);
            lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            Real byi = lw*b0_x2f(m,k,j,ie) + rw*b0_x2f(m,k,j+1,ie);
            lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            Real bzi = lw*b0_x3f(m,k,j,ie) + rw*b0_x3f(m,k+1,j,ie);
            if (cs_bc_) {
              const Real cg = ccell_bc_(m,k,j), sg = scell_bc_(m,k,j);
              byi = (byi + cg*bzi)/sg;
            }
            emag_i = 0.5*(SQR(bxi)+SQR(byi)+SQR(bzi));
          }
          Real e_i = EintFromCons(u0_, m, k, j, ie, (cs_bc_ ? ccell_bc_(m,k,j) : 0.0),
                                  cs_bc_, use_etotgrav, phicc0(m,k,j,ie), emag_i);
          // the pressure and temperature floors SingleC2P applies, in its order, so that
          // the extrapolation never starts from a state the EOS cannot represent.  The
          // temperature is solved WITHOUT a warm start -- the only guess that is the same
          // in a restarted run and in the run that wrote the file.
          Real t_i = -1.0, p_i = 0.0, g1_i = 0.0;
          bool tclamp_i = false;
          if (e_i > 0.0) {
            eos.TemperaturePressureGamma1(rho_i, e_i, -1.0, t_i, p_i, g1_i, tclamp_i);
          }
          if (!(e_i > 0.0) || p_i < eos.pfloor) {
            e_i = eos.EnergyFromPressure(rho_i, eos.pfloor, t_i);
          }
          if (t_i < eos.tfloor) {
            e_i = eos.EnergyFromTemperature(rho_i, eos.tfloor);
            t_i = eos.tfloor;
          }
          Real phi_i = phicc0(m,k,j,ie);
          Real q0_i = log(e_i);
          Real factor_i = rho_i/e_i*igm1;
          if (top_out_) {
            // problem/top_bc = outflow1: zero-gradient copy, one-way v_r.  rho_i and e_i
            // are the floored last-cell values above (the same extraction the hse branch
            // starts from); the tangential momenta are copied (same rho), the radial one
            // is clipped at zero so the ghost never pushes gas into the domain.
            const Real m1o = fmax(u0_(m,IM1,k,j,ie), 0.0);
            const Real m2o = u0_(m,IM2,k,j,ie);
            const Real m3o = u0_(m,IM3,k,j,ie);
            Real ke_o;
            if (cs_bc_) {
              const Real cg = ccell_bc_(m,k,j);
              ke_o = 0.5*(m1o*m1o + (m2o*m2o - 2.0*cg*m2o*m3o + m3o*m3o)/(1.0 - cg*cg))
                     /rho_i;
            } else {
              ke_o = 0.5*(SQR(m1o)+SQR(m2o)+SQR(m3o))/rho_i;
            }
            for (int i=0; i<ng; ++i) {
              const int ig = ie+i+1;
              u0_(m,IDN,k,j,ig) = rho_i;
              u0_(m,IM1,k,j,ig) = m1o;
              u0_(m,IM2,k,j,ig) = m2o;
              u0_(m,IM3,k,j,ig) = m3o;
              u0_(m,IEN,k,j,ig) = e_i + ke_o;
              if (use_etotgrav) u0_(m,IEN,k,j,ig) += rho_i*phicc0(m,k,j,ig);
              if (is_mhd_) u0_(m,IEN,k,j,ig) += 0.5*(SQR(bcc0(m,IBX,k,j,ig))
                                                     + SQR(bcc0(m,IBY,k,j,ig))
                                                     + SQR(bcc0(m,IBZ,k,j,ig)));
            }
          }
          for (int i=0; i<ng && !top_out_; ++i) {
            Real dM1mag = 0.0;
            if (is_mhd_) {
              // bcc0 in these ghost cells was built by usrboundaryx1_bcc_outer above.
              if (bc_outer_maxwell_) {
              // seamexact-1103: at the OUTERMOST tangential ghost row next to an
              // INTERIOR MeshBlock end (cs_seam_blockend = extend) the one-sided
              // difference above is layout dependent -- one block per panel side
              // differences centrally there; the owning neighbour's value, exchanged
              // (HseMaxwellBlockend), replaces it
              if (mxe_ && ((j == n2-1 && mxhj(m)) || (k == n3-1 && mxhk(m))) &&
                  !((j == n2-1) && (k == n3-1))) {
                dM1mag = mxq(m,i,k,j,ie);
              } else {
                dM1mag = mxdiv(m,k,j,i);
              }
              dM1mag += z_ov_rE(m,k,j,(ie+i+1)) * 0.5*SQR(bcc0(m,IBX,k,j,(ie+i+1)));
              }
            }
            // Hydrostatic extrapolation into the ghost zone, at fixed temperature. The
            // closed form below is the ideal-gas isothermal background; for a general EOS
            // the same statement is WBAdvance's isothermal branch, which integrates
            // dln d/dPhi = -d/(p chi_rho) instead of assuming p = (gamma-1)e.
            Real dphi_i = phicc0(m,k,j,(ie+i+1))-phi_i;
            // The magnetic force enters as an EFFECTIVE GRAVITY, not as an additive
            // shift of the extrapolated energy. dM1mag/rho is an acceleration, so
            // dphi -> dphi*(1 - dM1mag/(rho |g|)) is the same to first order but stays
            // inside the exponential (and inside WBAdvance's EOS-consistent
            // integration), so a magnetic force comparable to gravity changes the SCALE
            // HEIGHT instead of swinging the answer linearly through zero.
            //
            // The clamp is what makes it usable at high field. gmag -> 1 is the
            // force-free limit where the atmosphere stops falling off, and gmag > 1 is
            // net outward, which an outward-decaying ghost cannot represent; gmag < -1
            // would compress the ghost without bound. Both ends are held back, so the
            // ghost degrades to "very extended" rather than to nonsense.
            if (bc_outer_maxwell_) {
              // local gravity, so the clamp means the same fraction of the real weight
              Real gmag = dM1mag/(rho_i*fabs(GravAccAt(grav_acc, ap,
                                                       x1v_(m,ie+i+1), grav_pmass)));
              gmag = fmin(fmax(gmag, -1.0), 0.9);
              dphi_i *= (1.0 - gmag);
            }
            Real e0_hyd, rho0_hyd;
            if (eos.IsGeneral()) {
              rho0_hyd = rho_i;
              e0_hyd = e_i;
              Real t_hyd = -1.0;   // WBAdvance's temperature hand-off; unused here
              // Warm start from cell ie's own temperature, solved above from u0 (NOT
              // from the wtemp cache, which a restart does not have). This runs per ghost
              // cell per stage, and the isothermal branch inverts twice.
              WBAdvance(eos, 1, rho_i, e_i, dphi_i, rho0_hyd, e0_hyd, t_hyd, t_i);
            } else {
              e0_hyd = exp(q0_i - factor_i * dphi_i);
              rho0_hyd = e0_hyd/e_i*rho_i;
            }
            Real e0_ip = e0_hyd;
            if (e0_ip < 0.0) e0_ip = e_i;
            // density and energy now come from the SAME hydrostatic solve, so they are
            // thermodynamically consistent by construction; the old code had to rescale
            // the density by the energy's relative shift because the magnetic term was
            // bolted on afterwards
            Real rho0_ip = rho0_hyd;
//            rho0_ip = rho_i;
//            e0_ip = e_i;
            u0_(m,IDN,k,j,(ie+i+1)) = rho0_ip;
            u0_(m,IM2,k,j,(ie+i+1)) = u0_(m,IM2,k,j,ie)/rho_i*rho0_ip;
            u0_(m,IM3,k,j,(ie+i+1)) = u0_(m,IM3,k,j,ie)/rho_i*rho0_ip;
            Real mom = u0_(m,IM1,k,j,ie)/rho_i*rho0_ip; // fmax(0.0,u0_(m,IM1,k,j,ie)/rho_i*rho0_ip); // 
            u0_(m,IM1,k,j,(ie+i+1)) = mom;
            // Kinetic energy of the ghost cell.  On the CUBED SPHERE IM2/IM3 are
            // COVARIANT components on a non-orthogonal tangent basis, so
            // KE = 0.5 m_i G^{ij} m_j / rho with G^{-1} = [[1,-c],[-c,1]]/(1-c^2) --
            // the contraction GnomonicEquiangleRaiseVelMHD uses.  The plain sum of
            // squares (kept, bit for bit, for every other grid) drops the -2c m2 m3
            // cross term, up to 50 % of the tangential KE at a panel corner where
            // |c| = 1/2, and ConsToPrim then reads that error back out of E as
            // internal energy.  Same bug class as the rotation source that used to
            // fall through into the beta-plane.
            Real ke_g;
            {
              const Real m1g = u0_(m,IM1,k,j,(ie+i+1));
              const Real m2g = u0_(m,IM2,k,j,(ie+i+1));
              const Real m3g = u0_(m,IM3,k,j,(ie+i+1));
              if (cs_bc_) {
                const Real cg = ccell_bc_(m,k,j);
                ke_g = 0.5*(m1g*m1g + (m2g*m2g - 2.0*cg*m2g*m3g + m3g*m3g)/(1.0 - cg*cg))
                       /rho0_ip;
              } else {
                ke_g = 0.5*(SQR(m1g)+SQR(m2g)+SQR(m3g))/rho0_ip;
              }
            }
            u0_(m,IEN,k,j,(ie+i+1)) = e0_ip + ke_g;
            if (use_etotgrav) u0_(m,IEN,k,j,(ie+i+1)) += rho0_ip*phicc0(m,k,j,(ie+i+1));
            if (is_mhd_) u0_(m,IEN,k,j,(ie+i+1)) +=  0.5*(SQR(bcc0(m,IBX,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
          }
        }
    });

////    par_for("usrboundaryx2", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n1-1),
////    KOKKOS_LAMBDA(int m, int n, int k, int i) {
////        if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
////          for (int j=0; j<ng; ++j) {
////            if (n==(IVY)) {
////              u0_(m,n,k,js-j-1,i) = -u0_(m,n,k,js+j,i);
////            } else {
////              u0_(m,n,k,js-j-1,i) = u0_(m,n,k,js+j,i);
////            }
////          }
////        }
////        if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
////          for (int j=0; j<ng; ++j) {
////            if (n==(IVY)) {
////              u0_(m,n,k,je+j+1,i) = -u0_(m,n,k,je-j,i);
////            } else {
////              u0_(m,n,k,je+j+1,i) = u0_(m,n,k,je-j,i);
////            }
////          }
////        }
////    });
//    par_for("usrboundaryx2", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n1-1),
//    KOKKOS_LAMBDA(int m, int k, int i) {
//        if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//          for (int j=0; j<ng; ++j) {
//            u0_(m,IDN,k,js-j-1,i) = u0_(m,IDN,k,js,i);
//            u0_(m,IM1,k,js-j-1,i) = u0_(m,IM1,k,js,i);
//            u0_(m,IM3,k,js-j-1,i) = u0_(m,IM3,k,js,i);
//            u0_(m,IM2,k,js-j-1,i) = fmin(0.0,u0_(m,IM2,k,js,i));
//            u0_(m,IEN,k,js-j-1,i) = u0_(m,IEN,k,js,i)-0.5*SQR(u0_(m,IM2,k,js,i))+0.5*SQR(u0_(m,IM2,k,js-j-1,i));
//          }
//        }
//        if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//          for (int j=0; j<ng; ++j) {
//            u0_(m,IDN,k,je+j+1,i) = u0_(m,IDN,k,je,i);
//            u0_(m,IM1,k,je+j+1,i) = u0_(m,IM1,k,je,i);
//            u0_(m,IM3,k,je+j+1,i) = u0_(m,IM3,k,je,i);
//            u0_(m,IM2,k,je+j+1,i) = fmax(0.0,u0_(m,IM2,k,je,i));
//            u0_(m,IEN,k,je+j+1,i) = u0_(m,IEN,k,je,i)-0.5*SQR(u0_(m,IM2,k,je,i))+0.5*SQR(u0_(m,IM2,k,je+j+1,i));
//          }
//        }
//    });
//    if (pmbp->pmhd != nullptr) {
//      par_for("usrboundaryx2_bfield", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n1-1),
//      KOKKOS_LAMBDA(int m, int k, int i) {
//          if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
//                int jin = js+j;
//                int jex = js-j-1;
//                Real fac1 = -area1(m,k,jin,i)/area1(m,k,jex,i);
//                Real fac1p = -area1(m,k,jin,i+1)/area1(m,k,jex,i+1);
//                Real fac3 = -area3(m,k,jin,i)/area3(m,k,jex,i);
//                Real fac3p = -area3(m,k+1,jin,i)/area3(m,k+1,jex,i);
//                b0_x1f(m,k,jex,i) = b0_x1f(m,k,jin,i)*fac1;
//                if (i == n1-1) {b0_x1f(m,k,jex,i+1) = b0_x1f(m,k,jin,i+1)*fac1p;}
//                b0_x3f(m,k,jex,i) = b0_x3f(m,k,jin,i)*fac3;
//                if (k == n3-1) {b0_x3f(m,k+1,jex,i) = b0_x3f(m,k+1,jin,i)*fac3p;}
//                Real fac2 = area2(m,k,jin+1,i)/area2(m,k,jex,i);
//                b0_x2f(m,k,jex,i) = b0_x2f(m,k,jin+1,i)*fac2;
//              }
////              for (int j=0; j<ng; ++j) {
////                int jex = js-j-1;
////                Real div_rest = b0_x2f(m,k,jex+1,i)*area2(m,k,jex+1,i) + (b0_x1f(m,k,jex,i+1)*area1(m,k,jex,i+1)-b0_x1f(m,k,jex,i)*area1(m,k,jex,i)) + (b0_x3f(m,k+1,jex,i)*area3(m,k+1,jex,i)-b0_x3f(m,k,jex,i)*area3(m,k,jex,i));
////                b0_x2f(m,k,jex,i) = div_rest/area2(m,k,jex,i);
////              }
//          }
//          if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
//                int jin = je-j;
//                int jex = je+j+1;
//                Real fac1 = -area1(m,k,jin,i)/area1(m,k,jex,i);
//                Real fac1p = -area1(m,k,jin,i+1)/area1(m,k,jex,i+1);
//                Real fac3 = -area3(m,k,jin,i)/area3(m,k,jex,i);
//                Real fac3p = -area3(m,k+1,jin,i)/area3(m,k+1,jex,i);
//                b0_x1f(m,k,jex,i) = b0_x1f(m,k,jin,i)*fac1;
//                if (i == n1-1) {b0_x1f(m,k,jex,i+1) = b0_x1f(m,k,jin,i+1)*fac1p;}
//                b0_x3f(m,k,jex,i) = b0_x3f(m,k,jin,i)*fac3;
//                if (k == n3-1) {b0_x3f(m,k+1,jex,i) = b0_x3f(m,k+1,jin,i)*fac3p;}
//                Real fac2 = area2(m,k,jin,i)/area2(m,k,jex+1,i);
//                b0_x2f(m,k,jex+1,i) = b0_x2f(m,k,jin,i)*fac2;
//              }
////              for (int j=0; j<ng; ++j) {
////                int jex = je+j+1;
////                Real div_rest = b0_x2f(m,k,jex,i)*area2(m,k,jex,i) - (b0_x1f(m,k,jex,i+1)*area1(m,k,jex,i+1)-b0_x1f(m,k,jex,i)*area1(m,k,jex,i)) - (b0_x3f(m,k+1,jex,i)*area3(m,k+1,jex,i)-b0_x3f(m,k,jex,i)*area3(m,k,jex,i));
////                b0_x2f(m,k,jex+1,i) = div_rest/area2(m,k,jex+1,i);
////              }
//          }
//      });
//        par_for("usrboundaryx2_bfieldc", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n1-1),
//        KOKKOS_LAMBDA(int m, int k, int i) {
//            if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//                Real lw, rw;
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,js,i) = lw*b0_x1f(m,k,js,i) + rw*b0_x1f(m,k,js,i+1);
//                lw = (x2f_(m,js+1)-x2v_(m,js))/(x2f_(m,js+1)-x2f_(m,js));
//                rw = (x2v_(m,js)-x2f_(m,js))/(x2f_(m,js+1)-x2f_(m,js));
//                bcc0(m,IBY,k,js,i) = lw*b0_x2f(m,k,js,i) + rw*b0_x2f(m,k,js+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,js,i) = lw*b0_x3f(m,k,js,i) + rw*b0_x3f(m,k+1,js,i);
//
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,je,i) = lw*b0_x1f(m,k,je,i) + rw*b0_x1f(m,k,je,i+1);
//                lw = (x2f_(m,je+1)-x2v_(m,je))/(x2f_(m,je+1)-x2f_(m,je));
//                rw = (x2v_(m,je)-x2f_(m,je))/(x2f_(m,je+1)-x2f_(m,je));
//                bcc0(m,IBY,k,je,i) = lw*b0_x2f(m,k,je,i) + rw*b0_x2f(m,k,je+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,je,i) = lw*b0_x3f(m,k,je,i) + rw*b0_x3f(m,k+1,je,i);
//
//              for (int j=0; j<ng; ++j) {
////                Real lw, rw;
////                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                bcc0(m,IBX,k,js+j,i) = lw*b0_x1f(m,k,js+j,i) + rw*b0_x1f(m,k,js+j,i+1);
////                lw = (x2f_(m,js+j+1)-x2v_(m,js+j))/(x2f_(m,js+j+1)-x2f_(m,js+j));
////                rw = (x2v_(m,js+j)-x2f_(m,js+j))/(x2f_(m,js+j+1)-x2f_(m,js+j));
////                bcc0(m,IBY,k,js+j,i) = lw*b0_x2f(m,k,js+j,i) + rw*b0_x2f(m,k,js+j+1,i);
////                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                bcc0(m,IBZ,k,js+j,i) = lw*b0_x3f(m,k,js+j,i) + rw*b0_x3f(m,k+1,js+j,i);
//                u0_(m,IEN,k,js-j-1,i) -= 0.5*(SQR(bcc0(m,IBX,k,js,i))+SQR(bcc0(m,IBY,k,js,i))+SQR(bcc0(m,IBZ,k,js,i)));
//
//                Real lw, rw;
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,js-j-1,i) = lw*b0_x1f(m,k,js-j-1,i) + rw*b0_x1f(m,k,js-j-1,i+1);
//                lw = (x2f_(m,js-j-1+1)-x2v_(m,js-j-1))/(x2f_(m,js-j-1+1)-x2f_(m,js-j-1));
//                rw = (x2v_(m,js-j-1)-x2f_(m,js-j-1))/(x2f_(m,js-j-1+1)-x2f_(m,js-j-1));
//                bcc0(m,IBY,k,js-j-1,i) = lw*b0_x2f(m,k,js-j-1,i) + rw*b0_x2f(m,k,js-j-1+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,js-j-1,i) = lw*b0_x3f(m,k,js-j-1,i) + rw*b0_x3f(m,k+1,js-j-1,i);
//                u0_(m,IEN,k,js-j-1,i) += 0.5*(SQR(bcc0(m,IBX,k,js-j-1,i))+SQR(bcc0(m,IBY,k,js-j-1,i))+SQR(bcc0(m,IBZ,k,js-j-1,i)));
//              }
//            }
//            if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
////                  Real lw, rw;
////                  lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                  rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                  bcc0(m,IBX,k,je-j,i) = lw*b0_x1f(m,k,je-j,i) + rw*b0_x1f(m,k,je-j,i+1);
////                  lw = (x2f_(m,je-j+1)-x2v_(m,je-j))/(x2f_(m,je-j+1)-x2f_(m,je-j));
////                  rw = (x2v_(m,je-j)-x2f_(m,je-j))/(x2f_(m,je-j+1)-x2f_(m,je-j));
////                  bcc0(m,IBY,k,je-j,i) = lw*b0_x2f(m,k,je-j,i) + rw*b0_x2f(m,k,je-j+1,i);
////                  lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                  rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                  bcc0(m,IBZ,k,je-j,i) = lw*b0_x3f(m,k,je-j,i) + rw*b0_x3f(m,k+1,je-j,i);
//                  u0_(m,IEN,k,je+j+1,i) -= 0.5*(SQR(bcc0(m,IBX,k,je,i))+SQR(bcc0(m,IBY,k,je,i))+SQR(bcc0(m,IBZ,k,je,i)));
//
//                  Real lw, rw;
//                  lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                  rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                  bcc0(m,IBX,k,je+j+1,i) = lw*b0_x1f(m,k,je+j+1,i) + rw*b0_x1f(m,k,je+j+1,i+1);
//                  lw = (x2f_(m,je+j+1+1)-x2v_(m,je+j+1))/(x2f_(m,je+j+1+1)-x2f_(m,je+j+1));
//                  rw = (x2v_(m,je+j+1)-x2f_(m,je+j+1))/(x2f_(m,je+j+1+1)-x2f_(m,je+j+1));
//                  bcc0(m,IBY,k,je+j+1,i) = lw*b0_x2f(m,k,je+j+1,i) + rw*b0_x2f(m,k,je+j+1+1,i);
//                  lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                  rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                  bcc0(m,IBZ,k,je+j+1,i) = lw*b0_x3f(m,k,je+j+1,i) + rw*b0_x3f(m,k+1,je+j+1,i);
//                  u0_(m,IEN,k,je+j+1,i) += 0.5*(SQR(bcc0(m,IBX,k,je+j+1,i))+SQR(bcc0(m,IBY,k,je+j+1,i))+SQR(bcc0(m,IBZ,k,je+j+1,i)));
//              }
//            }
//        });
//    }
  return;
}
