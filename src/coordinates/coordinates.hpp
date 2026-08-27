#ifndef COORDINATES_COORDINATES_HPP_
#define COORDINATES_COORDINATES_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file coordinates.hpp
//! \brief implemention of light-weight coordinates class.  Provides data structure that
//! stores array of RegionSizes over (# of MeshBlocks), and inline functions for
//! computing positions.  In GR, also provides inline metric functions (currently only
//! Cartesian Kerr-Schild)

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/grid_stretch.hpp"

// forward declarations
struct EOS_Data;

// Enumerator for the excision method
enum class ExcisionScheme {
  fixed,
  lapse
};

//----------------------------------------------------------------------------------------
//! \struct CoordData
//! \brief container for Coordinate variables and functions needed inside kernels. Storing
//! everything in a container makes them easier to capture, and pass to inline functions,
//! inside kernels.

struct CoordData {
  // following data is only used in GR calculations to compute metric
  bool is_minkowski;               // flag to specify Minkowski (flat) space
  Real bh_spin;                    // needed for GR metric
  bool bh_excise;                  // flag to specify excision
  Real rexcise;                    // excision radius (SKS)
  Real dexcise;                    // rest-mass density inside excised region
  Real pexcise;                    // pressure inside excised region
  Real flux_excise_r;              // reduce to first-order inside this radius
  ExcisionScheme excision_scheme;  // excision method
  Real excise_lapse;               // if excision_scheme = lapse, excise under this lapse
};

//----------------------------------------------------------------------------------------
//! \class Coordinates
//! \brief data and functions for coordinates

class Coordinates {
 public:
  explicit Coordinates(ParameterInput *pin, MeshBlockPack *ppack);
  ~Coordinates() {}

  // flags to denote relativistic dynamics in these coordinates
  bool is_special_relativistic = false;
  bool is_general_relativistic = false;
  bool is_dynamical_relativistic = false;

  // data needed to compute metric in GR
  CoordData coord_data;

  // excision masks
  DvceArray4D<bool> excision_floor;  // cell-centered mask for C2P flooring about horizon
  DvceArray4D<bool> excision_flux;   // cell-centered mask for FOFC about horizon
    
  DvceArray4D<Real> volume;
  DvceFaceFld4D<Real> area;
  DvceArray4D<Real> dx1;
  DvceArray4D<Real> dx2;
  DvceArray4D<Real> dx3;
  DvceArray2D<Real> x1v;
  DvceArray2D<Real> x2v;
  DvceArray2D<Real> x3v;
  DvceArray2D<Real> xx1f;
  DvceArray2D<Real> xx2f;
  DvceArray2D<Real> xx3f;
  DvceEdgeFld4D<Real> dxedge;
  DvceFaceFld4D<Real> dxface;
  DvceEdgeFld4D<Real> areaedge;
  // Gnomonic-equiangle geometry. These depend on the two PANEL-TANGENTIAL angles only --
  // xi = x2 and eta = x3 -- so they are indexed (m,k,j) with no radial extent. The
  // _xi arrays live on x2 faces (staggered in j), the _eta arrays on x3 faces (in k).
  DvceArray3D<Real> sin_cell;
  DvceArray3D<Real> cos_cell;
  DvceArray3D<Real> sin_face_xi;
  DvceArray3D<Real> cos_face_xi;
  DvceArray3D<Real> sin_face_eta;
  DvceArray3D<Real> cos_face_eta;
  DvceArray4D<Real> x_ov_rD;
  DvceArray4D<Real> y_ov_rC;
  DvceArray4D<Real> z_ov_rE;
    
  void CoordSphericalPolar();
  void SrcTermsSphericalPolarHydro(const DvceArray5D<Real> &w0,
      const DvceArray4D<Real> &pwb, const DvceFaceFld5D<Real> uflx,
      const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
  void SrcTermsSphericalPolarMHD(const DvceArray5D<Real> &w0,
      const DvceArray5D<Real> &bcc0, const DvceArray4D<Real> &pwb,
      const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt,
      DvceArray5D<Real> &u0);
  void CoordGnomonicEquiangle();
  // The gnomonic tangent basis is NOT orthogonal, so the conserved momentum (covariant,
  // which is what GnomonicEquiangleFluxX* produces) and the primitive velocity
  // (contravariant, which is what GnomonicEquianglePrimFaceX* and the gnomonic source
  // terms consume) are different objects. These two convert between them with the metric
  // g = [[1, cos],[cos, 1]] on the unit-normalised basis, plus g_33 = 1 radially.
  void GnomonicEquiangleRaiseVel(const DvceArray5D<Real> &u0, DvceArray5D<Real> &w0,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku);
  void GnomonicEquiangleLowerMom(const DvceArray5D<Real> &w0, DvceArray5D<Real> &u0,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku);
  void SrcTermsGnomonicEquiangle(const DvceArray5D<Real> &w0,
      const DvceArray5D<Real> &wder, const DvceFaceFld5D<Real> uflx,
      const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
  void SrcTermsGnomonicEquiangleMHD(const DvceArray5D<Real> &w0,
      const DvceArray5D<Real> &bcc0, const DvceArray5D<Real> &wder,
      const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt,
      DvceArray5D<Real> &u0);
  void SrcTermsGnomonicEquiangleImpl(const DvceArray5D<Real> &w0,
      const DvceArray5D<Real> &bcc0, const bool is_mhd,
      const DvceArray5D<Real> &wder, const DvceFaceFld5D<Real> uflx,
      const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
    
    // GNOMONIC EQUIANGLE, x1 = RADIAL, xi = x2, eta = x3.
    //
    // The two angular basis vectors e_xi, e_eta are unit vectors separated by an angle
    // theta with cos(theta) = cos_cell, so the angular metric is g = [[1,c],[c,1]] and is
    // NOT the identity. The PrimFace* routines rotate the CONTRAVARIANT velocity held in
    // w0 into the locally ORTHONORMAL frame the Riemann solver needs; the Flux* routines
    // invert that rotation on the returned flux and then LOWER the index, so what the
    // update accumulates in u0(IM1,IM2,IM3) is the COVARIANT momentum. The radial
    // component is orthogonal to both angles and is never touched by these rotations.
    // See Coordinates::GnomonicEquiangleRaiseVel for the matching inversion.

    // --- x1 sweep: RADIAL face. Its normal is r-hat, so the normal component IVX passes
    // through untouched and only the angular pair is put in an orthonormal frame. The
    // coefficients depend on (k,j) alone, so they are loop-invariant in i.
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquianglePrimFaceX1(TeamMember_t const &member, const int m, const int k,
         const int j, const int il, const int iu, ScrArray2D<Real> &ql,
         ScrArray2D<Real> &qr) {
        const Real sin_theta = sin_cell(m,k,j);
        const Real cos_theta = cos_cell(m,k,j);
        par_for_inner(member, il, iu, [&](const int i) {
          ql(IVY,i+1) += cos_theta * ql(IVZ,i+1);
          ql(IVZ,i+1) *= sin_theta;
          qr(IVY,i)   += cos_theta * qr(IVZ,i);
          qr(IVZ,i)   *= sin_theta;
        });
      return;
    };

    // --- x2 sweep: XI face, normal component is IVY, tangential is IVZ
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquianglePrimFaceX2(TeamMember_t const &member, const int m, const int k,
         const int j, const int il, const int iu, ScrArray2D<Real> &ql_jp1,
         ScrArray2D<Real> &qr_j) {
        const Real sin_jp1 = sin_face_xi(m,k,j+1);
        const Real cos_jp1 = cos_face_xi(m,k,j+1);
        const Real sin_j = sin_face_xi(m,k,j);
        const Real cos_j = cos_face_xi(m,k,j);
        par_for_inner(member, il, iu, [&](const int i) {
          ql_jp1(IVZ,i) += cos_jp1 * ql_jp1(IVY,i);
          ql_jp1(IVY,i) *= sin_jp1;
          qr_j(IVZ,i)   += cos_j * qr_j(IVY,i);
          qr_j(IVY,i)   *= sin_j;
        });
      return;
    };

    // --- x3 sweep: ETA face, normal component is IVZ, tangential is IVY
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquianglePrimFaceX3(TeamMember_t const &member, const int m, const int k,
         const int j, const int il, const int iu, ScrArray2D<Real> &ql_kp1,
         ScrArray2D<Real> &qr_k) {
        const Real sin_kp1 = sin_face_eta(m,k+1,j);
        const Real cos_kp1 = cos_face_eta(m,k+1,j);
        const Real sin_k = sin_face_eta(m,k,j);
        const Real cos_k = cos_face_eta(m,k,j);
        par_for_inner(member, il, iu, [&](const int i) {
          ql_kp1(IVY,i) += cos_kp1 * ql_kp1(IVZ,i);
          ql_kp1(IVZ,i) *= sin_kp1;
          qr_k(IVY,i)   += cos_k * qr_k(IVZ,i);
          qr_k(IVZ,i)   *= sin_k;
        });
      return;
    };

    // --- MAGNETIC FIELD into the same orthonormal frames ------------------------------
    //
    // The face-centred field is stored as the FLUX DENSITY through its own face, i.e.
    // b.x2f = B.nhat_xi and b.x3f = B.nhat_eta, which is what makes  sum(area*b)  the
    // physical flux and lets the curvilinear CT update in mhd_ct.cpp use area/dxedge
    // unchanged. bcc is the average of two opposing faces, so it carries the same
    // convention: bcc = (B.rhat, B.nhat_xi, B.nhat_eta) = (B^r, s*B^xi, s*B^eta).
    //
    // That differs from w0, which holds the CONTRAVARIANT velocity, by exactly one
    // factor of s -- so B is put into the sweep's orthonormal frame by dividing by s and
    // then applying the same rotation the velocity gets. In every sweep the NORMAL
    // component comes from b.x*f and is already correct, and the frame's two orthonormal
    // face-parallel axes are the two edge directions, so only ONE of the three slots the
    // Riemann solver reads is actually wrong. Slot ordering matches wl/wr: 0,1,2 = r,
    // xi, eta.
    //
    // x1 sweep, frame {rhat, e_xi, (e_eta - c e_xi)/s}: slot 1 must be B.e_xi.
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFaceBX1(TeamMember_t const &member, const int m, const int k,
         const int j, const int il, const int iu, ScrArray2D<Real> &bl,
         ScrArray2D<Real> &br) {
        const Real sin_theta = sin_cell(m,k,j);
        const Real cos_theta = cos_cell(m,k,j);
        par_for_inner(member, il, iu, [&](const int i) {
          bl(1,i+1) = (bl(1,i+1) + cos_theta*bl(2,i+1))/sin_theta;
          br(1,i)   = (br(1,i)   + cos_theta*br(2,i))/sin_theta;
        });
      return;
    };

    // x2 sweep, frame {nhat_xi, e_eta, rhat}: slot 2 must be B.e_eta.
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFaceBX2(TeamMember_t const &member, const int m, const int k,
         const int j, const int il, const int iu, ScrArray2D<Real> &bl_jp1,
         ScrArray2D<Real> &br_j) {
        const Real sin_jp1 = sin_face_xi(m,k,j+1);
        const Real cos_jp1 = cos_face_xi(m,k,j+1);
        const Real sin_j = sin_face_xi(m,k,j);
        const Real cos_j = cos_face_xi(m,k,j);
        par_for_inner(member, il, iu, [&](const int i) {
          bl_jp1(2,i) = (cos_jp1*bl_jp1(1,i) + bl_jp1(2,i))/sin_jp1;
          br_j(2,i)   = (cos_j*br_j(1,i)     + br_j(2,i))/sin_j;
        });
      return;
    };

    // x3 sweep, frame {nhat_eta, rhat, e_xi}: slot 1 must be B.e_xi.
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFaceBX3(TeamMember_t const &member, const int m, const int k,
         const int j, const int il, const int iu, ScrArray2D<Real> &bl_kp1,
         ScrArray2D<Real> &br_k) {
        const Real sin_kp1 = sin_face_eta(m,k+1,j);
        const Real cos_kp1 = cos_face_eta(m,k+1,j);
        const Real sin_k = sin_face_eta(m,k,j);
        const Real cos_k = cos_face_eta(m,k,j);
        par_for_inner(member, il, iu, [&](const int i) {
          bl_kp1(1,i) = (bl_kp1(1,i) + cos_kp1*bl_kp1(2,i))/sin_kp1;
          br_k(1,i)   = (br_k(1,i)   + cos_k*br_k(2,i))/sin_k;
        });
      return;
    };

    // --- EMF returned by the x1 sweep -------------------------------------------------
    // CT integrates E around a face, so e2/e3 must be E projected on the EDGE directions
    // e_xi and e_eta, whose lengths dxedge.x2e/x3e it multiplies them by. Of the three
    // sweeps only x1 returns something else: its orthonormal frame is
    // {rhat, e_xi, (e_eta - c e_xi)/s}, whose third axis is NOT the eta edge, so
    //   E.e_eta = c*(E.e_xi) + s*(E.f2)  =  c*e21 + s*e31.
    // The x2 and x3 frames are {nhat, e_eta, rhat} and {nhat, rhat, e_xi}: both of their
    // face-parallel axes ARE edge directions, so those EMFs need no rotation.
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleEmfX1(TeamMember_t const &member, const int m, const int k,
         const int j, const int il, const int iu, DvceArray4D<Real> e31,
         DvceArray4D<Real> e21) {
        const Real sin_theta = sin_cell(m,k,j);
        const Real cos_theta = cos_cell(m,k,j);
        par_for_inner(member, il, iu, [&](const int i) {
          e31(m,k,j,i) = cos_theta*e21(m,k,j,i) + sin_theta*e31(m,k,j,i);
        });
      return;
    };

    // --- x1 flux: radial. IM2 comes back unchanged, which is correct: the orthonormal
    // frame used above is {e_xi, n-hat} with n-hat perpendicular to e_xi, so the first
    // component already IS the covariant m_2.
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFluxX1(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
                                 DvceArray5D<Real> flx) {
        const Real sin_theta = sin_cell(m,k,j);
        const Real cos_theta = cos_cell(m,k,j);
        par_for_inner(member, il, iu, [&](const int i) {
          Real fb = flx(m,IM3,k,j,i)/sin_theta;
          Real fa = flx(m,IM2,k,j,i) - fb*cos_theta;
          flx(m,IM2,k,j,i) = fa + fb*cos_theta;
          flx(m,IM3,k,j,i) = fb + fa*cos_theta;
        });
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFluxX2(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
                                 DvceArray5D<Real> flx) {
        const Real sin_theta = sin_face_xi(m,k,j);
        const Real cos_theta = cos_face_xi(m,k,j);
        const Real T22 = 1.0/sin_theta;
        const Real T32 = -cos_theta/sin_theta;
        par_for_inner(member, il, iu, [&](const int i) {
          Real fb = flx(m,IM3,k,j,i) + T32 * flx(m,IM2,k,j,i);
          Real fa = T22 * flx(m,IM2,k,j,i);
          flx(m,IM2,k,j,i) = fa + fb*cos_theta;
          flx(m,IM3,k,j,i) = fb + fa*cos_theta;
        });
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFluxX3(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
                                 DvceArray5D<Real> flx) {
        const Real sin_theta = sin_face_eta(m,k,j);
        const Real cos_theta = cos_face_eta(m,k,j);
        const Real T23 = -cos_theta/sin_theta;
        const Real T33 = 1.0/sin_theta;
        par_for_inner(member, il, iu, [&](const int i) {
          Real fa = flx(m,IM2,k,j,i) + T23 * flx(m,IM3,k,j,i);
          Real fb = T33 * flx(m,IM3,k,j,i);
          flx(m,IM2,k,j,i) = fa + fb*cos_theta;
          flx(m,IM3,k,j,i) = fb + fa*cos_theta;
        });
      return;
    };

  // functions
  void CoordSrcTerms(const DvceArray5D<Real> &w0, const EOS_Data &eos, const Real dt,
                     DvceArray5D<Real> &u0);
  void CoordSrcTerms(const DvceArray5D<Real> &w0, const DvceArray5D<Real> &bcc,
                     const EOS_Data &eos, const Real dt, DvceArray5D<Real> &u0);
  void SetExcisionMasks(DvceArray4D<bool> &floor, DvceArray4D<bool> &flux);

  void UpdateExcisionMasks();

 private:
  MeshBlockPack* pmy_pack;
};

#endif // COORDINATES_COORDINATES_HPP_
