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
  DvceArray3D<Real> sin_cell;
  DvceArray3D<Real> cos_cell;
  DvceArray3D<Real> sin_face1;
  DvceArray3D<Real> cos_face1;
  DvceArray3D<Real> sin_face2;
  DvceArray3D<Real> cos_face2;
  DvceArray4D<Real> x_ov_rD;
  DvceArray4D<Real> y_ov_rC;
  DvceArray4D<Real> z_ov_rE;
    
  void CoordSphericalPolar();
  void SrcTermsSphericalPolarHydro(const DvceArray5D<Real> &w0, const DvceArray5D<Real> &w0wb, const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
  void SrcTermsSphericalPolarMHD(const DvceArray5D<Real> &w0, const DvceArray5D<Real> &bcc0, const DvceArray5D<Real> &w0wb, const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
  KOKKOS_INLINE_FUNCTION
  void StretchR(const Real a, const Real r0, const Real r1, Real &r) {
    Real xi = (r-r0)/(r1-r0);
    Real denom = 1.0 - exp(-a);
    r = r0 + (r1 - r0)*(1.0 - exp(-a*xi))/denom;
//    r = r0*pow(r1/r0,xi);
  };
  KOKKOS_INLINE_FUNCTION
  void StretchTheta(const Real a, Real &t) {
    Real xi = t/M_PI;
    t = M_PI/2.0*(1.0+sinh(a*(2.0*xi-1.0))/sinh(a));
  };
    
  void CoordGnomonicEquiangle();
  void SrcTermsGnomonicEquiangle(const DvceArray5D<Real> &w0, const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
    
    KOKKOS_INLINE_FUNCTION
    void GnomonicEquianglePrimFaceX1(TeamMember_t const &member, const int m, const int j, const int il, const int iu,
         ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real sin_theta = sin_face1(m,j,i+1);
          Real cos_theta = cos_face1(m,j,i+1);
          Real T11 = sin_theta;
          Real T21 = cos_theta;
          ql(IVY,i+1) += T21 * ql(IVX,i+1);
          ql(IVX,i+1) *= T11;
          sin_theta = sin_face1(m,j,i);
          cos_theta = cos_face1(m,j,i);
          T11 = sin_theta;
          T21 = cos_theta;
          qr(IVY,i)   += T21 * qr(IVX,i);
          qr(IVX,i)   *= T11;
        });
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void GnomonicEquianglePrimFaceX2(TeamMember_t const &member, const int m, const int j, const int il, const int iu,
         ScrArray2D<Real> &ql_jp1, ScrArray2D<Real> &qr_j) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real sin_theta = sin_face2(m,j+1,i);
          Real cos_theta = cos_face2(m,j+1,i);
          Real T12 = cos_theta;
          Real T22 = sin_theta;
          ql_jp1(IVX,i) += T12 * ql_jp1(IVY,i);
          ql_jp1(IVY,i) *= T22;
          sin_theta = sin_face2(m,j,i);
          cos_theta = cos_face2(m,j,i);
          T12 = cos_theta;
          T22 = sin_theta;
          qr_j(IVX,i)   += T12 * qr_j(IVY,i);
          qr_j(IVY,i)   *= T22;
        });
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void GnomonicEquianglePrimFaceX3(TeamMember_t const &member, const int m, const int j, const int il, const int iu,
         ScrArray2D<Real> &ql_kp1, ScrArray2D<Real> &qr_k) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real sin_theta = sin_cell(m,j,i);
          Real cos_theta = cos_cell(m,j,i);
          ql_kp1(IVX,i) += cos_theta * ql_kp1(IVY,i);
          ql_kp1(IVY,i) *= sin_theta;
          qr_k(IVX,i)   += cos_theta * qr_k(IVY,i);
          qr_k(IVY,i)   *= sin_theta;
        });
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFluxX1(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
                                 DvceArray5D<Real> flx) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real sin_theta = sin_face1(m,j,i);
          Real cos_theta = cos_face1(m,j,i);
          Real T11 = 1.0/sin_theta;
          Real T21 = -1.0/sin_theta*cos_theta;
          Real fy = flx(m,IM2,k,j,i) + T21 * flx(m,IM1,k,j,i);
          Real fx = T11 * flx(m,IM1,k,j,i);
          flx(m,IM1,k,j,i) = fx + fy*cos_theta;
          flx(m,IM2,k,j,i) = fy + fx*cos_theta;
        });
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFluxX2(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
                                 DvceArray5D<Real> flx) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real sin_theta = sin_face2(m,j,i);
          Real cos_theta = cos_face2(m,j,i);
          Real T12 = -1.0/sin_theta*cos_theta;
          Real T22 = 1.0/sin_theta;
          Real fx = flx(m,IM1,k,j,i) + T12 * flx(m,IM2,k,j,i);
          Real fy = T22 * flx(m,IM2,k,j,i);
          flx(m,IM1,k,j,i) = fx + fy*cos_theta;
          flx(m,IM2,k,j,i) = fy + fx*cos_theta;
        });
      return;
    };

    KOKKOS_INLINE_FUNCTION
    void GnomonicEquiangleFluxX3(TeamMember_t const &member, const int m, const int k, const int j, const int il, const int iu,
                                 DvceArray5D<Real> flx) {
        par_for_inner(member, il, iu, [&](const int i) {
          Real sin_theta = sin_cell(m,j,i);
          Real cos_theta = cos_cell(m,j,i);
          Real fy = flx(m,IM2,k,j,i)/sin_theta;
          Real fx = flx(m,IM1,k,j,i) - fy*cos_theta;
          flx(m,IM1,k,j,i) = fx + fy*cos_theta;
          flx(m,IM2,k,j,i) = fy + fx*cos_theta;
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
