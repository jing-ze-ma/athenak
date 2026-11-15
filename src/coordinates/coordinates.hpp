#ifndef COORDINATES_COORDINATES_HPP_
#define COORDINATES_COORDINATES_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
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
#include "coordinates/gnomonic_kernels.hpp"

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
  // following data is only used in GR calculations to compute metric (initialized,
  // since kernels copy the whole struct in every run)
  bool is_minkowski = false;                            // flag to specify Minkowski space
  Real bh_spin = 0.0;                                   // needed for GR metric
  bool bh_excise = false;                               // flag to specify excision
  Real rexcise = 0.0;                                   // excision radius (SKS)
  Real dexcise = 0.0;                                   // density inside excised region
  Real pexcise = 0.0;                                   // pressure inside excised region
  Real flux_excise_r = 0.0;                             // first-order inside this radius
  ExcisionScheme excision_scheme = ExcisionScheme::fixed;  // excision method
  Real excise_lapse = 0.0;                // if excision_scheme = lapse, excise under it
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

  // geometry of the curvilinear grids (cubed sphere, spherical polar)
  DvceArray4D<Real> volume;
  DvceFaceFld4D<Real> area;
  DvceArray4D<Real> dx1;
  DvceArray4D<Real> dx2;
  DvceArray4D<Real> dx3;
  DvceArray2D<Real> x1v;
  DvceArray2D<Real> x2v;
  // x2 positions for the angular reconstruction: x2v, or with sp_x2_periodic_image the
  // periodic-image ghost centroids at a theta-periodic x2min/x2max (BuildX2vRecon)
  DvceArray2D<Real> x2v_rec;
  bool sp_x2_periodic_image = false;
  void BuildX2vRecon();
  DvceArray2D<Real> x3v;
  DvceArray2D<Real> xx1f;
  DvceArray2D<Real> xx2f;
  DvceArray2D<Real> xx3f;
  DvceEdgeFld4D<Real> dxedge;
  DvceFaceFld4D<Real> dxface;
  DvceEdgeFld4D<Real> areaedge;
  // gnomonic-equiangle trig of the angle between the tangent basis vectors, indexed
  // (m,k,j): cell centres, xi faces (staggered in j), eta faces (staggered in k)
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
  // spherical-polar phi (<mesh>/sp_phi_global) and theta (<mesh>/layout_invariant) faces
  // and centres from the global cell index instead of the block's own edges, so the
  // coordinates are the same for every MeshBlock layout (default <mesh>/layout_invariant)
  bool sp_phi_global = false;
  bool sp_theta_global = false;
  void SrcTermsSphericalPolarHydro(const DvceArray5D<Real> &w0,
      const DvceArray4D<Real> &pwb, const DvceFaceFld5D<Real> uflx,
      const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
  void SrcTermsSphericalPolarMHD(const DvceArray5D<Real> &w0,
      const DvceArray5D<Real> &bcc0, const DvceArray4D<Real> &pwb,
      const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt,
      DvceArray5D<Real> &u0);
  void CoordGnomonicEquiangle();
  // conversions between covariant conserved momentum and contravariant primitive
  // velocity on the non-orthogonal gnomonic basis, g = [[1, cos],[cos, 1]], g_33 = 1
  void GnomonicEquiangleRaiseVel(DvceArray5D<Real> &u0, DvceArray5D<Real> &w0,
                                 const EOS_Data &eos_data,
                                 const DvceArray5D<Real> &wder,
                                 const DvceArray4D<Real> &wtemp,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku);
  void GnomonicEquiangleRaiseVelMHD(DvceArray5D<Real> &u0,
                                    const DvceFaceFld4D<Real> &b0,
                                    DvceArray5D<Real> &bcc0, DvceArray5D<Real> &w0,
                                    const EOS_Data &eos_data,
                                    const DvceArray5D<Real> &wder,
                                    const DvceArray4D<Real> &wtemp,
                                    const int il, const int iu, const int jl,
                                    const int ju, const int kl, const int ku);
  // cell-centred field part of GnomonicEquiangleRaiseVelMHD only (<mhd>/rkg_lean)
  void GnomonicCellCenteredB(const DvceFaceFld4D<Real> &b0, DvceArray5D<Real> &bcc0,
                             const int il, const int iu, const int jl, const int ju,
                             const int kl, const int ku);
  // metric part of the two routines above with the EOS cache frozen, for the first
  // conversion of a restarted run under a general EOS (gnomonic_raisevel_frozen.cpp)
  void GnomonicEquiangleRaiseVelFrozen(DvceArray5D<Real> &u0, DvceArray5D<Real> &w0,
                                       const int il, const int iu, const int jl,
                                       const int ju, const int kl, const int ku);
  void GnomonicEquiangleRaiseVelMHDFrozen(DvceArray5D<Real> &u0,
                                          const DvceFaceFld4D<Real> &b0,
                                          DvceArray5D<Real> &bcc0, DvceArray5D<Real> &w0,
                                          const int il, const int iu, const int jl,
                                          const int ju, const int kl, const int ku);
  void GnomonicEquiangleLowerMom(const DvceArray5D<Real> &w0, DvceArray5D<Real> &u0,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku);
  // geometric source terms; pwb is the static well-balanced background pressure (used
  // only with use_wellbalance_static)
  void SrcTermsGnomonicEquiangle(const DvceArray5D<Real> &w0,
      const DvceArray5D<Real> &wder, const DvceArray4D<Real> &pwb,
      const DvceFaceFld5D<Real> uflx,
      const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
  void SrcTermsGnomonicEquiangleMHD(const DvceArray5D<Real> &w0,
      const DvceArray5D<Real> &bcc0, const DvceArray5D<Real> &wder,
      const DvceArray4D<Real> &pwb,
      const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt,
      DvceArray5D<Real> &u0);
  void SrcTermsGnomonicEquiangleImpl(const DvceArray5D<Real> &w0,
      const DvceArray5D<Real> &bcc0, const bool is_mhd,
      const DvceArray5D<Real> &wder, const DvceArray4D<Real> &pwb,
      const bool use_wb_static, const DvceFaceFld5D<Real> uflx,
      const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);

    // trig pairs captured by the cubed-sphere flux kernels (gnomonic_kernels.hpp): one
    // pair per sweep (cell centres for x1, xi faces for x2, eta faces for x3)
    GnomonicTrig GnomonicTrigCell() const {
      GnomonicTrig gt; gt.sn = sin_cell; gt.cs = cos_cell; return gt;
    }
    GnomonicTrig GnomonicTrigFaceXi() const {
      GnomonicTrig gt; gt.sn = sin_face_xi; gt.cs = cos_face_xi; return gt;
    }
    GnomonicTrig GnomonicTrigFaceEta() const {
      GnomonicTrig gt; gt.sn = sin_face_eta; gt.cs = cos_face_eta; return gt;
    }

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
