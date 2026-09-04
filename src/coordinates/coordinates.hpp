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
  // The face-sum (well-balanced) geometric source, shared by the cubed sphere and the
  // spherical-polar grid; the grid enters only through the geometry cache.
  void SrcTermsCurvilinearWB(const DvceArray5D<Real> &w0,
       const DvceArray5D<Real> &bcc0, const bool is_mhd,
       const DvceArray5D<Real> &wder, const EOS_Data &eos_data, const Real bdt,
       DvceArray5D<Real> &u0, const DvceArray4D<Real> &pwb, const bool use_wb_static);
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
    // The eight device helpers that used to live here are now FREE functions in
    // coordinates/gnomonic_kernels.hpp, taking the GnomonicTrig below. Reaching them
    // through `pmy_pack->pcoord->` made every flux kernel capture `this` and dereference
    // a host pointer on the device; see the note in that header.
    // One pair per sweep -- see the note in gnomonic_kernels.hpp on why a kernel must not
    // be handed all six.
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
  // CUBED-SPHERE DIAGNOSTIC: drop the MAGNETIC terms from the gnomonic geometric
  // source, leaving the hydro ones. The scheme is then inconsistent for MHD -- this is
  // a probe for whether the low-beta instability's feedback lives in that source.
  bool cs_diag_no_magsrc = false;

  // WELL-BALANCED geometric source (<mhd>/cs_wellbalanced_src on the cubed sphere,
  // <mhd>/sp_wellbalanced_src on the spherical-polar grid).  See the long note at
  // SrcTermsCurvilinearWB.
  bool cs_wellbalanced_src = false;
  bool sp_wellbalanced_src = false;
  // CARTESIAN momentum update in the two POLAR cell rows (<mhd>/sp_cart_polar_momentum):
  // the momentum balance there is formed from the face fluxes rotated to Cartesian
  // vectors with a vector-area correction, instead of local components plus the
  // cot(theta) source.  See SrcTermsSphericalPolarCartRows.
  bool sp_cart_polar_momentum = false;
  // CARTESIAN momentum update on the WHOLE cubed sphere (<mhd>/cs_cart_momentum): no
  // geometric source at all; see SrcTermsGnomonicCartMomentum.
  bool cs_cart_momentum = false;
  void SrcTermsGnomonicCartMomentum(const DvceArray5D<Real> &w0,
       const DvceArray5D<Real> &bcc0, const bool is_mhd,
       const DvceArray5D<Real> &wder, const DvceFaceFld5D<Real> uflx,
       const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
  void SrcTermsSphericalPolarCartRows(const DvceArray5D<Real> &w0,
       const DvceArray5D<Real> &bcc0, const bool is_mhd,
       const DvceArray5D<Real> &wder, const DvceFaceFld5D<Real> uflx,
       const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0);
  // Its per-(MeshBlock,k,j) geometry cache: the cell triad and the four tangential faces'
  // normals and triads depend on the angles alone, so they are built once per mesh (on
  // every call under AMR) instead of ten panel-map evaluations per cell per stage.
  static constexpr int NWBGEOM = 60;
  DvceArray4D<Real> wb_geom;   // (nmb, ncells3, ncells2, NWBGEOM)
  void BuildWBGeometry();

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
