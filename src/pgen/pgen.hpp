#ifndef PGEN_PGEN_HPP_
#define PGEN_PGEN_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pgen.hpp
//  \brief definitions for ProblemGenerator class

#include <functional>
#include <memory>
#include <vector>

#include "geodesic-grid/spherical_grid.hpp"
#include "parameter_input.hpp"

// --- PGEN STATE IN THE RESTART FILE.  A problem generator that carries RELAXED internal
// state -- red_giant's shell-mean MLT profile fmlt1d, which follows the shell flux
// deficit over problem/mlt_relax_time -- has to store it, or every restart re-seeds that
// state from scratch and hits the star with an impulse (measured: the radial kinetic
// energy jumped 2.5-4.6x in the one history interval containing the restart,
// bench/wt_he4/tests_r6/README.md section 3).  The block is OPTIONAL and MARKED: a pgen
// that leaves `pgen_rst_write_func` null writes nothing at all, so every restart file
// this code has ever written is unchanged and still readable.  The reader peeks eight
// bytes; when they are not the marker they are the variable data size that follows in
// the old layout, and nothing has been lost.
using PgenRestartStateFnPtr = std::vector<char> (*)();
// eight bytes that cannot be a plausible IOWrapperSizeT (the data size that occupies this
// position in a file without the block is O(1e6), i.e. five leading zero bytes)
constexpr char kPgenRstMagic[8] = {'P', 'G', 'E', 'N', 'S', 'T', '0', '1'};

using ProblemFinalizeFnPtr = void (*)(ParameterInput *pin, Mesh *pm);
using UserBoundaryFnPtr = void (*)(Mesh* pm);
using UserSrctermFnPtr = void (*)(Mesh* pm, const Real bdt);
using UserRefinementFnPtr = void (*)(MeshBlockPack* pmbp);
using UserHistoryFnPtr = void (*)(HistoryData *pdata, Mesh *pm);
// called ONCE PER CYCLE, after Mesh::NewTimeStep in Driver::Execute, for problem
// generators that want a per-cycle diagnostic. Not part of any task list: it sees the
// state at the end of the cycle, with the new dt already computed.
using UserCycleFnPtr = void (*)(Mesh *pm);

struct HotJupiterParam {
  // Initialised, because these are only filled when <problem>/hot_jupiter is true and the
  // Ohmic resistivity reads Rgas unconditionally on its ideal-gas branch. Left
  // uninitialised, a non-hot-Jupiter problem divided by whatever was on the stack.
  Real Teq = 0.0, omega = 0.0, grav = 0.0, ap = 0.0, Rgas = 0.0, met = 0.0, bbot = 0.0;
  // false = constant gravity g = grav (the historical thin-shell behaviour). true = a
  // point mass with the same surface gravity at ap, g(r) = grav*(ap/r)^2. See
  // GravAccAt/GravPotAt in deep_hot_jupiter_rt.cpp.
  bool grav_point_mass = false;
  // false = no stellar tide (the historical behaviour). true = add the host star's tidal
  // acceleration Omega^2 (2x, -y, -z) in the corotating frame. See TideAccR/TideAccT/
  // TideAccP in deep_hot_jupiter_rt.cpp.
  bool stellar_tide = false;
  // true = the centrifugal potential -Omega^2 R^2/2 is part of the potential the initial
  // column, the ghost column and the well-balanced background are built on (a barotropic
  // equilibrium of the TOTAL potential), and its radial part leaves the explicit source.
  // The theta part stays an explicit source. See TotPotAt/ZEffFromPot.
  bool rot_potential = false;
};

//----------------------------------------------------------------------------------------
//! \class ProblemGenerator

class ProblemGenerator {
 public:
  // constructor for new problems
  ProblemGenerator(ParameterInput *pin, Mesh *pmesh);
  // constructor for restarts
  ProblemGenerator(ParameterInput *pin, Mesh *pmesh, IOWrapper resfile,
                   bool single_file_per_rank=false);
  ~ProblemGenerator() = default;
    
  bool hot_jupiter;
  HotJupiterParam hot_jupiter_param;

  // true if user BCs are specified on any face
  bool user_bcs;

  // true if user srcterms are specified
  bool user_srcs;

  // true if user history outputs are specified
  bool user_hist;

  // vector of SphericalGrid objects for analysis
  std::vector<std::unique_ptr<SphericalGrid>> spherical_grids;

  // function pointer for final work after main loop (e.g. compute errors).  Called by
  // Driver::Finalize()
  ProblemFinalizeFnPtr pgen_final_func=nullptr;
  // function pointer for user-enrolled BCs.  Called in ApplyPhysicalBCs in task list
  UserBoundaryFnPtr user_bcs_func=nullptr;
  UserSrctermFnPtr user_srcs_func=nullptr;
  // an OPERATOR-SPLIT source applied Strang-wise around the whole time integrator --
  // half the cycle dt in "before_timeintegrator" and half in "after_timeintegrator".
  // Left null the two tasks are no-ops, so every existing run is bitwise unchanged.
  UserSrctermFnPtr user_split_func=nullptr;
  // when true the split source is applied ONCE per cycle with the FULL dt, in
  // "after_timeintegrator" only -- the "before_timeintegrator" call is skipped.
  bool user_split_once=false;
  UserRefinementFnPtr user_ref_func=nullptr;
  UserHistoryFnPtr user_hist_func=nullptr;
  // the two halves of the optional pgen state block (see kPgenRstMagic).  A pgen enrols
  // `pgen_rst_write_func` to have its bytes stored by every restart output, and consumes
  // `pgen_rststate` -- the bytes a restart read back, empty when the file carried none --
  // from its own restart path.
  PgenRestartStateFnPtr pgen_rst_write_func=nullptr;
  std::vector<char> pgen_rststate;
  // called once per cycle after Mesh::NewTimeStep in Driver::Execute
  UserCycleFnPtr user_cycle_func=nullptr;
  // predefined problem generator functions (default test suite)
  void CallProblemGenerator(ParameterInput *pin, bool is_restart);
  void Advection(ParameterInput *pin, const bool restart);
  void AlfvenWave(ParameterInput *pin, const bool restart);
  void BalsaraVortex(ParameterInput *pin, const bool restart);
  void BondiAccretion(ParameterInput *pin, const bool restart);
  void CShock(ParameterInput *pin, const bool restart);
  void Diffusion(ParameterInput *pin, const bool restart);
  void RadDiff2D(ParameterInput *pin, const bool restart);
  void LinearWave(ParameterInput *pin, const bool restart);
  void LWImplode(ParameterInput *pin, const bool restart);
  void Monopole(ParameterInput *pin, const bool restart);
  void MRI3d(ParameterInput *pin, const bool restart);
  void GreshoMHD(ParameterInput *pin, const bool restart);
  void OrszagTang(ParameterInput *pin, const bool restart);
  void ShockTube(ParameterInput *pin, const bool restart);
  void Shwave(ParameterInput *pin, const bool restart);
  void WbAtmosphere(ParameterInput *pin, const bool restart);
  void SphericalCollapse(ParameterInput *pin, const bool restart);
  void RadiationLinearWave(ParameterInput *pin, const bool restart);
  void RadiationBeam(ParameterInput *pin, const bool restart);
  void RadiationM1Tests(ParameterInput *pin, const bool restart);
  void RadiationM1Tests2(ParameterInput *pin, const bool restart);
  void Z4cBoostedPuncture(ParameterInput *pin, const bool restart);
  void Z4cLinearWave(ParameterInput *pin, const bool restart);
  void SelfGravity(ParameterInput *pin, const bool restart);
  void BinaryGravity(ParameterInput *pin, const bool restart);
  void BECollapse(ParameterInput *pin, const bool restart);

  // Generic error output function (using difference u0-u1)
  void OutputErrors(ParameterInput *pin, Mesh *pm);

  // template for user-specified problem generator
  void UserProblem(ParameterInput *pin, const bool restart);

 private:
  bool single_file_per_rank; // for restart file naming
  Mesh* pmy_mesh_;
};

#endif // PGEN_PGEN_HPP_
