#ifndef BVALS_BVALS_HPP_
#define BVALS_BVALS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file bvals.hpp
//! \brief defines classes for handling boundary values for both particles as well as all
//! types of Mesh variables. For Mesh variables, methods for cell-centered and
//! face-centered fields are currently implemented, based on derived classes from the
//! generic MeshBoundaryValue class.  A separate ParticlesBoundaryValues class is
//! implemented for partciles.

// identifiers for all 6 faces of a MeshBlock
enum BoundaryFace {undef=-1, inner_x1, outer_x1, inner_x2, outer_x2, inner_x3, outer_x3};

// identifiers for boundary conditions
enum class BoundaryFlag {undef=-1,block, panel, polar, reflect, inflow, outflow, diode, user, periodic,
                         shear_periodic, vacuum};

#include <algorithm>
#include <vector>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "tasklist/task_list.hpp"
//#include "particles/particles.hpp"

// Forward declarations
class MeshBlockPack;
namespace particles {
class Particles;
}

//----------------------------------------------------------------------------------------
//! \fn bool IsCubeVertexCorner(nghbr, mbpanel, m, n)
//! \brief Is x2x3-edge buffer index n (40..47) of MeshBlock m a cubed-sphere CUBE VERTEX?
//!
//! True when BOTH flanking FACE neighbours are on a different panel, which happens
//! only at a corner of the cube. Only THREE panels meet there, so no fourth block sits
//! diagonally across it and the generic pairing is meaningless: it points at the WRONG
//! CORNER of the right panel. A reciprocity audit of the whole neighbour table found
//! these to be the ONLY non-reciprocal slots -- every face and non-vertex diagonal pairs
//! correctly -- and non-reciprocal is fatal under MPI, because a send is tagged with the
//! RECEIVER's (lid, dest): the vertex slot is a posted receive that nobody satisfies, and
//! its own send collides with the legitimate sender for the slot it targets. That is a
//! guaranteed hang in ClearRecv, and it is why a cubed-sphere run on more than one rank
//! never got past cycle 0.
//!
//! So the exchange is skipped there, on the send and receive sides together. Nothing is
//! lost: FillPanelCornersCC / FillPanelCornersFC overwrite exactly that ng x ng corner
//! block afterwards by extrapolating from the two flanking face halos, and
//! flux_correct_fc's AveragePanelCornerEMF does the same for the corner EMF. The
//! predicate is local and symmetric, so every surviving slot has exactly one sender.
//!
//! Written as a template so it can be called with either the device or the host mirror.

template <class NghbrView, class PanelView>
KOKKOS_INLINE_FUNCTION
bool IsCubeVertexCorner(const NghbrView &nghbr, const PanelView &mbpanel,
                        const int m, const int n) {
  if (n < 40 || n >= 48) return false;
  const int c = (n - 40)/2;                        // 0..3 over (sign x2, sign x3)
  const int nj_id = (c & 1) ? 12 : 8;              // flanking x2 face
  const int nk_id = (c & 2) ? 28 : 24;             // flanking x3 face
  const int mp = mbpanel(m);
  return (nghbr(m,nj_id).gid >= 0 && nghbr(m,nj_id).panel != mp &&
          nghbr(m,nk_id).gid >= 0 && nghbr(m,nk_id).panel != mp);
}

//----------------------------------------------------------------------------------------
//! \fn bool IsSkippedPanelDiagonal(nghbr, mbpanel, mblev, m, n)
//! \brief Is buffer n of MeshBlock m a cubed-sphere EDGE or CORNER slot whose neighbour
//! is on another panel AND at another refinement level?
//!
//! Such a slot cannot be paired reliably. The cross-panel corner and edge pairing was
//! only ever made reciprocal for the CUBE VERTEX (IsCubeVertexCorner); the generic
//! diagonal was left alone because with one MeshBlock per panel almost every cross-panel
//! diagonal IS a vertex. Refine a panel and that stops being true: many ordinary
//! cross-panel diagonals appear, at the same level and at mixed levels, and a
//! reciprocity audit shows they are NOT reciprocal. Since a send is tagged with the
//! RECEIVER's (lid, dest), a non-reciprocal slot is a posted receive nobody satisfies --
//! a guaranteed hang in ClearRecv, which is exactly what a refined cubed sphere did on
//! more than one rank.
//!
//! Gated on `multilevel` so that the UNREFINED cubed sphere, which is validated and
//! whose non-vertex diagonals do pair, keeps its existing behaviour untouched.
//!
//! FACES are excluded here on purpose -- a face at mixed level across a seam is a level
//! boundary lying on a seam, which prolongation and restriction cannot do at all, and
//! CheckCubedSphereRefinement refuses it at startup instead of silently skipping it.
//!
//! Skipping is safe for the same reason it is at a cube vertex: a DIMENSIONALLY SPLIT
//! sweep reconstructs along one axis at a time, reading a straight line of cells, so it
//! touches face ghosts only and never these diagonal blocks. FillPanelCornersCC also
//! rewrites the tangential panel-corner block afterwards from the flanking face halos.
//!
//! Note it is the SPLITTING that makes this safe, not the order of reconstruction: PPM
//! or WENO still read along one axis, just a longer line. What DOES read diagonally is
//! the corner EMF of the CT update (see the note on FillPanelCornersCC, where O(1)
//! garbage there cost 1.3% of the magnetic energy on cycle 1) and the cross-derivatives
//! of viscosity and conduction. Both are moot today -- MHD with refinement, and
//! viscosity and conduction at all, are startup FATALs on the cubed sphere -- so nothing
//! that survives those guards looks at these cells. **Whoever lifts one of those guards
//! must deal with this**: the ghosts would be stale, and the run would complete looking
//! plausible rather than fail.
//!
//! The predicate is local and symmetric -- both sides see the same panels and levels --
//! so every surviving slot still has exactly one sender.

template <class NghbrView, class PanelView>
KOKKOS_INLINE_FUNCTION
bool IsSkippedPanelDiagonal(const NghbrView &nghbr, const PanelView &mbpanel,
                            const bool multilevel, const int m, const int n) {
  if (!multilevel) return false;
  const bool is_face = (n < 16) || (n >= 24 && n < 32);
  if (is_face) return false;
  return (nghbr(m,n).gid >= 0 && nghbr(m,n).panel != mbpanel(m));
}

//----------------------------------------------------------------------------------------
//! \fn int CreateBvals_MPI_Tag(int lid, int bufid)
//! \brief calculate an MPI tag for boundary buffer communications.  Note maximum size of
//! lid that can be encoded is set by (NUM_BITS_LID) macro defined in athena.hpp.
//! The convention in AthenaK is lid and bufid are both for the *receiving* process.
static int CreateBvals_MPI_Tag(int lid, int bufid) {
  return (bufid << (NUM_BITS_LID)) | lid;
}

//----------------------------------------------------------------------------------------
//! \struct BufferIndcs
//! \brief indices for range of cells packed/unpacked into boundary buffers

struct MeshBufferIndcs {
  int bis,bie,bjs,bje,bks,bke;  // start/end buffer ("b") indices in each dir
  MeshBufferIndcs() :
    bis(0), bie(0), bjs(0), bje(0), bks(0), bke(0) {}
};

//----------------------------------------------------------------------------------------
//! \struct MeshBoundaryBuffer
//! \brief container for index ranges, storage, and flags for boundary buffers

struct MeshBoundaryBuffer {
  // fixed-length-3 arrays used to store indices of each buffer for cell-centered vars, or
  // each component of a face-centered vector field ([0,1,2] --> [x1f, x2f, x3f]). For
  // cell-centered variables only first [0] component of index arrays are needed.
  MeshBufferIndcs isame[3];  // indices for pack/unpack when dest/src at same level
  MeshBufferIndcs icoar[3];  // indices for pack/unpack when dest/src at coarser level
  MeshBufferIndcs ifine[3];  // indices for pack/unpack when dest/src at finer level
  MeshBufferIndcs iprol[3];  // indices for prolongation (only used for receives)
  MeshBufferIndcs iflux_same[3];  // indices for pack/unpack for flux correction
  MeshBufferIndcs iflux_coar[3];  // indices for pack/unpack for flux correction
  // With Z4c higher-order prolongation/rstriction, must also send coarse data between
  // MeshBlocks at the same level, which requires an additional indices array
  MeshBufferIndcs isame_z4c;  // indices for pack/unpack with z4c when dst/src at same lvl

  // Maximum number of data elements (bie-bis+1) across 3 components of above
  int isame_ndat, isame_z4c_ndat, icoar_ndat, ifine_ndat, iflxs_ndat, iflxc_ndat;

  // 2D Views that store buffer data on device, dimensioned (nmb, ndata)
  DvceArray2D<Real> vars, flux;

#if MPI_PARALLEL_ENABLED
  // vectors of length (number of MBs) to hold MPI requests
  // Using STL vector causes problems with some GPU compilers, so just use plain C array
  MPI_Request *vars_req, *flux_req;
#endif

  // function to allocate memory for buffers for variables and their fluxes
  // Must only be called after BufferIndcs above are initialized
  void AllocateBuffers(int nmb, int nvars, bool is_z4c) {
    // With Z4c, buffers may contain BOTH same and coarse data
    if (is_z4c) {
      int nmax = std::max(isame_z4c_ndat, std::max(icoar_ndat, ifine_ndat) );
      Kokkos::realloc(vars, nmb, (nvars*nmax));
    } else {
      int nmax = std::max(isame_ndat, std::max(icoar_ndat, ifine_ndat) );
      Kokkos::realloc(vars, nmb, (nvars*nmax));
    }
    int nmax = std::max(iflxs_ndat, iflxc_ndat);
    Kokkos::realloc(flux, nmb, (nvars*nmax));
  }
};

// Forward declarations
class MeshBlockPack;

//----------------------------------------------------------------------------------------
//! \class MeshBoundaryValues
//  \brief Abstract base class for boundary values for different kinds of Mesh variables

class MeshBoundaryValues {
 public:
  MeshBoundaryValues(MeshBlockPack *ppack, ParameterInput *pin, bool z4c);
  ~MeshBoundaryValues();

  // data for all 56 buffers in most general 3D case. Not all elements used in most cases.
  // However each MeshBoundaryBuffer is lightweight, so the convenience of fixed array
  // sizes and index values for array elements outweighs cost of extra memory.
  MeshBoundaryBuffer sendbuf[56], recvbuf[56];

  // constant inflow states at each face, initialized in problem generator
  DualArray2D<Real> u_in, b_in, i_in;

#if MPI_PARALLEL_ENABLED
  // unique MPI communicators for each case (variables/fluxes)
  MPI_Comm comm_vars, comm_flux;
#endif

  //functions
  virtual void InitSendIndices(MeshBoundaryBuffer &buf,int x,int y,int z,int a,int b)=0;
  virtual void InitRecvIndices(MeshBoundaryBuffer &buf,int x,int y,int z,int a,int b)=0;
  void InitializeBuffers(const int nvar);

  TaskStatus InitRecv(const int nvar);
  virtual TaskStatus InitFluxRecv(const int nvar)=0;
  TaskStatus ClearRecv();
  TaskStatus ClearSend();
  TaskStatus ClearFluxRecv();
  TaskStatus ClearFluxSend();

  // BCs associated with various physics modules
  static void HydroBCs(MeshBlockPack *pp, DualArray2D<Real> uin, DvceArray5D<Real> u0);
  static void BFieldBCs(MeshBlockPack *pp, DualArray2D<Real> bin, DvceFaceFld4D<Real> b0);
  static void RadiationBCs(MeshBlockPack *pp,DualArray2D<Real> iin,DvceArray5D<Real> i0);
  static void Z4cBCs(MeshBlockPack *pp, DualArray2D<Real> uin, DvceArray5D<Real> u0,
                     DvceArray5D<Real> coarse_u0);
  static void PolarAzimuthalAverageBxBy(MeshBlockPack *pp, DvceFaceFld4D<Real> b0);

 protected:
  // must use pointer to MBPack and not parent physics module since parent can be one of
  // many types (Hydro, MHD, Radiation, Z4c, etc.)
  MeshBlockPack* pmy_pack;
  bool is_z4c_;   // flag to denote if this BoundaryValues is for Z4c module
};

//----------------------------------------------------------------------------------------
//! \class BoundaryValuesCC
//  \brief Derived class implementing boundary values for cell-centered variables

class MeshBoundaryValuesCC : public MeshBoundaryValues {
 public:
  MeshBoundaryValuesCC(MeshBlockPack *ppack, ParameterInput *pin, bool z4c);

  //functions
  void InitSendIndices(MeshBoundaryBuffer &b,int o1,int o2,int o3,int f1,int f2) override;
  void InitRecvIndices(MeshBoundaryBuffer &b,int o1,int o2,int o3,int f1,int f2) override;
  TaskStatus InitFluxRecv(const int nvar) override;

  // functions to communicate CC data
  TaskStatus PackAndSendCC(DvceArray5D<Real> &a, DvceArray5D<Real> &ca);
  TaskStatus RecvAndUnpackCC(DvceArray5D<Real> &a, DvceArray5D<Real> &ca);
  void FillPanelCornersCC(DvceArray5D<Real> &a);
  // functions to communicate fluxes of CC data
  TaskStatus PackAndSendFluxCC(DvceFaceFld5D<Real> &flx);
  TaskStatus RecvAndUnpackFluxCC(DvceFaceFld5D<Real> &flx);
  // functions to make the flux through a cubed-sphere PANEL SEAM single-valued
  TaskStatus InitFluxSeamRecv(const int nvar);
  TaskStatus PackAndSendFluxSeamCC(DvceFaceFld5D<Real> &flx);
  TaskStatus RecvAndUnpackFluxSeamCC(DvceFaceFld5D<Real> &flx);

  // functions to prolongate conserved and primitive CC variables
  void FillCoarseInBndryCC(DvceArray5D<Real> &a, DvceArray5D<Real> &ca,
       bool is_z4c=false);
  void ProlongateCC(DvceArray5D<Real> &a, DvceArray5D<Real> &ca, bool is_z4c=false);
  void ConsToPrimCoarseBndry(const DvceArray5D<Real> &cons, DvceArray5D<Real> &prim);
  void PrimToConsFineBndry(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons);
  void ConsToPrimCoarseBndry(const DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                             DvceArray5D<Real> &prim);
  void PrimToConsFineBndry(const DvceArray5D<Real> &prim, const DvceFaceFld4D<Real> &b,
                           DvceArray5D<Real> &cons);
};

//----------------------------------------------------------------------------------------
//! \class BoundaryValuesFC
//  \brief Derived class implementing boundary values for face-centered vector fields

class MeshBoundaryValuesFC : public MeshBoundaryValues {
 public:
  MeshBoundaryValuesFC(MeshBlockPack *ppack, ParameterInput *pin);

  //functions
  void InitSendIndices(MeshBoundaryBuffer &b,int o1,int o2,int o3,int f1,int f2) override;
  void InitRecvIndices(MeshBoundaryBuffer &b,int o1,int o2,int o3,int f1,int f2) override;
  TaskStatus InitFluxRecv(const int nvar) override;

  TaskStatus PackAndSendFC(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);
  TaskStatus RecvAndUnpackFC(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);
  void FillCoarseInBndryFC(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);
  void FillPanelCornersFC(DvceFaceFld4D<Real> &b);
  void ProlongateFC(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);

  TaskStatus PackAndSendFluxFC(DvceEdgeFld4D<Real> &flx);
  TaskStatus RecvAndUnpackFluxFC(DvceEdgeFld4D<Real> &flx);
  void SumBoundaryFluxes(DvceEdgeFld4D<Real> &flx, const bool same_level,
                         DvceArray2D<int> &nflx);
  void ZeroFluxesAtBoundaryWithFiner(DvceEdgeFld4D<Real> &flx, DvceArray2D<int> &nflx);
  void AverageBoundaryFluxes(DvceEdgeFld4D<Real> &flx, DvceArray2D<int> &nflx);
  void SavePanelCornerEMF(DvceEdgeFld4D<Real> &flx, DvceArray3D<Real> &save);
  void AveragePanelCornerEMF(DvceEdgeFld4D<Real> &flx, DvceArray3D<Real> &save);
};

//----------------------------------------------------------------------------------------
//! \struct ParticleLocationData
//! \brief data describing location of data for particles communicated with MPI

struct ParticleLocationData {
  int prtcl_indx;   // index in particle array
  int dest_gid;     // GID of target MeshBlock
  int dest_rank;    // rank of target MeshBlock
};

// Custom operators to sort ParticleLocationData array by dest_rank or prtcl_indx
struct {
  bool operator()(ParticleLocationData a, ParticleLocationData b)
    const { return a.dest_rank < b.dest_rank; }
} SortByRank;
struct {
  bool operator()(ParticleLocationData a, ParticleLocationData b)
    const { return a.prtcl_indx < b.prtcl_indx; }
} SortByIndex;

//----------------------------------------------------------------------------------------
//! \struct ParticleMessageData
//! \brief Data describing MPI messages containing particles

struct ParticleMessageData {
  int sendrank;  // rank of sender
  int recvrank;  // rank of receiver
  int nprtcls;   // number of particles in message
  ParticleMessageData(int a, int b, int c) :
    sendrank(a), recvrank(b), nprtcls(c) {}
};

//----------------------------------------------------------------------------------------
//! \class ParticlesBoundaryValues
//  \brief Defines boundary values class for particles

namespace particles {
class ParticlesBoundaryValues {
 public:
  ParticlesBoundaryValues(particles::Particles *ppart, ParameterInput *pin);
  ~ParticlesBoundaryValues();

  int nprtcl_send, nprtcl_recv;
  DualArray1D<ParticleLocationData> sendlist;

  // Data needed to count number of messages and particles to send between ranks
  int nsends; // number of MPI sends to neighboring ranks on this rank
  int nrecvs; // number of MPI recvs from neighboring ranks on this rank
  std::vector<int> nsends_eachrank;                // length nranks
  std::vector<ParticleMessageData> sends_thisrank; // length nsends
  std::vector<ParticleMessageData> recvs_thisrank; // length nrecvs
  std::vector<ParticleMessageData> sends_allranks; // length ncounts summed over ranks

#if MPI_PARALLEL_ENABLED
  DvceArray1D<Real> prtcl_rsendbuf, prtcl_rrecvbuf;
  DvceArray1D<int>  prtcl_isendbuf, prtcl_irecvbuf;
  std::vector<MPI_Request> rrecv_req, rsend_req;  // vectors of requests for Reals
  std::vector<MPI_Request> irecv_req, isend_req;  // vectors of requests for ints
  MPI_Comm mpi_comm_part;                       // unique MPI communicators for particles
#endif

  //functions
  TaskStatus SetNewPrtclGID();
  TaskStatus CountSendsAndRecvs();
  TaskStatus InitPrtclRecv();
  TaskStatus ClearPrtclRecv();
  TaskStatus PackAndSendPrtcls();
  TaskStatus ClearPrtclSend();
  TaskStatus RecvAndUnpackPrtcls();

 protected:
  particles::Particles* pmy_part;
};
} // namespace particles

#endif // BVALS_BVALS_HPP_
