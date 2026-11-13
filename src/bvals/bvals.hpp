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
enum class BoundaryFlag {undef=-1,block, panel, polar, reflect, inflow, outflow, diode,
                         user, periodic, shear_periodic, vacuum,
                         mg_zerograd, mg_zerofixed, mg_multipole};

#include <algorithm>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "tasklist/task_list.hpp"
//#include "particles/particles.hpp"

//----------------------------------------------------------------------------------------
//! \struct RankPackedVarEntry
//! \brief metadata for one (MeshBlock,neighbor) var-payload in a rank-packed message

struct RankPackedVarEntry {
  int m;
  int n;
  int lid;
  int dn;
  int data_size;
  int offset;
};

//----------------------------------------------------------------------------------------
//! \struct RankPackedVarMessage
//! \brief metadata for one aggregated MPI vars message between ranks

struct RankPackedVarMessage {
  int rank;
  int nentries;
  int entry_offset;
  int hdr_offset;
  int offset;
  int data_size;
};

// Forward declarations
class MeshBlockPack;
namespace particles {
class Particles;
}

//----------------------------------------------------------------------------------------
//! \fn bool IsCubeVertexCorner(nghbr, mbpanel, m, n)
//! \brief Is buffer index n of MeshBlock m -- an x2x3 edge (40..47) or a corner
//! (48..55) -- a cubed-sphere CUBE VERTEX?
//!
//! True when BOTH flanking FACE neighbours are on a different panel, which happens
//! only at a corner of the cube. Only THREE panels meet there, so no fourth block sits
//! diagonally across it and the generic pairing is meaningless: it points at the WRONG
//! CORNER of the right panel. A reciprocity audit of the whole neighbour table finds
//! these to be the ONLY non-reciprocal slots -- every face and every diagonal pairs
//! correctly, at one level and across a level boundary alike, once the cross-panel
//! child selection in meshblock.cpp goes through the panel map (see PanelEdgeMap) --
//! and non-reciprocal is fatal under MPI, because a send is tagged with the
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
  // x2x3 edges [40..47] AND the three-dimensional corners [48..55]. The corners were
  // once outside this predicate and never showed, because with ONE MeshBlock spanning
  // the radial direction a corner slot has no neighbour at all -- x1 ends on a physical
  // boundary. Split the radial direction, as any refined mesh does, and the cube-vertex
  // CORNER appears for the first time. It is non-reciprocal for exactly the same reason
  // as the edge and hangs ClearRecv the same way. See the note above.
  if (n < 40 || n >= 56) return false;
  int sj, sk;
  if (n < 48) {
    const int c = (n - 40)/2;                      // bit0 = sign x2, bit1 = sign x3
    sj = c & 1;  sk = (c >> 1) & 1;
  } else {
    const int c = n - 48;                          // bit0 = x1, bit1 = x2, bit2 = x3
    sj = (c >> 1) & 1;  sk = (c >> 2) & 1;
  }
  const int nj_id = sj ? 12 : 8;                   // flanking x2 face
  const int nk_id = sk ? 28 : 24;                  // flanking x3 face
  const int mp = mbpanel(m);
  return (nghbr(m,nj_id).gid >= 0 && nghbr(m,nj_id).panel != mp &&
          nghbr(m,nk_id).gid >= 0 && nghbr(m,nk_id).panel != mp);
}

//----------------------------------------------------------------------------------------
//! \fn void CsVertexFofcDrop(nghbr, mbpanel, m, k, j, js, je, ks, ke, ng, drop2, drop3)
//! \brief cornerfix-1103 <mesh>/cs_vertex_3cell: the FOFC trial state of a seam GHOST cell
//! (row k outside [ks,ke] or column j outside [js,je]) within ng of a cube-vertex corner
//! along its seam takes the along-seam flux divergence from faces whose reconstruction
//! stencils reach into the made-up corner block: drop2 / drop3 = leave the x2 / x3 terms
//! out.  Shared by Hydro::FOFC, MHD::FOFC and FofcTrialRemoveGrav.
template <class NghbrView, class PanelView>
KOKKOS_INLINE_FUNCTION
void CsVertexFofcDrop(const NghbrView &nghbr, const PanelView &mbpanel, const int m,
                      const int k, const int j, const int js, const int je, const int ks,
                      const int ke, const int ng, bool &drop2, bool &drop3) {
  drop2 = false; drop3 = false;
  const bool gk = (k < ks) || (k > ke), gj = (j < js) || (j > je);
  if (gk && !gj) {
    const int sk = (k > ke) ? 1 : 0;
    drop2 = (j < js + ng && IsCubeVertexCorner(nghbr, mbpanel, m, 40 + 2*(0 + 2*sk))) ||
            (j > je - ng && IsCubeVertexCorner(nghbr, mbpanel, m, 40 + 2*(1 + 2*sk)));
  } else if (gj && !gk) {
    const int sj = (j > je) ? 1 : 0;
    drop3 = (k < ks + ng && IsCubeVertexCorner(nghbr, mbpanel, m, 40 + 2*(sj + 0))) ||
            (k > ke - ng && IsCubeVertexCorner(nghbr, mbpanel, m, 40 + 2*(sj + 2)));
  }
}


//----------------------------------------------------------------------------------------
//! \fn bool IsX2X3DiagSlot(int n)
//! \brief Is buffer slot n an x2x3 EDGE (40..47) or a three-dimensional CORNER (48..55)?
//!
//! These are exactly the buffers that are ghost in BOTH x2 and x3, i.e. the DIAGONAL
//! transverse neighbours.  A 5-point cross stencil (k+-1 at fixed j, or j+-1 at fixed k)
//! never reads them, so an exchange dedicated to such a stencil can drop them entirely:
//! see MeshBoundaryValues::skip_x2x3_diag.  The slot
//! layout is the one nghbr_index.hpp builds (0-7 x1 faces, 8-15 x2 faces, 16-23 x1x2
//! edges, 24-31 x3 faces, 32-39 x3x1 edges, 40-47 x2x3 edges, 48-55 corners); in 1D/2D
//! nnghbr is 8/24 and no slot reaches 40, so the predicate is inert there.
KOKKOS_INLINE_FUNCTION
bool IsX2X3DiagSlot(const int n) {
  return (n >= 40 && n < 56);
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
  KOKKOS_INLINE_FUNCTION
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

//----------------------------------------------------------------------------------------
//! \struct MeshBufferData
//! \brief raw (pointer, row stride) handle on one (nmb, ndata) LayoutRight buffer View.
//! (m,i) is the element the View's own (m,i) addresses, so reads and writes through
//! it are bit-for-bit those of the View.

struct MeshBufferData {
  Real *ptr = nullptr;
  size_t s0 = 0;
  KOKKOS_INLINE_FUNCTION
  Real &operator()(const int m, const int i) const {
    return ptr[static_cast<size_t>(m)*s0 + static_cast<size_t>(i)];
  }
};

//----------------------------------------------------------------------------------------
//! \struct PackSegFC
//! \brief halofc-1102: one active (block m, buffer n, receiver component v) triple of the
//! face-field pack, with everything the old per-team setup derived from the neighbour
//! table (see MeshBoundaryValuesFC::BuildPackSegmentsFC).  Plain ints only.

struct PackSegFC {
  int m, n, v, vv;
  int il, jl, kl, ni, nj, nk, ndat;
  int dm, dn, local;
  int mode;               // 0 plain from b, 1 plain from cb, 2 seam / pole transform
  int aj, bj, ak, bk, sj, sk, signvar;
  int cs_xform, cs_coar, cs_slotidx, cs_gkind, cs_skind, cs_seam;
  int cs_psrc, cs_pdst;    // fcorder3-1103: source / destination panel
};

//----------------------------------------------------------------------------------------
//! \struct SubSegFC
//! \brief r1mspeed-1102: the sample box of one PackSegFC that a pack actually walks (the
//! whole segment, or only the source layer adjacent to the boundary for a depth-1 pack),
//! and, for a cross-panel seam segment, the box of along-seam source values (srcval) its
//! resample reads, evaluated once per pack into a scratch array instead of three times
//! per sample.  pvoff < 0: no precompute (the sample calls srcval itself).

struct SubSegFC {
  int q;                   // fc_seg index
  int il, ni, jl, nj, kl, nk;   // sample box (source-loop indices, as PackSegFC)
  int pvoff, ptlo, pnt, pnlo, pnn, pil, pni;   // srcval box: tangent x normal x i
};

//----------------------------------------------------------------------------------------
//! \struct MeshBufferDv
//! \brief DEVICE image of one MeshBoundaryBuffer: every index range and size a kernel
//! reads, and (pointer, stride) handles on its vars/flux Views.
//!
//! WHY: the 56-element host arrays sendbuf/recvbuf are ~31.5 kB EACH (19 index triples,
//! 2 Views and 2 pointers per element).  A kernel that captured them by value (every
//! pack/unpack, flux-correction and prolongation kernel) carried a 32-64 kB functor, over
//! Kokkos HIP's 32 kB constant-memory limit, so each launch went through
//! hip_parallel_launch_global_memory, whose functor staging begins with a
//! hipStreamSynchronize (kokkos/core/src/HIP/Kokkos_HIP_Instance.cpp, l. 322): one full
//! host-device sync per halo pack and per unpack.  Kernels now capture ONE View of these
//! (MeshBoundaryValues::SendBufDv/RecvBufDv) and launch through BvalsTeamFor below; the
//! arithmetic and every address are unchanged.

struct MeshBufferDv {
  MeshBufferIndcs isame[3], icoar[3], ifine[3], iprol[3], iflux_same[3], iflux_coar[3];
  MeshBufferIndcs isame_z4c;
  int isame_ndat, isame_z4c_ndat, icoar_ndat, ifine_ndat, iflxs_ndat, iflxc_ndat;
  MeshBufferData vars, flux;
};

//----------------------------------------------------------------------------------------
//! \fn void BvalsTeamFor(name, policy, functor)
//! \brief launch a boundary-value team kernel with the functor passed as a KERNEL
//! ARGUMENT (Kokkos HintLightWeight: local-memory launch for a driver < 4 kB), i.e. with
//! no host synchronisation at all.  The default mechanism for a 0.5-32 kB functor is the
//! constant-memory path, which waits on the previous constant-memory launch
//! (hip_event_synchronize on constantMemReusable) -- another hidden host stall.  The
//! static_assert keeps a functor that has grown from silently falling back to the
//! global-memory (hipStreamSynchronize) path, which HintLightWeight selects above 4 kB.

template <class F, class P = Kokkos::TeamPolicy<>>
inline void BvalsTeamFor(const char *name, const P &policy, const F &functor) {
  static_assert(sizeof(F) <= 3072, "boundary kernel functor too large for a kernel-"
                "argument launch: capture Views, not host structs (see bvals.hpp)");
#if defined(KOKKOS_ENABLE_CUDA) && (ATHENA_CUDA_BODY_GLOBAL > 0)
  // dhjcuda-1112: a large functor in device memory, its address the kernel argument
  // (athena.hpp AthBodyToDevice: no per-thread local copy of the closure; bitwise)
  if constexpr (sizeof(F) >= ATHENA_CUDA_BODY_GLOBAL) {
    const F *fp_ = AthBodyToDevice(DevExeSpace(), functor);
    Kokkos::parallel_for(name, Kokkos::Experimental::require(policy,
                         Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                         KOKKOS_LAMBDA(TeamMember_t tmember) {(*fp_)(tmember);});
  } else {
    Kokkos::parallel_for(name, Kokkos::Experimental::require(policy,
                         Kokkos::Experimental::WorkItemProperty::HintLightWeight), functor);
  }
#else
  Kokkos::parallel_for(name, Kokkos::Experimental::require(policy,
                       Kokkos::Experimental::WorkItemProperty::HintLightWeight), functor);
#endif
}

// d6fast5-1112: the team policy of the cell-centred pack / unpack (SendBuff / RecvBuff).
// On HIP a TeamPolicy<> kernel is compiled for up to 1024 threads per block, i.e. at most
// 128 VGPRs on gfx942; the cs / pole pack then spilled 816 B per lane to scratch
// (rocprofv3, HG small cs deck: SendBuff 0.86 s of 10.0 s GPU time).  LaunchBounds<256, 2>
// lifts the cap to 256 VGPRs; the copies (and so the result) are unchanged.
#if defined(KOKKOS_ENABLE_HIP)
using BvCcTeamPolicy = Kokkos::TeamPolicy<DevExeSpace, Kokkos::LaunchBounds<256, 2>>;
#else
using BvCcTeamPolicy = Kokkos::TeamPolicy<>;
#endif

// Forward declarations
class MeshBlockPack;

//----------------------------------------------------------------------------------------
//! \class MeshBoundaryValues
//  \brief Abstract base class for boundary values for different kinds of Mesh variables

// layoutbug-1103: <mesh>/cs_seam_blockend = clamp | extend, true for extend on a cubed
// sphere with more than one MeshBlock per panel side (bvals.cpp).  FRESH-start default
// extend there (prod1103, USER 11-03; recipe::NewDefaults); a restart whose input lacks
// the key keeps clamp.
bool CsSeamBlockendExtend(MeshBlockPack *ppack, ParameterInput *pin);

class MeshBoundaryValues {
 public:
  MeshBoundaryValues(MeshBlockPack *ppack, ParameterInput *pin, bool z4c);
  ~MeshBoundaryValues();

  // data for all 56 buffers in most general 3D case. Not all elements used in most cases.
  // However each MeshBoundaryBuffer is lightweight, so the convenience of fixed array
  // sizes and index values for array elements outweighs cost of extra memory.
  MeshBoundaryBuffer sendbuf[56], recvbuf[56];

  // device images of sendbuf/recvbuf for kernels (see MeshBufferDv).  Rebuilt by
  // SyncBufferDv() at the end of InitializeBuffers(), and again whenever an accessor sees
  // that a buffer View has been reallocated since.
  DvceArray1D<MeshBufferDv> SendBufDv() {
    if (BufferDvStale()) {
      SyncBufferDv();
    }
    return sbuf_dv_;
  }
  DvceArray1D<MeshBufferDv> RecvBufDv() {
    if (BufferDvStale()) {
      SyncBufferDv();
    }
    return rbuf_dv_;
  }
  void SyncBufferDv();

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

  // vmask (default 0 = every variable): post the receives for the variables whose bit
  // is set only (uniform grids: the selected variables travel as one compact block at
  // the start of each buffer, see BvVarList).  Must match the PackAndSendCC call on the
  // sending side.  Used by <mhd>/rkg_halo_ien_only.
  TaskStatus InitRecv(const int nvar, const int vmask = 0);
  virtual TaskStatus InitFluxRecv(const int nvar)=0;
  TaskStatus ClearRecv();
  TaskStatus ClearSend();
  TaskStatus ClearFluxRecv();
  TaskStatus ClearFluxSend();

  // BCs associated with various physics modules
  // `phi`: pass the module's phicc0 when <hydro|mhd>/etotgrav is on (u0's energy then
  // carries rho*Phi); left empty, the BCs are the plain copies they always were
  static void HydroBCs(MeshBlockPack *pp, DualArray2D<Real> uin, DvceArray5D<Real> u0,
                       DvceArray4D<Real> phi = DvceArray4D<Real>());
  static void BFieldBCs(MeshBlockPack *pp, DualArray2D<Real> bin, DvceFaceFld4D<Real> b0);
  static void RadiationBCs(MeshBlockPack *pp,DualArray2D<Real> iin,DvceArray5D<Real> i0);
  static void RadM1BCs(MeshBlockPack *pp, DualArray2D<Real> uin, DvceArray5D<Real> u0,
                       Real cl, Real e_floor);
  static void Z4cBCs(MeshBlockPack *pp, DualArray2D<Real> uin, DvceArray5D<Real> u0,
                     DvceArray5D<Real> coarse_u0);

 protected:
  // must use pointer to MBPack and not parent physics module since parent can be one of
  // many types (Hydro, MHD, Radiation, Z4c, etc.)
  MeshBlockPack* pmy_pack;
  bool is_z4c_;   // flag to denote if this BoundaryValues is for Z4c module

 private:
  DvceArray1D<MeshBufferDv> sbuf_dv_, rbuf_dv_;
  const Real *dv_ptr_[4][56] = {};   // vars/flux data pointers the images were built from
  bool BufferDvStale() const;

 public:
  // skip_x2x3_diag (default false): drop every DIAGONAL transverse buffer -- the x2x3
  // edges (40..47) and the corners (48..55) -- from this object's exchanges, on the
  // send, the receive, the MPI_Irecv/Isend, the completion test and the waits together,
  // so the two sides skip the identical slots and nothing is ever left unmatched.  Only
  // legal for a caller whose stencil is a 5-POINT CROSS and therefore never reads a
  // cell that is ghost in x2 and x3 at once.  Set on the RKL1 transverse-conduction
  // object alone (Conduction::pbval_tr); every other MeshBoundaryValues leaves it false
  // and is bitwise untouched.
  bool skip_x2x3_diag = false;
  // layoutbug-1103: <mesh>/cs_seam_blockend = extend (fresh-start default with more than
  // one MeshBlock per panel side, prod1103; clamp = the old halo, restarts w/o key) on
  // the fluid's conserved-state / face-field exchange (set by Hydro / MHD).  After the
  // exchange completes it is repeated once with the along-seam resample stencil allowed
  // into the sender's own along-seam ghosts at an INTERIOR block end, so the seam halo
  // does not depend on the MeshBlock layout (see RecvAndUnpackCC / RecvAndUnpackFC).
  bool cs_blockend_ext_ = false;
  bool cs_ext_pass_ = false;     // true while the second (extended) pass runs
};

//----------------------------------------------------------------------------------------
//! \struct BvVarList
//! \brief The variables a cell-centred exchange moves (rkghalo-1101): v[0..n-1] in
//! increasing order.  vmask = 0 means all nvar variables (v[vl] = vl).  At most 32.

struct BvVarList {
  int n;
  int v[32];
};
inline BvVarList MakeBvVarList(const int nvar, const int vmask) {
  BvVarList l;
  l.n = 0;
  for (int v=0; v<nvar && v<32; ++v) {
    if (vmask == 0 || ((vmask >> v) & 1)) {l.v[l.n++] = v;}
  }
  for (int q=l.n; q<32; ++q) {l.v[q] = 0;}
  return l;
}

//----------------------------------------------------------------------------------------
//! \class BoundaryValuesCC
//  \brief Derived class implementing boundary values for cell-centered variables

class MeshBoundaryValuesCC : public MeshBoundaryValues {
 public:
  MeshBoundaryValuesCC(MeshBlockPack *ppack, ParameterInput *pin, bool z4c);
  // <mesh>/cs_seam_rho_guard (set by Hydro/MHD on their conserved-state exchange only,
  // where variable 0 is the density; 0 = off): see PackAndSendCC.
  Real cs_rho_guard = 0.0;
  bool cs_lin_resample = false;   // <mesh>/cs_seam_resample = linear, same scope
  // wball-1104 <hydro|mhd>/wb_seam (opt-in, with wb_dirs = all on the cubed sphere; set by
  // Hydro/MHD on their conserved-state exchange only): the along-seam resample of the
  // DENSITY and ENERGY takes each stencil cell's state WALKED hydrostatically (isothermal,
  // or the local polytrope from the cell's inward seam-normal neighbour) to the potential
  // of the ghost point (a 5-point Lagrange interpolation of phicc0 along the source row),
  // then interpolates those: an equilibrium ghost is the equilibrium state at its own
  // potential to O(h^5) instead of the O(h^2) of resampling the stratified state itself.
  bool wb_seam = false;
  DvceArray4D<Real> wbs_phi;   // the fluid's phicc0
  EOS_Data wbs_eos;
  int wbs_opt = 1;             // WBAdvance branch: 1 isothermal, 3 local polytrope
  bool wbs_etg = false;        // E carries rho Phi (<hydro|mhd>/etotgrav)
  int wbs_ien = 4;             // the energy index of the exchanged state
  DvceArray5D<Real> wbs_bcc;   // MHD: the cell field (E carries |B|^2/2); empty for hydro
  // cs_noresample (default false): this object's cross-panel ghosts are the plain signed
  // index copy of the neighbour's own cells (no along-seam resample), i.e. the mirror
  // images of this block's cells across the seam.  Set by the implicit M1 on its scratch
  // exchanges (stage CS1, rad_m1_sph.cpp), whose seam faces pair cell to cell; the
  // global diagnostic env CS_NORESAMP_CC does the same for every object.
  bool cs_noresample = false;
  // r1mspeed-1102: >= 0 = the *FluxSeamCC exchange (averaging modes only) packs and
  // reconciles only this variable.  Set by an RKG super-stage around its flux exchange
  // (the super-stage update reads only the energy flux); -1 everywhere else.
  int seam_flux_var = -1;
  // csedge-1019: the *FluxSeamCC exchange moves RAW face values (no outward sign, the
  // receiver's array is overwritten by the neighbour's) -- the hydro's shared seam face
  // states (<mesh>/cs_seam_flux = shared, Hydro::CsSeamFaceStates)
  bool seam_raw = false;
  // layoutaudit-1103 (seam_raw only): the *FluxSeamCC exchange runs over the SAME-PANEL,
  // same-level, non-polar x2/x3 FACE neighbours (internal MeshBlock faces, incl. periodic
  // ones) instead of the cross-panel seams, with the identity index map and no basis
  // transform (one chart).  Hydro <hydro>/merge_lat_blockend = exchange.
  bool seam_same_panel = false;
  // ucthshare-1102 (seam_raw only): which variables of the exchanged face record are
  // tangential COVARIANT vector pairs (re-expressed on the receiver's basis by
  // cubed_sphere::ApplyMomentumXform).  Default: IM2/IM3 only (the hydro record).
  // seam_xform_b >= 0 adds the pair (seam_xform_b, seam_xform_b+1) (the MHD field);
  // seam_xform_stride > 0 repeats the pattern every stride variables (a record of two
  // sets: the active-radial-row states and the radial-ghost-row states); seam_xform_off
  // moves every variable verbatim (the MHD UCT seam packets, Cartesian vectors).
  int seam_xform_b = -1;
  int seam_xform_stride = 0;
  bool seam_xform_off = false;
  // cs_perm_pairs (default false): the tangential vector pairs of this object are
  // face-normal components; a seam permutes them with the axis signs instead of the
  // covariant metric transform (bvals_cc.cpp).  Set by the implicit M1 on pbval_th.
  bool cs_perm_pairs = false;
  // cs-s4-1017 cs_hrq0 >= 0 (default -1): slots cs_hrq0..+3 hold the HALF-RANGE ratios
  // (h+_2, h-_2, h+_3, h-_3) about this block's x2 / x3 axes (rad_m1 vgd_hr on the cubed
  // sphere).  Across a seam they take the axis permutation of the index map and,
  // where the destination axis runs reversed, h+ <-> -h- (no resample: the mirror
  // cell).  Exact for the seam-normal axis (the plane family is common to both charts).
  int cs_hrq0 = -1;
  // implhalo-1031 exch_depth > 0 (default 0 = ng): this object's SAME-level buffers carry
  // only the exch_depth ghost layers next to the active range (send and receive strips,
  // the cube-vertex widening); set before InitializeBuffers.  For exchanges whose readers
  // never touch the deeper ghost layers (the implicit M1 scratch arrays: poisoned-ghost
  // gate), the same values in the layers that are read with fewer bytes on the wire.
  int exch_depth = 0;

  //functions
  void InitSendIndices(MeshBoundaryBuffer &b,int o1,int o2,int o3,int f1,int f2) override;
  void InitRecvIndices(MeshBoundaryBuffer &b,int o1,int o2,int o3,int f1,int f2) override;
  TaskStatus InitFluxRecv(const int nvar) override;

  // functions to communicate CC data
  // iwl/iwu: an optional x1 INDEX WINDOW.  When >= 0, only the cells i in [iwl,iwu]
  // are packed/unpacked, and ONLY for the buffers whose i-range spans the whole x1
  // extent of the block (the x2/x3 transverse halos), which are exactly the buffers
  // whose sender and receiver index ranges are the SAME absolute x1 indices -- so both
  // sides skip the identical buffer entries and the layout is untouched.  Buffer
  // entries outside the window keep whatever they held; the caller must not read the
  // ghosts they feed.  Used by the RKL1 transverse-conduction substages, where the
  // stencil only touches a shrinking radial window of planes.
  // vmask (default 0 = every variable): pack/send/unpack (and corner-fill) only the
  // variables whose bit is set, in increasing order, as a COMPACT block of variables in
  // the buffer (selected variable number vl takes chunk vl; with every variable vl = v
  // and the layout is the usual one).  The other variables' ghosts are left as they
  // are.  Uniform grids only (multilevel: FATAL).  Used by <mhd>/rkg_halo_ien_only.
  TaskStatus PackAndSendCC(DvceArray5D<Real> &a, DvceArray5D<Real> &ca,
                           int iwl = -1, int iwu = -1, int vmask = 0);
  TaskStatus RecvAndUnpackCC(DvceArray5D<Real> &a, DvceArray5D<Real> &ca,
                             int iwl = -1, int iwu = -1, int vmask = 0);
  void FillPanelCornersCC(DvceArray5D<Real> &a, bool coarse = false, int vmask = 0);

  // COMPONENT-ROLE TABLE of this exchange.  The polar sign flip (sp) and the seam
  // basis transform (cs) act on the TANGENTIAL VECTOR PAIRS of the exchanged array and
  // leave every other slot alone.  By default the pair is (IVY, IVZ) = slots (2, 3),
  // whatever the array holds -- right for the fluid momenta/velocities, wrong for an
  // array whose slots 2,3 are scalars.  SetVectorPairs(nvar, {{a,b},...}) replaces the
  // default for THIS object: slot a is the x2 (theta / xi) component, slot b its x3
  // (phi / eta) partner, and every slot not listed is a SCALAR.  An empty list marks
  // every slot scalar.  Objects that never call it are bitwise unchanged.
  void SetVectorPairs(const int nvar, const std::vector<std::pair<int,int>> &pairs);
  bool vrole_set_ = false;
  // per slot: -1 = scalar; else 2*partner + (0 if the slot is the x2 member, 1 if x3)
  DvceArray1D<int> vrole_;
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
  // ucthshare-1102 <mhd>/cs_seam_emf_avg = false (diagnostic): RecvAndUnpackFluxFC leaves
  // every edge EMF as this block computed it (no same-level sum/average, no cube-vertex
  // average).  Single-block-per-panel cubed sphere only (the gate of the shared seam EMF).
  bool skip_emf_reconcile = false;

  //functions
  void InitSendIndices(MeshBoundaryBuffer &b,int o1,int o2,int o3,int f1,int f2) override;
  void InitRecvIndices(MeshBoundaryBuffer &b,int o1,int o2,int o3,int f1,int f2) override;
  TaskStatus InitFluxRecv(const int nvar) override;

  TaskStatus PackAndSendFC(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);
  TaskStatus RecvAndUnpackFC(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);
  void FillCoarseInBndryFC(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);
  void FillPanelCornersFC(DvceFaceFld4D<Real> &b, bool coarse = false);
  void ProlongateFC(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);

  TaskStatus PackAndSendFluxFC(DvceEdgeFld4D<Real> &flx);
  TaskStatus RecvAndUnpackFluxFC(DvceEdgeFld4D<Real> &flx);
  void SumBoundaryFluxes(DvceEdgeFld4D<Real> &flx, const bool same_level,
                         DvceArray2D<int> &nflx);
  void ZeroFluxesAtBoundaryWithFiner(DvceEdgeFld4D<Real> &flx, DvceArray2D<int> &nflx);
  void AverageBoundaryFluxes(DvceEdgeFld4D<Real> &flx, DvceArray2D<int> &nflx);
  void SavePanelCornerEMF(DvceEdgeFld4D<Real> &flx, DvceArray3D<Real> &save);
  void AveragePanelCornerEMF(DvceEdgeFld4D<Real> &flx, DvceArray3D<Real> &save);

  // CUBED-SPHERE SEAM GEOMETRY TABLES.  Everything the seam packer computes from
  // (xi,eta) alone -- the field transform's geometry, the along-seam resample's stencil
  // and the shear term's sigma -- is precomputed once per (block, destination panel,
  // staggering, kk, jj) and read back per cell.  See BuildCubedSphereSeamTables().
  bool cs_seam_tbl_ready = false;
  int cs_seam_nslot = 0;
  DualArray2D<int> cs_seam_slot;    // (nmb,6): destination panel -> slot, -1 if none
  DvceArray6D<Real> cs_seam_geom;   // (nmb,nslot,nkind,nk3,nk2,15): SeamFieldXform
  DvceArray6D<Real> cs_seam_stnc;   // (nmb,nkind,2,nk3,nk2,7): bs,wm,w0,wp,fpos,sden,sig
  DvceArray6D<Real> cs_seam_stnc_x; // layoutbug-1103: the same, cs_seam_blockend = extend
  DvceArray2D<int> cs_ext_jk;       // layoutbug-1103: (nmb,4) interior-end widening
  // seamexact-1103: the partner-component interpolants of the extended pass at the
  // sender's along-seam ghost positions, taken from the same-panel neighbour that owns
  // them (CsSeamExtPartner).  Only with > 2 MeshBlocks per panel side (cs_ext_q_), where
  // a block end lies off the panel midline: there the 4-point partner interpolation of
  // the outermost along-seam ghost would need one face beyond this block's ghosts.
  // cs_ext_q (nmb,4,n3,n2,n1): 0 = x3f at x2-face centres (slot j = face j), 1 = the
  // same with slot j = face j-1, 2 = x2f at x3-face centres (slot k = face k), 3 = the
  // same with slot k = face k-1 (the -1 slots carry the high-side outermost ghost face).
  bool cs_ext_q_ = false;
  DvceArray5D<Real> cs_ext_q;
  MeshBoundaryValuesCC *pbq_ = nullptr;
  void CsSeamExtPartner(DvceFaceFld4D<Real> &b);
  void BuildCubedSphereSeamTables();

  // halofc-1102: FLAT, LOAD-BALANCED cubed-sphere pack.  Every (block, buffer,
  // component) triple's index ranges, destination and seam parameters -- the per-team
  // setup of the old team kernel, which only depends on the (static) neighbour table --
  // are resolved ONCE into fc_seg, and the samples are walked by a RangePolicy in chunks
  // of FC_PACK_CHUNK, one thread per sample.  See PackAndSendFCFlat in bvals_fc.cpp.
  bool fc_seg_ready = false;
  int fc_nchunk_plain = 0, fc_nchunk_seam = 0;
  DvceArray1D<PackSegFC> fc_seg;          // the active (m,n,v) triples
  DvceArray1D<int> fc_chunk_seg;          // chunk -> fc_seg index (plain chunks first)
  DvceArray1D<int> fc_chunk_beg;          // chunk -> first flat sample of that triple
  // r1mspeed-1102: two walk sets over fc_seg, [0] all ng layers, [1] depth 1 (only the
  // source layer adjacent to each boundary; <mhd>/rkg_halo_d1, see resistivity.hpp).
  // Each set: SubSegFC boxes and chunk lists for the plain samples, the seam samples and
  // the seam srcval precompute (all chunk lists index fc_sub[set]).
  DvceArray1D<SubSegFC> fc_sub[2];
  DvceArray1D<int> fc_cs[2], fc_cb[2];    // chunk -> (sub index, first sample)
  int fc_np[2] = {0, 0}, fc_ns[2] = {0, 0}, fc_npre[2] = {0, 0};   // chunk counts
  DvceArray1D<Real> fc_preval;            // srcval scratch (max over the sets)
  // 1 = the next PackAndSendFC packs (and RecvAndUnpackFC unpacks) depth 1 only; reset
  // by the caller.  0 everywhere else.
  int pack_depth = 0;
  int last_pack_depth_ = 0;   // layoutbug-1103: pack_depth of the last pack
  // r1mspeed-1102: RecvAndUnpackFluxFC work arrays, kept (they were allocated, memset and
  // freed on every call: 2 cudaMalloc/cudaFree pairs and 4 fences per E exchange)
  DvceArray2D<int> nflx_buf;
  DvceArray3D<Real> cnr_e1_buf;
  void BuildPackSegmentsFC();
  void PackAndSendFCFlat(DvceFaceFld4D<Real> &b, DvceFaceFld4D<Real> &cb);
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
