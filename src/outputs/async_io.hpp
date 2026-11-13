#ifndef OUTPUTS_ASYNC_IO_HPP_
#define OUTPUTS_ASYNC_IO_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file async_io.hpp
//! \brief asyncio-1113: background writing of restart and binary output files.
//!
//! Opt-in by the environment variable ATHENA_ASYNC_OUTPUT=1 (no input key, so that the
//! parameter dump inside a restart file stays byte-identical).  With it set, an IOWrapper
//! marked SetDeferred() records every write (offset, bytes) instead of doing it; Close()
//! hands the record to ONE helper thread, which writes the bytes with POSIX pwrite at the
//! same offsets MPI-IO / fwrite would have used (same file contents, byte for byte) into
//! "<name>.part" and renames it to <name> once every writer of the file has finished:
//!  - one file per rank (or no MPI): the helper renames its own file;
//!  - one shared file (MPI-IO path): each rank's main thread notices its helper is done
//!    (Progress(), once per cycle) and joins a nonblocking MPI_Iallreduce of the error
//!    flags; when it completes, rank 0 renames.  The helper thread never calls MPI.
//! A crash mid-write therefore leaves only a .part file; the previous restart is intact.
//! The helper writes from host memory that the caller keeps unchanged until the job is
//! done: small writes (headers, < 1 MiB) are copied when recorded; large writes reference
//! the caller's host staging arrays (restart: the RestartOutput outarray_* members, which
//! RestartOutput::LoadOutputData refills only after WaitTag(this)) or a buffer the job
//! adopts (binary output: AdoptBuffer).

#include <string>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace async_io {

struct Job;

bool Enabled();  // ATHENA_ASYNC_OUTPUT=1 (read once)

#if MPI_PARALLEL_ENABLED
// called by IOWrapper::Open on the main thread of every rank of comm (collective when
// the file is shared)
Job *Open(const void *tag, const char *fname, bool shared, MPI_Comm comm);
#else
Job *Open(const void *tag, const char *fname, bool shared);
#endif
// record a write of n bytes at offset off; copy=true copies the bytes now
void Record(Job *job, const void *p, std::uint64_t off, std::uint64_t n, bool copy);
// the job frees p (new char[]) after writing it
void Adopt(Job *job, char *p);
// hand the job to the helper thread (IOWrapper::Close)
void Submit(Job *job);

void WaitTag(const void *tag);  // until every job of tag has been written locally
void Progress();                // nonblocking: finish (rename) completed jobs, in order
void WaitAll();                 // blocking: everything written and renamed

}  // namespace async_io

#endif  // OUTPUTS_ASYNC_IO_HPP_
