//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file async_io.cpp
//! \brief asyncio-1113: background writing of restart and binary output files (see
//! async_io.hpp).  Opt-in by ATHENA_ASYNC_OUTPUT=1.

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <atomic>
#include <cerrno>
#include <chrono>  // NOLINT(build/c++11)
#include <condition_variable>  // NOLINT(build/c++11)
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <memory>
#include <mutex>  // NOLINT(build/c++11)
#include <string>
#include <thread>  // NOLINT(build/c++11)
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "async_io.hpp"

namespace async_io {

using Clock = std::chrono::steady_clock;

struct Op {
  std::uint64_t off;
  const char *p;
  std::uint64_t n;
};

struct Job {
  const void *tag = nullptr;
  std::string fname, tname;
  bool shared = false;   // one file written by every rank of comm (the MPI-IO path)
  bool root = true;      // the rank that renames a shared file
  int fd = -1;
  std::vector<Op> ops;
  std::vector<std::unique_ptr<char[]>> owned;
  std::uint64_t bytes = 0;
  int err = 0;              // errno of the first failure (helper thread)
  std::string errwhat;
  std::atomic<bool> written{false};
  Clock::time_point t_submit;
  double secs = 0.0;        // submit -> written (helper thread)
  double wsecs = 0.0;       // time spent writing (helper thread)
#if MPI_PARALLEL_ENABLED
  MPI_Comm comm = MPI_COMM_NULL;
  bool posted = false;
  int sflag = 0, rflag = 0;
  MPI_Request req = MPI_REQUEST_NULL;
#endif
};

namespace {

// all state on the heap and never freed: a fatal std::exit with the helper running must
// not destroy a mutex or a joinable std::thread under it
struct State {
  std::mutex mu;
  std::condition_variable cv_todo, cv_done;
  std::deque<Job*> todo;      // submitted, not yet taken by the helper
  std::deque<Job*> pending;   // main thread only: submitted, not yet finished (renamed)
  bool started = false;
};
State *S = nullptr;

State *GetState() {
  if (S == nullptr) S = new State;
  return S;
}

void Fail(Job *j, const char *what) {
  if (j->err == 0) {
    j->err = (errno != 0) ? errno : EIO;
    j->errwhat = what;
  }
}

void HelperLoop() {
  State *s = S;
  for (;;) {
    Job *j;
    {
      std::unique_lock<std::mutex> lk(s->mu);
      s->cv_todo.wait(lk, [s] { return !s->todo.empty(); });
      j = s->todo.front();
      s->todo.pop_front();
    }
    auto t0 = Clock::now();
    for (const Op &op : j->ops) {
      std::uint64_t done = 0;
      while (done < op.n && j->err == 0) {
        std::uint64_t chunk = op.n - done;
        if (chunk > (std::uint64_t(1) << 30)) chunk = (std::uint64_t(1) << 30);
        ssize_t w = pwrite(j->fd, op.p + done, static_cast<size_t>(chunk),
                           static_cast<off_t>(op.off + done));
        if (w < 0) {
          if (errno == EINTR) continue;
          Fail(j, "pwrite");
        } else if (w == 0) {
          errno = EIO;
          Fail(j, "pwrite");
        } else {
          done += static_cast<std::uint64_t>(w);
        }
      }
      if (j->err != 0) break;
    }
    if (close(j->fd) != 0) Fail(j, "close");
    j->fd = -1;
    j->ops.clear();
    j->owned.clear();   // release the adopted / copied buffers now
    if (!j->shared && j->err == 0) {
      if (std::rename(j->tname.c_str(), j->fname.c_str()) != 0) Fail(j, "rename");
    }
    auto t1 = Clock::now();
    j->wsecs = std::chrono::duration<double>(t1 - t0).count();
    j->secs = std::chrono::duration<double>(t1 - j->t_submit).count();
    {
      std::lock_guard<std::mutex> lk(s->mu);
      j->written.store(true);
    }
    s->cv_done.notify_all();
  }
}

// main thread: the job is written locally (and, if shared, globally); rename + report
void Finish(Job *j, int global_err) {
  if (j->shared && j->root && global_err == 0) {
    if (std::rename(j->tname.c_str(), j->fname.c_str()) != 0) {
      std::cout << "### ERROR async_output: rename " << j->tname << " -> " << j->fname
                << " failed: " << std::strerror(errno) << std::endl;
    }
  }
  if (j->err != 0) {
    std::cout << "### ERROR async_output (rank " << global_variable::my_rank << "): "
              << j->errwhat << " on " << j->tname << " failed: "
              << std::strerror(j->err) << "; the file is left as " << j->tname
              << std::endl;
  } else if (global_err != 0 && j->root) {
    std::cout << "### ERROR async_output: a rank failed writing " << j->tname
              << "; not renamed" << std::endl;
  }
  if (global_variable::my_rank == 0 && j->err == 0 && global_err == 0) {
    std::cout << "async_output: " << j->fname << " written in background, "
              << j->bytes << " B on rank 0, write " << j->wsecs << " s, submit->done "
              << j->secs << " s" << std::endl;
  }
#if MPI_PARALLEL_ENABLED
  if (j->comm != MPI_COMM_NULL) MPI_Comm_free(&(j->comm));
#endif
  delete j;
}

void WaitWritten(Job *j) {
  std::unique_lock<std::mutex> lk(S->mu);
  S->cv_done.wait(lk, [j] { return j->written.load(); });
}

}  // namespace

bool Enabled() {
  static int en = -1;
  if (en < 0) {
    const char *e = std::getenv("ATHENA_ASYNC_OUTPUT");
    en = (e != nullptr && std::atoi(e) != 0) ? 1 : 0;
  }
  return en == 1;
}

#if MPI_PARALLEL_ENABLED
Job *Open(const void *tag, const char *fname, bool shared, MPI_Comm comm) {
#else
Job *Open(const void *tag, const char *fname, bool shared) {
#endif
  Job *j = new Job;
  j->tag = tag;
  j->fname = fname;
  j->tname = j->fname + ".part";
  j->shared = shared;
  bool create = true;
#if MPI_PARALLEL_ENABLED
  if (shared) {
    int crank = 0;
    MPI_Comm_rank(comm, &crank);
    MPI_Comm_dup(comm, &(j->comm));
    j->root = (crank == 0);
    create = j->root;
  }
#endif
  if (create) {
    unlink(j->tname.c_str());
    j->fd = open(j->tname.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
  }
#if MPI_PARALLEL_ENABLED
  if (shared) {
    int ok = (j->fd >= 0 || !create) ? 1 : 0, allok = 0;
    MPI_Allreduce(&ok, &allok, 1, MPI_INT, MPI_MIN, j->comm);  // also orders the create
    if (allok == 1 && !create) j->fd = open(j->tname.c_str(), O_WRONLY);
  }
#endif
  if (j->fd < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "File '" << j->tname << "' could not be opened: "
              << std::strerror(errno) << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return j;
}

void Record(Job *j, const void *p, std::uint64_t off, std::uint64_t n, bool copy) {
  if (n == 0) return;
  const char *src = static_cast<const char*>(p);
  if (copy) {
    std::unique_ptr<char[]> b(new char[n]);
    std::memcpy(b.get(), src, n);
    src = b.get();
    j->owned.push_back(std::move(b));
  }
  j->ops.push_back(Op{off, src, n});
  j->bytes += n;
}

void Adopt(Job *j, char *p) {
  j->owned.emplace_back(p);
}

void Submit(Job *j) {
  State *s = GetState();
  j->t_submit = Clock::now();
  s->pending.push_back(j);
  {
    std::lock_guard<std::mutex> lk(s->mu);
    s->todo.push_back(j);
    if (!s->started) {
      s->started = true;
      std::thread *t = new std::thread(HelperLoop);  // never joined, never deleted
      t->detach();
    }
  }
  s->cv_todo.notify_one();
}

void Progress() {
  if (S == nullptr) return;
  while (!S->pending.empty()) {
    Job *j = S->pending.front();
    if (!j->written.load()) break;
    int gerr = j->err;
#if MPI_PARALLEL_ENABLED
    if (j->shared) {
      if (!j->posted) {
        j->sflag = (j->err != 0) ? 1 : 0;
        MPI_Iallreduce(&(j->sflag), &(j->rflag), 1, MPI_INT, MPI_MAX, j->comm,
                       &(j->req));
        j->posted = true;
      }
      int flag = 0;
      MPI_Test(&(j->req), &flag, MPI_STATUS_IGNORE);
      if (!flag) break;
      gerr = j->rflag;
    }
#endif
    S->pending.pop_front();
    Finish(j, gerr);
  }
}

void WaitTag(const void *tag) {
  if (S == nullptr) return;
  for (Job *j : S->pending) {
    if (j->tag == tag) WaitWritten(j);
  }
  Progress();
}

void WaitAll() {
  if (S == nullptr) return;
  while (!S->pending.empty()) {
    Job *j = S->pending.front();
    WaitWritten(j);
    int gerr = j->err;
#if MPI_PARALLEL_ENABLED
    if (j->shared) {
      if (!j->posted) {
        j->sflag = (j->err != 0) ? 1 : 0;
        MPI_Iallreduce(&(j->sflag), &(j->rflag), 1, MPI_INT, MPI_MAX, j->comm,
                       &(j->req));
        j->posted = true;
      }
      MPI_Wait(&(j->req), MPI_STATUS_IGNORE);
      gerr = j->rflag;
    }
#endif
    S->pending.pop_front();
    Finish(j, gerr);
  }
}

}  // namespace async_io
