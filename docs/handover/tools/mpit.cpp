// mpit.cpp: where the host time outside kernels goes. MPI calls (PMPI wrappers) and Kokkos fences / copies, each keyed
// by the LAST Kokkos kernel launched before it (= the call site). No fence is added (requires_global_fencing = false),
// no source change, results untouched. Build: mpicxx -O2 -shared -fPIC -o libmpit.so mpit.cpp
// Run:  srun ... env LD_PRELOAD=<dir>/libmpit.so KOKKOS_TOOLS_LIBS=<dir>/libmpit.so MT_OUT=<dir>/mt <athena> ...
// One text file <MT_OUT>.<rank>.txt per rank at MPI_Finalize: seconds, calls, bytes, kind, key.
#include <mpi.h>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
struct Kokkos_Tools_ToolSettings { bool requires_global_fencing; bool padding[255]; };
struct SpaceHandle { char name[64]; };
using clk = std::chrono::steady_clock;
struct Acc { double t = 0, b = 0; long n = 0; };
using Key = std::pair<std::string, std::string>;
static std::map<Key, Acc>& A() { static auto* a = new std::map<Key, Acc>; return *a; }
static std::string& last() { static auto* s = new std::string("(none)"); return *s; }
static clk::time_point t0 = clk::now();
static int in_copy = 0;   // (a copy's own fences count as the copy)
static void add(const char* k, double t, double b) { auto& r = A()[Key(k, last())]; r.t += t; r.b += b; r.n++; }
struct Tm { const char* k; double b; clk::time_point s = clk::now(); Tm(const char* k_, double b_ = 0) : k(k_), b(b_) {}
  ~Tm() { add(k, std::chrono::duration<double>(clk::now() - s).count(), b); } };
static double bytes(int n, MPI_Datatype d) { int z = 0; PMPI_Type_size(d, &z); return double(n)*z; }
struct FS { std::string name; clk::time_point s; };
static std::map<uint64_t, FS>& F() { static auto* f = new std::map<uint64_t, FS>; return *f; }
static uint64_t fid = 1;
static clk::time_point cs;
extern "C" {
void kokkosp_request_tool_settings(const uint32_t, Kokkos_Tools_ToolSettings* s) { s->requires_global_fencing = false; }
void kokkosp_init_library(const int, const uint64_t, const uint32_t, void*) {}
void kokkosp_finalize_library() {}
void kokkosp_begin_parallel_for(const char* n, uint32_t, uint64_t* k) { last() = n; *k = 0; }
void kokkosp_end_parallel_for(uint64_t) {}
void kokkosp_begin_parallel_reduce(const char* n, uint32_t, uint64_t* k) { last() = n; *k = 0; }
void kokkosp_end_parallel_reduce(uint64_t) {}
void kokkosp_begin_parallel_scan(const char* n, uint32_t, uint64_t* k) { last() = n; *k = 0; }
void kokkosp_end_parallel_scan(uint64_t) {}
void kokkosp_begin_fence(const char* n, uint32_t, uint64_t* k) { *k = fid++; F()[*k] = {n, clk::now()}; }
void kokkosp_end_fence(uint64_t k) { auto it = F().find(k); if (it == F().end()) return;
  if (in_copy) { F().erase(it); return; }
  std::string nm = "fence: " + it->second.name.substr(0, 50);
  add(nm.c_str(), std::chrono::duration<double>(clk::now() - it->second.s).count(), 0); F().erase(it); }
void kokkosp_begin_deep_copy(SpaceHandle, const char*, const void*, SpaceHandle, const char*, const void*, uint64_t) { in_copy = 1; cs = clk::now(); }
void kokkosp_end_deep_copy() { in_copy = 0; add("deep_copy", std::chrono::duration<double>(clk::now() - cs).count(), 0); }
// MPI
int MPI_Allreduce(const void* s, void* r, int n, MPI_Datatype d, MPI_Op o, MPI_Comm c) { Tm t("Allreduce", bytes(n, d)); return PMPI_Allreduce(s, r, n, d, o, c); }
int MPI_Iallreduce(const void* s, void* r, int n, MPI_Datatype d, MPI_Op o, MPI_Comm c, MPI_Request* q) { Tm t("Iallreduce", bytes(n, d)); return PMPI_Iallreduce(s, r, n, d, o, c, q); }
int MPI_Reduce(const void* s, void* r, int n, MPI_Datatype d, MPI_Op o, int R, MPI_Comm c) { Tm t("Reduce", bytes(n, d)); return PMPI_Reduce(s, r, n, d, o, R, c); }
int MPI_Barrier(MPI_Comm c) { Tm t("Barrier"); return PMPI_Barrier(c); }
int MPI_Bcast(void* b, int n, MPI_Datatype d, int R, MPI_Comm c) { Tm t("Bcast", bytes(n, d)); return PMPI_Bcast(b, n, d, R, c); }
int MPI_Waitall(int n, MPI_Request* q, MPI_Status* s) { Tm t("Waitall"); return PMPI_Waitall(n, q, s); }
int MPI_Wait(MPI_Request* q, MPI_Status* s) { Tm t("Wait"); return PMPI_Wait(q, s); }
int MPI_Waitany(int n, MPI_Request* q, int* i, MPI_Status* s) { Tm t("Waitany"); return PMPI_Waitany(n, q, i, s); }
int MPI_Waitsome(int n, MPI_Request* q, int* m, int* i, MPI_Status* s) { Tm t("Waitsome"); return PMPI_Waitsome(n, q, m, i, s); }
int MPI_Test(MPI_Request* q, int* f, MPI_Status* s) { Tm t("Test"); return PMPI_Test(q, f, s); }
int MPI_Testall(int n, MPI_Request* q, int* f, MPI_Status* s) { Tm t("Testall"); return PMPI_Testall(n, q, f, s); }
int MPI_Isend(const void* b, int n, MPI_Datatype d, int p, int g, MPI_Comm c, MPI_Request* q) { Tm t("Isend", bytes(n, d)); return PMPI_Isend(b, n, d, p, g, c, q); }
int MPI_Irecv(void* b, int n, MPI_Datatype d, int p, int g, MPI_Comm c, MPI_Request* q) { Tm t("Irecv", bytes(n, d)); return PMPI_Irecv(b, n, d, p, g, c, q); }
int MPI_Start(MPI_Request* q) { Tm t("Start"); return PMPI_Start(q); }
int MPI_Startall(int n, MPI_Request* q) { Tm t("Startall"); return PMPI_Startall(n, q); }
int MPI_Finalize() {
  int r = 0; PMPI_Comm_rank(MPI_COMM_WORLD, &r);
  const char* p = getenv("MT_OUT"); std::string f = std::string(p ? p : "mt") + "." + std::to_string(r) + ".txt";
  if (FILE* o = fopen(f.c_str(), "w")) {
    fprintf(o, "# wall %.4f s\n", std::chrono::duration<double>(clk::now() - t0).count());
    for (auto& e : A()) fprintf(o, "%14.6f %10ld %16.0f\t%s\t%s\n", e.second.t, e.second.n, e.second.b, e.first.first.c_str(), e.first.second.c_str());
    fclose(o); }
  return PMPI_Finalize(); }
}
