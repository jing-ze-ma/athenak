// minimal Kokkos tool: fenced wall time per kernel name (parallel_for/reduce/scan) and
// per deep_copy, written at finalize to $KP_OUT.<rank> (rank from OMPI/PMI/SLURM env)
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
struct Kokkos_Tools_ToolSettings { bool requires_global_fencing; bool padding[255]; };
namespace {
using clk = std::chrono::steady_clock;
struct Rec { double t = 0.0; uint64_t n = 0; };
std::map<std::string, Rec> g_rec;
std::vector<std::pair<std::string, clk::time_point>> g_live;   // indexed by kID
clk::time_point g_t0;
uint64_t begin(const char *name) {
  g_live.emplace_back(name, clk::now());
  return g_live.size() - 1;
}
void end(uint64_t id) {
  auto &e = g_live[id];
  Rec &r = g_rec[e.first];
  r.t += std::chrono::duration<double>(clk::now() - e.second).count();
  r.n += 1;
  if (id + 1 == g_live.size()) g_live.pop_back();
}
}  // namespace
extern "C" {
void kokkosp_request_tool_settings(const uint32_t, Kokkos_Tools_ToolSettings *s) {
  s->requires_global_fencing = true;
}
void kokkosp_init_library(const int, const uint64_t, const uint32_t, void *) {
  g_t0 = clk::now();
  g_live.reserve(1024);
}
void kokkosp_begin_parallel_for(const char *n, const uint32_t, uint64_t *k) {*k = begin(n);}
void kokkosp_end_parallel_for(const uint64_t k) {end(k);}
void kokkosp_begin_parallel_reduce(const char *n, const uint32_t, uint64_t *k) {*k = begin(n);}
void kokkosp_end_parallel_reduce(const uint64_t k) {end(k);}
void kokkosp_begin_parallel_scan(const char *n, const uint32_t, uint64_t *k) {*k = begin(n);}
void kokkosp_end_parallel_scan(const uint64_t k) {end(k);}
static clk::time_point g_dc;
void kokkosp_begin_deep_copy(void *, const char *, const void *, void *, const char *,
                             const void *, uint64_t) {g_dc = clk::now();}
void kokkosp_end_deep_copy() {
  Rec &r = g_rec["<deep_copy>"];
  r.t += std::chrono::duration<double>(clk::now() - g_dc).count();
  r.n += 1;
}
void kokkosp_finalize_library() {
  const char *rk = nullptr;
  for (const char *v : {"OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK", "SLURM_PROCID"}) {
    if ((rk = std::getenv(v)) != nullptr) break;
  }
  const char *out = std::getenv("KP_OUT");
  std::string fn = std::string(out ? out : "kp") + "." + (rk ? rk : "0");
  FILE *f = std::fopen(fn.c_str(), "w");
  if (!f) return;
  std::vector<std::pair<std::string, Rec>> v(g_rec.begin(), g_rec.end());
  std::sort(v.begin(), v.end(), [](auto &a, auto &b) {return a.second.t > b.second.t;});
  std::fprintf(f, "# total %.6f s\n", std::chrono::duration<double>(clk::now() - g_t0).count());
  for (auto &e : v) std::fprintf(f, "%14.6f %10llu %s\n", e.second.t,
                                 (unsigned long long)e.second.n, e.first.c_str());
  std::fclose(f);
}
}
