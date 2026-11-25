// kt.cpp: minimal Kokkos kernel timer (kp_kernel_timer-like). Global fencing on, so each
// kernel's wall = its device time. KT_OUT=<prefix> writes <prefix>.<pid>.txt at finalize:
// total_s count kind name, sorted by time. Also times deep_copy and fences.
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>
struct Kokkos_Tools_ToolSettings { bool requires_global_fencing; bool padding[255]; };
struct KokkosPDeviceInfo { uint32_t deviceID; };
struct SpaceHandle { char name[64]; };
using clk = std::chrono::steady_clock;
struct Acc { double t = 0; long n = 0; };
static std::map<std::string, Acc> acc;
static std::vector<std::pair<std::string, clk::time_point>> open_k;  // by kID index
static std::map<uint64_t, std::pair<std::string, clk::time_point>> live;
static uint64_t next_id = 1;
static clk::time_point t_init;
static uint64_t beg(const char* kind, const char* name) {
  uint64_t id = next_id++; live[id] = {std::string(kind) + " " + name, clk::now()}; return id; }
static void fin(uint64_t id) {
  auto it = live.find(id); if (it == live.end()) return;
  Acc& a = acc[it->second.first];
  a.t += std::chrono::duration<double>(clk::now() - it->second.second).count(); a.n++; live.erase(it); }
extern "C" {
void kokkosp_request_tool_settings(const uint32_t, Kokkos_Tools_ToolSettings* s) {
  s->requires_global_fencing = true; }
void kokkosp_init_library(const int, const uint64_t, const uint32_t, void*) { t_init = clk::now(); }
void kokkosp_begin_parallel_for(const char* n, uint32_t, uint64_t* k) { *k = beg("F", n); }
void kokkosp_end_parallel_for(uint64_t k) { fin(k); }
void kokkosp_begin_parallel_reduce(const char* n, uint32_t, uint64_t* k) { *k = beg("R", n); }
void kokkosp_end_parallel_reduce(uint64_t k) { fin(k); }
void kokkosp_begin_parallel_scan(const char* n, uint32_t, uint64_t* k) { *k = beg("S", n); }
void kokkosp_end_parallel_scan(uint64_t k) { fin(k); }
void kokkosp_begin_fence(const char* n, uint32_t, uint64_t* k) { *k = beg("X", n); }
void kokkosp_end_fence(uint64_t k) { fin(k); }
static uint64_t dc_id = 0;
void kokkosp_begin_deep_copy(SpaceHandle, const char* dn, const void*, SpaceHandle, const char* sn,
                             const void*, uint64_t) {
  dc_id = beg("C", (std::string(sn) + "->" + dn).c_str()); }
void kokkosp_end_deep_copy() { fin(dc_id); }
void kokkosp_finalize_library() {
  const char* p = getenv("KT_OUT"); std::string f = std::string(p ? p : "kt") + "." +
      std::to_string(getpid()) + ".txt";
  std::vector<std::pair<std::string, Acc>> v(acc.begin(), acc.end());
  std::sort(v.begin(), v.end(), [](auto& a, auto& b) { return a.second.t > b.second.t; });
  FILE* o = fopen(f.c_str(), "w");
  fprintf(o, "# wall %.3f s\n", std::chrono::duration<double>(clk::now() - t_init).count());
  for (auto& e : v) fprintf(o, "%12.4f %9ld %s\n", e.second.t, e.second.n, e.first.c_str());
  fclose(o); }
}
