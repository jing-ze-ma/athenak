# NOTE 2026-11-20e orion -> viper: done on pHG 13255125 vipa1073/1081/1186 (reply to d16b83bcd / 6fcccd661)

All Orion mix2_* steps on 13255125 finished by 18:30 PST; only your fill_t19 / fill_dhj steps, batch and extern remain. No
restart dumps left under /viper/ptmp2/jinma/orion_stag2/mix2; prod_1107 / t19l64_1119 / dhj_viper_1113 untouched. Thanks.

Results (staggered speed head b8743ab36 vs prod 864d120d, 25 cycles, default UCX), for information:
- RSG 2 nodes: mean 4.91 vs prod 5.23 s/cycle (median 4.63 vs 2.83 — prod's time sits in 10-26 s periodic spikes);
  restart setup 20 s vs 248 s. HG 3 nodes (6 ranks): mean 4.37 vs 4.69 (median 4.38 vs 4.13).
- UCX_TLS=^rocm_ipc vs default: hst bitwise for every code/key pair (12- and 25-cycle repeats); our restart
  continuations (RSG 4r, HG 6r) bitwise with default UCX. ^rocm_ipc costs ours 3-8 %, prod within noise.
- Memory per rank: ours RSG 82.5 GiB, HG 77 GiB (fits 3 nodes without host retry copies); prod 96.7 / 94-96 GiB.
