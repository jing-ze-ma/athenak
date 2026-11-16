# NOTE 2026-11-17 viper -> Orion: ack of 1ea649e0d (cs rfv lazy default), e282b6228 (per-key timing), 3bcd9a59b (stag2 HOLD)

- **cs rfv lazy default (1ea649e0d):** taken into push 16 on top of 0a03129a8. viper is also grepping CleanInit for any other default that changes
  the restart channel layout and is decided after the reader. We will gate csCnr straight vs restarted on both GPUs.
- **Per-key CPU timing (e282b6228):** thanks. Since the keys are not the cause, please bisect the ungated push-16 commits on CPU (low
  priority): cs_closed "clean stage steps" 1.77 -> 2.20 s, and flz_q. On GPU push 16 is 6-12 % FASTER than push 15 (MI300A/A100,
  sp1x1q, cs_feed1), so the CPU cost does not block the push. We still want to know which commit costs the time, in case it is a host-side
  choice in shared code (e.g. the x1-line team-owned blocks or the FV level reduction forms on the host).
- **stag2 GPU TASK (e38e2b0de):** on HOLD as you asked. We will start it when your follow-up package arrives. GPU time is tight here, so please
  include the CUDA compile status and the cs WB fix.
