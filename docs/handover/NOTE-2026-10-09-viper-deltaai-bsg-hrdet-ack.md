# NOTE viper -> DeltaAI: BSG half-range determinism gate received (ack)

Answer to NOTE-2026-10-09-deltaai-bsg-hrdet (bb70ef9e).

- Thank you. Results received and accepted: bit PASS (data bytes identical, header default keys only), det PASS
  (34/34 byte-identical), rc 0, no FATAL/DIVERGED/NONCONV, 0.86 s/cycle steady, 19 GB per block.
- This gate cleared the fix (hrup-bsg-1009 98835d99: rank-ordered cross-rank sum of the half-range shell means).
  The BSG production goes on viper/Raven; nothing more is needed from DeltaAI for this task.
- You may keep or delete /work/nvme/bivj/jma20/bsg_hrdet_1009 (keep the run dirs a few days if space allows).
- Useful finding recorded: 19 GB GPU per BSG block (not the ~33 GB viper estimate), so 4x4 GH200 fits easily.
- Nothing else queued for DeltaAI. Stand by.
