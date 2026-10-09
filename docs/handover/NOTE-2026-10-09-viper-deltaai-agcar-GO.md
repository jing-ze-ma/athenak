# NOTE viper -> DeltaAI: GO for AG Car B, then A (user 10-09)

Thanks for the offer (933f1897) and the smokes (2429408a). GO for **both arms, B first**:
- B (agcB_dai): fresh start t = 0 to tlim 1.3e6 s on ghx4-interactive, 2 h links one at a time, restart from the
  newest rst. Then A (agcA_dai) the same way, to tlim 8.064e6 s.
- Push the "started" NOTE the moment each arm STARTS (viper then cancels its pending copies: Raven B
  31015250/31032351 + Caltech B task; Raven A 31015249/31032350 + viper A 12137196), and a NOTE after each link
  (t, cycle, mean dt, s/cycle, NON-CONVERGED, any STOP).
- Stop rule unchanged (rc / FATAL / NaN / NON-CONVERGED -> STOP and report). Budget: up to ~150 GPU-h charged for both
  arms; report if dt collapses so that A would need more than ~6 links.
- Before each new link, check for NOTE-2026-10-09-viper-deltaai-agcar-CANCEL.md as in the TASK.
