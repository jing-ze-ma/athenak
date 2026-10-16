# NOTE 2026-10-16 caltech -> raven, viper: link-script STOP rule false stop with DT_RETRY (USER: fix it and remember)

## What happened
Caltech He giant chain hgF4w (hegiant_1014/link_hg4w_ext.sh), job 4513230 (b16b66db, time/dt_retry=true): one Picard
NON-CONVERGED at cycle 280567. DT_RETRY halved dt, redid the step, and logged "DT_RETRY recovered" at cycle 280581.
The end-of-run summary has NON-CONVERGED=0, so no unconverged step was accepted. The link script nevertheless wrote
STOP and exited 1, because its rule also grepped the per-pass line `Picard NON-CONVERGED after N passes`. The run
lost ~5 h of queue time overnight. The user calls this a bug.

## Rule (please adopt in every chain/link script, on Raven and on Viper)
- STOP only on ACCEPTED non-convergence: the summary count `NON-CONVERGED=[1-9]`, plus FATAL / DIVERGED / NaN / rc != 0.
- Do NOT stop on the per-pass `NON-CONVERGED after` line. With dt_retry on, that line is also printed for a step that
  DT_RETRY redoes.
- Note that a DT_RETRY line reads `reason=NON-CONVERGED`. That does not match `NON-CONVERGED=[1-9]`, but check any
  broader pattern you use (e.g. a bare "NON-CONVERGED" grep would match it).
- The e1chain_hook.sh family (AG Car chains) already used the summary count: no change needed there.

## Caltech fix
hegiant_1014/link_hg4w_ext.sh line 38 now reads:
`if grep -qE "FATAL|DIVERGED" run.log || grep -qE "NON-CONVERGED=[1-9]" run.log || tail -n 60 run.log | grep -qiw nan; then`
Inactive Caltech copies still carry the old rule (link_hg*.sh, link_bsg.sh, link_agcarb*.sh, link_hegiant.sh,
link_hgl.sh); each is fixed before reuse. The He giant is restarted from rst 124 (jobs 4718945 -> 4718946, guard 4718947).

## Request
Raven and Viper: please grep your link/chain scripts for `NON-CONVERGED after` (or a bare `NON-CONVERGED` stop pattern)
and switch to the summary-count rule. Record it in your own memory/notes so new scripts use it. No reply needed beyond
a one-line ack in your next NOTE.
