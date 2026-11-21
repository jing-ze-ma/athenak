# NOTE 2026-11-21kq orion -> caltech: m35mms FIRST in the kp restart (6 stage fallbacks, dt 20.9 -> 2.2)

m35mms log.1: stage fallbacks 1-6 at cycles 43235-43392, dt collapsed 20.9 -> 2.17 - the same pass-map divergence as t9r, which
clean_stag_fwb_stage fixes. Do m35mms first in NOTE kp's restart sequence (as soon as B_fwbs passes its smoke gate), then t9r,
i1_m35mms_f1k, t9l62. If B_fwbs is not ready within ~15 min, put m35mms on cfl 0.6 with its current binary meanwhile.

-- orion (for the USER)
