# NOTE 2026-11-21km orion -> caltech: USER decision — rotating-ray EMA change_form = linear (HG/RSG and any rotating run)

USER (2026-11-21): rad_m1/vet_rotate_avg_ema_change_form = linear, as the default for rotating ray sets AND in the HG/RSG production keys.
Why: the v20.4 default "step" (adopted 11-20 only to match the K1 restart files) breaks the EMA ring's 2nd order in time
(sp ladder E/F/T order 1.66/1.39/1.63 with step vs 1.95/1.85/1.96 with linear; fixed rays 2.0). Orion key files hg_latest.txt and
rsg_latest.txt now name linear; code default flip = commit 28c7c7766 (branch clean-stag2-lean-v207-chlin, not yet in a release).
Action: at the NEXT (re)start of any rotating-ray HG/RSG/R136 run, name rad_m1/vet_rotate_avg_ema_change_form=linear on the command line.
A restart from a step file switching to linear is a one-time change of the change-weight form (clamp(c/eps-1,0,1) -> min(1,c/eps));
no other key changes. Do not restart running production just for this. (A cs order confirmation arm is finishing on Orion.)

-- orion (for the USER)
