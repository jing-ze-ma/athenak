# NOTE caltech -> raven: please check the BSG radial flux F_r in 2-D on viper (USER request 10-15)

USER (10-15): "ask raven to check also on viper the 2d plot for Fr in bsg".

## Why
Our 3-D movies colour the tau_R = 2/3 photosphere by T_eff = (F_r(r_ph)/sigma)^1/4 (local radial flux at the photosphere,
F <= 0 clipped to a floor). On the AG Car LBV run this map is speckled/blocky: F_r at r_ph is noisy and partly <= 0 in the
convective, inflated envelope. Caltech now applies one 3x3 box pass to F_r r^2 at r_ph for AG Car only. On the BSG
(Linhao IC) run the same unsmoothed map looks clean. We want an independent look at F_r on viper to confirm whether the BSG
flux near the photosphere is genuinely smooth, and whether some smoothing or a different T_eff definition should be used for all runs.

## Ask (diagnostic only, no runs needed)
1. From the viper BSG run(s) (latest dumps), plot 2-D slices of F_r (lab-frame radial radiation flux; also F_r r^2 / L_*-normalised)
   - (r, theta) at one phi and (theta, phi) at r = r_ph (tau_R = 2/3 per column), linear colour, diverging map around 0;
   - overlay the tau_R = 2/3 and tau_R = 10 contours.
2. Report: fraction of photosphere columns with F_r <= 0; p2/50/98 of F_r(r_ph) r^2 / (L_*/4 pi); radial range over which F_r
   changes sign; whether the convective (F_conv) share at tau ~ 1 is significant.
3. If viper has its own T_eff surface maps: how was T_eff defined (local F_r at r_ph, escaping flux at the top scaled by r^2,
   any smoothing)?
4. Optional: the same for the He giant and AG Car states if available on viper, for comparison.

Push the plots (PNG, < 45 MB each) and a short NOTE back to bsg-files-1009.

Caltech references: renderer /resnick/groups/carnegie_poc/jingze/vis3d/render3d.py (teff(), RUNS['bsg_li'] Teff 12.35-17.37 kK,
Tsrc gas); BSG run bsgli_1014/pagerun.
