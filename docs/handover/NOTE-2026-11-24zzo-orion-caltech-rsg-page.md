# NOTE 2026-11-24zzo orion -> Caltech: USER: update the RSG page (URL) with the current RSG run

The USER asks to update the RSG page. That's the RSG artifact page Caltech maintains ("RSG E Caltech", https://claude.ai/artifact/LBASg2iK2ArcLY1V2yBCAR, per TASK-2026-11-07 page builders; or whichever RSG page you now publish).
Please:
1. Rebuild it from the CURRENT RSG run on 7333072: the clean scheme on rc3, continued from old-scheme rsg.01707 on 11-24, now c ~793500, t ~7.43e8 s. Include the continuity across the switch: L_top / L_bot / L_in and R_ph vs time across t of rsg.01707. Use L_bot comoving (f6d4e119a's meaning) and note the 1.2 % frame term in the old column. Keep the 3-D / slice views the page had.
2. Republish to the SAME URL (update in place). If the page can't be updated there, publish a new one.
3. Post the URL (and what changed) in a NOTE, so I can give it to the USER.
4. Low priority vs production: use idle CPU (srun --overlap --gres=none inside a holder, or RSG's own post-processing); stop if production s/cyc slows > 2 %.

-- orion (coordinator)
