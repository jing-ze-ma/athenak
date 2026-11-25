# NOTE 2026-11-26p orion -> Caltech: movies for the documentation (addendum to 26o, the USER's wish, 11:57 PST)

"or if you can put movies in the doc that's even better" (USER). So, AFTER the palette pick, for the figures of 26o also one short movie per star (bsg, m35mms, m35zams, t9l62, t9mlt) and, if it costs little, one montage:
- Frames: the photosphere map (v_r or T) and the slice, from existing dumps at their cadence (m35mms 1800 s, others as they are); no new runs. 8-12 s at 24-30 fps; a time stamp in d; the same colormap and limits as the figures, limits fixed over the whole movie; no deep stable layer velocity (as 26o).
- Format: H.264 MP4 (yuv420p, faststart) + a poster PNG of one frame; width 960-1280; target <= 4 MB each, hard limit 8 MB (the docs repository must stay small); a <= 2 MB WebM only if your ffmpeg has it. No audio.
- State per movie: time range, number of dumps, wall time of the run on one node of 4 H200, and the caption line. Start-up transients may be shortened but say so ("first 2 d not shown").
- Deliver to docs/handover/caltech_results/DOCFIG/ as small files; if the MP4 exceeds the git comfort size, tell me the total first. The palette sheet (26o item 1) comes first.

-- orion (coordinator)
