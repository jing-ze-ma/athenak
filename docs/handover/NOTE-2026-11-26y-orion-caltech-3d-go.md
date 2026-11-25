# NOTE 2026-11-26y orion -> Caltech: 3-D videos: the USER's go (13:37 PST: "3d video looks good")

The test is approved. Please render all five (bsg, m35mms, m35zams, t9l62, t9mlt) as planned in your 26w, from existing dumps on free CPU cores only (no production node, no new run):
- bsg as the test. Camera for the t9 wedges (90 x 180 deg) and m35zams (3.7 deg wide) set per star so that the wedge fills the frame; show one frame of each star before the full movie and say the camera; if a frame looks wrong, fix it before rendering all.
- m35mms and m35zams: the USER found the 2-D m35 movies too fast, so the 3-D ones run slowly from the start (every dump, 12 dumps per s, about 40 s, cross-faded and labelled "frames interpolated", or a shorter excerpt of the last 2-3 d at 12 dumps per s without interpolation: choose the one that looks better and say which).
- Limits and deep-layer rules as 26w; frames as 1080p H.264, <= 6 MB each (<= 8 MB m35), a poster each.
- Push on docfig-movies-1126 (new files mov3d_<star>.mp4 and mov3d_<star>_poster.png, the five 2-D ones unchanged). Report sizes and the cost (core-hours) and any frame you are unsure of. t9l62 / t9mlt: render what exists now and say the last dump; re-render when the runs end (not before).

-- orion (coordinator)
