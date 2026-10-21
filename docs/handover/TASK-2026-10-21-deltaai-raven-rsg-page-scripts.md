# TASK 2026-10-21 deltaai -> raven: please push your RSG cubed-sphere page pipeline (USER: plots + video "like the golden url" for the DeltaAI RSG E copy)

1. REQUEST: push the scripts that build your "RSG g25 Production" page (claude.ai/artifact/TdGw9wAzb191iuy5x229TD; footer
   says upd_page.sh) = the cubed-sphere versions of reduce.py / mk_figs.py / render (movie) / build_page.py + any config, e.g. to
   docs/handover/rsg-page/. hg145_page (wedge) is all I have; it does not handle 6 panels. If you also have the panel->Cartesian
   assembly for slices / 3-D (block gid -> panel = gid // 4 here, 2x2 blocks of 64x64 per panel), include it. I will run it on
   the DeltaAI RSG E clip copy (bin dir /work/hdd/bivj/jma20/rsg_E_1021/run/bin, basename rsg) and publish a page in the
   Golden RSG Runs style (time series, radial profiles, slices, movies).
2. STATUS link 1 (3425380): clean, ~3.0 s/cycle incl. I/O (2.25 s without rst writes), dt 1e3 s, ~14 sim d per wall hour.
   Latest: elapsed=1.477469e+03 cycle=450 time=4.500000e+05 dt=1.000000e+03
3. DISK (FYI, USER decides): rst = 30.8 GB each every 0.5 d (~every 2 min wall); prune keeps 4 -> ~150 GB steady. Dumps
   hydro_w 1.08 GB + m1 0.86 GB per day + m1_vet 0.86 GB per 5 d -> ~29 GB per wall hour. /work/hdd/bivj 248 of 1000 GB used ->
   quota in ~1 day of continuous running (a 48 h batch link would write ~1.4 TB). Deck unchanged per your NOTE; options are
   moving dumps to /projects/bivj (1 TB free) and/or thinning old dumps. Tell me if you have a preference.
