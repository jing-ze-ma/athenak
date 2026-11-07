3-D movies of the GLOBAL cubed-sphere runs (cs3d.py, 10-21)
=============================================================
Runs: rsgE_glob = RSG E (rsg_E_1020/prod63, rsg.*), hggl_glob = He giant global (hgglob_1021/prod, hegiant.*). READ ONLY on the runs.

Refresh (incremental, returns at once; one expansion CPU job, 8 cores / 64 GB; skipped if one is already queued):
    bash /resnick/groups/carnegie_poc/jingze/vis3d/rsgE_glob/refresh.sh
    bash /resnick/groups/carnegie_poc/jingze/vis3d/hggl_glob/refresh.sh
  VIEWS=single (or layer) bash .../refresh.sh   renders only that view (default both).
  The job (render_cs.sh -> venv_vtk/pvpy.sh cs3d.py RUN --view V) reduces only NEW dumps (cache RUN/red/r_NNNNN.npz,
  ~35-50 s per dump: tau_R = 2/3 crossing, T_eff, v_r on cut faces and 2 inner shells), renders missing frames
  (~1.5-3 s each; frames are re-rendered when the newest dump changes the colour ranges / view radius), encodes
  H.264 yuv420p +faststart (CRF raised until < 8 MB) and copies mp4 + poster into the page publish/movies dir.

Outputs:  RUN/out/{rsgE_3d,hggl_3d}{,_layers}.mp4, *_poster.png, *.json (clims, frames, fps)
          copied to rsgE_page/out/publish/movies/ and hgglob_page/out/publish/movies/
Views:    single = 90 x 90 deg wedge (octant lon 0..90, lat 0..90) cut out of the corrugated tau_R = 2/3 photosphere;
          layer  = photosphere cut lon +-90 x lat 0..90, v_r sphere at shells[0] cut lon +-67.5, v_r sphere at shells[1]
                   cut lon +-45 (RSG 0.7 / 0.4 R_ph,IC; He giant 0.75 / 0.35 R_ph,IC).
Style:    black bg, photosphere T_eff = (F_r(r_ph)/sigma)^1/4 (RSG 'hot' from 0.22, He giant 'Blues'), cut faces / shells v_r
          coolwarm, clims = 1-99 % of the newest dump (v_r symmetric), 1000x760, bars 24/22 pt, two-line title.
Test frames: pvpy.sh cs3d.py RUN --view V --only N --tag _x  -> RUN/test/ (no encode). Run in sbatch (login 8 GB cap).
