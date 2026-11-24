M35FAST Caltech m35mms arms from rst m35mms.00010.rst (cycle 21945, t 1.728e5), 1 node 4 H200, deck = restart (no -i):
  env ATHENA_DT_REPORT=25 (+ the production env), then:  $BIN -r <rst dir>/m35mms.00010.rst -d <run dir> $(cat Dx_args.txt)
  Dx_args.txt names every <hydro>, <rad_m1>, <time> key (and mesh layout keys) of the IA12.1 tile deck (rst-stored-keys rule).
  D0  : BIN = IA12.1 fa59bffd6, key off, 2000 cycles (nlim 23945)                              -- reference + lever B late cost on H200
  D1  : BIN = caltech-m35fast tip, hydro/rad_signal_speed_diff=true, 2000 cycles                -- lever A
  D1h : same binary, key on, cfl 0.15, 1000 cycles (nlim 22945)                                  -- dt-halving yardstick
  D2  : same binary, key on + clean_gas_secant=max_stage2 + zalesak_converge=true, 2000 cycles    -- levers A+B
  If the restart FATALs on mesh/layout_invariant_mode=tile (IA11 rst written without tile), drop the 3 mesh/ keys and say so.
Return: hst (user+hydro), run.log grep "DTREP|DTREG|CLEAN WARNING|FATAL|^CLEAN cycle|elapsed=", hydro_w bins (output2), s/cycle.
