# NOTE Caltech -> Raven (10-15 06:10 PDT): agcE1r and hgF4w move to fin51-r-fix b16b66db at their next link (USER)

- After the re-gate (7b54337a: He giant all bitwise; E1 C = L2 = F = F_all bitwise, P != B only from the first natural retry),
  USER chose to switch both production runs to b16b66db (md5 2df1a7744c95857e42e4904ccf32bc18). run.cfg of
  rol2u_1015/e1/chain and hegiant_1014/wide4/run now point to it (previous cfg saved as run.cfg.5ab8c87e); the running links keep
  5ab8c87e (run.cfg is read at link start only): agcE1r 4553459 ends ~09:03 -> 4553460 uses b16b66db; hgF4w 4513229 (started 05:56,
  12 h slot, ROL2U from day 30, ~0.3 s/cycle) -> 4513230 / twins use b16b66db. Same ROL2U keys.
