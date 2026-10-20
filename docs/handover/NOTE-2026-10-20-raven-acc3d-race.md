# NOTE 2026-10-20 raven -> caltech: fresh accretor (3896d529 + coriolis_lz) is a RACE with Raven run12 (USER)

- Raven run12 (31291344, 2 nodes x 4 A100, est start ~16 h) runs the same fresh setup: fin51-r-fix 3896d529 accretor A100 nofma
  (d7cbe61d), deck v11amwtr + <problem>/coriolis_lz = true (md5 26ed974d).
- Rule: first to START wins; the other is cancelled only AFTER the winner has started (Raven tells you, or you tell Raven in a NOTE the
  moment your new link 1 starts). Your old chain 5088058: cancel when either new copy has started (dir kept).
- RSG E: same rule as before (first start of Raven/Viper/Caltech wins; DeltaAI stood down).
