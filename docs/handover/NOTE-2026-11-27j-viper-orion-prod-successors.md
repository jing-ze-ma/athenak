# NOTE 2026-11-27j viper -> orion: successor holders queued for the three viper old-scheme productions (USER, 11-25 15:40 PST)

The USER: "queue successors for prod runs". Submitted with the same sbatch files, names and sizes; job ids appended
to HOLDERS in the cfgs (backups *.bak_succ1126; written as temp + mv, the running oprun2.sh re-sources each loop):
- RSG   pRSG            13701779 (4n, holder_rsg.sbatch)  rsg2.cfg   HOLDERS="13514310 13585068 13701779"
- HG    orion_hghold    13701780 (6n, holder6.sbatch)     hg2.cfg    HOLDERS="13360518 13587214 13701780"
- HG262 orion_hg262hold 13701781 (6n, holder6.sbatch)     hg262.cfg  HOLDERS="13255126 13587215 13701781"
All PENDING (Resources). The current holders end 11-25 17:19 (HG262), 18:35 (HG), 18:39 (RSG) PST. No successors
for pr136f 13610529 / pool 13610530 (lent to you idle; say if you want them renewed).

-- viper
