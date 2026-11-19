# NOTE 2026-11-19n viper -> Orion (USER): Raven pool holders are yours (supersedes the pool line of 9c692380a)

USER (~00:55 PST): "don't run dhj on raven pools. also give raven pools to orion claude".

Yours from now:
- Raven pool 32490420 (ravg1106, 4 A100, ~22.5 h left): idle, only keepalive tasks.
- Raven pool 32490421 (1 node, starts ~01:07 PST).

Use them as you like. Their batch scripts are viper's pool_1107 holder.sbatch, which runs tasks from
/raven/ptmp/jinma/pool_1107/queue. Viper has emptied that queue of its tasks:
- the dhj handback copies are withdrawn;
- the R136 HO task is parked and its requeue gate is off (QUEUE_OK_r136.raven.off_1119).

You can put your own tasks there, or attach steps with srun --jobid --overlap. Please don't cancel the holders.

Still viper's on Raven:
- pool 32435128 (ravg1099, ~5.5 h left): W189 va1e7 dhj step prod.w189va1e7, running since 11-18.
- pHG262 32307710 / 32307711: HG production (started 00:41 PST).

Unchanged: orion_integ 32574067 / 32574068 are yours.

-- viper
