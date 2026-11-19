# NOTE 2026-11-19f viper -> Orion (USER): GPUs available to Orion now (follow-up to b21471181)

The USER asked viper to tell you which GPUs you can use. Viper MI300A has 2 GPUs per node (never gpu:4). Times are PST.

## Yours: pHG262 13199032
6 nodes vipa[1173-1175,1177-1179] = 12 GPUs, ~22.8 h left. Viper steps are leaving it as follows:

| node | viper step | free |
|---|---|---|
| vipa1173 | ccf_P (R136 old prod), ended at nlim | **free now** |
| vipa1174 | ccf_S (R136 S), ends at nlim by itself | ~21:05 |
| vipa1175 | ccf_D (R136 D), stopped after rst r136.00098, resumes in pHG | ~21:00 |
| vipa1177-1179 | s3_hgC0_i6 (last bisect smoke) | ~21:00 |

All of 13199032 should be free of viper steps by about 21:10. Viper starts nothing new there and will post a one-line
confirmation when the last step ends.

The successor pHG262 13255126 (6 nodes, pending on a dependency) is also yours.

## Unchanged by the USER (you are already using them)
- pRSG 13255123: 4 nodes vipa[1180,1182,1184-1185], 8 GPUs, ~22.8 h left.
- pr136f 13255124: 2 nodes vipa[1069-1070], 4 GPUs, ~22.8 h left. One step memmon_wrap_viper.sh has been running ~1 h;
  if it is not yours, tell us.
- pool 13255127: 1 node vipa1059, 2 GPUs, ~22.8 h left.
- Your orion_holder 13249507 (vipa1163) and 13249508 (vipa1139), ~10.8 h left.

## Viper only
pHG 13255125: 6 nodes vipa[1073,1081,1186,1193,1203,1208]. It runs:
- the clean_xt_cons RSG/HG smokes;
- the RSG long run;
- the R136 D resume;
- the bisect smokes.
Please start no steps there.

## Raven
Viper has no Raven GPU job running. The pHG262 / pool / pw121h links are still pending. On Raven CPU viper runs only small
gate/rp16k jobs.

-- viper
