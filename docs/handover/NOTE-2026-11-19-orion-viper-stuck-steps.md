# NOTE 2026-11-19 orion -> viper (USER): please clear stuck Orion steps in pHG262 13199032

The USER asks viper to clear these steps, or to escalate to MPCDF support if you cannot. They are all Orion's own test steps
inside pHG262 13199032. Each was cancelled with scancel but has not exited, ~45-105 min later:

| step | name | nodes | what happened |
|---|---|---|---|
| 13199032.40 | orion_vg_S3B_rsg684506_ck__rs | vipa[1177-1179] | hung in its final restart dump (collective MPI-IO, ~124 GB to /viper/ptmp2); even `stat` on its rst/ dir hung |
| 13199032.55 | orion_vg_S3B_rsg684506_ck__r1 | vipa[1177-1179] | follow-on arm; cancelled |
| 13199032.54 | orion_vg_S3B_rsg684506_head__rs | vipa[1173-1175] | hung in its final restart dump, same pattern |
| 13199032.59 | orion_vg_build | vipa1173 | build step; cancelled |

Effects:
- vipa1177 cannot launch tasks at all ("task 0 launch failed: Unspecified error").
- New steps that ask for real memory (e.g. --mem=16000M builds) fail to launch on vipa1173.
- Only vipa1178/1179 work reliably.

Likely cause: a /viper/ptmp2 stall during two large collective MPI-IO restart writes (6 ranks on 3 nodes each). Our new
independent/streamed writer is meant to avoid this.

Options, as you see fit:
- (a) clear the steps (root/admin) or escalate to MPCDF support;
- (b) tell us if ending pHG262 13199032 so that its successor 13255126 starts on fresh nodes is acceptable to you. We will not
  touch the holder ourselves.

Orion's queue now skips vipa1177 (vgpu/skipnodes).

-- orion
