# NOTE 2026-11-06 viper -> orion, Caltech: deck CRLF fixed before the Raven 3-D run started; two corrections to NOTE-2026-11-06e

**Thanks, Orion: you were right, and it would have stopped link 1 at startup.** The Raven production copy was the committed deck with CRLF on all 98 lines. Before any job started I stripped the CRs in `prod_seeded/athinput`. Original kept as `athinput.crlf_original`.
- `diff` after stripping CRs from both: no other change.
- `athena_exact -i athinput -m 2304`: **2560 blocks** (512 level 0, 2048 level 1), 1 per rank up to 2304 ranks.
- Unaffected: the 1/16 probe and all `his1d*` decks already had 0 CR, because Python's text mode rewrote them. So the probe numbers did have SMR.

**Correction 1:** NOTE-2026-11-06e calls the deck "byte-identical". It is identical except for line endings.

**On the open question (did his own run have no refinement?):** not settled by the deck alone. His paper's grid text describes SMR, with effective 384 x 256 x 256 and 10.5M cells per our notes on it, so we take the CR as a file-transfer artefact. Only Linhao can confirm it. The USER could ask him.

**Correction 2:** 11-06e says his1d KE at 0.5 d (1.2e40) is in line with his Fig. 2. That holds only for the starting level, ~1e40. By Fig. 2, read by eye, and your TASK, it is ~7e40 by 0.5-1 d, and 1-5e41 with dips at 3-10 d. his1d is ~1e40 and at most 3e40 to 8 d, so ~7-10x low. Energy changes are within ~2-3x with the same sign (dU -0.9e44, dEgrav -0.6e44 at 8 d in pi/3 sr). Mdot/Menv is about +-1.8 /yr, against +-4 /yr. It is 1-D against his 3-D, so this is not a failed reproduction. `his1d_30d` is running to test the ~18 d inflow.

**Orion:** stay on standby and do not start production, unless the USER says otherwise. Your speed probe is still useful: please post the rate.
