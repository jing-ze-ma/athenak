# NOTE 2026-11-17 Orion -> viper (USER): stag2 floor / work fixes pushed (answer to TASK-2026-11-17-viper-orion-push-stag2-floor-fixes.md)

## Branches on the fork (pushed now)
- **clean-cc-stag2-gpu4 = 056fda78f**: Orion's local gpu3 head (g58), the newest gated staggered head. It contains:
  - fcfc0929e (BE step's E pays fref_wacc): the exact commit.
  - ac5bc4dca = dccd3e031 (stage step's E pays the reference work of its hydro stages): the same patch rebased (git patch-id --stable identical, 291f8411b41f).
  - 1d247bebe (write-back eint reject with a round-off margin, your item 1).
  - 4ba11d71b = 32b045052 (time2_auto_be = 0), then **3fb92f134 REVERTS it** (see the follow-ups below).
  - 32bd5aa46 = ce55710e3 (BE closure: the predicted closure state relaxes the local exchange also without K1; patch-id identical).
  - kokkos gitlink 160000.
- Exact hashes for porting:
  - **clean-cc-stag2-ref-dccd3e031 = dccd3e031** (tip of the rsthand3 line).
  - **clean-cc-stag2-ref-32b045052 = 32b045052** (tip of the abe0 line).
  - fcfc0929e is reachable from gpu4.

## Follow-ups to those commits
1. **32b045052 is reverted on staggered (3fb92f134).**
   - With the BE closure fix ce55710e3, the default auto-BE path holds dt: R136 r136.00009, 150 cycles, 4 A100, mean dt 5.412 with auto_be 16 or with 0 (old production 5.40).
   - Production decks name auto_be 16, so staggered keeps production's default.
   - The earlier dt loss under auto_be came from the BE closure skipping the exchange relax: un-relaxed hot gas in the rays gave a BE face-F spike.
   - Suggestion for your item 3: port ce55710e3 first, then check whether the cell layout still needs auto_be 0. If not, keep 16.
2. fcfc0929e / dccd3e031: no later fixes on staggered. One related fix is local and still gating: 21852f122 (FCT reads the owner's Phi at a MeshBlock's high face). It matters only with time2_face_fct on: there the limiter at the theta wrap and seams created energy, R136 booked ledger +7.5e-8 of <L_top>T.
3. Item 2 (clean_stage_floor default false): consistent with staggered, which has no stage floor. Staggered instead rejects a write-back that would cool gas below e(tfloor) (4e6456619).

## Status
- gpu4 is the last fully gated staggered head (trimmed CPU gate + one A100 round).
- Newer local candidates (g61/g62) are not landed. A star-row regression is being bisected.
- Ready for the fresh CPU regression against push 18 (incl. sfl/, sflc/ and the t / dt columns) when clean-cc-physfix2-viper is up.

-- Orion
