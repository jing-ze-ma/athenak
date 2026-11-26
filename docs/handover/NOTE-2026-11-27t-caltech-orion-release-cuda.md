# NOTE 2026-11-27t caltech -> orion: CUDA suite of the release code caltech-rel-82231a337 (your 27t)

## Builds (worktree wt_orion1121/rel82; CUDA 12.9 + g++13; 02-03 (suite) and 02-17 (star_vet), 28 cores each)
- BUILD rel82 rc=0 wall=299 s TUs=280
- BUILD_OK rel82 md5=e80fcc9c65def95726b5dfd80f9b861a prov=82231a337f4e5f5a5a4c314dfaaeda8e88b2ae14
- GATE FAIL rel82: athena_rel82_82231a337 md5 e80fcc9c node hpc-sm-02-17 pass1= log date 11-25 17:05:55 end: 
- BUILD suite rc=0 wall=291 s
- BUILD_OK suite md5=ee84ed64f900c6eb485736a344307211 prov=82231a337f4e5f5a5a4c314dfaaeda8e88b2ae14
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).

## 3. Star restart: m35zams.00010.rst (hard link of the finished run's restart, c223416), +2 cycles, 4 ranks on 02-16, own keys; release 82231a337 against b5ade8753 (star_vet binaries)
- rel82_82231a337 rc=0 Terminating on cycle limit CW 0 FATAL 0
- next5_b5ade8753 rc=0 Terminating on cycle limit CW 0 FATAL 0
- m35zams.hydro_w.00433.bin: payload md5 5a184d18333e / 5a184d18333e -> IDENTICAL
- m35zams.vet.00433.bin: payload md5 7f264e050072 / 7f264e050072 -> IDENTICAL
- m35zams.vet_face.00019.bin: payload md5 88abb8c857c9 / 88abb8c857c9 -> IDENTICAL
- m35zams.vet_vet.00019.bin: payload md5 f7f480942c33 / f7f480942c33 -> IDENTICAL
- user.hst: 2 / 2 rows, printed-bitwise
- hydro.hst: 2 / 2 rows, printed-bitwise

## Gate (27u): the old tiny gate deck FATALs on the removed key, as you expected
- GATE FAIL rel82 (16b039e9a) = "<rad_vet/clean_face_hctau> = false: the key was removed, the code runs true" (rad_vet_args.hpp FATAL), before pass 1. Gate deck without that line (b1_1124/gate_relt.athinput; no other removed key FATALed), 1 rank, same arguments, on both binaries:
- GATE2 rel82_82231a337 md5 e80fcc9c: pass1 3.720e-11 end: Terminating on cycle limit FATAL 0 CW 0; hst md5 73cc8651512b
- GATE2 next5_b5ade8753 md5 706572b0: pass1 3.720e-11 end: Terminating on cycle limit FATAL 0 CW 0; hst md5 73cc8651512b
- GATE2 DONE
- So the release and b5ade8753 agree on the gate to the printed digits and in the whole user.hst (md5 73cc8651512b each); 3.720e-11 replaces the hctau = false reference 3.686e-11.
