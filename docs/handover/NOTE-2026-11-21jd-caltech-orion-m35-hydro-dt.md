# NOTE 2026-11-21jd Caltech -> orion: m35mms on B_p1 passed the old FATAL cycle; hydro dt dips from floor cells

- m35mms (cfl 0.9): c45505 "floor_pay = rad_column left 1.446214e+33 of the fluid floors energy unpaid (created)", then hydro dt COLLAPSE 21.7 -> 1.16e-3 at c45506 (rank 2); c45702 dt 6.0 -> 1.0; CW at c45729 (t2s 2 stage fallback 1, recovered), hydro dt 20.3 -> 0.54 at c45730. Now past the old FATAL cycle 45706: c45750, dt 2.7 and recovering, FATAL 0, dt_retry 0. CW 2 total, both > 100 cycles after the restart, so no cfl fallback applied.
- i1_m35mms_f1k (cfl 0.9): CW 1; hydro dt dips 24 -> 1.5 (c43435) and 24 -> 4.3-5.9 (c43716-43777), recovered to dt 22.9 at c43800.
- t9l62 (cfl 0.6 since 16:41): CW 0, dt 100. t9r (cfl 0.6 + 1000 K): CW 0, dt 36.
The m35mms dips are hydro-limited (like t9r c26512: one floor cell), the 1e-3 one right after a floor_pay energy creation. Holders end 17:46.

-- Caltech
