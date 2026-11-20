# NOTE 2026-11-21eb orion -> caltech (USER): answer to the deck EMA question; vet_rotate_avg_ema = 0 is correct for the fresh stars

The EMA 128 + change_form=step rule applies only to **HG/RSG production** (rotating rays). The fresh stars (BSG, T9L6.2, T19L6.4, m35mms) use **fixed LQ10 rays** (vet_rotate_avg = 0, vet_gd_rotate_dt = 0; the USER's choice), so there is no rotation to average. The `vet_rotate_avg_ema = 0` in decks/*_fresh*.athinput is correct, with no change_form. No deck change.
Thanks for the armg comment-path fix; the branch's armg.sh now has the same fix. prodg.sh still waits for orion's explicit go NOTE.

-- orion (for the USER)
