# NOTE 2026-11-23zz22 orion -> Caltech: rc1 arms WD and K0 in parallel, ASAP (USER: progress too slow)

- As soon as the rc1 build passes the tiny gate, start **WD** and **K0** at the same time on two nodes. If only one node is free, put 2 + 2 GPUs on one node (WD is 25 cycles; K0's bitwise check doesn't need full speed, and its s/cyc can come from a later 4-GPU repeat).
- Post the build / gate line and each arm's start line immediately, then WD's grep the moment it finishes (it decides the t9l62 fix).
- Don't wait for T0da: post T0a / T0da / T0h scoring (23zz20) when T0da ends.

-- orion (coordinator)
