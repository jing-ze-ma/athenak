# NOTE 2026-11-24zzzze6 orion -> Caltech: shorten the queue (USER: progress is slow)

1. **No node is reserved for a t9 MLT wedge:** it runs on viper. 02-09 is free after its items.
2. The **rel-b suite is the release gate**; the rel-a suite may be stopped when the rel-b suite starts, if both would wait for the same node (rel-b contains rel-a's code; report what rel-a ran so far).
3. **Run HR as soon as two nodes are free together** (4 min), even between suite tests. Its GPUs don't wait for a build.
4. Speed F: run m35mms, M35ZAMS and t9l62 quiet base / tip arms first (two repeats only for t9l62 and m35mms), nsys after. Then G, H. The nsys profiles are last.
5. cs arm L2 on the next free node (5 min), then the momentum counter. Any node idle for more than 10 minutes: post it.

-- orion (coordinator)
