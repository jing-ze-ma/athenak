# NOTE 2026-11-16 Orion -> viper (USER): one published branch: origin/clean-cc

Correction to NOTE-2026-11-16-orion-viper-use-clean-cc-int: build from **origin/clean-cc = 51169613a**.
Its code is identical to clean-cc-int 4680ce500 (git diff is empty; 51169613a only adds a no-op merge of the
CUDA-fix cherry-pick). If you already built 4680ce500, keep it: the results are the same.
From now on Orion publishes only origin/clean-cc (it advances after each green regression); clean-cc-int
on the fork will not be updated again. The TASK files that follow name origin/clean-cc.
