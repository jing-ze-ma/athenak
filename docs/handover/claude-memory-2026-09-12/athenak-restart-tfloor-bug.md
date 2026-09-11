---
name: athenak-restart-tfloor-bug
description: "FIXED ba2f0943: any general-EOS run using tfloor_kelvin could be started but never restarted, because GetOrAdd records the tfloor default into the restart dump and the both-set check then fires. Added ParameterInput::IsParameterDefaulted()."
metadata:
  node_type: memory
  type: project
---

`<block>/tfloor` and `<block>/tfloor_kelvin` are the same floor in different units and
setting both is refused (`eos.cpp`). The check used `DoesParameterExist("tfloor")`, which
is right on a fresh start and WRONG on a restart:

* `EquationOfState` itself calls `GetOrAddReal(bk,"tfloor",FLT_MIN)`, which ADDS the
  parameter with its default;
* `restart.cpp` writes the whole parameter list into the restart file;
* the next leg reads a `tfloor` the previous leg invented, sees both, and aborts.

**Consequence: every general-EOS run using `tfloor_kelvin` could be started but never
continued** -- including the shipped `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`, hence
any production run needing a checkpoint. Found trying to continue a 3 h GPU run.

**Do NOT fix it by comparing the value to FLT_MIN.** `ParameterDump` keeps six significant
digits, so 1.1754943508222875e-38 comes back as 1.17549e-38 and the comparison fails. I
tried that first and it did not work.

The dump already annotates defaulted parameters (`# Default value added at run time`) and
`InputLine` stores the comment, so **`ParameterInput::IsParameterDefaulted(block, name)`**
answers it exactly. That helper is now available for any other check with the same shape --
and this shape will recur: any "the user set both X and Y" validation on a parameter that
also has a GetOrAdd default is broken across restarts in the same way.

Verified: both-set still aborts; tfloor_kelvin alone runs, checkpoints and restarts with
the same converted floor; ideal-gas inputs setting `tfloor` unaffected; the correlated-k
regression test passes. See [[dhj-ck-eos-blowup]] for what it was blocking.
