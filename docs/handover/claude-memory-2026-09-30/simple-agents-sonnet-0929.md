---
name: simple-agents-sonnet-0929
description: user 09-29: use Sonnet (the "sonnet" model = Sonnet 5; there is no Sonnet 5.5) for SIMPLE agents from now on; Opus 5.5 for everything else
metadata:
  type: feedback
---
Simple, mechanical, verifiable agent tasks -> Agent(subagent_type "worker", model "sonnet"): run a given build
(build_inc_viper.sh <target> <commit>) and report md5; submit a prepared job script and report ids; collect numbers
from finished runs with an existing script into a table; file housekeeping (copy/list/delete named files); status
checks. Everything with judgement (physics gates, noise tests, code changes, merges, diagnosis, analysis design)
stays on Opus 5.5 (worker default). Verify every Sonnet result as usual.
**Why:** user 09-29 "try to use sonnet 5.5 as simple agents from now on" (earlier the same night the user said no;
this supersedes that and the Caltech "no haiku/sonnet" line for viper). Sonnet 5.5 does not exist here: a model "sonnet" agent reported itself as "Claude Sonnet 5, model ID claude-sonnet-5" (checked 09-29). Re-check the version by asking a sonnet agent if the user mentions a newer one.
**How to apply:** pass model: "sonnet" on the Agent call for simple tasks; say in the brief exactly what to run and what
to report; no decisions left to the agent. Related: [[delegate-to-sonnet-when-adequate]], [[model-setup-opus55-medium]].
