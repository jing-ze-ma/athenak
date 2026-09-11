---
name: confounded-tests-rejected
description: "The user rejects experiments whose result would be confounded, and calls a halt when hypotheses are churning faster than they are being settled"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ee41f6f0-a3cb-4361-83ad-7778d2fbc04d
  modified: 2026-08-23T09:26:35.996Z
---

On 2026-08-23 I proposed testing a Gamma_1 / reconstruction hypothesis by rerunning with
donor-cell reconstruction. The user replied only "I don't think that helps" -- and was
right: donor cell changes the whole scheme's diffusivity, and dissipation ALONE was
already known to delay that failure (the failure times 0.370 < 0.485 < 0.584 < 1.072
tracked increasing dissipation). A surviving donor-cell run would have proved nothing.
Shortly after, having watched several of my hypotheses die the same day I proposed them,
they said "stop and write up what's established".

**Why:** in this project a single test costs hours of GPU allocation, and the failure mode
being chased is marginal, so almost any change that adds numerical dissipation delays it.
That makes confounded tests actively expensive: they consume the allocation AND produce a
result that cannot discriminate. The user tracks this and will veto the design rather than
the topic.

**How to apply:** before proposing an experiment, state what each possible outcome would
prove, and check that the manipulated variable is the ONLY thing that changed. Prefer
tests that hold diffusivity, resolution and the initial state fixed -- e.g. "evaluate
Gamma_1 at the interface instead of reconstructing the stored value" (same order, same
dissipation) over "switch to a more diffusive reconstruction". Flag known confounds
yourself, in the proposal, rather than after the result comes back. And when a run of
hypotheses has been refuted without converging, offer to consolidate instead of proposing
the next one -- see [[dhj-ck-eos-blowup]], where nine died before the user called it.
