---
name: use-fork-not-origin
description: "Fetch and push only from the user's fork jing-ze-ma/athenak, never from origin (IAS-Astrophysics/athenak)"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 727372a4-3557-4a71-b294-5104e4fbe37f
  modified: 2026-08-16T13:13:25.363Z
---

For git fetch/pull/push in this repo, use the user's fork `https://github.com/jing-ze-ma/athenak.git` only. Do not fetch from or push to `origin`, which points at the upstream `IAS-Astrophysics/athenak`.

**Why:** The user's own work (e.g. the `general-eos` branch) lives on the fork; upstream is the shared IAS repo where their branches don't exist and where pushes are not wanted.

**How to apply:** Before any fetch/push, target the fork URL (or a fork-named remote) explicitly rather than relying on `origin`. The `fork` remote is the SSH one; use `git push fork <branch>`. If a branch the user names is missing on the fork, ask rather than falling back to upstream.

**Careful:** the `origin` URL in `.git/config` has a GitHub personal access token embedded in it in cleartext, so `git remote -v` prints a live credential. Do not run it unredacted (pipe through `sed -E 's#//[^@]*@#//***@#'`) and never paste its output anywhere. Flagged to the user 2026-08-21; if they have since rotated the token and switched `origin` to SSH, drop this paragraph.
