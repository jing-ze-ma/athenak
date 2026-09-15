---
name: athenak-fork-git-setup
description: "Git remotes and auth for AthenaK — origin is the user's fork over SSH, upstream is IAS-Astrophysics over https"
metadata:
  node_type: memory
  type: project
  originSessionId: d0224e66-3587-4359-a545-9608682177c3
  modified: 2026-08-12T04:34:36.095Z
---

The AthenaK checkout at `/orion/u/jinma/ATHENAK/athenak`:

- `origin` = `git@github.com:jing-ze-ma/athenak.git` (the user's fork) — **all git operations go
  here**
- `upstream` = `https://github.com/IAS-Astrophysics/athenak.git` (the shared repo, read only in
  practice)

The user's real work lives on fork branches, not on local uncommitted state — when asked to "get
the files from GitHub", check the fork's branches first, not `main`.

## Auth (settled 2026-08-12)

Pushing over https does not work on this machine: no `gh`, no credential helper, no `~/.netrc`,
and the user's PAT was rejected 401. **`origin` now uses SSH** with a dedicated key
`~/.ssh/github_ed25519` (a `github.com` block in `~/.ssh/config` pins `IdentitiesOnly`, so the
unrelated FAU-cluster `id_ed25519` is not offered). The public half is added on GitHub as a
**deploy key on the fork with write access**, not an account key — `ssh -T git@github.com`
answers "Hi jing-ze-ma/athenak!" rather than naming the user. It pushes to that ONE repo; any
other repo would need its own key or an account-level key.

Git identity is set repo-locally to `Jingze Ma <jingzema0211@gmail.com>` (it was previously
unset, so commits were authored `jinma@orion02.bc.mpcdf.mpg.de`).

**The `upstream` URL used to carry a GitHub PAT in plaintext** in `.git/config` (`ghp_IT5C...`,
a classic token, so possibly `repo` scope on the SHARED upstream repo). Removed from the config
2026-08-12. Filesystem exposure was contained — `/orion/u/jinma` is `drwx------`, so no other
cluster user could reach it — but the value was printed into a session transcript. If a token is
ever needed again it goes in a credential helper, never in a remote URL.

## Open items (as of 2026-08-12)

- **The old PAT is still NOT REVOKED.** The user said "I'll do it later" — github.com/settings
  /tokens, "Tokens (classic)" tab, delete. Nothing depends on it: `origin` is SSH and
  `IAS-Astrophysics/athenak` is public so `fetch upstream` works anonymously. Worth a nudge, not
  a blocker.
- Local branch `general-eos-backup-preauthor` still exists — the safety net from the author
  rewrite. The push is verified, so it can be deleted whenever the user says so.

**The user's terminal can paste IN but cannot copy OUT.** Anything they need to transfer to a
browser must go to a file (`/u/jinma/github_key.txt` was used for the public key) or to a
published artifact with a copy button — do not ask them to copy from terminal output.

See [[run-directory-untouchable]] and [[general-eos-project]].
