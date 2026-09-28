# NOTE for the Caltech and DeltaAI sessions: remove a revoked GitHub token (2026-09-28)

On viper, the `origin` remote of the AthenaK checkout had a GitHub classic personal access token embedded in its
URL (`https://jing-ze-ma:ghp_...@github.com/IAS-Astrophysics/athenak.git`). The user has **deleted that token on
GitHub** (09-28), and viper's remote is clean now. Do the same cleanup on your machine, once, at the start of
your next session. Report what you found to the user.

## What to do

1. **Look for any GitHub token on this machine** (prints file names only, never the token):
   ```bash
   for d in $(find ~ -maxdepth 4 -name .git -type d 2>/dev/null); do
     git -C "$d/.." remote -v 2>/dev/null | grep -q 'ghp_\|github_pat_' && echo "TOKEN IN REMOTE: $d"
   done
   grep -l 'ghp_\|github_pat_' ~/.git-credentials ~/.gitconfig ~/.config/gh/hosts.yml ~/.netrc \
       ~/.bashrc ~/.bash_profile ~/.profile 2>/dev/null
   grep -c 'ghp_\|github_pat_' ~/.bash_history 2>/dev/null
   env | grep -l 'GH_TOKEN\|GITHUB_TOKEN' >/dev/null 2>&1 && echo "token in environment"
   ```
   Also check the git config of every AthenaK clone and worktree: `git config --show-origin -l | grep -i url`.
2. **Remove what you find:**
   - **A remote URL:** `git remote set-url origin https://github.com/IAS-Astrophysics/athenak.git` (upstream,
     read-only, no token needed). Pushes go to the fork over SSH:
     `git@github.com:jing-ze-ma/athenak.git`.
   - **A line in `~/.git-credentials`, `~/.netrc` or a shell rc file:** delete that line. Do not print it.
   - **Shell history:** tell the user how many lines match, and let the user decide whether to scrub them.
3. **Do not create a new token.** Use an SSH key for pushes: Caltech already uses `~/.ssh/id_ed25519`. On
   DeltaAI, ask the user to add a key if pushing is needed.
   - Never put a token in a remote URL, a file, a commit, a job script or an environment variable in a
     committed script.
4. **Report to the user:**
   - which clones or files had a token and what you changed;
   - that `git remote -v` in every AthenaK clone now shows no credentials.

If a fetch or push fails with an authentication error mentioning the old URL, that is this revoked token.
Switch to the URLs above.
