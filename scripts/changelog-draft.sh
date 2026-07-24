#!/usr/bin/env bash
#
# Draft an [Unreleased] CHANGELOG block from the PRs merged since the last release tag.
#
#   scripts/changelog-draft.sh
#
# It groups those PRs by conventional-commit type (feat -> Added, fix -> Fixed, perf/refactor ->
# Changed) and prints a Keep-a-Changelog block, each entry carrying its PR number. This is a STARTING
# POINT, not the final text: curate the prose, merge related entries, and drop internal ones before
# you release. Prints to stdout — review it, then paste the parts you want under `## [Unreleased]` in
# CHANGELOG.md. (The changelog stays hand-written on purpose; this just removes the blank page.)
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

PY="$(command -v python3 || command -v python)" || { echo "changelog-draft: python not found" >&2; exit 1; }
command -v gh >/dev/null || { echo "changelog-draft: the GitHub CLI 'gh' is required (brew install gh)" >&2; exit 1; }

# the newest vX.Y.Z tag, and the moment it was cut — PRs merged after it are what this release adds
LAST_TAG="$(git tag --list 'v*' --sort=-v:refname | grep -E '^v[0-9]+\.[0-9]+\.[0-9]+$' | head -1)"
[[ -n "$LAST_TAG" ]] || { echo "changelog-draft: no vX.Y.Z tag found — is this the first release?" >&2; exit 1; }
SINCE="$(git log -1 --format=%cI "$LAST_TAG")"

# pass the PR list through the environment, not a pipe: the python program arrives on stdin via the
# heredoc, so stdin is already taken — a `gh | python - <<PY` would silently feed python an empty stdin.
PRS_JSON="$(gh pr list --state merged --base main --limit 300 --json number,title,mergedAt)" \
  "$PY" - "$LAST_TAG" "$SINCE" <<'PY'
import json, os, re, sys

last_tag, since = sys.argv[1], sys.argv[2]
prs = [p for p in json.loads(os.environ["PRS_JSON"]) if p["mergedAt"] > since]

# conventional type -> Keep-a-Changelog section. Types not here (docs, chore, test, ci, style, build)
# are not user-facing, so they are listed separately for you to pull in only what matters.
SECTION = {"feat": "Added", "fix": "Fixed", "perf": "Changed", "refactor": "Changed"}
buckets = {"Added": [], "Changed": [], "Fixed": []}
skipped = []
head = re.compile(r"^(\w+)(?:\([^)]*\))?!?:\s*(.*)$")   # type(scope)!: description

for p in sorted(prs, key=lambda p: p["number"]):
    m = head.match(p["title"])
    typ = m.group(1).lower() if m else ""
    desc = (m.group(2) if m else p["title"]).strip().rstrip(".")
    entry = f"- {desc[:1].upper()}{desc[1:]} (#{p['number']})" if desc else f"- {p['title']} (#{p['number']})"
    (buckets[SECTION[typ]].append(entry) if typ in SECTION
     else skipped.append(f"#  - {p['title']} (#{p['number']})"))

if not prs:
    print(f"# no PRs merged since {last_tag} — nothing to draft")
    raise SystemExit

print(f"# DRAFT — {len(prs)} PRs merged since {last_tag}. Curate before releasing; then paste under ## [Unreleased].\n")
for sec in ("Added", "Changed", "Fixed"):
    if buckets[sec]:
        print(f"### {sec}")
        print("\n".join(buckets[sec]))
        print()
if skipped:
    print(f"# {len(skipped)} non-user-facing PRs skipped (docs/chore/test/…) — uncomment any that matter:")
    print("\n".join(skipped))
PY
