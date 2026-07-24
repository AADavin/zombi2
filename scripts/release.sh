#!/usr/bin/env bash
#
# One-command release for ZOMBI2.
#
#   scripts/release.sh patch     # bump the last number   (0.4.0 -> 0.4.1) — a fix release
#   scripts/release.sh minor     # bump the middle number (0.4.0 -> 0.5.0) — a feature release
#   scripts/release.sh major     # bump the first number  (0.4.0 -> 1.0.0)
#   scripts/release.sh 0.4.0     # or an explicit version, if you must
#
# The version is COMPUTED from the current one (no typing a version, no typos). What it does, on an
# up-to-date `main` with a clean working tree:
#   1. sets the single-source version (`__version__` in zombi2/__init__.py — hatchling reads it),
#   2. rolls the CHANGELOG's [Unreleased] entries into a dated `[X.Y.Z]` section,
#   3. commits `release: bump version to X.Y.Z`, tags `vX.Y.Z`, and pushes both,
#   4. publishes the GitHub Release — which is what triggers the PyPI upload
#      (.github/workflows/release.yml, gated on `release: published`).
#
# It stops and asks before step 3/4 (the irreversible, outward part). Pass --yes to skip the prompt.
set -euo pipefail

YES=0
ARGS=()
for a in "$@"; do [[ "$a" == "--yes" ]] && YES=1 || ARGS+=("$a"); done
BUMP="${ARGS[0]:-}"

die() { echo "release: $*" >&2; exit 1; }
[[ -n "$BUMP" ]] || die "usage: scripts/release.sh (patch|minor|major | X.Y.Z) [--yes]"
cd "$(git rev-parse --show-toplevel)"

PY="$(command -v python3 || command -v python)" || die "python not found"
command -v gh >/dev/null || die "the GitHub CLI 'gh' is required (brew install gh)"

# --- read the current version and compute the next -----------------------------------------------
# The version is derived from zombi2/__init__.py, not typed, so a keyword can never fat-finger it.
CURRENT="$("$PY" - <<'PY'
import re, pathlib
m = re.search(r'^__version__ = "(.*)"$', pathlib.Path("zombi2/__init__.py").read_text(), re.M)
print(m.group(1) if m else "")
PY
)"
[[ "$CURRENT" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
    || die "could not read a valid current version from zombi2/__init__.py (got '$CURRENT')"

case "$BUMP" in
    patch|minor|major)
        VERSION="$("$PY" - "$CURRENT" "$BUMP" <<'PY'
import sys
x, y, z = (int(n) for n in sys.argv[1].split("."))
kind = sys.argv[2]
if kind == "patch":   x, y, z = x, y, z + 1
elif kind == "minor": x, y, z = x, y + 1, 0
else:                 x, y, z = x + 1, 0, 0
print(f"{x}.{y}.{z}")
PY
)" ;;
    *) [[ "$BUMP" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
           || die "not a bump keyword (patch|minor|major) or an X.Y.Z version: '$BUMP'"
       VERSION="$BUMP" ;;
esac
TAG="v$VERSION"
DATE="$("$PY" -c 'import datetime; print(datetime.date.today().isoformat())')"

# --- preconditions: fail before touching anything ------------------------------------------------
[[ "$(git branch --show-current)" == "main" ]] || die "not on main (checkout main first)"
git diff --quiet && git diff --cached --quiet || die "working tree is not clean"
git fetch --quiet origin
[[ "$(git rev-parse HEAD)" == "$(git rev-parse origin/main)" ]] \
    || die "main is not in sync with origin/main — 'git pull' first"
git rev-parse "$TAG" >/dev/null 2>&1 && die "tag $TAG already exists"
[[ -f CHANGELOG.md ]] || die "CHANGELOG.md not found"
grep -q '^## \[Unreleased\]' CHANGELOG.md || die "CHANGELOG.md has no '## [Unreleased]' section"

echo "ZOMBI2 release  $CURRENT -> $VERSION   (tag $TAG, $DATE)"
echo "  · set __version__ in zombi2/__init__.py"
echo "  · roll CHANGELOG [Unreleased] -> [$VERSION] - $DATE"
echo "  · commit + tag $TAG + push origin main $TAG"
echo "  · gh release create $TAG   <-- PUBLISHES to PyPI"
if [[ "$YES" != 1 ]]; then
    read -r -p "Proceed? [y/N] " ok
    [[ "$ok" == y || "$ok" == Y ]] || die "aborted"
fi

# --- 1. bump the single-source version -----------------------------------------------------------
"$PY" - "$VERSION" <<'PY'
import re, sys, pathlib
v = sys.argv[1]
p = pathlib.Path("zombi2/__init__.py")
s = p.read_text()
new = re.sub(r'^__version__ = ".*"$', f'__version__ = "{v}"', s, count=1, flags=re.M)
if new == s:
    raise SystemExit("release: could not find the __version__ line in zombi2/__init__.py")
p.write_text(new)
PY

# --- 2. roll the CHANGELOG, and extract this version's notes --------------------------------------
NOTES="$("$PY" - "$VERSION" "$DATE" <<'PY'
import sys, pathlib
v, d = sys.argv[1], sys.argv[2]
p = pathlib.Path("CHANGELOG.md")
lines = p.read_text().splitlines()
# insert the dated header just after [Unreleased]
out, notes, grab = [], [], False
for ln in lines:
    out.append(ln)
    if ln.strip() == "## [Unreleased]":
        out.append("")
        out.append(f"## [{v}] - {d}")
        grab = True
        continue
    if grab:
        if ln.startswith("## ["):     # next section: stop collecting notes
            grab = False
        else:
            notes.append(ln)
p.write_text("\n".join(out) + "\n")
print("\n".join(notes).strip())
PY
)"
[[ -n "$NOTES" ]] || die "the [Unreleased] section is empty — nothing to release"

# --- 3. commit, tag, push ------------------------------------------------------------------------
git add zombi2/__init__.py CHANGELOG.md
git commit -m "release: bump version to $VERSION"
git tag -a "$TAG" -m "$TAG"
git push origin main "$TAG"

# --- 4. publish the GitHub Release (triggers the PyPI upload) -------------------------------------
printf '%s\n' "$NOTES" | gh release create "$TAG" --title "$TAG" --notes-file -

echo
echo "Released $TAG. The PyPI publish is running:"
echo "  https://github.com/AADavin/zombi2/actions/workflows/release.yml"
