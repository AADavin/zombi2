"""The gallery's example code calls the live ``zombi2`` API — check it still resolves.

The `extant()` -> `extant_leaves()` rename broke every gallery figure that colours extinct lineages,
and nothing caught it: the examples call the library, but no test ran them, so a published example
tracebacked for anyone who copied it.

Rendering them in CI is the wrong guard. The gallery draws with Phylustrator, and its examples target
Phylustrator's *development* API — no ``pip``-installable Phylustrator renders them — so a render job
would couple this repo's CI to a sibling package's release cadence and be red for reasons that have
nothing to do with ``zombi2``.

So this mocks Phylustrator out and runs each example's render function. The example's ``zombi2`` calls
then execute for real — ``simulate_species_tree(...)``, ``ct.extant_leaves()``, ``res.node_values`` —
while every drawing call is a no-op on the mock. A mock cannot raise ``AttributeError`` or
``TypeError``; a real ``zombi2`` object whose method was renamed or whose signature changed does. So
those two exception types escaping a render *are* a zombi2 API break — deterministically, with no
rendering and no Phylustrator installed. It is the lean canary the bug needed.

Scope: the species / traits / joining examples, which call zombi2's **Python API** directly. The
genome and sequence examples reach zombi2 through its **CLI** (``h.ordered_run()`` and friends shell
out to ``zombi2 ... run``) and read the results back through Phylustrator — that path is guarded by
``test_cli.py``, and mocking Phylustrator underneath it would only feed mock objects into the glue and
manufacture failures. One of those helpers also downloads a genome, which must not run in CI.
"""

from __future__ import annotations

import importlib
import pathlib
import sys
from unittest.mock import MagicMock

import pytest

GALLERY = pathlib.Path(__file__).resolve().parent.parent / "gallery"
# The examples that call the zombi2 Python API directly (see the module docstring for why not genomes
# / sequences). This is where the extant()->extant_leaves() rename bit, and where a mock cleanly
# isolates the zombi2 calls.
_MODULES = ("species", "traits", "joining")
# A zombi2 rename surfaces as one of these on a real result object; a mocked Phylustrator call cannot
# raise them. Anything else out of a render (a FileNotFoundError from a composite reading back a PNG
# the mock never wrote, say) is the drawing layer, not our concern.
_ZOMBI2_BREAK = (AttributeError, TypeError)


class _Draw(MagicMock):
    """A stand-in for a Phylustrator figure: every attribute and call returns another one, ``+``
    composes layers, and it iterates empty — so a gallery line that loops over a drawing return does
    not raise a spurious ``TypeError``. Only the zombi2 calls are under test."""

    def __add__(self, other):
        return self

    def __iter__(self):
        return iter(())


@pytest.fixture
def gallery_examples(monkeypatch, tmp_path):
    # Phylustrator and matplotlib are drawing, not simulation — mock them (and their submodules) in
    # sys.modules before the gallery imports them, so `import phylustrator as ph` and helpers' matplotlib
    # imports bind the mocks. That keeps this test dependency-free: it needs only zombi2 + pytest (the
    # [dev] extra), so it rides in the existing test job with no Phylustrator, matplotlib, or rendering.
    # helpers itself stays real — it defines the Example/EXAMPLES the test iterates and the compute
    # helpers (LTT, node_values, ...) that run on the real zombi2 objects.
    for name in ("phylustrator", "phylustrator.trees", "phylustrator.genomes", "phylustrator.zombi",
                 "matplotlib", "matplotlib.pyplot", "matplotlib.image", "matplotlib.cm",
                 "matplotlib.colors"):
        monkeypatch.setitem(sys.modules, name, _Draw())
    monkeypatch.syspath_prepend(str(GALLERY))
    monkeypatch.chdir(tmp_path)                       # any file a render writes lands here, not the repo

    examples = []
    for name in (*_MODULES, "helpers"):
        sys.modules.pop(name, None)                  # fresh import so the mock is the one bound
    for name in _MODULES:
        module = importlib.import_module(name)
        examples += [(name, ex) for ex in module.EXAMPLES]
    return examples


def _break(render, out) -> str | None:
    """Run one render with Phylustrator mocked; return a message if a zombi2 call broke, else None."""
    try:
        render(out)
    except _ZOMBI2_BREAK as e:
        return f"{type(e).__name__}: {e}"
    except Exception:
        return None                                  # drawing / file I/O on the mock — not a zombi2 break
    return None


def test_every_gallery_example_calls_a_current_zombi2_api(gallery_examples, tmp_path):
    assert len(gallery_examples) >= 9, "the gallery shrank — is the import finding the examples?"

    failures = [f"{mod}.{ex.id}: {msg}"
                for mod, ex in gallery_examples
                if (msg := _break(ex.render, str(tmp_path / f"{ex.id}.png"))) is not None]
    assert not failures, ("gallery examples call a zombi2 API that no longer exists — the published "
                          "figures would traceback for anyone who copied them:\n  " + "\n  ".join(failures))


def test_the_guard_catches_a_renamed_zombi2_method(gallery_examples, tmp_path, monkeypatch):
    """A guard on the guard: rename a method out from under the examples and confirm it is flagged, so
    a green result means the checks ran rather than that nothing did."""
    from zombi2.tree import Tree
    assert hasattr(Tree, "extant_leaves")            # the method the examples rely on today
    monkeypatch.delattr(Tree, "extant_leaves")       # simulate the next accidental rename

    caught = [f"{mod}.{ex.id}"
              for mod, ex in gallery_examples
              if (msg := _break(ex.render, str(tmp_path / f"{ex.id}.png"))) and "extant_leaves" in msg]
    assert caught, "removing Tree.extant_leaves should have tripped at least one gallery example"
