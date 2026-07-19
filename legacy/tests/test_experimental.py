"""The experimental staging namespace: its warning machinery, and its isolation
from the stable core surface."""
from __future__ import annotations

import warnings

import pytest

import zombi2.experimental as ex


def test_warn_experimental_emits_experimental_warning():
    ex._warned.discard("DemoModel")
    with pytest.warns(ex.ExperimentalWarning, match="DemoModel"):
        ex.warn_experimental("DemoModel")


def test_warn_experimental_is_quiet_after_the_first_call():
    ex._warned.discard("DemoModel")
    with pytest.warns(ex.ExperimentalWarning):
        ex.warn_experimental("DemoModel")
    with warnings.catch_warnings():
        warnings.simplefilter("error")     # any warning would now raise
        ex.warn_experimental("DemoModel")  # already warned -> stays silent


def test_experimental_is_not_reexported_into_the_core():
    import zombi2

    # the staging namespace is deliberately kept out of the stable top-level API
    assert "experimental" not in getattr(zombi2, "__all__", [])
    assert not hasattr(zombi2, "warn_experimental")
