"""ZOMBI2 command-line interface.

Split into one module per subcommand (``species``, ``genomes``, ``traits``, ``sequences``,
``coevolve``, ``tools``, ``experimental``) over a shared :mod:`~zombi2.cli.framework`; :mod:`main`
assembles the parser and dispatches. ``main`` is the console entry point (``zombi2 = zombi2.cli:main``).

A few private helpers are re-exported here because the test-suite imports them from ``zombi2.cli``.
"""
from zombi2.cli.main import main

# Backward-compatible re-exports for `from zombi2.cli import <helper>` used by tests.
from zombi2.cli.framework import _int_or_float  # noqa: F401
from zombi2.cli.genomes import _extension_from_mean_length  # noqa: F401

__all__ = ["main"]
