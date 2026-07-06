# zombi2_core

The compiled Rust engine behind [ZOMBI2](https://github.com/AADavin/zombi2).

It runs ZOMBI2's **built-in** gene-family model (the order-free
`UnorderedGenome` with shared duplication/transfer/loss/origination rates),
which has no pure-Python fallback so that a given `seed` is reproducible against
a single engine. Every other model in ZOMBI2 — species trees, traits,
sequences, and the flexible genome rate models — runs in pure Python.

You normally do not install this package directly. It is a dependency of
`zombi2`:

```bash
pip install zombi2
```

Binary wheels are published for Linux, macOS, and Windows and work on
CPython 3.10 and newer (built against the stable ABI). If no wheel matches your
platform, `pip` builds it from source, which requires a
[Rust toolchain](https://rustup.rs).

Licensed under GPL-3.0-or-later, the same as ZOMBI2.
