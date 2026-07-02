# ZOMBI2

A Python library for simulating **species-tree** and **gene-family** evolution in
phylogenetics. A ground-up redesign of [ZOMBI](https://github.com/AADavin/Zombi).

- The **species tree** is simulated *backward* in time (reconstructed birth–death,
  conditioned on the number of extant species).
- **Gene families** are simulated *forward* in time along that fixed species tree, with
  duplication, transfer, loss and origination events.

The architecture is deliberately interface-first so that additions — genome
representations (an ordered-chromosome model with inversions/transpositions ships today),
rate models, gene-family coupling (a Potts model), ghost lineages — arrive as new
subclasses rather than a rewrite.

**Documentation:** full guides and API reference live in `docs/` (build with
`pip install -e ".[docs]" && mkdocs serve`). Start with `docs/index.md`,
`docs/quickstart.md`, and the user guide under `docs/guide/`.

## Quick start

Two symmetric steps — simulate the species tree, then the gene families along it.
Models are objects; simple cases have a keyword shorthand.

```python
import zombi2 as z

# 1. species tree (backward). Yule(birth) == BirthDeath(birth, death=0).
tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3),
                               n_tips=20, age=5.0, seed=1)

# 2a. gene families where every family shares the same D/T/L (shorthand):
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_size=40, seed=42)

# 2b. ...or every family draws its OWN D/T/L from distributions (ZOMBI-1 style):
genomes = z.simulate_genomes(tree, z.FamilySampledRates(
    duplication=z.Gamma(2, 0.06), transfer=z.Exponential(0.08),
    loss=z.Gamma(2, 0.07), origination=0.5), initial_size=40, seed=42)

print(tree.to_newick())
print(genomes.profiles.matrix)          # families x extant-species copy numbers
complete, extant = genomes.gene_trees()["1"]   # per-family reconstructed gene trees
genomes.write("out/")                   # trees, event tables, transfers, profiles
```

`genomes.write("out/")` produces ZOMBI-1-style output: `species_tree.nwk`,
per-family event tables (`gene_family_events/`), reconstructed **complete** and
**extant** gene trees (`gene_trees/`), `Transfers.tsv`, `Gene_family_summary.tsv`, and
the presence/copy-number matrices (`Profiles.tsv` / `Presence.tsv`).

**More models.** Species trees also support **episodic/skyline** rates and incomplete
sampling (`z.EpisodicBirthDeath(birth=[...], death=[...], shifts=[...], sampling_fraction=ρ)`);
gene families also support **genome-wise** rates (`z.GenomeWiseRates`, constant per-genome
rather than per-copy) and an ordered-chromosome genome with inversions/transpositions
(`z.OrderedGenome`). See the notebooks in `examples/` and the guides in `docs/`.

**Transfer mechanics.** Pass a `TransferModel` to control how a transfer resolves:

```python
z.simulate_genomes(tree, transfer=0.3, ..., transfers=z.TransferModel(
    replacement=0.2,        # 0 = additive (adds a copy); 1 = replacement (adds + removes one, net 0)
    distance_decay=2.0,     # None = uniform recipient; larger = more phylogenetically local transfers
    allow_self=True,        # a self-transfer = a duplication (lets you run transfer/loss only)
))
```

**Bounding gene-family growth.** Both duplication *and* transfer create copies, so a family
can grow like `e^{(d−l)t}` without bound. The principled control is `max_family_size` — a
hard ceiling across all events, given either as an absolute integer or (recommended) a
**fraction of the number of species**; an over-cap additive transfer is turned into a
replacement:

```python
z.simulate_genomes(tree, duplication=0.5, transfer=0.2, loss=0.1, origination=0.3,
                   max_family_size=0.5)   # cap = round(0.5 * N_species)
```

A softer, duplication-only alternative is `carrying_capacity=K` on the rate model (logistic
density dependence — family size settles around `K`).

Rate models are the flexible knob. `z.UniformRates(...)` and `z.FamilySampledRates(...)`
ship today; new models (genome-wise rates, gene-family coupling / Potts) are subclasses
of `z.RateModel` that need no change to the simulator. Distribution arguments accept a
built-in (`z.Gamma`, `z.Exponential`, `z.LogNormal`, `z.Uniform`, `z.Fixed`), any
`scipy.stats` frozen distribution, or a callable `rng -> float`.

**Large-scale profiles (optional Rust engine).** When you only need the presence/copy-number
matrix — not the event log or gene trees — an optional native engine runs the forward
Gillespie in Rust over per-family counts (~50× faster; the 10 000-tip profile matrix in
~0.4 s). Build it once (`pip install maturin && cd rust && maturin build --release -i python3
&& pip install --force-reinstall rust/target/wheels/*.whl`), then:

```python
if z.rust_available():
    profiles = z.simulate_profiles_fast(tree, duplication=0.05, transfer=0.03, loss=0.1,
                                        origination=0.5, initial_size=200,
                                        max_family_size=0.3, seed=42)  # -> ProfileMatrix
```

Need the full event log and gene trees, just faster? `z.simulate_genomes_fast(...)` tracks gene
lineages in Rust and returns a complete `Genomes` (with `.event_log`, `.gene_trees()`,
`.write()`) — a drop-in for `simulate_genomes` (~3× at 10k tips). And for large datasets on
disk, `z.simulate_and_write_fast(tree, "out/", ...)` simulates, reconstructs gene trees, **and
writes the whole ZOMBI-1 output in Rust** — ~10× vs Python simulate + write. All three fast
paths cover the built-in `UnorderedGenome` + `UniformRates` model; the pure-Python
`simulate_genomes` stays the default. See `docs/guide/rust-fast-path.md`.

### CLI

A thin wrapper over the library with three subcommands, mirroring the two-step design —
simulate a species tree, simulate gene families along one, or do both in one run:

```bash
# species tree only (backward birth–death) -> out/species_tree.nwk
zombi2 species --birth 1 --death 0.3 --tips 5000 --age 5 --seed 1 -o out/

# gene families along a supplied Newick tree (your own, or one from `species`)
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --max-family-size 0.5 --seed 42 -o out/

# species tree, then gene families along it, in one run
zombi2 all --birth 1 --death 0.2 --tips 20 --age 5 \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/
```

`species` and `all` build the tree; `genomes` reads one from `--tree` (any Newick file),
so `zombi2 species … -o out/` then `zombi2 genomes --tree out/species_tree.nwk … -o out/`
is the split form of `zombi2 all`. `species` writes `species_tree.nwk`; `genomes` and `all`
write the full ZOMBI-1-style output described above.

| Option | Commands | Meaning |
| --- | --- | --- |
| `--birth` / `--death` | `species`, `all` | speciation / extinction rate (`--death` defaults to 0 = Yule) |
| `--tips` / `--age` | `species`, `all` | number of extant species N / tree age |
| `--age-type {crown,stem}` | `species`, `all` | interpret `--age` as crown (default) or stem age |
| `--tree` / `-t` | `genomes` | input species tree in Newick format |
| `--dup` `--trans` `--loss` `--orig` | `genomes`, `all` | per-copy duplication / transfer / loss / origination rates |
| `--initial-size` | `genomes`, `all` | number of gene families seeded at the root (default 20) |
| `--max-family-size` | `genomes`, `all` | growth cap — integer = absolute, decimal = fraction of N (e.g. `0.5`) |
| `--fast` | `genomes`, `all` | use the Rust engine (same full output, much faster; see below) |
| `--profiles-only` | `genomes`, `all` | with `--fast`, write only the profile matrices |
| `--seed` | all | RNG seed for reproducibility |
| `-o` / `--out` | all | output directory |

**`--fast` (Rust).** Routes `genomes`/`all` through the optional Rust engine
(`simulate_and_write_fast`): it simulates, reconstructs the gene trees, and writes the **full
ZOMBI-1 output** — the same files as the default — entirely in Rust (~10× faster at scale).
Add `--profiles-only` to write just `species_tree.nwk` + `Profiles.tsv`/`Presence.tsv` (the
even-faster `simulate_profiles_fast`, no event log or gene trees). Requires the compiled
`zombi2_core` extension (see the Rust section above); without it, `--fast` exits with a build
hint.

```bash
zombi2 all --birth 1 --tips 5000 --age 5 --dup 0.2 --loss 0.25 --orig 0.5 --fast -o out/
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --fast --profiles-only -o out/
```

Run `zombi2 <command> --help` for the full list. The CLI covers the common uniform-rate
case; for family-sampled or genome-wise rates, transfer mechanics, ordered genomes, or
replicate parallelism, use the Python API above.

## Development

```bash
uv pip install -e ".[dev]"   # or: pip install -e ".[dev]"
pytest
```
