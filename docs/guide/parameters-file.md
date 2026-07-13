# Parameters files

Every simulation command takes its settings as command-line flags. For a reproducible or shared
configuration you can instead collect those settings in a **`--params` file** — a small
[TOML](https://toml.io) file whose keys are the command's option names. It is the modern successor
to ZOMBI1's `*Parameters.tsv` files.

`species`, `genomes` and `sequence` accept `--params FILE`.

```bash
zombi2 species --params species.toml -o run/
```

## How it applies

The file supplies **defaults**; any flag you also pass on the command line **overrides** the file.
So the file captures the model, and you can still tweak one knob per run:

```bash
# use species.toml but override the seed for this replicate
zombi2 species --params species.toml --seed 99 -o run/rep99
```

The required **I/O paths stay on the command line** — `-o`/`--out`, `-t`/`--tree`, `--genomes`.
The parameters file is for the *model*, not for where files go.

## The format

Keys are the command's long option names (write them with hyphens or underscores — both work), and
values are native TOML scalars or arrays:

```toml
# species.toml
birth = 1.0
death = 0.3
tips  = 20
age   = 5.0
seed  = 42
```

```bash
zombi2 species --params species.toml -o run/
# identical to:
zombi2 species --birth 1.0 --death 0.3 --tips 20 --age 5.0 --seed 42 -o run/
```

A few conventions:

- **List-valued options** take a TOML array — e.g. `write = ["profiles", "trees"]` for
  `zombi2 genomes`. For an option that *can* take several values but you only want one (the
  episodic `--birth`/`--death`, say), a bare scalar works too: `birth = 1.0` is read as `[1.0]`.
- **Unknown keys are an error**, so a typo is caught rather than silently ignored.

## One file for a whole pipeline

Give the file a table per subcommand and reuse it across the run — each command reads its own
`[section]`:

```toml
# pipeline.toml
[species]
birth = 1.0
death = 0.3
tips  = 20
age   = 5.0
seed  = 42

[genomes]
dup   = 0.2
trans = 0.1
loss  = 0.25
orig  = 0.5
write = ["profiles", "trees"]
seed  = 42
```

```bash
zombi2 species  --params pipeline.toml -o run/
zombi2 genomes  --params pipeline.toml -t run/species_tree.nwk -o run/
```

(A file with no `[section]` tables is read as a flat set of parameters for whichever command you
pass it to.)
