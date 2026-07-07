# ZOMBI 1 — genome-step (`Gm`) rate semantics & benchmarking notes

Everything below is derived from reading the ZOMBI 1 source (a local ZOMBI v1
checkout, pointed at by `$ZOMBI1_DIR`; nothing in that repo was modified).

## 1. Parameter-file syntax (`GenomeParameters.tsv` / `SpeciesTreeParameters.tsv`)

Parsing lives in `AuxiliarFunctions.py`:

- `read_parameters(file)` (line 51) reads `KEY<TAB>VALUE` or `KEY<space>VALUE`
  lines, skipping `#` comments and blanks, into a `dict[str, str]`.
- `obtain_value(value)` (line 107) interprets the `PREFIX:PAYLOAD` mini-language
  used by most rate values. The prefix picks a *distribution*, and a value is
  **drawn from it every time `obtain_value` is called**:
  - `f:X`   → fixed value `X` (e.g. `f:5` = always 5).
  - `n:m;s` → Normal(m, s), abs.
  - `l:m;s` → LogNormal(m, s), abs.
  - `u:a;b` → Uniform(a, b), abs.
  - `g:p`   → **Geometric(p)** — returns an *integer ≥ 1*. Used by the
    `*_EXTENSION` params: `g:1` ⇒ Geometric(1) ⇒ almost always 1 gene per event.
  - `e:λ`   → Exponential(λ).
- `prepare_genome_parameters(parameters)` (line 238) coerces a fixed set of keys
  to `int` (`PROFILES`, `EVENTS_PER_BRANCH`, `GENE_TREES`, `RECONCILED_TREES`,
  `VERBOSE`, `MIN_GENOME_SIZE`, `SEED`, `SCALE_TREE`, …) and leaves the rate
  strings (`DUPLICATION`, `TRANSFER`, `LOSS`, `INVERSION`, `TRANSPOSITION`,
  `ORIGINATION`, extensions) as raw strings — they are parsed later by
  `obtain_value`. So `DUPLICATION f:5` is stored as the string `"f:5"` and only
  turned into a number when a family's rate is generated.

## 2. What the D/T/L/O rates *mean* in `Gm` mode

`Gm` = "genome mode m" = **per-gene-family DTL rates** (`GenomeSimulator.run_m`,
line 591). The flow:

1. `fill_genome(family_rates=True)` (line 374) builds the root genome of
   `INITIAL_GENOME_SIZE` genes. Each gene is created by `make_origination(...,
   family_mode=True)` (line 1347), which calls `generate_new_rates()` (line 802).
   `generate_new_rates` just calls `obtain_value` on `DUPLICATION`, `TRANSFER`,
   `LOSS` — so **each gene family gets its own D, T, L drawn from those
   distributions**. With `f:X` (fixed) every family gets exactly X.
2. The Gillespie clock is `get_time_to_next_event_family_mode()` (line 1320):

       total_rate = Σ_lineages Σ_genes (family.D + family.T + family.L)
                  + n_lineages · (INVERSION + TRANSPOSITION + ORIGINATION)
       Δt ~ Exponential(1 / total_rate)

   So **D/T/L are per-gene, per-time-unit rates**, summed over *every gene in
   every currently-alive lineage*. **INVERSION / TRANSPOSITION / ORIGINATION are
   per-lineage, per-time-unit rates** (they scale with the number of live
   lineages, not the genome size).
3. `evolve_genomes_m(time)` (line 1015) picks a lineage ∝ its total rate, then an
   event ∝ (D_sum, T_sum, L_sum, I, P, O):
   - **D** duplicates a gene (same family → that family's D/T/L now apply to one
     more gene → the family's contribution to `total_rate` grows).
   - **T** transfers a gene to another live lineage. If `REPLACEMENT_TRANSFER=0`
     the transfer is **additive** (recipient genome *grows*); if `=1` it replaces
     a gene (no net growth). This is the single biggest knob for runaway growth.
   - **L** deletes a gene (guarded by `MIN_GENOME_SIZE`: losses that would push a
     genome below the floor are skipped).
   - **O** originates a brand-new family (fresh D/T/L from `generate_new_rates`)
     and inserts one gene → genome and family count both grow.

### Why the DEFAULTS explode
Default `GenomeParameters.tsv`: `D=5, T=12, L=4, INITIAL_GENOME_SIZE=100`,
`REPLACEMENT_TRANSFER=0.5`. Per gene, birth ≈ D + (additive fraction of T) ≈
5 + 6 = 11 ≫ death L = 4. Every gene reproduces ~2.75× faster than it dies, so
the genome (and hence `total_rate`) grows without bound and `run_m` never
reaches the tips. **Verified**: on a 20-tip tree, `Gm` with the defaults did not
finish in 45 s (killed), whereas the benchmark regime below finished in 0.69 s.

## 3. Units & the branch-length (tree-scale) interaction

- Rates are **events per gene (or per lineage) per unit of branch length** on the
  species tree. `run_m` walks species-tree events in tree-time (`Events.tsv`
  timestamps) and races the genome Gillespie clock against them, so the *same
  rate* produces more events on a taller tree.
- `SCALE_RATES True` (`GenomeSimulator.__init__`, line 36) only matters when
  `RATE_FILE` supplies *empirical* rates — it divides them by the crown length.
  With `RATE_FILE False` (our case) the D/T/L strings are used as-is, so **the
  effective per-gene rate equals the number in the file** and the total number of
  genome events scales with the tree-time integral of (genes × rate).
- **Consequence for matching ZOMBI2.** ZOMBI2's comparable regime is
  `BirthDeath(λ=1.0, μ=0.3)`, **fixed tree age 2.0**, `D=0.2 T=0.1 L=0.25
  O=0.5`, size 20. ZOMBI 1 with `STOPPING_RULE=1` grows the tree until it has *N*
  extant tips, so under a Yule process with λ=1.0 the tree **height ≈ ln(N)**
  (≈ 4.6 at N=100, ≈ 9.2 at N=10 000) rather than a fixed 2.0. This is an
  **unavoidable difference**: matching the tip count and the per-event rates
  (task-approved) means the ZOMBI 1 tree is deeper, so per tip it does somewhat
  *more* genome work than ZOMBI2 at age 2.0. For an order-of-magnitude
  wall-clock-vs-size ceiling this is fine and is noted here for honesty.

## 4. Controlling the number of extant tips (`T` mode)

From `SpeciesTreeParameters.tsv` + `SpeciesTreeSimulator.run()` (line 37) and
`prepare_species_tree_parameters` (`AuxiliarFunctions.py` line 175):

- `STOPPING_RULE` — `0` = stop at `TOTAL_TIME`; **`1` = stop when the number of
  *alive* lineages equals `TOTAL_LINEAGES`** (this is the clean tip-count knob).
- `TOTAL_LINEAGES` — target number of alive lineages when `STOPPING_RULE=1`.
- `TOTAL_TIME` — stop time when `STOPPING_RULE=0`.
- `SPECIATION` / `EXTINCTION` — birth/death rates (`f:λ`, etc.).
- `MIN_LINEAGES` / `MAX_LINEAGES` — abort/retry guards; the tree build retries up
  to 100 times if it dies out or blows past `MAX_LINEAGES`.
- `TURNOVER` — only used by the `Tp` (lineage-profile) mode.

**Tip-count strategy used by the harness:** `STOPPING_RULE=1`,
`TOTAL_LINEAGES=N`, `SPECIATION=f:1.0`, **`EXTINCTION=f:0`**. With zero extinction
the process is pure-birth, every stopped lineage is *extant*, so
`T/ExtantTree.nwk` has **exactly N leaves** and there are no failed-tree retries
(important at N=10 000+). `MAX_LINEAGES` is set very high so the stopping rule,
not the guard, ends the run.

## 5. Chosen non-exploding, ZOMBI2-comparable regime (`GenomeParameters_bench.tsv`)

    DUPLICATION  f:0.2      (per gene, per time)   <- ZOMBI2 D
    TRANSFER     f:0.1      (per gene, per time)   <- ZOMBI2 T
    LOSS         f:0.25     (per gene, per time)   <- ZOMBI2 L
    ORIGINATION  f:0.5      (per lineage, per time)<- ZOMBI2 O
    INVERSION    f:0
    TRANSPOSITION f:0
    INITIAL_GENOME_SIZE 20                          <- ZOMBI2 initial size
    REPLACEMENT_TRANSFER 1   (transfers replace, so they do NOT grow genomes)
    MIN_GENOME_SIZE 5
    PROFILES 1, GENE_TREES 0, RECONCILED_TREES 0, EVENTS_PER_BRANCH 0, VERBOSE 0

Rationale / stability: with `REPLACEMENT_TRANSFER=1`, per-gene birth = D = 0.2 <
per-gene death = L = 0.25, so families do **not** run away (honours the
"loss ≥ duplication + transfer" guidance — here L ≥ D and transfer is
size-neutral). Origination adds families at O=0.5 per lineage; genome size stays
in the low tens. Gene-tree / reconciliation / events-per-branch writing is turned
OFF so the measured time is the genome *simulation* step, not disk-heavy
post-processing — the fair thing to compare against ZOMBI2's genome sim.

Differences from ZOMBI2 that cannot be removed (documented for the figure):
- ZOMBI 1 tree **height grows as ln(N)** (Yule, λ=1) vs ZOMBI2's fixed age 2.0.
- ZOMBI 1 draws one fixed D/T/L per family (here identical for all families
  because we use `f:`), matching ZOMBI2's per-family constant-rate spirit.
- ZOMBI2's μ=0.3 (some extinction) vs our μ=0 (clean tip count). Extinction only
  changes the *complete* tree; the genome sim cost is dominated by extant tips,
  which we match exactly.

## 6. Measured scaling & practical ceiling (180s hard timeout)

| target N | actual extant tips | Gm seconds | families | status  |
|---------:|-------------------:|-----------:|---------:|:--------|
| 20       | 20                 | 0.69       | 29       | ok      |
| 100      | 100                | 1.05       | 58       | ok      |
| 300      | 300                | 3.64       | 166      | ok      |
| 1000     | 1000               | 49.6       | 550      | ok      |
| 1200     | 1200               | 71.4       | 610      | ok      |
| 1500     | 1500               | >180       | —        | timeout |
| 3000     | 3000               | >180       | —        | timeout |

Cost is strongly **super-linear** in tip count (300→1000 tips ≈ 14× the time for
≈ 3.3× the tips), because the Gillespie total rate scales with (genes × lineages)
and the Yule tree height grows as ln(N), so total genome work grows much faster
than N. **Practical ceiling ≈ 1200–1400 extant tips** under a couple-of-minutes
budget (1200 tips finishes in ~71 s; 1500 tips does not finish in 180 s). ZOMBI 1
is pure Python with no vectorisation, so this is the honest order-of-magnitude
limit for the genome step.

## 7. Files
- `harness.py` — the timing harness (see its docstring).
- `GenomeParameters_bench.tsv` — the regime above.
- `SpeciesTreeParameters_template.tsv` — T-mode template (`__NTIPS__`, `__SEED__`).
- `GenomeParameters_default_reference.tsv` — copy of ZOMBI 1 defaults (the
  exploding regime), kept for reference only.
- `results.json` — raw measurement lines.
