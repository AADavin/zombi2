# Examples

Ready-to-run ZOMBI2 pipelines. Each parameters file in `parameters/` drives a full
run; the data files it needs are **not** committed (they are large and public), so
you fetch them once into `data/` with the commands below.

## `parameters/ecoli.toml` — a full *E. coli* genome down a species tree

A three-step pipeline: a 20-tip dated species tree, the whole *E. coli* K-12 MG1655
genome evolving along it at the **nucleotide** resolution (real annotation + real DNA),
and DNA evolved down every gene tree. Event rates are set deliberately high so the
duplications, losses, inversions and transfers are plentiful and easy to see in the logs.

### 1. Fetch the initial genome (once)

The run starts from the real *E. coli* K-12 MG1655 assembly
([GCF_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)).
From this `examples/` folder:

```bash
mkdir -p data
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"
curl -sL "$BASE/GCF_000005845.2_ASM584v2_genomic.gff.gz" | gunzip > data/ecoli.gff   # annotation (~2.4 MB)
curl -sL "$BASE/GCF_000005845.2_ASM584v2_genomic.fna.gz" | gunzip > data/ecoli.fna   # DNA        (~4.6 MB)
```

`ecoli.gff` gives the gene layout; `ecoli.fna` is the DNA, so the real MG1655 sequence
descends to the tips. Drop `ecoli.fna` (and the `fasta =` line in the TOML) and the run
still works — the letters then come from the substitution model instead of real E. coli.

### 2. Run the pipeline

Name the same run directory each time; every step reads the same parameters file:

```bash
zombi2 species   ecoli_run/ --params parameters/ecoli.toml   # dated tree, 20 tips   (~instant)
zombi2 genomes   ecoli_run/ --params parameters/ecoli.toml   # full genome           (~2 s)
zombi2 sequences ecoli_run/ --params parameters/ecoli.toml   # DNA down each gene tree (~20 s)
```

The `sequences` step assembles a full genome for **every node** of the tree, so it writes
a few hundred MB by default. To keep only the assembled genomes, add
`write = ["genomes"]` to the `[sequences]` table.

### A note on the rates

Every rate in the `[genomes]` table is **per lineage, per unit of branch time** — that is
the scope at the nucleotide resolution (unlike the per-copy rates of the unordered/ordered
resolutions). The number says how often a lineage does the event; the matching
`*-length` says how many base pairs it touches. See the manual, *Genomes III: nucleotide*.
