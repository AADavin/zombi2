# References

Internal bibliography for the ZOMBI2 models — the source for the future PDF manual's
bibliography and a lookup while writing the docs. **Every entry was verified against Crossref /
the DOI record** (not written from memory) on 2026-07-04; the one item that could not be
confirmed is marked **⚠ UNVERIFIED**. Each line gives a suggested BibTeX `key`, the DOI, and the
ZOMBI2 feature it supports. This page is intentionally kept out of the site navigation.

*Verification method:* each work was matched via the Crossref REST API
(`api.crossref.org/works?query.bibliographic=…`) and its DOI record; author list, year, title and
venue were confirmed against that record. Where a starting guess differed from the record, the
verified values were used.

## The ZOMBI lineage

- **Davín, Tricou, Tannier, de Vienne & Szöllősi (2020)** — *Zombi: a phylogenetic simulator of
  trees, genomes and sequences that accounts for dead lineages.* Bioinformatics 36(4):1286–1288.
  [10.1093/bioinformatics/btz710](https://doi.org/10.1093/bioinformatics/btz710). `davin2020zombi`
  · the original ZOMBI simulator.
- **Davín, Woodcroft, Soo, Morel, Murali, Schrempf, Clark, Álvarez-Carretero, Boussau, Moody,
  Szánthó, Richy, Pisani, Hemp, Fischer, Donoghue, Spang, Hugenholtz, Williams & Szöllősi (2025)**
  — *A geological timescale for bacterial evolution and oxygen adaptation.* Science
  388(6742):eadp1853. [10.1126/science.adp1853](https://doi.org/10.1126/science.adp1853).
  `davin2025goe` · the GOE-based bacterial-tree dating that motivates trait-linked gene families.

## Species trees & the reconstructed process

- **Nee, May & Harvey (1994)** — *The reconstructed evolutionary process.* Phil. Trans. R. Soc. B
  344(1309):305–311. [10.1098/rstb.1994.0068](https://doi.org/10.1098/rstb.1994.0068).
  `nee1994reconstructed` · the reconstructed birth–death process.
- **Hartmann, Wong & Stadler (2010)** — *Sampling Trees from Evolutionary Models.* Syst. Biol.
  59(4):465–476. [10.1093/sysbio/syq026](https://doi.org/10.1093/sysbio/syq026).
  `hartmann2010sampling` · the backward, tip-conditioned species-tree sampler.
- **Stadler (2009)** — *On incomplete sampling under birth–death models and connections to the
  sampling-based coalescent.* J. Theor. Biol. 261(1):58–66.
  [10.1016/j.jtbi.2009.07.018](https://doi.org/10.1016/j.jtbi.2009.07.018). `stadler2009incomplete`
  · constant-rate birth–death and incomplete sampling (TreeSim).
- **Lambert & Stadler (2013)** — *Birth–death models and coalescent point processes: the shape and
  probability of reconstructed phylogenies.* Theor. Popul. Biol. 90:113–128.
  [10.1016/j.tpb.2013.10.002](https://doi.org/10.1016/j.tpb.2013.10.002). `lambert2013birthdeath`
  · the coalescent-point-process view behind the backward sampler.

## Diversification models

- **Stadler & Bonhoeffer (2013)** — *Uncovering epidemiological dynamics in heterogeneous host
  populations using phylogenetic methods.* Phil. Trans. R. Soc. B 368(1614):20120198.
  [10.1098/rstb.2012.0198](https://doi.org/10.1098/rstb.2012.0198). `stadler2013skyline`
  · the birth–death skyline (episodic rates).
- **Stadler, Kühnert, Bonhoeffer & Drummond (2013)** — *Birth–death skyline plot reveals temporal
  changes of epidemic spread in HIV and hepatitis C virus (HCV).* PNAS 110(1):228–233.
  [10.1073/pnas.1207965110](https://doi.org/10.1073/pnas.1207965110). `stadler2013birthdeathskyline`
  · the skyline plot / rate-shift estimation.
- **Höhna, May & Moore (2016)** — *TESS: an R package for efficiently simulating phylogenetic trees
  and performing Bayesian inference of lineage diversification rates.* Bioinformatics 32(5):789–791.
  [10.1093/bioinformatics/btv651](https://doi.org/10.1093/bioinformatics/btv651). `hohna2016tess`
  · episodic BD simulation/inference; the mass-extinction pulse formulation.
- **Stadler (2011)** — *Mammalian phylogeny reveals recent diversification rate shifts.* PNAS
  108(15):6187–6192. [10.1073/pnas.1016876108](https://doi.org/10.1073/pnas.1016876108).
  `stadler2011mammalian` · instantaneous mass-extinction survival pulses (TreeSim).
- **Rabosky & Lovette (2008)** — *Density-dependent diversification in North American wood
  warblers.* Proc. R. Soc. B 275(1649):2363–2371.
  [10.1098/rspb.2008.0630](https://doi.org/10.1098/rspb.2008.0630). `rabosky2008densitydependent`
  · diversity-dependent diversification.
- **Etienne, Haegeman, Stadler, Aze, Pearson, Purvis & Phillimore (2012)** —
  *Diversity-dependence brings molecular phylogenies closer to agreement with the fossil record.*
  Proc. R. Soc. B 279(1732):1300–1309.
  [10.1098/rspb.2011.1439](https://doi.org/10.1098/rspb.2011.1439). `etienne2012diversitydependence`
  · the diversity-dependent likelihood (DDD).
- **Maliet, Hartig & Morlon (2019)** — *A model with many small shifts for estimating
  species-specific diversification rates.* Nat. Ecol. Evol. 3(7):1086–1092.
  [10.1038/s41559-019-0908-0](https://doi.org/10.1038/s41559-019-0908-0). `maliet2019clads`
  · ClaDS (per-lineage diversification-rate shifts).
- **Gavryushkina, Welch, Stadler & Drummond (2014)** — *Bayesian Inference of Sampled Ancestor
  Trees for Epidemiology and Fossil Calibration.* PLoS Comput. Biol. 10(12):e1003919.
  [10.1371/journal.pcbi.1003919](https://doi.org/10.1371/journal.pcbi.1003919).
  `gavryushkina2014sampledancestor` · sampled ancestors / fossilized birth–death.
- **Heath, Huelsenbeck & Stadler (2014)** — *The fossilized birth–death process for coherent
  calibration of divergence-time estimates.* PNAS 111(29):E2957–E2966.
  [10.1073/pnas.1319091111](https://doi.org/10.1073/pnas.1319091111). `heath2014fossilized`
  · the fossilized birth–death (FBD) process.
- **Andréoletti, Zwaans, Warnock, Aguirre-Fernández, Barido-Sottani, Gupta, Stadler & Manceau
  (2022)** — *The Occurrence Birth–Death Process for Combined-Evidence Analysis in Macroevolution
  and Epidemiology.* Syst. Biol. 71(6):1440–1452.
  [10.1093/sysbio/syac037](https://doi.org/10.1093/sysbio/syac037). `andreoletti2022occurrence`
  · the occurrence birth–death process (roadmap).
- **Louca & Pennell (2020)** — *Extant timetrees are consistent with a myriad of diversification
  histories.* Nature 580(7804):502–505.
  [10.1038/s41586-020-2176-1](https://doi.org/10.1038/s41586-020-2176-1). `louca2020identifiability`
  · identifiability limits of time-varying birth–death (a caveat).

## Ghost lineages

- **Szöllősi, Tannier, Lartillot & Daubin (2013)** — *Lateral Gene Transfer from the Dead.* Syst.
  Biol. 62(3):386–397. [10.1093/sysbio/syt003](https://doi.org/10.1093/sysbio/syt003).
  `szollosi2013lgtdead` · why extinct/unsampled ("ghost") lineages matter for gene transfer.

## Gene families & phylogenetic profiles

- **Pellegrini, Marcotte, Thompson, Eisenberg & Yeates (1999)** — *Assigning protein functions by
  comparative genome analysis: protein phylogenetic profiles.* PNAS 96(8):4285–4288.
  [10.1073/pnas.96.8.4285](https://doi.org/10.1073/pnas.96.8.4285). `pellegrini1999profiles`
  · the founding phylogenetic-profiles idea (functional co-occurrence).

## Trait evolution

- **Felsenstein (1985)** — *Phylogenies and the Comparative Method.* Am. Nat. 125(1):1–15.
  [10.1086/284325](https://doi.org/10.1086/284325). `felsenstein1985comparative`
  · Brownian motion / independent contrasts.
- **Lewis (2001)** — *A Likelihood Approach to Estimating Phylogeny from Discrete Morphological
  Character Data.* Syst. Biol. 50(6):913–925.
  [10.1080/106351501753462876](https://doi.org/10.1080/106351501753462876). `lewis2001mk`
  · the Mk model (`--model mk`).
- **Pagel (1994)** — *Detecting correlated evolution on phylogenies: a general method for the
  comparative analysis of discrete characters.* Proc. R. Soc. B 255(1342):37–45.
  [10.1098/rspb.1994.0006](https://doi.org/10.1098/rspb.1994.0006). `pagel1994correlated`
  · correlated evolution of binary characters (`CorrelatedBinary`).
- **Pagel (1999)** — *Inferring the historical patterns of biological evolution.* Nature
  401(6756):877–884. [10.1038/44766](https://doi.org/10.1038/44766). `pagel1999inferring`
  · the λ/κ/δ branch transforms.
- **Felsenstein (2012)** — *A Comparative Method for Both Discrete and Continuous Characters Using
  the Threshold Model.* Am. Nat. 179(2):145–156.
  [10.1086/663681](https://doi.org/10.1086/663681). `felsenstein2012threshold`
  · the threshold (liability) model (`--model threshold`).
- **Huelsenbeck, Nielsen & Bollback (2003)** — *Stochastic Mapping of Morphological Characters.*
  Syst. Biol. 52(2):131–158.
  [10.1080/10635150390192780](https://doi.org/10.1080/10635150390192780). `huelsenbeck2003stochastic`
  · the stochastic character map (per-branch discrete histories).
- **Nielsen (2002)** — *Mapping Mutations on Phylogenies.* Syst. Biol. 51(5):729–739.
  [10.1080/10635150290102393](https://doi.org/10.1080/10635150290102393). `nielsen2002mapping`
  · the mutational-mapping precursor to stochastic maps.
- **Beaulieu, O'Meara & Donoghue (2013)** — *Identifying Hidden Rate Changes in the Evolution of a
  Binary Morphological Character.* Syst. Biol. 62(5):725–737.
  [10.1093/sysbio/syt034](https://doi.org/10.1093/sysbio/syt034). `beaulieu2013hidden`
  · hidden-rates model / corHMM (`HiddenStateMk`).
- **Hansen (1997)** — *Stabilizing Selection and the Comparative Analysis of Adaptation.* Evolution
  51(5):1341–1351. [10.1111/j.1558-5646.1997.tb01457.x](https://doi.org/10.1111/j.1558-5646.1997.tb01457.x).
  `hansen1997stabilizing` · the Ornstein–Uhlenbeck model (`--model ou`).
- **Butler & King (2004)** — *Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive
  Evolution.* Am. Nat. 164(6):683–695. [10.1086/426002](https://doi.org/10.1086/426002).
  `butler2004phylogenetic` · OU adaptive evolution (ouch); multi-optimum OU.
- **Beaulieu, Jhwueng, Boettiger & O'Meara (2012)** — *Modeling Stabilizing Selection: Expanding
  the Ornstein–Uhlenbeck Model of Adaptive Evolution.* Evolution 66(8):2369–2383.
  [10.1111/j.1558-5646.2012.01619.x](https://doi.org/10.1111/j.1558-5646.2012.01619.x). `beaulieu2012ouwie`
  · OUwie / variable-regime OU (`MultiOptimumOU`).
- **Harmon et al. (2010)** — *Early Bursts of Body Size and Shape Evolution Are Rare in Comparative
  Data.* Evolution 64(8):2385–2396.
  [10.1111/j.1558-5646.2010.01025.x](https://doi.org/10.1111/j.1558-5646.2010.01025.x). `harmon2010earlyburst`
  · the Early Burst / ACDC model (`--model eb`).
- **Clavel, Escarguel & Merceron (2015)** — *mvMORPH: an R package for fitting multivariate
  evolutionary models to morphometric data.* Methods Ecol. Evol. 6(11):1311–1319.
  [10.1111/2041-210X.12420](https://doi.org/10.1111/2041-210X.12420). `clavel2015mvmorph`
  · multivariate/correlated continuous traits (`MultivariateBrownian`/`MultivariateOU`).
- **Ree & Smith (2008)** — *Maximum Likelihood Inference of Geographic Range Evolution by Dispersal,
  Local Extinction, and Cladogenesis.* Syst. Biol. 57(1):4–14.
  [10.1080/10635150701883881](https://doi.org/10.1080/10635150701883881). `ree2008dec`
  · the DEC biogeography model (`--model dec`).

## Coevolution — state-dependent diversification (SSE)

- **Maddison, Midford & Otto (2007)** — *Estimating a Binary Character's Effect on Speciation and
  Extinction.* Syst. Biol. 56(5):701–710.
  [10.1080/10635150701607033](https://doi.org/10.1080/10635150701607033). `maddison2007bisse`
  · BiSSE.
- **FitzJohn, Maddison & Otto (2009)** — *Estimating Trait-Dependent Speciation and Extinction
  Rates from Incompletely Resolved Phylogenies.* Syst. Biol. 58(6):595–611.
  [10.1093/sysbio/syp067](https://doi.org/10.1093/sysbio/syp067). `fitzjohn2009trait`
  · the multi-state SSE likelihood (MuSSE precursor).
- **FitzJohn (2012)** — *Diversitree: comparative phylogenetic analyses of diversification in R.*
  Methods Ecol. Evol. 3(6):1084–1092.
  [10.1111/j.2041-210X.2012.00234.x](https://doi.org/10.1111/j.2041-210X.2012.00234.x). `fitzjohn2012diversitree`
  · diversitree / MuSSE.
- **FitzJohn (2010)** — *Quantitative Traits and Diversification.* Syst. Biol. 59(6):619–633.
  [10.1093/sysbio/syq053](https://doi.org/10.1093/sysbio/syq053). `fitzjohn2010quasse`
  · QuaSSE.
- **Beaulieu & O'Meara (2016)** — *Detecting Hidden Diversification Shifts in Models of
  Trait-Dependent Speciation and Extinction.* Syst. Biol. 65(4):583–601.
  [10.1093/sysbio/syw022](https://doi.org/10.1093/sysbio/syw022). `beaulieu2016hisse`
  · HiSSE.
- **Maddison & FitzJohn (2015)** — *The Unsolved Challenge to Phylogenetic Correlation Tests for
  Categorical Characters.* Syst. Biol. 64(1):127–136.
  [10.1093/sysbio/syu070](https://doi.org/10.1093/sysbio/syu070). `maddison2015unreplicated`
  · the "unreplicated evolution" caveat for SSE / correlation tests.

## Sequence evolution & rate variation

- **Thorne, Kishino & Painter (1998)** — *Estimating the rate of evolution of the rate of molecular
  evolution.* Mol. Biol. Evol. 15(12):1647–1657.
  [10.1093/oxfordjournals.molbev.a025892](https://doi.org/10.1093/oxfordjournals.molbev.a025892).
  `thorne1998autocorrelated` · the autocorrelated relaxed molecular clock.
- **⚠ UNVERIFIED — the GTDB discrete-bin within-branch clock.** `zombi2/rate_variation.py` describes
  its model as "the discrete-bin model from the GTDB archaea study" but stores no citation, and it
  could not be pinned to a single confirmed paper. Most likely **Davín et al. (2025)** *Science*
  (`davin2025goe`, same group) — **needs your confirmation**. Other candidates checked and *not*
  matched: Moody et al. 2022 (*eLife* 11:e66695, [10.7554/eLife.66695](https://doi.org/10.7554/eLife.66695));
  Wang & Luo 2025 (*Syst. Biol.* 74(4):639–655, [10.1093/sysbio/syae071](https://doi.org/10.1093/sysbio/syae071),
  a different group). · supports: `RateVariation` (discrete-bin lineage clock).

## Inference — Approximate Bayesian Computation

- **Beaumont, Zhang & Balding (2002)** — *Approximate Bayesian Computation in Population Genetics.*
  Genetics 162(4):2025–2035.
  [10.1093/genetics/162.4.2025](https://doi.org/10.1093/genetics/162.4.2025). `beaumont2002abc`
  · the ABC regression adjustment (`--regression-adjust`).
- **Toni, Welch, Strelkowa, Ipsen & Stumpf (2009)** — *Approximate Bayesian computation scheme for
  parameter inference and model selection in dynamical systems.* J. R. Soc. Interface 6(31):187–202.
  [10.1098/rsif.2008.0172](https://doi.org/10.1098/rsif.2008.0172). `toni2009abcsmc`
  · ABC-SMC (`--smc`).

---

*Status: 42 verified against Crossref/DOI, 1 flagged UNVERIFIED (the GTDB rate-variation clock).
Open confirmation for you: which paper the GTDB discrete-bin clock comes from. The same 42 entries
are in [`references.bib`](references.bib) (BibTeX, same keys) for the LaTeX/pandoc PDF build —
validated with `biber --tool`.*
