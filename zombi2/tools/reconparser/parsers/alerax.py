"""Parser for AleRax reconciliation output files.

AleRax is the successor to ALE (Amalgamated Likelihood Estimation),
providing scalable gene tree / species tree reconciliation with
support for duplications, transfers, and losses (DTL model).

This module provides two classes:
- AleRaxFamily: per-family parser (lazy, loads only what you ask for)
- AleRaxRun: run-level parser that wraps the full output directory
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re

import ete3
import pandas as pd


class AleRaxFamily:
    """
    Lazy parser for a single gene family within an AleRax run.

    Nothing is loaded at construction time. Data is parsed only when
    you call the corresponding getter, and then cached so the file
    is read at most once.

    Parameters
    ----------
    family_name : str
        The gene family identifier (e.g., "K00192").
    output_dir : str or Path
        Path to the AleRax output directory (the one containing
        ``reconciliations/``, ``species_trees/``, etc.).

    Examples
    --------
    >>> fam = AleRaxFamily("K00192", "output_UNIFORM_GLOBAL")
    >>> tree = fam.get_consensus_tree()
    >>> transfers = fam.get_transfers()
    >>> counts = fam.get_event_counts(sample=0)
    """

    def __init__(self, family_name: str, output_dir: str | Path):
        self.family_name = family_name
        self.output_dir = Path(output_dir)

        # Derived paths
        self._all_dir = self.output_dir / "reconciliations" / "all"
        self._summaries_dir = self.output_dir / "reconciliations" / "summaries"

        # Caches (None = not yet loaded)
        self._consensus_tree: Optional[ete3.Tree] = None
        self._sampled_gene_trees: Optional[List[ete3.Tree]] = None
        self._rec_uml_tree: Optional[ete3.Tree] = None
        self._summary_transfers: Optional[pd.DataFrame] = None
        self._summary_perspecies: Optional[pd.DataFrame] = None
        self._event_counts: Dict[int, Dict[str, int]] = {}
        self._sample_transfers: Dict[int, pd.DataFrame] = {}
        self._sample_perspecies: Dict[int, pd.DataFrame] = {}
        self._num_samples: Optional[int] = None

    # ------------------------------------------------------------------
    #  Path helpers
    # ------------------------------------------------------------------

    def _newick_path(self) -> Path:
        return self._all_dir / f"{self.family_name}.newick"

    def _rec_uml_path(self) -> Path:
        return self._all_dir / f"{self.family_name}.rec_uml"

    def _consensus_path(self) -> Path:
        # AleRax names these with a majority-rule threshold, e.g. _consensus_50
        candidates = sorted(self._summaries_dir.glob(
            f"{self.family_name}_consensus_*.newick"
        ))
        if candidates:
            return candidates[0]
        # Fallback to exact name
        return self._summaries_dir / f"{self.family_name}_consensus_50.newick"

    def _summary_transfers_path(self) -> Path:
        return self._summaries_dir / f"{self.family_name}_transfers.txt"

    def _summary_perspecies_path(self) -> Path:
        return self._summaries_dir / f"{self.family_name}_perspecies_eventcount.txt"

    def _sample_eventcount_path(self, sample: int) -> Path:
        return self._all_dir / f"{self.family_name}_eventcount_{sample}.txt"

    def _sample_transfers_path(self, sample: int) -> Path:
        return self._all_dir / f"{self.family_name}_transfers_{sample}.txt"

    def _sample_perspecies_path(self, sample: int) -> Path:
        return self._all_dir / f"{self.family_name}_perspecies_eventcount_{sample}.txt"

    def _xml_path(self, sample: int) -> Path:
        return self._all_dir / f"{self.family_name}_{sample}.xml"

    # ------------------------------------------------------------------
    #  Number of samples
    # ------------------------------------------------------------------

    def get_num_samples(self) -> int:
        """
        Return the number of sampled gene trees for this family.

        Determined by counting lines in the .newick file.

        Returns
        -------
        int
            Number of sampled reconciled gene trees.
        """
        if self._num_samples is not None:
            return self._num_samples

        nwk = self._newick_path()
        if not nwk.exists():
            raise FileNotFoundError(f"Newick file not found: {nwk}")

        with open(nwk, "r") as f:
            self._num_samples = sum(1 for line in f if line.strip())

        return self._num_samples

    # ------------------------------------------------------------------
    #  Consensus tree (from summaries/)
    # ------------------------------------------------------------------

    def get_consensus_tree(self) -> ete3.Tree:
        """
        Return the majority-rule consensus gene tree.

        The tree has support values (clade frequencies across samples)
        as internal-node support.

        Returns
        -------
        ete3.Tree
            Consensus gene tree with support values.

        Raises
        ------
        FileNotFoundError
            If no consensus tree file is found.
        """
        if self._consensus_tree is not None:
            return self._consensus_tree

        path = self._consensus_path()
        if not path.exists():
            raise FileNotFoundError(f"Consensus tree not found: {path}")

        with open(path, "r") as f:
            tree_str = f.read().strip()

        self._consensus_tree = ete3.Tree(tree_str, format=1)
        return self._consensus_tree

    # ------------------------------------------------------------------
    #  Sampled gene trees (from all/*.newick)
    # ------------------------------------------------------------------

    def get_sampled_gene_trees(self) -> List[ete3.Tree]:
        """
        Return all sampled gene trees for this family.

        Internal node names carry event-type labels (S, T, D).

        .. warning::
            This loads *all* samples into memory. For large families
            consider :meth:`get_sampled_gene_tree` to load one at a time.

        Returns
        -------
        List[ete3.Tree]
            All sampled gene trees.
        """
        if self._sampled_gene_trees is not None:
            return self._sampled_gene_trees

        nwk = self._newick_path()
        if not nwk.exists():
            raise FileNotFoundError(f"Newick file not found: {nwk}")

        trees: List[ete3.Tree] = []
        with open(nwk, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    trees.append(ete3.Tree(line, format=1))

        self._sampled_gene_trees = trees
        return self._sampled_gene_trees

    def get_sampled_gene_tree(self, sample: int) -> ete3.Tree:
        """
        Return a single sampled gene tree by index.

        This avoids loading all trees into memory.

        Parameters
        ----------
        sample : int
            Zero-based sample index.

        Returns
        -------
        ete3.Tree
            The requested sampled gene tree.
        """
        # If all trees are already cached, just index
        if self._sampled_gene_trees is not None:
            return self._sampled_gene_trees[sample]

        nwk = self._newick_path()
        if not nwk.exists():
            raise FileNotFoundError(f"Newick file not found: {nwk}")

        with open(nwk, "r") as f:
            for i, line in enumerate(f):
                if i == sample:
                    return ete3.Tree(line.strip(), format=1)

        raise IndexError(
            f"Sample index {sample} out of range for family {self.family_name}"
        )

    # ------------------------------------------------------------------
    #  Reconciled gene tree with annotations (from all/*.rec_uml)
    # ------------------------------------------------------------------

    def get_reconciled_gene_tree(self) -> ete3.Tree:
        """
        Return the annotated reconciled gene tree (.rec_uml format).

        Node names contain reconciliation annotations:

        - ``.T@donor->recipient`` for transfers
        - ``.D@species`` for duplications
        - ``.S`` for speciations
        - Species-node identifiers for the mapping

        Returns
        -------
        ete3.Tree
            Annotated reconciled gene tree.
        """
        if self._rec_uml_tree is not None:
            return self._rec_uml_tree

        path = self._rec_uml_path()
        if not path.exists():
            raise FileNotFoundError(f"rec_uml file not found: {path}")

        with open(path, "r") as f:
            tree_str = f.read().strip()

        self._rec_uml_tree = ete3.Tree(tree_str, format=1)
        return self._rec_uml_tree

    # ------------------------------------------------------------------
    #  Summary transfers (from summaries/)
    # ------------------------------------------------------------------

    def get_transfers(self) -> pd.DataFrame:
        """
        Return averaged transfer events across all samples.

        Returns
        -------
        pd.DataFrame
            Columns: ``from``, ``to``, ``freq``.
            Sorted by descending frequency.
        """
        if self._summary_transfers is not None:
            return self._summary_transfers

        path = self._summary_transfers_path()
        if not path.exists():
            raise FileNotFoundError(f"Summary transfers not found: {path}")

        self._summary_transfers = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["from", "to", "freq"],
            dtype={"from": str, "to": str, "freq": float},
            engine="python",
        )
        self._summary_transfers.sort_values(
            "freq", ascending=False, inplace=True, ignore_index=True
        )
        return self._summary_transfers

    def get_transfers_as_dict_list(self) -> List[Dict[str, float | str]]:
        """
        Convenience: transfers as a list of dicts.

        Returns
        -------
        List[Dict]
            Each dict has keys ``from``, ``to``, ``freq``.
        """
        return self.get_transfers().to_dict("records")

    # ------------------------------------------------------------------
    #  Summary per-species event counts (from summaries/)
    # ------------------------------------------------------------------

    def get_perspecies_events(self) -> pd.DataFrame:
        """
        Return per-species event counts averaged across samples.

        Returns
        -------
        pd.DataFrame
            Columns: ``species_label``, ``speciations``, ``duplications``,
            ``losses``, ``transfers``, ``presence``, ``origination``,
            ``copies``, ``singletons``.
        """
        if self._summary_perspecies is not None:
            return self._summary_perspecies

        path = self._summary_perspecies_path()
        if not path.exists():
            raise FileNotFoundError(
                f"Summary per-species events not found: {path}"
            )

        self._summary_perspecies = pd.read_csv(
            path, skipinitialspace=True
        )
        return self._summary_perspecies

    # ------------------------------------------------------------------
    #  Per-sample event counts (from all/*_eventcount_N.txt)
    # ------------------------------------------------------------------

    def get_event_counts(self, sample: int) -> Dict[str, int]:
        """
        Return event counts for a single sampled reconciliation.

        Parameters
        ----------
        sample : int
            Zero-based sample index.

        Returns
        -------
        Dict[str, int]
            Keys: ``S`` (speciations), ``SL`` (speciation-losses),
            ``D`` (duplications), ``DL`` (duplication-losses),
            ``T`` (transfers), ``TL`` (transfer-losses),
            ``L`` (losses), ``Leaf`` (leaves).
        """
        if sample in self._event_counts:
            return self._event_counts[sample]

        path = self._sample_eventcount_path(sample)
        if not path.exists():
            raise FileNotFoundError(f"Event count file not found: {path}")

        counts: Dict[str, int] = {}
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if ":" in line:
                    key, val = line.split(":", 1)
                    counts[key.strip()] = int(val.strip())

        self._event_counts[sample] = counts
        return counts

    def get_all_event_counts(self) -> pd.DataFrame:
        """
        Return event counts across all samples as a DataFrame.

        Returns
        -------
        pd.DataFrame
            One row per sample. Columns: ``S``, ``SL``, ``D``, ``DL``,
            ``T``, ``TL``, ``L``, ``Leaf``.
        """
        n = self.get_num_samples()
        rows = [self.get_event_counts(i) for i in range(n)]
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    #  Per-sample transfers (from all/*_transfers_N.txt)
    # ------------------------------------------------------------------

    def get_sample_transfers(self, sample: int) -> pd.DataFrame:
        """
        Return transfer events for a single sampled reconciliation.

        Parameters
        ----------
        sample : int
            Zero-based sample index.

        Returns
        -------
        pd.DataFrame
            Columns: ``from``, ``to``, ``count``.
        """
        if sample in self._sample_transfers:
            return self._sample_transfers[sample]

        path = self._sample_transfers_path(sample)
        if not path.exists():
            raise FileNotFoundError(
                f"Sample transfers file not found: {path}"
            )

        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["from", "to", "count"],
            dtype={"from": str, "to": str, "count": float},
            engine="python",
        )
        self._sample_transfers[sample] = df
        return df

    # ------------------------------------------------------------------
    #  Per-sample per-species events (from all/*_perspecies_eventcount_N.txt)
    # ------------------------------------------------------------------

    def get_sample_perspecies_events(self, sample: int) -> pd.DataFrame:
        """
        Return per-species event counts for a single sample.

        Parameters
        ----------
        sample : int
            Zero-based sample index.

        Returns
        -------
        pd.DataFrame
            Same schema as :meth:`get_perspecies_events`.
        """
        if sample in self._sample_perspecies:
            return self._sample_perspecies[sample]

        path = self._sample_perspecies_path(sample)
        if not path.exists():
            raise FileNotFoundError(
                f"Sample per-species events not found: {path}"
            )

        df = pd.read_csv(path, skipinitialspace=True)
        self._sample_perspecies[sample] = df
        return df

    # ------------------------------------------------------------------
    #  File existence check
    # ------------------------------------------------------------------

    def files_exist(self) -> Dict[str, bool]:
        """
        Check which output files exist for this family.

        Returns
        -------
        Dict[str, bool]
            Mapping of file type to existence.
        """
        return {
            "newick": self._newick_path().exists(),
            "rec_uml": self._rec_uml_path().exists(),
            "consensus_tree": self._consensus_path().exists(),
            "summary_transfers": self._summary_transfers_path().exists(),
            "summary_perspecies": self._summary_perspecies_path().exists(),
        }

    def __repr__(self) -> str:
        return (
            f"AleRaxFamily('{self.family_name}', "
            f"output_dir='{self.output_dir}')"
        )


# ======================================================================
#  Run-level parser
# ======================================================================


class AleRaxRun:
    """
    Lazy parser for a complete AleRax run directory.

    Provides access to run-level data (species trees, model parameters,
    likelihoods, global transfers, global per-species events, origins)
    as well as per-family data via :class:`AleRaxFamily` objects.

    Nothing is loaded at construction time — data is parsed on demand.

    Parameters
    ----------
    output_dir : str or Path
        Path to the AleRax output directory (the one containing
        ``reconciliations/``, ``species_trees/``, ``alerax.log``, etc.).

    Examples
    --------
    >>> run = AleRaxRun("output_UNIFORM_GLOBAL")
    >>> tree = run.get_species_tree()
    >>> params = run.get_model_parameters()
    >>> fam = run.get_family("K00192")
    >>> fam.get_transfers()
    """

    def __init__(self, output_dir: str | Path):
        self.output_dir = Path(output_dir)
        if not self.output_dir.is_dir():
            raise NotADirectoryError(
                f"AleRax output directory not found: {self.output_dir}"
            )

        # Caches
        self._families: Optional[List[str]] = None
        self._family_parsers: Dict[str, AleRaxFamily] = {}
        self._species_tree: Optional[ete3.Tree] = None
        self._starting_species_tree: Optional[ete3.Tree] = None
        self._model_parameters: Optional[pd.DataFrame] = None
        self._per_fam_likelihoods: Optional[pd.DataFrame] = None
        self._global_transfers: Optional[pd.DataFrame] = None
        self._global_perspecies: Optional[pd.DataFrame] = None
        self._run_info: Optional[Dict] = None
        self._fraction_missing: Optional[pd.DataFrame] = None
        self._per_species_coverage: Optional[pd.DataFrame] = None
        self._ccp_dimensions: Optional[pd.DataFrame] = None

    # ------------------------------------------------------------------
    #  Path helpers
    # ------------------------------------------------------------------

    def _species_tree_path(self, which: str = "inferred") -> Path:
        name = (
            "inferred_species_tree.newick"
            if which == "inferred"
            else "starting_species_tree.newick"
        )
        return self.output_dir / "species_trees" / name

    # ------------------------------------------------------------------
    #  Family discovery
    # ------------------------------------------------------------------

    def get_family_names(self) -> List[str]:
        """
        Return sorted list of gene family names in this run.

        Families are discovered from .newick files in
        ``reconciliations/all/``.

        Returns
        -------
        List[str]
            Sorted family identifiers (e.g., ``['K00192', 'K00193', ...]``).
        """
        if self._families is not None:
            return self._families

        all_dir = self.output_dir / "reconciliations" / "all"
        if not all_dir.is_dir():
            raise FileNotFoundError(
                f"reconciliations/all/ not found in {self.output_dir}"
            )

        self._families = sorted(
            p.stem for p in all_dir.glob("*.newick")
        )
        return self._families

    def get_family(self, family_name: str) -> AleRaxFamily:
        """
        Return a lazy :class:`AleRaxFamily` parser for one gene family.

        Parameters
        ----------
        family_name : str
            Gene family identifier (e.g., ``"K00192"``).

        Returns
        -------
        AleRaxFamily
        """
        if family_name not in self._family_parsers:
            self._family_parsers[family_name] = AleRaxFamily(
                family_name, self.output_dir
            )
        return self._family_parsers[family_name]

    def get_families(self) -> Dict[str, AleRaxFamily]:
        """
        Return a dict of all :class:`AleRaxFamily` parsers keyed by name.

        Returns
        -------
        Dict[str, AleRaxFamily]
        """
        return {name: self.get_family(name) for name in self.get_family_names()}

    # ------------------------------------------------------------------
    #  Species trees
    # ------------------------------------------------------------------

    def get_species_tree(self) -> ete3.Tree:
        """
        Return the inferred (optimised) species tree.

        Returns
        -------
        ete3.Tree
        """
        if self._species_tree is not None:
            return self._species_tree

        path = self._species_tree_path("inferred")
        if not path.exists():
            raise FileNotFoundError(f"Inferred species tree not found: {path}")

        with open(path, "r") as f:
            self._species_tree = ete3.Tree(f.read().strip(), format=1)
        return self._species_tree

    def get_starting_species_tree(self) -> ete3.Tree:
        """
        Return the starting (input) species tree.

        Returns
        -------
        ete3.Tree
        """
        if self._starting_species_tree is not None:
            return self._starting_species_tree

        path = self._species_tree_path("starting")
        if not path.exists():
            raise FileNotFoundError(f"Starting species tree not found: {path}")

        with open(path, "r") as f:
            self._starting_species_tree = ete3.Tree(f.read().strip(), format=1)
        return self._starting_species_tree

    # ------------------------------------------------------------------
    #  Model parameters
    # ------------------------------------------------------------------

    def get_model_parameters(self) -> pd.DataFrame:
        """
        Return optimised DTL rate parameters for each gene.

        Returns
        -------
        pd.DataFrame
            Columns: ``gene``, ``dup_rate``, ``loss_rate``,
            ``transfer_rate``.
        """
        if self._model_parameters is not None:
            return self._model_parameters

        path = self.output_dir / "model_parameters" / "model_parameters.txt"
        if not path.exists():
            raise FileNotFoundError(f"Model parameters not found: {path}")

        self._model_parameters = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["gene", "dup_rate", "loss_rate", "transfer_rate"],
            dtype={
                "gene": str,
                "dup_rate": float,
                "loss_rate": float,
                "transfer_rate": float,
            },
            engine="python",
        )
        return self._model_parameters

    # ------------------------------------------------------------------
    #  Per-family likelihoods
    # ------------------------------------------------------------------

    def get_per_family_likelihoods(self) -> pd.DataFrame:
        """
        Return per-family log-likelihoods.

        Returns
        -------
        pd.DataFrame
            Columns: ``family``, ``log_likelihood``.
            Sorted by ascending likelihood (worst first).
        """
        if self._per_fam_likelihoods is not None:
            return self._per_fam_likelihoods

        path = self.output_dir / "per_fam_likelihoods.txt"
        if not path.exists():
            raise FileNotFoundError(
                f"Per-family likelihoods not found: {path}"
            )

        self._per_fam_likelihoods = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["family", "log_likelihood"],
            dtype={"family": str, "log_likelihood": float},
            engine="python",
        )
        return self._per_fam_likelihoods

    def get_total_log_likelihood(self) -> float:
        """
        Return the sum of per-family log-likelihoods.

        Returns
        -------
        float
        """
        return float(self.get_per_family_likelihoods()["log_likelihood"].sum())

    # ------------------------------------------------------------------
    #  Global transfers (reconciliations/transfers.txt)
    # ------------------------------------------------------------------

    def get_transfers(self) -> pd.DataFrame:
        """
        Return transfer events aggregated across all families and samples.

        Returns
        -------
        pd.DataFrame
            Columns: ``from``, ``to``, ``score``.
        """
        if self._global_transfers is not None:
            return self._global_transfers

        path = self.output_dir / "reconciliations" / "transfers.txt"
        if not path.exists():
            raise FileNotFoundError(f"Global transfers not found: {path}")

        self._global_transfers = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["from", "to", "score"],
            dtype={"from": str, "to": str, "score": float},
            engine="python",
        )
        self._global_transfers.sort_values(
            "score", ascending=False, inplace=True, ignore_index=True
        )
        return self._global_transfers

    # ------------------------------------------------------------------
    #  Global per-species events (reconciliations/perspecies_eventcount.txt)
    # ------------------------------------------------------------------

    def get_perspecies_events(self) -> pd.DataFrame:
        """
        Return per-species event counts aggregated over all families.

        Returns
        -------
        pd.DataFrame
            Columns: ``species_label``, ``speciations``, ``duplications``,
            ``losses``, ``transfers``, ``presence``, ``origination``,
            ``copies``, ``singletons``.
        """
        if self._global_perspecies is not None:
            return self._global_perspecies

        path = (
            self.output_dir / "reconciliations" / "perspecies_eventcount.txt"
        )
        if not path.exists():
            raise FileNotFoundError(
                f"Global per-species event counts not found: {path}"
            )

        self._global_perspecies = pd.read_csv(path, skipinitialspace=True)
        return self._global_perspecies

    # ------------------------------------------------------------------
    #  Origins (reconciliations/origins/)
    # ------------------------------------------------------------------

    def get_origin(self, node_name: str) -> Dict[str, float]:
        """
        Return origination probabilities for a given species-tree node.

        Parameters
        ----------
        node_name : str
            Node identifier (e.g., ``"Node_a1001_a1048_0"``).

        Returns
        -------
        Dict[str, float]
            Mapping of source identifiers (including ``"vertical"``)
            to origination probability.
        """
        path = self.output_dir / "reconciliations" / "origins" / f"{node_name}.txt"
        if not path.exists():
            raise FileNotFoundError(f"Origin file not found: {path}")

        result: Dict[str, float] = {}
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if ":" in line:
                    key, val = line.split(":", 1)
                    result[key.strip()] = float(val.strip())
                elif "," in line:
                    key, val = line.split(",", 1)
                    result[key.strip()] = float(val.strip())

        return result

    def get_origin_node_names(self) -> List[str]:
        """
        Return sorted list of node names that have origin files.

        Returns
        -------
        List[str]
        """
        origins_dir = self.output_dir / "reconciliations" / "origins"
        if not origins_dir.is_dir():
            raise FileNotFoundError(
                f"Origins directory not found: {origins_dir}"
            )
        return sorted(p.stem for p in origins_dir.glob("*.txt"))

    # ------------------------------------------------------------------
    #  Fraction missing / coverage
    # ------------------------------------------------------------------

    def get_fraction_missing(self) -> pd.DataFrame:
        """
        Return per-species fraction of missing gene families.

        Returns
        -------
        pd.DataFrame
            Columns: ``species``, ``fraction_missing``.
        """
        if self._fraction_missing is not None:
            return self._fraction_missing

        path = self.output_dir / "fractionMissing.txt"
        if not path.exists():
            raise FileNotFoundError(f"fractionMissing.txt not found: {path}")

        rows = []
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("SPECIES"):
                    continue
                parts = line.split(":")
                if len(parts) == 2:
                    rows.append({
                        "species": parts[0].strip(),
                        "fraction_missing": float(parts[1].strip()),
                    })

        self._fraction_missing = pd.DataFrame(rows)
        return self._fraction_missing

    def get_per_species_coverage(self) -> pd.DataFrame:
        """
        Return per-species family coverage ratio.

        Returns
        -------
        pd.DataFrame
            Columns: ``species``, ``coverage``.
        """
        if self._per_species_coverage is not None:
            return self._per_species_coverage

        path = self.output_dir / "perSpeciesCoverage.txt"
        if not path.exists():
            raise FileNotFoundError(
                f"perSpeciesCoverage.txt not found: {path}"
            )

        self._per_species_coverage = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["species", "coverage"],
            dtype={"species": str, "coverage": float},
            engine="python",
        )
        return self._per_species_coverage

    # ------------------------------------------------------------------
    #  CCP dimensions
    # ------------------------------------------------------------------

    def get_ccp_dimensions(self) -> pd.DataFrame:
        """
        Return CCP dimension info for each family.

        Returns
        -------
        pd.DataFrame
            Columns: ``family``, ``num_nodes``, ``num_clades``.
        """
        if self._ccp_dimensions is not None:
            return self._ccp_dimensions

        path = self.output_dir / "ccpdim.txt"
        if not path.exists():
            raise FileNotFoundError(f"ccpdim.txt not found: {path}")

        self._ccp_dimensions = pd.read_csv(
            path,
            header=None,
            names=["family", "num_nodes", "num_clades"],
            dtype={"family": str, "num_nodes": int, "num_clades": int},
        )
        return self._ccp_dimensions

    # ------------------------------------------------------------------
    #  Run info (parsed from alerax.log)
    # ------------------------------------------------------------------

    def get_run_info(self) -> Dict:
        """
        Parse run metadata from ``alerax.log``.

        Returns
        -------
        Dict
            Keys include ``version``, ``command``, ``rec_model``,
            ``parametrization``, ``transfer_constraint``,
            ``num_families``, ``num_species``, ``total_genes``,
            ``gene_tree_samples``, ``origination``.
        """
        if self._run_info is not None:
            return self._run_info

        path = self.output_dir / "alerax.log"
        if not path.exists():
            raise FileNotFoundError(f"alerax.log not found: {path}")

        info: Dict = {}
        with open(path, "r") as f:
            text = f.read()

        # Version
        m = re.search(r"AleRax (v[\d.]+)", text)
        if m:
            info["version"] = m.group(1)

        # Command
        m = re.search(r"AleRax was called as follow:\n(.+)", text)
        if m:
            info["command"] = m.group(1).strip()

        # Reconciliation model
        m = re.search(r"Reconciliation model:\s*(\S+)", text)
        if m:
            info["rec_model"] = m.group(1)

        # Parametrization
        m = re.search(r"Model parametrization:\s*(.+)", text)
        if m:
            info["parametrization"] = m.group(1).strip()

        # Transfer constraint
        m = re.search(r"Transfer constraints:\s*(.+)", text)
        if m:
            info["transfer_constraint"] = m.group(1).strip()

        # Num families
        m = re.search(r"Number of gene families:\s*(\d+)", text)
        if m:
            info["num_families"] = int(m.group(1))

        # Num species
        m = re.search(r"Number of species:\s*(\d+)", text)
        if m:
            info["num_species"] = int(m.group(1))

        # Total genes
        m = re.search(r"Total number of genes:\s*(\d+)", text)
        if m:
            info["total_genes"] = int(m.group(1))

        # Gene tree samples
        m = re.search(r"gene trees to sample:\s*(\d+)", text)
        if m:
            info["gene_tree_samples"] = int(m.group(1))

        # Origination
        m = re.search(r"Origination strategy:\s*(.+)", text)
        if m:
            info["origination"] = m.group(1).strip()

        self._run_info = info
        return self._run_info

    # ------------------------------------------------------------------
    #  File existence check
    # ------------------------------------------------------------------

    def files_exist(self) -> Dict[str, bool]:
        """
        Check which run-level output files/directories exist.

        Returns
        -------
        Dict[str, bool]
        """
        return {
            "alerax.log": (self.output_dir / "alerax.log").exists(),
            "species_trees": (self.output_dir / "species_trees").is_dir(),
            "model_parameters": (
                self.output_dir / "model_parameters" / "model_parameters.txt"
            ).exists(),
            "per_fam_likelihoods": (
                self.output_dir / "per_fam_likelihoods.txt"
            ).exists(),
            "reconciliations": (
                self.output_dir / "reconciliations"
            ).is_dir(),
            "global_transfers": (
                self.output_dir / "reconciliations" / "transfers.txt"
            ).exists(),
            "global_perspecies": (
                self.output_dir
                / "reconciliations"
                / "perspecies_eventcount.txt"
            ).exists(),
            "origins": (
                self.output_dir / "reconciliations" / "origins"
            ).is_dir(),
            "ccpdim": (self.output_dir / "ccpdim.txt").exists(),
            "fractionMissing": (
                self.output_dir / "fractionMissing.txt"
            ).exists(),
            "perSpeciesCoverage": (
                self.output_dir / "perSpeciesCoverage.txt"
            ).exists(),
        }

    def __repr__(self) -> str:
        try:
            n = len(self.get_family_names())
        except FileNotFoundError:
            n = "?"
        return f"AleRaxRun('{self.output_dir}', families={n})"
