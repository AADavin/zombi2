"""Parser for ALE (Amalgamated Likelihood Estimation) output files."""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import ete3
import pandas as pd
import re


class ALEParser:
    """
    Parser for ALE reconciliation output files.

    ALE produces several output files:
    - .ucons_tree: Consensus gene tree in Newick format
    - .uTs: Transfer events with frequencies (tab-separated)
    - .uml_rec: Maximum likelihood reconciliation with extensive data:
        * Species tree
        * Reconciled gene trees (variable number, default 100)
        * ML rates (Duplications, Transfers, Losses)
        * Log-likelihood
        * Per-branch statistics table

    Parameters
    ----------
    base_path : str or Path
        Path to ALE output files. Can be:
        - Base path without extension (e.g., "results.ale")
        - Path to a specific file (e.g., "results.ale.ucons_tree")

    Examples
    --------
    >>> parser = ALEParser("results.ale")
    >>> gene_tree = parser.get_consensus_tree()
    >>> transfers = parser.get_transfers()
    >>> ml_rates = parser.get_ml_rates()
    >>> logl = parser.get_log_likelihood()
    >>> branch_stats = parser.get_branch_statistics()
    """

    def __init__(self, base_path: str | Path):
        """Initialize ALE parser with base path to output files."""
        self.base_path = Path(base_path)

        # If a specific file was provided, extract the base path
        base_str = str(self.base_path)
        for ext in ['.ucons_tree', '.uTs', '.uml_rec']:
            if base_str.endswith(ext):
                self.base_path = Path(base_str[:-len(ext)])
                break

        # Define file paths
        self.ucons_tree_path = Path(str(self.base_path) + '.ucons_tree')
        self.uts_path = Path(str(self.base_path) + '.uTs')
        self.uml_rec_path = Path(str(self.base_path) + '.uml_rec')

        # Cache for parsed data
        self._consensus_tree: Optional[ete3.Tree] = None
        self._transfers: Optional[pd.DataFrame] = None
        self._reconciled_tree: Optional[ete3.Tree] = None
        self._reconciled_gene_trees: Optional[List[ete3.Tree]] = None
        self._ml_rates: Optional[Dict[str, float]] = None
        self._log_likelihood: Optional[float] = None
        self._branch_statistics: Optional[pd.DataFrame] = None
        self._summary_statistics: Optional[Dict[str, float]] = None

        # Cache for full uml_rec parse
        self._uml_rec_parsed: bool = False

    def get_consensus_tree(self) -> ete3.Tree:
        """
        Parse and return the consensus gene tree.

        Returns
        -------
        ete3.Tree
            The consensus gene tree with support values

        Raises
        ------
        FileNotFoundError
            If .ucons_tree file doesn't exist
        """
        if self._consensus_tree is not None:
            return self._consensus_tree

        if not self.ucons_tree_path.exists():
            raise FileNotFoundError(
                f"Consensus tree file not found: {self.ucons_tree_path}"
            )

        with open(self.ucons_tree_path, 'r') as f:
            # Skip header comment line
            lines = f.readlines()
            tree_line = None
            for line in lines:
                if line.strip() and not line.startswith('#'):
                    tree_line = line.strip()
                    break

            if tree_line is None:
                raise ValueError("No tree found in consensus tree file")

            # Parse the Newick tree (format=1 for standard format with support values)
            self._consensus_tree = ete3.Tree(tree_line, format=1)

        return self._consensus_tree

    def get_transfers(self) -> pd.DataFrame:
        """
        Parse and return transfer events with frequencies.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: 'from', 'to', 'freq'
            - 'from': Source node/leaf name
            - 'to': Target node/leaf name
            - 'freq': Transfer frequency (float between 0 and 1)

        Raises
        ------
        FileNotFoundError
            If .uTs file doesn't exist
        """
        if self._transfers is not None:
            return self._transfers

        if not self.uts_path.exists():
            raise FileNotFoundError(
                f"Transfer file not found: {self.uts_path}"
            )

        # Read the tab-separated file, skipping the header comment
        self._transfers = pd.read_csv(
            self.uts_path,
            sep='\t',
            comment='#',
            names=['from', 'to', 'freq'],
            skipinitialspace=True,
            dtype={'from': str, 'to': str, 'freq': float}
        )

        # Reset index to ensure clean integer index
        self._transfers = self._transfers.reset_index(drop=True)

        return self._transfers

    def get_transfers_as_dict_list(self) -> List[Dict[str, float | str]]:
        """
        Get transfers as a list of dictionaries.

        Convenient format for visualization libraries that accept
        transfer data as list of dicts.

        Returns
        -------
        List[Dict]
            List of transfer dictionaries with keys: 'from', 'to', 'freq'

        Examples
        --------
        >>> parser = ALEParser("results.ale")
        >>> transfers = parser.get_transfers_as_dict_list()
        >>> # Each transfer is a dict: {'from': 'node1', 'to': 'node2', 'freq': 0.5}
        """
        df = self.get_transfers()
        return df.to_dict('records')

    def _parse_uml_rec_file(self):
        """
        Parse the complete .uml_rec file and cache all data.

        This is called internally by the various get_* methods to avoid
        parsing the file multiple times.
        """
        if self._uml_rec_parsed:
            return

        if not self.uml_rec_path.exists():
            raise FileNotFoundError(
                f"Reconciliation file not found: {self.uml_rec_path}"
            )

        with open(self.uml_rec_path, 'r') as f:
            lines = f.readlines()

        # Parse different sections
        line_idx = 0

        # Skip header comment
        while line_idx < len(lines) and (not lines[line_idx].strip() or lines[line_idx].startswith('#')):
            line_idx += 1

        # Line 3: Species tree (S:)
        if line_idx < len(lines) and lines[line_idx].startswith('S:'):
            tree_str = lines[line_idx].split('S:', 1)[1].strip()
            self._reconciled_tree = ete3.Tree(tree_str, format=1)
            line_idx += 1

        # Skip empty lines
        while line_idx < len(lines) and not lines[line_idx].strip():
            line_idx += 1

        # Parse metadata section (lines 6-11 approximately)
        # Line: "Input ale from:\t<filename>"
        while line_idx < len(lines) and not lines[line_idx].strip().startswith('('):
            line = lines[line_idx].strip()

            # Parse log-likelihood
            if line.startswith('>logl:'):
                logl_str = line.split(':', 1)[1].strip()
                self._log_likelihood = float(logl_str)

            # Parse ML rates header (skip it)
            elif line.startswith('rate of'):
                pass

            # Parse ML rates values
            elif line.startswith('ML'):
                parts = line.split('\t')
                if len(parts) >= 4:
                    self._ml_rates = {
                        'duplications': float(parts[1]),
                        'transfers': float(parts[2]),
                        'losses': float(parts[3])
                    }

            # Parse number of reconciled gene trees
            elif 'reconciled G-s:' in line:
                match = re.search(r'(\d+)\s+reconciled', line)
                if match:
                    num_trees = int(match.group(1))
                    # Note: We'll parse the trees below

            line_idx += 1

        # Parse reconciled gene trees (start with "(")
        reconciled_trees = []
        while line_idx < len(lines) and lines[line_idx].strip().startswith('('):
            tree_str = lines[line_idx].strip()
            try:
                # These trees have reconciliation annotations, format=1
                tree = ete3.Tree(tree_str, format=1)
                reconciled_trees.append(tree)
            except:
                # If parsing fails, skip this tree
                pass
            line_idx += 1

        self._reconciled_gene_trees = reconciled_trees

        # Skip empty lines and look for summary statistics
        while line_idx < len(lines) and not lines[line_idx].strip().startswith('#'):
            line_idx += 1

        # Parse summary statistics (lines with "# of\tDuplications\tTransfers...")
        if line_idx < len(lines) and '# of' in lines[line_idx]:
            line_idx += 1  # Skip header
            if line_idx < len(lines) and lines[line_idx].strip().startswith('Total'):
                parts = lines[line_idx].strip().split('\t')
                if len(parts) >= 5:
                    self._summary_statistics = {
                        'total_duplications': float(parts[1]),
                        'total_transfers': float(parts[2]),
                        'total_losses': float(parts[3]),
                        'total_speciations': float(parts[4])
                    }
                line_idx += 1

        # Skip to branch statistics header line
        # ALE v0.4 header: "# of\tDuplications\tTransfers\tLosses\tOriginations\tcopies"
        # ALE v1.0 header: "# of\tDuplications\tTransfers\tLosses\tOriginations\tcopies\tsingletons\textinction_prob\tpresence\tLL"
        while line_idx < len(lines) and not lines[line_idx].strip().startswith('S_'):
            # Detect column header to determine format
            if line_idx < len(lines) and lines[line_idx].startswith('# of'):
                header_line = lines[line_idx].strip()
                self._branch_has_singletons = 'singletons' in header_line.lower()
            line_idx += 1

        # Parse branch statistics table
        branch_data = []
        while line_idx < len(lines):
            line = lines[line_idx].strip()
            if line.startswith('S_terminal_branch') or line.startswith('S_internal_branch'):
                parts = line.split('\t')
                if len(parts) >= 7:
                    # Clean branch_id: ALE v1.0 uses "name(idx)" for terminals
                    # e.g. "a10(0)" → "a10"
                    branch_id = parts[1]
                    paren_match = re.match(r'^(.+?)\(\d+\)$', branch_id)
                    if paren_match:
                        branch_id = paren_match.group(1)

                    row = {
                        'branch_type': parts[0],
                        'branch_id': branch_id,
                        'duplications': float(parts[2]),
                        'transfers': float(parts[3]),
                        'losses': float(parts[4]),
                        'originations': float(parts[5]),
                        'copies': float(parts[6]),
                    }

                    # ALE v1.0 extra columns
                    if len(parts) >= 11:
                        row['singletons'] = float(parts[7])
                        row['extinction_prob'] = float(parts[8])
                        row['presence'] = float(parts[9])
                        row['branch_LL'] = float(parts[10])

                    branch_data.append(row)
            line_idx += 1

        if branch_data:
            self._branch_statistics = pd.DataFrame(branch_data)

        self._uml_rec_parsed = True

    def get_reconciled_tree(self) -> ete3.Tree:
        """
        Parse and return the reconciled species tree.

        Returns
        -------
        ete3.Tree
            The reconciled species tree

        Raises
        ------
        FileNotFoundError
            If .uml_rec file doesn't exist
        """
        if self._reconciled_tree is None:
            self._parse_uml_rec_file()

        if self._reconciled_tree is None:
            raise ValueError("No species tree found in reconciliation file")

        return self._reconciled_tree

    def get_reconciled_gene_trees(self) -> List[ete3.Tree]:
        """
        Parse and return all reconciled gene trees.

        The .uml_rec file contains multiple reconciled gene trees
        (default 100, but this can vary). Each tree contains reconciliation
        annotations showing duplications (D@) and transfers (T@).

        Returns
        -------
        List[ete3.Tree]
            List of reconciled gene trees with annotations

        Raises
        ------
        FileNotFoundError
            If .uml_rec file doesn't exist
        """
        if self._reconciled_gene_trees is None:
            self._parse_uml_rec_file()

        return self._reconciled_gene_trees

    def get_ml_rates(self) -> Dict[str, float]:
        """
        Parse and return the ML rates for Duplications, Transfers, and Losses.

        Returns
        -------
        Dict[str, float]
            Dictionary with keys: 'duplications', 'transfers', 'losses'
            Each value is the maximum likelihood rate estimate

        Raises
        ------
        FileNotFoundError
            If .uml_rec file doesn't exist
        """
        if self._ml_rates is None:
            self._parse_uml_rec_file()

        if self._ml_rates is None:
            raise ValueError("ML rates not found in reconciliation file")

        return self._ml_rates

    def get_log_likelihood(self) -> float:
        """
        Parse and return the log-likelihood of the reconciliation.

        Returns
        -------
        float
            The log-likelihood value

        Raises
        ------
        FileNotFoundError
            If .uml_rec file doesn't exist
        """
        if self._log_likelihood is None:
            self._parse_uml_rec_file()

        if self._log_likelihood is None:
            raise ValueError("Log-likelihood not found in reconciliation file")

        return self._log_likelihood

    def get_summary_statistics(self) -> Dict[str, float]:
        """
        Parse and return summary statistics of reconciliation events.

        Returns
        -------
        Dict[str, float]
            Dictionary with keys:
            - 'total_duplications': Total number of duplication events
            - 'total_transfers': Total number of transfer events
            - 'total_losses': Total number of loss events
            - 'total_speciations': Total number of speciation events

        Raises
        ------
        FileNotFoundError
            If .uml_rec file doesn't exist
        """
        if self._summary_statistics is None:
            self._parse_uml_rec_file()

        if self._summary_statistics is None:
            raise ValueError("Summary statistics not found in reconciliation file")

        return self._summary_statistics

    def get_branch_statistics(self) -> pd.DataFrame:
        """
        Parse and return per-branch reconciliation statistics.

        Supports both ALE v0.4 (7 columns) and ALE v1.0 (11 columns).
        Terminal branch IDs are cleaned: "a10(0)" → "a10".

        Returns
        -------
        pd.DataFrame
            DataFrame with columns (always present):
            - 'branch_type': 'S_terminal_branch' or 'S_internal_branch'
            - 'branch_id': Branch identifier (leaf name or internal node ID)
            - 'duplications': Number of duplication events on this branch
            - 'transfers': Number of transfer events on this branch
            - 'losses': Number of loss events on this branch
            - 'originations': Number of origination events on this branch
            - 'copies': Number of gene copies on this branch

            Additional columns (ALE v1.0 only):
            - 'singletons': Number of singleton events on this branch
            - 'extinction_prob': Extinction probability on this branch
            - 'presence': Presence probability on this branch
            - 'branch_LL': Per-branch log-likelihood contribution

        Raises
        ------
        FileNotFoundError
            If .uml_rec file doesn't exist
        """
        if self._branch_statistics is None:
            self._parse_uml_rec_file()

        if self._branch_statistics is None:
            raise ValueError("Branch statistics not found in reconciliation file")

        return self._branch_statistics

    def get_file_paths(self) -> Dict[str, Path]:
        """
        Get paths to all ALE output files.

        Returns
        -------
        Dict[str, Path]
            Dictionary mapping file types to their paths
        """
        return {
            'consensus_tree': self.ucons_tree_path,
            'transfers': self.uts_path,
            'reconciliation': self.uml_rec_path
        }

    def files_exist(self) -> Dict[str, bool]:
        """
        Check which ALE output files exist.

        Returns
        -------
        Dict[str, bool]
            Dictionary mapping file types to whether they exist
        """
        return {
            'consensus_tree': self.ucons_tree_path.exists(),
            'transfers': self.uts_path.exists(),
            'reconciliation': self.uml_rec_path.exists()
        }

    def get_all_data(self) -> Dict:
        """
        Parse and return all available data from ALE output files.

        This is a convenience method that extracts everything in one call.

        Returns
        -------
        Dict
            Dictionary containing all parsed data:
            - 'consensus_tree': ete3.Tree or None
            - 'species_tree': ete3.Tree or None
            - 'reconciled_gene_trees': List[ete3.Tree] or None
            - 'transfers': pd.DataFrame or None
            - 'ml_rates': Dict[str, float] or None
            - 'log_likelihood': float or None
            - 'summary_statistics': Dict[str, float] or None
            - 'branch_statistics': pd.DataFrame or None
        """
        data = {}

        # Try to get consensus tree
        try:
            data['consensus_tree'] = self.get_consensus_tree()
        except (FileNotFoundError, ValueError):
            data['consensus_tree'] = None

        # Try to get transfers
        try:
            data['transfers'] = self.get_transfers()
        except (FileNotFoundError, ValueError):
            data['transfers'] = None

        # Try to get reconciliation data
        try:
            data['species_tree'] = self.get_reconciled_tree()
            data['reconciled_gene_trees'] = self.get_reconciled_gene_trees()
            data['ml_rates'] = self.get_ml_rates()
            data['log_likelihood'] = self.get_log_likelihood()
            data['summary_statistics'] = self.get_summary_statistics()
            data['branch_statistics'] = self.get_branch_statistics()
        except (FileNotFoundError, ValueError):
            data['species_tree'] = None
            data['reconciled_gene_trees'] = None
            data['ml_rates'] = None
            data['log_likelihood'] = None
            data['summary_statistics'] = None
            data['branch_statistics'] = None

        return data
