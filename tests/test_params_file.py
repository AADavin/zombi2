"""``--params``: a TOML parameters file applied as CLI defaults (CLI flags override)."""

import pytest

from zombi2._params import load_params_file
from zombi2.cli import main


# --- the loader ------------------------------------------------------------- #
def test_load_flat_with_hyphen_and_underscore(tmp_path):
    f = tmp_path / "p.toml"
    f.write_text("birth = 1.0\ntips = 12\nmax-family-size = 0.5\n")
    out = load_params_file(str(f), {"birth", "tips", "max_family_size"}, "species")
    assert out == {"birth": 1.0, "tips": 12, "max_family_size": 0.5}  # hyphen -> underscore dest


def test_load_scopes_to_command_section(tmp_path):
    f = tmp_path / "p.toml"
    f.write_text("[species]\nbirth = 2.0\n\n[genomes]\ndup = 0.3\n")
    assert load_params_file(str(f), {"birth"}, "species") == {"birth": 2.0}
    assert load_params_file(str(f), {"dup"}, "genomes") == {"dup": 0.3}


def test_load_unknown_key_raises(tmp_path):
    f = tmp_path / "p.toml"
    f.write_text("birth = 1.0\nnope = 5\n")
    with pytest.raises(ValueError, match="nope"):
        load_params_file(str(f), {"birth"}, "species")


# --- through the CLI -------------------------------------------------------- #
_SPECIES = "birth = 1.0\ndeath = 0.3\ntips = 12\nage = 4.0\nseed = 7\n"


def test_params_file_equals_equivalent_flags(tmp_path):
    f = tmp_path / "s.toml"
    f.write_text(_SPECIES)
    a, b = tmp_path / "A", tmp_path / "B"
    assert main(["species", "--params", str(f), "-o", str(a)]) == 0
    assert main(["species", "--birth", "1.0", "--death", "0.3", "--tips", "12",
                 "--age", "4.0", "--seed", "7", "-o", str(b)]) == 0
    # a scalar `birth = 1.0` is wrapped to match the nargs='+' --birth, and the run is identical
    assert (a / "species_tree.nwk").read_text() == (b / "species_tree.nwk").read_text()


def test_command_line_flag_overrides_file(tmp_path):
    f = tmp_path / "s.toml"
    f.write_text(_SPECIES)
    a, c = tmp_path / "A", tmp_path / "C"
    main(["species", "--params", str(f), "-o", str(a)])
    main(["species", "--params", str(f), "--seed", "999", "-o", str(c)])
    assert (a / "species_tree.nwk").read_text() != (c / "species_tree.nwk").read_text()


def test_pipeline_file_with_sections_and_list_value(tmp_path):
    f = tmp_path / "pipeline.toml"
    f.write_text("[species]\nbirth = 1.0\ndeath = 0.3\ntips = 8\nage = 3.0\nseed = 5\n\n"
                 "[genomes]\ndup = 0.2\nloss = 0.25\norig = 0.5\n"
                 'write = ["profiles", "trees"]\nseed = 5\n')
    sp, g = tmp_path / "S", tmp_path / "G"
    assert main(["species", "--params", str(f), "-o", str(sp)]) == 0
    assert main(["genomes", "--params", str(f), "-t", str(sp / "species_tree.nwk"),
                 "-o", str(g)]) == 0
    assert (g / "Profiles.tsv").exists() and (g / "gene_trees").exists()  # the list write applied


def test_unknown_key_is_a_clean_cli_error(tmp_path):
    f = tmp_path / "s.toml"
    f.write_text("birth = 1.0\nbogus = 5\n")
    with pytest.raises(SystemExit):  # subp.error -> exit 2, not a traceback
        main(["species", "--params", str(f), "-o", str(tmp_path / "X")])


def test_missing_params_file_is_a_clean_cli_error(tmp_path):
    with pytest.raises(SystemExit):
        main(["species", "--params", str(tmp_path / "nope.toml"), "-o", str(tmp_path / "X")])
