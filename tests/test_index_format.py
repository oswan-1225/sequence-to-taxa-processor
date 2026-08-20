"""Tests for the on-disk k-mer index format (save_index / load_index).

The point of the format is not storage, it is refusal. Every failure it
catches is one that otherwise produces a complete, plausible-looking, wrong
result rather than an error:

  - a pre-canonical (v1) index queried with canonical k-mers misses nearly
    every lookup and yields a CSV of nulls
  - an index built at one k and queried at another matches almost nothing

Both used to run to completion and print a success message.
"""
import pickle

import pytest

from classifier_functions import (
    INDEX_FORMAT_VERSION,
    build_kmer_index,
    load_index,
    save_index,
)

GENOMES = {"species_a": "ACGTACGTACGT", "species_b": "TTTTGGGGTTTT"}


def test_save_then_load_round_trips(tmp_path):
    index = build_kmer_index(GENOMES, 4)
    path = str(tmp_path / "index.pkl")

    save_index(path, index, 4)

    assert load_index(path, 4) == index


def test_saved_index_records_its_metadata(tmp_path):
    """The header has to be readable without this module, since diagnosing a
    bad index is exactly the situation where you cannot trust the loader."""
    path = str(tmp_path / "index.pkl")
    save_index(path, build_kmer_index(GENOMES, 4), 4)

    with open(path, "rb") as f:
        payload = pickle.load(f)

    assert payload["format_version"] == INDEX_FORMAT_VERSION
    assert payload["k"] == 4
    assert payload["canonical"] is True
    assert "index" in payload


def test_load_rejects_old_bare_dict_index(tmp_path):
    """A v1 index is a bare {kmer: set(species)} dict with no header.

    This is the case that matters most: two such files are sitting in the
    repo's data/ and results/ folders, and reusing one via --index would
    silently produce garbage.
    """
    path = str(tmp_path / "old.pkl")
    with open(path, "wb") as f:
        pickle.dump(build_kmer_index(GENOMES, 4), f)

    with pytest.raises(ValueError, match="old-format"):
        load_index(path, 4)


def test_old_format_error_names_the_fix(tmp_path):
    """An error a user cannot act on is only marginally better than silence."""
    path = str(tmp_path / "old.pkl")
    with open(path, "wb") as f:
        pickle.dump(build_kmer_index(GENOMES, 4), f)

    with pytest.raises(ValueError) as excinfo:
        load_index(path, 4)

    assert "build_reference.py" in str(excinfo.value)


def test_load_rejects_future_format_version(tmp_path):
    """Refuse to guess at a format written by newer code."""
    path = str(tmp_path / "future.pkl")
    with open(path, "wb") as f:
        pickle.dump({
            "format_version": INDEX_FORMAT_VERSION + 1,
            "k": 4,
            "canonical": True,
            "index": {},
        }, f)

    with pytest.raises(ValueError, match="format v"):
        load_index(path, 4)


def test_load_rejects_k_mismatch(tmp_path):
    """A pre-existing silent failure, closed as a side effect of the header.

    The CLIs take --k independently of --index and nothing checked that they
    agreed, so classifying at k=5 against a k=4 index ran happily and
    reported almost everything as unclassified.
    """
    path = str(tmp_path / "index.pkl")
    save_index(path, build_kmer_index(GENOMES, 4), 4)

    with pytest.raises(ValueError, match="built with k=4"):
        load_index(path, 5)


def test_k_mismatch_error_reports_both_values(tmp_path):
    path = str(tmp_path / "index.pkl")
    save_index(path, build_kmer_index(GENOMES, 4), 4)

    with pytest.raises(ValueError) as excinfo:
        load_index(path, 21)

    message = str(excinfo.value)
    assert "k=4" in message and "k=21" in message
