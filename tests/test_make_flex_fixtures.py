"""Tests for the fixture generator `make_flex_fixtures.py`."""

import filecmp
import os

import numpy as np
import pandas as pd
import scipy.io

import make_flex_fixtures

FLEX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "flex")


def _all_files(root):
    out = []
    for d, _, files in os.walk(root):
        out += [os.path.relpath(os.path.join(d, f), root) for f in files]
    return sorted(out)


def test_committed_fixtures_equal_generator_output(tmp_path):
    make_flex_fixtures.make_all(str(tmp_path))
    committed = _all_files(FLEX_DIR)
    assert committed == _all_files(str(tmp_path))
    _, mismatch, errors = filecmp.cmpfiles(
        FLEX_DIR, str(tmp_path), committed, shallow=False
    )
    assert mismatch == [] and errors == [], f"changed files: {mismatch + errors}"


def test_fixtures_are_small():
    size = sum(os.path.getsize(os.path.join(FLEX_DIR, f)) for f in _all_files(FLEX_DIR))
    assert size < 50 * 1024


def _probe_set(path):
    return pd.read_csv(path, comment="#", dtype=str)


def test_flex_probe_matrix_summed_by_gene_equals_gene_matrix():
    d = os.path.join(FLEX_DIR, "flex_probe", "alevin")
    gene = scipy.io.mmread(os.path.join(d, "quants_mat.mtx")).toarray()
    probe = scipy.io.mmread(os.path.join(d, "probe_quants_mat.mtx")).toarray()
    genes = pd.read_table(os.path.join(d, "quants_mat_cols.txt"), header=None)[0]
    probes = pd.read_table(os.path.join(d, "probe_quants_mat_cols.txt"), header=None)[0]
    csv = _probe_set(os.path.join(FLEX_DIR, "flex_probe", "probe_set.csv"))
    p2g = csv.set_index("probe_id")["gene_id"]
    summed = np.zeros_like(gene)
    for j, p in enumerate(probes):
        summed[:, genes.tolist().index(p2g[p])] += probe[:, j]
    assert np.array_equal(summed, gene)


def test_probe_columns_are_the_included_csv_probes_in_another_order():
    d = os.path.join(FLEX_DIR, "flex_probe")
    probes = pd.read_table(
        os.path.join(d, "alevin", "probe_quants_mat_cols.txt"), header=None
    )[0].tolist()
    csv = _probe_set(os.path.join(d, "probe_set.csv"))
    included = csv["probe_id"][csv["included"] == "TRUE"].tolist()
    assert sorted(probes) == sorted(included)
    assert probes != included
    assert (csv["included"] == "FALSE").sum() == 2


def test_old_format_probe_set_has_same_probes_without_names():
    new = _probe_set(os.path.join(FLEX_DIR, "flex_probe", "probe_set.csv"))
    old = _probe_set(os.path.join(FLEX_DIR, "probe_set_v1.csv"))
    assert old.columns.tolist() == ["gene_id", "probe_seq", "probe_id", "included"]
    assert old.equals(new[old.columns.tolist()])


def test_flex_bad_featuredump_differs_from_rows_file():
    d = os.path.join(FLEX_DIR, "flex_bad_featuredump")
    fd = pd.read_table(os.path.join(d, "featureDump.txt"), dtype=str)
    rows = pd.read_table(
        os.path.join(d, "alevin", "quants_mat_rows.txt"), header=None, dtype=str
    )[0]
    labels = fd["sample_name"] + "_" + fd["CB"]
    assert sorted(labels) == sorted(rows)
    assert (labels != rows).sum() == 2
