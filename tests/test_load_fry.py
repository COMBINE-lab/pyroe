"""Offline tests for `pyroe.load_fry` on the fixtures in tests/data/flex."""

import json
import os
import shutil

import anndata
import numpy as np
import pytest
import scipy.io

from pyroe import load_fry
from pyroe.load_fry import PROBE_VAR_COLUMNS

FLEX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "flex")
LEGACY = os.path.join(FLEX_DIR, "legacy_usa")
FLEX_GENE = os.path.join(FLEX_DIR, "flex_gene")
FLEX_SINGLE = os.path.join(FLEX_DIR, "flex_single")
FLEX_PROBE = os.path.join(FLEX_DIR, "flex_probe")
PROBE_SET = os.path.join(FLEX_PROBE, "probe_set.csv")
PROBE_SET_V1 = os.path.join(FLEX_DIR, "probe_set_v1.csv")


def _copy_fixture(tmp_path, name):
    dst = str(tmp_path / name)
    shutil.copytree(os.path.join(FLEX_DIR, name), dst)
    return dst


def _read_lines(path):
    with open(path) as f:
        return f.read().splitlines()


def _write_lines(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _gene_names(d):
    return dict(line.split("\t") for line in _read_lines(os.path.join(d, "gene_id_to_name.tsv")))


def _csv_rows(path):
    """The data rows of a probe-set CSV as dicts, keyed by probe_id."""
    lines = [line for line in _read_lines(path) if not line.startswith("#")]
    head = lines[0].split(",")
    return {r["probe_id"]: r for r in (dict(zip(head, x.split(","))) for x in lines[1:])}


@pytest.fixture(scope="module")
def flex_gene():
    return load_fry(FLEX_GENE, quiet=True)


@pytest.fixture(scope="module")
def flex_probe():
    """flex_probe loaded with the probe_set.csv that load_fry finds in frydir."""
    return load_fry(FLEX_PROBE, feature_level="probe", quiet=True)


class TestOldOutputs:
    def test_legacy_usa_has_no_new_fields(self):
        adata = load_fry(LEGACY, quiet=True)
        assert adata.obs.columns.tolist() == ["barcodes"]
        assert adata.var.columns.tolist() == ["gene_ids", "gene_names"]
        assert list(adata.uns.keys()) == []

    def test_legacy_usa_scrna_counts(self):
        adata = load_fry(LEGACY, quiet=True)
        m = scipy.io.mmread(os.path.join(LEGACY, "alevin", "quants_mat.mtx")).toarray()
        s, u, a = np.split(m, 3, axis=1)
        assert np.array_equal(adata.X.toarray(), s + a)
        assert np.array_equal(adata.layers["unspliced"].toarray(), u)

    def test_single_sample_featuredump_adds_no_columns(self):
        adata = load_fry(FLEX_SINGLE, quiet=True)
        assert adata.obs.columns.tolist() == ["barcodes"]

    def test_single_sample_featuredump_is_not_checked(self, tmp_path):
        d = _copy_fixture(tmp_path, "flex_single")
        path = os.path.join(d, "featureDump.txt")
        lines = _read_lines(path)
        _write_lines(path, [lines[0], lines[2], lines[1]] + lines[3:-1])
        adata = load_fry(d, quiet=True)
        assert adata.obs.columns.tolist() == ["barcodes"]


class TestFlexGene:
    def test_obs_names_are_row_labels(self, flex_gene):
        rows = _read_lines(os.path.join(FLEX_GENE, "alevin", "quants_mat_rows.txt"))
        assert flex_gene.obs_names.tolist() == rows
        assert flex_gene.obs["barcodes"].tolist() == rows

    def test_obs_gets_only_cell_barcode_and_sample_name(self, flex_gene):
        assert flex_gene.obs.columns.tolist() == ["barcodes", "cell_barcode", "sample_name"]

    def test_sample_name_and_cell_barcode_rebuild_row_label(self, flex_gene):
        obs = flex_gene.obs
        rebuilt = (obs["sample_name"] + "_" + obs["cell_barcode"]).tolist()
        assert rebuilt == flex_gene.obs_names.tolist()

    def test_sample_names_keep_underscores_and_number_like_names(self, flex_gene):
        assert sorted(flex_gene.obs["sample_name"].unique()) == ["1e5", "S_1", "S_1_A"]

    def test_cells_per_sample_match_sample_info(self, flex_gene):
        with open(os.path.join(FLEX_GENE, "sample_info.json")) as f:
            info = json.load(f)
        expected = {s["name"]: s["num_cells"] for s in info["samples"]}
        assert flex_gene.obs["sample_name"].value_counts().to_dict() == expected

    def test_same_bare_barcode_in_two_samples_keeps_rows_apart(self, flex_gene):
        assert flex_gene.obs["cell_barcode"].duplicated(keep=False).sum() == 2
        assert flex_gene.obs_names.is_unique

    def test_x_equals_mtx(self, flex_gene):
        m = scipy.io.mmread(os.path.join(FLEX_GENE, "alevin", "quants_mat.mtx"))
        assert np.array_equal(flex_gene.X.toarray(), m.toarray())

    def test_var_is_left_join_of_gene_names(self, flex_gene):
        cols = _read_lines(os.path.join(FLEX_GENE, "alevin", "quants_mat_cols.txt"))
        names = _gene_names(FLEX_GENE)
        assert len(names) > len(cols)
        assert flex_gene.var_names.tolist() == cols
        assert flex_gene.var["gene_ids"].tolist() == cols
        assert flex_gene.var["gene_names"].tolist() == [names[g] for g in cols]

    def test_nonzero_removes_zero_count_gene(self):
        adata = load_fry(FLEX_GENE, nonzero=True, quiet=True)
        assert adata.n_vars == 5
        assert "GENE000060" not in adata.var_names


class TestSampleSplit:
    def test_bad_featuredump_row_order_raises(self):
        with pytest.raises(ValueError, match="same row order"):
            load_fry(os.path.join(FLEX_DIR, "flex_bad_featuredump"), quiet=True)

    def test_featuredump_missing_line_raises(self, tmp_path):
        d = _copy_fixture(tmp_path, "flex_gene")
        path = os.path.join(d, "featureDump.txt")
        _write_lines(path, _read_lines(path)[:-1])
        with pytest.raises(ValueError, match="rows"):
            load_fry(d, quiet=True)

    def test_without_featuredump_split_on_last_underscore(self, flex_gene):
        adata = load_fry(os.path.join(FLEX_DIR, "flex_no_featuredump"), quiet=True)
        assert adata.obs.equals(flex_gene.obs)

    def test_without_featuredump_bad_label_raises(self, tmp_path):
        d = _copy_fixture(tmp_path, "flex_no_featuredump")
        path = os.path.join(d, "alevin", "quants_mat_rows.txt")
        lines = _read_lines(path)
        lines[0] = lines[0].replace("_", "-")
        _write_lines(path, lines)
        with pytest.raises(ValueError, match="Cannot split"):
            load_fry(d, quiet=True)

    def test_without_featuredump_or_multiplex_files_no_split(self, tmp_path):
        d = _copy_fixture(tmp_path, "flex_no_featuredump")
        os.remove(os.path.join(d, "sample_info.json"))
        adata = load_fry(d, quiet=True)
        assert adata.obs.columns.tolist() == ["barcodes"]


class TestFlexProbe:
    def test_var_names_are_probe_columns_in_matrix_order(self, flex_probe):
        cols = _read_lines(os.path.join(FLEX_PROBE, "alevin", "probe_quants_mat_cols.txt"))
        assert flex_probe.var_names.tolist() == cols
        assert flex_probe.var["probe_id"].tolist() == cols

    def test_var_columns_come_from_probe_set_csv(self, flex_probe):
        rows = _csv_rows(PROBE_SET)
        var = flex_probe.var
        assert var.columns.tolist() == PROBE_VAR_COLUMNS
        assert var["gene_ids"].tolist() == [rows[p]["gene_id"] for p in var.index]
        assert var["gene_names"].tolist() == [rows[p]["gene_name"] for p in var.index]
        assert var["region"].tolist() == [rows[p]["region"] for p in var.index]

    def test_excluded_probes_are_not_in_var(self, flex_probe):
        rows = _csv_rows(PROBE_SET)
        excluded = [p for p, r in rows.items() if r["included"] == "FALSE"]
        assert len(excluded) == 2
        assert not set(excluded) & set(flex_probe.var_names)

    def test_explicit_argument_is_used(self, tmp_path, flex_probe):
        d = _copy_fixture(tmp_path, "flex_probe")
        os.remove(os.path.join(d, "probe_set.csv"))
        adata = load_fry(d, feature_level="probe", probe_ref=PROBE_SET, quiet=True)
        assert adata.var.equals(flex_probe.var)

    def test_old_format_csv_takes_gene_names_from_gene_id_to_name(self, flex_probe):
        adata = load_fry(FLEX_PROBE, feature_level="probe", probe_ref=PROBE_SET_V1, quiet=True)
        assert adata.var.columns.tolist() == ["probe_id", "gene_ids", "gene_names"]
        assert adata.var["gene_names"].equals(flex_probe.var["gene_names"])

    def test_empty_gene_name_and_region_are_nan(self, tmp_path):
        ref = str(tmp_path / "probe_set.csv")
        lines = _read_lines(PROBE_SET)
        i = next(k for k, line in enumerate(lines) if line.startswith("GENE000010,"))
        lines[i] = ",".join(lines[i].split(",")[:4] + ["", ""])
        _write_lines(ref, lines)
        adata = load_fry(FLEX_PROBE, feature_level="probe", probe_ref=ref, quiet=True)
        row = adata.var[adata.var["gene_ids"] == "GENE000010"]
        assert row["gene_names"].isna().all() and row["region"].isna().all()

    def test_no_csv_found_raises_with_path(self, tmp_path):
        d = _copy_fixture(tmp_path, "flex_probe")
        os.remove(os.path.join(d, "probe_set.csv"))
        with pytest.raises(FileNotFoundError, match="probe-set CSV") as err:
            load_fry(d, feature_level="probe", quiet=True)
        assert os.path.join(d, "probe_set.csv") in str(err.value)

    def test_missing_argument_file_raises(self, tmp_path):
        path = str(tmp_path / "missing.csv")
        with pytest.raises(FileNotFoundError, match="missing.csv"):
            load_fry(FLEX_PROBE, feature_level="probe", probe_ref=path, quiet=True)

    def test_matrix_probe_not_in_csv_raises(self, tmp_path):
        ref = str(tmp_path / "probe_set.csv")
        _write_lines(ref, [x for x in _read_lines(PROBE_SET) if not x.startswith("GENE000040,")])
        with pytest.raises(ValueError, match="not in"):
            load_fry(FLEX_PROBE, feature_level="probe", probe_ref=ref, quiet=True)

    def test_duplicate_probe_in_csv_raises(self, tmp_path):
        ref = str(tmp_path / "probe_set.csv")
        lines = _read_lines(PROBE_SET)
        _write_lines(ref, lines + lines[-1:])
        with pytest.raises(ValueError):
            load_fry(FLEX_PROBE, feature_level="probe", probe_ref=ref, quiet=True)

    def test_unknown_region_value_raises(self, tmp_path):
        ref = str(tmp_path / "probe_set.csv")
        _write_lines(ref, [x.replace(",unspliced,", ",U,") for x in _read_lines(PROBE_SET)])
        with pytest.raises(ValueError, match="spliced"):
            load_fry(FLEX_PROBE, feature_level="probe", probe_ref=ref, quiet=True)

    def test_csv_without_gene_id_raises(self, tmp_path):
        ref = str(tmp_path / "probe_set.csv")
        _write_lines(ref, [x.replace("gene_id,", "gene,", 1) for x in _read_lines(PROBE_SET)])
        with pytest.raises(ValueError, match="gene_id"):
            load_fry(FLEX_PROBE, feature_level="probe", probe_ref=ref, quiet=True)

    def test_probe_sum_by_gene_equals_gene_matrix(self, flex_probe):
        gene = load_fry(FLEX_PROBE, quiet=True)
        genes = gene.var_names.tolist()
        p2g = np.zeros((flex_probe.n_vars, gene.n_vars))
        for j, g in enumerate(flex_probe.var["gene_ids"]):
            p2g[j, genes.index(g)] = 1
        assert np.array_equal(flex_probe.X.toarray() @ p2g, gene.X.toarray())

    def test_obs_equal_gene_obs(self, flex_probe, flex_gene):
        assert flex_probe.obs.equals(flex_gene.obs)

    def test_x_equals_probe_mtx(self, flex_probe):
        m = scipy.io.mmread(os.path.join(FLEX_PROBE, "alevin", "probe_quants_mat.mtx"))
        assert np.array_equal(flex_probe.X.toarray(), m.toarray())

    def test_nonzero_removes_zero_count_probe(self):
        adata = load_fry(FLEX_PROBE, feature_level="probe", nonzero=True, quiet=True)
        assert adata.n_vars == 6
        assert "GENE000060" not in adata.var["gene_ids"].tolist()

    @pytest.mark.parametrize("name", ["flex_gene", "legacy_usa"])
    def test_probe_level_without_probe_mtx_raises(self, name):
        with pytest.raises(ValueError, match="probe_mtx"):
            load_fry(
                os.path.join(FLEX_DIR, name),
                feature_level="probe",
                probe_ref=PROBE_SET,
                quiet=True,
            )

    def test_unknown_feature_level_raises(self):
        with pytest.raises(ValueError, match="feature_level"):
            load_fry(FLEX_PROBE, feature_level="transcript", quiet=True)


class TestH5adRoundTrip:
    @pytest.mark.parametrize("level", ["gene", "probe"])
    def test_flex_anndata_writes_and_reads_h5ad(self, tmp_path, level):
        adata = load_fry(FLEX_PROBE, feature_level=level, quiet=True)
        path = str(tmp_path / "flex.h5ad")
        adata.write_h5ad(path)
        back = anndata.read_h5ad(path)
        assert np.array_equal(back.X.toarray(), adata.X.toarray())
        assert back.obs.equals(adata.obs)
        assert back.var.equals(adata.var)
