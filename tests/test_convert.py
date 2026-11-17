"""Offline tests for `pyroe.convert` and the `pyroe convert` command."""

import argparse
import os
import subprocess
import sys

import anndata
import numpy as np
import pytest
import scipy.io

from pyroe.convert import validate_convert_args
from pyroe.load_fry import PREDEFINED_FORMATS

HERE = os.path.dirname(os.path.abspath(__file__))
LEGACY = os.path.join(HERE, "data", "flex", "legacy_usa")
PYROE_BIN = os.path.join(os.path.dirname(HERE), "bin", "pyroe")


def _args(output_structure):
    return argparse.Namespace(
        output_structure=output_structure, output_format="h5ad", geneid_to_name=None
    )


def _legacy_blocks():
    m = scipy.io.mmread(os.path.join(LEGACY, "alevin", "quants_mat.mtx")).toarray()
    s, u, a = np.split(m, 3, axis=1)
    return s, u, a


class TestValidateConvertArgs:
    @pytest.mark.parametrize("name", ["scRNA", "S+A", "snRNA", "all", "U+S+A", "raw"])
    def test_preset_name_becomes_format_dict(self, name):
        args = _args(name)
        validate_convert_args(args)
        assert args.output_structure == PREDEFINED_FORMATS[name.lower()]

    def test_dict_string_becomes_format_dict(self):
        args = _args("{'X': ['U'], 'sa': ['S', 'A']}")
        validate_convert_args(args)
        assert args.output_structure == {"X": ["U"], "sa": ["S", "A"]}

    def test_bad_structure_exits(self):
        with pytest.raises(SystemExit):
            validate_convert_args(_args("not a format"))


class TestConvertCommand:
    def _run(self, tmp_path, *extra):
        out = str(tmp_path / "out.h5ad")
        subprocess.run(
            [sys.executable, PYROE_BIN, "convert", LEGACY, out, *extra],
            check=True,
            capture_output=True,
        )
        return anndata.read_h5ad(out)

    def test_default_output_structure_is_scrna(self, tmp_path):
        adata = self._run(tmp_path)
        s, u, a = _legacy_blocks()
        assert np.array_equal(adata.X.toarray(), s + a)
        assert np.array_equal(adata.layers["unspliced"].toarray(), u)

    def test_sum_preset_has_no_extra_layer(self, tmp_path):
        adata = self._run(tmp_path, "--output-structure", "U+S+A")
        s, u, a = _legacy_blocks()
        assert np.array_equal(adata.X.toarray(), s + u + a)
        assert list(adata.layers.keys()) == []
