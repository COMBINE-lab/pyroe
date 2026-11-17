"""Offline tests for `pyroe.load_fry` on the fixtures in tests/data/flex."""

import os

import numpy as np
import pytest
import scipy.io

from pyroe import load_fry
from pyroe.load_fry import PREDEFINED_FORMATS, process_output_format

FLEX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "flex")
LEGACY = os.path.join(FLEX_DIR, "legacy_usa")


@pytest.fixture(scope="module")
def legacy_blocks():
    """The S, U and A blocks of the legacy USA matrix, as dense arrays."""
    m = scipy.io.mmread(os.path.join(LEGACY, "alevin", "quants_mat.mtx")).toarray()
    s, u, a = np.split(m, 3, axis=1)
    return {"S": s, "U": u, "A": a}


def _sum_blocks(blocks, names):
    return sum(blocks[n] for n in names)


class TestOutputFormatPresets:
    @pytest.mark.parametrize("name", ["S+A", "U+S+A", "s+a", "u+s+a"])
    def test_process_output_format_accepts_sum_presets(self, name):
        assert process_output_format(name, True) == PREDEFINED_FORMATS[name.lower()]

    @pytest.mark.parametrize("name", ["scRNA", "SCRNA", "snRNA", "All", "Velocity", "RAW"])
    def test_process_output_format_ignores_case(self, name):
        assert process_output_format(name, True) == PREDEFINED_FORMATS[name.lower()]

    def test_process_output_format_rejects_unknown_preset(self):
        with pytest.raises(ValueError, match="Invalid output_format"):
            process_output_format("S+U", True)

    def test_predefined_format_keys_are_lower_case(self):
        assert all(k == k.lower() for k in PREDEFINED_FORMATS)

    @pytest.mark.parametrize("name", sorted(PREDEFINED_FORMATS))
    def test_load_fry_preset_layers_are_block_sums(self, name, legacy_blocks):
        adata = load_fry(LEGACY, output_format=name, quiet=True)
        spec = PREDEFINED_FORMATS[name]
        assert np.array_equal(adata.X.toarray(), _sum_blocks(legacy_blocks, spec["X"]))
        assert set(adata.layers.keys()) == set(spec) - {"X"}
        for layer, names in spec.items():
            if layer != "X":
                assert np.array_equal(
                    adata.layers[layer].toarray(), _sum_blocks(legacy_blocks, names)
                ), layer
