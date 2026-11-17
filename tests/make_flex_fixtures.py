"""Write small synthetic alevin-fry output directories for the loader tests.

The files copy the on-disk format of alevin-fry 0.19.0 (``quant`` and
``quant-probe``) and of the older layout with ``meta_info.json``.
All values come from a seeded random generator, so the output is the same
on each run.  The test ``test_make_flex_fixtures.py`` makes sure that the
files in ``tests/data/flex`` are equal to the output of this script.

Usage::

    python tests/make_flex_fixtures.py [OUTDIR]

OUTDIR is ``tests/data/flex`` if you do not give it.

Fixtures (each one is an ``af_quant`` directory, except ``probe_set_v1.csv``):

* ``legacy_usa``: old layout, ``meta_info.json``, USA mode, bare barcodes,
  no ``featureDump.txt``.
* ``flex_gene``: multiplex ``quant-probe`` run without ``--probe-mtx``.
  3 samples with ``_`` in the names, interleaved rows, one empty row,
  ``MeanByMax = NaN``, one all-zero gene, extra ids in
  ``gene_id_to_name.tsv``.
* ``flex_probe``: ``flex_gene`` plus the ``--probe-mtx`` files and the 10x
  probe-set CSV ``probe_set.csv`` (format 3.0), which a future simpleaf
  version will copy into ``af_quant`` (simpleaf 0.31 does not).
  The CSV has 2 probes with ``included=FALSE`` that are not matrix columns.
  The probe columns are not in the CSV order.
  The gene matrix is the probe matrix summed by gene.
* ``flex_no_featuredump``: ``flex_gene`` without ``featureDump.txt``.
* ``flex_bad_featuredump``: ``flex_gene`` with two ``featureDump.txt``
  lines in the wrong order.
* ``flex_single``: single-sample run, bare barcodes, ``featureDump.txt``
  without the ``sample_name`` column.
* ``probe_set_v1.csv``: the same probes in an older 10x probe-set CSV
  format, without the ``gene_name`` and ``region`` columns.
"""

import hashlib
import json
import os
import shutil
import sys

import numpy as np

SEED = 20261117
BC_LEN = 16

# gene id, gene symbol, number of probes
GENES = [
    ("GENE000010", "Alpha", 1),
    ("GENE000020", "Beta", 2),
    ("GENE000030", "Gamma", 1),
    ("GENE000040", "Delta", 1),
    ("GENE000050", "Epsilon", 1),
    ("GENE000060", "Zeta", 1),
]
# the last gene has no counts (for the `nonzero` tests)
ZERO_GENE = "GENE000060"
# ids that the index knows but that are not matrix columns
EXTRA_GENES = [("GENE000015", "AlphaBis"), ("GENE000070", "Eta")]
# 10x probe region of the included probes
REGIONS = ["spliced", "spliced", "unspliced", "spliced", "unspliced", "spliced", "spliced"]
# probes with included=FALSE: gene id, gene symbol, region
EXCLUDED_PROBES = [
    ("GENE000015", "AlphaBis", "unspliced"),
    ("GENE000030", "Gamma", "spliced"),
]
PROBE_SET_HEADER = [
    "#panel_name=Synthetic Test Probe Set",
    "#panel_type=predesigned",
    "#reference_genome=TEST1",
    "#reference_version=2026-A",
]
# the probe column order (RAD header order) as indices into the t2g order
PROBE_COL_ORDER = [5, 2, 0, 6, 1, 4, 3]

# sample names: "S_1" is a prefix of "S_1_A", and "1e5" looks like a number
SAMPLES = ["S_1", "S_1_A", "1e5"]
CELLS_PER_SAMPLE = [4, 3, 3]
# row order of the matrix: the samples are interleaved
ROW_SAMPLE_ORDER = [0, 1, 0, 2, 1, 0, 2, 1, 2, 0]
# 0-based row with no counts
EMPTY_ROW = 6

TG_MAP_PLACEHOLDER = "/path/to/probe_index/probe_t2g.tsv"
QUANT_CMD = "alevin-fry quant-probe -i af_quant -o af_quant -t 4 -m probe_t2g.tsv -r cr-like --use-mtx"


def fmt(v):
    """Print a number as Rust `Display` prints an f32 (`1`, `0.5`, `NaN`)."""
    if np.isnan(v):
        return "NaN"
    return np.format_float_positional(np.float32(v), trim="-")


def barcode(rng):
    return "".join(rng.choice(list("ACGT"), size=BC_LEN))


def probe_id(gid, sym, i):
    h = hashlib.sha1(f"{gid}{i}".encode()).hexdigest()[:7]
    return f"{gid}|{sym}|{h}"


def probe_ids():
    """Return (probe_id, gene_id, region) of the included probes, in CSV order."""
    out = []
    k = 0
    for gid, sym, n in GENES:
        for i in range(n):
            out.append((probe_id(gid, sym, i), gid, REGIONS[k]))
            k += 1
    return out


def probe_set_csv(fmt):
    """Return the lines of a 10x probe-set CSV, sorted by gene id as 10x does.

    `fmt` "3.0" has the gene_name and region columns; "1.0" does not.
    """
    symbol = {g: s for g, s, _ in GENES + [(g, s, 0) for g, s in EXTRA_GENES]}
    rows = [(p, g, "TRUE", r) for p, g, r in probe_ids()]
    rows += [(probe_id(g, s, 9), g, "FALSE", r) for g, s, r in EXCLUDED_PROBES]
    lines = [f"#probe_set_file_format={fmt}"] + PROBE_SET_HEADER
    if fmt == "3.0":
        lines.append("gene_id,probe_seq,probe_id,included,region,gene_name")
    else:
        lines.append("gene_id,probe_seq,probe_id,included")
    for p, g, inc, r in sorted(rows, key=lambda x: (x[1], x[0])):
        digest = hashlib.sha512(p.encode()).digest()
        seq = "".join("ACGT"[b % 4] for b in digest[:50])
        line = f"{g},{seq},{p},{inc}"
        lines.append(line + (f",{r},{symbol[g]}" if fmt == "3.0" else ""))
    return lines


def write_text(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=4)
        f.write("\n")


def write_mtx(path, mat, rng):
    """Write `mat` (rows = cells) as the sprs MatrixMarket writer does.

    The size line pads nnz to 20 characters, and the entries are not sorted.
    """
    rows, cols = np.nonzero(mat)
    order = rng.permutation(len(rows))
    lines = [
        "%%MatrixMarket matrix coordinate real general",
        "% written by sprs",
        f"{mat.shape[0]} {mat.shape[1]} {len(rows):<20}",
    ]
    for i in order:
        lines.append(f"{rows[i] + 1} {cols[i] + 1} {fmt(mat[rows[i], cols[i]])}")
    write_text(path, lines)


def feature_dump(labels, cbs, samples, gene_mat, rng):
    """Return the featureDump.txt lines; `samples` is None for one sample."""
    head = ["CB"] + (["sample_name"] if samples is not None else [])
    head += [
        "CorrectedReads",
        "MappedReads",
        "DeduplicatedReads",
        "MappingRate",
        "DedupRate",
        "MeanByMax",
        "NumGenesExpressed",
        "NumGenesOverMean",
    ]
    lines = ["\t".join(head)]
    for i in range(len(labels)):
        row = gene_mat[i]
        umis = row.sum()
        mapped = int(umis) + int(rng.integers(1, 20))
        corrected = mapped + int(rng.integers(0, 5))
        expressed = row[row > 0]
        if len(expressed) > 0:
            mean = expressed.mean()
            mean_by_max = mean / expressed.max()
            over_mean = int((expressed > mean).sum())
        else:
            mean_by_max = float("nan")
            over_mean = 0
        vals = [cbs[i]] + ([samples[i]] if samples is not None else [])
        vals += [
            str(corrected),
            str(mapped),
            fmt(umis),
            fmt(mapped / corrected),
            fmt(umis / mapped),
            fmt(mean_by_max),
            str(len(expressed)),
            str(over_mean),
        ]
        lines.append("\t".join(vals))
    return lines


def quant_json(num_cells, num_genes, feature_axis, num_probes, probe_nnz):
    probe_level = feature_axis == "probe"
    probe_mtx = probe_nnz is not None
    return {
        "alt_resolved_cell_numbers": [],
        "cmd": QUANT_CMD + (" --probe-mtx" if probe_mtx else ""),
        "dump_eq": False,
        # chunk numbers in the collated RAD, not matrix rows
        "empty_resolved_cell_numbers": [7],
        "feature_axis": feature_axis,
        "molecule_table": None,
        "num_genes": num_genes,
        "num_probes": num_probes,
        "num_quantified_cells": num_cells,
        "num_same_gene_multimapping_reads": 1 if probe_level else None,
        "num_tiny_cell_resolved": 0,
        "probe_mtx": probe_mtx,
        "probe_nnz": probe_nnz,
        "quant_options": {
            "cmdline": QUANT_CMD + (" --probe-mtx" if probe_mtx else ""),
            "crlike_umi_edit": 0,
            "dump_eq": False,
            "dump_molecules": False,
            "filter_list": None,
            "init_uniform": False,
            "input_dir": "af_quant",
            "large_graph_thresh": 0,
            "num_bootstraps": 0,
            "num_threads": 4,
            "output_dir": "af_quant",
            "probe_level": probe_level,
            "probe_mtx": probe_mtx,
            "pug_exact_umi": False,
            "resolution": "CellRangerLike",
            "sa_model": "WinnerTakeAll",
            "small_thresh": 100,
            "summary_stat": False,
            "tg_map": TG_MAP_PLACEHOLDER,
            "version": "0.19.0",
        },
        "resolution_strategy": "CellRangerLike",
        "tiny_cell_resolved_cell_numbers": [],
        "usa_mode": False,
        "version_str": "0.19.0",
    }


def sample_info(row_samples):
    samples = []
    for k, name in enumerate(SAMPLES):
        n = int(sum(s == name for s in row_samples))
        samples.append(
            {
                "barcode": f"0x{k + 1:05x}",
                "collatable_reads": 100 * n,
                "name": name,
                "num_cells": n,
                "num_reads": 100 * n,
            }
        )
    return {
        "num_barcodes": 2,
        "num_samples": len(SAMPLES),
        "samples": samples,
        "total_cells": len(row_samples),
    }


def gene_id_to_name_lines():
    pairs = [(g, s) for g, s, _ in GENES] + EXTRA_GENES
    return [f"{g}\t{s}" for g, s in sorted(pairs)]


def write_flex_gene(d, labels, cbs, samples, gene_mat, num_probes, rng):
    write_text(os.path.join(d, "alevin", "quants_mat_rows.txt"), labels)
    write_text(os.path.join(d, "alevin", "quants_mat_cols.txt"), [g for g, _, _ in GENES])
    write_mtx(os.path.join(d, "alevin", "quants_mat.mtx"), gene_mat, rng)
    write_text(os.path.join(d, "gene_id_to_name.tsv"), gene_id_to_name_lines())
    write_json(os.path.join(d, "sample_info.json"), sample_info(samples))
    write_text(
        os.path.join(d, "featureDump.txt"),
        feature_dump(labels, cbs, samples, gene_mat, rng),
    )
    write_json(
        os.path.join(d, "quant.json"),
        quant_json(len(labels), len(GENES), "probe", num_probes, None),
    )


def make_flex(outdir, rng):
    probes = probe_ids()
    n_cells = len(ROW_SAMPLE_ORDER)

    # one barcode per (sample, cell); the first cell of S_1 and S_1_A
    # share the same bare barcode, as real multiplex data can
    sample_cbs = [[barcode(rng) for _ in range(n)] for n in CELLS_PER_SAMPLE]
    sample_cbs[1][0] = sample_cbs[0][0]
    next_cell = [0] * len(SAMPLES)
    samples, cbs = [], []
    for s in ROW_SAMPLE_ORDER:
        samples.append(SAMPLES[s])
        cbs.append(sample_cbs[s][next_cell[s]])
        next_cell[s] += 1
    labels = [f"{s}_{cb}" for s, cb in zip(samples, cbs)]

    # probe matrix in t2g order, then reorder the columns to the RAD order
    probe_t2g_mat = rng.integers(0, 4, size=(n_cells, len(probes)))
    probe_t2g_mat[rng.random(probe_t2g_mat.shape) < 0.4] = 0
    probe_t2g_mat[EMPTY_ROW, :] = 0
    for j, (_, gid, _) in enumerate(probes):
        if gid == ZERO_GENE:
            probe_t2g_mat[:, j] = 0
    probe_mat = probe_t2g_mat[:, PROBE_COL_ORDER]
    probe_cols = [probes[j][0] for j in PROBE_COL_ORDER]

    # gene matrix = probe matrix summed through the t2g
    gene_index = {g: i for i, (g, _, _) in enumerate(GENES)}
    p2g = np.zeros((len(probes), len(GENES)), dtype=int)
    for j, (_, gid, _) in enumerate(probes):
        p2g[j, gene_index[gid]] = 1
    gene_mat = probe_t2g_mat @ p2g

    write_text(os.path.join(outdir, "probe_set_v1.csv"), probe_set_csv("1.0"))

    gene_dir = os.path.join(outdir, "flex_gene")
    write_flex_gene(gene_dir, labels, cbs, samples, gene_mat, len(probes), rng)

    # the other fixtures are copies of flex_gene with one change each
    d = os.path.join(outdir, "flex_probe")
    shutil.copytree(gene_dir, d, dirs_exist_ok=True)
    write_mtx(os.path.join(d, "alevin", "probe_quants_mat.mtx"), probe_mat, rng)
    shutil.copyfile(
        os.path.join(d, "alevin", "quants_mat_rows.txt"),
        os.path.join(d, "alevin", "probe_quants_mat_rows.txt"),
    )
    write_text(os.path.join(d, "alevin", "probe_quants_mat_cols.txt"), probe_cols)
    write_text(os.path.join(d, "probe_set.csv"), probe_set_csv("3.0"))
    write_json(
        os.path.join(d, "quant.json"),
        quant_json(
            n_cells, len(GENES), "probe", len(probes), int(np.count_nonzero(probe_mat))
        ),
    )

    d = os.path.join(outdir, "flex_no_featuredump")
    shutil.copytree(gene_dir, d, dirs_exist_ok=True)
    os.remove(os.path.join(d, "featureDump.txt"))

    # swap two featureDump lines of different samples
    d = os.path.join(outdir, "flex_bad_featuredump")
    shutil.copytree(gene_dir, d, dirs_exist_ok=True)
    fd_path = os.path.join(d, "featureDump.txt")
    with open(fd_path) as f:
        lines = f.read().splitlines()
    lines[1], lines[2] = lines[2], lines[1]
    write_text(fd_path, lines)


def make_flex_single(outdir, rng):
    n_cells = 4
    cbs = [barcode(rng) for _ in range(n_cells)]
    gene_mat = rng.integers(0, 4, size=(n_cells, len(GENES)))
    d = os.path.join(outdir, "flex_single")
    write_text(os.path.join(d, "alevin", "quants_mat_rows.txt"), cbs)
    write_text(os.path.join(d, "alevin", "quants_mat_cols.txt"), [g for g, _, _ in GENES])
    write_mtx(os.path.join(d, "alevin", "quants_mat.mtx"), gene_mat, rng)
    write_text(os.path.join(d, "gene_id_to_name.tsv"), gene_id_to_name_lines())
    write_text(
        os.path.join(d, "featureDump.txt"),
        feature_dump(cbs, cbs, None, gene_mat, rng),
    )
    qj = quant_json(n_cells, len(GENES), "gene", None, None)
    qj["cmd"] = qj["quant_options"]["cmdline"] = qj["cmd"].replace(
        "quant-probe", "quant"
    )
    qj["quant_options"]["probe_level"] = False
    write_json(os.path.join(d, "quant.json"), qj)


def make_legacy_usa(outdir, rng):
    n_cells = 5
    ng = len(GENES)
    cbs = [barcode(rng) for _ in range(n_cells)]
    # column blocks: S (gene id), U (gene-U), A (gene-A)
    mat = rng.integers(0, 4, size=(n_cells, 3 * ng))
    mat[rng.random(mat.shape) < 0.4] = 0
    ids = [g for g, _, _ in GENES]
    d = os.path.join(outdir, "legacy_usa")
    write_text(os.path.join(d, "alevin", "quants_mat_rows.txt"), cbs)
    write_text(
        os.path.join(d, "alevin", "quants_mat_cols.txt"),
        ids + [f"{g}-U" for g in ids] + [f"{g}-A" for g in ids],
    )
    write_mtx(os.path.join(d, "alevin", "quants_mat.mtx"), mat, rng)
    write_text(os.path.join(d, "gene_id_to_name.tsv"), gene_id_to_name_lines())
    cmd = "alevin-fry quant -i af_quant -o af_quant -t 4 -m t2g_3col.tsv -r cr-like --use-mtx"
    write_json(
        os.path.join(d, "meta_info.json"),
        {
            "alt_resolved_cell_numbers": [],
            "cmd": cmd,
            "dump_eq": False,
            "empty_resolved_cell_numbers": [],
            "num_genes": 3 * ng,
            "num_quantified_cells": n_cells,
            "resolution_strategy": "CellRangerLike",
            "usa_mode": True,
            "version_str": "0.4.0",
        },
    )


def make_all(outdir):
    rng = np.random.default_rng(SEED)
    make_legacy_usa(outdir, rng)
    make_flex(outdir, rng)
    make_flex_single(outdir, rng)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "data", "flex")
    make_all(out)
