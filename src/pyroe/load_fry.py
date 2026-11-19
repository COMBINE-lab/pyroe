from .pyroe_utils import say


def load_fry(
    frydir,
    output_format="scRNA",
    aux_columns=["X", "Y"],
    gene_id_to_name=None,
    nonzero=False,
    quiet=False,
    feature_level="gene",
    probe_ref=None,
):
    """
    load alevin-fry quantification result into an AnnData object

    Required Parameters
    ----------
    frydir : `str`
        The path to a output directory returned by alevin-fry quant command. \\
        The directory containing the alevin-fry quantification (i.e. the the quant.json file & alevin subdirectory).

    Optional Parameters
    ----------
    output_format : `str` or `dict`
        A string represents one of the pre-defined output formats, which are "scRNA", "S+A", "snRNA", "all", "U+S+A" and "velocity". \\
        If a customized format of the returned `AnnData` is needed, one can pass a Dictionary.\\
        See Notes section for details.

    aux_columns : `list[str]` (default: `["X", "Y"]`)
        A list of strings contains the column names of the auxiliary information in the barcodes file starting from the second column. \\
        The first column is assumed to be the barcodes and is named as "barcodes". \\
        Extra auxiliary columns in the barcodes file without a specified name will be ignored.

    gene_id_to_name : `str` or `None` (default: `None`)
        The path to a file that contains the mapping from gene names to gene ids. \\
        It is only needed if \\
            1. you are not using the simpleaf pipeline (`simpleaf index` + `simpleaf quant`), \\
            2. you have such a file, and,
            3. you want to add this information to the coldata of your anndata.
        If you do, please ensure it is a tab-separated, two-column file without a header, and the first column is the gene ids and the second column is the gene names.

    nonzero : `bool` (default: `False`)
        True if cells with non-zero expression value across all genes should be filtered in each layer.
        False if unexpressed genes should be kept.

    quiet : `bool` (default: `False`)
        True if function should be quiet.
        False if messages (including error messages) should be printed out.

    feature_level : `str` (default: `"gene"`)
        "gene" loads the cell x gene matrix `alevin/quants_mat.mtx`.
        "probe" loads the cell x probe matrix `alevin/probe_quants_mat.mtx`, which
        `alevin-fry quant-probe --probe-mtx` (`simpleaf multiplex-quant --probe-mtx`)
        writes. The function raises a `ValueError` if `quant.json` does not have
        `"probe_mtx": true`.

    probe_ref : `str` or `None` (default: `None`)
        Only for `feature_level="probe"`. The 10x probe-set CSV of the run, which
        maps each probe to its gene. If it is `None`, the function uses
        `probe_set.csv` in `frydir`. The function raises a `FileNotFoundError`
        if the file does not exist.
        The function skips the lines at the top that start with `#`. The CSV must
        have a `probe_id` column. The probes are joined to the CSV by `probe_id`,
        because the probe column order of the matrix is not the order of the CSV.
        The function raises a `ValueError` if a matrix probe is not in the CSV.

    Notes
    ----------
    The `output_format` argument takes either a dictionary that defines the customized format or
    a string that represents one of the pre-defined format of the returned `AnnData` object.

    Each of the pre-defined formats contains a `X` field and some optional extra `AnnData.layers`
    obtained from the submatrices representing unspliced (U), spliced (S) and ambiguous (A) counts
    returned by alevin-fry.

    The following formats are defined:

    * "scRNA" and "S+A": \\
        This format is recommended for single cell RNA-sequencing experiments.
        It returns a `X` field that contains the S+A count of each gene in each cell,
        and a `unspliced` field that contains the U count of each gene.

    * "snRNA", "all" and "U+S+A": \\
        These three formats are the same. They return a `X` field that contains the U+S+A
        count of each gene in each cell without any extra layers.
        It is recommended for single-nucleus RNA-sequencing experiments.
        CellRanger 7 returns this format for both single-cell and single-nucleus experiments.

    * "raw": \\
        This format uses the S count matrix as the `X` field and put the U, S, and A counts into three
        separate layers, which are "unspliced", "spliced" and "ambiguous".

    * "velocity": \\
        This format is the same as "scRNA", except it contains two extra layers: the "spliced" layer,
        which contains the S+A counts, and the "unspliced" layer, which contains the U counts.

    A custom output format can be defined using a Dictionary specifying the desired format of the output `Anndata` object.
    If the input is not a USA mode quantification directory, this parameter is ignored
    and the count matrix is returned in the `X` field of the returned `AnnData` object.  If the input
    quantification directory contains a USA mode quantification, then there are 3 sub-matrices that can
    be referenced in the dictionary; 'U', 'S', 'A' containing, respectively, unspliced, spliced and
    ambiguous counts.  The dictionary should have entries of the form `key` (str) : `value` (list[str]).
    The following constraints apply : there should be one key-value pair with the key `X`, the resulting
    counts will be returned in the `X` field of the AnnData object. There can be an arbitrary number
    of other key-value pairs, but each will be returned as a layer of the resulting AnnData object.
    Within the key-value pairs, the key refers to the layer name that will be given to the combined
    count matrix upon output, and the value should be a subset of `['U', 'S', 'A']` that defines
    which sub-matrices should be summed.  For example:
    `{'X' : ['S', 'A'], 'unspliced' : ['U']}`
    will result in a return AnnData object where the X field has a matrix in which each entry
    corresponds to the summed spliced and ambiguous counts for each gene in each cell, and there
    is an additional "unspliced" layer, whose counts are taken directly from the unspliced sub-matrix.

    Returns:
    ----------
        An AnnData object with X and layers corresponding to the requested `output_format`.

        * For a multiplex (Flex) run, the row label (`obs_names`, `obs["barcodes"]`) is
          `<sample_name>_<cell_barcode>`, and `obs["sample_name"]` and
          `obs["cell_barcode"]` give its two parts. They come from the `sample_name`
          and `CB` columns of `featureDump.txt`, whose line i describes row i; the
          function raises a `ValueError` if a line does not match its row. If
          `featureDump.txt` does not exist, but `collation_manifest.bin` or
          `sample_info.json` shows a multiplex run, the function gives a warning and
          does not add the two columns. Other outputs get no new columns.
        * With `feature_level="probe"`, `var_names` are the probe ids in
          `probe_quants_mat_cols.txt`, and `var` is the probe-set CSV row of each
          probe: all CSV columns, with their CSV names (for example `gene_id`,
          `probe_seq`, `probe_id`, `included`, `region`, `gene_name`), as
          `pandas.read_csv` reads them with its default settings: for example
          `included` is a bool, and an empty value is `NaN`.

    """
    import json
    import os
    import pandas as pd
    import scanpy

    # since alevin-fry 0.4.1 the generic "meta_info.json"
    # has been replaced by a more informative name for each
    # sub-command. For quantification, it is "quant.json".
    # we check for both files here, in order.
    meta_info_files = ["quant.json", "meta_info.json"]

    fpath = os.path.sep.join([frydir, meta_info_files[0]])
    # first, check for the new file, if we don't find it, check
    # for the old one.
    if not os.path.exists(fpath):
        say(
            quiet,
            f"Did not find a {meta_info_files[0]} file, checking for older {meta_info_files[1]}.",
        )

        fpath = os.path.sep.join([frydir, meta_info_files[1]])
        # if we don't find the old one either, then return None
        if not os.path.exists(fpath):
            raise IOError(
                "The profvided `frydir` doesn't contain required meta info file; cannot proceed."
            )

    # if we got here then we had a valid json file, so
    # use it to get the number of genes, and if we are
    # in USA mode or not.
    meta_info = json.load(open(fpath))
    ng = meta_info["num_genes"]
    usa_mode = meta_info["usa_mode"]

    if feature_level not in ["gene", "probe"]:
        raise ValueError(
            f"feature_level must be 'gene' or 'probe', not '{feature_level}'."
        )
    probe_level = feature_level == "probe"
    if probe_level and meta_info.get("probe_mtx") is not True:
        raise ValueError(
            f"feature_level='probe' needs a probe matrix, but {fpath} does not have "
            "'probe_mtx': true. Run alevin-fry quant-probe (or simpleaf "
            "multiplex-quant) with --probe-mtx."
        )
    mat_prefix = "probe_quants_mat" if probe_level else "quants_mat"

    say(quiet, f"USA mode: {usa_mode}")

    # if we are in USA mode
    if usa_mode:
        # preparation
        # each gene has 3 splicing statuses, so the actual number of distinct
        # genes is ng/3.
        ng = int(ng / 3)
        output_assays = process_output_format(output_format, quiet)
    else:
        say(
            quiet,
            "Processing input in standard mode, the count matrix will be stored in field 'X'.",
        )
        if output_format != "scRNA":
            say(quiet, "Output_format will be ignored.")

    # if we have a gene name to id mapping, use it
    # Otherwise, we see if there is a default mapping file
    # We first hope the user gives one
    gene_id_to_name_path = gene_id_to_name
    # make sure the file exists
    if gene_id_to_name_path is not None:
        if not os.path.exists(gene_id_to_name_path):
            say(
                quiet,
                f"The provided gene id to name mapping file {gene_id_to_name_path} does not exist; ignored.",
            )
            gene_id_to_name_path = None

    # we then check if we can find a default mapping file
    if gene_id_to_name_path is None:
        default_gene_id_to_name_path = os.path.sep.join([frydir, "gene_id_to_name.tsv"])
        if os.path.exists(default_gene_id_to_name_path):
            gene_id_to_name_path = default_gene_id_to_name_path
            say(
                quiet,
                f"Using simpleaf gene id to name mapping file: {gene_id_to_name_path}",
            )

    if probe_level:
        afg_df = _read_probe_var(frydir, _find_probe_ref(frydir, probe_ref))
    else:
        # read the gene ids
        afg_df = pd.read_table(
            os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"]),
            names=["gene_ids"],
            nrows=ng
        )

    # read the file if we find it
    if gene_id_to_name_path is not None and not probe_level:
        gene_id_to_name_df = pd.read_table(
            gene_id_to_name_path, names=["gene_ids", "gene_names"]
        )
        afg_df = pd.merge(afg_df, gene_id_to_name_df, on="gene_ids", how="left")

    afg_df = afg_df.set_index("probe_id" if probe_level else "gene_ids", drop=False)

    # read the barcodes file
    abc_df = pd.read_table(
        os.path.sep.join([frydir, "alevin", f"{mat_prefix}_rows.txt"]), header=None
    )

    # if column names and num columns don't match, match them
    # we have only barcode
    if abc_df.shape[1] == 1:
        columns = ["barcodes"]
    else:
        say(quiet, f"Found {abc_df.shape[1] - 1} auxiliary columns in barcodes file.")
        columns = ["barcodes"] + aux_columns

    if abc_df.shape[1] != len(columns):
        ncol = min(abc_df.shape[1], len(columns))

        say(
            quiet,
            f"Number of auxiliary columns in barcodes file does not match the provided column names. Using the first {ncol} columns in the barcode file with names: {columns[:ncol]}.",
        )

        abc_df.columns = columns[:ncol]
        columns = columns[:ncol]

    abc_df = abc_df.set_axis(columns, axis=1)
    abc_df = abc_df.set_index(columns[0], drop=False)
    abc_df = _add_sample_info(frydir, abc_df, quiet)

    say(quiet, "Reading the quantification matrix.")
    x = _read_mtx(os.path.sep.join([frydir, "alevin", f"{mat_prefix}.mtx"]))
    # if we're not in USA mode, just combine this info into
    # an AnnData object
    if not usa_mode:
        say(quiet, "Constructing the AnnData object.")
        af = scanpy.AnnData(x.T, var=abc_df, obs=afg_df)
        af = af.T

    else:  # USA mode
        say(quiet, "Converting USA mode quantification into specified layers.")
        # otherwise, combine the sub-matrices into the output object as
        # specified by `output_assays`
        rd = {"S": range(0, ng), "U": range(ng, 2 * ng), "A": range(2 * ng, 3 * ng)}
        xcounts = output_assays["X"]
        o = x[:, rd[xcounts[0]]]
        for wc in xcounts[1:]:
            o += x[:, rd[wc]]
        af = scanpy.AnnData(o.T, var=abc_df, obs=afg_df)
        af = af.T

        # now, if there are other layers requested, populate those
        for other_layer in output_assays.keys() - "X":
            xcounts = output_assays[other_layer]
            o = x[:, rd[xcounts[0]]]
            for wc in xcounts[1:]:
                o += x[:, rd[wc]]
            af.layers[other_layer] = o

    if nonzero:
        import numpy as np

        not_zero_genes = af.X.sum(axis=0).A1 > 0
        if usa_mode:
            for other_layer in output_assays.keys() - "X":
                not_zero_genes = np.logical_or(
                    not_zero_genes, af.layers[other_layer].sum(axis=0).A1 > 0
                )

        af = af[:, not_zero_genes]

        say(quiet, f"Filtered {np.sum(~not_zero_genes)} non-expressed genes.")

    say(quiet, "Done.")
    return af


def _read_mtx(path, block_bytes=8 << 20):
    """
    Read a Matrix Market file into a `scipy.sparse.csr_matrix`.

    The result is the same as `scanpy.read_mtx(path).X`: float32 data,
    int32 indices and indptr (int64 if the matrix is too large for int32),
    sorted columns in each row. A coordinate file with two entries for the
    same row and column is refused (alevin-fry does not write such files).
    An array (dense) file is read with `scipy.io.mmread`, as before.

    The function does not keep the full COO matrix in memory. It reads the
    file two times, one block of whole lines at a time, and parses each block
    with `scipy.io.mmread`. The first pass counts the entries in each row.
    The second pass writes each entry into its row of the result, with a
    cursor for each row. Thus the peak memory is approximately the size of
    the result plus one block, for any order of the entries in the file.
    """
    import io
    import itertools
    import os
    import numpy as np
    import scipy
    from packaging.version import Version
    from scipy.io import mmread
    from scipy.sparse import csr_matrix

    class Parts:
        # a stream that gives one part for each read() call, without a copy
        def __init__(self, parts):
            self.parts = [p for p in reversed(parts) if p]

        def read(self, size=-1):
            return self.parts.pop() if self.parts else b""

    # scipy >= 1.12 reads with fast_matrix_market; the older mmread is
    # written in Python
    fmm = Version(scipy.__version__).release >= (1, 12)  # also 1.12.0rc1

    def stream(parts):
        # fast_matrix_market uses only read(); the older mmread also needs
        # readline(), so join the parts
        return Parts(parts) if fmm else io.BytesIO(b"".join(parts))

    with open(path, "rb") as f:
        banner = f.readline()
        if b" coordinate " not in banner:  # a dense array: read it as before
            return csr_matrix(mmread(path).astype(np.float32))
        size, line = f.readline(), 2  # line: the line number of size
        while size.startswith(b"%") or not size.strip():
            if not size:
                raise ValueError(f"{path} has no size line (rows, columns, entries).")
            size, line = f.readline(), line + 1
        body = f.tell()
    nrows, ncols, nlines = (int(v) for v in size.split())
    # at least 8 blocks, so that a small file also needs less memory than
    # the full COO matrix
    block_bytes = min(block_bytes, max(1 << 16, os.path.getsize(path) // 8))

    # the body lines that mmread skips, as today: lines of only these bytes,
    # and with the older mmread also lines that start with "%"
    skip = b" \t\r" if fmm else b" \t\r\v\f"
    blank = np.zeros(256, dtype=bool)
    blank[list(skip)] = True

    def entry_lines(text, tail, exact):
        # the number of line ends in text + tail, and the number of lines
        # that are not blank (exact=True) or of all lines (exact=False);
        # the Matrix Market format allows blank lines after the first line
        a = np.frombuffer(text, np.uint8)
        is_end = a == 10
        ends = int(np.count_nonzero(is_end))
        n = ends
        if exact:  # go back from each line end over the blank bytes
            end = np.flatnonzero(is_end)
            first = np.r_[0, end[:-1] + 1]  # the first byte of each line
            pos = end - 1
            todo = np.flatnonzero(pos >= first)
            n -= end.size - todo.size  # empty lines
            while todo.size:
                todo = todo[blank[a[pos[todo]]]]
                pos[todo] -= 1
                done = pos[todo] < first[todo]
                n -= int(np.count_nonzero(done))
                todo = todo[~done]
            if not fmm:
                n -= int(np.count_nonzero(a[first] == 37))  # "%"
        last = text[text.rfind(b"\n") + 1:] + tail  # the line that tail ends
        data = bool(last.strip(skip + b"\n")) and (fmm or last[:1] != b"%")
        return n + data, ends + tail.endswith(b"\n")

    def parse(n, text, tail, first_line):
        if not n:
            return None
        head = banner + b"%d %d %d\n" % (nrows, ncols, n)
        try:
            return mmread(stream([head, text, tail]))
        except ValueError as e:
            raise ValueError(
                f"{path}: cannot read the entry lines from file line {first_line} "
                f"on (line 3 in the message below is file line {first_line}): {e}"
            ) from e

    counts = []  # (entry lines, line ends) of each block, from pass 1

    def blocks(first_pass):
        # yield the (row, col, value) arrays of each block of whole lines
        first_line = line + 1
        with open(path, "rb") as f:
            f.seek(body)
            for i in itertools.count():
                text = f.read(block_bytes)
                if not text:
                    break
                tail = f.readline()  # the rest of the last line
                if not first_pass:
                    m = parse(counts[i][0], text, tail, first_line)
                else:
                    n, ends = entry_lines(text, tail, exact=False)
                    try:
                        m = parse(n, text, tail, first_line)
                    except ValueError:
                        # blank lines make the first count too large
                        k = entry_lines(text, tail, exact=True)[0]
                        if k == n:
                            raise
                        n = k
                        m = parse(n, text, tail, first_line)
                    counts.append((n, ends))
                first_line += counts[i][1]
                if m is not None:
                    yield m.row, m.col, m.data

    def runs(row):
        # the start, row and length of each run of equal row values; start
        # and length are int32, because a block has less than 2**31 entries
        start = np.flatnonzero(row[1:] != row[:-1]).astype(np.int32)
        start += 1
        start = np.concatenate((np.zeros(1, np.int32), start))
        return start, row[start], np.diff(start, append=np.int32(row.size))

    # pass 1: the number of entries in each row
    row_count = np.zeros(nrows, dtype=np.int64)
    for r, _, _ in blocks(True):
        _, rr, ll = runs(r)
        np.add.at(row_count, rr, ll)
    if sum(n for n, _ in counts) != nlines:
        raise ValueError(
            f"{path} has {sum(n for n, _ in counts)} entry lines, "
            f"but its header says {nlines}."
        )
    nnz = int(row_count.sum())
    idx_dtype = np.int32 if max(nnz, ncols) < 2**31 else np.int64
    indptr = np.zeros(nrows + 1, dtype=np.int64)
    np.cumsum(row_count, out=indptr[1:])
    del row_count

    # pass 2: write each entry into its row; cursor is the next free place
    # of each row
    indices = np.empty(nnz, dtype=idx_dtype)
    data = np.empty(nnz, dtype=np.float32)
    cursor = indptr[:-1].copy()
    for r, col, val in blocks(False):
        start, rr, ll = runs(r)
        dest = cursor[rr]
        # a row with more than one run in this block: put its runs one after
        # the other, in file order (stable sort of the runs by row)
        o = np.argsort(rr, kind="stable")
        sorted_len = ll[o]
        before = np.cumsum(sorted_len)  # in sorted order: the lengths of
        before -= sorted_len  # the earlier runs of the block
        sorted_row = rr[o]
        new_row = np.ones(o.size, dtype=bool)
        np.not_equal(sorted_row[1:], sorted_row[:-1], out=new_row[1:])
        before -= np.maximum.accumulate(np.where(new_row, before, 0))
        dest[o] += before  # the lengths of the earlier runs of the same row
        del o, sorted_len, before, sorted_row, new_row
        np.add.at(cursor, rr, ll)
        dest -= start
        pos = np.repeat(dest, ll)
        del start, rr, ll, dest
        pos += np.arange(pos.size)
        indices[pos] = col
        data[pos] = val
        del pos

    x = csr_matrix(
        (data, indices, indptr.astype(idx_dtype)), shape=(nrows, ncols), copy=False
    )
    x.sum_duplicates()  # sort the rows that came in more than one run, in place
    if x.nnz != nnz:
        raise ValueError(f"{path} has more than one entry for a row and column.")
    return x


def _add_sample_info(frydir, abc_df, quiet):
    """
    For a multiplex run, add `cell_barcode` and `sample_name` to the row table
    `abc_df`. The row label of a multiplex run is `<sample_name>_<CB>`.

    If featureDump.txt has a `sample_name` column, use its `CB` and
    `sample_name` columns: line i describes row i. If featureDump.txt does not
    exist in a multiplex run, warn and add nothing. Other outputs do not change.
    """
    import os
    import warnings
    import numpy as np
    import pandas as pd

    labels = abc_df["barcodes"].astype(str).to_numpy()
    fd_path = os.path.join(frydir, "featureDump.txt")

    if os.path.exists(fd_path):
        with open(fd_path) as f:
            if "sample_name" not in f.readline().rstrip("\r\n").split("\t"):
                return abc_df
        # keep_default_na=False: a sample name such as "NA" stays a string
        fd = pd.read_table(
            fd_path, usecols=["CB", "sample_name"], dtype=str, keep_default_na=False
        )
        if fd.shape[0] != len(labels):
            raise ValueError(
                f"{fd_path} has {fd.shape[0]} cells, but the count matrix has {len(labels)} rows."
            )
        expected = (fd["sample_name"] + "_" + fd["CB"]).to_numpy()
        bad = np.flatnonzero(expected != labels)
        if bad.size > 0:
            raise ValueError(
                f"{fd_path} does not have the same row order as the rows file: "
                f"{bad.size} rows are different. The first is row {bad[0] + 1}: "
                f"'{expected[bad[0]]}' in featureDump.txt, '{labels[bad[0]]}' in the rows file."
            )
        abc_df["cell_barcode"] = fd["CB"].to_numpy()
        abc_df["sample_name"] = fd["sample_name"].to_numpy()
        return abc_df

    multiplex_files = ["collation_manifest.bin", "sample_info.json"]
    if any(os.path.exists(os.path.join(frydir, f)) for f in multiplex_files):
        warnings.warn(
            f"{fd_path} is missing in a multiplex run, so obs['sample_name'] and "
            "obs['cell_barcode'] are not added."
        )

    return abc_df


def _find_probe_ref(frydir, probe_ref):
    """Return the path of the 10x probe-set CSV (see `load_fry`)."""
    import os

    path = probe_ref if probe_ref is not None else os.path.join(frydir, "probe_set.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Cannot find the probe-set CSV: {path} does not exist. Give the 10x "
            "probe-set CSV of the run with the probe_ref argument."
        )
    return path


def _read_probe_var(frydir, probe_ref_path):
    """
    Return the rows of the 10x probe-set CSV for the probes of the probe
    matrix, in the column order of `probe_quants_mat.mtx`. All CSV columns are
    kept, as pandas reads them with its default settings.
    """
    import os
    import pandas as pd

    with open(probe_ref_path) as f:
        n_meta = 0
        for line in f:
            if not line.startswith("#"):
                break
            n_meta += 1
    csv = pd.read_csv(probe_ref_path, skiprows=n_meta)

    probes = pd.read_table(
        os.path.join(frydir, "alevin", "probe_quants_mat_cols.txt"),
        names=["probe_id"],
        dtype=str,
        keep_default_na=False,
    )
    # join by name: the matrix column order is not the order of the CSV
    var = pd.merge(
        probes, csv, on="probe_id", how="left", validate="one_to_one", indicator=True
    )
    missing = var["probe_id"][var["_merge"] == "left_only"]
    if len(missing) > 0:
        raise ValueError(
            f"{len(missing)} probes of the probe matrix are not in {probe_ref_path}, "
            f"for example '{missing.iloc[0]}'."
        )
    return var.drop(columns="_merge")[csv.columns]


def process_output_format(output_format, quiet):
    # make sure output_format isn't empty
    if not output_format:
        raise ValueError("output_format cannot be empty")

    if isinstance(output_format, (str, dict)):
        if isinstance(output_format, str):
            predefined_format = {
                "scrna": {"X": ["S", "A"], "unspliced": ["U"]},
                "S+A": {"X": ["S", "A"]},
                "snrna": {"X": ["U", "S", "A"]},
                "all": {"X": ["U", "S", "A"]},
                "U+S+A": {"X": ["U", "S", "A"]},
                "velocity": {
                    "X": ["S", "A"],
                    "spliced": ["S", "A"],
                    "unspliced": ["U"],
                },
                "raw": {
                    "X": ["S"],
                    "spliced": ["S"],
                    "unspliced": ["U"],
                    "ambiguous": ["A"],
                },
            }

            output_format = output_format.lower()
            if output_format not in predefined_format.keys():
                # invalid output_format string
                say(quiet, "A undefined Provided output_format string provided.")
                say(quiet, "See function help message for details.")
                raise ValueError("Invalid output_format.")
            say(quiet, f"Using pre-defined output format: {output_format}")
            say(
                quiet,
                f"Will populate output field X with sum of counts from {predefined_format[output_format]['X']}.",
            )
            for k, v in predefined_format[output_format].items():
                if k != "X":
                    say(quiet, f"Will combine {v} into output layer {k}.")

            return predefined_format[output_format]
        else:
            say(quiet, "Processing user-defined output format.")
            # make sure the X is there
            if "X" not in output_format.keys():
                raise ValueError(
                    'In USA mode some sub-matrices must be assigned to the "X" (default) output.'
                )
            say(
                quiet,
                f"Will populate output field X with sum of counts frorm {output_format['X']}.",
            )

            for k, v in output_format.items():
                if not v:
                    # empty list
                    raise ValueError(
                        f"The element list of key '{k}' in output_format is empty. Please remove it."
                    )

                # v contains Non-USA element
                if len(set(v) - set(["U", "S", "A"])) != 0:
                    # invalid value
                    raise ValueError(
                        f"Found non-USA element in output_format element list '{v}' for key '{k}'; cannot proceed."
                    )
                if k != "X":
                    say(quiet, f"Will combine {v} into output layer {k}.")

            return output_format
    else:
        raise ValueError(
            "Provided invalid output_format. See function help message for details"
        )
