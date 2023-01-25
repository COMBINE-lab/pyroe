try:
    import scanpy
except ModuleNotFoundError:
    print(
        "scanpy must be installed to enable the load_fry() function. Use `conda install -c scanpy ` or `pip install scanpy` to install it."
    )
    import sys

    sys.exit(1)


def load_fry(frydir, output_format="scRNA", nonzero=False, quiet=False):
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

    nonzero : `bool` (default: `False`)
        True if cells with non-zero expression value across all genes should be filtered in each layer.
        False if unexpressed genes should be kept.

    quiet : `bool` (default: `False`)
        True if function should be quiet.
        False if messages (including error messages) should be printed out.

    Notes
    ----------
    The `output_format` argument takes either a dictionary that defines the customized format or
    a string that represents one of the pre-defined format of the returned `AnnData` object.

    Each of the pre-defined formats contains a `X` field and some optional extra `AnnData.layers`
    obtained from the submatrices representing unspliced (U), spliced (S) and ambiguous (A) counts
    returned by alevin-fry.

    The following formats are defined:

    * "scRNA": \\
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

    """
    import json
    import os
    import pandas as pd

    # since alevin-fry 0.4.1 the generic "meta_info.json"
    # has been replaced by a more informative name for each
    # sub-command. For quantification, it is "quant.json".
    # we check for both files here, in order.
    meta_info_files = ["quant.json", "meta_info.json"]

    fpath = os.path.sep.join([frydir, meta_info_files[0]])
    # first, check for the new file, if we don't find it, check
    # for the old one.
    if not os.path.exists(fpath):
        if not quiet:
            print(
                f"Did not find a {meta_info_files[0]} file, checking for older {meta_info_files[1]}."
            )
        fpath = os.path.sep.join([frydir, meta_info_files[1]])
        # if we don't find the old one either, then return None
        if not os.path.exists(fpath):
            raise IOError(f"Found no {meta_info_files[1]} file either; cannot proceed.")

    # if we got here then we had a valid json file, so
    # use it to get the number of genes, and if we are
    # in USA mode or not.
    meta_info = json.load(open(fpath))
    ng = meta_info["num_genes"]
    usa_mode = meta_info["usa_mode"]
    if not quiet:
        print(f"USA mode: {usa_mode}")

    # if we are in USA mode
    if usa_mode:
        # preparation
        # each gene has 3 splicing statuses, so the actual number of distinct
        # genes is ng/3.
        ng = int(ng / 3)
        output_assays = process_output_format(output_format, quiet)
    elif not quiet:
        print(
            "Processing input in standard mode, the count matrix will be stored in field 'X'."
        )
        if output_format != "scRNA":
            print("Output_format will be ignored.")

    # read the actual input matrix
    af_raw = scanpy.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
    afg = [
        line.rstrip()
        for line in open(
            os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])
        ).readlines()
    ][:ng]
    # read the gene ids
    afg_df = pd.DataFrame(afg, columns=["gene_ids"])
    afg_df = afg_df.set_index("gene_ids")
    # and the barcodes
    abc = [
        line.rstrip()
        for line in open(
            os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])
        ).readlines()
    ]
    abc_df = pd.DataFrame(abc, columns=["barcodes"])
    abc_df.index = abc_df["barcodes"]

    x = af_raw.X
    # if we're not in USA mode, just combine this info into
    # an AnnData object
    if not usa_mode:
        af = scanpy.AnnData(x.T, var=abc_df, obs=afg_df)
        af = af.T

    else:  # USA mode
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

        if not quiet:
            print(f"Filtered {np.sum(~not_zero_genes)} non-expressed genes.")

    return af


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
                if not quiet:
                    print("A undefined Provided output_format string provided.")
                    print("See function help message for details.")
                raise ValueError("Invalid output_format.")
            if not quiet:
                print("Using pre-defined output format:", output_format)
                print(
                    f"Will populate output field X with sum of counts frorm {predefined_format[output_format]['X']}."
                )
                for (k, v) in predefined_format[output_format].items():
                    if k != "X":
                        print(f"Will combine {v} into output layer {k}.")

            return predefined_format[output_format]
        else:
            if not quiet:
                print("Processing user-defined output format.")
            # make sure the X is there
            if "X" not in output_format.keys():
                raise ValueError(
                    'In USA mode some sub-matrices must be assigned to the "X" (default) output.'
                )
            if not quiet:
                print(
                    f"Will populate output field X with sum of counts frorm {output_format['X']}."
                )

            for (k, v) in output_format.items():
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
                if not quiet and (k != "X"):
                    print(f"Will combine {v} into output layer {k}.")

            return output_format
    else:
        raise ValueError(
            "Provided invalid output_format. See function help message for details"
        )
