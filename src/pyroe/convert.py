import logging
import json

from .load_fry import load_fry
from .load_fry import process_output_format
from .pyroe_utils import output_formats


def validate_convert_args(args):
    """
    Perform validation of the arguments for the covert sub-command.
    This will e.g. make sure that the input exists, that the output
    format is supported, and that the `output-structure` option makes
    sense as parsed.
    """
    import ast
    import sys
    import pathlib

    output_fmt = output_formats()
    if args.output_format not in output_fmt:
        print(f"The output format {args.output_format} was invalid")
        sys.exit(1)

    if args.geneid_to_name is not None:
        p = pathlib.Path(args.geneid_to_name)
        if not (p.is_file() or p.is_fifo()):
            print(f"The path {args.geneid_to_name} doesn't point to a valid file")
            sys.exit(1)

    # if we don't have one of these, then attempt to covert
    # the string to a dictionary.
    out_struct = {}
    if args.output_structure in ["scRNA", "snRNA", "raw", "velocity"]:
        args.ouput_structure = process_output_format(args.output_structure, True)
    else:
        try:
            out_struct = ast.literal_eval(args.output_structure)
            args.output_structure = process_output_format(out_struct, True)
        except Exception:
            print(
                f"Could not parse {args.output_structure} argument to --output-structure"
            )
            sys.exit(1)


def get_id_to_name_map(id_to_name_file):
    d = {}
    with open(id_to_name_file) as ifile:
        for line in ifile:
            toks = line.rstrip().split("\t")
            d[toks[0]] = toks[1]
    return d


def convert_id_to_name(adata, id_to_name):
    unmapped = adata.var_names[adata.var_names.map(id_to_name).isna()].tolist()
    # drop the unmappable names and covnert the rest
    adata = adata[:, ~adata.var_names.map(id_to_name).isna()]
    # make the names unique
    adata.var_names = adata.var_names.map(id_to_name)
    adata.var_names_make_unique()
    return (adata, unmapped)


def convert(args):
    # first make sure that the input is such that
    # the conversion makes sense
    validate_convert_args(args)
    # offload the work of loading the input to `load_fry`
    A = load_fry(args.quant_dir, output_format=args.output_structure, quiet=True)

    id_to_name = args.geneid_to_name
    if id_to_name is not None:
        id_name_map = get_id_to_name_map(id_to_name)
        A, unmapped = convert_id_to_name(A, id_name_map)
        if len(unmapped) > 0:
            logging.info(f"There were {len(unmapped)} gene ids without a mapped name.")
            uout = f"{args.output}_unmapped_ids.json"
            logging.info(f"Writing them to {uout}.")
            udict = {"unmapped_geneids": unmapped}
            with open(uout, "w") as f:
                json.dump(udict, f, indent=4, sort_keys=True)

    # write the output in the requested format
    output_fn = {
        "h5ad": A.write,
        "loom": A.write_loom,
        "csvs": A.write_csvs,
        "zarr": A.write_zarr,
    }

    if args.output_format in output_fn:
        output_fn[args.output_format](args.output)
