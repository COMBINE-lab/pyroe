from .load_fry import load_fry
from .load_fry import process_output_format
from .pyroe_utils import output_formats


def validate_convert_args(args):
    """
    Perform validation of the arguments for the covert sub-command.
    This will e.g. make sure that the input exists, that the ouput
    format is supported, and that the `output-structure` option makes
    sense as parsed.
    """
    import ast
    import sys

    output_fmt = output_formats()
    if args.output_format not in output_fmt:
        print(f"The output format {args.output_format} was invalid")

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


def convert(args):
    # first make sure that the input is such that
    # the conversion makes sense
    validate_convert_args(args)
    # offload the work of loading the input to `load_fry`
    A = load_fry(args.quant_dir, output_format=args.output_structure, quiet=True)

    # write the output in the requested format
    output_fn = {
        "h5ad": A.write,
        "loom": A.write_loom,
        "csvs": A.write_csvs,
        "zarr": A.write_zarr,
    }

    if args.output_format in output_fn:
        output_fn[args.output_format](args.output)
