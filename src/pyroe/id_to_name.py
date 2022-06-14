import pyranges
import logging
import pathlib
import sys


def id_to_name(params):
    format_map = {
        ".gtf": pyranges.read_gtf,
        ".gff": pyranges.read_gff3,
        ".gff3": pyranges.read_gff3,
    }
    annot_reader = None
    if params.format is None:
        p = pathlib.Path(params.gtf_file)
        suffs = [z.lower() for z in p.suffixes]

        z = None
        if len(suffs) == 1:
            z = suffs[0]
        elif len(suffs) == 2 and suffs[-1] == ".gz":
            z = suffs[-2]

        if z in format_map.keys():
            annot_reader = format_map[z]
        else:
            logging.error(
                "Could not determine format of annotation file. Please provide it explicitly."
            )
            sys.exit(1)
    else:
        fmt = params.format.lower()
        if fmt not in ["gtf", "gff3"]:
            logging.error(
                f'Format must be either "gtf" or "gff3", but {fmt} was provided.'
            )
            sys.exit(1)
        annot_reader = format_map[fmt]

    a = annot_reader(params.gtf_file)
    # only bother looking at `gene` features
    a_genes = a[(a.Feature == "gene")]
    id_name = {}
    for k, df in a_genes:
        cdf = df[["gene_id", "gene_name"]].to_dict(orient="records")
        id_name.update({d["gene_id"]: d["gene_name"] for d in cdf})

    logging.info(f"generated mappings for {len(id_name)} gene ids.")
    logging.info(f"writing output to {params.output}")

    with open(params.output, "w") as ofile:
        for k, v in id_name.items():
            ofile.write(f"{k}\t{v}\n")
