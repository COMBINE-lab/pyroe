import pyranges
import logging
import pathlib
import sys


def id_to_name(params):
    format_map = {
        "gtf": pyranges.read_gtf,
        "gff": pyranges.read_gff3,
        "gff3": pyranges.read_gff3,
    }
    annot_reader = None
    if params.format is None:
        p = pathlib.Path(params.gtf_file)
        suffs = [z.lower().strip('.') for z in p.suffixes]

        z = None
        if len(suffs) >= 1:
            # look at the final suffix
            z = suffs[-1]
            # if the final suffix is gz and there are
            # suffixes preceding it, check the penultimate
            # one and use that
            if z == "gz" and len(suffs) > 1:
                z = suffs[-2]

        if z is None or z not in format_map.keys():
            logging.error(
                "Could not determine format of annotation file. Please provide it explicitly."
            )
            sys.exit(1)
        else:
            annot_reader = format_map[z]
    else:
        fmt = params.format.lower()
        if fmt not in format_map.keys():
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
