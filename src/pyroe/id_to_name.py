import pyranges
import logging


def id_to_name(params):
    a = pyranges.read_gtf(params.gtf_file)
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
