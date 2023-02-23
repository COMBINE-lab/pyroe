import os
import warnings
import subprocess
import shutil
import pyranges as pr
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from packaging.version import parse as parse_version
import logging

bed_required_fields = [
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Name",
    "Gene",
    "splice_status",
]


def append_extra(extra_infile, out_fa, out_t2g3col, id2name_path, col_status):
    """
    extra_infile : str
        path to the extra FASTA format input file from which to read the extra records
    out_fa : str
        path to the output FASTA file being created (will be *appended* to)
    out_t2g3col : str
        path to the output t2g 3-column file (will be *appended* to)
    id2name_path : str
        path to the output (existing) gene id to gene name tsv (will be *appended* to)
    col_status : char
        splicing status to use for these records (one of 'S' or 'U')
    """

    # write extra sequences to the t2g file
    with open(extra_infile) as extra_in_fa:
        with open(out_fa, "a") as splici_fa:
            with open(out_t2g3col, "a") as t2g:
                with open(id2name_path, "a") as id_to_name:
                    for title, sequence in SimpleFastaParser(extra_in_fa):
                        tid = title.split()[0]
                        splici_fa.write(f">{tid}\n")
                        splici_fa.write(f"{sequence}\n")
                        t2g.write(f"{tid}\t{tid}\t{col_status}\n")
                        id_to_name.write(f"{tid}\t{tid}\n")


def dedup_sequences(output_dir, in_fa, out_fa):
    """
    This function deduplicates the fasta sequences in `in_fa`, writes the
    deduplicated sequences in fasta format to `out_fa`, and also produces a
    2-column tsv file in `ouput_dir`/duplicate_entries.tsv recording the
    entries that have been removed and the retained entry that they match.
    Note: It is allowed (and in some cases intended) that `in_fa` == `out_fa`.
    In this case, the input fasta file will be overwritten.

    output_dir : str
        The output directory where duplicate_entries.tsv will be written
    in_fa : str
        The path to the input fasta file to be deduplicated
    out_fa : str
        The path to where the ouput deduplicated fasta file should be written
    """
    record_representatives = {}

    # read from the input fasta file and track
    # the duplicate sequences
    for record in SeqIO.parse(in_fa, "fasta"):
        seq = str(record.seq)
        name = record.id.split()[0]
        # if we haven't seen it yet, this is the representative
        # otherwise append the name of the duplicate
        if seq in record_representatives:
            record_representatives[seq].append(name)
        else:
            record_representatives[seq] = [name]

    # write out the duplicate entries in the same format
    # the salmon / pufferfish indexer uses
    dup_file_name = os.path.sep.join([output_dir, "duplicate_entries.tsv"])
    with open(dup_file_name, "w") as dup_file:
        dup_file.write("RetainedRef\tDuplicateRef\n")
        for k, vu in record_representatives.items():
            if len(vu) > 1:
                v = sorted(vu)
                for dup in v[1:]:
                    dup_file.write(f"{v[0]}\t{dup}\n")

    # write the deduplicated output to file
    # Note: it might be that out_file == in_file and we
    # are overwriting the input.
    with open(out_fa, "w") as ofile:
        for k, vu in record_representatives.items():
            v = sorted(vu)
            rec = SeqRecord(Seq(k), id=v[0], description="")
            SeqIO.write(rec, ofile, "fasta")


def check_gr(gr, output_dir):
    """
    This function checks the validity of a PyRanges object and generates a clean GTF file in the expected format if there is any invalid record in the input PyRanges object.
    It applies the following rules:
    1. Each non-gene record has to have a valid transcript_id. If this is not satisfied, it returns an error. Only the records with a valid transcript_id will be written to the clean_gtf.gtf.
    2. For gene_id and gene_name metadata field,
            - If these two fields are entirely missing in the GTF file, An error will be returned. At the same time, in the clean_gtf.gtf, the two fields will be imputed using the transcript_id fields.
            - If one of these two fields is completely missing, a warning will be generated, and the missing field will be imputed using the other one.
            - if some records have missing gene_id and/or gene_name, a warning will be printed, and the missing values will be imputed by the following rules: For records miss gene_id or gene_name, impute the missing one using the other one; If both are missing, impute them using transcript_id, which cannot be missing.
    3. If there is no "transcript" or "gene" feature record, a warning will be printed. Moreover, those missing records will be imputed using the "exon" feature records: The Start and End site of the gene/transcript will be imputed as the bounds of their corresponding exons.
    4. If the boundaries defined in the transcripts'/genes' feature records do not match those implied by their exons' feature records, report a warning but still use transcripts'/genes' feature records to extract unspliced sequences. To be specific, if some but not all transcripts/genes have their corresponding transcripts'/genes' feature records, or the Start and/or End site defined in the transcript/gene feature records do not match the corresponding exons' bounds, then the existing transcripts'/genes' feature records will be used to extract unspliced transcripts. At the same time, in the clean_gtf.gtf, all genes/transcripts that appeared in the exon feature records will have their corresponding transcripts'/genes' feature records, in which the boundaries match the corresponding exons' bounds.

    Args:
        gr (`PyRanges`): A stranded PyRanges object

    Return:
        `PyRanges`: A PyRanges with missing values imputed if possible.

        Besides, a clean_gtf.gtf file will be generated in the output directory
        if there is any invalid records in the input GTF file.
    """

    # split gene type records with others
    # we don't use gene records in splici construction
    clean_gtf_path = os.path.join(output_dir, "clean_gtf.gtf")
    clean_gr = pr.PyRanges()

    # If required fields are missing, quit
    if "transcript_id" not in gr.columns:
        logging.critical(
            " The input GTF file doesn't contain transcript_id metadata field; Cannot proceed."
        )

    if "Feature" not in gr.columns:
        logging.critical(
            " The input GTF file doesn't contain feature field; Cannot proceed."
        )

    if "gene_id" not in gr.columns:
        # use gene_name as gene_id if exists, return an error otherwise
        if "gene_name" not in gr.columns:
            logging.error(
                " The input GTF file doesn't contain gene_id and gene_name metadata field; Cannot proceed."
            )
        else:
            logging.warning(
                " The gene_id field does not exist; Imputing using gene_name."
            )
            gene_id = pd.Series(data=gr.gene_name, name="gene_id")
            gr = gr.insert(gene_id)

    if "gene_name" not in gr.columns:
        logging.warning("The gene_name field does not exist; Imputing using gene_id.")
        gene_name = pd.Series(data=gr.gene_id, name="gene_name")
        gr = gr.insert(gene_name)

    # keep only the fields we need
    gr = gr[
        [
            "Chromosome",
            "Feature",
            "Start",
            "End",
            "Strand",
            "gene_id",
            "gene_name",
            "transcript_id",
        ]
    ]

    # If there is any NaN in transcript_id field, quit
    # gene features don't have transcript_id, ignore

    if gr[gr.Feature == "exon"].transcript_id.isnull().any():
        # first, write an clean GTF if needed
        clean_gr = gr[np.logical_and(gr.transcript_id.notnull(), gr.Feature != "gene")]
        clean_gr.to_gtf(clean_gtf_path)
        # Then, raise a value error
        logging.critical(
            "".join(
                [
                    " Found missing value in exons' transcript ID; Cannot proceed."
                    f" An clean GTF file without missing transcript_id records is written to {clean_gtf_path}.",
                    " If needed, rerun using the clean GTF file",
                ]
            )
        )

    # Impute missing gene_id and gene_name values
    # define an object
    num_nan = gr.df[["gene_id", "gene_name"]].isnull().sum(axis=1)

    if num_nan.any():
        # create intermediate df
        gene_df = gr.df[["gene_id", "gene_name"]]

        # Firstly, write the problematic records to a file before imputation
        problematic_gtf_path = os.path.join(
            output_dir, "missing_gene_id_or_name_records.gtf"
        )
        gr[num_nan != 0].to_gtf(problematic_gtf_path)
        missing_record_msg = f" Found records with missing gene_id/gene_name field. These records are reported in {problematic_gtf_path}."

        # impute double missing using transcript_id
        double_missing = num_nan == 2
        if double_missing.sum():
            gene_df.loc[num_nan == 2, "gene_id"] = gr.transcript_id[num_nan == 2]
            gene_df.loc[num_nan == 2, "gene_name"] = gr.transcript_id[num_nan == 2]
            double_missing_msg = f" Imputed {(num_nan == 2).sum()} records with missing 'gene_id' and 'gene_name' using transcript_id."
        else:
            double_missing_msg = ""

        # If one field is missing, impute using the other
        # missing only gene_id
        gene_id_missing = gene_df["gene_id"].isnull()
        if gene_id_missing.sum():
            gene_df.loc[gene_df["gene_id"].isnull(), "gene_id"] = gene_df.loc[
                gene_df["gene_id"].isnull(), "gene_name"
            ]
            gene_id_missing_msg = (
                f" Imputed {gene_id_missing.sum()} missing gene_id using gene_name."
            )
        else:
            gene_id_missing_msg = ""

        gene_name_missing = gene_df["gene_name"].isnull()
        if gene_name_missing.sum():
            gene_df.loc[gene_df["gene_name"].isnull(), "gene_name"] = gene_df.loc[
                gene_df["gene_name"].isnull(), "gene_id"
            ]
            gene_name_missing_msg = (
                f" Imputed {gene_name_missing.sum()} missing gene_name using gene_id."
            )
        else:
            gene_name_missing_msg = ""

        # write the warning message
        logging.warning(
            "".join(
                [
                    missing_record_msg,
                    double_missing_msg,
                    gene_id_missing_msg,
                    gene_name_missing_msg,
                ]
            )
        )

        # replace the old gene_id and gene_name fields by imputed.
        gr = gr.drop(["gene_id", "gene_name"])
        gr = gr.insert(gene_df)

        # Then, records all exon records and gene records
        clean_gr = pr.concat([clean_gr, gr[gr.Feature == "exon"]])

    # check if the transcripts and genes are well defined
    # first, we get the transcript annotation from exons and from the transcript feature records
    # from GTF
    transcript_gr = gr[gr.Feature == "transcript"].sort(["transcript_id"])
    gene_gr = gr[gr.Feature == "gene"].sort(["gene_id"])

    # from exons
    transcript_bound_from_exons = (
        gr[gr.Feature == "exon"]
        .boundaries(group_by=["transcript_id", "gene_id", "gene_name"])
        .sort(["transcript_id"])
    )
    transcript_bound_from_exons.Feature = "transcript"

    gene_bound_from_exons = (
        gr[gr.Feature == "exon"]
        .boundaries(group_by=["gene_id", "gene_name"])
        .sort(["gene_id"])
    )
    gene_bound_from_exons.Feature = "gene"

    # define clean transcript and gene gr
    clean_transcript_gr = pr.PyRanges()
    clean_gene_gr = pr.PyRanges()

    # If there is no transcript or annotation, we
    # 1. give a warning,
    # 2. using exons' bounds as the bound of transcripts/genes to extract unspliced sequences
    # 3. write those bounds to the clean gtf file

    if transcript_gr.empty or gene_gr.empty:
        if transcript_gr.empty:
            logging.warning(
                "".join(
                    [
                        " The given GTF file doesn't have transcript feature records;",
                        " Imputing using exon feature records.",
                    ]
                )
            )

            transcript_gr = transcript_bound_from_exons
            clean_transcript_gr = transcript_bound_from_exons
            gr = pr.concat([gr, transcript_bound_from_exons])

        if gene_gr.empty:
            logging.warning(
                "".join(
                    [
                        " The given GTF file doesn't have gene feature records;",
                        " Imputing using exon feature records.",
                    ]
                )
            )

            gene_gr = gene_bound_from_exons
            clean_gene_gr = gene_bound_from_exons
            gr = pr.concat([gr, gene_bound_from_exons])
    else:
        # If some of them are missing, we report a warning
        # and use the transcripts' bounds in the original GTF file to extract introns,
        # but impute the missing annotations in the clean GTF file.
        # We will say in the warning message that if the users want to
        # use the annotations we generated,
        # they should rerun pyroe using the clean GTF file.

        # if some transcripts don't have exons
        if transcript_gr.length > transcript_bound_from_exons.length:
            # pyranges will ignore it anyway. Here I filter them out manually.
            transcript_gr = transcript_gr[
                transcript_gr.transcript_id.isin(
                    transcript_bound_from_exons.transcript_id
                )
            ]

            # complain
            logging.warning(" Found transcript(s) without exons; Ignored.")
            clean_transcript_gr = transcript_bound_from_exons

        # if some transcripts don't have features
        elif transcript_gr.length < transcript_bound_from_exons.length:
            clean_transcript_gr = transcript_bound_from_exons

            # complain
            logging.warning(
                "".join(
                    [
                        " Found transcripts without corresponding transcript feature record;",
                        " Those transcripts were not used to extract unspliced sequences.",
                    ]
                )
            )

        # if some genes don't have exons
        if gene_gr.length > gene_bound_from_exons.length:
            gene_gr = gene_gr[gene_gr.gene_id.isin(gene_bound_from_exons.gene_id)]

            # complain
            logging.warning("".join([" Found gene(s) without exons; Ignored."]))
            clean_gene_gr = gene_bound_from_exons

        # if some genes don't have features
        elif gene_gr.length < gene_bound_from_exons.length:
            clean_gene_gr = gene_bound_from_exons

            # complain
            logging.warning(
                "".join(
                    [
                        " Found genes without corresponding gene feature record.",
                        " Those genes were not used to extract unspliced sequences.",
                    ]
                )
            )

        # If the transcripts'/genes' bounds defined in the original GTF file
        # and those found manually (using exons' bounds) are different,
        # we report a warning and extract unspliced sequences
        # using the transcripts'/genes' annotation in the original GTF file,
        # but use manually defined transcript annotations (from their exons' bounds)
        # in the clean GTF file.

        # transcripts
        intersecting_txs = set(transcript_gr.transcript_id).intersection(
            set(transcript_bound_from_exons.transcript_id)
        )

        if not transcript_gr[transcript_gr.transcript_id.isin(intersecting_txs)][
            ["Start", "End"]
        ].df.equals(
            transcript_bound_from_exons[
                transcript_bound_from_exons.transcript_id.isin(intersecting_txs)
            ][["Start", "End"]].df
        ):
            clean_transcript_gr = transcript_bound_from_exons

            logging.warning(
                "".join(
                    [
                        " Found transcripts whose boundaries defined in their transcript feature record do not match their exons' bounds.",
                        " However, those boundaries were still used to extract unspliced sequences.",
                    ]
                )
            )

        # genes
        intersecting_gs = set(gene_gr.gene_id).intersection(
            set(gene_bound_from_exons.gene_id)
        )

        if not gene_gr[gene_gr.gene_id.isin(intersecting_gs)][
            ["Start", "End"]
        ].df.equals(
            gene_bound_from_exons[gene_bound_from_exons.gene_id.isin(intersecting_gs)][
                ["Start", "End"]
            ].df
        ):
            clean_gene_gr = gene_bound_from_exons

            logging.warning(
                "".join(
                    [
                        " Found genes whose boundaries defined in the gene feature records do not equal to their exons' bounds.",
                        " However, those boundaries were still used to extract unspliced sequences.",
                    ]
                )
            )

    # if clean_gr is not empty, write it
    clean_gr = pr.concat([clean_gr, clean_transcript_gr, clean_gene_gr])
    if not clean_gr.empty:
        clean_gr.to_gtf(clean_gtf_path)
        clean_gtf_msg = f" A clean GTF file with all issues fixed is generated at {clean_gtf_path}. If needed, please rerun using this clean GTF file."
        logging.warning(clean_gtf_msg)

    # return imputed gr
    return gr


def check_bedtools_version(bt_path):
    try:
        vstr = (
            subprocess.run([bt_path, "--version"], capture_output=True)
            .stdout.decode()
            .strip()
            .split("v")[1]
        )
        found_ver = parse_version(vstr)
        req_ver = parse_version("2.30.0")
        return found_ver >= req_ver
    except Exception as err:
        # in this case couldn't even run subprocess
        logging.warning(f" Cannot check bedtools version. The error message was: {err}")
        return False


def make_splici_txome(
    genome_path,
    gtf_path,
    read_length,
    output_dir,
    flank_trim_length=5,
    filename_prefix="splici",
    extra_spliced=None,
    extra_unspliced=None,
    dedup_seqs=False,
    no_bt=False,
    bt_path="bedtools",
    no_flanking_merge=False,
):
    """
    Construct the splici (spliced + introns) transcriptome for alevin-fry.

    Required Parameters
    ----------
    genome_path : str
        The path to a genome fasta file.

    gtf_path : str
        The path to a gtf file.

    read_length : int
        The read length of the single-cell experiment being processed.

    output_dir : str
        The output directory, where the splici reference files will
        be written.

    Optional Parameters
    ----------

    flank_trim_length : int (default: `5`)
        The flank trimming length.
        The final flank length is obtained by subtracting
        the flank_trim_length from the read_length.

    filename_prefix : str (default: `splici`)
        The file name prefix of the generated output files.
        The derived flank length will be automatically
        appended to the provided prefix.

    extra_spliced : str
        A path to a fasta file. The records in this fasta file will be
        regarded as spliced transcripts.

    extra_unspliced : str
        The path to a fasta file. The records in this fasta file will be
        regarded as introns.

    dedup_seqs : bool  (default: `False`)
        If True, the repeated sequences in the splici reference will be
        deduplicated.

    no_bt : bool (default: `False`)
        If true, biopython, instead of bedtools, will be used for
        generating splici reference files.

    bt_path : str
        The path to bedtools v2.30.0 or greater if it is not in the environment PATH.

    no_flanking_merge : bool (default: `False`)
        If true, overlapping introns caused by the added flanking length will not be merged.

    Returns
    -------
    Nothing will be returned. The splici reference files will be written
    to disk.

    Notes
    -----
    * The input GTF file will be processed before extracting unspliced sequences. If pyroe finds invalid records, a `clean_gtf.gtf` file will be generated in the specified output directory.  **Note** : The features extracted in the spliced + intronic transcriptome will not necessarily be those present in the `clean_gtf.gtf` file — as this command will prefer the input in the user-provided file wherever possible.  More specifically:
            * If the required metadata fields contain missing values, pyroe will impute them if possible, or return an error if not.
            * **Pyroe will always extract unspliced sequences according to the boundaries defined in the transcript/gene feature records unless there is no transcript/gene feature record in the GTF file.** In this case, pyroe imputes all transcripts/genes boundaries as the bounds of the corresponding exons to extract unspliced sequences.
            * If the transcript/gene feature records do not match their exon feature records, pyroe will still use transcript/gene feature records, but correct those transcript/gene feature records in the `celan_grf.gtf` according to exon feature records.
    * If using bedtools, a temp.bed and a temp.fa will be created and
        then deleted. These two files encode the introns of each gene
        and the exons of each transcript of each gene.

    """
    # Preparation

    # check flanking length
    flank_length = read_length - flank_trim_length

    if flank_length < 0:
        logging.critical(" Flank trim length cannot be larger than read length!")

    # check fasta file
    if not os.path.isfile(genome_path):
        logging.critical(" Cannot open the input fasta file!")

    # check gtf file
    if not os.path.isfile(gtf_path):
        logging.critical(" Cannot open the input gtf file!")

    # check bedtools
    if not no_bt:
        # check at the provided path
        if not check_bedtools_version(bt_path):
            # if it's not ok at the provided path, check
            # the standard system path
            if bt_path == "bedtools":
                # in this case, there's nowhere else to check
                # so give up on bedtools
                logging.warning(
                    "".join(
                        [
                            " Bedtools in the environemnt PATH is either",
                            " older than v.2.30.0 or doesn't exist.",
                            " Biopython will be used to extract sequences.",
                        ]
                    )
                )
                no_bt = True
            else:
                logging.warning(
                    " Bedtools specified by bt_path is either",
                    " older than v.2.30.0 or doesn't exist.",
                    " Trying to find bedtools in the environmental PATH.",
                )
                # if it's not ok at the standard system path
                # fallback to biopython
                if not check_bedtools_version("bedtools"):
                    logging.warning(
                        "".join(
                            [
                                " Bedtools in the environemnt PATH is either",
                                " older than v.2.30.0 or doesn't exist.",
                                " Biopython will be used to extract sequences.",
                            ]
                        )
                    )
                    no_bt = True
                # found it at the system path
                else:
                    bt_path = "bedtools"
                    logging.warning(" Using bedtools in the environmental PATH.")

    # create out folder and temp folder inside
    # create output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # specify output file names
    filename_prefix = filename_prefix + "_fl" + str(flank_length)
    out_fa = os.path.join(output_dir, filename_prefix + ".fa")
    out_t2g3col = os.path.join(output_dir, filename_prefix + "_t2g_3col.tsv")
    id2name_path = os.path.join(output_dir, "gene_id_to_name.tsv")

    # load gtf
    try:
        gr = pr.read_gtf(gtf_path)
    except Exception as err:
        # in this case couldn't even read the GTF file
        logging.error(
            "".join(
                [
                    " PyRanges failed to parse the input GTF file.",
                    " Please check the PyRanges documentation for the expected GTF format constraints at",
                    " https://pyranges.readthedocs.io/en/latest/autoapi/pyranges/readers/index.html?highlight=read_gtf#pyranges.readers.read_gtf .",
                    f" The error message was: {str(err)}",
                ]
            ),
            exc_info=True,
        )

    # check the validity of gr
    gr = check_gr(gr, output_dir)

    # write gene id to name tsv file
    gr.df[["gene_id", "gene_name"]].drop_duplicates().to_csv(
        id2name_path, sep="\t", header=False, index=False
    )

    # get introns
    # the introns() function uses inplace=True argument from pandas,
    # which will trigger an FutureWarning.
    warnings.simplefilter(action="ignore", category=FutureWarning)
    introns = gr.features.introns(by="transcript")
    warnings.simplefilter(action="default", category=FutureWarning)

    introns.Name = introns.gene_id

    if no_flanking_merge:
        introns = introns.merge(strand=True, by=["Name"], slack=0)

    introns = introns.extend(flank_length)

    if not no_flanking_merge:
        introns = introns.merge(strand=True, by=["Name"], slack=0)

    introns.Gene = introns.Name
    introns.Name = [
        "-I".join(map(str, z))
        for z in zip(
            introns.Name,
            introns.Name.groupby(introns.Name)
            .cumcount()
            .astype(str)
            .replace(
                "0",
                "",
            )
            .values,
        )
    ]

    # trim outbounded introns
    with open(genome_path) as fasta_file:
        chromsize = {
            title.split()[0]: len(sequence)
            for title, sequence in SimpleFastaParser(fasta_file)
        }

    # in case the genome and gene annotaitons do not match
    # try it and raise value error if this fails
    try:
        introns = pr.gf.genome_bounds(introns, chromsize, clip=True)
    except Exception as err:
        logging.error(
            "".join(
                [
                    " Failed to refine intron bounds using genome bounds.",
                    " Please check if the input genome FASTA file and GTF file match each other, especially the chromosome names.",
                    f" The error message was: {str(err)}",
                ]
            ),
            exc_info=True,
        )  # deduplicate introns
    if dedup_seqs:
        introns.drop_duplicate_positions()
    # add splice status for introns
    introns.splice_status = "U"

    introns = introns[bed_required_fields]

    # get exons
    exons = gr[gr.Feature == "exon"]

    exons.Name = exons.transcript_id
    exons.Gene = exons.gene_id
    exons = exons.sort(["Name", "Start", "End"])
    # add splice status for exons
    exons.splice_status = "S"
    # keep only required fields
    exons = exons[bed_required_fields]

    # concat spliced transcripts and introns as splici
    splici = pr.concat([exons, introns])

    # write to files
    # t2g_3col.tsv
    splici.df[["Name", "Gene", "splice_status"]].drop_duplicates().to_csv(
        out_t2g3col, sep="\t", header=False, index=False
    )
    # print(splici.head())
    tid2strand = dict(zip(splici.Name, splici.Strand))

    # splici fasta
    if not no_bt:
        try:

            # create temp folder
            temp_dir = os.path.join(output_dir, "temp")
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
            temp_fa = os.path.join(temp_dir, "temp.fa")
            temp_bed = os.path.join(temp_dir, "temp.bed")

            # write bed file
            splici.to_bed(temp_bed, keep=True)

            # run bedtools, ignore strand for now
            bt_r = subprocess.run(
                " ".join(
                    [
                        bt_path,
                        "getfasta",
                        "-fi",
                        genome_path,
                        "-fo",
                        temp_fa,
                        "-bed",
                        temp_bed,
                        # "-s",
                        "-nameOnly",
                    ]
                ),
                shell=True,
                capture_output=True,
            )

            # check return code
            if bt_r.returncode != 0:
                logging.exception(" Bedtools failed.", exc_info=True)

            # parse temp fasta file to concat exons of each transcript
            ei_parser = SeqIO.parse(temp_fa, "fasta")
            prev_rec = next(ei_parser)
            # prev_rec.id = prev_rec.id.split("(")[0]
            prev_rec.description = ""
            with open(out_fa, "w") as out_handle:

                for seq_record in ei_parser:
                    # seq_record.id = seq_record.id.split("(")[0]
                    seq_record.description = ""
                    if seq_record.id == prev_rec.id:
                        prev_rec += seq_record

                    else:
                        if tid2strand[prev_rec.id] == "-":
                            prev_rec = prev_rec.reverse_complement(
                                id=True, description=True
                            )
                        SeqIO.write(prev_rec, out_handle, "fasta")
                        prev_rec = seq_record
                # Don't forget our last customer
                if tid2strand[prev_rec.id] == "-":
                    prev_rec = prev_rec.reverse_complement(id=True, description=True)
                SeqIO.write(prev_rec, out_handle, "fasta")
            shutil.rmtree(temp_dir, ignore_errors=True)
        except Exception as err:
            no_bt = True
            logging.warning(
                f" Bedtools failed; Using biopython instead. The error message was: \n{err}"
            )
            shutil.rmtree(temp_dir, ignore_errors=True)

    if no_bt:
        with open(out_fa, "w") as out_handle:
            # read fasta, process a chromosome at a time
            for seq_record in SeqIO.parse(genome_path, "fasta"):
                # get all records on that chromosome
                chr_records = introns[introns.Chromosome == seq_record.id].df
                if not chr_records.empty:
                    chr_records.Strand = chr_records.Strand.replace(
                        ["+", "-"], [+1, -1]
                    )
                    # init seq list
                    intron_seqs = []
                    # for each intron record
                    for (idx, intron_record) in chr_records.iterrows():
                        # create Seqeture object for extracting sequence from chromosome
                        intron_feature = SeqFeature(
                            FeatureLocation(intron_record.Start, intron_record.End),
                            type="intron",
                            id=intron_record.Name,
                        )
                        intron_feature.strand = intron_record.Strand

                        # append the intron sequence to the seq list, specify name as well
                        intron_seqs.append(
                            SeqRecord(
                                intron_feature.extract(seq_record).seq,
                                id=intron_record.Name,
                                description="",
                            )
                        )
                    # finally, write all intron sequences at once.
                    SeqIO.write(intron_seqs, out_handle, "fasta")

                # Then, process spliced transcripts
                chr_records = exons[exons.Chromosome == seq_record.id].df
                if not chr_records.empty:
                    txp_seqs = []
                    # as spliced txps are the concat of all exon sequences, fist get the sequence of each exon separately,then sum them up.
                    for (tid, exon_records) in chr_records.groupby("Name"):

                        # init exon seq list
                        exon_seqs = []
                        # get the sequence of each exon
                        for (idx, exon_record) in exon_records.iterrows():
                            # create SeqFeature object for the exon record
                            # ignore strand for now, get reverse complement later if needed
                            exon_feature = SeqFeature(
                                FeatureLocation(exon_record.Start, exon_record.End),
                                type="exon",
                            )
                            # extract exon sequence from chromosome and append to exon seq list
                            exon_seqs.append(
                                SeqRecord(
                                    exon_feature.extract(seq_record).seq,
                                    id=tid,
                                    description="",
                                )
                            )
                        # append the txp sequence to spliced txp seq list
                        # consider strand
                        if tid2strand[tid] == "-":
                            txp_seqs.append(
                                sum(exon_seqs, Seq("")).reverse_complement(
                                    id=True, description=True
                                )
                            )
                        else:
                            txp_seqs.append(sum(exon_seqs, Seq("")))
                    # write all spliced transcript serquence at once.
                    SeqIO.write(txp_seqs, out_handle, "fasta")

    # append extra spliced transcript onto splici
    if extra_spliced is not None:
        append_extra(extra_spliced, out_fa, out_t2g3col, id2name_path, "S")

    # append extra unspliced transcript onto splici
    if extra_unspliced is not None:
        append_extra(extra_unspliced, out_fa, out_t2g3col, id2name_path, "U")

    if dedup_seqs:
        # Note: out_fa is intentionally passed as both
        # the input and output file name parameters because
        # we want to overwirte the duplicate fasta with
        # the deduplicated fasta.
        dedup_sequences(output_dir, out_fa, out_fa)


def make_spliceu_txome(
    genome_path,
    gtf_path,
    output_dir,
    filename_prefix,
    extra_spliced=None,
    extra_unspliced=None,
    dedup_seqs=False,
    no_bt=False,
    bt_path="bedtools",
):
    """
    Construct the spliceu (spliced + unspliced) transcriptome for alevin-fry.

    Required Parameters
    ----------
    genome_path : str
        The path to a genome fasta file.

    gtf_path : str
        The path to a gtf file.

    output_dir : str
        The output directory, where the spliceu reference files will
        be written.

    Optional Parameters
    ----------

    filename_prefix : str (default: `spliceu`)
        The file name prefix of the generated output files.
        The derived flank length will be automatically
        appended to the provided prefix.

    extra_spliced : str
        A path to a fasta file. The records in this fasta file will be
        regarded as spliced transcripts.

    extra_unspliced : str
        The path to a fasta file. The records in this fasta file will be
        regarded as introns.

    dedup_seqs : bool  (default: `False`)
        If True, the repeated sequences in the spliceu reference will be
        deduplicated.

    no_bt : bool (default: `False`)
        If true, biopython, instead of bedtools, will be used for
        generating spliceu reference files.

    bt_path : str
        The path to bedtools v2.30.0 or greater if it is not in the environment PATH.

    Returns
    -------
    Nothing will be returned. The spliceu reference files will be written
    to disk.

    Notes
    -----
    * The input GTF file will be processed before extracting unspliced sequences. If pyroe finds invalid records, a `clean_gtf.gtf` file will be generated in the specified output directory.  **Note** : The features extracted in the spliced + unspliced transcriptome will not necessarily be those present in the `clean_gtf.gtf` file — as this command will prefer the input in the user-provided file wherever possible.  More specifically:
            * If the required metadata fields contain missing values, pyroe will impute them if possible, or return an error if not.
            * **Pyroe will always extract unspliced sequences according to the boundaries defined in the transcript/gene feature records unless there is no transcript/gene feature record in the GTF file.** In this case, pyroe imputes all transcripts/genes boundaries as the bounds of the corresponding exons to extract unspliced sequences.
            * If the transcript/gene feature records do not match their exon feature records, pyroe will still use transcript/gene feature records, but correct those transcript/gene feature records in the `celan_grf.gtf` according to exon feature records.
    * If using bedtools, a temp.bed and a temp.fa will be created and
        then deleted. These two files encode the introns of each gene
        and the exons of each transcript of each gene.

    """
    # Preparation
    # check fasta file
    if not os.path.isfile(genome_path):
        logging.critical(" Cannot open the input fasta file!")

    # check gtf file
    if not os.path.isfile(gtf_path):
        logging.critical(" Cannot open the input gtf file!")

    # check bedtools
    if not no_bt:
        # check at the provided path
        if not check_bedtools_version(bt_path):
            # if it's not ok at the provided path, check
            # the standard system path
            if bt_path == "bedtools":
                # in this case, there's nowhere else to check
                # so give up on bedtools
                logging.warning(
                    "".join(
                        [
                            " Bedtools in the environemnt PATH is either",
                            " older than v.2.30.0 or doesn't exist.",
                            " Biopython will be used.",
                        ]
                    )
                )
                no_bt = True
            else:
                logging.warning(
                    "".join(
                        [
                            "Bedtools specified by bt_path is either",
                            "older than v.2.30.0 or doesn't exist.",
                            "Trying to find bedtools in the environmental PATH.",
                        ]
                    )
                )
                # if it's not ok at the standard system path
                # fallback to biopython
                if not check_bedtools_version("bedtools"):
                    logging.warning(
                        "".join(
                            [
                                " Bedtools in the environemnt PATH is either",
                                " older than v.2.30.0 or doesn't exist.",
                                " Biopython will be used.",
                            ]
                        )
                    )
                    no_bt = True
                # found it at the system path
                else:
                    bt_path = "bedtools"
                    logging.warning(" Using bedtools in the environmental PATH.")

    # create out folder and temp folder inside
    # create output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # specify output file names
    out_fa = os.path.join(output_dir, filename_prefix + ".fa")
    out_t2g3col = os.path.join(output_dir, filename_prefix + "_t2g_3col.tsv")
    out_t2g = os.path.join(output_dir, filename_prefix + "_t2g.tsv")
    out_g2g = os.path.join(output_dir, filename_prefix + "_g2g.tsv")
    id2name_path = os.path.join(output_dir, "gene_id_to_name.tsv")

    # load gtf
    try:
        gr = pr.read_gtf(gtf_path)
    except ValueError as err:
        # in this case couldn't even run subprocess
        logging.error(
            "".join(
                [
                    " PyRanges failed to parse the input GTF file.",
                    " Please check the PyRanges documentation for",
                    " the expected GTF format constraints at",
                    " https://pyranges.readthedocs.io/en/latest/autoapi/pyranges/readers/index.html?highlight=read_gtf#pyranges.readers.read_gtf .",
                    f" The error message was: {str(err)}",
                ]
            ),
            exc_info=True,
        )

    # check the validity of gr
    gr = check_gr(gr, output_dir)

    # write gene id to name tsv file
    gr.df[["gene_id", "gene_name"]].drop_duplicates().to_csv(
        id2name_path, sep="\t", header=False, index=False
    )

    # get unspliced
    unspliced = gr[gr.Feature == "gene"]
    # unspliced = gr.boundaries("gene_id")
    unspliced.Name = unspliced.gene_id + "-I"
    unspliced.Gene = unspliced.gene_id
    # add splice status for unspliced
    unspliced.splice_status = "U"
    # keep only required fields
    unspliced = unspliced[bed_required_fields]

    # get exons
    exons = gr[gr.Feature == "exon"]

    exons.Name = exons.transcript_id
    exons.Gene = exons.gene_id

    exons = exons.sort(["Name", "Start", "End"])
    # add splice status for exons
    exons.splice_status = "S"
    # keep only required fields
    exons = exons[bed_required_fields]

    # concat spliced transcripts and unspliced as spliceu
    spliceu = pr.concat([exons, unspliced])

    # write to files
    # t2g_3col.tsv
    t2g_3col = spliceu.df[["Name", "Gene", "splice_status"]].drop_duplicates()
    t2g_3col.to_csv(out_t2g3col, sep="\t", header=False, index=False)

    # t2g.csv
    t2g_3col[["Name", "Gene"]].to_csv(out_t2g, sep="\t", header=False, index=False)

    # g2g.csv
    t2g_3col[["Gene", "Gene"]].to_csv(out_g2g, sep="\t", header=False, index=False)

    tid2strand = dict(zip(spliceu.Name, spliceu.Strand))

    # spliceu fasta
    if not no_bt:
        try:
            # create temp folder
            temp_dir = os.path.join(output_dir, "temp")
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
            temp_fa = os.path.join(temp_dir, "temp.fa")
            temp_bed = os.path.join(temp_dir, "temp.bed")

            # write bed file
            spliceu.to_bed(temp_bed, keep=True)

            # run bedtools, ignore strand for now
            bt_r = subprocess.run(
                " ".join(
                    [
                        bt_path,
                        "getfasta",
                        "-fi",
                        genome_path,
                        "-fo",
                        temp_fa,
                        "-bed",
                        temp_bed,
                        # "-s",
                        "-nameOnly",
                    ]
                ),
                shell=True,
                capture_output=True,
            )

            # check return code
            if bt_r.returncode != 0:
                logging.exception("Bedtools failed.", exc_info=True)

            # parse temp fasta file to concat exons of each transcript
            ei_parser = SeqIO.parse(temp_fa, "fasta")
            prev_rec = next(ei_parser)
            # prev_rec.id = prev_rec.id.split("(")[0]
            prev_rec.description = ""
            with open(out_fa, "w") as out_handle:

                for seq_record in ei_parser:
                    # seq_record.id = seq_record.id.split("(")[0]
                    seq_record.description = ""
                    if seq_record.id == prev_rec.id:
                        prev_rec += seq_record

                    else:
                        if tid2strand[prev_rec.id] == "-":
                            prev_rec = prev_rec.reverse_complement(
                                id=True, description=True
                            )
                        SeqIO.write(prev_rec, out_handle, "fasta")
                        prev_rec = seq_record
                # Don't forget our last customer
                if tid2strand[prev_rec.id] == "-":
                    prev_rec = prev_rec.reverse_complement(id=True, description=True)
                SeqIO.write(prev_rec, out_handle, "fasta")
            shutil.rmtree(temp_dir, ignore_errors=True)
        except Exception as err:
            no_bt = True
            logging.warning(
                f" Bedtools failed; Using biopython instead. The error message was: {err}"
            )
            shutil.rmtree(temp_dir, ignore_errors=True)

    if no_bt:
        with open(out_fa, "w") as out_handle:
            # read fasta, process a chromosome at a time
            for seq_record in SeqIO.parse(genome_path, "fasta"):
                # get all records on that chromosome
                chr_records = unspliced[unspliced.Chromosome == seq_record.id].df
                if not chr_records.empty:
                    chr_records.Strand = chr_records.Strand.replace(
                        ["+", "-"], [+1, -1]
                    )
                    # init seq list
                    unspliced_seqs = []
                    # for each unspliced record
                    for (idx, unspliced_record) in chr_records.iterrows():
                        # create Seqeture object for extracting sequence from chromosome
                        unspliced_feature = SeqFeature(
                            FeatureLocation(
                                unspliced_record.Start, unspliced_record.End
                            ),
                            type="unspliced",
                            id=unspliced_record.Name,
                        )
                        unspliced_feature.strand = unspliced_record.Strand
                        # append the unspliced sequence to the seq list, specify name as well
                        unspliced_seqs.append(
                            SeqRecord(
                                unspliced_feature.extract(seq_record).seq,
                                id=unspliced_record.Name,
                                description="",
                            )
                        )
                    # finally, write all unspliced sequences at once.
                    SeqIO.write(unspliced_seqs, out_handle, "fasta")

                # Then, process spliced transcripts
                chr_records = exons[exons.Chromosome == seq_record.id].df
                if not chr_records.empty:
                    txp_seqs = []
                    # as spliced txps are the concat of all exon sequences, fist get the sequence of each exon separately,then sum them up.
                    for (tid, exon_records) in chr_records.groupby("Name"):

                        # init exon seq list
                        exon_seqs = []
                        # get the sequence of each exon
                        for (idx, exon_record) in exon_records.iterrows():
                            # create SeqFeature object for the exon record
                            # ignore strand for now, get reverse complement later if needed
                            exon_feature = SeqFeature(
                                FeatureLocation(exon_record.Start, exon_record.End),
                                type="exon",
                            )
                            # extract exon sequence from chromosome and append to exon seq list
                            exon_seqs.append(
                                SeqRecord(
                                    exon_feature.extract(seq_record).seq,
                                    id=tid,
                                    description="",
                                )
                            )
                        # append the txp sequence to spliced txp seq list
                        # consider strand
                        if tid2strand[tid] == "-":
                            txp_seqs.append(
                                sum(exon_seqs, Seq("")).reverse_complement(
                                    id=True, description=True
                                )
                            )
                        else:
                            txp_seqs.append(sum(exon_seqs, Seq("")))
                    # write all spliced transcript serquence at once.
                    SeqIO.write(txp_seqs, out_handle, "fasta")

    # append extra spliced transcript onto spliceu
    if extra_spliced is not None:
        append_extra(extra_spliced, out_fa, out_t2g3col, id2name_path, "S")

    # append extra unspliced transcript onto spliceu
    if extra_unspliced is not None:
        append_extra(extra_unspliced, out_fa, out_t2g3col, id2name_path, "U")

    if dedup_seqs:
        # Note: out_fa is intentionally passed as both
        # the input and output file name parameters because
        # we want to overwirte the duplicate fasta with
        # the deduplicated fasta.
        dedup_sequences(output_dir, out_fa, out_fa)
