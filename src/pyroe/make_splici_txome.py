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

    from Bio.SeqIO.FastaIO import SimpleFastaParser

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
    import os
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

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


def check_gr(gr, output_dir, write_clean_gtf):
    """
    This function checks the validity of a PyRanges object.
    It will follow these rules:
    1. For gene type records, gene_id and gene_name cannot be both missing.
    2. All other types of records have to have a valid (not NaN) transcript_id.
    3. For all types of records other than gene type,
        - If the gene_id and gene_name are both missing, then the transcript_id
    will be used to impute them.
        - If one of gene_id and gene_name is missing, then the other will be
    used to impute the missing one.

    Args:
        gr (pyranges): Stranded PyRanges object
    """

    import pandas as pd
    import os
    import pyranges as pr
    import warnings

    # split gene type records with others
    # we don't use gene records in splici construction
    gene_gr = gr[gr.Feature == "gene"]
    gr = gr[gr.Feature != "gene"]
    clean_gtf_path = os.path.join(output_dir, "clean_gtf.gtf")

    # If required fields are missing, quit
    if "transcript_id" not in gr.columns:
        raise ValueError(
            "The input GTF file doesn't contain transcript_id field; Cannot proceed."
        )

    if "gene_id" not in gr.columns:
        # use gene_name as gene_id if exists, return an error otherwise
        if "gene_name" not in gr.columns:
            raise ValueError(
                "The input GTF file doesn't contain gene_id and gene_name field; Cannot proceed."
            )
        else:
            warnings.warn("gene_id field does not exist, use gene_name as gene_id.")
            gene_id = pd.Series(data=gr.gene_name, name="gene_id")
            gr = gr.insert(gene_id)

    if "gene_name" not in gr.columns:
        warnings.warn("gene_name field does not exist, use gene_id as gene_name.")
        gene_name = pd.Series(data=gr.gene_id, name="gene_name")
        gr = gr.insert(gene_name)

    # If there is any NaN in transcript_id field, quit
    if sum(gr.transcript_id.isnull()):

        # first, write an clean GTF if needed
        if write_clean_gtf:
            gr = pr.concat([gene_gr, gr[gr.transcript_id.notnull()]])
            gr.to_gtf(clean_gtf_path)
            # gr[gr.transcript_id.notnull()].to_gtf(clean_gtf_path)
            clean_gtf_msg = f"An clean GTF file is written to {clean_gtf_path}."
        else:
            clean_gtf_msg = "Set the write_clean_gtf flag if a clean GTF without the invalid records is needed."
        # Then, raise a value error
        raise ValueError(
            f"Found NaN value in exons' transcript ID; Cannot proceed.\n{clean_gtf_msg}"
        )

    # Impute missing gene_id and gene_name values
    # define an object
    num_nan = gr.df[["gene_id", "gene_name"]].isnull().sum(axis=1)

    if num_nan.sum():

        # create intermediate df
        gene_df = gr.df[["gene_id", "gene_name"]]

        # Firstly, write the problematic records to a file before imputation
        problematic_gtf_path = os.path.join(
            output_dir, "missing_gene_id_or_name_records.gtf"
        )
        gr[num_nan > 0].to_gtf(problematic_gtf_path)
        missing_record_msg = f"\nFound records with missing gene_id/gene_name field.\nThese records are reported in {problematic_gtf_path}."

        # impute using transcript_id
        double_missing = num_nan == 2
        if double_missing.sum():
            gene_df.loc[num_nan == 2, "gene_id"] = gr.transcript_id[num_nan == 2]
            gene_df.loc[num_nan == 2, "gene_name"] = gr.transcript_id[num_nan == 2]
            double_missing_msg = f"\n  - Found {(num_nan == 2).sum()} records missing gene_id and gene_name, imputing using transcript_id."
        else:
            double_missing_msg = ""

        # If one field is missing, impute using the other
        # missing only gene_id
        gene_id_missing = gene_df["gene_id"].isnull()
        if gene_id_missing.sum():
            gene_df.loc[gene_df["gene_id"].isnull(), "gene_id"] = gene_df.loc[
                gene_df["gene_id"].isnull(), "gene_name"
            ]
            gene_id_missing_msg = f"\n  - Found {gene_id_missing.sum()} records missing gene_id, imputing using gene_name."
        else:
            gene_id_missing_msg = ""

        gene_name_missing = gene_df["gene_name"].isnull()
        if gene_name_missing.sum():
            gene_df.loc[gene_df["gene_name"].isnull(), "gene_name"] = gene_df.loc[
                gene_df["gene_name"].isnull(), "gene_id"
            ]
            gene_name_missing_msg = f"\n  - Found {gene_name_missing.sum()} records missing gene_name, imputing using gene_id."
        else:
            gene_name_missing_msg = ""

        # write the warning message
        warnings.warn(
            "".join(
                [
                    missing_record_msg,
                    double_missing_msg,
                    gene_id_missing_msg,
                    gene_name_missing_msg,
                ]
            )
        )

        # replace the old gene_id and gene_name fields using imputed one.
        gr = gr.drop(["gene_id", "gene_name"])
        gr = gr.insert(gene_df)
        # if gene_gr is used in the future, then concat them.
        if write_clean_gtf:
            clean_gr = pr.concat([gene_gr, gr])
            clean_gr.to_gtf(clean_gtf_path)
            clean_gtf_msg = f"An clean GTF file is written to {clean_gtf_path}."
            print(clean_gtf_msg)

            # gr = pr.concat([gr, gene_gr])

    # return imputed gr
    return gr


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
    write_clean_gtf=False,
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

    write_clean_gtf : bool (default: `False`)
        If true, when the input GTF contains invalid records, a clean GTF
        file `clean_gtf.gtf` with these invalid records removed will be
        exported to the output dir.
        An invalid record is an exon record without an transcript ID
        in the `transcript_id` field.


    Returns
    -------
    Nothing will be returned. The splici reference files will be written
    to disk.

    Notes
    -----
    * If using bedtools, a temp.bed and a temp.fa will be created and
        then deleted. These two files encode the introns of each gene
        and the exons of each transcript of each gene.

    """

    import pyranges as pr

    import warnings

    import os

    import subprocess
    import shutil

    from packaging.version import parse as parse_version

    from Bio import SeqIO
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    # Preparation

    # check flanking length
    flank_length = read_length - flank_trim_length

    if flank_length < 0:
        raise ValueError("Flank trim length cannot be larger than read length!")

    # check fasta file
    if not os.path.isfile(genome_path):
        raise IOError("Cannot open the input fasta file!")

    # check gtf file
    if not os.path.isfile(gtf_path):
        raise IOError("Cannot open the input gtf file!")

    def check_bedtools_version(bt_check_path):
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
        except subprocess.CalledProcessError as err:
            # in this case couldn't even run subprocess
            warnings.warn(f"Cannot check bedtools version.\n{err}")
            return False

    # check bedtools
    if not no_bt:
        # check at the provided path
        if not check_bedtools_version(bt_path):
            # if it's not ok at the provided path, check
            # the standard system path
            if bt_path == "bedtools":
                # in this case, there's nowhere else to check
                # so give up on bedtools
                print(
                    "bedtools in the environemnt PATH is either",
                    "older than v.2.30.0 or doesn't exist.",
                    "\nBiopython will be used.",
                )
                no_bt = True
            else:
                print(
                    "bedtools specified by bt_path is either",
                    "older than v.2.30.0 or doesn't exist.",
                    "\nTry finding bedtools in the environmental PATH.",
                )
                # if it's not ok at the standard system path
                # fallback to biopython
                if not check_bedtools_version("bedtools"):
                    print(
                        "bedtools in the environemnt PATH is either",
                        "older than v.2.30.0 or doesn't exist.",
                        "\nBiopython will be used.",
                    )
                    no_bt = True
                # found it at the system path
                else:
                    bt_path = "bedtools"
                    print("Using bedtools in the environmental PATH.")

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
    except ValueError:
        # in this case couldn't even run subprocess
        raise RuntimeError(
            "PyRanges failed to parse the input GTF file. Please check the PyRanges documentation for the expected GTF format constraints.\nhttps://pyranges.readthedocs.io/en/latest/autoapi/pyranges/readers/index.html?highlight=read_gtf#pyranges.readers.read_gtf"
        )

    gr = pr.read_gtf(gtf_path)

    # check the validity of gr
    gr = check_gr(gr, output_dir, write_clean_gtf)

    # write gene id to name tsv file
    gr.df[["gene_id", "gene_name"]].drop_duplicates().to_csv(
        id2name_path, sep="\t", header=False, index=False
    )

    # get introns
    # the introns() function uses inplace=True argument from pandas,
    # which will trigger an FutureWarning.
    warnings.simplefilter(action="ignore", category=FutureWarning)
    introns = gr.features.introns(by="transcript")
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
    introns = pr.gf.genome_bounds(introns, chromsize, clip=True)

    # deduplicate introns
    if dedup_seqs:
        introns.drop_duplicate_positions()

    # add splice status for introns
    introns.splice_status = "U"

    # get exons
    exons = gr[gr.Feature == "exon"]

    exons.Name = exons.transcript_id
    exons.Gene = exons.gene_id
    exons = exons.drop(exons.columns[~exons.columns.isin(introns.columns)].tolist())
    exons = exons.sort(["Name", "Start", "End"])
    # add splice status for exons
    exons.splice_status = "S"

    # concat spliced transcripts and introns as splici
    splici = pr.concat([exons, introns])

    # splici = splici.sort(["Name", "Start", "End", "Gene"])

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
                raise ValueError("Bedtools failed.")

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
        except subprocess.CalledProcessError as err:
            no_bt = True
            warnings.warn(f"Bedtools failed. Use biopython instead.\n{err}")
            shutil.rmtree(temp_dir, ignore_errors=True)

    if no_bt:

        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        from Bio.SeqRecord import SeqRecord

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
                            strand=intron_record.Strand,
                            id=intron_record.Name,
                        )
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
