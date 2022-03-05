# import itertools
# import numpy as np
import pandas as pd
import pyranges as pr
import warnings
import os
import subprocess
import shutil

def make_splici_txome(
    gtf_path,
    genome_path,
    read_length,
    output_dir,
    flank_trim_length=5,
    file_name_prefix = "splici",
    extra_spliced = None,
    extra_unspliced = None,
    dedup_seqs = False,
    use_bt = True,
    bt_path = "bedtools"
):
    """
    Construct the splici (spliced + introns) transcriptome for alevin-fry.

    Required Parameters
    ----------
    gtf_path : str
        The path to a gtf file.

    genome_path : str
        The path to a genome fasta file.

    read_length : int
        The read length of the single-cell experiment being processed.
    
    Optional Parameters
    ----------
    output_dir : str
        The output directory, where the splici reference files will 
        be written.

    flank_trim_length : int
        The flank trimming length. 
        The final flank length is obtained by subtracting
        the flank_trim_length from the read_length.

    file_name_prefix : str
        The file name prefix of the generated output files. 
        The derived flank length will be automatically
        appended to the provided prefix.

    extra_spliced : str
        A path to a fasta file. The records in this fasta file will be 
        regarded as spliced transcripts.

    extra_unspliced : str
        The path to a fasta file. The records in this fasta file will be 
        regarded as introns.

    dedup_seqs : bool
        If True, the repeated sequences in the splici reference will be
        deduplicated.
    use_bt : bool
        If true, bedtools will be used for extracting sequences from 
        genome file.

    bt_path : str
        The path to bedtools if it is not in the environment PATH. 
        
    Returns
    -------
    Nothing will be returned. The splici reference files will be written 
    to disk. 

    Notes
    -----
    * If using bedtools, a temp.bed and a temp.fa will be created and
        then deleted. These two files encodes the introns of each gene 
        and the exons of each transcript of each gene.

    """
    # Preparation
    ## check flanking length
    flank_length = read_length - flank_trim_length
    if flank_length < 0:
        raise ValueError("Flank trim length cannot be larger than read length!")

    ## check fasta file
    if not os.path.isfile(genome_path):
        raise IOError("Cannot open the input fasta file!")
    
    ## check gtf file
    if not os.path.isfile(gtf_path):
        raise IOError("Cannot open the input gtf file!")

    ## check bedtools
    if use_bt:
        try:
            subprocess.check_call([bt_path, "--help"], 
                           stdout=subprocess.DEVNULL, 
                           stderr=subprocess.DEVNULL
                          )
        except:
            use_bt = False
            warnings.warn("Bedtools is not available. Biopython will be used to write splici fasta file.")

    ## create out folder and temp folder inside
    ### create output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ### create temp folder
    temp_dir = os.path.join(output_dir, "temp")
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    ## specify output file names
    file_name_prefix = file_name_prefix + "_fl" + str(flank_length)
    out_fa = os.path.join(output_dir, file_name_prefix + ".fa")
    out_t2g3col = os.path.join(output_dir, file_name_prefix + "_t2g_3col.tsv")
    temp_fa = os.path.join(temp_dir, "temp.fa")
    temp_bed = os.path.join(temp_dir, "temp.bed")

    # load gtf
    gr = pr.read_gtf(gtf_path)

    # get introns
    introns = gr.features.introns(by="transcript", nb_cpu=5)
    introns.Name = introns.gene_id
    introns_merged = introns.merge(strand=True, by="Name", slack=0)
    introns_merged.Gene = introns_merged.Name
    introns_merged_extended = introns_merged.extend(flank_len)
    introns_merged_extended.Name = ["-I".join(map(str, z)) for tid, size in introns_merged_extended.Name.value_counts().items() for z in zip([tid] * size, list(range(1, size + 1)))]

    ## trim outbounded introns
    with open(genome_path) as fasta_file:
        chromsize =  {title.split()[0]:len(sequence) for title, sequence in SimpleFastaParser(fasta_file)}
    introns_merged_extended = pr.gf.genome_bounds(introns_merged_extended, chromsize, clip=True)

    ## deduplicate introns
    if dedup_seqs:
        introns_merged_extended.drop_duplicate_positions()

    # get exons
    exons = gr[gr.Feature == "exon"]
    exons.Name = exons.transcript_id
    exons.Gene = exons.gene_id
    exons = exons.drop(exons.columns[~exons.columns.isin(introns_merged_extended.columns)].tolist())
    
    # concat spliced transcripts and introns as splici
    
    splici = pr.concat([exons, introns_merged_extended])
    splici = splici.sort(["Name", "Start", "End", "Gene"])
    
    # write to files
    ## t2g_3col.tsv
    splici.df[["Name", "Gene"]].drop_duplicates().to_csv(out_t2g3col, sep="\t", header=False, index=False)

    # splici fasta
    if use_bt:
        try:
            # write bed file
            splici.to_bed(temp_bed, keep=True)

            # run bedtools
            os.system(" ".join([bt_path, "getfasta",
                                  "-fi", genome_path,
                                  "-fo", temp_fa, 
                                  "-bed", temp_bed,
                                  "-s", "-nameOnly"]))

            # parse temp fasta file to concat exons of each transcript
            ei_parser = SeqIO.parse(temp_fa, "fasta")
            prev_rec = next(ei_parser)
            prev_rec.id = prev_rec.id.split("(")[0]
            prev_rec.description = ""
            with open(out_fa, "w") as out_handle:
                for seq_record in ei_parser:
                    seq_record.id = seq_record.id.split("(")[0]
                    seq_record.description = ""
                    if seq_record.id == prev_rec.id:
                        prev_rec += seq_record
                    else:
                        SeqIO.write(prev_rec, out_handle, "fasta")
                        prev_rec = seq_record
                SeqIO.write(prev_rec, out_handle, "fasta")
            shutil.rmtree(temp_dir, ignore_errors=True)
        except:
            use_bt = False
            warnings.warn("Bedtools failed, use Biopython instead.a")

    if not use_bt:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqIO.FastaIO import SimpleFastaParser
        with open(out_fa, "w") as out_handle:
            # read fasta, process a chromosome at a time
            for seq_record in SeqIO.parse(genome_path, "fasta"):
                # get all records on that chromosome
                chr_records = introns_merged_extended[introns_merged_extended.Chromosome == seq_record.id].df
                if not chr_records.empty:
                    chr_records.Strand.replace(['+', '-'],[+1, -1], inplace=True)
                    # init seq list
                    intron_seqs = []
                    # for each intron record
                    for (idx, intron_record) in chr_records.iterrows():
                        # create Seqeture object for extracting sequence from chromosome
                        intron_feature = SeqFeature(FeatureLocation(intron_record.Start, intron_record.End), type="intron", strand=intron_record.Strand, id=intron_record.Name)
                        # append the intron sequence to the seq list, specify name as well
                        intron_seqs.append(SeqRecord(intron_feature.extract(seq_record).seq, id=intron_record.Name, description=""))
                    # finally, write all intron sequences at once.
                    SeqIO.write(intron_seqs, out_handle, "fasta")

                # Then, process spliced transcripts
                chr_records = exons[exons.Chromosome == seq_record.id].df
                if not chr_records.empty:
                    chr_records.Strand.replace(['+', '-'],[+1, -1], inplace=True)
                    txp_seqs = []
                    # as spliced txps are the concat of all exon sequences, fist get the sequence of each exon separately,then sum them up.
                    for (tid, exon_records) in chr_records.groupby("Name"):
                        # init exon seq list
                        exon_seqs = []
                        # get the sequence of each exon
                        for (idx, exon_record) in exon_records.iterrows():
                            # create SeqFeature object for the exon record
                            exon_feature = SeqFeature(FeatureLocation(exon_record.Start, exon_record.End), type="exon", strand=exon_record.Strand, id=exon_record.Name)
                            # extract exon sequence from chromosome and append to exon seq list
                            exon_seqs.append(SeqRecord(exon_feature.extract(seq_record).seq, id=exon_record.Name, description=""))
                        # append the txp sequence to spliced txp seq list
                        txp_seqs.append(sum(exon_seqs, Seq("")))
                    # write all spliced transcript serquence at once.
                    SeqIO.write(txp_seqs, out_handle, "fasta")
    
