def make_splici_txome(
    genome_path,
    gtf_path,
    read_length,
    output_dir,
    flank_trim_length=5,
    filename_prefix = "splici",
    extra_spliced = None,
    extra_unspliced = None,
    dedup_seqs = False,
    no_bt = False,
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

    filename_prefix : str
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
    no_bt : bool
        If true, biopython, instead of bedtools, will not be used for 
        extracting sequences from genome file.

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
    
    import pandas as pd
    import pyranges as pr
    import warnings
    import os
    import subprocess
    import shutil

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqIO.FastaIO import SimpleFastaParser

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
    if not no_bt:
        try:
            # try default setting 
            [btv_major, btv_minor, btv_patch] = subprocess.run([bt_path, "--version"], 
                        capture_output=True
                        ).stdout.decode().strip().split("v")[1].split(".")
            # Check version
            if int(btv_major) < 2 & int(btv_minor) < 30:
                raise ValueError("Old bedtools found.")
        except:
            # if bt_path is user-defined, try bedtools in the default path
            if bt_path != "bedtools":
                print("bedtools specified by bt_path is either",
                        "older than v.2.30.0 or doesn't exist.",
                        "\nTry finding bedtools in the environmental PATH.")
                try:
                    # try bedtools in default path
                    [btv_major, btv_minor, btv_patch] = subprocess.run(["bedtools", "--version"], 
                                capture_output=True
                            ).stdout.decode().strip().split("v")[1].split(".")
                    # check version
                    if int(btv_major) < 2 & int(btv_minor) < 30:
                        raise Exception("Please update bedtools to at least v2.30.0!")
                except:
                    # If bedtools in the default path doesn't work
                    # use biopython to write fasta file.
                    print("bedtools in the environemnt PATH is either",
                            "older than v.2.30.0 or doesn't exist.",
                            "\nBiopython will be used.")
                    no_bt = True
                else:
                    bt_path = "bedtools"
                    print("Use bedtools in the environmental PATH.")

    ## create out folder and temp folder inside
    ### create output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ### create temp folder
    temp_dir = os.path.join(output_dir, "temp")
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    ## specify output file names
    filename_prefix = filename_prefix + "_fl" + str(flank_length)
    out_fa = os.path.join(output_dir, filename_prefix + ".fa")
    out_t2g3col = os.path.join(output_dir, filename_prefix + "_t2g_3col.tsv")
    temp_fa = os.path.join(temp_dir, "temp.fa")
    temp_bed = os.path.join(temp_dir, "temp.bed")

    # load gtf
    gr = pr.read_gtf(gtf_path)

    # get introns
    introns = gr.features.introns(by="transcript", nb_cpu=5)
    introns.Name = introns.gene_id
    introns_merged = introns.merge(strand=True, by=["Name"], slack=0)
    introns_merged.Gene = introns_merged.Name
    introns_merged_extended = introns_merged.extend(flank_length)
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
    exons = exons.sort(["Name", "Start", "End"])
    
    # concat spliced transcripts and introns as splici
    splici = pr.concat([exons, introns_merged_extended])
    # splici = splici.sort(["Name", "Start", "End", "Gene"])
    
    # write to files
    ## t2g_3col.tsv
    splici.df[["Name", "Gene"]].drop_duplicates().to_csv(out_t2g3col, sep="\t", header=False, index=False)
    # print(splici.head())
    tid2strand = dict(zip(splici.Name, splici.Strand))

    # splici fasta
    if not no_bt:
        try:
            # write bed file
            splici.to_bed(temp_bed, keep=True)

            # run bedtools, ignore strand for now
            bt_r = subprocess.run(" ".join([bt_path, "getfasta",
                                            "-fi", genome_path,
                                            "-fo", temp_fa, 
                                            "-bed", temp_bed,
                                            # "-s", 
                                            "-nameOnly"]),
                                    shell=True,
                                    capture_output=True)
            # check return code
            if bt_r.returncode != 0:
                raise ValueError("Bedtools failed.")
            
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
                        if tid2strand[prev_rec.id] == "-":
                            prev_rec = prev_rec.reverse_complement(id=True, description=True)
                        SeqIO.write(prev_rec, out_handle, "fasta")
                        prev_rec = seq_record
                # Don't forget our last customer
                if tid2strand[prev_rec.id] == "-":
                    prev_rec = prev_rec.reverse_complement(id=True, description=True)
                SeqIO.write(prev_rec, out_handle, "fasta")
            shutil.rmtree(temp_dir, ignore_errors=True)
        except:
            no_bt = True
            warnings.warn("Bedtools failed, use biopython instead.")

    if no_bt:
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
                    chr_records.Strand = chr_records.Strand.replace(['+', '-'],[+1, -1])
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
                    txp_seqs = []
                    # as spliced txps are the concat of all exon sequences, fist get the sequence of each exon separately,then sum them up.
                    for (tid, exon_records) in chr_records.groupby("Name"):
                        # init exon seq list
                        exon_seqs = []
                        # get the sequence of each exon
                        for (idx, exon_record) in exon_records.iterrows():
                            # create SeqFeature object for the exon record
                            # ignore strand for now, get reverse complement later if needed
                            exon_feature = SeqFeature(FeatureLocation(exon_record.Start, exon_record.End), type="exon")
                            # extract exon sequence from chromosome and append to exon seq list
                            exon_seqs.append(SeqRecord(exon_feature.extract(seq_record).seq, id=tid, description=""))
                        # append the txp sequence to spliced txp seq list
                        # consider strand
                        if tid2strand[tid] == "-":
                            txp_seqs.append(sum(exon_seqs, Seq("")).reverse_complement(id=True, description=True))
                        else:
                            txp_seqs.append(sum(exon_seqs, Seq("")))
                    # write all spliced transcript serquence at once.
                    SeqIO.write(txp_seqs, out_handle, "fasta")

        # append extra spliced transcript onto splici
    if extra_spliced is not None:
        ## trim outbounded introns
        with open(extra_spliced) as extra_spliced_fa:
            with open(out_fa, "a") as splici_fa:
                with open(out_t2g3col, "a") as t2g:
                    for title, sequence in SimpleFastaParser(extra_spliced_fa):
                        tid = title.split()[0]
                        # splici_fa.writelines("\n")
                        splici_fa.write("\n>"+tid)
                        splici_fa.write("\n"+sequence)
                        t2g.write("\n"+tid+"\t"+tid)
                        
    # append extra unspliced transcript onto splici
    if extra_unspliced is not None:
        ## trim outbounded introns
        with open(extra_unspliced) as extra_unspliced_fa:
            with open("test_splici_fl145.fa", "a") as splici_fa:
                with open("test_splici_t2g.tsv", "a") as t2g:
                    for title, sequence in SimpleFastaParser(extra_unspliced_fa):
                        tid = title.split()[0]
                        # splici_fa.writelines("\n")
                        splici_fa.write("\n>"+tid+"-I")
                        splici_fa.write("\n"+sequence)
                        t2g.write("\n"+tid+"-I\t"+tid)
                    

                
            



if __name__ == "__main__":
    import argparse

    # Create the parser
    parser = argparse.ArgumentParser(description='The pyroe package provides useful functions for preparing input files required by alevin-fry.',
                                        prog='pyroe')
    subparsers = parser.add_subparsers(title='subcommands',
                                        description='valid subcommands',
                                        help='additional help')
    parser_makeSplici = subparsers.add_parser('make-splici', help='Make splici reference')
    parser_makeSplici.add_argument('genome_path', metavar='genome-path', type=str, help='The path to a gtf file.')
    parser_makeSplici.add_argument('gtf_path', metavar='gtf-path', type=str, help='The path to a gtf file.')
    parser_makeSplici.add_argument('read_length', metavar='read-length', type=int, help='Read length (determines flank size).')
    parser_makeSplici.add_argument('output_dir', metavar='output-dir', type=str, help='Output directory where splici reference information will be written.')
    parser_makeSplici.add_argument('--filename-prefix', type=str, default="splici", help='The file name prefix of the generated output files.')
    parser_makeSplici.add_argument('--flank-trim-length', type=int, default=5, help='Determines the amount subtracted from the read length to get the flank length.')
    parser_makeSplici.add_argument('--extra-spliced', type=str, help='The path to an extra spliced sequence fasta file.')
    parser_makeSplici.add_argument('--extra-unspliced', type=str, help='The path to an extra unspliced sequence fasta file.')
    parser_makeSplici.add_argument('--bt-path', type=str, default="bedtools", help='The path to bedtools.')
    parser_makeSplici.add_argument('--dedup-seqs', action='store_true', help='a flag indicates whether to deduplicate identical sequences.')
    parser_makeSplici.add_argument('--no-bt', action='store_true', help='A flag indicates whether to disable bedtools.')

    # Execute the parse_args() method
    args = parser.parse_args()
    
    make_splici_txome(
    genome_path=args.genome_path,
    gtf_path=args.gtf_path,
    read_length=args.read_length,
    output_dir=args.output_dir,
    flank_trim_length=args.flank_trim_length,
    filename_prefix=args.filename_prefix,
    extra_spliced=args.extra_spliced,
    extra_unspliced=args.extra_unspliced,
    dedup_seqs=args.dedup_seqs,
    no_bt=args.no_bt,
    bt_path=args.bt_path
)
