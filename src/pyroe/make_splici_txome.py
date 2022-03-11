
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
    bt_path = "bedtools",
    no_flanking_merge = False
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
        If true, biopython, instead of bedtools, will be used for 
        generating splici reference files.

    bt_path : str
        The path to bedtools v2.30.0 or greater if it is not in the environment PATH. 

    no_flanking_merge : bool
        If true, overlapping introns caused by the added flanking length will not be merged.


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

    from packaging.version import parse as parse_version

    from Bio import SeqIO
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


    def check_bedtools_version(bt_check_path):
        try:
            vstr = subprocess.run([bt_path, "--version"],
                                  capture_output=True).stdout.decode().strip().split("v")[1]
            found_ver = parse_version(vstr)
            req_ver = parse_version("2.30.0")
            return found_ver >= req_ver
        except:
            # in this case couldn't even run subprocess
            return False

    ## check bedtools
    if not no_bt:
        # check at the provided path
        if not check_bedtools_version(bt_path):
            # if it's not ok at the provided path, check
            # the standard system path
            if bt_path == "bedtools":
                # in this case, there's nowhere else to check
                # so give up on bedtools
                print("bedtools in the environemnt PATH is either",
                      "older than v.2.30.0 or doesn't exist.",
                      "\nBiopython will be used.")
                no_bt = True
            else:
                print("bedtools specified by bt_path is either",
                        "older than v.2.30.0 or doesn't exist.",
                        "\nTry finding bedtools in the environmental PATH.")
                # if it's not ok at the standard system path
                # fallback to biopython
                if not check_bedtools_version("bedtools"):
                    print("bedtools in the environemnt PATH is either",
                            "older than v.2.30.0 or doesn't exist.",
                            "\nBiopython will be used.")
                    no_bt = True
                # found it at the system path
                else:
                    bt_path = "bedtools"
                    print("Using bedtools in the environmental PATH.")


    ## create out folder and temp folder inside
    ### create output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    ## specify output file names
    filename_prefix = filename_prefix + "_fl" + str(flank_length)
    out_fa = os.path.join(output_dir, filename_prefix + ".fa")
    out_t2g3col = os.path.join(output_dir, filename_prefix + "_t2g_3col.tsv")

    # load gtf
    gr = pr.read_gtf(gtf_path)

    # get introns
    # the introns() function uses inplace=True argument from pandas,
    # which will trigger an FutureWarning. 
    warnings.simplefilter(action='ignore', category=FutureWarning)
    introns = gr.features.introns(by="transcript")
    introns.Name = introns.gene_id

    if no_flanking_merge:
        introns = introns.merge(strand=True, by=["Name"], slack=1)

    introns = introns.extend(flank_length)

    if not no_flanking_merge:
        introns = introns.merge(strand=True, by=["Name"], slack=1)

    introns.Gene = introns.Name
    introns.Name = ["-I".join(map(str, z)) for tid, size in introns.Name.value_counts().items() for z in zip([tid] * size, list(range(1, size + 1)))]

    ## trim outbounded introns
    with open(genome_path) as fasta_file:
        chromsize =  {title.split()[0]:len(sequence) for title, sequence in SimpleFastaParser(fasta_file)}
    introns = pr.gf.genome_bounds(introns, chromsize, clip=True)

    ## deduplicate introns
    if dedup_seqs:
        introns.drop_duplicate_positions()
    
    ## add splice status for introns
    introns.splice_status = "U"



    # get exons
    exons = gr[gr.Feature == "exon"]

    exons.Name = exons.transcript_id
    exons.Gene = exons.gene_id
    exons = exons.drop(exons.columns[~exons.columns.isin(introns.columns)].tolist())
    exons = exons.sort(["Name", "Start", "End"])
    ## add splice status for exons
    exons.splice_status = "S"

    # concat spliced transcripts and introns as splici
    splici = pr.concat([exons, introns])


    # splici = splici.sort(["Name", "Start", "End", "Gene"])
    
    # write to files
    ## t2g_3col.tsv
    splici.df[["Name", "Gene", "splice_status"]].drop_duplicates().to_csv(out_t2g3col, sep="\t", header=False, index=False)
    # print(splici.head())
    tid2strand = dict(zip(splici.Name, splici.Strand))

    # splici fasta
    if not no_bt:
        try:

            ### create temp folder
            temp_dir = os.path.join(output_dir, "temp")
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
            temp_fa = os.path.join(temp_dir, "temp.fa")
            temp_bed = os.path.join(temp_dir, "temp.bed")

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
                    

                
            



