from pyroe.make_splici_txome import make_splici_txome
import tempfile
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def test_make_splici():
    location = os.path.dirname(os.path.realpath(__file__))
    test_file_dir = os.path.join(location, "data", "test_make_splici")
    genome_path = os.path.join(test_file_dir, "small_example_genome.fa")
    gtf_path = os.path.join(test_file_dir, "small_example.gtf")
    filename_prefix = "splici"
    extra_spliced = os.path.join(test_file_dir, "extra_spliced.txt")
    extra_unspliced = os.path.join(test_file_dir, "extra_unspliced.txt")
    output_dir = tempfile.TemporaryDirectory()
    read_length = 5
    flank_trim_length = 2
    make_splici_txome(
        genome_path=genome_path,
        gtf_path=gtf_path,
        read_length=read_length,
        extra_spliced=extra_spliced,
        extra_unspliced=extra_unspliced,
        flank_trim_length=flank_trim_length,
        output_dir=output_dir.name,
        filename_prefix=filename_prefix,
        no_bt=True,
    )

    out_fa = os.path.join(output_dir.name, "splici_fl3.fa")
    out_t2g = os.path.join(output_dir.name, "splici_fl3_t2g_3col.tsv")
    genome = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    # gene_annotation = pr.read_gtf(gtf_path).df
    splici_seqs = SeqIO.to_dict(SeqIO.parse(out_fa, "fasta"))
    t2g_colnames = ["txp_name", "gene_name", "splice_status"]
    t2g_3col = pd.read_csv(
        out_t2g, sep="\t", header=None, names=t2g_colnames
    ).set_index(t2g_colnames[0])

    # Define expected sequences
    # For each gene,
    # we have 3 spliced txps. Indexing starting from 1
    # tx1 ([1,2], [36, 45], [71,80]), with intron regions ([3,35], [46,70])
    # tx2 ([46,55], [91, 100]), with intron ([56,90])
    # txp3 ([121,130], [156, 160], [191,200]) with introns ([131,155],[161,190])
    # By adding flanking length, the intron regions become
    #  ([0,38], [43,73])
    #  ([53,93])
    #  ([128,158],[158,193])
    # after collapsing, we will left three introns for each gene
    # I1: [1,38]
    # I2: [43,93]
    # I3: [128,193]

    # python subset works as [)

    chr_1_seq = genome["chr1"].seq
    chr_2_seq = genome["chr2"].seq
    # check spliced sequences
    # tx1_1
    extracted_tx1_1 = splici_seqs["tx1.1"].seq
    expected_tx1_1 = sum([chr_1_seq[0:2], chr_1_seq[35:45], chr_1_seq[70:80]], Seq(""))
    assert extracted_tx1_1 == expected_tx1_1

    # tx1_2
    extracted_tx1_2 = splici_seqs["tx1.2"].seq
    expected_tx1_2 = sum([chr_1_seq[45:55], chr_1_seq[90:100]], Seq(""))
    assert extracted_tx1_2 == expected_tx1_2

    # tx1_3
    extracted_tx1_3 = splici_seqs["tx1.3"].seq
    expected_tx1_3 = sum(
        [chr_1_seq[120:130], chr_1_seq[155:160], chr_1_seq[190:200]], Seq("")
    )
    assert extracted_tx1_3 == expected_tx1_3

    # tx2_1
    extracted_tx2_1 = splici_seqs["tx2.1"].seq
    expected_tx2_1 = sum(
        [chr_2_seq[0:2], chr_2_seq[35:45], chr_2_seq[70:80]], Seq("")
    ).reverse_complement()
    assert extracted_tx2_1 == expected_tx2_1

    # tx2_2
    extracted_tx2_2 = splici_seqs["tx2.2"].seq
    expected_tx2_2 = sum(
        [chr_2_seq[45:55], chr_2_seq[90:100]], Seq("")
    ).reverse_complement()
    assert extracted_tx2_2 == expected_tx2_2

    # tx2_3
    extracted_tx2_3 = splici_seqs["tx2.3"].seq
    expected_tx2_3 = sum(
        [chr_2_seq[120:130], chr_2_seq[155:160], chr_2_seq[190:200]], Seq("")
    ).reverse_complement()
    assert extracted_tx2_3 == expected_tx2_3

    # check intron seqs
    # I1: [1,38]
    # I2: [43,93]
    # I3: [128,193]

    # g1_I [0,38]

    extract_g1_I = splici_seqs["g1-I"].seq
    expected_g1_I = sum([chr_1_seq[0:38]], Seq(""))
    assert extract_g1_I == expected_g1_I

    # g1_I1 [43,93]
    extract_g1_I1 = splici_seqs["g1-I1"].seq
    expected_g1_I1 = sum([chr_1_seq[42:93]], Seq(""))
    assert extract_g1_I1 == expected_g1_I1

    # g1_I2 [128,193]
    extract_g1_I2 = splici_seqs["g1-I2"].seq
    expected_g1_I2 = sum([chr_1_seq[127:193]], Seq(""))
    assert extract_g1_I2 == expected_g1_I2

    # g2_I [0,38]

    extract_g2_I = splici_seqs["g2-I"].seq
    expected_g2_I = sum([chr_2_seq[0:38]], Seq("")).reverse_complement()
    assert extract_g2_I == expected_g2_I

    # g2_I1 [43,93]
    extract_g2_I1 = splici_seqs["g2-I1"].seq
    expected_g2_I1 = sum([chr_2_seq[42:93]], Seq("")).reverse_complement()
    assert extract_g2_I1 == expected_g2_I1

    # g2_I2 [128,193]
    extract_g2_I2 = splici_seqs["g2-I2"].seq
    expected_g2_I2 = chr_2_seq[127:193].reverse_complement()
    assert extract_g2_I2 == expected_g2_I2

    # check t2g_3col
    t2g_dict = t2g_3col[[t2g_colnames[1]]].to_dict()["gene_name"]
    t2s_dict = t2g_3col[[t2g_colnames[2]]].to_dict()["splice_status"]

    # check txp name to gene name mapping
    # for g1
    assert t2g_dict["tx1.1"] == "g1"
    assert t2g_dict["tx1.2"] == "g1"
    assert t2g_dict["tx1.3"] == "g1"
    assert t2g_dict["g1-I"] == "g1"
    assert t2g_dict["g1-I1"] == "g1"
    assert t2g_dict["g1-I2"] == "g1"

    # for g2
    assert t2g_dict["tx2.1"] == "g2"
    assert t2g_dict["tx2.2"] == "g2"
    assert t2g_dict["tx2.3"] == "g2"
    assert t2g_dict["g2-I"] == "g2"
    assert t2g_dict["g2-I1"] == "g2"
    assert t2g_dict["g2-I2"] == "g2"

    # check txp name to splice status mapping
    assert t2s_dict["tx1.1"] == "S"
    assert t2s_dict["tx1.2"] == "S"
    assert t2s_dict["tx1.3"] == "S"
    assert t2s_dict["g1-I"] == "U"
    assert t2s_dict["g1-I1"] == "U"
    assert t2s_dict["g1-I2"] == "U"

    # for S
    assert t2s_dict["tx2.1"] == "S"
    assert t2s_dict["tx2.2"] == "S"
    assert t2s_dict["tx2.3"] == "S"
    assert t2s_dict["g2-I"] == "U"
    assert t2s_dict["g2-I1"] == "U"
    assert t2s_dict["g2-I2"] == "U"
