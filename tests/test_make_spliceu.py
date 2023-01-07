from pyroe.make_txome import make_spliceu_txome
import tempfile
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def test_make_spliceu():
    location = os.path.dirname(os.path.realpath(__file__))
    test_file_dir = os.path.join(location, "data", "test_make_splici")
    genome_path = os.path.join(test_file_dir, "small_example_genome.fa")
    gtf_path = os.path.join(test_file_dir, "small_example.gtf")
    filename_prefix = "spliceu"
    extra_spliced = os.path.join(test_file_dir, "extra_spliced.txt")
    extra_unspliced = os.path.join(test_file_dir, "extra_unspliced.txt")
    output_dir = tempfile.TemporaryDirectory()
    make_spliceu_txome(
        genome_path=genome_path,
        gtf_path=gtf_path,
        extra_spliced=extra_spliced,
        extra_unspliced=extra_unspliced,
        output_dir=output_dir.name,
        filename_prefix=filename_prefix,
        no_bt=True,
    )

    out_fa = os.path.join(output_dir.name, "".join([filename_prefix, ".fa"]))
    out_t2g = os.path.join(output_dir.name, "".join([filename_prefix, "_t2g_3col.tsv"]))
    genome = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    spliceu_seqs = SeqIO.to_dict(SeqIO.parse(out_fa, "fasta"))
    t2g_colnames = ["txp_name", "gene_name", "splice_status"]
    t2g_3col = pd.read_csv(
        out_t2g, sep="\t", header=None, names=t2g_colnames
    ).set_index(t2g_colnames[0])

    # Define expected sequences
    # For each gene,
    # we have 3 spliced txps. Indexing starting from 1
    # tx1 ([1,2], [36, 45], [71,80])
    # tx2 ([46,55], [91, 100])
    # txp3 ([121,130], [156, 160], [191,200])
    # unspliced [1,200]

    chr_1_seq = genome["chr1"].seq
    chr_2_seq = genome["chr2"].seq
    # check spliced sequences
    # tx1_1
    extracted_tx1_1 = spliceu_seqs["tx1.1"].seq
    expected_tx1_1 = sum([chr_1_seq[0:2], chr_1_seq[35:45], chr_1_seq[70:80]], Seq(""))
    assert extracted_tx1_1 == expected_tx1_1

    # tx1_2
    extracted_tx1_2 = spliceu_seqs["tx1.2"].seq
    expected_tx1_2 = sum([chr_1_seq[45:55], chr_1_seq[90:100]], Seq(""))
    assert extracted_tx1_2 == expected_tx1_2

    # tx1_3
    extracted_tx1_3 = spliceu_seqs["tx1.3"].seq
    expected_tx1_3 = sum(
        [chr_1_seq[120:130], chr_1_seq[155:160], chr_1_seq[190:200]], Seq("")
    )
    assert extracted_tx1_3 == expected_tx1_3

    # tx2_1
    extracted_tx2_1 = spliceu_seqs["tx2.1"].seq
    expected_tx2_1 = sum(
        [chr_2_seq[0:2], chr_2_seq[35:45], chr_2_seq[70:80]], Seq("")
    ).reverse_complement()
    assert extracted_tx2_1 == expected_tx2_1

    # tx2_2
    extracted_tx2_2 = spliceu_seqs["tx2.2"].seq
    expected_tx2_2 = sum(
        [chr_2_seq[45:55], chr_2_seq[90:100]], Seq("")
    ).reverse_complement()
    assert extracted_tx2_2 == expected_tx2_2

    # tx2_3
    extracted_tx2_3 = spliceu_seqs["tx2.3"].seq
    expected_tx2_3 = sum(
        [chr_2_seq[120:130], chr_2_seq[155:160], chr_2_seq[190:200]], Seq("")
    ).reverse_complement()
    assert extracted_tx2_3 == expected_tx2_3

    # check intron seqs
    # g1_I [0,200]

    extract_g1_I = spliceu_seqs["g1-I"].seq
    expected_g1_I = sum([chr_1_seq[0:200]], Seq(""))
    assert extract_g1_I == expected_g1_I

    # g2_I [0,200]

    extract_g2_I = spliceu_seqs["g2-I"].seq
    expected_g2_I = sum([chr_2_seq[0:200]], Seq("")).reverse_complement()
    assert extract_g2_I == expected_g2_I

    # check t2g_3col
    t2g_dict = t2g_3col[[t2g_colnames[1]]].to_dict()["gene_name"]
    t2s_dict = t2g_3col[[t2g_colnames[2]]].to_dict()["splice_status"]

    # check txp name to gene name mapping
    # for g1
    assert t2g_dict["tx1.1"] == "g1"
    assert t2g_dict["tx1.2"] == "g1"
    assert t2g_dict["tx1.3"] == "g1"
    assert t2g_dict["g1-I"] == "g1"

    # for g2
    assert t2g_dict["tx2.1"] == "g2"
    assert t2g_dict["tx2.2"] == "g2"
    assert t2g_dict["tx2.3"] == "g2"
    assert t2g_dict["g2-I"] == "g2"

    # check txp name to splice status mapping
    assert t2s_dict["tx1.1"] == "S"
    assert t2s_dict["tx1.2"] == "S"
    assert t2s_dict["tx1.3"] == "S"
    assert t2s_dict["g1-I"] == "U"

    # for S
    assert t2s_dict["tx2.1"] == "S"
    assert t2s_dict["tx2.2"] == "S"
    assert t2s_dict["tx2.3"] == "S"
    assert t2s_dict["g2-I"] == "U"
