# pyroe

## Introduction

[`Alevin-fry`](https://github.com/COMBINE-lab/alevin-fry) is a fast, accurate, and memory frugal quantification tool for preprocessing single-cell RNA-sequencing data. Detailed information can be found in the alevin-fry [pre-print](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2), and [paper](https://www.nature.com/articles/s41592-022-01408-3).

The `pyroe` package provides useful functions for preparing input files required by `alevin-fry`, which consists of

1. preparing the *splici* reference for the `USA` mode of `alevin-fry`, which will export a unspliced, a spliced, and an ambiguous molecule count for each gene within each cell.

## Installation
The `pyroe` package can be accessed from its [github repository](https://github.com/COMBINE-lab/pyroe), installed via [`pip`](https://pip.pypa.io/en/stable/). To install the `pyroe` package via `pip` use the command:

```
pip install pyroe
```

## Preparing a splici index for quantification with alevin-fry

The USA mode in alevin-fry requires a special index reference, which is called the *splici* reference. The *splici* reference contains the spliced transcripts plus the intronic sequences of each gene. The `make_splici_txome()` function is designed to make the *splici* reference by taking a genome FASTA file and a gene annotation GTF file as the input. Details about the *splici* can be found in Section S2 of the supplementary file of the [alevin-fry paper](https://www.nature.com/articles/s41592-022-01408-3). To run pyroe, you also need to specify the read length argument `read_length` of the experiment you are working on and the flank trimming length `flank_trim_length`. A final flank length will be computed as the difference between the read_length and flank trimming length and will be attached to the ends of each intron to absorb the intron-exon junctional reads.

Following is an example of calling the `pyroe` to make the *splici* index reference. The final flank length is calculated as the difference between the read length and the flank_trim_length, i.e., 5-2=3. This function allows you to add extra spliced and unspliced sequences to the *splici* index, which will be useful when some unannotated sequences, such as mitochondrial genes, are important for your experiment. **Note** : to make `pyroe` work more quickly, it is recommended to have the latest version of [`bedtools`](https://bedtools.readthedocs.io/en/latest/) ([Aaron R. Quinlan and Ira M. Hall, 2010](https://doi.org/10.1093/bioinformatics/btq033)) installed.

```
pyroe make-splici extdata/small_example_genome.fa extdata/small_example.gtf 5 splici_txome \
      --flank-trim-length 2 --filename-prefix transcriptome_splici --dedup-seqs
```

The `pyroe` program writes two files to your specified output directory `output_dir`. They are 
- A FASTA file that stores the extracted splici sequences.
- A three columns' transcript-name-to-gene-name file that stores the name of each transcript in the splici index reference, their corresponding gene name, and the splicing status (`S` for spliced and `U` for unspliced) of those transcripts.

## Full usage

```
usage: pyroe make-splici [-h] [--filename-prefix FILENAME_PREFIX]
                         [--flank-trim-length FLANK_TRIM_LENGTH]
                         [--extra-spliced EXTRA_SPLICED]
                         [--extra-unspliced EXTRA_UNSPLICED]
                         [--bt-path BT_PATH] [--dedup-seqs] [--no-bt]
                         [--no-flanking-merge]
                         genome-path gtf-path read-length output-dir

positional arguments:
  genome-path           The path to a gtf file.
  gtf-path              The path to a gtf file.
  read-length           The read length of the single-cell experiment 
                          being processed (determines flank size).
  output-dir            The output directory where splici reference 
                          files will be written.

optional arguments:
  -h, --help            show this help message and exit
  --filename-prefix FILENAME_PREFIX
                        The file name prefix of the generated output files.
  --flank-trim-length FLANK_TRIM_LENGTH
                        Determines the amount subtracted from the read length
                        to get the flank length.
  --extra-spliced EXTRA_SPLICED
                        The path to an extra spliced sequence fasta file.
  --extra-unspliced EXTRA_UNSPLICED
                        The path to an extra unspliced sequence fasta file.
  --bt-path BT_PATH     The path to bedtools v2.30.0 or greater.
  --dedup-seqs          A flag indicates whether identical sequences will be
                          deduplicated.
  --no-bt               A flag indicates whether to disable bedtools.
  --no-flanking-merge   A flag indicates whether introns will be merged after
                          adding flanking length.
```

### the *splici* index

The *splici* index of a given species consists of the transcriptome of the species, i.e., the spliced transcripts, and the intronic sequences of the species. Within a gene, if the flanked intronic sequences overlap with each other, the overlapped intronic sequences will be collapsed as a single intronic sequence to make sure each base will appear only once in the intronic sequences. For more detailed information, please check the section S2 in the supplementary file of [alevin-fry manuscript](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2).





