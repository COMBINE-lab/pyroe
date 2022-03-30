# pyroe

## Introduction

[`Alevin-fry`](https://github.com/COMBINE-lab/alevin-fry) is a fast, accurate, and memory frugal quantification tool for preprocessing single-cell RNA-sequencing data. Detailed information can be found in the alevin-fry [pre-print](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2), and [paper](https://www.nature.com/articles/s41592-022-01408-3).

The `pyroe` package provides useful functions for analyzing single-cell or single-nucleus RNA-sequencing data using `alevin-fry`, which consists of

1. preparing the *splici* reference for the `USA` mode of `alevin-fry`, which will export a unspliced, a spliced, and an ambiguous molecule count for each gene within each cell.
2. processing the quantification result of `alevin-fry` into python as an [`AnnData`](https://anndata.readthedocs.io/en/latest/) object. 

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

### Full usage

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

## Processing alevin-fry quantification result

The quantification result of alevin-fry can be loaded into python by the `load_fry()` function. This function takes a output directory returned by `alevin-fry quant` command as the minimum input, and load the quantification result as an `AnnData` object. When processing USA mode result, it assumes that the data comes from a single-cell RNA-sequencing experiment. If one wants to process single-nucleus RNA-sequencing data or prepare the single-cell data for RNA-velocity analysis, the `output_format` argument should be set as `snRNA` or `velocity` correspondingly. One can also define customized output format, see the Full usage section for detail.

### Full usage

load alevin-fry quantification result into an AnnData object

#### Required Parameters

frydir : `str`
    The path to a output directory returned by alevin-fry quant command. \\
    The directory containing the alevin-fry quantification (i.e. the the quant.json file & alevin subdirectory).

#### Optional Parameters

output_format : `str` or `dict`
    A string represents one of the pre-defined output formats, which are "scRNA", "snRNA" and "velocity". \\
    If a customized format of the returned `AnnData` is needed, one can pass a Dictionary.\\
    See Notes section for details.

quiet : `bool` (default: `True`)
    True if function should be quiet.
    False if messages (including error messages) should be printed out. 
    
nonzero : `bool` (default: `False`)
    True if cells with non-zero expression value across all genes should be filtered in each layer.
    False if unexpressed genes should be kept.

#### Notes

The `output_format` argument takes either a dictionary that defines the customized format or 
a string that represents one of the pre-defined format of the returned `AnnData` object.

Each of the pre-defined formats contains a `X` field and some optional extra `AnnData.layers` 
obtained from the submatrices representing unspliced (U), spliced (S) and ambiguous (A) counts 
returned by alevin-fry. 

The following formats are defined:

* "scRNA": \
    This format is recommended for single cell RNA-sequencing experiments. 
    It returns a `X` field that contains the S+A count of each gene in each cell without any extra layers.

* "snRNA": \
    This format is recommended for single nucleus RNA-sequencing experiments. 
    It returns a `X` field that contains the U+S+A count of each gene in each cell without any extra layers.

* "raw": \
    This format uses the S count matrix as the `X` field and put the U, S, and A counts into three 
    separate layers, which are "unspliced", "spliced" and "ambiguous".

* "velocity": \
    This format is the same as "scRNA", except it contains two extra layers: the "spliced" layer, 
    which contains the S+A counts, and the "unspliced" layer, which contains the U counts.

A custom output format can be defined using a Dictionary specifying the desired format of the output `Anndata` object.  
If the input is not a USA mode quantification directory, this parameter is ignored
and the count matrix is returned in the `X` field of the returned `AnnData` object.  If the input
quantification directory contains a USA mode quantification, then there are 3 sub-matrices that can 
be referenced in the dictionary; 'U', 'S', 'A' containing, respectively, unspliced, spliced and 
ambiguous counts.  The dictionary should have entries of the form `key` (str) : `value` (list[str]).
The following constraints apply : there should be one key-value pair with the key `X`, the resulting
counts will be returned in the `X` field of the AnnData object. There can be an arbitrary number
of other key-value pairs, but each will be returned as a layer of the resulting AnnData object.
Within the key-value pairs, the key refers to the layer name that will be given to the combined 
count matrix upon output, and the value should be a subset of `['U', 'S', 'A']` that defines 
which sub-matrices should be summed.  For example:
`{'X' : ['S', 'A'], 'unspliced' : ['U']}`
will result in a return AnnData object where the X field has a matrix in which each entry 
corresponds to the summed spliced and ambiguous counts for each gene in each cell, and there
is an additional "unspliced" layer, whose counts are taken directly from the unspliced sub-matrix.

#### Returns
An AnnData object with X and layers corresponding to the requested `output_format`.
        
