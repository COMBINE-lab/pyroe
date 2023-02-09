#################################################################################
Preparing an expanded transcriptomic reference for quantification with alevin-fry
#################################################################################

The USA mode in alevin-fry requires an expanded index reference, in which sequences represent spliced and unspliced transcripts. Pyroe provides CLI programs and python functions to build the pre-defined expanded references, the spliced + intronic (*splici*) reference, which includes the spliced transcripts plus the (merged and collapsed) intronic sequences of each gene and the spliced + unspliced (*spliceu*) reference, which consists of the spliced transcripts plus the unspliced transcript (genes' entire genomic interval) of each gene. The ``make_splici_txome()`` and ``make_spliceu_txome()`` python functions are designed to make the *splici* and *spliceu* reference by taking a genome FASTA file and a gene annotation GTF file as the input.

Preparing a *spliced+intronic* transcriptome reference
-------------------------------------------------------

The *splici* index reference of a given species consists of the transcriptome of the species, i.e., the spliced transcripts and the intronic sequences of the species. Within a gene, if the flanked intronic sequences overlap with each other, the overlapped intronic sequences will be collapsed as a single intronic sequence to make sure each base will appear only once in the intronic sequences of each gene.

To prepare a *splici* reference using pyroe, in addition to a genome FASTA file and a gene annotation GTF file, you also need to specify the read length argument ``read_length`` of the experiment you are working on. Besides, you can provide a flank trimming length ``flank_trim_length``. Then, a final flank length will be computed as the difference between the ``read_length`` and flank trimming length and will be attached to the ends of each intron to absorb the intron-exon junctional reads. Details about *splici* can be found in the supplementary section S2 of the `alevin-fry paper <https://www.nature.com/articles/s41592-022-01408-3>`_. 

Following is an example of calling `pyroe` in the command line to make the *splici* index reference. The final flank length is calculated as the difference between the read length and the flank_trim_length, i.e., 5-2=3. This program allows you to add extra spliced and unspliced sequences to the *splici* index reference, which will be useful when some unannotated sequences, such as mitochondrial genes, are important for your experiment. **Note** : To make `pyroe` work more quickly, it is recommended to have the latest version of `bedtools <https://bedtools.readthedocs.io/en/latest/>`_ (`Aaron R. Quinlan and Ira M. Hall, 2010 <https://doi.org/10.1093/bioinformatics/btq033>`_) installed.

.. code:: bash

  pyroe make-spliced+intronic extdata/small_example_genome.fa extdata/small_example.gtf 5 splici_txome \
  --flank-trim-length 2 --filename-prefix splici --dedup-seqs


The `pyroe make-spliced+intronic` program writes three files to your specified output directory `splici_txome`. They are 

* A FASTA file that stores the extracted splici sequences.
* A three-column transcript-name-to-gene-name file that stores the name of each reference sequence in the splici index reference, their corresponding gene name, and the splicing status (`S` for spliced and `U` for unspliced) of those transcripts.
* A two-column TSV file that maps gene ids (used as the keys in eventual alevin-fry output) to gene names. This can later be used with the ``pyroe convert`` command line program to convert gene ids to gene names in the count matrix.

Full usage
^^^^^^^^^^

.. code:: bash 

  usage: pyroe make-spliced+intronic [-h] [--filename-prefix FILENAME_PREFIX]
                                    [--flank-trim-length FLANK_TRIM_LENGTH]
                                    [--extra-spliced EXTRA_SPLICED]
                                    [--extra-unspliced EXTRA_UNSPLICED]
                                    [--bt-path BT_PATH] [--no-bt]
                                    [--dedup-seqs] [--no-flanking-merge]
                                    [--write-clean-gtf]
                                    genome-path gtf-path read-length output-dir

  positional arguments:
    genome-path           The path to a genome fasta file.
    gtf-path              The path to a gtf file.
    read-length           The read length of the single-cell experiment being
                          processed (determines flank size).
    output-dir            The output directory where splici reference files will
                          be written.

  options:
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
    --bt-path BT_PATH     The path to beadtools v2.30.0 or greater.
    --no-bt               A flag indicates whether bedtools will be used for
                          generating splici reference files.
    --dedup-seqs          A flag indicates whether identical sequences will be
                          deduplicated.
    --no-flanking-merge   A flag indicates whether flank lengths will be
                          considered when merging introns.
    --write-clean-gtf     A flag indicates whether a clean gtf will be written
                          if encountered invalid records.

Preparing a *spliced+unspliced* index reference
-----------------------------------------------

Recently, `He et al., 2023 <https://www.biorxiv.org/content/10.1101/2023.01.04.522742>`_ introduced the spliced + unspliced (*spliceu*) index in alevin-fry. This requires the spliced + unspliced transcriptome, where the unspliced transcripts of each gene represent the entire genomic interval of that gene. Details about the *spliceu* can be found in `the preprint <https://www.biorxiv.org/content/10.1101/2023.01.04.522742>`_. To make the spliceu reference using pyroe, one can call the ``make_spliceu_txome()`` python function or ``pyroe make-spliced+unspliced`` or its alias ``pyroe make-spliceu`` from the command line. The following example shows the shell command of building a spliceu reference from a given reference set in the directory ``spliceu_txome``.

.. code:: bash

  pyroe make-spliced+unspliced extdata/small_example_genome.fa extdata/small_example.gtf spliceu_txome \
  --filename-prefix spliceu

Full usage
^^^^^^^^^^

.. code::

  usage: pyroe make-spliced+unspliced [-h] [--filename-prefix FILENAME_PREFIX]
                                      [--extra-spliced EXTRA_SPLICED] [--extra-unspliced EXTRA_UNSPLICED]
                                      [--bt-path BT_PATH] [--no-bt] [--dedup-seqs]
                                      genome-path gtf-path output-dir

  positional arguments:
    genome-path           The path to a genome fasta file.
    gtf-path              The path to a gtf file.
    output-dir            The output directory where Spliceu reference files will be written.

  options:
    -h, --help            show this help message and exit
    --filename-prefix FILENAME_PREFIX
                          The file name prefix of the generated output files.
    --extra-spliced EXTRA_SPLICED
                          The path to an extra spliced sequence fasta file.
    --extra-unspliced EXTRA_UNSPLICED
                          The path to an extra unspliced sequence fasta file.
    --bt-path BT_PATH     The path to bedtools v2.30.0 or greater.
    --no-bt               A flag indicates whether bedtools will be used for generating Spliceu reference
                          files.
    --dedup-seqs          A flag indicates whether identical sequences will be deduplicated.



Notes on the input gene annotation GTF files for building an expanded reference
----------------------------------------------------------------------------------
Pyroe builds expanded transcriptome references, the spliced + intronic (*splici*) and the spliced + unspliced (*spliceu*) transcriptome reference, based on a genome build FASTA file and a gene annotation GTF file.

The input GTF file will be processed before extracting unspliced sequences. If pyroe finds invalid records, a ``clean_gtf.gtf`` file will be generated in the specified output directory.  **Note** : The features extracted in the spliced + unspliced transcriptome will not necessarily be those present in the ``clean_gtf.gtf`` file — as this command will prefer the input in the user-provided file wherever possible. One can rerun pyroe using the ``clean_gtf.gtf`` file if needed. More specifically:

#. The non-gene level records, those whose ``feature`` field value is not "gene, " must have a valid ``transcript_id``. If this is not satisfied, pyroe returns an error and writes only the records with a valid ``transcript_id`` to the ``clean_gtf.gtf`` file. One can rerun pyroe using the `clean_gtf.gtf` file to ignore those invalid records if needed.

#. For ``gene_id`` and ``gene_name`` metadata field, 

    * If these two fields are entirely missing in the GTF file, An error will be returned. At the same time, in the ``clean_gtf.gtf``, the two fields will be imputed using the ``transcript_id`` field.
    * If one of ``gene_id`` and ``gene_name`` is completely missing, pyroe will print a warning, impute the missing field using the other one, and move to the next step with the imputed data.
    * if some records have missing ``gene_id``, ``gene_name``, or both, pyroe will print a warning and move to the next step after imputing the missing values by the following rules: For records missing ``gene_id`` or ``gene_name``, pyroe imputes the missing one using the other one; If both are missing, pyroe imputes both of them using its ``transcript_id``, which cannot be missing. 
  
#. If the GTF file does not contain transcript or gene level records, those whose ``feature`` field value is "transcript" or "gene", pyroe will print a warning and impute those missing records using the exon level records of transcripts and genes, in which the ``Start`` and ``End`` fields will be imputed as the bounds of the corresponding exons.

#. If the boundaries of transcripts/genes defined in the "transcript" or "gene" level records -- those whose ``feature`` field value is either "transcript" or "gene" -- do not match those implied by their exons' feature records, or the transcript/gene level records of some transcripts/genes' are missing, pyroe will report a warning, fix all those gene/transcript level records using their exon level records and write them to the ``clean_gtf.gtf`` file, but still extract unspliced sequences based on the existing transcript/gene level records.

