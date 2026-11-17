Processing alevin-fry quantification result
-------------------------------------------

The quantification result of alevin-fry can be loaded into python by the ``load_fry()`` function. This function takes a output directory returned by ``alevin-fry quant`` command as the minimum input, and load the quantification result as an ``AnnData`` object. When processing USA mode result, it assumes that the data comes from a single-cell RNA-sequencing experiment. If one wants to process single-nucleus RNA-sequencing data or prepare the single-cell data for RNA-velocity analysis, the ``output_format`` argument should be set as ``snRNA`` or ``velocity`` correspondingly. One can also define customized output format, see the Full Usage section for detail.

``load_fry()`` also loads the output of a 10x Flex run (``simpleaf multiplex-quant``, which runs ``alevin-fry quant-probe``). See the section `Flex and other multiplex output`_.

.. code:: python

  from pyroe import load_fry

  # Flex: cell x gene matrix, with the sample of each cell in obs["sample_name"]
  adata = load_fry("simpleaf_multiplex_quant/af_quant")

  # Flex with --probe-mtx: cell x probe matrix, with the 10x probe-set CSV
  probes = load_fry(
      "simpleaf_multiplex_quant/af_quant",
      feature_level="probe",
      probe_ref="Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.csv",
  )

``load_fry()`` full Usage
-------------------------

load alevin-fry quantification result into an AnnData object

Required Parameters
^^^^^^^^^^^^^^^^^^^

frydir : ``str``
    The path to a output directory returned by alevin-fry quant command. \
    The directory containing the alevin-fry quantification (i.e. the the quant.json file & alevin subdirectory).


Optional Parameters
^^^^^^^^^^^^^^^^^^^

output_format : ``str`` or ``dict``
    A string represents one of the pre-defined output formats, which are "scRNA", "snRNA" and "velocity". \
    If a customized format of the returned `AnnData` is needed, one can pass a Dictionary. \
    See Notes section for details.

quiet : ``bool`` (default: ``True``)
    True if function should be quiet.\
    False if messages (including error messages) should be printed out. 
    
nonzero : ``bool`` (default: ``False``)
    True if cells with non-zero expression value across all genes should be filtered in each layer.\
    False if unexpressed genes should be kept.

feature_level : ``str`` (default: ``"gene"``)
    ``"gene"`` loads the cell x gene matrix ``alevin/quants_mat.mtx``. \
    ``"probe"`` loads the cell x probe matrix ``alevin/probe_quants_mat.mtx``, which ``alevin-fry quant-probe --probe-mtx`` (``simpleaf multiplex-quant --probe-mtx``) writes. \
    The function raises a ``ValueError`` if ``quant.json`` does not have ``"probe_mtx": true``.

probe_ref : ``str`` or ``None`` (default: ``None``)
    Only for ``feature_level="probe"``. The 10x probe-set CSV of the run. See `Flex and other multiplex output`_. \
    If it is ``None``, the function uses ``probe_set.csv`` in ``frydir``. \
    The function raises a ``FileNotFoundError`` if the file does not exist.

`load_fry` Notes
^^^^^^^^^^^^^^^^

The ``output_format`` argument takes either a dictionary that defines the customized format or 
a string that represents one of the pre-defined format of the returned ``AnnData`` object.

Each of the pre-defined formats contains a ``X`` field and some optional extra ``AnnData.layers`` 
obtained from the submatrices representing unspliced (U), spliced (S) and ambiguous (A) counts 
returned by alevin-fry. 

The following formats are defined:

* "scRNA": \
    This format is recommended for single cell RNA-sequencing experiments. \
    It returns a `X` field that contains the S+A count of each gene in each cell , and an extra `unspliced` field that contains the U count of each gene in each cell.

* "snRNA", "U+S+A", "all": \
    These formats are recommended for single nucleus RNA-sequencing experiments. Furthermore, these formats match the behaviors of Cell Ranger 7, which by default includes all intronic reads in the output gene count matrix for both single-cell and single-nucleus experiments.\
    These formats return a `X` field that contains the U+S+A count of each gene in each cell without any extra layers.

* "raw": \
    This format uses the S count matrix as the `X` field and put the U, S, and A counts into three separate layers, which are `unspliced`, `spliced` and `ambiguous`.

* "S+A": \
    This format uses the  U + S counts as the `X` field without any extra layers.

* "velocity": \
    This format is the same as "scRNA", except it contains a `spliced` layer, 
    which contains the S+A counts.

A custom output format can be defined using a Dictionary specifying the desired format of the output ``Anndata`` object.  
If the input is not a USA mode quantification directory, this parameter is ignored
and the count matrix is returned in the `X` field of the returned ``AnnData`` object.  If the input
quantification directory contains a USA mode quantification, then there are 3 sub-matrices that can 
be referenced in the dictionary; 'U', 'S', 'A' containing, respectively, unspliced, spliced and 
ambiguous counts.  The dictionary should have entries of the form ``key`` (str) : ``value`` (list[str]).
The following constraints apply : there should be one key-value pair with the key ``X``, the resulting
counts will be returned in the ``X`` field of the AnnData object. There can be an arbitrary number
of other key-value pairs, but each will be returned as a layer of the resulting AnnData object.
Within the key-value pairs, the key refers to the layer name that will be given to the combined 
count matrix upon output, and the value should be a subset of ``['U', 'S', 'A']`` that defines 
which sub-matrices should be summed.  For example:
``{'X' : ['S', 'A'], 'unspliced' : ['U']}``
will result in a return AnnData object where the X field has a matrix in which each entry 
corresponds to the summed spliced and ambiguous counts for each gene in each cell, and there
is an additional "unspliced" layer, whose counts are taken directly from the unspliced sub-matrix.

Returns
^^^^^^^

An AnnData object with X and layers corresponding to the requested ``output_format``.

Flex and other multiplex output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a multiplex run, one matrix holds the cells of all samples.
The row label (``obs_names`` and ``obs["barcodes"]``) is ``<sample_name>_<cell_barcode>``.
The same cell barcode can occur in two samples, so the row label is the cell id.
Sample names can contain ``_``.
The rows of one sample are not together, and the row order can change between runs.

For a multiplex run, ``obs["sample_name"]`` and ``obs["cell_barcode"]`` give the two parts of the row label.
They come from the ``sample_name`` and ``CB`` columns of ``featureDump.txt``, whose line *i* describes row *i* of the matrix.
The function raises a ``ValueError`` if the number or the order of the lines does not match the rows.
If ``featureDump.txt`` does not exist, and ``collation_manifest.bin`` or ``sample_info.json`` is in ``frydir``, the function splits each row label on the last ``_``.
The barcode contains only ``A``, ``C``, ``G`` and ``T``, so this split is correct.
An output without a ``sample_name`` column in ``featureDump.txt`` (for example a single-sample run) gets no new columns.

The ``gene_id_to_name.tsv`` of a Flex run can have more gene ids than the matrix columns.
The function keeps only the matrix columns (a left join).

With ``--probe-mtx``, alevin-fry also writes the cell x probe matrix ``alevin/probe_quants_mat.mtx``, which ``feature_level="probe"`` loads.
It has the same rows as the gene matrix, in the same order, and the sum of the probe counts of each gene is equal to the gene matrix.
The probe reference is the original 10x probe-set CSV.
The function uses the ``probe_ref`` argument, or else ``probe_set.csv`` in ``frydir``.

* The function skips the lines at the top that start with ``#``.
* The CSV must have the columns ``probe_id`` and ``gene_id``. The columns ``gene_name`` and ``region`` are optional.
* The probe columns are not in the order of the CSV, so the function joins the probes to the CSV by ``probe_id``. It raises a ``ValueError`` if a matrix probe is not in the CSV, or if a probe id occurs two times in the CSV.
* ``var_names`` and ``var["probe_id"]`` are the probe ids. ``var["gene_ids"]``, ``var["gene_names"]`` and ``var["region"]`` come from the ``gene_id``, ``gene_name`` and ``region`` columns. An empty value is ``NaN``.
* A 10x "unspliced" probe does not span a splice junction; it does not mean pre-mRNA.
* If the CSV has no ``gene_name`` column (older formats), the gene names come from the gene id to name file.
