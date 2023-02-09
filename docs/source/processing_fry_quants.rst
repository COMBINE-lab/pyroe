Processing alevin-fry quantification result
-------------------------------------------

The quantification result of alevin-fry can be loaded into python by the ``load_fry()`` function. This function takes a output directory returned by ``alevin-fry quant`` command as the minimum input, and load the quantification result as an ``AnnData`` object. When processing USA mode result, it assumes that the data comes from a single-cell RNA-sequencing experiment. If one wants to process single-nucleus RNA-sequencing data or prepare the single-cell data for RNA-velocity analysis, the ``output_format`` argument should be set as ``snRNA`` or ``velocity`` correspondingly. One can also define customized output format, see the Full Usage section for detail.

``load_fry()`` full Usage
-------------------------

load alevin-fry quantification result into an AnnData object

Required Parameters
^^^^^^^^^^^^^^^^^^^

frydir : ``str``
    The path to a output directory returned by alevin-fry quant command. \\
    The directory containing the alevin-fry quantification (i.e. the the quant.json file & alevin subdirectory).


Optional Parameters
^^^^^^^^^^^^^^^^^^^

output_format : ``str`` or ``dict``
    A string represents one of the pre-defined output formats, which are "scRNA", "snRNA" and "velocity". \\
    If a customized format of the returned `AnnData` is needed, one can pass a Dictionary.\\
    See Notes section for details.

quiet : ``bool`` (default: ``True``)
    True if function should be quiet.
    False if messages (including error messages) should be printed out. 
    
nonzero : ``bool`` (default: ``False``)
    True if cells with non-zero expression value across all genes should be filtered in each layer.
    False if unexpressed genes should be kept.

`load_fry` Notes
^^^^^^^^^^^^^^^^

The ``output_format`` argument takes either a dictionary that defines the customized format or 
a string that represents one of the pre-defined format of the returned ``AnnData`` object.

Each of the pre-defined formats contains a ``X`` field and some optional extra ``AnnData.layers`` 
obtained from the submatrices representing unspliced (U), spliced (S) and ambiguous (A) counts 
returned by alevin-fry. 

The following formats are defined:

* "scRNA": \
    This format is recommended for single cell RNA-sequencing experiments. 
    It returns a `X` field that contains the S+A count of each gene in each cell , and an extra `unspliced` field that contains the U count of each gene in each cell.

* "snRNA", "U+S+A", "all": \
    These formats are recommended for single nucleus RNA-sequencing experiments. Furthermore, these formats match the behaviors of Cell Ranger 7, which by default includes all intronic reads in the output gene count matrix for both single-cell and single-nucleus experiments.
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
