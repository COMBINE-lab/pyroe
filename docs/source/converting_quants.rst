Converting quantification results
=================================

The ``convert`` sub-command of ``pyroe`` can convert the output of `alevin-fry` into several common formats, such as 
the native `AnnData` format (``h5ad``).  Further, when performing this conversion, it can organize the unspliced, 
spliced, and ambiguous counts as desired by the user.

The sub-command takes as input a quantification directory produced by ``alevin-fry``, and an output location.
Additionally, the user should pass in command line parameters to describe the desired output structure, and
output format. The output structure defines how the ``U``, ``S``, and ``A`` layers of the input quantification should
be represented in the converted matrix.  The syntax for this flag exactly mimics the ``output_format`` argument of
the ``load_fry`` function, which you can read about `here <https://pyroe.readthedocs.io/en/latest/building_splici_index.html#load-fry-notes>`_.
Note that, if you pass in a custom output structure, you should enclose your format description in quotes.  For
example, to output to an object where the "main" layer (``X``) contains the sum of ``U``, ``S``, and ``A``, and where
there is an additional layer named `unspliced` having just the unspliced counts, you would pass
``--output-structure '{ "X" : ["U", "S", "A"], "unspliced" : ["U"]}'``. 

If you do not explicitly provide an ``--output-format``, the default of ``h5ad`` will be used.

``convert`` command full usage
------------------------------

.. code:: bash

  usage: pyroe convert [-h] [--output-structure OUTPUT_STRUCTURE] [--output-format OUTPUT_FORMAT] [--geneid-to-name GENEID_TO_NAME] quant_dir output

  positional arguments:
    quant_dir             The input quantification directory containing the matrix to be converted.
    output                The output name where the quantification matrix should be written. For `csvs` output format, this will be a directory. For all others, it will be a file.

  optional arguments:
    -h, --help            show this help message and exit
    --output-structure OUTPUT_STRUCTURE
                          The structure that U,S and A counts should occupy in the output matrix.
    --output-format OUTPUT_FORMAT
                          The format in which the output should be written, one of {'zarr', 'loom', 'csvs', 'h5ad'}.
    --geneid-to-name GENEID_TO_NAME
                          A 2 column tab-separated list of gene ID to gene name mappings. Providing this file will project gene IDs to gene names in the output.