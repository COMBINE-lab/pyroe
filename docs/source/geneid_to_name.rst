Generating a gene id to gene name mapping
=========================================

It is often useful to perform analyses with gene *names* rather than gene *identifiers*. The `convert <https://pyroe.readthedocs.io/en/latest/converting_quants.html>`_ command of ``pyroe`` allows you to specify an id to name mapping so that the converted output matrix will be labeled with gene names rather than identifiers.  However, you must provide it with a 2-column tab-separated file mapping IDs to names.  This command can help you with that task.

The ``id-to-name`` command takes as input 2 parameters, an annotation file and the location where an output file should be written. The annotation file can be either a `GTF <https://mblab.wustl.edu/GTF22.html>`_ or `GFF3 <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>`_ file (which, optionally, can be gzipped).  The program will attempt to figure out the format automatically from the file suffix(es); if it cannot, you may also provide the parameter ``--format`` taking the argument either ``gtf`` or ``gff3`` to specify the format of the input.  The ``id-to-name`` command will use the `PyRanges <https://pubmed.ncbi.nlm.nih.gov/31373614/>`_ package to parse the file and extract the gene id to gene name mapping, and will write this mapping at the provided output path.

**Note**: The ``id-to-name`` sub-command assumes that the input annotation file has a field named ``gene_id`` and another named ``gene_name``.  Please make sure the input has these fields to produce a proper output mapping. We may allow specifying a custom ``gene_id`` identifier and ``gene_name`` identifier in the future.


``id-to-name`` full usage
-------------------------

.. code:: bash

	usage: pyroe id-to-name [-h] [--format FORMAT] gtf_file output

	positional arguments:
	gtf_file         The GTF input file.
	output           The path to where the output tsv file will be written.

	optional arguments:
	-h, --help       show this help message and exit
	--format FORMAT  The input format of the file (must be either GTF or GFF3). This will be inferred from the filename, but if that fails it can be provided explicitly.