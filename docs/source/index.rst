Welcome to the documentation for pyroe
============================================

What is pyroe?
===================

The pyroe package provides useful functions for analyzing single-cell or single-nucleus RNA-sequencing data using `alevin-fry`.
The main purpose of `pyroe` is to provide the function `load_fry` to load `alevin-fry` and `simpleaf` quantification results into Python as an `anndata <http://anndata.readthedocs.io/>`_ object, so as to perform downstream analysis provided by `scanpy <https://scanpy.readthedocs.io/en/stable/index.html>`_. Moreover, `pyroe` also provides functions to fetch the pre-computed quantification results from the `quantaf <https://combine-lab.github.io/quantaf/>`_ database. 

In previous versions (before v 0.10.0), pyroe also provided the functions to construct the augmented transcriptome references. Since `simpleaf` version 0.14.0, `roers <https://github.com/COMBINE-lab/roers>`_, instead of pyroe, became the default augmented reference constructor for `alevin-fry` and `simpleaf`. If you would like to use the deprecated functions to construct the augmented references, please install an older version of pyroe. Notice that old versions of pyroe are compatitble with pandas version less than 2.0.0. So, we suggest you to install the old versions of pyroe in a conda environment with a isolated environment, so as to not affect the other packages in your system.  
 
**To note that** although pyroe is available on bioconda and can be easily installed, if you encounter any problem during installation, you can define the `load_fry` function locally in your python script by copying the function definition defined `here <https://github.com/COMBINE-lab/pyroe/blob/main/src/pyroe/load_fry.py>`_. The only dependency of `load_fry` is `scanpy <https://scanpy.readthedocs.io/en/stable/installation.html>`_. 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installing
   processing_fry_quants
   converting_quants
   fetching_processed_quants
   building_splici_index
   geneid_to_name
   LICENSE.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
