Welcome to the documentation for pyroe
============================================

What is pyroe?
===================

The pyroe package provides useful functions for analyzing single-cell or single-nucleus RNA-sequencing data using `alevin-fry`. Since `simpleaf` version 0.14.0, `roers <https://github.com/COMBINE-lab/roers>`_, instead of pyroe, became as the default augmented reference constructor for `alevin-fry` and `simpleaf`. Now, the main purpose of `pyroe` is to provide the function `load_fry` to load `alevin-fry` quantification results into Python as an `anndata <http://anndata.readthedocs.io/>`_ object, so as to be compatible with `scanpy <https://scanpy.readthedocs.io/en/stable/index.html>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installing
   building_splici_index
   processing_fry_quants
   fetching_processed_quants
   geneid_to_name
   converting_quants
   LICENSE.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
