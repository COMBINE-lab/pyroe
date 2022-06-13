Fetching and loading preprocessed quantification results
========================================================

The raw data for many single-cell and single-nucleus RNA-seq experiments is publicly available.  However, certain datasets are used *again and again*, to demonstrate data processing in tutorials, as benchmark datasets for novel methods (e.g. for clustering, dimensionality reduction, cell type identification, etc.).  In particular, 10x Genomics hosts various publicly available datasets generated using their technology and processed via their Cell Ranger software `on their website for download <https://www.10xgenomics.com/resources/datasets>`_.

We have created a `Nextflow <https://www.nextflow.io>`_-based ``alevin-fry`` workflow that one can use to easily quantify single-cell RNA-sequencing data in a single workflow.  The pipeline can be found `here <https://github.com/COMBINE-lab/10x-requant>`_.  To test out this initial pipeline, we have begun to reprocess the publicly-available datasets collected from the 10x website. We have focused the initial effort on standard single-cell and single-nucleus gene-expression data generated using the Chromium v2 and v3 chemistries, but hope to expand the pipeline to more complex protocols soon (e.g. feature barcoding experiments) and process those data as well.  We note that these more complex protocols can already be processed with ``alevin-fry`` (see the `alevin-fry tutorials <https://combine-lab.github.io/alevin-fry-tutorials/>`_), but these have just not yet been incorprated into the automated Nextflow-based workflow linked above.

We provide two python functions:

* ``fetch_processed_quant()`` can fetch the quantification result of one or more available datasets according to the provided ``dataset_ids`` vector, and store them to a local folder. 
* ``load_processed_quant()`` can fetch the quantification result of one or more available dataset as ``fetch_processed_quant()``, and load them into python as ``AnnData`` objects. We also provide a CLI for fetching quantification results.


.. code:: bash 

  pyroe fetch-quant 1 3 6


or start python, and run

.. code:: python 

  import pyroe

  # fetch, decompress and load the quantification result of dastset #1, 3 and 6
  pq_dict = pyroe.load_processed_quant([1,3,6])

  # get the ProcessedQuant class object for dataset #1 and #3
  pq_ds1 = pq_dict[1]
  pq_ds3 = pq_dict[3]

  # get the dataset name
  pq_ds1.dataset_name
  pq_ds3.dataset_name


  # get the path to the quantification result
  pq_ds1.quant_path
  pq_ds3.quant_path

  # get the AnnData
  pq_ds1.anndata
  pq_da3.anndata



Full usage
----------

.. code:: bash 

  usage: pyroe fetch-quant [-h] [--fetch_dir FETCH_DIR] [--force] [--delete_tar]
                           [--quiet]
                           dataset-ids [dataset-ids ...]

  positional arguments:
    dataset-ids           The ids of the datasets to fetch

  optional arguments:
    -h, --help            show this help message and exit
    --fetch_dir FETCH_DIR
                          The path to a directory for storing fetched datasets.
    --force               A flag indicates whether existing datasets will be redownloaded by force.
    --delete_tar          A flag indicates whether fetched tar files will be deleted.
    --quiet               A flag indicates whether help messaged should not be printed.

  1. 500 Human PBMCs, 3' LT v3.1, Chromium Controller
  2. 500 Human PBMCs, 3' LT v3.1, Chromium X
  3. 1k PBMCs from a Healthy Donor (v3 chemistry)
  4. 10k PBMCs from a Healthy Donor (v3 chemistry)
  5. 10k Human PBMCs, 3' v3.1, Chromium X
  6. 10k Human PBMCs, 3' v3.1, Chromium Controller
  7. 10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, Single Indexed
  8. 10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, Dual Indexed
  9. 20k Human PBMCs, 3' HT v3.1, Chromium X
  10. PBMCs from EDTA-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)
  11. PBMCs from Heparin-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)
  12. PBMCs from ACD-A Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)
  13. PBMCs from Citrate-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)
  14. PBMCs from Citrate-Treated Cell Preparation Tubes (3' v3.1 Chemistry)
  15. PBMCs from a Healthy Donor: Whole Transcriptome Analysis
  16. Whole Blood RBC Lysis for PBMCs and Neutrophils, Granulocytes, 3'
  17. Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 5)
  18. Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 1)
  19. Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Chromium Connect (channel 5)
  20. Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Chromium Connect (channel 1)
  21. Hodgkin's Lymphoma, Dissociated Tumor: Whole Transcriptome Analysis
  22. 200 Sorted Cells from Human Glioblastoma Multiforme, 3’ LT v3.1
  23. 750 Sorted Cells from Human Invasive Ductal Carcinoma, 3’ LT v3.1
  24. 2k Sorted Cells from Human Glioblastoma Multiforme, 3’ v3.1
  25. 7.5k Sorted Cells from Human Invasive Ductal Carcinoma, 3’ v3.1
  26. Human Glioblastoma Multiforme: 3’v3 Whole Transcriptome Analysis
  27. 1k Brain Cells from an E18 Mouse (v3 chemistry)
  28. 10k Brain Cells from an E18 Mouse (v3 chemistry)
  29. 1k Heart Cells from an E18 mouse (v3 chemistry)
  30. 10k Heart Cells from an E18 mouse (v3 chemistry)
  31. 10k Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells, Single Indexed
  32. 10k Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells, Dual Indexed
  33. 1k PBMCs from a Healthy Donor (v2 chemistry)
  34. 1k Brain Cells from an E18 Mouse (v2 chemistry)
  35. 1k Heart Cells from an E18 mouse (v2 chemistry)


The ProcessedQuant class
------------------------

To store the information of a dataset, we provide the ``ProcessedQuant`` class, which can be simply instantiated using a dataset id, for example, ``ProcessedQuant(2)`` will return an instance of the ``ProcessedQuant`` class containing the detail of dataset #2, 500 Human PBMCs, 3' LT v3.1, Chromium X. This class contains methods for fetching, decompressing and loading the quantification result of the corresponding dataset. After getting an instance of the class, i.e., running ``pq = ProcessedQuant(dataset_id)``, one can run the following commands to fetch, decompress and/or load the quantification result of the dataset:

* ``pq.fetch_quant()`` fetches the compressed quantification result of the corresponding dataset into a local directory and stores the path in its ``tar_path`` attribute.
* ``pq.decompress_quant()`` decompresses the fetched quantification result into a local directory and stores the path in its ``quant_path`` attribute.
* ``pq.load_quant()`` loads the decompressed quantification result into python as an `AnnData` object and stores the object in its ``anndata`` attribute.

Besides, we have some helper function for printing and loading the information of the available datasets:

* ``ProcessedQuant.get_available_dataset_df()`` returns the detail of available datasets as a pandas dataframe.
* ``ProcessedQuant.print_available_datasets()`` prints the index and name of the available datasets.
