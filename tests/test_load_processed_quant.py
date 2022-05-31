# This file contains the testing functions for ProcessedQuant class methods.
# Here we will use dataset #22, which has only 200 expected cells
import tempfile
import os
import pandas as pd
from pyroe import ProcessedQuant, load_processed_quant
import scanpy


# test fetch_quant and its parameters
def test_fetch_quant():
    ds_22 = ProcessedQuant.get_available_dataset_df().iloc[21, :]
    output_dir = tempfile.TemporaryDirectory()
    ds = load_processed_quant([22], fetch_dir=output_dir.name).get(22)
    tar_dir = os.path.join(output_dir.name, "quant_tar")
    expected_tar_path = os.path.join(tar_dir, f"{ds.dataset_id}.tar")
    quant_dir = os.path.join(output_dir.name)
    expected_quant_path = os.path.join(
        quant_dir,
        str(ds.dataset_id),
        "_".join([ds.fastq_MD5sum, "fry_unfilt_quant_usa_cr-like"]),
    )

    # check fetched file existence
    assert ds.dataset_id == ds_22["dataset_id"]
    assert ds.chemistry == ds_22["chemistry"]
    assert ds.reference == ds_22["reference"]
    assert ds.dataset_name == ds_22["dataset_name"]
    assert ds.dataset_url == ds_22["dataset_url"]
    assert ds.fastq_url == ds_22["fastq_url"]
    assert ds.fastq_MD5sum == ds_22["fastq_MD5sum"]
    assert ds.delete_fastq == ds_22["delete_fastq"]
    assert pd.isna(ds.feature_barcode_csv_url)
    assert pd.isna(ds.multiplexing_library_csv_url)
    assert ds.quant_tar_url == ds_22["quant_tar_url"]
    assert os.path.exists(expected_tar_path)
    assert os.path.exists(expected_quant_path)
    assert isinstance(ds.anndata, scanpy.AnnData)
