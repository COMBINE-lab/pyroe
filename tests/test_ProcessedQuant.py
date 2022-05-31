# This file contains the testing functions for ProcessedQuant class methods.
# Here we will use dataset #22, which has only 200 expected cells
import tempfile
import os
import pandas as pd
from pyroe.ProcessedQuant import ProcessedQuant
import scanpy


# test class instantiation
def test_init():
    ds_22 = ProcessedQuant.get_available_dataset_df().iloc[21, :]
    ds = ProcessedQuant(22)
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
    assert ds.quant_path is None
    assert ds.tar_path is None
    assert ds.anndata is None


# test fetch_quant and its parameters
def test_fetch_quant():
    output_dir = tempfile.TemporaryDirectory()
    file_name = "test_ds22_quant_tar"
    ds = ProcessedQuant(22)
    ds.fetch_quant(tar_dir=output_dir.name, file_name=file_name)
    expected_tar_path = os.path.join(output_dir.name, f"{file_name}.tar")

    # check fetched file existence
    assert os.path.exists(expected_tar_path)

    # check whether force works
    ds.fetch_quant(tar_dir=output_dir.name, file_name=file_name, force=True)
    assert os.path.exists(expected_tar_path)


def test_decompress_quant():
    # setup
    output_dir = tempfile.TemporaryDirectory()
    tar_dir = os.path.join(output_dir.name, "quant_tar")
    file_name = "test_ds22_quant_tar"
    ds = ProcessedQuant(22)
    ds.fetch_quant(tar_dir=tar_dir, file_name=file_name)

    quant_dir = os.path.join(output_dir.name, "processed_quant")
    quant_path_name = "test_ds22_quant_dir"
    expected_quant_path = os.path.join(
        quant_dir,
        quant_path_name,
        "_".join([ds.fastq_MD5sum, "fry_unfilt_quant_usa_cr-like"]),
    )

    # check whether decompress_quant works
    ds.decompress_quant(quant_dir=quant_dir, quant_path_name=quant_path_name)
    assert ds.quant_path == expected_quant_path
    assert os.path.exists(expected_quant_path)

    # check force
    ds.decompress_quant(
        quant_dir=quant_dir, quant_path_name=quant_path_name, force=True
    )
    assert os.path.exists(expected_quant_path)


def test_load_quant():
    # setup
    output_dir = tempfile.TemporaryDirectory()
    tar_dir = os.path.join(output_dir.name, "quant_tar")
    file_name = "test_ds22_quant_tar"
    ds = ProcessedQuant(22)
    ds.fetch_quant(tar_dir=tar_dir, file_name=file_name)
    quant_dir = os.path.join(output_dir.name, "processed_quant")
    quant_path_name = "test_ds22_quant_dir"
    ds.decompress_quant(quant_dir=quant_dir, quant_path_name=quant_path_name)

    # check whether load_quant works
    ds.load_quant()
    assert isinstance(ds.anndata, scanpy.AnnData)

    # check force
    ds.load_quant(force=True)
    assert isinstance(ds.anndata, scanpy.AnnData)

    # check another output_format

    ds.load_quant(output_format="raw", force=True)
    observed_layer_names = sorted(list(ds.anndata.layers.keys()))
    expected_layer_names = ["ambiguous", "spliced", "unspliced"]
    assert observed_layer_names == expected_layer_names

    # check nonzero
    ds.load_quant(output_format="raw", force=True, nonzero=True)
    assert isinstance(ds.anndata, scanpy.AnnData)
