from .pyroe_utils import say
import pandas as pd
import os
import shutil
import urllib.request
import tarfile
from .load_fry import load_fry


class ProcessedQuant:
    """
    A class stores the information of the quantification
    result of a processed dataset
    """

    def get_available_dataset_df():
        """
        get the dataframe in which each row contains
        the information of an available dataset that
        can be fetched.
        """

        # load available dataset sheet
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, "data", "available_datasets.tsv")
        available_datasets = pd.read_csv(my_file, sep="\t")

        return available_datasets

    def print_available_datasets():
        """
        Print the index and name of the available datasets.
        """
        available_datasets = ProcessedQuant.get_available_dataset_df()
        epilog = "\n".join(
            [
                "".join([f"{idx+1}", ". ", dataset_name])
                for (idx, dataset_name) in zip(
                    range(available_datasets.shape[0]),
                    available_datasets["dataset_name"].tolist(),
                )
            ]
        )
        epilog = "  \n".join(["Index of the available datasets:", epilog])
        print(epilog)

    def __init__(self, dataset_id: int):
        available_datasets = ProcessedQuant.get_available_dataset_df()
        if dataset_id < 0 or dataset_id >= available_datasets.shape[0]:
            raise ValueError(
                "Invalid dataset_id, run",
                "ProcessedQuant.print_available_datasets()",
                "to get available dataset ids.",
            )

        # get the info of the queried dataset id, python is zero based.
        available_dataset = available_datasets.iloc[dataset_id - 1, :]
        self.dataset_id = available_dataset["dataset_id"]
        self.chemistry = available_dataset["chemistry"]
        self.reference = available_dataset["reference"]
        self.dataset_name = available_dataset["dataset_name"]
        self.dataset_url = available_dataset["dataset_url"]
        self.fastq_url = available_dataset["fastq_url"]
        self.fastq_MD5sum = available_dataset["fastq_MD5sum"]
        self.delete_fastq = available_dataset["delete_fastq"]
        self.feature_barcode_csv_url = available_dataset["feature_barcode_csv_url"]
        self.multiplexing_library_csv_url = available_dataset[
            "multiplexing_library_csv_url"
        ]
        self.quant_tar_url = available_dataset["quant_tar_url"]
        self.quant_path = None
        self.tar_path = None
        self.anndata = None

    def fetch_quant(
        self, tar_dir="quant_tar", file_name=None, force=False, quiet=False
    ):
        """
        Fetch processed quantification to a local directory.\\
        The path to the fetched tar file will be sotred
        as the `ProcessedQuant.tar_path` attribute.

        Parameters
        ----------
        tar_dir: `str` (default: `quant_tar`)
            The directory for saving the fetched tar file.

        file_name: `str` (default: dataset id)
            Customized file name of the fetched tar file.
            Default is the dataset id.

        force: `bool` (default: `False`)
            If `True`, any existing tar file will be overwritten.

        quiet: `bool` (default: `False`)
            If `True`, help messaged will be printed out.
        """

        self.check_validity()

        say(quiet, f"Fetching the quant result of dataset #{self.dataset_id}")

        # check whether tar file exist,
        # download it if needed
        if self.tar_path is not None:
            if os.path.exists(self.tar_path) and (not force):
                say(
                    quiet,
                    "  - The tar_path attribute is not None and the path exists:",
                )
                say(quiet, f"    {self.tar_path}")
                say(quiet, "  - Pass force=True to fetch it again\n")
                return

        # folder for (temporarily) storing tar files.
        if not os.path.exists(tar_dir):
            os.makedirs(tar_dir)

        # process file_name
        if file_name is None:
            file_name = "".join([f"{self.dataset_id}", ".tar"])
        elif not file_name.endswith(".tar"):
            file_name = "".join([f"{file_name}", ".tar"])

        # update tar_path
        tar_path = os.path.join(tar_dir, file_name)

        if os.path.exists(tar_path):
            if force:
                say(quiet, "  - Overwriting the existing tar file:")
                say(quiet, f"    {tar_path}")
            else:
                say(quiet, "  - Use the existing file as tar_path:")
                say(quiet, f"    {tar_path}")
                say(quiet, "  - Pass force=True to overwrite it")
                self.tar_path = tar_path
                return

        # download tar file
        urllib.request.urlretrieve(self.quant_tar_url, tar_path)
        self.tar_path = tar_path
        say(quiet, "  - Fetched quant tar is saved as:")
        say(quiet, f"    {self.tar_path}")

    def decompress_quant(
        self,
        quant_dir="processed_quant",
        quant_path_name=None,
        force=False,
        quiet=False,
    ):
        """
        Decompress the fetched quantification to a local directory.\\
        The path to the decompressed quantification result will be sotred
        as the `ProcessedQuant.quant_path` attribute.
        Parameters
        ----------
        quant_dir: `str` (default: `processed_quant`)
            The directory for saving decompressed quantification result folder.

        quant_path_name: `str` (default: dataset id)
            Customized folder name of the quantification result folder.
            Default is the dataset id.

        force: `bool` (default: `False`)
            If `True`, existing tar file will be overwritten.

        quiet: `bool` (default: `False`)
            If `True`, help messaged will be printed out.
        """
        # make sure class is valid
        self.check_validity()

        # make sure tar file is valid
        if self.tar_path is None:
            raise ValueError(
                "tar_path attribute is None, run ProcessedQuant.fetch_quant() method to fetch the tar file."
            )

        say(
            quiet,
            f"Decompressing the quant result of dataset #{self.dataset_id} using:\n  {self.tar_path}",
        )

        # if quant_path is not None, return unless force=TRUE
        if self.quant_path is not None:
            if os.path.exists(self.tar_path) and (not force):
                say(
                    quiet,
                    "  - The quant_path attribute is not None and the path exists:",
                )
                say(quiet, f"    {self.quant_path}")
                say(quiet, "  - pass force=True to decompress it again")
                return

        # check expected output dir
        if quant_path_name is None:
            quant_path_name = self.dataset_id

        quant_parent_dir = os.path.join(quant_dir, f"{quant_path_name}")

        if os.path.exists(quant_parent_dir):
            if force:
                say(quiet, "  - Removing existing quant folder:")
                say(quiet, f"    {quant_parent_dir}")
                shutil.rmtree(quant_parent_dir)
            else:
                say(quiet, "  - Use the existing directory as quant_path:")
                say(quiet, f"    {quant_parent_dir}")
                say(quiet, "  - pass force=True to overwrite it")
                self.quant_path = os.path.join(
                    quant_parent_dir, next(os.walk(quant_parent_dir))[1][0]
                )
                return

        # decompress the tar file
        tf = tarfile.open(self.tar_path)
        tf.extractall(quant_parent_dir)
        self.quant_path = os.path.join(
            quant_parent_dir, next(os.walk(quant_parent_dir))[1][0]
        )
        say(quiet, "  - Decompressed quant result is saved as:")
        say(quiet, f"    {self.quant_path}")

    def load_quant(
        self, output_format="scRNA", force=False, nonzero=False, quiet=False
    ):
        """
        Load the quantification result as the `ProcessedQuant.anndata` attribute.\\

        Parameters
        ----------
        output_format: `str` or `dict` (default: `scRNA`)
            A string represents one of the pre-defined output formats, which are "scRNA", "snRNA" and "velocity". \\
            If a customized format of the returned `AnnData` is needed, one can pass a dictionary.\\
            See [load_fry](https://github.com/COMBINE-lab/pyroe/blob/main/src/pyroe/load_fry.py) for details.

        nonzero: `bool` (default: `False`)
            If `True`, the genes that have zero expression across all cells will be removed.

        quiet: `bool` (default: `False`)
            If `True`, help messaged will not be printed out.
        """
        self.check_validity()

        # make sure quant dir is valid
        if self.quant_path is None:
            raise ValueError(
                "The quant_path attribute is None, run ProcessedQuant.fetch_quant() and then ProcessedQuant.decompress_quant() to generate it."
            )

        if not os.path.exists(self.quant_path):
            raise ValueError(
                "The quant_path attribute is invalid, run ProcessedQuant.fetch_quant() and then ProcessedQuant.decompress_quant() to regenerate it."
            )

        if (self.anndata is not None) and (not force):
            say(quiet, "  - The anndata attribute is not None.")
            say(quiet, "  - pass force=True to update it")
            return

        say(quiet, f"Loading dataset #{self.dataset_id} from:")
        say(quiet, f"  {self.quant_path}")

        self.anndata = load_fry(
            frydir=self.quant_path,
            output_format=output_format,
            nonzero=nonzero,
            quiet=quiet,
        )

    def FDL(
        dataset_id: int,
        tar_dir="quant_tar",
        tar_file_name=None,
        quant_dir="processed_quant",
        quant_path_name=None,
        output_format="scRNA",
        nonzero=False,
        force=False,
        quiet=False,
    ):
        """
        Call `ProcessedQuant.fetch_quant()`, ProcessedQuant.decompress_quant() and ProcessedQuant.load_quant() in turn
        for a dataset to generate a complete ProcessedQuant object.

        Parameters
        -----------------------
        dataset_id: `int`
            The id of an available dataset

        tar_dir: `str` (default: `quant_tar`)
            The directory for saving the fetched tar file.

        tar_file_name: `str` (default: dataset id)
            Customized file name of the fetched tar file.
            Default is the dataset id.

        quant_dir: `str` (default: `processed_quant`)
            The directory for saving decompressed quantification result folder.

        quant_path_name: `str` (default: dataset id)
            Customized folder name of the quantification result folder.
            Default is the dataset id.
        output_format: `str` or `dict` (default: `scRNA`)
            A string represents one of the pre-defined output formats, which are "scRNA", "snRNA" and "velocity". \\
            If a customized format of the returned `AnnData` is needed, one can pass a Dictionary.\\
            See [load_fry](https://github.com/COMBINE-lab/pyroe/blob/main/src/pyroe/load_fry.py) for details.

        nonzero: `bool` (default: `False`)
            If `True`, existing tar file will be overwritten.

        force: `bool` (default: `False`)
            If `True`, existing tar file will be overwritten.

        quiet: `bool` (default: `False`)
            If `True`, help messaged will be printed out.
        """
        processed_quant = ProcessedQuant(dataset_id)

        # fetch it
        processed_quant.fetch_quant(
            tar_dir=tar_dir, file_name=tar_file_name, force=force, quiet=quiet
        )

        # decompress it
        processed_quant.decompress_quant(
            quant_dir=quant_dir,
            quant_path_name=quant_path_name,
            force=force,
            quiet=quiet,
        )

        # load it
        processed_quant.load_quant(
            output_format=output_format, force=force, nonzero=nonzero, quiet=quiet
        )

        return processed_quant

    def check_validity(self):
        if (
            self.quant_tar_url is None
            or self.dataset_id is None
            or self.chemistry is None
            or self.reference is None
            or self.dataset_name is None
            or self.dataset_url is None
            or self.fastq_url is None
            or self.fastq_MD5sum is None
            or self.delete_fastq is None
            or self.feature_barcode_csv_url is None
            or self.multiplexing_library_csv_url is None
            or self.quant_tar_url is None
        ):
            raise ValueError(
                "Incomplete class object, use",
                "ProcessedQuant(dataset_id)",
                "to instantiate it.",
            )
