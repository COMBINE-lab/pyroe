from .pyroe_utils import say

def fetch_processed_quant(
    dataset_ids = [],
    fetch_dir = "processed_quant",
    force = False,
    delete_tar = True,
    quiet = False
):
    """
    Download the quantification result of the preprocessed 10x datasets.

    Required Parameters
    ----------
    dataset_ids : `int` or `list`
        The list of the id of some available datasets.

    Optional Parameters
    ----------
    fetch_dir : `str` (default: `processed_quant`)
        The path to a directory for storing fetched datasets.
    
    force : `bool` (default: `False`)
        True if existing datasets should be re-downloaded.
        
    delete_tar : `bool` (default: `True`)
        True if intermediate tar files should be deleted.
        If False, they will be stored in the datasets_tar
        folder under the fetch_dir.
        
    quiet : `bool` (default: `False`)
        True if function should be quiet.
        False if messages (including error messages) should be printed out. 

    Returns
    -------
    If an empty dataset_ids list is given, a dataframe 
    containing the information of all available datasets
    will be returned. If one or more dataset ids are 
    provided as dataset_ids, a dictionary of str paths 
    will be returned, in which the keys are the dataset ids,
    and the values are the path to the folder of the 
    corresponding quantification result.

    Notes
    -----
    10x Genomics provides many publicly available single-cell
    RNA-sequencing experiments on their 
    [website](https://www.10xgenomics.com/resources/datasets).
    To avoid reinventing wheels, we processed these datasets
    using a nextflow-based 
    [alevin-fry workflow](https://github.com/COMBINE-lab/10x-requant) 
    and made the quantification results available for free downloading. 
    Currently, the available datasets include (Notice that dataset id starts form one, not zero):
    
    1. [500 Human PBMCs, 3' LT v3.1, Chromium Controller](https://www.10xgenomics.com/resources/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-controller-3-1-low-6-1-0): [link to the quant result](https://umd.box.com/shared/static/tg919re5gd4klua39z3zemcg9ya422am.tar)
    1. [500 Human PBMCs, 3' LT v3.1, Chromium X](https://www.10xgenomics.com/resources/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-x-3-1-low-6-1-0): [link to the quant result](https://umd.box.com/shared/static/lrl68q2lz0ltsvs89iazbr302p50wnqj.tar)
    1. [1k PBMCs from a Healthy Donor (v3 chemistry)](https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/wrn19wsmkem1jyc9seqpe4pxto5zimwa.tar)
    1. [10k PBMCs from a Healthy Donor (v3 chemistry)](https://www.10xgenomics.com/resources/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/01j9574g1yd93noz2pqlqjfrdhx0m1ff.tar)
    1. [10k Human PBMCs, 3' v3.1, Chromium X](https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-ht-v3-1-chromium-x-3-1-high): [link to the quant result](https://umd.box.com/shared/static/jvvzacmo98vxfnoimg4dgi52lifhl2aa.tar)
    1. [10k Human PBMCs, 3' v3.1, Chromium Controller](https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high): [link to the quant result](https://umd.box.com/shared/static/5dzu2tw8nz9tijt8lgmelll6sbaaomh4.tar)
    1. [10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, Single Indexed](https://www.10xgenomics.com/resources/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-single-indexed-3-1-standard-4-0-0): [link to the quant result](https://umd.box.com/shared/static/iol9bxiv740xq6m29p2fzcoe8volsi7i.tar)
    1. [10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, Dual Indexed](https://www.10xgenomics.com/resources/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-dual-indexed-3-1-standard-4-0-0): [link to the quant result](https://umd.box.com/shared/static/5dzu2tw8nz9tijt8lgmelll6sbaaomh4.tar)
    1. [20k Human PBMCs, 3' HT v3.1, Chromium X](https://www.10xgenomics.com/resources/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0): [link to the quant result](https://umd.box.com/shared/static/c609sk8w6cbn4w0tcwofz4qcyjp67506.tar)
    1. [PBMCs from EDTA-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)](https://www.10xgenomics.com/resources/datasets/pbmcs-3p_edta_sepmate-3-1-standard): [link to the quant result](https://umd.box.com/shared/static/imedrs558dx4tzxy9uhhxvy0dmjlhjsh.tar)
    1. [PBMCs from Heparin-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)](https://www.10xgenomics.com/resources/datasets/pbmcs-3p_heparin_sepmate-3-1-standard): [link to the quant result](https://umd.box.com/shared/static/e8gqxali0lwy2nashh5rmmoc6bgj92xm.tar)
    1. [PBMCs from ACD-A Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)](https://www.10xgenomics.com/resources/datasets/pbmcs-3p_acda_sepmate-3-1-standard): [link to the quant result](https://umd.box.com/shared/static/w1kdz3vifqi4ixtqkuwqgc2mpkkiehky.tar)
    1. [PBMCs from Citrate-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)](https://www.10xgenomics.com/resources/datasets/pbmcs-3p_citrate_sepmate-3-1-standard): [link to the quant result](https://umd.box.com/shared/static/cs0s6e2u0j7d8uc36xsdo6922c7dle6y.tar)
    1. [PBMCs from Citrate-Treated Cell Preparation Tubes (3' v3.1 Chemistry)](https://www.10xgenomics.com/resources/datasets/pbmcs-3p_citrate_cpt-3-1-standard): [link to the quant result](https://umd.box.com/shared/static/2tqrzreghvi6nxe94oob1ei1vi4458br.tar)
    1. [PBMCs from a Healthy Donor: Whole Transcriptome Analysis](https://www.10xgenomics.com/resources/datasets/pbm-cs-from-a-healthy-donor-whole-transcriptome-analysis-3-1-standard-4-0-0): [link to the quant result](https://umd.box.com/shared/static/dk0hmj5mpqjq56afkr5jibavy9e3yil8.tar)
    1. [Whole Blood RBC Lysis for PBMCs and Neutrophils, Granulocytes, 3'](https://www.10xgenomics.com/resources/datasets/whole-blood-rbc-lysis-for-pbmcs-neutrophils-granulocytes-3-3-1-standard): [link to the quant result](https://umd.box.com/shared/static/0gnwx7d9hbdmptyi0ddz6mfa79d1l8be.tar)
    1. [Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 5)](https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-manual-channel-5-3-1-standard-3-1-0): [link to the quant result](https://umd.box.com/shared/static/tn884ctombnj214abt8rp77p7kih5i02.tar)
    1. [Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 1)](https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-manual-channel-1-3-1-standard-3-1-0): [link to the quant result](https://umd.box.com/shared/static/0jcgdgy8woj30oarkwhybk8fly7gb7v8.tar)
    1. [Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Chromium Connect (channel 5)](https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-chromium-connect-channel-5-3-1-standard-3-1-0): [link to the quant result](https://umd.box.com/shared/static/kybks0ncf609xhcwvhv7z743zrmvlg94.tar)
    1. [Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Chromium Connect (channel 1)](https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-chromium-connect-channel-1-3-1-standard-3-1-0): [link to the quant result](https://umd.box.com/shared/static/vtuexhbqiyvfob7qdpvsxl1nbqlo074f.tar)
    1. [Hodgkin's Lymphoma, Dissociated Tumor: Whole Transcriptome Analysis](https://www.10xgenomics.com/resources/datasets/hodgkins-lymphoma-dissociated-tumor-whole-transcriptome-analysis-3-1-standard-4-0-0): [link to the quant result](https://umd.box.com/shared/static/qis4ovf34wvq12n2uabdiem6w355qry7.tar)
    1. [200 Sorted Cells from Human Glioblastoma Multiforme, 3’ LT v3.1](https://www.10xgenomics.com/resources/datasets/200-sorted-cells-from-human-glioblastoma-multiforme-3-lt-v-3-1-3-1-low-6-0-0): [link to the quant result](https://umd.box.com/shared/static/2xf9xf8m1n5vbvmpo1vshwigs7f7o5jd.tar)
    1. [750 Sorted Cells from Human Invasive Ductal Carcinoma, 3’ LT v3.1](https://www.10xgenomics.com/resources/datasets/750-sorted-cells-from-human-invasive-ductal-carcinoma-3-lt-v-3-1-3-1-low-6-0-0): [link to the quant result](https://umd.box.com/shared/static/3txnreehxoj2plyypfs6fkibnnbo72h4.tar)
    1. [2k Sorted Cells from Human Glioblastoma Multiforme, 3’ v3.1](https://www.10xgenomics.com/resources/datasets/2-k-sorted-cells-from-human-glioblastoma-multiforme-3-v-3-1-3-1-standard-6-0-0): [link to the quant result](https://umd.box.com/shared/static/n0vpgbdwbnnqdw1h9of2ykk7ive9p6pt.tar)
    1. [7.5k Sorted Cells from Human Invasive Ductal Carcinoma, 3’ v3.1](https://www.10xgenomics.com/resources/datasets/7-5-k-sorted-cells-from-human-invasive-ductal-carcinoma-3-v-3-1-3-1-standard-6-0-0): [link to the quant result](https://umd.box.com/shared/static/aly78r6bppqf01npbqfopc3epmp17weu.tar)
    1. [Human Glioblastoma Multiforme: 3’v3 Whole Transcriptome Analysis](https://www.10xgenomics.com/resources/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0): [link to the quant result](https://umd.box.com/shared/static/suf8pt3avv4rchxfw0bqrshslzieygef.tar)
    1. [1k Brain Cells from an E18 Mouse (v3 chemistry)](https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/4w5eiq3qafbru5ocler39j5j28bvgz98.tar)
    1. [10k Brain Cells from an E18 Mouse (v3 chemistry)](https://www.10xgenomics.com/resources/datasets/10-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/tym9m73frtp13vo15jhit9uwuk3mtfdq.tar)
    1. [1k Heart Cells from an E18 mouse (v3 chemistry)](https://www.10xgenomics.com/resources/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/d838oy3udjvtzjo7tsdiao7u6sazabeg.tar)
    1. [10k Heart Cells from an E18 mouse (v3 chemistry)](https://www.10xgenomics.com/resources/datasets/10-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/efinlf6p8weich13kv3bzrlndsx963v4.tar)
    1. [10k Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells, Single Indexed](https://www.10xgenomics.com/resources/datasets/10-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-single-indexed-3-1-standard-4-0-0): [link to the quant result](https://umd.box.com/shared/static/mr0yolo83rjdcdqgu6om4q133fpime8r.tar)
    1. [10k Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells, Dual Indexed](https://www.10xgenomics.com/resources/datasets/10-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-dual-indexed-3-1-standard-4-0-0): [link to the quant result](https://umd.box.com/shared/static/mr7raea3v5ccn4dchemwhcimpz7t1cwl.tar)
    1. [1k PBMCs from a Healthy Donor (v2 chemistry)](https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-2-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/xeya5zr0t0wg0t8c20zu0pdhclxywx3c.tar)
    1. [1k Brain Cells from an E18 Mouse (v2 chemistry)](https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-v-2-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/a53twm69uo2xf6778asuvw2aft7wkur5.tar)
    1. [1k Heart Cells from an E18 mouse (v2 chemistry)](https://www.10xgenomics.com/resources/datasets/1-k-heart-cells-from-an-e-18-mouse-v-2-chemistry-3-standard-3-0-0): [link to the quant result](https://umd.box.com/shared/static/p4ieuzimfgrjfsr9rzhrn48kved4ha7m.tar)

    To obtain the information of the available datasets as 
    a dataframe, one can run `load_processed_quant()`
    """
    
    import pandas as pd
    import os
    import shutil
    import urllib.request 
    import tarfile

    # load available dataset sheet
    location = os.path.dirname(os.path.realpath(__file__))
    my_file = os.path.join(location, 'data', 'available_datasets.tsv')
    # # my_file = os.path.join('data', 'available_datasets.tsv')
    available_datasets = pd.read_csv(my_file, sep="\t")

    # if no dataset is provided, just return the available dataset dataframe
    if len(dataset_ids) == 0:
        return available_datasets

    # check the validity of dataset_ids
    n_ds = available_datasets.shape[0]
    invalid_ids = []
    for idx, dataset_id in enumerate(dataset_ids):
        if (type(dataset_id) == int):
            if dataset_id > n_ds & dataset_id < 1:
                print(f"Found invalid dataset id '{dataset_id}', ignored.")
                invalid_ids.append(idx)
        else:
            print(f"Found invalid dataset id '{dataset_id}', ignored.")
            invalid_ids.append(idx)
    for i in reversed(invalid_ids):
        del dataset_ids[i]

    # if no id left, return an error
    if not dataset_ids:
        raise ValueError(f"No valid dataset id found, can not proceed")

    # download the quantification tar file for each queried dataset.
    quant_dir_list = []
    # folder for (temporarily) storing tar files.
    tar_dir = os.path.join(fetch_dir, "datasets_tar")
    if not os.path.exists(tar_dir):
        os.makedirs(tar_dir)

    # download the quantification tar file for each queried dataset.
    for dataset_id in dataset_ids:
        say(quiet, f"Processing dataset #{dataset_id}")
        quant_parent_dir = os.path.join(fetch_dir, 
                                        f"{dataset_id}")
        # python index is zero-based
        dataset_id -= 1
        tar_file = os.path.join(tar_dir,
                            "".join([available_datasets.iloc[dataset_id,5], ".tar"]))

        if os.path.exists(quant_parent_dir):
            say(quiet, f"  - output dir exists:\n    {quant_parent_dir}")
            
            if force:
                say(quiet, "  - force re-processing")
                shutil.rmtree(quant_parent_dir)
            else:
                say(quiet, "  - use the existing quant result")
                quant_dir_list.append(os.path.join(quant_parent_dir, next(os.walk(quant_parent_dir))[1][0]))
                continue

        say(quiet, "  - Downloading quant result")
        url = available_datasets.iloc[dataset_id, 9]
        urllib.request.urlretrieve(url, tar_file)
        # decompress the downloaded tar files
        say(quiet, "  - Decompressing quant result")
        tf = tarfile.open(tar_file)
        tf.extractall(quant_parent_dir)
        quant_dir_list.append(os.path.join(quant_parent_dir, next(os.walk(quant_parent_dir))[1][0]))

    # delete tar if needed
    if delete_tar:
        say(quiet, "Removing downloaded tar files")
    shutil.rmtree(tar_dir)

    say(quiet, "Done")
    return dict(zip(dataset_ids, quant_dir_list))
