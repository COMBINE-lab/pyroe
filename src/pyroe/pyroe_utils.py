def say(quiet, words):
    if not quiet:
        print(words)


def output_formats():
    return set(["h5ad", "loom", "csvs", "zarr"])


def check_dataset_ids(n_ds, dataset_ids):
    # check the validity of dataset_ids
    invalid_ids = []
    for idx, dataset_id in enumerate(dataset_ids):
        if type(dataset_id) == int:
            if dataset_id > n_ds or dataset_id < 1:
                print(f"Found invalid dataset id '{dataset_id}', ignored.")
                invalid_ids.append(idx)
        else:
            print(f"Found invalid dataset id '{dataset_id}', ignored.")
            invalid_ids.append(idx)

    if invalid_ids:
        for i in reversed(invalid_ids):
            del dataset_ids[i]

    return dataset_ids
