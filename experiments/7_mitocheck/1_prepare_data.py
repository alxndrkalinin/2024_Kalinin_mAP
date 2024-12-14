import pathlib

import pandas as pd

from pycytominer import feature_select


def split_data(pycytominer_output: pd.DataFrame, dataset: str = "CP_and_DP"):
    all_cols = pycytominer_output.columns.tolist()

    # get DP,CP, or both features from all columns depending on desired dataset
    if dataset == "CP":
        feature_cols = [col for col in all_cols if "CP__" in col]
    elif dataset == "DP":
        feature_cols = [col for col in all_cols if "DP__" in col]
    elif dataset == "CP_and_DP":
        feature_cols = [col for col in all_cols if "P__" in col]

    # metadata columns is all columns except feature columns
    metadata_cols = [col for col in all_cols if "P__" not in col]

    metadata_dataframe = pycytominer_output[metadata_cols]
    feature_data = pycytominer_output[feature_cols].values

    return metadata_dataframe, feature_data


training_singlecell_data = pathlib.Path("data/raw/training_data__no_ic.csv.gz").resolve(
    strict=True
)
neg_control_data = pathlib.Path(
    "data/raw/normalized_data/negative_control_data.csv.gz"
).resolve(strict=True)

map_out_dir = pathlib.Path("data/processed/mAP_scores/")
map_out_dir.mkdir(parents=True, exist_ok=True)

training_sc_data = pd.read_csv(training_singlecell_data).drop("Unnamed: 0", axis=1)
neg_control_sc_data = pd.read_csv(neg_control_data)
neg_control_sc_data.insert(0, "Mitocheck_Phenotypic_Class", "neg_control")

training_sc_data.insert(1, "Metadata_is_control", 0)
neg_control_sc_data.insert(1, "Metadata_is_control", 1)

training_sc_data = training_sc_data.drop("Metadata_Object_Outline", axis=1)

print("control shape:", neg_control_sc_data.shape)
print("training shape:", training_sc_data.shape)

fs_ops = ["variance_threshold", "correlation_threshold", "drop_na_columns", "blocklist"]

cp_cols = [
    colname for colname in training_sc_data.columns if colname.startswith("CP__")
]

train_meta, train_features = split_data(training_sc_data, dataset="CP")
cp_data = pd.concat([train_meta, pd.DataFrame(train_features)], axis=1)
cp_data.columns = train_meta.columns.tolist() + cp_cols

pycytm_cp_training_feats_df = feature_select(cp_data, features=cp_cols)
pycytm_cp_training_feats_df = pycytm_cp_training_feats_df[
    [
        cols
        for cols in pycytm_cp_training_feats_df.columns.tolist()
        if cols.startswith("CP__")
    ]
]
del cp_data

training_sc_data = training_sc_data[
    [col for col in training_sc_data.columns.tolist() if not col.startswith("CP__")]
]
training_sc_data = pd.concat([training_sc_data, pycytm_cp_training_feats_df], axis=1)
print(f"{training_sc_data.shape=}")

cp_cols = [
    colname for colname in neg_control_sc_data.columns if colname.startswith("CP__")
]

neg_control_meta, neg_control_features = split_data(neg_control_sc_data, dataset="CP")
cp_data = pd.concat([neg_control_meta, pd.DataFrame(neg_control_features)], axis=1)
cp_data.columns = neg_control_meta.columns.tolist() + cp_cols

pycytm_cp_training_feats_df = feature_select(cp_data, features=cp_cols)
pycytm_cp_training_feats_df = pycytm_cp_training_feats_df[
    [
        cols
        for cols in pycytm_cp_training_feats_df.columns.tolist()
        if cols.startswith("CP__")
    ]
]
del cp_data

neg_control_sc_data = neg_control_sc_data[
    [col for col in neg_control_sc_data.columns.tolist() if not col.startswith("CP__")]
]
neg_control_sc_data = pd.concat(
    [neg_control_sc_data, pycytm_cp_training_feats_df], axis=1
)
print(f"{neg_control_sc_data.shape=}")

neg_control_sc_data_subset = neg_control_sc_data.sample(frac=0.005, random_state=0)
neg_control_sc_data_subset.to_parquet(
    "data/processed/neg_control_sc_fs_subset.parquet", index=False
)
print(f"{neg_control_sc_data_subset.shape=}")
training_sc_data.to_parquet("data/processed/training_sc_fs.parquet", index=False)
print(f"{training_sc_data.shape=}")
