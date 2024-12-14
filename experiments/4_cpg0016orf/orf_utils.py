from pathlib import Path

import pandas as pd


MORPHMAP_MAP_URL = "https://github.com/jump-cellpainting/2024_Chandrasekaran_Morphmap/raw/5e6da3cae56443584a89cd9205a2e4190d8ebbaf/03.retrieve-annotations/output/phenotypic-activity-consistency.parquet"


def load_orf_map_results(data_dir="data", file_name="orf_map_results.csv"):
    file_path = Path(data_dir) / file_name
    if not file_path.exists():
        df = pd.read_parquet(MORPHMAP_MAP_URL)
        df = df.query("Modality == 'ORF'").reset_index(drop=True)
        df.to_csv(file_path, index=False)
    else:
        df = pd.read_csv(file_path)
    return df
