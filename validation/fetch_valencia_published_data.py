#!/usr/bin/env python3
"""Stage VALENCIA's published dataset locally for the P2 validation gate.

Downloads the ~8 MB published composition table (with the paper's Val_CST/Val_subCST
labels) into .validation_data/ so test_cst_validation_gate.py can run. The data is
VALENCIA's (ravel-lab/VALENCIA, MIT) and is intentionally NOT committed to this repo.

    python validation/fetch_valencia_published_data.py
"""

from pathlib import Path
from urllib.request import urlretrieve

REPO_ROOT = Path(__file__).resolve().parents[1]
DEST = REPO_ROOT / ".validation_data"
RAW = "https://raw.githubusercontent.com/ravel-lab/VALENCIA/master"
FILES = {
    "all_samples_composition.csv":
        f"{RAW}/Publication_materials/Data_and_metadata/all_samples_taxonomic_composition_data.csv",
    "cst_centroids_012920.csv": f"{RAW}/CST_centroids_012920.csv",
}


def main() -> None:
    DEST.mkdir(exist_ok=True)
    for name, url in FILES.items():
        out = DEST / name
        if out.exists():
            print(f"have {name}")
            continue
        print(f"downloading {name} ...")
        urlretrieve(url, out)
    print(f"staged in {DEST}")


if __name__ == "__main__":
    main()
