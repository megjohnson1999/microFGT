"""P2 validation gate: centroid CST must reproduce VALENCIA on the published data.

Target (constraint B): >=99.9% subCST agreement vs the paper's own labels on the 13k
published samples. This needs VALENCIA's large published dataset
(``all_samples_taxonomic_composition_data.csv``, ~8 MB), which is NOT committed. Stage it
under ``.validation_data/`` (the helper script in ``validation/`` downloads it); the test
skips cleanly when it's absent so CI stays light, but the gate is fully reproducible.
"""

from pathlib import Path

import pandas as pd
import pytest

from microfgt.cst import classify_centroid

REPO_ROOT = Path(__file__).resolve().parents[2]
COMPOSITION = REPO_ROOT / ".validation_data" / "all_samples_composition.csv"

# Non-taxon leading columns in the published composition file.
_META_COLS = [
    "Sample_number_for_SRA", "Subject_number",
    "HC_CST", "HC_subCST", "Val_CST", "Val_subCST", "total_reads",
]


@pytest.mark.skipif(not COMPOSITION.exists(), reason="published dataset not staged (see validation/)")
def test_centroid_reproduces_valencia_on_published_data():
    df = pd.read_csv(COMPOSITION)
    taxa = [c for c in df.columns if c not in _META_COLS]

    counts = df[taxa]
    counts.index = df["Sample_number_for_SRA"].astype(str)
    rc = pd.Series(df["total_reads"].to_numpy(), index=counts.index)

    # Pin the 2020 centroids: this gate validates our port against genuine VALENCIA output on the
    # published dataset, which is 2020-named (see design/method_log.md M4). The default is 2024.
    out = classify_centroid(counts, read_count=rc, reference="2020")

    truth_subcst = df["Val_subCST"].astype(str).to_numpy()
    truth_cst = df["Val_CST"].astype(str).to_numpy()
    subcst_agree = (out["subCST"].to_numpy() == truth_subcst).mean()
    cst_agree = (out["CST"].to_numpy() == truth_cst).mean()

    # Report exact numbers regardless of pass/fail.
    print(f"\nsubCST agreement vs Val_subCST: {subcst_agree:.5f} ({len(df)} samples)")
    print(f"CST    agreement vs Val_CST:    {cst_agree:.5f}")

    assert subcst_agree >= 0.999, f"subCST agreement {subcst_agree:.5f} < 0.999"
