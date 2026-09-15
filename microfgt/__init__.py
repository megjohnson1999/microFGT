"""microFGT — a one-stop tool for female genital tract (FGT) microbiome analysis.

Layer 1 (the integrated object) is a MuData; the import layer (this package's `io`
module) turns real tool outputs into it. See `microfgt/README.md` for the build phases.
"""

__version__ = "0.0.1"

from microfgt.characterize import describe_composition
from microfgt.cst import classify_cst
from microfgt.io import (
    build_mudata,
    collapse_to_taxon,
    existing_cst,
    import_phyloseq,
    import_speciateit,
    import_valencia,
    import_virgo,
)

__all__ = [
    "__version__",
    "build_mudata",
    "classify_cst",
    "collapse_to_taxon",
    "describe_composition",
    "existing_cst",
    "import_phyloseq",
    "import_speciateit",
    "import_valencia",
    "import_virgo",
]
