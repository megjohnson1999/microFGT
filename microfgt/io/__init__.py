"""Import layer — real tool outputs in, integrated object out.

Each importer does only the FGT-specific reshape (the "glue" the user should never
have to write); the container is a generic ``anndata.AnnData`` / ``mudata.MuData``.
Formats are grounded in ``prototype/real_fixtures/FORMATS.md`` and the tools' own code,
NOT guessed.
"""

from microfgt.io.integrate import build_mudata
from microfgt.io.phyloseq import existing_cst, import_phyloseq
from microfgt.io.speciateit import collapse_to_taxon, import_speciateit
from microfgt.io.valencia import import_valencia
from microfgt.io.virgo import (
    collapse_virgo2_to_taxon,
    import_virgo,
    import_virgo2,
    import_virgo2_taxonomy,
)
from microfgt.io.vista import import_mgcst

__all__ = [
    "build_mudata",
    "collapse_to_taxon",
    "collapse_virgo2_to_taxon",
    "existing_cst",
    "import_mgcst",
    "import_phyloseq",
    "import_speciateit",
    "import_valencia",
    "import_virgo",
    "import_virgo2",
    "import_virgo2_taxonomy",
]
