"""Analysis layer — two tiers over the integrated object.

* **Primitives** (``transforms`` / ``diversity`` / ``diffabund``): validated commodity stats
  bought from scikit-bio, each mutating one modality in place (a layer, an ``obs`` column, a
  distance matrix). The compute wheels, not reinvented.
* **Hypothesis-test verbs** (``compare_alpha`` / ``compare_beta`` / …): the tests a researcher
  actually runs — "is this community feature associated with the variable I picked, adjusting
  for these covariates?" Each takes the object + which ``obs`` variables play which role + a
  subset, and returns a uniform :class:`~microfgt.analysis.results.AnalysisResult` (tidy table
  + headline stats + a declarative plot spec). That uniform return is the contract Layer 3
  (the dashboard) binds to, so the surface stays a thin presenter of the same calls a power
  user makes.
"""

from microfgt.analysis.abundance import compare_abundance
from microfgt.analysis.association import associate
from microfgt.analysis.diffabund import differential_abundance
from microfgt.analysis.dispatch import VERBS, run_verb
from microfgt.analysis.diversity import alpha_diversity, beta_diversity, ordinate
from microfgt.analysis.hypothesis import compare_alpha, compare_beta
from microfgt.analysis.results import AnalysisResult
from microfgt.analysis.transforms import clr_transform, relative_abundance

__all__ = [
    "relative_abundance",
    "clr_transform",
    "alpha_diversity",
    "beta_diversity",
    "ordinate",
    "differential_abundance",
    # hypothesis-test verbs (return AnalysisResult)
    "AnalysisResult",
    "compare_alpha",
    "compare_beta",
    "associate",
    "compare_abundance",
    "run_verb",
    "VERBS",
]
