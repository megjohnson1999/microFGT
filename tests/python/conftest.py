"""Shared fixture paths for the importer tests.

Real tool-output fixtures live in ``prototype/real_fixtures/`` (genuine VIRGO/VALENCIA
outputs + speciateIT inputs). Synthetic speciateIT classification fixtures (no genuine
output ships) live next to these tests under ``data/``.
"""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REAL_FIXTURES = REPO_ROOT / "prototype" / "real_fixtures"
TEST_DATA = Path(__file__).resolve().parent / "data"


@pytest.fixture
def real_fixtures() -> Path:
    return REAL_FIXTURES


@pytest.fixture
def test_data() -> Path:
    return TEST_DATA
