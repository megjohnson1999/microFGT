"""discover_pairs handles both Illumina _R1/_R2 and ENA/SRA _1/_2 naming."""

from microfgt.orchestrate.cutadapt import discover_pairs


def _touch(d, *names):
    for n in names:
        (d / n).write_text("@r\nACGT\n+\nIIII\n")


def test_illumina_r1_r2(tmp_path):
    _touch(tmp_path, "sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz")
    pairs = discover_pairs(tmp_path)
    assert [(s, r1.name, r2.name) for s, r1, r2 in pairs] == [
        ("sampleA", "sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz")
    ]


def test_illumina_lane_suffix(tmp_path):
    # The classic _R1_001 suffix must still resolve.
    _touch(tmp_path, "S1_L001_R1_001.fastq.gz", "S1_L001_R2_001.fastq.gz")
    (sample, r1, r2), = discover_pairs(tmp_path)
    assert sample == "S1_L001" and r2 is not None and r2.name == "S1_L001_R2_001.fastq.gz"


def test_ena_underscore_1_2(tmp_path):
    _touch(tmp_path, "ERR4421550_1.fastq.gz", "ERR4421550_2.fastq.gz")
    (sample, r1, r2), = discover_pairs(tmp_path)
    assert sample == "ERR4421550"
    assert r1.name == "ERR4421550_1.fastq.gz" and r2.name == "ERR4421550_2.fastq.gz"


def test_ena_sample_name_containing_underscore_one(tmp_path):
    # A sample whose own name ends in _1 must not be mangled — only the final marker counts.
    _touch(tmp_path, "cohort_1_1.fastq.gz", "cohort_1_2.fastq.gz")
    (sample, r1, r2), = discover_pairs(tmp_path)
    assert sample == "cohort_1"
    assert r2 is not None and r2.name == "cohort_1_2.fastq.gz"


def test_r1_preferred_over_underscore_one_in_mixed_dir(tmp_path):
    # If both conventions are present, _R1/_R2 wins and _1/_2 is not also emitted.
    _touch(tmp_path, "a_R1.fastq.gz", "a_R2.fastq.gz", "b_1.fastq.gz", "b_2.fastq.gz")
    samples = [s for s, _, _ in discover_pairs(tmp_path)]
    assert samples == ["a"]


def test_missing_mate_reports_none(tmp_path):
    _touch(tmp_path, "ERR999_1.fastq.gz")   # no _2
    (sample, r1, r2), = discover_pairs(tmp_path)
    assert sample == "ERR999" and r2 is None
