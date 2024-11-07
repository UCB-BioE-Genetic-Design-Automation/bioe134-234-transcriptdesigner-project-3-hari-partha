import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

@pytest.fixture
def internal_rbs_checker():
    checker = InternalRBSChecker()
    return checker

def test_no_internal_rbs(internal_rbs_checker):
    dna_sequence = "ATGCGTACGTTAGCGATTTGCGGCCGTTGCCG"
    result, problematic_sequence = internal_rbs_checker.run(dna_sequence)
    print(f"Result: {result}, Sequence: {dna_sequence}")
    assert result is True
    assert problematic_sequence is None

def test_internal_rbs_found(internal_rbs_checker):
    dna_sequence = "AAAGGAGGTAGGGGTGATGAAA"
    result, problematic_sequence = internal_rbs_checker.run(dna_sequence)
    print(f"Result: {result}, Problematic Sequence: {problematic_sequence}")
    assert result is False
    assert problematic_sequence == "AGGAGGTGATG"

def test_internal_rbs_multiple_sites(internal_rbs_checker):
    dna_sequence = "ATGGGAGGCGTACGGTGTTAGGTGGAGGTAGTGATG"
    result, problematic_sequence = internal_rbs_checker.run(dna_sequence)
    print(f"Result: {result}, Problematic Sequence: {problematic_sequence}")
    assert result is False
    assert problematic_sequence == "GGAGGTAGTG"

def test_case_insensitive_detection(internal_rbs_checker):
    dna_sequence = "ggaggttGCCGTAGTGAAA"
    result, problematic_sequence = internal_rbs_checker.run(dna_sequence)
    print(f"Result: {result}, Problematic Sequence: {problematic_sequence}")
    assert result is False
    assert problematic_sequence == "GGAGTTGCC"

def test_edge_case_no_shine_dalgarno(internal_rbs_checker):
    dna_sequence = "TTTCCCGTGGGCACTGAGCACTG"
    result, problematic_sequence = internal_rbs_checker.run(dna_sequence)
    print(f"Result: {result}, Sequence: {dna_sequence}")
    assert result is True
    assert problematic_sequence is None
