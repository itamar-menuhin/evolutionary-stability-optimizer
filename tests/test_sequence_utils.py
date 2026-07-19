from eso.sequence_utils import reverse_complement_seq, parse_region


def test_reverse_complement_seq():
    assert reverse_complement_seq("ACGT") == "ACGT"
    assert reverse_complement_seq("AATTCCGG") == "CCGGAATT"


def test_parse_region_empty():
    assert parse_region('') == ()
    assert parse_region('None') == ()


def test_parse_region_single():
    assert parse_region('1-9') == [(0, 9)]


def test_parse_region_multiple():
    assert parse_region('1-9, 21-29') == [(0, 9), (20, 29)]


def test_parse_region_malformed():
    assert parse_region('not-a-region') == 'error'
