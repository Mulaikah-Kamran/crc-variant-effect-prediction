"""
Unit tests for src/pipeline_utils.py — the functions actually used by the
pipeline scripts (03, 04, 05 import from this module directly, not copies
of it), so these tests exercise production logic.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from pipeline_utils import (
    parse_protein_change,
    validate_variant,
    get_window,
    tukey_hinges_iqr,
    ESM2_CONTEXT_LIMIT,
)


# ---- parse_protein_change ----

def test_parse_simple_missense():
    assert parse_protein_change("R175H") == ("R", 175, "H")


def test_parse_single_digit_position():
    assert parse_protein_change("K9R") == ("K", 9, "R")


def test_parse_rejects_nonsense():
    assert parse_protein_change("R1443*") is None


def test_parse_rejects_frameshift_notation():
    assert parse_protein_change("A2fs*12") is None


def test_parse_rejects_empty_string():
    assert parse_protein_change("") is None
    assert parse_protein_change(None) is None


def test_parse_rejects_indel_notation():
    assert parse_protein_change("A2del") is None


# ---- validate_variant ----

def test_validate_passes_correct_variant():
    seq = "MEDYTKIEKIGEGTYGVVYKG"
    # position 1 is 'M'
    status, reason = validate_variant(seq, 1, "M", "V")
    assert status == "pass"
    assert reason is None


def test_validate_catches_wt_mismatch():
    seq = "MEDYTKIEKIGEGTYGVVYKG"
    status, reason = validate_variant(seq, 1, "A", "V")  # actual is M, not A
    assert status == "excluded_wt_mismatch"
    assert "M" in reason


def test_validate_catches_out_of_bounds():
    seq = "MEDY"
    status, reason = validate_variant(seq, 100, "M", "V")
    assert status == "excluded_out_of_bounds"


def test_validate_catches_invalid_residue():
    seq = "MEDY"
    status, reason = validate_variant(seq, 1, "M", "X")  # X not a standard AA here
    assert status == "excluded_invalid_residue"


def test_validate_position_is_one_indexed():
    seq = "ABCDE"
    # position 5 (1-indexed) should be 'E'
    status, _ = validate_variant(seq, 5, "E", "A")
    assert status == "pass"
    status, _ = validate_variant(seq, 5, "D", "A")  # off-by-one check
    assert status == "excluded_wt_mismatch"


# ---- get_window ----

def test_get_window_no_window_needed_when_short():
    seq = "A" * 500
    windowed, pos = get_window(seq, 250, window_size=1023)
    assert windowed == seq
    assert pos == 250


def test_get_window_applies_when_sequence_exceeds_limit():
    seq = "A" * 1068  # PIK3CA's actual length
    windowed, pos = get_window(seq, 545, window_size=1023)  # a middle position
    assert len(windowed) == 1023
    assert windowed[pos - 1] == "A"


def test_get_window_shifts_inward_near_n_terminus():
    seq = "".join(chr(65 + (i % 26)) for i in range(1068))  # distinct-ish sequence
    pos = 13  # near the very start, like PIK3CA's I69N... use a small position
    windowed, new_pos = get_window(seq, pos, window_size=1023)
    assert len(windowed) == 1023
    assert windowed[new_pos - 1] == seq[pos - 1]
    # window must start exactly at sequence position 0 when shifted this far left
    assert seq.startswith(windowed[:20])


def test_get_window_shifts_inward_near_c_terminus():
    seq = "".join(chr(65 + (i % 26)) for i in range(1068))
    pos = 1060  # near the very end
    windowed, new_pos = get_window(seq, pos, window_size=1023)
    assert len(windowed) == 1023
    assert windowed[new_pos - 1] == seq[pos - 1]
    assert seq.endswith(windowed[-20:])


def test_get_window_never_exceeds_window_size_across_full_range():
    seq = "".join(chr(65 + (i % 26)) for i in range(1068))
    for pos in range(1, 1069, 37):  # sample across the whole protein
        windowed, new_pos = get_window(seq, pos, window_size=1023)
        assert len(windowed) <= 1023
        assert windowed[new_pos - 1] == seq[pos - 1]


def test_context_limit_matches_confirmed_esm2_value():
    # This is the value confirmed empirically against the model config
    # (see docs/methodology.md, Section 5) - not a made-up constant.
    assert ESM2_CONTEXT_LIMIT == 1022


# ---- tukey_hinges_iqr ----

def test_iqr_known_even_sample():
    # [1,2,3,4]: lower half [1,2] -> median 1.5; upper half [3,4] -> median 3.5
    q1, q3 = tukey_hinges_iqr([1, 2, 3, 4])
    assert q1 == 1.5
    assert q3 == 3.5


def test_iqr_known_odd_sample():
    # [1,2,3,4,5]: median is 5 excluded from both halves under this convention
    # lower half [1,2] -> 1.5; upper half [4,5] -> 4.5
    q1, q3 = tukey_hinges_iqr([1, 2, 3, 4, 5])
    assert q1 == 1.5
    assert q3 == 4.5


def test_iqr_too_small_returns_none():
    assert tukey_hinges_iqr([5]) is None
    assert tukey_hinges_iqr([]) is None


def test_iqr_order_independent():
    a = tukey_hinges_iqr([5, 1, 3, 2, 4])
    b = tukey_hinges_iqr([1, 2, 3, 4, 5])
    assert a == b
