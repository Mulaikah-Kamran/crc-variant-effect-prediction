"""
Shared functions used by the pipeline scripts and directly exercised by
tests/test_pipeline.py — kept here once rather than duplicated so a test
against this module is a test against what actually ran.
"""
import re
import statistics as stats

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
PROTEIN_CHANGE_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

ESM2_CONTEXT_LIMIT = 1022  # confirmed empirically from model config, all ESM2 sizes


def parse_protein_change(pc):
    """Parse a simple missense notation like 'R175H' into (wt, pos, mt).
    Returns None for anything else (nonsense, frameshift, indels, etc.)."""
    if not pc:
        return None
    m = PROTEIN_CHANGE_RE.match(pc.strip())
    if not m:
        return None
    wt, pos, mt = m.group(1), int(m.group(2)), m.group(3)
    return wt, pos, mt


def validate_variant(sequence, position, wt_residue, mutant_residue):
    """
    Validate a candidate variant against a protein sequence.
    Returns (qc_status, exclusion_reason_or_None).
    """
    if wt_residue not in VALID_AA or mutant_residue not in VALID_AA:
        return "excluded_invalid_residue", f"non-standard amino acid code: wt={wt_residue}, mt={mutant_residue}"

    seq_len = len(sequence)
    if not (1 <= position <= seq_len):
        return "excluded_out_of_bounds", f"position {position} outside sequence length {seq_len}"

    actual = sequence[position - 1]
    if actual != wt_residue:
        return "excluded_wt_mismatch", (
            f"record states WT={wt_residue} at position {position}, "
            f"but sequence has '{actual}' there"
        )

    return "pass", None


def get_window(seq, pos, window_size=1023):
    """
    Return (windowed_sequence, new_position_1indexed) for scoring a protein
    that exceeds ESM2's context limit. Centers a window of `window_size`
    residues on `pos`, shifting inward near either terminus so the window
    never runs past the sequence bounds. If the full sequence already fits,
    returns it unchanged.
    """
    L = len(seq)
    if L <= window_size:
        return seq, pos

    half = window_size // 2
    start = pos - 1 - half
    end = start + window_size

    if start < 0:
        end -= start
        start = 0
    if end > L:
        start -= (end - L)
        end = L
    start = max(start, 0)

    return seq[start:end], pos - start


def tukey_hinges_iqr(values):
    """Q1/Q3 via Tukey's hinges (median-of-halves). See docs/methodology.md
    for why this convention is named explicitly rather than left implicit."""
    if len(values) < 2:
        return None
    sorted_v = sorted(values)
    q1 = stats.median(sorted_v[: len(sorted_v) // 2])
    upper_half = sorted_v[(len(sorted_v) + 1) // 2:]
    q3 = stats.median(upper_half) if upper_half else q1
    return q1, q3
