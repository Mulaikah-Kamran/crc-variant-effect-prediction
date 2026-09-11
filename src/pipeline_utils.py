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


# ---- BLOSUM62 classical substitution-matrix baseline ----
#
# Standard BLOSUM62 (Henikoff & Henikoff, 1992), the version distributed with
# NCBI BLAST. Defined here as the upper triangle only (20 diagonal + 190
# off-diagonal pairs) and mirrored programmatically below, so symmetry is
# guaranteed by construction rather than relying on 400 hand-typed cells
# being individually correct. tests/test_pipeline.py checks the result
# against several independently-known reference values and against full
# matrix symmetry.
_BLOSUM62_ORDER = "ARNDCQEGHILKMFPSTWYV"

_BLOSUM62_UPPER = {
    "A": {"A": 4, "R": -1, "N": -2, "D": -2, "C": 0, "Q": -1, "E": -1, "G": 0, "H": -2, "I": -1,
          "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 0, "W": -3, "Y": -2, "V": 0},
    "R": {"R": 5, "N": 0, "D": -2, "C": -3, "Q": 1, "E": 0, "G": -2, "H": 0, "I": -3, "L": -2,
          "K": 2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3},
    "N": {"N": 6, "D": 1, "C": -3, "Q": 0, "E": 0, "G": 0, "H": 1, "I": -3, "L": -3, "K": 0,
          "M": -2, "F": -3, "P": -2, "S": 1, "T": 0, "W": -4, "Y": -2, "V": -3},
    "D": {"D": 6, "C": -3, "Q": 0, "E": 2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3,
          "F": -3, "P": -1, "S": 0, "T": -1, "W": -4, "Y": -3, "V": -3},
    "C": {"C": 9, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2,
          "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1},
    "Q": {"Q": 5, "E": 2, "G": -2, "H": 0, "I": -3, "L": -2, "K": 1, "M": 0, "F": -3, "P": -1,
          "S": 0, "T": -1, "W": -2, "Y": -1, "V": -2},
    "E": {"E": 5, "G": -2, "H": 0, "I": -3, "L": -3, "K": 1, "M": -2, "F": -3, "P": -1, "S": 0,
          "T": -1, "W": -3, "Y": -2, "V": -2},
    "G": {"G": 6, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S": 0, "T": -2,
          "W": -2, "Y": -3, "V": -3},
    "H": {"H": 8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2,
          "Y": 2, "V": -3},
    "I": {"I": 4, "L": 2, "K": -3, "M": 1, "F": 0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V": 3},
    "L": {"L": 4, "K": -2, "M": 2, "F": 0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V": 1},
    "K": {"K": 5, "M": -1, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "M": {"M": 5, "F": 0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V": 1},
    "F": {"F": 6, "P": -4, "S": -2, "T": -2, "W": 1, "Y": 3, "V": -1},
    "P": {"P": 7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"S": 4, "T": 1, "W": -3, "Y": -2, "V": -2},
    "T": {"T": 5, "W": -2, "Y": -2, "V": 0},
    "W": {"W": 11, "Y": 2, "V": -3},
    "Y": {"Y": 7, "V": -1},
    "V": {"V": 4},
}


def _build_blosum62_full():
    """Mirror the upper-triangle table into a full, symmetric 20x20 lookup."""
    full = {aa: {} for aa in _BLOSUM62_ORDER}
    for a in _BLOSUM62_ORDER:
        for b in _BLOSUM62_ORDER:
            if b in _BLOSUM62_UPPER.get(a, {}):
                score = _BLOSUM62_UPPER[a][b]
            elif a in _BLOSUM62_UPPER.get(b, {}):
                score = _BLOSUM62_UPPER[b][a]
            else:
                score = None
            full[a][b] = score
    return full


BLOSUM62 = _build_blosum62_full()


def blosum62_score(wt_residue, mutant_residue):
    """
    Look up the standard BLOSUM62 substitution score for a WT->mutant amino
    acid pair. Higher (less negative / more positive) values indicate a more
    evolutionarily tolerated substitution; this is a classical log-odds
    substitution score, not a probability, pathogenicity score, or
    cancer-specific measure. Raises ValueError for non-standard residues
    rather than silently returning a placeholder.
    """
    if wt_residue not in VALID_AA or mutant_residue not in VALID_AA:
        raise ValueError(f"Invalid amino acid code(s): wt={wt_residue!r}, mt={mutant_residue!r}")
    return BLOSUM62[wt_residue][mutant_residue]
