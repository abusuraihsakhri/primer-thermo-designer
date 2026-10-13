"""Regression cases for scientific boundaries and HMAC verification."""
import os
import pytest

os.environ.setdefault("AUDIT_SECRET_KEY", "test-audit-key-for-unit-tests-only")

from primer_designer import calc_tm, design_primer_pairs, parse_fasta, has_gc_clamp
from agents.base import AuditTrail


@pytest.mark.parametrize("primer_conc, salt_conc", [
    (0.0, 0.05), (-1.0, 0.05), (float("nan"), 0.05),
    (2.5e-7, 0.0), (2.5e-7, -0.05), (2.5e-7, float("inf")),
])
def test_tm_rejects_nonpositive_or_nonfinite_concentrations(primer_conc, salt_conc):
    with pytest.raises(ValueError, match="concentration"):
        calc_tm("ACGTGCATCGTAGCTAGCTG", primer_conc, salt_conc)


def test_fasta_rejects_multiple_templates():
    with pytest.raises(ValueError, match="Multiple FASTA"):
        parse_fasta(">first\nACGT\n>second\nTGCA")


def test_empty_gc_clamp():
    assert has_gc_clamp("") is False


@pytest.mark.parametrize("changes", [
    {"current_hash": "not-a-real-signature"},
    {"event_type": "FORGED_EVENT"},
    {"payload_hash": "tampered"},
])
def test_audit_trail_recomputes_hmac(changes):
    trail = AuditTrail("isolated-test-key")
    trail.log("actor", "tier", "event", {"x": 1})
    assert trail.verify_integrity()
    trail.logs[0].update(changes)
    assert trail.verify_integrity() is False


def test_adjacent_primers_can_span_exact_product():
    target = "ACGTGCATCGTAGCTAGCTG" + "CAGCTAGCTACGATGCACGT"
    pairs = design_primer_pairs(
        target, product_size_range=(40, 40), primer_len_range=(20, 20),
        tm_range=(-100, 120), gc_range=(0, 1), require_gc_clamp=False,
        three_prime_dg_threshold=-100, top_n=20,
    )
    assert any(p.forward.start == 0 and p.reverse.start == 20 and p.product_size == 40 for p in pairs)


def test_forward_candidate_at_last_valid_start_is_not_skipped():
    target = ("ATCG" * 5) + ("ACGTGCATCGTAGCTAGCTG" + "CAGCTAGCTACGATGCACGT")
    pairs = design_primer_pairs(
        target, product_size_range=(40, 40), primer_len_range=(20, 20),
        tm_range=(-100, 120), gc_range=(0, 1), require_gc_clamp=False,
        three_prime_dg_threshold=-100, top_n=1000,
    )
    assert any(p.forward.start == 20 and p.reverse.start == 40 for p in pairs)


@pytest.mark.parametrize("kwargs", [
    {"primer_len_range": (1, 20)},
    {"primer_len_range": (25, 18)},
    {"gc_range": (-0.1, 0.7)},
    {"tm_range": (65, 55)},
    {"primer_conc_nm": 0},
    {"na_conc_mm": float("nan")},
    {"top_n": 0},
])
def test_invalid_design_controls_raise(kwargs):
    with pytest.raises(ValueError):
        design_primer_pairs("ACGT" * 40, product_size_range=(40, 100), **kwargs)
