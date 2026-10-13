"""
Focused tests for security fixes and edge cases in Primer Thermo Designer.
"""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Set test audit key before importing agents
os.environ.setdefault("AUDIT_SECRET_KEY", "test-audit-key-for-unit-tests-only")

import pytest
from primer_designer import (
    calc_tm,
    validate_sequence,
    gc_fraction,
    reverse_complement,
    design_primer_pairs,
    find_hairpin,
    sliding_window_complementarity,
    calc_dg37,
    parse_fasta,
)
from agents.base import AuditTrail, PHIGuard, SecurityException, assert_no_phi


class TestAuditTrailSecurity:
    """Tests for HMAC-SHA256 audit trail security."""

    def test_audit_trail_requires_secret_key(self):
        """AuditTrail should raise SecurityException when no key is provided."""
        original_key = os.environ.pop("AUDIT_SECRET_KEY", None)
        try:
            with pytest.raises(SecurityException):
                AuditTrail()
        finally:
            if original_key:
                os.environ["AUDIT_SECRET_KEY"] = original_key

    def test_audit_trail_accepts_env_key(self):
        """AuditTrail should accept key from environment variable."""
        os.environ["AUDIT_SECRET_KEY"] = "test-key-123"
        trail = AuditTrail()
        assert trail.secret_key == b"test-key-123"

    def test_audit_trail_chained_integrity(self):
        """Audit trail should maintain cryptographic chain integrity."""
        os.environ["AUDIT_SECRET_KEY"] = "test-key-123"
        trail = AuditTrail()
        trail.log("test", "unit", "TEST_EVENT", {"data": "value1"})
        trail.log("test", "unit", "TEST_EVENT", {"data": "value2"})
        assert trail.verify_integrity() is True
        assert len(trail.get_trail()) == 2

    def test_audit_trail_detects_tampering(self):
        """Audit trail should detect tampered chain links."""
        os.environ["AUDIT_SECRET_KEY"] = "test-key-123"
        trail = AuditTrail()
        trail.log("test", "unit", "TEST_EVENT", {"data": "value1"})
        trail.log("test", "unit", "TEST_EVENT", {"data": "value2"})
        # Tamper with the chain: modify the first entry's hash
        trail.logs[0]["current_hash"] = "tampered_hash"
        # The second entry's prev_hash should no longer match
        assert trail.verify_integrity() is False

    def test_audit_trail_blocks_phi_in_logs(self):
        """Audit trail should reject log entries containing PHI."""
        os.environ["AUDIT_SECRET_KEY"] = "test-key-123"
        trail = AuditTrail()
        with pytest.raises(SecurityException):
            trail.log("test", "unit", "PHI_EVENT", {"patient": "MRN-123456"})


class TestPHIGuard:
    """Tests for PHI outbound guard."""

    def test_blocks_mrn(self):
        with pytest.raises(SecurityException):
            assert_no_phi("Patient MRN-123456")

    def test_blocks_ssn(self):
        with pytest.raises(SecurityException):
            assert_no_phi("SSN: 123-45-6789")

    def test_blocks_phone_number(self):
        with pytest.raises(SecurityException):
            assert_no_phi("Call 555-123-4567")

    def test_blocks_email(self):
        with pytest.raises(SecurityException):
            assert_no_phi("Email: user@example.com")

    def test_blocks_patient_name(self):
        with pytest.raises(SecurityException):
            assert_no_phi("Patient Name: John Smith")

    def test_allows_clean_text(self):
        assert_no_phi("Analytical specimen KEY-001 optimal")  # Should not raise

    def test_handles_empty_string(self):
        assert_no_phi("")  # Should not raise

    def test_handles_none(self):
        assert_no_phi(None)  # Should not raise


class TestSequenceValidation:
    """Tests for DNA sequence validation edge cases."""

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError, match="empty"):
            validate_sequence("")

    def test_whitespace_only_raises(self):
        with pytest.raises(ValueError, match="empty"):
            validate_sequence("   ")

    def test_lowercase_converted(self):
        result = validate_sequence("atcg")
        assert result == "ATCG"

    def test_mixed_case_converted(self):
        result = validate_sequence("AtCg")
        assert result == "ATCG"

    def test_whitespace_stripped(self):
        result = validate_sequence("  ATCG  ")
        assert result == "ATCG"

    def test_invalid_chars_rejected(self):
        with pytest.raises(ValueError, match="non-ACGT"):
            validate_sequence("ATCGX")

    def test_rna_uracil_rejected(self):
        with pytest.raises(ValueError, match="non-ACGT"):
            validate_sequence("AUCG")

    def test_numbers_rejected(self):
        with pytest.raises(ValueError, match="non-ACGT"):
            validate_sequence("ATCG123")


class TestFastaParser:
    """Tests for FASTA file parsing."""

    def test_simple_fasta(self):
        fasta = ">seq1\nATCG\nATCG"
        result = parse_fasta(fasta)
        assert result == "ATCGATCG"

    def test_multiple_sequences_rejected(self):
        fasta = ">seq1\nATCG\n>seq2\nGCTA"
        with pytest.raises(ValueError, match="Multiple FASTA"):
            parse_fasta(fasta)

    def test_raw_sequence_no_header(self):
        result = parse_fasta("ATCGATCG")
        assert result == "ATCGATCG"

    def test_empty_lines_ignored(self):
        fasta = ">seq1\n\nATCG\n\n"
        result = parse_fasta(fasta)
        assert result == "ATCG"


class TestCalcTmEdgeCases:
    """Tests for Tm calculation edge cases."""

    def test_gc_fraction_extremes(self):
        assert gc_fraction("AAAA") == 0.0
        assert gc_fraction("GCGC") == 1.0
        assert gc_fraction("") == 0.0

    def test_reverse_complement(self):
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"

    def test_self_complementary_palindrome(self):
        # Palindromic sequence
        result = reverse_complement("GAATTC")
        assert result == "GAATTC"

    def test_dg37_calculation(self):
        dg = calc_dg37("ATCG")
        assert isinstance(dg, float)

    def test_dg37_rejects_short(self):
        with pytest.raises(ValueError):
            calc_dg37("A")


class TestPrimerDesignEdgeCases:
    """Tests for primer pair design edge cases."""

    def test_no_valid_pairs_raises_for_short_target(self):
        """Short target that can't fit the product size raises ValueError."""
        with pytest.raises(ValueError):
            design_primer_pairs(
                "ACGT",
                product_size_range=(50, 100),
                primer_len_range=(18, 25),
            )

    def test_invalid_product_size_raises(self):
        with pytest.raises(ValueError):
            design_primer_pairs(
                "ACGTACGTACGT",
                product_size_range=(100, 50),  # min > max
            )

    def test_target_too_short_raises(self):
        with pytest.raises(ValueError):
            design_primer_pairs(
                "ACGT",
                product_size_range=(100, 200),
                primer_len_range=(18, 25),
            )

    def test_hairpin_no_stem(self):
        result = find_hairpin("ATATATATAT", min_stem=4, min_loop=3)
        assert result.stem_len == 0

    def test_complementarity_no_match(self):
        # AAAAAA and GGGGGG are NOT complementary (A pairs with T, G pairs with C)
        result = sliding_window_complementarity("AAAAAA", "GGGGGG")
        assert result.score == 0
        assert result.max_run == 0
