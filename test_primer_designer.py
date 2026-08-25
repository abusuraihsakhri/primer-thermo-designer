"""Tests for primer_designer.py — plain assert / unittest, stdlib only."""

import unittest

from primer_designer import (
    calc_tm,
    reverse_complement,
    gc_fraction,
    has_gc_clamp,
    longest_homopolymer_run,
    sliding_window_complementarity,
    find_hairpin,
    design_primer_pairs,
    validate_sequence,
)


class TestSequenceUtils(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATCG"), "CGAT")
        self.assertEqual(reverse_complement("GAATTC"), "GAATTC")  # palindrome

    def test_gc_fraction(self):
        self.assertAlmostEqual(gc_fraction("GGCC"), 1.0)
        self.assertAlmostEqual(gc_fraction("AATT"), 0.0)
        self.assertAlmostEqual(gc_fraction("ATGC"), 0.5)

    def test_validate_sequence_rejects_bad_chars(self):
        with self.assertRaises(ValueError):
            validate_sequence("ATCGX")
        with self.assertRaises(ValueError):
            validate_sequence("")


class TestGcClampAndPolyRun(unittest.TestCase):
    def test_gc_clamp_present(self):
        # last base G/C, last 5 bases contain 2 G/C -> within [1,3]
        self.assertTrue(has_gc_clamp("AAAAAATTATC"))

    def test_gc_clamp_absent_wrong_last_base(self):
        self.assertFalse(has_gc_clamp("AAAAAATTATA"))

    def test_gc_clamp_absent_too_many_gc(self):
        # last 5 bases all G/C -> exceeds max_gc=3
        self.assertFalse(has_gc_clamp("AAAAAGCGCG"))

    def test_homopolymer_run(self):
        self.assertEqual(longest_homopolymer_run("ATCGAAAATGC"), 4)
        self.assertEqual(longest_homopolymer_run("ATCGATCG"), 1)
        self.assertEqual(longest_homopolymer_run("GGGGGG"), 6)


class TestTmCalculation(unittest.TestCase):
    def test_tm_increases_with_gc_content(self):
        at_rich = "AAATATATATATATATATTT"       # 21 nt, 0% GC
        gc_rich = "GGCGCGCGCGCGCGCGCGCGC"       # 21 nt, ~95% GC
        tm_at = calc_tm(at_rich)
        tm_gc = calc_tm(gc_rich)
        self.assertGreater(tm_gc, tm_at)

    def test_tm_reasonable_range_for_typical_primer(self):
        primer = "ACGTGCATCGTAGCTAGCTG"  # 20 nt, mixed composition
        tm = calc_tm(primer)
        # A typical 20-mer at 250nM primer / 50mM Na+ should land in a
        # biologically plausible window.
        self.assertTrue(30.0 < tm < 90.0)

    def test_tm_decreases_with_lower_salt(self):
        primer = "ACGTGCATCGTAGCTAGCTG"
        tm_high_salt = calc_tm(primer, na_conc_m=0.10)
        tm_low_salt = calc_tm(primer, na_conc_m=0.01)
        self.assertGreater(tm_high_salt, tm_low_salt)

    def test_tm_increases_with_primer_concentration(self):
        primer = "ACGTGCATCGTAGCTAGCTG"
        tm_low_conc = calc_tm(primer, primer_conc_m=10e-9)
        tm_high_conc = calc_tm(primer, primer_conc_m=1000e-9)
        self.assertGreater(tm_high_conc, tm_low_conc)

    def test_tm_rejects_too_short_primer(self):
        with self.assertRaises(ValueError):
            calc_tm("A")


class TestComplementarityScoring(unittest.TestCase):
    def test_palindrome_shows_strong_self_complementarity(self):
        # EcoRI site: perfectly self-complementary (its own reverse complement)
        result = sliding_window_complementarity("GAATTC", "GAATTC")
        self.assertEqual(result.max_run, 6)
        self.assertEqual(result.score, 6)

    def test_low_complementarity_sequence(self):
        # A sequence with minimal self-complementarity should score low.
        result = sliding_window_complementarity("AAAAAA", "GGGGGG")
        self.assertEqual(result.max_run, 0)
        self.assertEqual(result.score, 0)

    def test_heterodimer_detects_3prime_complementarity(self):
        # seq_a ends in a run complementary to seq_b's start (antiparallel).
        seq_a = "TTTTTGGGGG"
        seq_b = "CCCCCAAAAA"
        result = sliding_window_complementarity(seq_a, seq_b)
        self.assertGreaterEqual(result.max_run, 5)


class TestHairpinDetection(unittest.TestCase):
    def test_detects_designed_hairpin(self):
        # stem "GGGGG" ... loop "AAAA" ... stem-complement "CCCCC"
        seq = "GGGGG" + "AAAA" + "CCCCC"
        hp = find_hairpin(seq, min_stem=3, min_loop=3)
        self.assertEqual(hp.stem_len, 5)
        self.assertEqual(hp.loop_len, 4)

    def test_no_hairpin_in_non_self_complementary_sequence(self):
        seq = "ATATATATATATATATATAT"
        hp = find_hairpin(seq, min_stem=4, min_loop=3)
        self.assertEqual(hp.stem_len, 0)


class TestPrimerPairDesign(unittest.TestCase):
    def setUp(self):
        # Build a synthetic target: a forward-primer-friendly region,
        # a filler region, and a region whose reverse complement is a
        # reverse-primer-friendly sequence, so a valid pair should exist.
        forward_region = "ACGTGCATCGTAGCTAGCTG"      # 20 nt, GC clamp at 3' (G)
        filler = "ATCG" * 28                          # 112 nt, low-complexity spacer
        reverse_region_top_strand = "CAGCTAGCTACGATGCACGT"  # 20 nt
        self.target = forward_region + filler + reverse_region_top_strand
        # Sanity: length should exceed the minimum product size used below.
        self.assertGreater(len(self.target), 130)

    def test_finds_at_least_one_valid_pair(self):
        pairs = design_primer_pairs(
            self.target,
            product_size_range=(50, 145),
            primer_len_range=(18, 22),
            tm_range=(30.0, 90.0),   # loose bounds so the synthetic target passes
            gc_range=(0.20, 0.80),
            top_n=5,
        )
        self.assertGreater(len(pairs), 0)

        for pair in pairs:
            # Structural correctness: primers are real substrings/revcomps
            # of the target, non-overlapping, and product size is in range.
            fwd_window = self.target[pair.forward.start:pair.forward.end]
            self.assertEqual(fwd_window, pair.forward.seq)

            rev_window = self.target[pair.reverse.start:pair.reverse.end]
            self.assertEqual(reverse_complement(rev_window), pair.reverse.seq)

            self.assertLess(pair.forward.end, pair.reverse.start + 1)
            self.assertTrue(50 <= pair.product_size <= 145)

    def test_no_pairs_when_product_size_too_small_to_fit_two_primers(self):
        # A product window smaller than two non-overlapping primers can
        # ever produce should yield no pairs, without raising an error.
        pairs = design_primer_pairs(
            self.target,
            product_size_range=(5, 10),
            primer_len_range=(18, 22),
            tm_range=(30.0, 90.0),
            gc_range=(0.20, 0.80),
            top_n=5,
        )
        self.assertEqual(pairs, [])


if __name__ == "__main__":
    unittest.main()
