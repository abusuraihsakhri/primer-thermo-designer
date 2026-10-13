"""
Primer Thermo Designer
=======================

Designs and validates PCR primer pairs from a target DNA sequence using
the SantaLucia (1998) unified nearest-neighbor (NN) thermodynamic
parameters for melting-temperature (Tm) prediction, combined with
sliding-window complementarity scoring for self-dimer, heterodimer and
hairpin screening (in the spirit of Primer3's core heuristics).

Reference:
    SantaLucia J Jr. "A unified view of polymer, dumbbell, and
    oligonucleotide DNA nearest-neighbor thermodynamics."
    Proc Natl Acad Sci USA. 1998;95(4):1460-1465.
"""

from __future__ import annotations

import argparse
import math
import sys
import heapq
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

# ---------------------------------------------------------------------------
# Nearest-neighbor thermodynamic parameters
# ---------------------------------------------------------------------------
# SantaLucia (1998) unified NN parameters.
# Keys are 5'->3' dinucleotide steps. dH in kcal/mol, dS in cal/(mol*K).
# The 16 possible dinucleotide steps collapse to 10 unique duplexes because
# a step and its reverse complement represent the same physical duplex.
NN_PARAMS = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

# Initiation parameters (SantaLucia 1998), applied once per primer end.
INIT_GC = (0.1, -2.8)   # dH, dS for a terminal G or C
INIT_AT = (2.3, 4.1)    # dH, dS for a terminal A or T

GAS_CONSTANT_R = 1.987  # cal / (mol * K)

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

VALID_BASES = set("ACGT")


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------
def validate_sequence(seq: str) -> str:
    """Uppercase and validate a DNA sequence contains only A/C/G/T."""
    if not isinstance(seq, str):
        raise ValueError("Sequence must be a string.")
    seq = seq.strip().upper()
    if not seq:
        raise ValueError("Sequence is empty.")
    bad = set(seq) - VALID_BASES
    if bad:
        raise ValueError(
            f"Sequence contains non-ACGT characters: {''.join(sorted(bad))}"
        )
    return seq


def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq)


def parse_fasta(text: str) -> str:
    """Read one FASTA record or a raw DNA sequence.

    Multiple records are rejected: joining them would create an artificial
    PCR template and may produce primers spanning record boundaries.
    """
    if not isinstance(text, str):
        raise ValueError("FASTA input must be text.")
    sequence = []
    record_count = 0
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            record_count += 1
            if record_count > 1:
                raise ValueError("Multiple FASTA records are not supported; select one target.")
            if sequence:
                raise ValueError("FASTA header must precede sequence lines.")
        else:
            sequence.append(line)
    return "".join(sequence)


# ---------------------------------------------------------------------------
# Tm calculation (SantaLucia 1998 unified NN parameters)
# ---------------------------------------------------------------------------
def calc_tm(
    primer: str,
    primer_conc_m: float = 250e-9,
    na_conc_m: float = 50e-3,
) -> float:
    """Predict primer melting temperature (deg C) using the SantaLucia
    (1998) unified nearest-neighbor thermodynamic parameters, with a
    salt correction applied to the total entropy and the two-state
    Tm equation.

    primer_conc_m: total strand (primer) concentration in mol/L.
    na_conc_m: monovalent cation ([Na+] equivalent) concentration in mol/L.
    """
    seq = validate_sequence(primer)
    n = len(seq)
    if n < 2:
        raise ValueError("Primer must be at least 2 nt long for NN Tm calculation.")
    for value, label in ((primer_conc_m, "Primer concentration"), (na_conc_m, "Salt concentration")):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be finite and greater than zero.")

    dH = 0.0
    dS = 0.0

    # Initiation contributions at both ends of the duplex.
    for end_base in (seq[0], seq[-1]):
        h, s = INIT_GC if end_base in "GC" else INIT_AT
        dH += h
        dS += s

    # Nearest-neighbor stacking contributions.
    for i in range(n - 1):
        step = seq[i:i + 2]
        h, s = NN_PARAMS[step]
        dH += h
        dS += s

    # Self-complementary (palindromic) correction and strand-count factor x
    # in the Tm equation: x = 1 for self-complementary duplexes (single
    # species at equilibrium), x = 4 for non-self-complementary duplexes
    # of two distinct, equimolar strands.
    if seq == reverse_complement(seq):
        dS += -1.4
        x = 1.0
    else:
        x = 4.0

    # Salt correction (SantaLucia 1998): entropy is adjusted per
    # phosphate (i.e. per NN step) as a function of [Na+].
    dS_adj = dS + 0.368 * (n - 1) * math.log(na_conc_m)

    ct = primer_conc_m
    denom = dS_adj + GAS_CONSTANT_R * math.log(ct / x)
    if denom == 0:
        raise ValueError("Degenerate Tm calculation (zero denominator).")

    tm_kelvin = (dH * 1000.0) / denom
    return tm_kelvin - 273.15


def calc_dg37(seq: str) -> float:
    """Free energy (kcal/mol) at 37C (310.15 K) for a short NN duplex
    window, ignoring salt correction (used for 3'-end stability screening,
    a purely relative/comparative metric)."""
    seq = validate_sequence(seq)
    n = len(seq)
    if n < 2:
        raise ValueError("Need at least 2 nt to compute NN free energy.")
    dH = 0.0
    dS = 0.0
    for end_base in (seq[0], seq[-1]):
        h, s = INIT_GC if end_base in "GC" else INIT_AT
        dH += h
        dS += s
    for i in range(n - 1):
        h, s = NN_PARAMS[seq[i:i + 2]]
        dH += h
        dS += s
    t = 310.15
    return dH - (t * dS / 1000.0)


# ---------------------------------------------------------------------------
# Sequence quality checks
# ---------------------------------------------------------------------------
def has_gc_clamp(seq: str, window: int = 5, min_gc: int = 1, max_gc: int = 3) -> bool:
    """A GC clamp requires the 3'-terminal base to be G or C, and the
    number of G/C bases within the last `window` nt to fall in a
    moderate range (too few = weak 3' anneal, too many = excess
    3' stability / mispriming risk)."""
    if not seq or window < 1:
        return False
    if len(seq) < window:
        window = len(seq)
    tail = seq[-window:]
    gc_count = sum(1 for b in tail if b in "GC")
    return seq[-1] in "GC" and (min_gc <= gc_count <= max_gc)


def longest_homopolymer_run(seq: str) -> int:
    """Length of the longest run of an identical base (poly-X run)."""
    longest = 1
    current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest


def three_prime_stability(seq: str, window: int = 5) -> float:
    """Free energy (kcal/mol) of the last `window` 3'-terminal bases.
    Strongly negative values indicate an overly stable 3' end, which
    raises the risk of mispriming."""
    window = min(window, len(seq))
    return calc_dg37(seq[-window:])


# ---------------------------------------------------------------------------
# Sliding-window complementarity scoring (dimers & hairpins)
# ---------------------------------------------------------------------------
@dataclass
class ComplementarityResult:
    score: int          # total complementary base pairs at best alignment
    max_run: int         # longest consecutive complementary run
    offset: int           # alignment offset at which the best score occurs


def complement_base(b: str) -> str:
    return COMPLEMENT[b]


def sliding_window_complementarity(seq_a: str, seq_b: str) -> ComplementarityResult:
    """Slide seq_b (read 3'->5', i.e. reversed) across seq_a and score
    Watson-Crick complementary base pairs at every possible antiparallel
    alignment offset. This models how two primer strands (or one primer
    folding back on itself) could anneal to each other, and is used for
    both self-dimer / heterodimer screening.

    Returns the alignment with the highest total complementary score
    (ties broken by longest consecutive run).
    """
    seq_a = validate_sequence(seq_a)
    seq_b = validate_sequence(seq_b)
    rev_b = seq_b[::-1]
    la, lb = len(seq_a), len(rev_b)

    best = ComplementarityResult(score=0, max_run=0, offset=0)
    for offset in range(-(lb - 1), la):
        score = 0
        run = 0
        max_run = 0
        for i in range(lb):
            a_idx = offset + i
            if 0 <= a_idx < la:
                if complement_base(seq_a[a_idx]) == rev_b[i]:
                    score += 1
                    run += 1
                    max_run = max(max_run, run)
                else:
                    run = 0
        if score > best.score or (score == best.score and max_run > best.max_run):
            best = ComplementarityResult(score=score, max_run=max_run, offset=offset)
    return best


@dataclass
class HairpinResult:
    stem_len: int
    loop_len: int
    left_start: int
    right_start: int


def find_hairpin(seq: str, min_stem: int = 3, min_loop: int = 3) -> HairpinResult:
    """Local self-complementarity scan for intramolecular hairpins.

    For every candidate outermost base pair (p, q) with p < q, extend
    inward while seq[p+k] is Watson-Crick complementary to seq[q-k],
    stopping when the loop constraint is violated or bases stop
    pairing. Returns the longest stem found that satisfies the minimum
    loop size requirement (stem_len == 0 if none found).
    """
    seq = validate_sequence(seq)
    n = len(seq)
    best = HairpinResult(stem_len=0, loop_len=0, left_start=-1, right_start=-1)

    min_span = min_loop + 2 * min_stem - 1
    for p in range(n):
        for q in range(p + min_span, n):
            k = 0
            max_k = (q - p + 1) // 2  # cannot extend past the midpoint
            while k < max_k and complement_base(seq[p + k]) == seq[q - k]:
                k += 1
            loop_len = (q - k) - (p + k) + 1
            if k >= min_stem and loop_len >= min_loop:
                if k > best.stem_len:
                    best = HairpinResult(
                        stem_len=k, loop_len=loop_len, left_start=p, right_start=q - k + 1
                    )
    return best


# ---------------------------------------------------------------------------
# Primer candidate generation and pairing
# ---------------------------------------------------------------------------
@dataclass
class Primer:
    seq: str
    start: int   # 0-based start in target (forward-strand coordinates)
    end: int     # exclusive end in target (forward-strand coordinates)
    strand: str  # 'forward' or 'reverse'
    tm: float
    gc_percent: float
    self_dimer: ComplementarityResult
    hairpin: HairpinResult
    warnings: List[str] = field(default_factory=list)


@dataclass
class PrimerPair:
    forward: Primer
    reverse: Primer
    product_size: int
    heterodimer: ComplementarityResult
    score: float
    warnings: List[str] = field(default_factory=list)


def evaluate_primer(
    seq: str,
    start: int,
    end: int,
    strand: str,
    *,
    primer_conc_m: float,
    na_conc_m: float,
    gc_min: float,
    gc_max: float,
    tm_min: float,
    tm_max: float,
    max_poly_run: int,
    require_gc_clamp: bool,
    three_prime_dg_threshold: float,
    self_dimer_run_threshold: int,
    hairpin_stem_threshold: int,
) -> Optional[Primer]:
    """Score a candidate primer sequence and apply hard filters. Returns
    None if the primer fails a hard constraint, otherwise a Primer with
    any soft warnings attached."""
    gc = gc_fraction(seq)
    if not (gc_min <= gc <= gc_max):
        return None

    if longest_homopolymer_run(seq) > max_poly_run:
        return None

    if require_gc_clamp and not has_gc_clamp(seq):
        return None

    dg3 = three_prime_stability(seq)
    if dg3 <= three_prime_dg_threshold:
        return None

    tm = calc_tm(seq, primer_conc_m=primer_conc_m, na_conc_m=na_conc_m)
    if not (tm_min <= tm <= tm_max):
        return None

    self_dimer = sliding_window_complementarity(seq, seq)
    hairpin = find_hairpin(seq)

    warnings: List[str] = []
    if self_dimer.max_run >= self_dimer_run_threshold:
        warnings.append(
            f"self-dimer risk: {self_dimer.max_run} consecutive complementary bases"
        )
    if hairpin.stem_len >= hairpin_stem_threshold:
        warnings.append(
            f"hairpin risk: {hairpin.stem_len}bp stem / {hairpin.loop_len}nt loop"
        )

    return Primer(
        seq=seq, start=start, end=end, strand=strand, tm=tm,
        gc_percent=gc * 100.0, self_dimer=self_dimer, hairpin=hairpin,
        warnings=warnings,
    )


def generate_candidates(
    target: str,
    strand: str,
    primer_len_range: Tuple[int, int],
    search_start: int,
    search_end: int,
    **eval_kwargs,
) -> List[Primer]:
    """Enumerate candidate primers of every allowed length at every
    position within [search_start, search_end) of the target."""
    candidates: List[Primer] = []
    min_len, max_len = primer_len_range
    n = len(target)
    search_end = min(search_end, n)
    for length in range(min_len, max_len + 1):
        for start in range(search_start, search_end - length + 1):
            end = start + length
            window = target[start:end]
            if strand == "forward":
                seq = window
            else:
                seq = reverse_complement(window)
            primer = evaluate_primer(seq, start, end, strand, **eval_kwargs)
            if primer is not None:
                candidates.append(primer)
    return candidates


def score_pair(
    forward: Primer,
    reverse: Primer,
    heterodimer: ComplementarityResult,
    target_tm: float,
) -> float:
    """Composite ranking score (lower is better) combining Tm matching
    between the two primers, deviation from the target Tm, and
    specificity penalties from dimer/hairpin complementarity scans."""
    tm_diff = abs(forward.tm - reverse.tm)
    tm_dev = abs(((forward.tm + reverse.tm) / 2.0) - target_tm)
    dimer_penalty = heterodimer.max_run
    hairpin_penalty = max(forward.hairpin.stem_len, reverse.hairpin.stem_len)
    self_dimer_penalty = max(forward.self_dimer.max_run, reverse.self_dimer.max_run)
    return (
        tm_diff * 2.0
        + tm_dev * 1.0
        + dimer_penalty * 1.5
        + hairpin_penalty * 1.0
        + self_dimer_penalty * 1.0
    )


def design_primer_pairs(
    target: str,
    product_size_range: Tuple[int, int],
    primer_len_range: Tuple[int, int] = (18, 25),
    tm_range: Tuple[float, float] = (57.0, 63.0),
    target_tm: float = 60.0,
    primer_conc_nm: float = 250.0,
    na_conc_mm: float = 50.0,
    gc_range: Tuple[float, float] = (0.30, 0.70),
    max_poly_run: int = 4,
    require_gc_clamp: bool = True,
    three_prime_dg_threshold: float = -9.0,
    self_dimer_run_threshold: int = 5,
    hairpin_stem_threshold: int = 4,
    heterodimer_run_threshold: int = 5,
    top_n: int = 5,
) -> List[PrimerPair]:
    """Enumerate candidate primer pairs across the target sequence within
    the requested amplicon-size bounds, screen each pair for dimer /
    hairpin risk, and return the top-N pairs ranked by combined Tm-match
    and specificity score."""
    target = validate_sequence(target)
    n = len(target)
    min_product, max_product = product_size_range
    min_len, max_len = primer_len_range
    if min_len < 2 or max_len < min_len:
        raise ValueError("Invalid primer_len_range; lengths must be at least 2 nt.")
    if not (math.isfinite(tm_range[0]) and math.isfinite(tm_range[1]) and
            tm_range[0] <= tm_range[1] and math.isfinite(target_tm)):
        raise ValueError("Invalid Tm range or target Tm.")
    if not (0 <= gc_range[0] <= gc_range[1] <= 1):
        raise ValueError("Invalid GC range; expected fractions between zero and one.")
    if max_poly_run < 1 or top_n < 1:
        raise ValueError("max_poly_run and top_n must be positive integers.")
    if not (math.isfinite(primer_conc_nm) and primer_conc_nm > 0 and
            math.isfinite(na_conc_mm) and na_conc_mm > 0):
        raise ValueError("Primer and salt concentrations must be finite and positive.")
    if not math.isfinite(three_prime_dg_threshold):
        raise ValueError("3' stability threshold must be finite.")
    if min_product < 2 or max_product < min_product:
        raise ValueError("Invalid product_size_range.")
    if n < min_product:
        raise ValueError(
            f"Target sequence ({n} nt) is shorter than the minimum product size "
            f"({min_product} nt)."
        )

    primer_conc_m = primer_conc_nm * 1e-9
    na_conc_m = na_conc_mm * 1e-3

    eval_kwargs = dict(
        primer_conc_m=primer_conc_m,
        na_conc_m=na_conc_m,
        gc_min=gc_range[0],
        gc_max=gc_range[1],
        tm_min=tm_range[0],
        tm_max=tm_range[1],
        max_poly_run=max_poly_run,
        require_gc_clamp=require_gc_clamp,
        three_prime_dg_threshold=three_prime_dg_threshold,
        self_dimer_run_threshold=self_dimer_run_threshold,
        hairpin_stem_threshold=hairpin_stem_threshold,
    )

    # Forward primers may start anywhere that still leaves room for the
    # minimum product size; reverse primers may end anywhere after the
    # minimum product size from the start of the target.
    fwd_candidates = generate_candidates(
        target, "forward", primer_len_range,
        search_start=0, search_end=min(n, n - min_product + max_len + 1),
        **eval_kwargs,
    )
    rev_candidates = generate_candidates(
        target, "reverse", primer_len_range,
        search_start=0, search_end=n,
        **eval_kwargs,
    )

    # Restrict pairing to reverse primers whose end coordinates can
    # yield an amplicon within the requested product range. Keep only
    # the best top_n matches to avoid quadratic memory growth.
    rev_candidates.sort(key=lambda candidate: candidate.end)
    rev_ends = [candidate.end for candidate in rev_candidates]
    # Protect both the browser and the API against prohibitively expensive
    # unrestricted searches. The estimate is conservative: overlapping
    # primer pairs are counted even though they are filtered later.
    candidate_evaluations = 0
    bounds = []
    for f in fwd_candidates:
        lower = bisect_left(rev_ends, f.start + min_product)
        upper = bisect_right(rev_ends, f.start + max_product)
        candidate_evaluations += upper - lower
        if candidate_evaluations > 500_000:
            raise ValueError(
                "Search exceeds the 500,000-pair screening budget. "
                "Narrow the amplicon, primer-length, GC or Tm ranges."
            )
        bounds.append((lower, upper))

    best_pairs: List[Tuple[float, int, PrimerPair]] = []
    index = 0
    for f, (lower, upper) in zip(fwd_candidates, bounds):
        for r in rev_candidates[lower:upper]:
            product_size = r.end - f.start
            if r.start < f.end:
                continue  # primers must not overlap (adjacent ends are allowed)
            if not (min_product <= product_size <= max_product):
                continue
            heterodimer = sliding_window_complementarity(f.seq, r.seq)
            warnings = []
            if heterodimer.max_run >= heterodimer_run_threshold:
                warnings.append(
                    f"heterodimer risk: {heterodimer.max_run} consecutive "
                    "complementary bases between primers"
                )
            score = score_pair(f, r, heterodimer, target_tm)
            pair = PrimerPair(
                forward=f, reverse=r, product_size=product_size,
                heterodimer=heterodimer, score=score, warnings=warnings,
            )
            entry = (-score, -index, pair)
            index += 1
            if len(best_pairs) < top_n:
                heapq.heappush(best_pairs, entry)
            elif entry > best_pairs[0]:
                heapq.heapreplace(best_pairs, entry)

    return [item[2] for item in sorted(best_pairs, key=lambda e: (-e[0], -e[1]))]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _read_sequence(args: argparse.Namespace) -> str:
    if args.fasta:
        with open(args.fasta, "r", encoding="utf-8") as fh:
            text = fh.read()
        return validate_sequence(parse_fasta(text))
    if args.sequence:
        return validate_sequence(args.sequence)
    raise SystemExit("error: provide --sequence or --fasta")


def cmd_design(args: argparse.Namespace) -> None:
    target = _read_sequence(args)
    pairs = design_primer_pairs(
        target,
        product_size_range=(args.product_min, args.product_max),
        primer_len_range=(args.primer_len_min, args.primer_len_max),
        tm_range=(args.tm_min, args.tm_max),
        target_tm=args.target_tm,
        primer_conc_nm=args.primer_conc,
        na_conc_mm=args.na_conc,
        gc_range=(args.gc_min / 100.0, args.gc_max / 100.0),
        max_poly_run=args.max_poly_run,
        require_gc_clamp=not args.no_gc_clamp,
        top_n=args.top,
    )

    if not pairs:
        print("No primer pairs found satisfying the given constraints.")
        return

    for i, p in enumerate(pairs, start=1):
        print(f"--- Pair {i} (score={p.score:.2f}, product={p.product_size} bp) ---")
        print(
            f"  Forward: {p.forward.seq}  "
            f"[{p.forward.start}-{p.forward.end}]  "
            f"Tm={p.forward.tm:.1f}C  GC={p.forward.gc_percent:.1f}%"
        )
        print(
            f"  Reverse: {p.reverse.seq}  "
            f"[{p.reverse.start}-{p.reverse.end}]  "
            f"Tm={p.reverse.tm:.1f}C  GC={p.reverse.gc_percent:.1f}%"
        )
        all_warnings = p.forward.warnings + p.reverse.warnings + p.warnings
        if all_warnings:
            for w in all_warnings:
                print(f"  Warning: {w}")
        print()


def cmd_check(args: argparse.Namespace) -> None:
    seq = validate_sequence(args.primer)
    tm = calc_tm(seq, primer_conc_m=args.primer_conc * 1e-9, na_conc_m=args.na_conc * 1e-3)
    gc = gc_fraction(seq) * 100.0
    clamp = has_gc_clamp(seq)
    run = longest_homopolymer_run(seq)
    dg3 = three_prime_stability(seq)
    self_dimer = sliding_window_complementarity(seq, seq)
    hairpin = find_hairpin(seq)

    print(f"Primer:            {seq}")
    print(f"Length:            {len(seq)} nt")
    print(f"Tm (NN, salt-adj): {tm:.2f} C")
    print(f"GC content:        {gc:.1f}%")
    print(f"GC clamp:          {'yes' if clamp else 'no'}")
    print(f"Longest poly-X run:{run} nt")
    print(f"3' end dG (37C):   {dg3:.2f} kcal/mol")
    print(
        f"Self-dimer:        score={self_dimer.score}, "
        f"max consecutive run={self_dimer.max_run}"
    )
    if hairpin.stem_len > 0:
        print(
            f"Hairpin:           stem={hairpin.stem_len}bp, "
            f"loop={hairpin.loop_len}nt"
        )
    else:
        print("Hairpin:           none detected")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="primer_designer",
        description=(
            "Design and validate PCR primer pairs using SantaLucia (1998) "
            "nearest-neighbor thermodynamics."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    design = sub.add_parser("design", help="Enumerate and rank primer pairs for a target sequence.")
    seq_group = design.add_mutually_exclusive_group(required=True)
    seq_group.add_argument("--sequence", help="Target DNA sequence (ACGT).")
    seq_group.add_argument("--fasta", help="Path to a FASTA file containing the target sequence.")
    design.add_argument("--product-min", type=int, default=100, help="Minimum amplicon size (bp). Default 100.")
    design.add_argument("--product-max", type=int, default=300, help="Maximum amplicon size (bp). Default 300.")
    design.add_argument("--primer-len-min", type=int, default=18, help="Minimum primer length (nt). Default 18.")
    design.add_argument("--primer-len-max", type=int, default=25, help="Maximum primer length (nt). Default 25.")
    design.add_argument("--tm-min", type=float, default=57.0, help="Minimum acceptable primer Tm (C). Default 57.")
    design.add_argument("--tm-max", type=float, default=63.0, help="Maximum acceptable primer Tm (C). Default 63.")
    design.add_argument("--target-tm", type=float, default=60.0, help="Ideal Tm used for ranking (C). Default 60.")
    design.add_argument("--primer-conc", type=float, default=250.0, help="Primer strand concentration (nM). Default 250.")
    design.add_argument("--na-conc", type=float, default=50.0, help="Monovalent salt concentration (mM). Default 50.")
    design.add_argument("--gc-min", type=float, default=30.0, help="Minimum primer GC%% . Default 30.")
    design.add_argument("--gc-max", type=float, default=70.0, help="Maximum primer GC%% . Default 70.")
    design.add_argument("--max-poly-run", type=int, default=4, help="Maximum allowed homopolymer run length. Default 4.")
    design.add_argument("--no-gc-clamp", action="store_true", help="Disable the 3' GC clamp requirement.")
    design.add_argument("--top", type=int, default=5, help="Number of top-ranked pairs to report. Default 5.")
    design.set_defaults(func=cmd_design)

    check = sub.add_parser("check", help="Evaluate a single primer sequence.")
    check.add_argument("primer", help="Primer sequence (ACGT).")
    check.add_argument("--primer-conc", type=float, default=250.0, help="Primer strand concentration (nM). Default 250.")
    check.add_argument("--na-conc", type=float, default=50.0, help="Monovalent salt concentration (mM). Default 50.")
    check.set_defaults(func=cmd_check)

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
