#!/usr/bin/env python3
"""
Bacterial Genome Synteny Agent & Comparative Microbial Genomics Toolkit.

Implements:
- Syntenic Collinear Block Detection & Inversion Identification
- Synteny Conservation Metrics: Conserved Adjacency Fraction, Breakpoint Distance
- Rank Correlation Metrics: Spearman's Rho & Kendall's Tau on Ortholog Order
- Average Nucleotide Identity (ANI) & Microbial Species Delimitation Cutoffs
- ASCII Dotplot Matrix Rendering (Collinear '\', Inversion '/', Unaligned '.')
- Inversion Breakpoint Junction Detection & PCR Primer Design (Tm, GC clamp)
- Pangenome Partitioning (Core >=99%, Soft-Core >=95%, Shell 15-95%, Cloud <15%)

Pure Python Standard Library (no external dependencies required).
"""

from __future__ import annotations
import math
import re
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Set, Tuple, Optional, Any


@dataclass
class SyntenyBlock:
    """Collinear or inverted syntenic block between reference and query genomes."""
    block_id: str
    ref_start: int
    ref_end: int
    qry_start: int
    qry_end: int
    strand: str = "+"  # '+' (collinear forward) or '-' (inverted reverse)
    identity_pct: float = 98.0
    ref_genes: List[str] = field(default_factory=list)
    qry_genes: List[str] = field(default_factory=list)

    @property
    def ref_length(self) -> int:
        return abs(self.ref_end - self.ref_start)

    @property
    def qry_length(self) -> int:
        return abs(self.qry_end - self.qry_start)

    @property
    def is_inverted(self) -> bool:
        if self.strand == "-":
            return True
        return (self.qry_end - self.qry_start) * (self.ref_end - self.ref_start) < 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "block_id": self.block_id,
            "ref_start": self.ref_start,
            "ref_end": self.ref_end,
            "qry_start": self.qry_start,
            "qry_end": self.qry_end,
            "strand": self.strand,
            "identity_pct": round(self.identity_pct, 2),
            "ref_length": self.ref_length,
            "qry_length": self.qry_length,
            "is_inverted": self.is_inverted,
            "ref_genes": self.ref_genes,
            "qry_genes": self.qry_genes,
        }


@dataclass
class SyntenyAnalysisResult:
    """Overall comparative synteny and ANI report."""
    reference_name: str
    query_name: str
    ref_genome_length: int
    qry_genome_length: int
    num_blocks: int
    collinear_blocks: int
    inverted_blocks: int
    ani_pct: float
    ref_coverage_pct: float
    qry_coverage_pct: float
    taxonomic_call: str
    conserved_adjacency_score: float  # 0.0 to 1.0
    spearman_rho: float              # -1.0 to +1.0
    kendall_tau: float               # -1.0 to +1.0
    breakpoint_count: int
    blocks: List[SyntenyBlock] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reference_name": self.reference_name,
            "query_name": self.query_name,
            "ref_genome_length": self.ref_genome_length,
            "qry_genome_length": self.qry_genome_length,
            "num_blocks": self.num_blocks,
            "collinear_blocks": self.collinear_blocks,
            "inverted_blocks": self.inverted_blocks,
            "ani_pct": round(self.ani_pct, 2),
            "ref_coverage_pct": round(self.ref_coverage_pct, 2),
            "qry_coverage_pct": round(self.qry_coverage_pct, 2),
            "taxonomic_call": self.taxonomic_call,
            "conserved_adjacency_score": round(self.conserved_adjacency_score, 4),
            "spearman_rho": round(self.spearman_rho, 4),
            "kendall_tau": round(self.kendall_tau, 4),
            "breakpoint_count": self.breakpoint_count,
            "blocks": [b.to_dict() for b in self.blocks],
        }


class SyntenyAnalyzer:
    """Core comparative genomics synteny engine."""

    @staticmethod
    def compute_ani(blocks: List[SyntenyBlock], ref_length: int, qry_length: int) -> Tuple[float, float, float, str]:
        """
        Compute weighted Average Nucleotide Identity (ANI) and coverage fractions.
        Returns: (ani_pct, ref_coverage_pct, qry_coverage_pct, taxonomic_call)
        """
        if not blocks or ref_length <= 0 or qry_length <= 0:
            return 0.0, 0.0, 0.0, "No alignment"

        total_ref_aligned = sum(b.ref_length for b in blocks)
        total_qry_aligned = sum(b.qry_length for b in blocks)

        if total_ref_aligned == 0:
            return 0.0, 0.0, 0.0, "No alignment"

        weighted_id = sum(b.ref_length * b.identity_pct for b in blocks) / total_ref_aligned
        ref_cov = min(100.0, (total_ref_aligned / ref_length) * 100.0)
        qry_cov = min(100.0, (total_qry_aligned / qry_length) * 100.0)

        # Taxonomic species boundary standards (Richter & Rossello-Mora 2009)
        if weighted_id >= 95.0 and ref_cov >= 60.0:
            tax_call = "Same Species (Conspecific)"
        elif weighted_id >= 93.0:
            tax_call = "Subspecies / Borderline Species Divergence"
        elif weighted_id >= 85.0:
            tax_call = "Different Species (Same Genus)"
        else:
            tax_call = "Divergent Taxon / Different Genus"

        return round(weighted_id, 2), round(ref_cov, 2), round(qry_cov, 2), tax_call

    @staticmethod
    def compute_spearman_rho(ref_order: List[int], qry_order: List[int]) -> float:
        """
        Spearman's rank correlation coefficient on ortholog gene positions:
        rho = 1 - (6 * sum(d_i^2)) / (n * (n^2 - 1))
        """
        n = len(ref_order)
        if n < 2 or len(qry_order) != n:
            return 1.0 if n <= 1 else 0.0
        d_sq_sum = sum((r - q) ** 2 for r, q in zip(ref_order, qry_order))
        denom = n * (n ** 2 - 1)
        if denom == 0:
            return 1.0
        rho = 1.0 - (6.0 * d_sq_sum) / denom
        return max(-1.0, min(1.0, rho))

    @staticmethod
    def compute_kendall_tau(ref_order: List[int], qry_order: List[int]) -> float:
        """
        Kendall's tau-b rank correlation between reference and query orders.
        tau = (concordant_pairs - discordant_pairs) / (n * (n - 1) / 2)
        """
        n = len(ref_order)
        if n < 2:
            return 1.0
        concordant = 0
        discordant = 0
        for i in range(n):
            for j in range(i + 1, n):
                r_diff = ref_order[i] - ref_order[j]
                q_diff = qry_order[i] - qry_order[j]
                prod = r_diff * q_diff
                if prod > 0:
                    concordant += 1
                elif prod < 0:
                    discordant += 1
        total_pairs = n * (n - 1) / 2.0
        if total_pairs == 0:
            return 1.0
        return (concordant - discordant) / total_pairs

    @staticmethod
    def compute_conserved_adjacency(ref_gene_list: List[str], qry_gene_list: List[str]) -> Tuple[float, int]:
        """
        Synteny Conservation Score (SCS) / Conserved Adjacency Fraction:
        Calculates fraction of adjacent gene pairs (A, B) in reference that remain
        adjacent (either forward (A, B) or inverted (B, A)) in query.
        Returns: (fraction_conserved, breakpoint_count)
        """
        if len(ref_gene_list) < 2 or len(qry_gene_list) < 2:
            return 1.0, 0

        # Build query adjacency set
        qry_adjacencies = set()
        for i in range(len(qry_gene_list) - 1):
            g1, g2 = qry_gene_list[i], qry_gene_list[i + 1]
            qry_adjacencies.add((g1, g2))
            qry_adjacencies.add((g2, g1))

        ref_total_pairs = len(ref_gene_list) - 1
        conserved_count = 0
        breakpoints = 0

        for i in range(len(ref_gene_list) - 1):
            g1, g2 = ref_gene_list[i], ref_gene_list[i + 1]
            if (g1, g2) in qry_adjacencies:
                conserved_count += 1
            else:
                breakpoints += 1

        score = conserved_count / ref_total_pairs if ref_total_pairs > 0 else 1.0
        return round(score, 4), breakpoints

    @classmethod
    def analyze(
        cls,
        reference_name: str,
        query_name: str,
        ref_length: int,
        qry_length: int,
        blocks: List[SyntenyBlock],
    ) -> SyntenyAnalysisResult:
        """Perform comprehensive synteny analysis."""
        ani_pct, ref_cov, qry_cov, tax_call = cls.compute_ani(blocks, ref_length, qry_length)

        collinear = sum(1 for b in blocks if not b.is_inverted)
        inverted = sum(1 for b in blocks if b.is_inverted)

        # Collect gene orders across blocks
        ref_genes = []
        qry_genes = []
        for b in sorted(blocks, key=lambda x: x.ref_start):
            ref_genes.extend(b.ref_genes)
        for b in sorted(blocks, key=lambda x: min(x.qry_start, x.qry_end)):
            qry_genes.extend(b.qry_genes)

        # Gene orders
        if ref_genes and qry_genes:
            common = [g for g in ref_genes if g in qry_genes]
            if len(common) >= 2:
                ref_rank = list(range(len(common)))
                qry_pos_map = {g: idx for idx, g in enumerate(qry_genes)}
                qry_rank = [qry_pos_map[g] for g in common]
                rho = cls.compute_spearman_rho(ref_rank, qry_rank)
                tau = cls.compute_kendall_tau(ref_rank, qry_rank)
                cas, bp_count = cls.compute_conserved_adjacency(common, [g for g in qry_genes if g in common])
            else:
                rho, tau, cas, bp_count = (1.0 if not inverted else 0.0, 1.0 if not inverted else 0.0, 1.0, 0)
        else:
            # Fallback based on block order
            sorted_by_ref = sorted(blocks, key=lambda x: x.ref_start)
            ref_rank = list(range(len(sorted_by_ref)))
            sorted_by_qry = sorted(blocks, key=lambda x: min(x.qry_start, x.qry_end))
            qry_pos_map = {b.block_id: idx for idx, b in enumerate(sorted_by_qry)}
            qry_rank = [qry_pos_map[b.block_id] for b in sorted_by_ref]
            rho = cls.compute_spearman_rho(ref_rank, qry_rank)
            tau = cls.compute_kendall_tau(ref_rank, qry_rank)
            block_ids_ref = [b.block_id for b in sorted_by_ref]
            block_ids_qry = [b.block_id for b in sorted_by_qry]
            cas, bp_count = cls.compute_conserved_adjacency(block_ids_ref, block_ids_qry)

        return SyntenyAnalysisResult(
            reference_name=reference_name,
            query_name=query_name,
            ref_genome_length=ref_length,
            qry_genome_length=qry_length,
            num_blocks=len(blocks),
            collinear_blocks=collinear,
            inverted_blocks=inverted,
            ani_pct=ani_pct,
            ref_coverage_pct=ref_cov,
            qry_coverage_pct=qry_cov,
            taxonomic_call=tax_call,
            conserved_adjacency_score=cas,
            spearman_rho=rho,
            kendall_tau=tau,
            breakpoint_count=bp_count,
            blocks=blocks,
        )


class DotplotRenderer:
    """ASCII 2D Comparative Synteny Dotplot Visualizer."""

    @staticmethod
    def render(blocks: List[SyntenyBlock], ref_length: int, qry_length: int, grid_size: int = 30) -> str:
        """
        Renders an ASCII dotplot where:
        '\\' represents collinear forward alignment
        '/' represents inverted reverse alignment
        '.' represents unaligned matrix space
        """
        grid = [["." for _ in range(grid_size)] for _ in range(grid_size)]

        for b in blocks:
            r0 = max(0, min(grid_size - 1, int((b.ref_start / ref_length) * grid_size)))
            r1 = max(0, min(grid_size - 1, int((b.ref_end / ref_length) * grid_size)))
            if r1 < r0:
                r0, r1 = r1, r0
            r1 = max(r1, r0 + 1)

            if b.is_inverted:
                q0 = max(0, min(grid_size - 1, int((abs(b.qry_end) / qry_length) * grid_size)))
                q1 = max(0, min(grid_size - 1, int((abs(b.qry_start) / qry_length) * grid_size)))
                if q1 < q0:
                    q0, q1 = q1, q0
                q1 = max(q1, q0 + 1)
            else:
                q0 = max(0, min(grid_size - 1, int((b.qry_start / qry_length) * grid_size)))
                q1 = max(0, min(grid_size - 1, int((b.qry_end / qry_length) * grid_size)))
                if q1 < q0:
                    q0, q1 = q1, q0
                q1 = max(q1, q0 + 1)

            r_span = list(range(r0, min(r1, grid_size)))
            q_span = list(range(q0, min(q1, grid_size)))
            steps = max(len(r_span), len(q_span))

            for i in range(steps):
                ri = r_span[min(i, len(r_span) - 1)]
                if not b.is_inverted:
                    qi = q_span[min(i, len(q_span) - 1)]
                    grid[ri][qi] = "\\"
                else:
                    qi = q_span[len(q_span) - 1 - min(i, len(q_span) - 1)]
                    grid[ri][qi] = "/"

        header = "      " + "".join(str(i // 10 % 10) for i in range(grid_size))
        lines = [header]
        for idx, row in enumerate(grid):
            lines.append(f"{idx:4d}  " + "".join(row))
        return "\n".join(lines)


class PrimerDesigner:
    """Breakpoint Junction PCR Primer Design for Inversion Validation."""

    @staticmethod
    def calculate_tm(seq: str) -> float:
        """Salt-adjusted nearest-neighbor / Wallace melting temperature."""
        seq_u = seq.upper()
        gc = seq_u.count("G") + seq_u.count("C")
        at = len(seq_u) - gc
        if len(seq_u) < 14:
            return float(2 * at + 4 * gc)
        return round(64.9 + 41.0 * (gc - 16.4) / len(seq_u), 2)

    @staticmethod
    def calculate_gc(seq: str) -> float:
        """Percentage of Guanine and Cytosine bases."""
        if not seq:
            return 0.0
        seq_u = seq.upper()
        gc = seq_u.count("G") + seq_u.count("C")
        return round((gc / len(seq_u)) * 100.0, 2)

    @classmethod
    def design_flanking_primers(
        cls,
        left_flank_seq: str,
        right_flank_seq: str,
        primer_len: int = 20,
        target_amplicon_bp: int = 400
    ) -> Dict[str, Any]:
        """
        Designs forward primer on left flank and reverse primer on right flank
        spanning an inversion or rearrangement breakpoint.
        """
        clean_left = re.sub(r'[^ACGTacgt]', '', left_flank_seq).upper()
        clean_right = re.sub(r'[^ACGTacgt]', '', right_flank_seq).upper()

        if len(clean_left) < primer_len or len(clean_right) < primer_len:
            raise ValueError(f"Flank sequences must be at least {primer_len} bp in length.")

        fwd_seq = clean_left[-primer_len:]
        # Reverse complement of 5' end of right flank
        rc_trans = str.maketrans("ACGT", "TGCA")
        rev_seq = clean_right[:primer_len].translate(rc_trans)[::-1]

        fwd_tm = cls.calculate_tm(fwd_seq)
        fwd_gc = cls.calculate_gc(fwd_seq)
        fwd_has_clamp = fwd_seq[-1] in ("G", "C")

        rev_tm = cls.calculate_tm(rev_seq)
        rev_gc = cls.calculate_gc(rev_seq)
        rev_has_clamp = rev_seq[-1] in ("G", "C")

        fwd_valid = (40.0 <= fwd_gc <= 60.0) and (55.0 <= fwd_tm <= 68.0)
        rev_valid = (40.0 <= rev_gc <= 60.0) and (55.0 <= rev_tm <= 68.0)
        tm_diff = abs(fwd_tm - rev_tm)

        return {
            "forward_primer": {
                "sequence": fwd_seq,
                "length": len(fwd_seq),
                "tm_c": fwd_tm,
                "gc_pct": fwd_gc,
                "gc_clamp": fwd_has_clamp,
                "valid": fwd_valid,
            },
            "reverse_primer": {
                "sequence": rev_seq,
                "length": len(rev_seq),
                "tm_c": rev_tm,
                "gc_pct": rev_gc,
                "gc_clamp": rev_has_clamp,
                "valid": rev_valid,
            },
            "tm_difference": round(tm_diff, 2),
            "amplicon_size_bp": target_amplicon_bp,
            "overall_qc_passed": fwd_valid and rev_valid and (tm_diff <= 3.0),
            "purpose": "PCR amplification across genomic breakpoint junction for structural verification",
        }


class PangenomeAnalyzer:
    """Pangenome Partitioning (Core, Soft Core, Shell, Cloud) across bacterial isolates."""

    @staticmethod
    def partition(gene_presence_matrix: Dict[str, Set[str]]) -> Dict[str, Any]:
        """
        Partitions gene families based on frequency of occurrence across strain collection:
        - Core Genome: Present in >= 99% of strains (or 100%)
        - Soft Core: Present in 95% - 99% of strains
        - Shell Genome: Present in 15% - 95% of strains (flexible accessory genome)
        - Cloud Genome: Present in < 15% of strains (strain-specific rare genes)
        """
        all_strains: Set[str] = set()
        for strains in gene_presence_matrix.values():
            all_strains |= strains

        n_strains = len(all_strains)
        if n_strains == 0:
            return {"total_families": 0, "n_strains": 0, "core": [], "soft_core": [], "shell": [], "cloud": []}

        core = []
        soft_core = []
        shell = []
        cloud = []

        for family_id, strains in sorted(gene_presence_matrix.items()):
            freq = len(strains) / n_strains
            if freq >= 0.99:
                core.append(family_id)
            elif freq >= 0.95:
                soft_core.append(family_id)
            elif freq >= 0.15:
                shell.append(family_id)
            else:
                cloud.append(family_id)

        total = len(gene_presence_matrix)
        core_pct = round((len(core) / total) * 100.0, 2) if total else 0.0
        shell_pct = round((len(shell) / total) * 100.0, 2) if total else 0.0
        cloud_pct = round((len(cloud) / total) * 100.0, 2) if total else 0.0

        return {
            "total_families": total,
            "n_strains": n_strains,
            "strains": sorted(list(all_strains)),
            "core_count": len(core),
            "soft_core_count": len(soft_core),
            "shell_count": len(shell),
            "cloud_count": len(cloud),
            "core_pct": core_pct,
            "shell_pct": shell_pct,
            "cloud_pct": cloud_pct,
            "core_families": core,
            "soft_core_families": soft_core,
            "shell_families": shell,
            "cloud_families": cloud,
            "pangenome_type": "Open (High Diversity)" if (cloud_pct + shell_pct) > 50.0 else "Closed (Conserved)",
        }


def get_curated_synteny_benchmarks() -> Dict[str, Dict[str, Any]]:
    """Curated bacterial synteny benchmark pairs."""
    return {
        "ecoli_k12_vs_o157": {
            "ref_name": "Escherichia coli K-12 MG1655",
            "qry_name": "Escherichia coli O157:H7 EDL933",
            "ref_length": 4641652,
            "qry_length": 5528445,
            "blocks": [
                SyntenyBlock("BLK_1", 0, 1_200_000, 0, 1_350_000, "+", 99.1,
                             ["dnaA", "dnaN", "gyrB", "recF", "gyrA"], ["dnaA", "dnaN", "gyrB", "recF", "gyrA"]),
                SyntenyBlock("BLK_2_INV", 1_200_000, 2_400_000, 3_800_000, 2_500_000, "-", 98.4,
                             ["lacZ", "lacY", "cysK", "trpA", "trpB"], ["trpB", "trpA", "cysK", "lacY", "lacZ"]),
                SyntenyBlock("BLK_3", 2_400_000, 3_800_000, 3_800_000, 5_000_000, "+", 98.9,
                             ["hisA", "hisB", "rpoB", "rpoC", "atpA"], ["hisA", "hisB", "rpoB", "rpoC", "atpA"]),
                SyntenyBlock("BLK_4", 3_800_000, 4_641_652, 5_000_000, 5_528_445, "+", 98.2,
                             ["ftsZ", "groEL", "dnaK", "tufA"], ["ftsZ", "groEL", "dnaK", "tufA"]),
            ]
        },
        "salmonella_lt2_vs_ct18": {
            "ref_name": "Salmonella enterica serovar Typhimurium LT2",
            "qry_name": "Salmonella enterica serovar Typhi CT18",
            "ref_length": 4857432,
            "qry_length": 4809037,
            "blocks": [
                SyntenyBlock("SB_1", 0, 1_500_000, 0, 1_490_000, "+", 99.4,
                             ["thrA", "thrB", "carA", "dnaJ"], ["thrA", "thrB", "carA", "dnaJ"]),
                SyntenyBlock("SB_2_INV", 1_500_000, 3_100_000, 3_120_000, 1_520_000, "-", 99.1,
                             ["invA", "hilA", "prgH", "sseA"], ["sseA", "prgH", "hilA", "invA"]),
                SyntenyBlock("SB_3", 3_100_000, 4_857_432, 3_120_000, 4_809_037, "+", 99.3,
                             ["viaB", "sopB", "fliC", "clpB"], ["viaB", "sopB", "fliC", "clpB"]),
            ]
        },
        "pseudomonas_pao1_vs_lesb58": {
            "ref_name": "Pseudomonas aeruginosa PAO1",
            "qry_name": "Pseudomonas aeruginosa LESB58",
            "ref_length": 6264404,
            "qry_length": 6601757,
            "blocks": [
                SyntenyBlock("PBLK_1", 0, 2_000_000, 0, 2_100_000, "+", 98.8,
                             ["dnaA", "oprF", "algD"], ["dnaA", "oprF", "algD"]),
                SyntenyBlock("PBLK_2", 2_000_000, 4_200_000, 2_100_000, 4_400_000, "+", 98.5,
                             ["toxA", "lasB", "pcrV"], ["toxA", "lasB", "pcrV"]),
                SyntenyBlock("PBLK_3_INV", 4_200_000, 5_500_000, 6_000_000, 4_700_000, "-", 97.9,
                             ["pyoS2", "mexA", "mexB"], ["mexB", "mexA", "pyoS2"]),
                SyntenyBlock("PBLK_4", 5_500_000, 6_264_404, 6_000_000, 6_601_757, "+", 98.1,
                             ["clpP", "rpoD"], ["clpP", "rpoD"]),
            ]
        }
    }


def parse_gene_list(val: str) -> List[str]:
    """Parse comma-separated or semicolon-separated gene lists from CSV field."""
    if not val or not val.strip():
        return []
    clean = val.strip().strip('"').strip("'")
    if not clean:
        return []
    delimiter = ";" if ";" in clean else ","
    return [g.strip() for g in clean.split(delimiter) if g.strip()]


def process_synteny_batch_csv(input_csv_path: str, output_csv_path: str) -> List[Dict[str, Any]]:
    """
    Process batch CSV file containing comparative bacterial genomics parameters.
    Supports either:
    1) Per-comparison pairs (grouped rows by pair_id / ref_name+qry_name or single row per block),
       aggregating blocks and computing full synteny metrics (ANI, conserved adjacency score,
       Spearman's rho, Kendall's tau, breakpoint count, collinear/inverted block counts).
    2) Detailed output with both block-level and genome-level synteny metrics.

    Writes output CSV with enriched comparative metrics and returns list of result summaries.
    """
    import csv

    with open(input_csv_path, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"Input CSV file '{input_csv_path}' is empty or invalid.")
        rows = list(reader)

    if not rows:
        raise ValueError(f"No records found in CSV file '{input_csv_path}'.")

    # Group rows by comparison pair (pair_id or reference_name+query_name)
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for idx, row in enumerate(rows):
        pair_id = row.get("pair_id") or row.get("comparison_id") or f"{row.get('reference_name', 'ref')}_vs_{row.get('query_name', 'qry')}"
        if not pair_id:
            pair_id = f"comparison_{idx + 1}"
        if pair_id not in grouped:
            grouped[pair_id] = []
        grouped[pair_id].append(row)

    results: List[Dict[str, Any]] = []

    for pair_id, pair_rows in grouped.items():
        first = pair_rows[0]
        ref_name = first.get("reference_name") or first.get("ref_name") or "Reference_Genome"
        qry_name = first.get("query_name") or first.get("qry_name") or "Query_Genome"

        try:
            ref_len = int(float(first.get("ref_genome_length") or first.get("ref_length") or 5_000_000))
        except (ValueError, TypeError):
            ref_len = 5_000_000

        try:
            qry_len = int(float(first.get("qry_genome_length") or first.get("qry_length") or 5_000_000))
        except (ValueError, TypeError):
            qry_len = 5_000_000

        blocks: List[SyntenyBlock] = []
        for b_idx, r in enumerate(pair_rows):
            block_id = r.get("block_id") or f"BLK_{b_idx + 1}"
            try:
                r_start = int(float(r.get("ref_start", 0)))
                r_end = int(float(r.get("ref_end", ref_len)))
                q_start = int(float(r.get("qry_start", 0)))
                q_end = int(float(r.get("qry_end", qry_len)))
            except (ValueError, TypeError):
                r_start, r_end, q_start, q_end = 0, ref_len, 0, qry_len

            strand = r.get("strand", "+").strip()
            if strand not in ("+", "-"):
                strand = "-" if "inv" in block_id.lower() or "inv" in r.get("event_type", "").lower() else "+"

            try:
                ident = float(r.get("identity_pct", 98.0))
            except (ValueError, TypeError):
                ident = 98.0

            ref_genes = parse_gene_list(r.get("ref_genes", ""))
            qry_genes = parse_gene_list(r.get("qry_genes", ""))

            # If no gene list is provided, infer from ortholog pair fields
            if not ref_genes and r.get("ref_gene"):
                ref_genes = [r["ref_gene"]]
            if not qry_genes and r.get("qry_gene"):
                qry_genes = [r["qry_gene"]]

            b = SyntenyBlock(
                block_id=block_id,
                ref_start=r_start,
                ref_end=r_end,
                qry_start=q_start,
                qry_end=q_end,
                strand=strand,
                identity_pct=ident,
                ref_genes=ref_genes,
                qry_genes=qry_genes,
            )
            blocks.append(b)

        analysis = SyntenyAnalyzer.analyze(ref_name, qry_name, ref_len, qry_len, blocks)

        # Detect genomic islands / horizontal gene transfer candidate regions
        # Heuristic: segments with >15% length disparity or reduced sequence identity (<90%)
        island_candidates = 0
        for b in blocks:
            len_ratio = min(b.ref_length, b.qry_length) / max(b.ref_length, b.qry_length) if max(b.ref_length, b.qry_length) > 0 else 1.0
            if b.identity_pct < 92.0 or len_ratio < 0.70:
                island_candidates += 1

        res_dict = {
            "pair_id": pair_id,
            "reference_name": analysis.reference_name,
            "query_name": analysis.query_name,
            "ref_genome_length": analysis.ref_genome_length,
            "qry_genome_length": analysis.qry_genome_length,
            "num_blocks": analysis.num_blocks,
            "collinear_blocks": analysis.collinear_blocks,
            "inverted_blocks": analysis.inverted_blocks,
            "ani_pct": analysis.ani_pct,
            "ref_coverage_pct": analysis.ref_coverage_pct,
            "qry_coverage_pct": analysis.qry_coverage_pct,
            "taxonomic_call": analysis.taxonomic_call,
            "conserved_adjacency_score": analysis.conserved_adjacency_score,
            "spearman_rho": analysis.spearman_rho,
            "kendall_tau": analysis.kendall_tau,
            "breakpoint_count": analysis.breakpoint_count,
            "genomic_island_candidates": island_candidates,
            "synteny_status": "CONSERVED_SYNTENY" if analysis.conserved_adjacency_score >= 0.70 and analysis.inverted_blocks == 0
                              else ("REARRANGED_INVERSIONS" if analysis.inverted_blocks > 0 else "PARTIALLY_CONSERVED"),
        }
        results.append(res_dict)

    fieldnames = [
        "pair_id",
        "reference_name",
        "query_name",
        "ref_genome_length",
        "qry_genome_length",
        "num_blocks",
        "collinear_blocks",
        "inverted_blocks",
        "ani_pct",
        "ref_coverage_pct",
        "qry_coverage_pct",
        "taxonomic_call",
        "conserved_adjacency_score",
        "spearman_rho",
        "kendall_tau",
        "breakpoint_count",
        "genomic_island_candidates",
        "synteny_status",
    ]

    with open(output_csv_path, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    return results

