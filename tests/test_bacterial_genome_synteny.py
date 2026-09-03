#!/usr/bin/env python3
"""
Unit Test Suite for Bacterial Genome Synteny Agent.
Covers:
- SyntenyBlock orientation and span calculations
- Weighted Average Nucleotide Identity (ANI) & taxonomic thresholds
- Gene order collinearity: Spearman's rank rho, Kendall's tau
- Conserved adjacency score & chromosomal breakpoint counting
- 2D ASCII Dotplot matrix rendering
- Inversion breakpoint PCR primer design (Tm, GC, GC clamp)
- Pangenome partitioning (Core, Soft Core, Shell, Cloud)
- Curated bacterial benchmark pairs
- CLI interface and JSON output parsing
"""

import json
import math
import os
import sys
import unittest
from io import StringIO
from unittest.mock import patch

from bacterial_genome_synteny import (
    SyntenyBlock,
    SyntenyAnalyzer,
    DotplotRenderer,
    PrimerDesigner,
    PangenomeAnalyzer,
    get_curated_synteny_benchmarks,
)
import cli


class TestSyntenyBlockAndANI(unittest.TestCase):
    def test_synteny_block_forward(self):
        b = SyntenyBlock("B1", 0, 1000, 0, 1000, "+", 99.0)
        self.assertEqual(b.ref_length, 1000)
        self.assertEqual(b.qry_length, 1000)
        self.assertFalse(b.is_inverted)

    def test_synteny_block_inverted_strand(self):
        b = SyntenyBlock("B2", 1000, 2000, 2000, 1000, "-", 98.0)
        self.assertEqual(b.ref_length, 1000)
        self.assertEqual(b.qry_length, 1000)
        self.assertTrue(b.is_inverted)

    def test_ani_perfect_match(self):
        blocks = [
            SyntenyBlock("B1", 0, 1_000_000, 0, 1_000_000, "+", 100.0),
            SyntenyBlock("B2", 1_000_000, 2_000_000, 1_000_000, 2_000_000, "+", 100.0),
        ]
        ani, ref_cov, qry_cov, call = SyntenyAnalyzer.compute_ani(blocks, 2_000_000, 2_000_000)
        self.assertEqual(ani, 100.0)
        self.assertEqual(ref_cov, 100.0)
        self.assertEqual(qry_cov, 100.0)
        self.assertEqual(call, "Same Species (Conspecific)")

    def test_ani_species_cutoff_thresholds(self):
        # 96.0% ANI -> Same Species
        b_same = [SyntenyBlock("B1", 0, 1_000_000, 0, 1_000_000, "+", 96.0)]
        ani1, _, _, call1 = SyntenyAnalyzer.compute_ani(b_same, 1_000_000, 1_000_000)
        self.assertEqual(call1, "Same Species (Conspecific)")

        # 93.5% ANI -> Subspecies / Borderline
        b_sub = [SyntenyBlock("B1", 0, 1_000_000, 0, 1_000_000, "+", 93.5)]
        ani2, _, _, call2 = SyntenyAnalyzer.compute_ani(b_sub, 1_000_000, 1_000_000)
        self.assertEqual(call2, "Subspecies / Borderline Species Divergence")

        # 88.0% ANI -> Different Species
        b_diff = [SyntenyBlock("B1", 0, 1_000_000, 0, 1_000_000, "+", 88.0)]
        ani3, _, _, call3 = SyntenyAnalyzer.compute_ani(b_diff, 1_000_000, 1_000_000)
        self.assertEqual(call3, "Different Species (Same Genus)")

        # 80.0% ANI -> Different Genus
        b_genus = [SyntenyBlock("B1", 0, 1_000_000, 0, 1_000_000, "+", 80.0)]
        ani4, _, _, call4 = SyntenyAnalyzer.compute_ani(b_genus, 1_000_000, 1_000_000)
        self.assertEqual(call4, "Divergent Taxon / Different Genus")

    def test_ani_empty_blocks(self):
        ani, ref_cov, qry_cov, call = SyntenyAnalyzer.compute_ani([], 1_000_000, 1_000_000)
        self.assertEqual(ani, 0.0)
        self.assertEqual(call, "No alignment")


class TestCollinearityAndSyntenyMetrics(unittest.TestCase):
    def test_spearman_rho_perfect_order(self):
        ref_order = [0, 1, 2, 3, 4]
        qry_order = [0, 1, 2, 3, 4]
        rho = SyntenyAnalyzer.compute_spearman_rho(ref_order, qry_order)
        self.assertAlmostEqual(rho, 1.0, places=4)

    def test_spearman_rho_inverted_order(self):
        ref_order = [0, 1, 2, 3, 4]
        qry_order = [4, 3, 2, 1, 0]
        rho = SyntenyAnalyzer.compute_spearman_rho(ref_order, qry_order)
        self.assertAlmostEqual(rho, -1.0, places=4)

    def test_kendall_tau_perfect_order(self):
        ref_order = [0, 1, 2, 3]
        qry_order = [0, 1, 2, 3]
        tau = SyntenyAnalyzer.compute_kendall_tau(ref_order, qry_order)
        self.assertAlmostEqual(tau, 1.0, places=4)

    def test_kendall_tau_inverted_order(self):
        ref_order = [0, 1, 2, 3]
        qry_order = [3, 2, 1, 0]
        tau = SyntenyAnalyzer.compute_kendall_tau(ref_order, qry_order)
        self.assertAlmostEqual(tau, -1.0, places=4)

    def test_conserved_adjacency_score(self):
        ref_genes = ["A", "B", "C", "D", "E"]
        qry_genes = ["A", "B", "C", "D", "E"]
        score, bp = SyntenyAnalyzer.compute_conserved_adjacency(ref_genes, qry_genes)
        self.assertEqual(score, 1.0)
        self.assertEqual(bp, 0)

    def test_conserved_adjacency_with_inversion(self):
        # Inversion of (C, D)
        ref_genes = ["A", "B", "C", "D", "E"]
        qry_genes = ["A", "B", "D", "C", "E"]
        # Pairs in ref: (A,B), (B,C), (C,D), (D,E)
        # Adjacencies in qry: (A,B), (B,D), (D,C)=(C,D), (C,E)
        # Conserved: (A,B) [yes], (B,C) [no], (C,D) [yes], (D,E) [no] -> 2 / 4 = 0.5
        score, bp = SyntenyAnalyzer.compute_conserved_adjacency(ref_genes, qry_genes)
        self.assertEqual(score, 0.5)
        self.assertEqual(bp, 2)


class TestDotplotAndPrimers(unittest.TestCase):
    def test_dotplot_renderer(self):
        blocks = [
            SyntenyBlock("B1", 0, 1000, 0, 1000, "+", 99.0),
            SyntenyBlock("B2", 1000, 2000, 2000, 1000, "-", 98.0),
        ]
        ascii_plot = DotplotRenderer.render(blocks, 2000, 2000, grid_size=10)
        self.assertIn("\\", ascii_plot)
        self.assertIn("/", ascii_plot)
        self.assertTrue(len(ascii_plot.split("\n")) >= 10)

    def test_primer_designer_tm_calculation(self):
        # 20-mer: 10 GC, 10 AT
        seq = "GCGCGCGCGCAAAAAATTTT"
        gc = PrimerDesigner.calculate_gc(seq)
        tm = PrimerDesigner.calculate_tm(seq)
        self.assertEqual(gc, 50.0)
        self.assertTrue(50.0 <= tm <= 65.0)

    def test_flanking_primers_design(self):
        lf = "ATCGATCGATCGATCGATCG" + "GGCACCGATTCAGGCACTAG"
        rf = "CTGTAGCGTTGGCATTGAAG" + "ATCGATCGATCGATCGATCG"
        res = PrimerDesigner.design_flanking_primers(lf, rf, primer_len=20)
        self.assertIn("forward_primer", res)
        self.assertIn("reverse_primer", res)
        self.assertEqual(len(res["forward_primer"]["sequence"]), 20)
        self.assertEqual(len(res["reverse_primer"]["sequence"]), 20)
        self.assertTrue(res["forward_primer"]["gc_clamp"])
        self.assertTrue(res["reverse_primer"]["gc_clamp"])

    def test_primer_short_sequence_error(self):
        with self.assertRaises(ValueError):
            PrimerDesigner.design_flanking_primers("ACTG", "ACTG", primer_len=20)


class TestPangenomePartitioning(unittest.TestCase):
    def test_pangenome_core_shell_cloud(self):
        matrix = {
            "core_gene_1": {"S1", "S2", "S3", "S4"},
            "core_gene_2": {"S1", "S2", "S3", "S4"},
            "shell_gene_1": {"S1", "S2"},
            "cloud_gene_1": {"S1"},  # 1/4 = 25% is shell under 15-95%
        }
        res = PangenomeAnalyzer.partition(matrix)
        self.assertEqual(res["n_strains"], 4)
        self.assertEqual(res["core_count"], 2)
        self.assertIn("core_gene_1", res["core_families"])
        self.assertIn("core_gene_2", res["core_families"])

    def test_pangenome_empty_matrix(self):
        res = PangenomeAnalyzer.partition({})
        self.assertEqual(res["total_families"], 0)
        self.assertEqual(res["n_strains"], 0)


class TestBenchmarksAndCLI(unittest.TestCase):
    def test_curated_benchmarks_available(self):
        benchmarks = get_curated_synteny_benchmarks()
        self.assertIn("ecoli_k12_vs_o157", benchmarks)
        self.assertIn("salmonella_lt2_vs_ct18", benchmarks)
        self.assertIn("pseudomonas_pao1_vs_lesb58", benchmarks)

    def test_full_analysis_ecoli(self):
        b = get_curated_synteny_benchmarks()["ecoli_k12_vs_o157"]
        res = SyntenyAnalyzer.analyze(b["ref_name"], b["qry_name"], b["ref_length"], b["qry_length"], b["blocks"])
        self.assertEqual(res.collinear_blocks, 3)
        self.assertEqual(res.inverted_blocks, 1)
        self.assertTrue(res.ani_pct >= 98.0)
        self.assertEqual(res.taxonomic_call, "Same Species (Conspecific)")

    def test_cli_list_benchmarks(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["--list-benchmarks"])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("ecoli_k12_vs_o157", output)

    def test_cli_benchmark_analysis(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["--benchmark", "ecoli_k12_vs_o157"])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("BACTERIAL GENOME SYNTENY REPORT", output)
            self.assertIn("Escherichia coli K-12 MG1655", output)

    def test_cli_benchmark_dotplot(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["--benchmark", "ecoli_k12_vs_o157", "--dotplot"])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("ASCII Comparative Synteny Dotplot", output)
            self.assertIn("\\", output)
            self.assertIn("/", output)

    def test_cli_benchmark_json(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["--benchmark", "salmonella_lt2_vs_ct18", "--json"])
            self.assertEqual(ret, 0)
            data = json.loads(fake_out.getvalue())
            self.assertEqual(data["reference_name"], "Salmonella enterica serovar Typhimurium LT2")
            self.assertTrue(data["ani_pct"] >= 99.0)

    def test_cli_pangenome(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["--pangenome"])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("Pangenome Partitioning Analysis", output)

    def test_cli_primers(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["--primers"])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("Inversion Breakpoint Validation PCR Primers", output)

    def test_cli_invalid_benchmark(self):
        with patch('sys.stderr', new=StringIO()) as fake_err:
            ret = cli.main(["--benchmark", "invalid_species_xyz"])
            self.assertEqual(ret, 1)


class TestBatchProcessingAndSubcommands(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_csv = os.path.join(self.temp_dir.name, "test_input.csv")
        self.output_csv = os.path.join(self.temp_dir.name, "test_output.csv")

        # Create test CSV
        with open(self.input_csv, "w", encoding="utf-8") as f:
            f.write(
                "pair_id,reference_name,query_name,ref_genome_length,qry_genome_length,block_id,ref_start,ref_end,qry_start,qry_end,strand,identity_pct,ref_genes,qry_genes,event_type\n"
                "TEST_PAIR_1,E. coli K12,E. coli O157,4600000,5500000,BLK1,0,2000000,0,2100000,+,99.2,dnaA;gyrB,dnaA;gyrB,collinear\n"
                "TEST_PAIR_1,E. coli K12,E. coli O157,4600000,5500000,BLK2,2000000,4600000,5000000,2400000,-,98.5,lacZ;trpA,trpA;lacZ,inversion\n"
                "TEST_PAIR_2,S. enterica LT2,S. enterica CT18,4800000,4800000,BLK1,0,4800000,0,4800000,+,99.5,invA;sseA,invA;sseA,collinear\n"
            )

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_process_synteny_batch_csv_function(self):
        from bacterial_genome_synteny import process_synteny_batch_csv
        res = process_synteny_batch_csv(self.input_csv, self.output_csv)
        self.assertEqual(len(res), 2)
        self.assertEqual(res[0]["pair_id"], "TEST_PAIR_1")
        self.assertEqual(res[0]["inverted_blocks"], 1)
        self.assertEqual(res[0]["collinear_blocks"], 1)
        self.assertEqual(res[1]["pair_id"], "TEST_PAIR_2")
        self.assertEqual(res[1]["inverted_blocks"], 0)
        self.assertTrue(os.path.exists(self.output_csv))

    def test_cli_batch_subcommand(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["batch", "--input", self.input_csv, "--output", self.output_csv])
            self.assertEqual(ret, 0)
            self.assertIn("BATCH COMPARATIVE SYNTENY PROCESSING COMPLETE", fake_out.getvalue())
            self.assertTrue(os.path.exists(self.output_csv))

    def test_cli_batch_subcommand_json(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["batch", "-i", self.input_csv, "-o", self.output_csv, "--json"])
            self.assertEqual(ret, 0)
            data = json.loads(fake_out.getvalue())
            self.assertEqual(data["status"], "SUCCESS")
            self.assertEqual(data["comparisons_processed"], 2)

    def test_cli_root_batch_flags(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["-i", self.input_csv, "-o", self.output_csv])
            self.assertEqual(ret, 0)
            self.assertIn("BATCH COMPARATIVE SYNTENY PROCESSING COMPLETE", fake_out.getvalue())

    def test_cli_benchmark_subcommand(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["benchmark", "ecoli_k12_vs_o157", "--dotplot"])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("BACTERIAL GENOME SYNTENY REPORT", output)
            self.assertIn("ASCII Comparative Synteny Dotplot", output)

    def test_cli_pangenome_subcommand(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["pangenome", "--json"])
            self.assertEqual(ret, 0)
            data = json.loads(fake_out.getvalue())
            self.assertIn("total_families", data)

    def test_cli_primers_subcommand(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["primers", "--json"])
            self.assertEqual(ret, 0)
            data = json.loads(fake_out.getvalue())
            self.assertTrue(data["overall_qc_passed"])

    def test_cli_batch_nonexistent_file(self):
        with patch('sys.stderr', new=StringIO()) as fake_err:
            ret = cli.main(["batch", "-i", "nonexistent_file_xyz.csv", "-o", self.output_csv])
            self.assertEqual(ret, 1)


if __name__ == "__main__":
    unittest.main()

