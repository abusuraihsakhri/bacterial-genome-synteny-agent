#!/usr/bin/env python3
"""
Command Line Interface for Bacterial Genome Synteny Agent & Comparative Genomics Toolkit.

Usage examples:
  python cli.py --benchmark ecoli_k12_vs_o157
  python cli.py --benchmark ecoli_k12_vs_o157 --dotplot
  python cli.py --pangenome
  python cli.py --primers
  python cli.py --list-benchmarks
  python cli.py --interactive
  python cli.py --benchmark salmonella_lt2_vs_ct18 --json
"""

import argparse
import json
import sys
from typing import List, Optional

from bacterial_genome_synteny import (
    SyntenyBlock,
    SyntenyAnalyzer,
    DotplotRenderer,
    PrimerDesigner,
    PangenomeAnalyzer,
    get_curated_synteny_benchmarks,
)


def format_synteny_summary(res, show_dotplot: bool = False) -> str:
    lines = []
    lines.append("=" * 80)
    lines.append(f"  BACTERIAL GENOME SYNTENY REPORT")
    lines.append(f"  Reference: {res.reference_name} ({res.ref_genome_length:,} bp)")
    lines.append(f"  Query:     {res.query_name} ({res.qry_genome_length:,} bp)")
    lines.append("-" * 80)
    lines.append(f"  Taxonomic Species Call:  [{res.taxonomic_call}]")
    lines.append(f"  Average Nucleotide Identity (ANI): {res.ani_pct:.2f}%")
    lines.append(f"  Genome Coverage: Ref: {res.ref_coverage_pct:.1f}% | Qry: {res.qry_coverage_pct:.1f}%")
    lines.append(f"  Syntenic Blocks: {res.num_blocks} total ({res.collinear_blocks} collinear, {res.inverted_blocks} inverted)")
    lines.append(f"  Gene Order Metrics:")
    lines.append(f"    * Conserved Adjacency Score: {res.conserved_adjacency_score:.4f}")
    lines.append(f"    * Spearman's Rank Rho:       {res.spearman_rho:.4f}")
    lines.append(f"    * Kendall's Rank Tau:        {res.kendall_tau:.4f}")
    lines.append(f"    * Breakpoint Count:          {res.breakpoint_count}")
    lines.append("-" * 80)
    lines.append("  Syntenic Alignment Blocks:")
    for b in res.blocks:
        strand_sym = "(+ Collinear)" if not b.is_inverted else "(- INVERTED)"
        lines.append(f"    * [{b.block_id:10s}] Ref: {b.ref_start:>9,}-{b.ref_end:>9,} ({b.ref_length:>8,} bp) <==> "
                     f"Qry: {b.qry_start:>9,}-{b.qry_end:>9,} ({b.qry_length:>8,} bp) | ID: {b.identity_pct:.1f}% {strand_sym}")
    if show_dotplot:
        lines.append("-" * 80)
        lines.append("  ASCII Comparative Synteny Dotplot ('\\' Collinear, '/' Inversion):")
        lines.append(DotplotRenderer.render(res.blocks, res.ref_genome_length, res.qry_genome_length, grid_size=30))
    lines.append("=" * 80)
    return "\n".join(lines)


def run_interactive():
    benchmarks = get_curated_synteny_benchmarks()
    print("\n" + "=" * 80)
    print("  BACTERIAL GENOME SYNTENY AGENT - INTERACTIVE SHELL")
    print("=" * 80)
    print("Type 'help' for commands, 'exit' or 'quit' to exit.\n")

    while True:
        try:
            cmd_line = input("synteny-agent> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nExiting.")
            break

        if not cmd_line:
            continue

        parts = cmd_line.split()
        cmd = parts[0].lower()

        if cmd in ("exit", "quit", "q"):
            print("Goodbye.")
            break
        elif cmd == "help":
            print("""
Available Commands:
  benchmarks                      List available comparative genomic benchmark pairs
  analyze <benchmark_id>          Run complete synteny & ANI analysis
  dotplot <benchmark_id>          Render 2D ASCII dotplot of collinear/inverted blocks
  pangenome                       Run pangenome core/shell/cloud partitioning
  primers                         Design breakpoint PCR validation primers
  help                            Show this menu
  exit / quit                     Exit the shell
            """)
        elif cmd == "benchmarks":
            print("\nAvailable Benchmark Genome Pairs:")
            for k, v in benchmarks.items():
                print(f"  [{k:25s}] {v['ref_name']} vs {v['qry_name']}")
            print()
        elif cmd == "analyze":
            if len(parts) < 2:
                print("Usage: analyze <benchmark_id>")
                continue
            b_id = parts[1]
            if b_id not in benchmarks:
                print(f"Benchmark '{b_id}' not found. Use 'benchmarks' to see options.")
                continue
            b_data = benchmarks[b_id]
            res = SyntenyAnalyzer.analyze(
                b_data["ref_name"], b_data["qry_name"],
                b_data["ref_length"], b_data["qry_length"],
                b_data["blocks"]
            )
            print(format_synteny_summary(res, show_dotplot=False))
        elif cmd == "dotplot":
            if len(parts) < 2:
                print("Usage: dotplot <benchmark_id>")
                continue
            b_id = parts[1]
            if b_id not in benchmarks:
                print(f"Benchmark '{b_id}' not found.")
                continue
            b_data = benchmarks[b_id]
            res = SyntenyAnalyzer.analyze(
                b_data["ref_name"], b_data["qry_name"],
                b_data["ref_length"], b_data["qry_length"],
                b_data["blocks"]
            )
            print(format_synteny_summary(res, show_dotplot=True))
        elif cmd == "pangenome":
            sample_matrix = {
                "dnaA": {"S1", "S2", "S3", "S4", "S5"},
                "rpoB": {"S1", "S2", "S3", "S4", "S5"},
                "gyrA": {"S1", "S2", "S3", "S4", "S5"},
                "recA": {"S1", "S2", "S3", "S4", "S5"},
                "lacZ": {"S1", "S2", "S3"},
                "cysK": {"S1", "S2", "S4"},
                "blaCTX-M-15": {"S2", "S3"},
                "colV_plasmid_tra": {"S4"},
                "phage_tail_integrase": {"S5"},
            }
            res = PangenomeAnalyzer.partition(sample_matrix)
            print("\nPangenome Partitioning (5 Strains, 9 Gene Families):")
            print(f"  Core ({res['core_count']}):       {', '.join(res['core_families'])}")
            print(f"  Shell ({res['shell_count']}):      {', '.join(res['shell_families'])}")
            print(f"  Cloud ({res['cloud_count']}):      {', '.join(res['cloud_families'])}")
            print(f"  Pangenome Trajectory: [{res['pangenome_type']}]\n")
        elif cmd == "primers":
            lf = "ACGTACGTACGTACGTACGTGGCACTGATTCAGGCATTAAGCTT"
            rf = "TTGACCATTGCCATTCGAAACGTACGTACGTACGTACGT"
            primers = PrimerDesigner.design_flanking_primers(lf, rf)
            print("\nBreakpoint Junction Validation Primers:")
            print(f"  Forward: {primers['forward_primer']['sequence']} (Tm: {primers['forward_primer']['tm_c']}C, GC: {primers['forward_primer']['gc_pct']}%)")
            print(f"  Reverse: {primers['reverse_primer']['sequence']} (Tm: {primers['reverse_primer']['tm_c']}C, GC: {primers['reverse_primer']['gc_pct']}%)")
            print(f"  QC Passed: {primers['overall_qc_passed']}\n")
        else:
            print(f"Unknown command '{cmd}'. Type 'help' for instructions.")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="bacterial-genome-synteny-agent",
        description="Bacterial Genome Synteny Agent & Comparative Microbial Genomics Toolkit",
    )

    parser.add_argument("--benchmark", "-b", help="Analyze curated benchmark genome pair (e.g., ecoli_k12_vs_o157)")
    parser.add_argument("--dotplot", action="store_true", help="Render ASCII comparative synteny dotplot")
    parser.add_argument("--pangenome", action="store_true", help="Run pangenome core/shell/cloud partitioning")
    parser.add_argument("--primers", action="store_true", help="Design PCR primers for breakpoint junction validation")
    parser.add_argument("--list-benchmarks", action="store_true", help="List all available curated benchmark pairs")
    parser.add_argument("--interactive", "-i", action="store_true", help="Launch interactive comparative genomics shell")
    parser.add_argument("--json", action="store_true", help="Output results in JSON format")

    args = parser.parse_args(argv)
    benchmarks = get_curated_synteny_benchmarks()

    if args.interactive:
        run_interactive()
        return 0

    if args.list_benchmarks:
        if args.json:
            print(json.dumps({k: {"ref": v["ref_name"], "qry": v["qry_name"], "blocks": len(v["blocks"])}
                             for k, v in benchmarks.items()}, indent=2))
        else:
            print("\nAvailable Curated Bacterial Synteny Benchmarks:")
            print("=" * 70)
            for k, v in benchmarks.items():
                print(f"  [{k:26s}] {v['ref_name']} <==> {v['qry_name']}")
            print("=" * 70)
        return 0

    if args.benchmark:
        if args.benchmark not in benchmarks:
            err = {"error": f"Benchmark '{args.benchmark}' not found. Choose from: {list(benchmarks.keys())}"}
            if args.json:
                print(json.dumps(err, indent=2))
            else:
                print(f"Error: {err['error']}", file=sys.stderr)
            return 1

        b_data = benchmarks[args.benchmark]
        res = SyntenyAnalyzer.analyze(
            b_data["ref_name"], b_data["qry_name"],
            b_data["ref_length"], b_data["qry_length"],
            b_data["blocks"]
        )

        if args.json:
            out_dict = res.to_dict()
            if args.dotplot:
                out_dict["dotplot_ascii"] = DotplotRenderer.render(res.blocks, res.ref_genome_length, res.qry_genome_length)
            print(json.dumps(out_dict, indent=2))
        else:
            print(format_synteny_summary(res, show_dotplot=args.dotplot))
        return 0

    if args.pangenome:
        sample_matrix = {
            "dnaA": {"Strain_A", "Strain_B", "Strain_C", "Strain_D", "Strain_E"},
            "rpoB": {"Strain_A", "Strain_B", "Strain_C", "Strain_D", "Strain_E"},
            "gyrA": {"Strain_A", "Strain_B", "Strain_C", "Strain_D", "Strain_E"},
            "lacZ": {"Strain_A", "Strain_B", "Strain_C"},
            "cysK": {"Strain_A", "Strain_B", "Strain_D"},
            "blaKPC_carbapenemase": {"Strain_B", "Strain_E"},
            "phage_tail_integrase": {"Strain_E"},
            "capsule_k1_polysaccharide": {"Strain_A"},
        }
        res = PangenomeAnalyzer.partition(sample_matrix)
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            print("\nPangenome Partitioning Analysis:")
            print("=" * 70)
            print(f"  Strains Evaluated:   {res['n_strains']}")
            print(f"  Total Gene Families: {res['total_families']}")
            print(f"  Core Genome (>=99%): {res['core_count']} ({res['core_pct']}%) -> {', '.join(res['core_families'])}")
            print(f"  Shell Genome (15-95%): {res['shell_count']} ({res['shell_pct']}%) -> {', '.join(res['shell_families'])}")
            print(f"  Cloud Genome (<15%): {res['cloud_count']} ({res['cloud_pct']}%) -> {', '.join(res['cloud_families'])}")
            print(f"  Trajectory:          {res['pangenome_type']}")
            print("=" * 70)
        return 0

    if args.primers:
        left_flank = "GATCGATCGATCGATCGATC" + "GGCACTGATTCAGGCATTAAGCTT"
        right_flank = "TTGACCATTGCCATTCGAA" + "GATCGATCGATCGATCGATC"
        p_res = PrimerDesigner.design_flanking_primers(left_flank, right_flank)
        if args.json:
            print(json.dumps(p_res, indent=2))
        else:
            print("\nInversion Breakpoint Validation PCR Primers:")
            print("=" * 70)
            fwd = p_res["forward_primer"]
            rev = p_res["reverse_primer"]
            print(f"  Forward Primer: 5'- {fwd['sequence']} -3'")
            print(f"    Length: {fwd['length']} bp | Tm: {fwd['tm_c']} C | GC: {fwd['gc_pct']}% | GC Clamp: {fwd['gc_clamp']}")
            print(f"  Reverse Primer: 5'- {rev['sequence']} -3'")
            print(f"    Length: {rev['length']} bp | Tm: {rev['tm_c']} C | GC: {rev['gc_pct']}% | GC Clamp: {rev['gc_clamp']}")
            print(f"  Tm Delta: {p_res['tm_difference']} C | QC Status: {'PASSED' if p_res['overall_qc_passed'] else 'FLAGGED'}")
            print(f"  Amplicon Size: ~{p_res['amplicon_size_bp']} bp")
            print("=" * 70)
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
