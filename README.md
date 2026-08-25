# Bacterial Genome Synteny Agent & Comparative Genomics Toolkit

A specialized comparative microbial genomics engine for analyzing chromosomal synteny blocks, large-scale inversions, Average Nucleotide Identity (ANI), gene order collinearity, breakpoint junction validation primers, and pangenome core/shell/cloud partitioning.

## Core Capabilities

- **Collinear & Inverted Synteny Block Mapping**:
  - Detects forward collinear alignments and reverse-strand inversions across bacterial chromosomes.
  - Computes genome-wide syntenic coverage and conservation statistics.
- **Microbial Species Delimitation & ANI Calculation**:
  - Calculates weighted Average Nucleotide Identity (ANI) over all aligned blocks.
  - Implements standard taxonomic species thresholds (Konstantinidis & Tiedje 2005, Richter & Rosselló-Móra 2009):
    - $\ge 95.0\%$: Same Species (Conspecific boundary)
    - $93.0\% - 94.9\%$: Subspecies / Borderline divergence
    - $85.0\% - 92.9\%$: Different Species (Con-generic)
    - $< 85.0\%$: Divergent Taxon / Different Genus
- **Gene Order Collinearity Metrics**:
  - **Synteny Conservation Score (SCS)**: Conserved adjacent gene pair fraction $\frac{\text{Conserved Adjacencies}}{N - 1}$.
  - **Spearman's Rank Correlation ($\rho$)**: Measures global collinearity of orthologous gene positions.
  - **Kendall's Rank Correlation ($\tau$)**: Quantifies concordant vs discordant pairwise gene orientations.
  - **Breakpoint Distance**: Exact count of chromosomal rearrangement boundaries.
- **ASCII Comparative Synteny Dotplot Visualizer**:
  - Renders 2D ASCII matrices highlighting forward collinear synteny (`\`), inversions (`/`), and unaligned regions (`.`).
- **Breakpoint Junction PCR Primer Designer**:
  - Automatically designs left/right flanking validation primers spanning rearrangement junctions.
  - Salt-adjusted nearest-neighbor melting temperature ($T_m$), GC content calculation, and 3' GC clamp checking.
- **Pangenome Partitioning**:
  - Partitions gene families into **Core Genome** ($\ge 99\%$), **Soft Core** ($95\% - 99\%$), **Shell Genome** ($15\% - 95\%$), and **Cloud Genome** ($< 15\%$).
  - Characterizes open vs closed pangenome trajectory.

---

## Benchmark Datasets

Includes curated, gold-standard bacterial comparison pairs:
1. **Escherichia coli K-12 MG1655 vs Escherichia coli O157:H7 EDL933**:
   - Captures large-scale inversion around replication terminus and O-island pathogenic insertions.
2. **Salmonella enterica serovar Typhimurium LT2 vs Typhi CT18**:
   - High ANI (>99%) with characteristic *Salmonella* pathogenicity island rearrangements and pseudogene accumulation.
3. **Pseudomonas aeruginosa PAO1 vs LESB58**:
   - Epidemic cystic fibrosis isolate genomic islands and inversion blocks.

---

## CLI Usage

### 1. Analyze Benchmark Synteny Pair
```bash
python cli.py --benchmark ecoli_k12_vs_o157
```

### 2. Render ASCII Synteny Dotplot
```bash
python cli.py --benchmark ecoli_k12_vs_o157 --dotplot
```

### 3. Pangenome Core/Shell/Cloud Partitioning
```bash
python cli.py --pangenome
```

### 4. Design Breakpoint Junction PCR Primers
```bash
python cli.py --primers
```

### 5. Structured JSON Output
```bash
python cli.py --benchmark salmonella_lt2_vs_ct18 --json
```

### 6. Interactive Comparative Genomics Shell
```bash
python cli.py --interactive
```

---

## Unit Testing

Run unit tests via Python's standard `unittest` framework:

```bash
python -m unittest test_bacterial_genome_synteny.py
```

All 26 test cases validate ANI calculations, rank correlations, dotplot generation, primer quality checks, pangenome partitions, and CLI commands.
