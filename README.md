# Bacterial Genome Synteny Agent

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Bacterial Genome Synteny Agent & Comparative Microbial Genomics Toolkit.

Implements:
- Syntenic Collinear Block Detection & Inversion Identification
- Synteny Conservation Metrics: Conserved Adjacency Fraction, Breakpoint Distance
- Rank Correlation Metrics: Spearman's Rho & Kendall's Tau on Ortholog Order
- Average Nucleotide Identity (ANI) & Microbial Species Delimitation Cutoffs
- ASCII Dotplot Matrix Rendering (Collinear '', Inversion '/', Unaligned '.')
- Inversion Breakpoint Junction Detection & PCR Primer Design (Tm, GC clamp)
- Pangenome Partitioning (Core >=99%, Soft-Core >=95%, Shell 15-95%, Cloud <15%)

Pure Python Standard Library (no external dependencies required).

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`SyntenyBlock`**: Collinear or inverted syntenic block between reference and query genomes.
- **`SyntenyAnalysisResult`**: Overall comparative synteny and ANI report.
- **`SyntenyAnalyzer`**: Core comparative genomics synteny engine.
- **`DotplotRenderer`**: ASCII 2D Comparative Synteny Dotplot Visualizer.
- **`PrimerDesigner`**: Breakpoint Junction PCR Primer Design for Inversion Validation.
- **`PangenomeAnalyzer`**: Pangenome Partitioning (Core, Soft Core, Shell, Cloud) across bacterial isolates.

---

## 📐 Mathematical Formulation & Logic

```text
  return (self.qry_end - self.qry_start) * (self.ref_end - self.ref_start) < 0
  return (concordant - discordant) / total_pairs
  Calculates fraction of adjacent gene pairs (A, B) in reference that remain
  score = conserved_count / ref_total_pairs if ref_total_pairs > 0 else 1.0
  fwd_tm = cls.calculate_tm(fwd_seq)
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --benchmark <value> --dotplot <value> --pangenome <value> --primers <value>
```

### Parameter Reference
- `--benchmark`: Specifies input measurement or parameter value.
- `--dotplot`: Specifies input measurement or parameter value.
- `--pangenome`: Specifies input measurement or parameter value.
- `--primers`: Specifies input measurement or parameter value.
- `--list-benchmarks`: Specifies input measurement or parameter value.
- `--interactive`: Specifies input measurement or parameter value.
- `--json`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `suite_name` | Parameter / observation metric | Required |
| `system_slug` | Parameter / observation metric | Required |
| `standard_reference` | Parameter / observation metric | Required |
| `test_cases` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t bacterial-genome-synteny-agent .
docker run -p 8000:8000 bacterial-genome-synteny-agent
```
