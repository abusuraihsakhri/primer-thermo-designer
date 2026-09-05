# Primer Thermo Designer

> **Domain:** PCR Primer Design & Thermodynamic Validation

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Primer Thermo Designer designs and validates PCR primer pairs from a target DNA sequence using
the SantaLucia (1998) unified nearest-neighbor (NN) thermodynamic
parameters for melting-temperature (Tm) prediction, combined with
sliding-window complementarity scoring for self-dimer, heterodimer and
hairpin screening (in the spirit of Primer3's core heuristics).

Reference:
    SantaLucia J Jr. "A unified view of polymer, dumbbell, and
    oligonucleotide DNA nearest-neighbor thermodynamics."
    Proc Natl Acad Sci USA. 1998;95(4):1460-1465.

---

## ⚙️ Key Capabilities

### Core Primer Design Engine (`primer_designer.py`)

- **SantaLucia Tm Calculation** — Nearest-neighbor thermodynamic model with salt correction
- **Primer Pair Generation** — Enumerates candidate pairs with configurable size/length/Tm ranges
- **Dimer Screening** — Sliding-window complementarity scoring for self-dimer and heterodimer detection
- **Hairpin Detection** — Local self-complementarity scan for intramolecular hairpin structures
- **GC Clamp Validation** — Ensures appropriate 3' end stability
- **Homopolymer Run Detection** — Screens for problematic repetitive sequences

### Enterprise Security & Audit (`agents/`)

- **Zero-PHI Outbound Interceptor** — Regex-based detection blocking SSNs, MRNs, phone numbers, and patient identifiers
- **Tamper-Evident HMAC-SHA256 Audit Trail** — Chained, cryptographically signed logs for every evaluation
- **FastAPI REST API** — OpenAPI 3.1 REST endpoints with Prometheus metrics (`/metrics`)

---

## 💻 Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/primer-thermo-designer.git
cd primer-thermo-designer

# Install dependencies
pip install fastapi uvicorn pydantic

# For development (includes pytest)
pip install -e ".[dev]"
```

---

## 🧪 CLI Usage

### 1. Design Primer Pairs

```bash
# Design primers for a target sequence
python primer_designer.py design --sequence ACGTGCATCGTAGCTAGCTGATCGATCGATCG --product-min 100 --product-max 300

# Using a FASTA file
python primer_designer.py design --fasta target.fasta --product-min 100 --product-max 300
```

### 2. Check a Single Primer

```bash
python primer_designer.py check ACGTGCATCGTAGCTAGCTG
```

### 3. Parameters

| Parameter | Description | Default |
|:----------|:------------|:--------|
| `--sequence` | Target DNA sequence (ACGT) | — |
| `--fasta` | Path to FASTA file | — |
| `--product-min` | Minimum amplicon size (bp) | 100 |
| `--product-max` | Maximum amplicon size (bp) | 300 |
| `--primer-len-min` | Minimum primer length (nt) | 18 |
| `--primer-len-max` | Maximum primer length (nt) | 25 |
| `--tm-min` | Minimum acceptable Tm (°C) | 57.0 |
| `--tm-max` | Maximum acceptable Tm (°C) | 63.0 |
| `--target-tm` | Ideal Tm for ranking (°C) | 60.0 |
| `--primer-conc` | Primer strand concentration (nM) | 250 |
| `--na-conc` | Monovalent salt concentration (mM) | 50 |
| `--gc-min` | Minimum primer GC% | 30 |
| `--gc-max` | Maximum primer GC% | 70 |
| `--max-poly-run` | Maximum homopolymer run length | 4 |
| `--no-gc-clamp` | Disable GC clamp requirement | False |
| `--top` | Number of top pairs to report | 5 |

---

## 🚀 REST API Server

### Start the Server

```bash
# Using the CLI
python cli.py serve --host 127.0.0.1 --port 8000

# Or directly
uvicorn agents.api:app --host 127.0.0.1 --port 8000
```

### Endpoints

| Method | Path | Description |
|:-------|:-----|:------------|
| GET | `/health` | Health check |
| GET | `/metrics` | Prometheus operational metrics |
| POST | `/api/audit` | Submit task for evaluation |
| POST | `/api/chat` | Supervisory chat query |
| GET | `/api/audit/logs` | Retrieve HMAC audit trail |

---

## 🐳 Docker Deployment

```bash
# Build and run with Docker Compose
cp .env.example .env
# Edit .env to set AUDIT_SECRET_KEY
docker-compose up --build

# Or use Docker directly
docker build -t primer-thermo-designer .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-secret-key primer-thermo-designer
```

---

## 🧪 Testing

```bash
# Run all tests
pytest -v

# Run specific test files
pytest test_primer_designer.py -v
pytest tests/ -v
```

### Test Coverage

- **test_primer_designer.py** — 19 tests for core primer design functionality
- **tests/test_primer_thermo_designer.py** — 5 tests for enterprise security features
- **tests/test_enrichment.py** — 2 tests for enrichment modules

**Total: 24 tests**

---

## 📁 Project Structure

```
primer-thermo-designer/
├── primer_designer.py      # Core primer design engine (SantaLucia NN)
├── cli.py                  # Main CLI interface
├── enrichment.py           # Enrichment feature modules
├── simulator.py            # High-throughput simulation/testing
├── agents/
│   ├── __init__.py
│   ├── api.py              # FastAPI REST server
│   ├── base.py             # PHI guard, HMAC audit trail
│   ├── models.py           # Pydantic data models
│   ├── supervisor.py       # Multi-worker orchestrator
│   ├── workers.py          # Specialized evaluation workers
│   ├── llm_factory.py      # LLM provider abstraction
│   ├── learning.py         # Bayesian calibration engine
│   ├── metrics.py          # Prometheus metrics exporter
│   └── streamer.py         # WebSocket telemetry
├── tests/                  # Pytest test suite
├── web/index.html          # Operations console UI
├── Dockerfile
├── docker-compose.yml
├── pyproject.toml
└── README.md
```

---

## 🔒 Security

- **No hardcoded secrets** — `AUDIT_SECRET_KEY` must be set via environment variable
- **PHI outbound guard** — Prevents accidental leakage of patient identifiers
- **HMAC-SHA256 audit trail** — Tamper-evident logging with chained hashes

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
