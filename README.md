# Primer Thermo Designer

PCR primer design and assessment using a Python implementation of SantaLucia (1998) nearest-neighbor DNA duplex thermodynamics.

## What it does

- Finds and ranks forward/reverse primer pairs by predicted melting temperature, amplicon size and sequence-complementarity heuristics
- Calculates salt-adjusted Tm, GC content, 3′ GC clamp, homopolymer runs and approximate 3′ stability
- Screens self-dimers, heterodimers and hairpins using contiguous complementarity heuristics
- Accepts a raw ACGT target or **one** FASTA record, with zero-based, half-open target coordinates
- Exports candidate pairs to CSV in the browser

**Scope:** This is a candidate-screening tool, not a substitute for Primer3, genome-wide off-target searches, Mg²⁺/dNTP-aware models, or wet-lab verification. The supplied Tm model approximates monovalent salt effects and perfect-match duplexes; thermodynamic scores and structural flags are heuristic.

## Browser application

The interface in `web/` runs the **same** `primer_designer.py` engine inside the browser using Pyodide (Python/WebAssembly). No application server or API credentials are needed. Sequence processing stays inside the browser; the first launch downloads the Pyodide runtime from jsDelivr. No sequence is submitted to that CDN, but loading runtime assets does make a network request.

To preview locally after cloning:

```bash
mkdir -p _site
cp web/index.html web/app.js web/style.css primer_designer.py _site/
python -m http.server --directory _site 8000
# Open http://localhost:8000
```

GitHub Actions includes a Pages deployment workflow for `master`. A public deployment link will be added after the deployed site has been verified.

## Python CLI

Requires Python 3.10 or later. The thermodynamics engine itself uses only the standard library.

```bash
git clone https://github.com/abusuraihsakhri/primer-thermo-designer.git
cd primer-thermo-designer
python -m pip install -e ".[dev]"

python primer_designer.py check ACGTGCATCGTAGCTAGCTG
python primer_designer.py design --fasta target.fasta --product-min 100 --product-max 300
python -m pytest -q
```

The installed `primer-design` console entry point provides the `design`/`check` commands. The existing `primer-thermo-designer` entry point continues to operate the auxiliary audit CLI.

Main defaults: 18–25 nt primer length, 57–63 °C Tm, 30–70% GC, 100–300 bp amplicon, 250 nM primer concentration and 50 mM monovalent salt. Use `python primer_designer.py design --help` for all options.

## Optional API and audit service

```bash
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_hex(32))')"
uvicorn agents.api:app --host 127.0.0.1 --port 8000
```

- `POST /api/primer/design`: rank primer pairs from JSON design parameters
- `POST /api/primer/check`: assess one primer sequence
- `GET /health`: service health
- Legacy `/api/audit`, `/api/chat`, `/api/audit/logs`, and `/metrics` remain available for compatibility; they are **not** PCR thermodynamics calculations.
- `/metrics` currently returns JSON, not Prometheus exposition text.
- The audit ledger is **in-memory** and HMAC-chained; it is not persistent or a regulatory-grade audit system. Its pattern-based PHI filter does not establish clinical privacy compliance.

Container operation is available through `Dockerfile` and `docker-compose.yml`; set `AUDIT_SECRET_KEY` using the provided `.env.example` before launching.

## Scientific reference

SantaLucia J Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proc Natl Acad Sci USA*. 1998;95(4):1460–1465. https://doi.org/10.1073/pnas.95.4.1460

## Testing and implementation

- Scientific engine: Python standard library
- Optional API: FastAPI, Pydantic, Uvicorn
- Browser UI: HTML, CSS, JavaScript, Pyodide 0.27.7
- CI: Python 3.10, 3.11 and 3.12; Pytest and JavaScript syntax checks
- Supported browsers: modern evergreen browsers with WebAssembly and network access to jsDelivr; other environments are not validated.

## License

[MIT](LICENSE)
