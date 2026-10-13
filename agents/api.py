"""
FastAPI REST API Server for Primer Thermo Designer.
"""
from dataclasses import asdict
from typing import Dict, Any, List, Tuple
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
from primer_designer import (
    calc_tm, gc_fraction, has_gc_clamp, longest_homopolymer_run,
    three_prime_stability, sliding_window_complementarity, find_hairpin,
    validate_sequence, parse_fasta, design_primer_pairs,
)
from .base import AuditLogger, PHIGuard
from .models import SystemTaskPayload, ConsensusDossier
from .supervisor import SystemSupervisor

supervisor = SystemSupervisor(model_provider="mock")

app = FastAPI(
    title="Primer Thermo Designer API",
    description="PCR primer design and sequence assessment, plus legacy audit endpoints",
    version="1.1.0",
)


class ChatRequest(BaseModel):
    query: str


@app.get("/health")
def health():
    return {"status": "HEALTHY", "service": "primer-thermo-designer", "domain": "PCR primer thermodynamics", "version": "1.1.0"}


@app.get("/metrics")
def metrics():
    return {
        "dossiers_processed_total": len(supervisor.dossier_registry),
        "audit_blocks_total": len(AuditLogger.get_trail()),
        "system_status": "NOMINAL_OPTIMAL"
    }


@app.post("/api/audit")
def api_audit(payload: SystemTaskPayload):
    dossier = supervisor.process_task(payload)
    return dossier.to_dict()


@app.post("/api/chat")
def api_chat(req: ChatRequest):
    try:
        ans = supervisor.query_supervisory_chat(req.query)
        return {"response": ans}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/audit/logs")
def api_audit_logs():
    return {"audit_trail": AuditLogger.get_trail(), "verified": AuditLogger.verify_integrity()}


class PrimerCheckRequest(BaseModel):
    sequence: str = Field(..., min_length=2, max_length=100)
    primer_conc_nm: float = Field(default=250.0, gt=0)
    na_conc_mm: float = Field(default=50.0, gt=0)


class PrimerDesignRequest(BaseModel):
    target: str = Field(..., min_length=4, max_length=5000)
    product_size_range: Tuple[int, int] = (100, 300)
    primer_len_range: Tuple[int, int] = (18, 25)
    tm_range: Tuple[float, float] = (57.0, 63.0)
    target_tm: float = 60.0
    primer_conc_nm: float = Field(default=250.0, gt=0)
    na_conc_mm: float = Field(default=50.0, gt=0)
    gc_range: Tuple[float, float] = (0.30, 0.70)
    max_poly_run: int = Field(default=4, ge=1)
    require_gc_clamp: bool = True
    top_n: int = Field(default=5, ge=1, le=50)


@app.post("/api/primer/check")
def check_primer(req: PrimerCheckRequest):
    try:
        seq = validate_sequence(req.sequence)
        return {
            "sequence": seq,
            "tm_c": calc_tm(seq, req.primer_conc_nm * 1e-9, req.na_conc_mm * 1e-3),
            "gc_percent": gc_fraction(seq) * 100,
            "gc_clamp": has_gc_clamp(seq),
            "longest_homopolymer": longest_homopolymer_run(seq),
            "three_prime_dg37_kcal_mol": three_prime_stability(seq),
            "self_dimer": asdict(sliding_window_complementarity(seq, seq)),
            "hairpin": asdict(find_hairpin(seq)),
        }
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.post("/api/primer/design")
def design_primers(req: PrimerDesignRequest):
    try:
        data = req.model_dump()
        data["target"] = validate_sequence(parse_fasta(data["target"]))
        return {"pairs": [asdict(pair) for pair in design_primer_pairs(**data)]}
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
