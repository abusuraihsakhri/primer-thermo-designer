"use strict";
let enginePromise = null;
let lastPairs = [];
const $ = id => document.getElementById(id);
const inputNumber = id => Number($(id).value);
const status = (id, message, error = false) => {
  const node = $(id);
  node.textContent = message;
  node.classList.toggle("error", error);
};
const fmt = (n, digits = 1) => Number(n).toFixed(digits);
function parseSequence(raw) {
  if (typeof raw !== "string") throw new Error("Sequence must be text.");
  const lines = raw.split(/\r?\n/).map(x => x.trim()).filter(Boolean);
  const headers = lines.filter(x => x.startsWith(">"));
  if (headers.length > 1) throw new Error("Multiple FASTA records are not supported. Select one sequence.");
  if (headers.length && !lines[0].startsWith(">")) throw new Error("FASTA header must precede sequence lines.");
  const sequence = lines.filter(x => !x.startsWith(">")).join("").toUpperCase();
  if (!sequence) throw new Error("Enter a DNA sequence.");
  if (/[^ACGT]/.test(sequence)) throw new Error("Only unambiguous A, C, G, and T bases are accepted.");
  if (sequence.length > 2000) throw new Error("Browser limit: 2,000 nt. Use the Python CLI for longer targets.");
  return sequence;
}
async function engine() {
  if (!enginePromise) {
    enginePromise = (async () => {
      if (typeof loadPyodide !== "function") throw new Error("Python runtime unavailable. Check your internet connection.");
      const runtime = await loadPyodide({indexURL:"https://cdn.jsdelivr.net/pyodide/v0.27.7/full/"});
      const source = await fetch("./primer_designer.py");
      if (!source.ok) throw new Error("Python source file not found in this deployment.");
      runtime.FS.writeFile("primer_designer.py", await source.text());
      runtime.runPython("from primer_designer import design_primer_pairs, calc_tm, gc_fraction, has_gc_clamp, longest_homopolymer_run, three_prime_stability, sliding_window_complementarity, find_hairpin, validate_sequence; from dataclasses import asdict; import json");
      return runtime;
    })().catch(err => { enginePromise = null; throw err; });
  }
  return enginePromise;
}
function execute(runtime, payload, code) {
  runtime.globals.set("web_payload_json", JSON.stringify(payload));
  try { return JSON.parse(runtime.runPython(code)); }
  finally { runtime.globals.delete("web_payload_json"); }
}
function addDatum(container, title, value, mono = false) {
  const cell = document.createElement("div");
  cell.className = "datum";
  const label = document.createElement("small");
  label.textContent = title;
  const data = document.createElement("span");
  data.className = "value" + (mono ? " mono" : "");
  data.textContent = String(value);
  cell.append(label, data);
  container.append(cell);
}
function appendResult(container, title, datas, warnings = []) {
  const card = document.createElement("article");
  card.className = "result-card";
  const heading = document.createElement("div");
  heading.className = "result-title";
  heading.textContent = title;
  const grid = document.createElement("div");
  grid.className = "data-grid";
  for (const [label, value, mono] of datas) addDatum(grid, label, value, mono);
  card.append(heading, grid);
  if (warnings.length) {
    const risk = document.createElement("p");
    risk.className = "risk";
    risk.textContent = "Heuristic warnings: " + warnings.join("; ");
    card.append(risk);
  }
  container.append(card);
}
function csvCell(value) {
  const text = String(value ?? "");
  return '"' + text.replace(/"/g, '""') + '"';
}
function exportCSV() {
  if (!lastPairs.length) return;
  const header = ["rank","forward","reverse","forward_start_0based","reverse_end_exclusive","product_bp","forward_tm_C","reverse_tm_C","forward_gc_percent","reverse_gc_percent","pair_score","warnings"];
  const lines = [header.map(csvCell).join(",")];
  lastPairs.forEach((pair, idx) => {
    const warnings = [...pair.forward.warnings,...pair.reverse.warnings,...pair.warnings];
    const row = [idx + 1,pair.forward.seq,pair.reverse.seq,pair.forward.start,pair.reverse.end,pair.product_size,pair.forward.tm,pair.reverse.tm,pair.forward.gc_percent,pair.reverse.gc_percent,pair.score,warnings.join("; ")];
    lines.push(row.map(csvCell).join(","));
  });
  const url = URL.createObjectURL(new Blob([lines.join("\r\n") + "\r\n"], {type:"text/csv;charset=utf-8"}));
  const link = document.createElement("a");
  link.href = url; link.download = "primer-pairs.csv"; link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
for (const button of document.querySelectorAll("[data-view]")) {
  button.addEventListener("click", () => {
    const chosen = button.dataset.view;
    $("design-view").classList.toggle("hidden", chosen !== "design");
    $("check-view").classList.toggle("hidden", chosen !== "check");
    for (const tab of document.querySelectorAll("[data-view]")) {
      tab.classList.toggle("active", tab === button);
      tab.setAttribute("aria-pressed", String(tab === button));
    }
  });
}
$("target").addEventListener("input", () => {
  $("sequence-length").textContent = $("target").value.replace(/\s/g, "").replace(/^>[^\n]*\n?/gm, "").length + " nt";
});
$("fasta-file").addEventListener("change", async event => {
  const file = event.target.files?.[0];
  if (!file) return;
  if (file.size > 200000) { status("design-status","File is too large (200 KB limit).",true); return; }
  $("target").value = await file.text();
  $("target").dispatchEvent(new Event("input"));
});
$("csv-button").addEventListener("click", exportCSV);
$("design-form").addEventListener("submit", async event => {
  event.preventDefault();
  lastPairs = []; $("csv-button").disabled = true; $("design-results").replaceChildren();
  const button = $("design-button"); button.disabled = true;
  try {
    const target = parseSequence($("target").value);
    const payload = {
      target,
      product_size_range:[inputNumber("product-min"),inputNumber("product-max")],
      primer_len_range:[inputNumber("length-min"),inputNumber("length-max")],
      tm_range:[inputNumber("tm-min"),inputNumber("tm-max")],
      gc_range:[inputNumber("gc-min")/100,inputNumber("gc-max")/100],
      target_tm:inputNumber("target-tm"),
      primer_conc_nm:inputNumber("concentration"),
      na_conc_mm:inputNumber("salt"),
      require_gc_clamp:$("gc-clamp").checked,
      top_n:inputNumber("top-n")
    };
    if (!Number.isInteger(payload.primer_len_range[0]) || !Number.isInteger(payload.primer_len_range[1])) throw new Error("Primer lengths must be integers.");
    if (payload.primer_len_range[1] > 35) throw new Error("Maximum browser primer length is 35 nt.");
    status("design-status","Loading Python and evaluating candidates…");
    const runtime = await engine();
    lastPairs = execute(runtime,payload,"json.dumps([asdict(p) for p in design_primer_pairs(**json.loads(web_payload_json))])");
    status("design-status",lastPairs.length ? lastPairs.length + " ranked pair(s). Coordinates use zero-based start and exclusive end." : "No primer pairs match these constraints. Consider widening the Tm/GC ranges.");
    $("csv-button").disabled = !lastPairs.length;
    lastPairs.forEach((pair,index) => appendResult($("design-results"),
      "Pair " + (index+1) + " · " + pair.product_size + " bp · score " + fmt(pair.score,2),
      [["Forward (5′ → 3′)",pair.forward.seq,true],["Reverse (5′ → 3′)",pair.reverse.seq,true],
       ["Forward Tm",fmt(pair.forward.tm)+" °C"],["Reverse Tm",fmt(pair.reverse.tm)+" °C"],
       ["Forward GC",fmt(pair.forward.gc_percent)+" %"],["Reverse GC",fmt(pair.reverse.gc_percent)+" %"],
       ["Forward start",pair.forward.start],["Reverse end (exclusive)",pair.reverse.end],
       ["Heterodimer longest run",pair.heterodimer.max_run + " bp"]],
      [...pair.forward.warnings,...pair.reverse.warnings,...pair.warnings]));
  } catch(error) { status("design-status",error.message || String(error),true); }
  finally { button.disabled = false; }
});
$("check-form").addEventListener("submit", async event => {
  event.preventDefault(); $("check-results").replaceChildren();
  const button = $("check-button"); button.disabled = true;
  try {
    const seq = parseSequence($("primer-input").value);
    if (seq.length < 2) throw new Error("Primer must be at least 2 nt.");
    status("check-status","Loading Python and assessing sequence…");
    const runtime = await engine();
    const data = execute(runtime,{seq,primer_conc_m:inputNumber("check-concentration")*1e-9,na_conc_m:inputNumber("check-salt")*1e-3},
      "json.dumps((lambda p: {'tm':calc_tm(p['seq'],p['primer_conc_m'],p['na_conc_m']), 'gc':gc_fraction(p['seq'])*100, 'clamp':has_gc_clamp(p['seq']), 'poly':longest_homopolymer_run(p['seq']), 'dg3':three_prime_stability(p['seq']), 'self':asdict(sliding_window_complementarity(p['seq'],p['seq'])), 'hairpin':asdict(find_hairpin(p['seq']))})(json.loads(web_payload_json)))");
    status("check-status","Thermodynamic and structure heuristic assessment completed.");
    appendResult($("check-results"),"Primer properties",[["Sequence",seq,true],["Length",seq.length+" nt"],
      ["Tm (salt-adjusted)",fmt(data.tm,2)+" °C"],["GC fraction",fmt(data.gc)+" %"],
      ["3′ GC clamp",data.clamp?"Present":"Absent"],["Longest homopolymer",data.poly+" nt"],
      ["3′ end ΔG at 37 °C",fmt(data.dg3,2)+" kcal/mol"],
      ["Self-dimer max run",data.self.max_run+" bp"],["Hairpin stem",data.hairpin.stem_len+" bp"]],
      [...(data.self.max_run>=5?["Potential self-dimer"]:[]),...(data.hairpin.stem_len>=4?["Potential hairpin"]:[])]);
  } catch(error) { status("check-status",error.message || String(error),true); }
  finally { button.disabled = false; }
});
