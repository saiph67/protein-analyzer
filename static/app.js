// Shared helpers for the Protein Analyzer pages (no build step, no external libraries).
const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];

function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

async function api(url, opts = {}) {
  const init = { ...opts };
  if (opts.json !== undefined) {
    init.method = init.method || "POST";
    init.headers = { "Content-Type": "application/json" };
    init.body = JSON.stringify(opts.json);
  }
  delete init.json;
  let r;
  try { r = await fetch(url, init); }
  catch (e) { throw new Error("Could not reach the local server. Is the app still running?"); }
  let data = null;
  try { data = await r.json(); } catch (e) { /* non-JSON */ }
  if (!r.ok) {
    const msg = (data && data.error) || `Request failed (HTTP ${r.status})`;
    const err = new Error(msg); err.status = r.status; err.source = data && data.source; throw err;
  }
  return data;
}

function showMsg(el, text, kind = "error") {
  el.innerHTML = text ? `<div class="msg ${kind}">${esc(text)}</div>` : "";
}

function busy(btn, on, label) {
  if (on) { btn.dataset.label = btn.textContent; btn.disabled = true; btn.innerHTML = `<span class="spinner"></span> ${esc(label || "Working")}`; }
  else { btn.disabled = false; btn.textContent = btn.dataset.label || btn.textContent; }
}

function download(filename, text, type = "text/plain") {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = filename; document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function toCsv(header, rows) {
  const q = v => { const s = String(v ?? ""); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
  return [header.map(q).join(","), ...rows.map(r => r.map(q).join(","))].join("\n") + "\n";
}

// ---- provenance: version, date, UniProt release, methods, sources on every export --------
const _prov = {};
async function provenanceFor(kind) {
  if (!_prov[kind]) {
    try { _prov[kind] = await api("/api/provenance?kind=" + encodeURIComponent(kind)); }
    catch (e) { _prov[kind] = { app: "Protein Analyzer", generated: new Date().toISOString(), methods: [], sources: [] }; }
  }
  return { ..._prov[kind], generated: new Date().toISOString() };
}
function provenanceLines(p, prefix = "# ") {
  return [`${p.app} · generated ${p.generated}` + (p.uniprot_release ? ` · UniProtKB release ${p.uniprot_release}` : ""),
    ...p.methods.map(m => "Method: " + m), ...p.sources.map(s => "Data: " + s),
    "Methods & citations: see the app's Methods page"].map(l => prefix + l);
}
// Download helpers that prepend/attach provenance.
async function downloadCsv(filename, kind, header, rows) {
  const p = await provenanceFor(kind);
  download(filename, provenanceLines(p).join("\n") + "\n" + toCsv(header, rows), "text/csv");
}
async function downloadJson(filename, kind, data) {
  const p = await provenanceFor(kind);
  download(filename, JSON.stringify({ provenance: p, data }, null, 2), "application/json");
}
async function downloadText(filename, kind, text, prefix = "# ") {
  const p = await provenanceFor(kind);
  download(filename, provenanceLines(p, prefix).join("\n") + "\n\n" + text);
}
// Printed pages get a provenance footer (filled just before printing).
addEventListener("beforeprint", () => {
  const el = document.getElementById("printprov"); if (!el) return;
  const p = _prov[document.body.dataset.kind || ""];
  el.textContent = `Printed ${new Date().toLocaleString()} from ${p ? p.app : "Protein Analyzer"}` +
    (p && p.uniprot_release ? ` · UniProtKB release ${p.uniprot_release}` : "") + ". Methods, data sources and citations: Methods & citations page.";
});

const fmt = (x, d = 2) => (x === null || x === undefined) ? "–" : Number(x).toLocaleString(undefined, { minimumFractionDigits: d, maximumFractionDigits: d });
const pct = (x, d = 1) => (x === null || x === undefined) ? "–" : (100 * x).toFixed(d) + "%";
const ACC_RE = /^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$/i;

// ---- theme ----------------------------------------------------------------
(function initTheme() {
  let t = null; try { t = localStorage.getItem("pa_theme"); } catch (e) {}
  if (t === "light" || t === "dark") document.documentElement.dataset.theme = t;
})();
function setTheme(t) {
  if (t === "auto") delete document.documentElement.dataset.theme; else document.documentElement.dataset.theme = t;
  try { localStorage.setItem("pa_theme", t); } catch (e) {}
  updateThemeButton();
  document.dispatchEvent(new Event("themechange"));
}
function updateThemeButton() {
  const t = document.documentElement.dataset.theme || "auto";
  $$("[data-theme-set]").forEach(b => b.setAttribute("aria-pressed", String(b.dataset.themeSet === t)));
}
matchMedia("(prefers-color-scheme: dark)").addEventListener?.("change", () => document.dispatchEvent(new Event("themechange")));
function isDark() {
  const t = document.documentElement.dataset.theme;
  return t ? t === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
}

// ---- glossary + tooltips ---------------------------------------------------
const GLOSSARY = {
  length: ["Length", "Number of amino acids (residues) in the chain."],
  mw: ["Molecular weight", "Mass of the protein in daltons (Da); 1 kDa = 1000 Da. An average amino acid weighs ~110 Da."],
  gravy: ["GRAVY", "Grand average of hydropathy: the mean Kyte–Doolittle score of all residues. Above 0 = overall water-repelling (often membrane proteins); below 0 = water-loving (often soluble)."],
  pi: ["Isoelectric point (pI)", "The pH at which the protein has no net electric charge. Below its pI a protein is positively charged, above it negatively charged."],
  charge: ["Charge at pH 7", "Estimated net charge at neutral pH, from the acidic (D, E, C, Y, C-terminus) and basic (K, R, H, N-terminus) groups."],
  aromaticity: ["Aromaticity", "Fraction of aromatic residues: phenylalanine (F), tryptophan (W) and tyrosine (Y)."],
  tm: ["Transmembrane (TM) segment", "A stretch that likely crosses a cell membrane. Predicted here as regions where the 19-residue average hydropathy is ≥ 1.6 (Kyte–Doolittle). A rough rule: adjacent helices can merge and signal peptides can look like TM segments."],
  hydropathy: ["Hydropathy plot", "Average water-repelling score along the sequence (sliding 19-residue window). Peaks above the dashed line are candidate membrane-spanning segments. Hover to read values; click a shaded region to highlight it in the sequence."],
  composition: ["Amino-acid composition", "How often each of the 20 standard amino acids occurs, as a percentage of the sequence."],
  identity: ["Identity", "Percentage of alignment columns where both proteins have exactly the same amino acid."],
  similarity: ["Similarity", "Percentage of columns with the same or a chemically similar amino acid (positive BLOSUM62 score). Always ≥ identity."],
  score: ["Alignment score", "Sum of BLOSUM62 substitution scores minus gap penalties. Higher = better alignment; only comparable between alignments of similar lengths."],
  gaps: ["Gaps", "Columns where one protein has an insertion relative to the other (shown as -)."],
  coverage: ["Coverage", "How much of each protein takes part in the alignment. Local alignments often cover only a shared region."],
  global: ["Global vs local", "Global (Needleman–Wunsch) aligns the proteins end to end; use for proteins of similar length and family. Local (Smith–Waterman) finds the best-matching region; use when proteins share only a domain."],
  blosum: ["BLOSUM62", "A table of how often each amino acid swaps for another in related proteins; used to score alignments."],
  kmer: ["3-mer similarity", "Fraction of 3-letter 'words' two sequences share (Jaccard index). A fast, rough similarity: unrelated proteins ≈ 0.02, close relatives > 0.3."],
  reviewed: ["Reviewed (Swiss-Prot)", "Entries checked by UniProt curators. Unreviewed (TrEMBL) entries are automatic annotations: far more numerous but less reliable."],
  go: ["Gene Ontology (GO)", "Standard vocabulary for what a protein does (molecular function), which processes it takes part in (biological process) and where it is in the cell (cellular component)."],
  ec: ["EC number", "Enzyme Commission number classifying the reaction an enzyme catalyses, e.g. 2.7.x.x = transferring phosphate groups (kinases)."],
  keyword: ["UniProt keywords", "Curated labels summarising an entry. The ML classes here are defined by keywords such as KW-0418 (Kinase)."],
  gc: ["GC content", "Share of G and C bases. It varies between species (about 25–75%) and is higher in genes than in the rest of many genomes."],
  orf: ["Open reading frame (ORF)", "A stretch of DNA that starts with a start codon (usually ATG) and runs without a stop codon, so it could encode a protein. Long ORFs are good gene candidates; very short ones often occur by chance."],
  cam: ["Carbamidomethylation", "Before mass spectrometry, cysteines are usually blocked with iodoacetamide, which adds 57.021 Da to each one. Tick this if your sample was treated that way."],
  coverage_ms: ["Sequence coverage", "Share of the protein covered by peptides you could detect (here: those within your length and mass limits). Higher coverage = more confident identification."],
  mono: ["Monoisotopic mass", "Mass using only the most common isotope of each element (12C, 1H, 14N, 16O, 32S): the first peak a high-resolution mass spectrometer sees. [M+H]⁺ etc. are the m/z values for ions carrying 1, 2 or 3 extra protons."],
  epsilon: ["Extinction coefficient (ε)", "How strongly 1 M of the protein absorbs 280 nm light through 1 cm. Mostly from tryptophan (5500 each) and tyrosine (1490 each), plus a little from disulfide bonds (125 each)."],
  ss: ["Secondary structure", "Local shapes of the chain: α-helices (spirals) and β-strands (extended segments that pair into sheets). Here they are read from the AlphaFold model's geometry, so they are a prediction of a prediction, but usually accurate for confident regions."],
  disorder: ["Disordered regions", "Parts with no fixed 3D shape that stay flexible, often for binding partners. AlphaFold confidence below 50 is a good indicator of disorder."],
  string: ["STRING interactions", "Known and predicted partners from STRING. Evidence combines lab experiments, curated databases, text mining, co-expression and genomic context. Score 0–1: above 0.7 = high confidence, 0.4–0.7 = medium. Lines show associations, not necessarily direct physical contact."],
  reactome: ["Reactome pathways", "Expert-curated maps of biological processes (signalling, metabolism, disease) that this protein takes part in. Mainly human, with predicted versions for other model organisms."],
  tmscore: ["TM-score", "Structural similarity from 0 to 1 that doesn't depend much on protein length. Above 0.5: same overall fold. Below 0.3: essentially unrelated shapes."],
  rmsd: ["RMSD", "Root-mean-square distance between matched C-alpha atoms after the best overlay, in ångströms (Å). Under ~2 Å is very close for the core of related proteins."],
  pfam: ["Pfam domains", "Pfam groups protein regions into families that share an evolutionary origin and usually a 3D fold and function, e.g. 'protein kinase domain'. Recognising a known domain is one of the most reliable ways to guess what a new protein does."],
  am: ["AlphaMissense", "A DeepMind model (2023) that predicts, for human proteins, how likely each single amino-acid change is to cause disease: 0 = likely benign, 1 = likely pathogenic. It mainly detects changes that break the protein, so it can miss 'gain-of-function' mutations. Sickle-cell E7V is a famous example: it scores as benign."],
  plddt: ["pLDDT (AlphaFold confidence)", "AlphaFold's per-residue confidence (0–100). Above 90: very high; 70–90: confident; 50–70: low; below 50: very low, often a disordered region rather than a wrong prediction."],
  features: ["Sequence annotations", "Regions and sites recorded by UniProt curators (from experiments or reliable inference). Shown next to this app's own TM prediction so you can compare them."],
  signal: ["Signal peptide", "A short N-terminal tag that sends the protein into the secretory pathway; usually cut off. Hydrophobic, so simple TM predictors often mistake it for a TM segment."],
  accuracy: ["Accuracy", "Fraction of hold-out proteins whose class was predicted correctly."],
  balanced_accuracy: ["Balanced accuracy", "Average of the per-class recall, so every class counts equally even if some have more proteins."],
  macro_f1: ["Macro-F1", "F1 (balance of precision and recall) computed per class, then averaged."],
  holdout: ["Time-split hold-out", "Proteins first published after the cutoff date, never seen in training. Mimics predicting genuinely new proteins."],
  low_sim: ["Low-similarity accuracy", "Accuracy only on hold-out proteins that share little with any training protein (3-mer similarity < 0.10). Near-copies of training proteins are easy, so this is the honest measure."],
  identity_bin: ["< 30% identity", "Accuracy on hold-out proteins less than 30% identical to every training protein: roughly the 'twilight zone' where simple homology lookup stops working."],
  redundancy: ["Redundancy-filtered accuracy", "Accuracy after removing every test protein that is ≥ 40% identical to any of its 5 most similar training proteins (checked by full alignment). This approximates clustering the data at 40% identity, the standard way to stop near-copies leaking between training and test sets."],
  baseline: ["Baselines", "Simple reference methods. A useful model must beat them, especially '1-NN homology transfer' (copy the label of the most similar training protein)."],
  confusion: ["Confusion matrix", "Rows are the true class, columns the predicted class. The diagonal holds correct predictions; off-diagonal cells show which classes get mixed up."],
  chance: ["Chance level", "Accuracy from guessing uniformly at random among the classes."],
  confidence: ["Confidence", "The model's probability for its top class. It is not a guarantee: check the similarity warning too."],
  nn: ["Nearest training protein", "The most similar protein in the training set (3-mer similarity). If it is very low, the model has never seen anything like this protein."],
  other: ["'other' class", "Proteins with none of the chosen function keywords. Lets the model answer 'none of these' instead of forcing every protein into a class."],
  esm: ["ESM-2", "A protein language model from Meta AI, pre-trained on millions of natural sequences to predict masked amino acids. Its internal representation (embedding) captures evolutionary and structural signals. The model never saw our labels; we only train a small classifier on top."],
  handcrafted: ["Hand-crafted features", "434 numbers computed directly from the sequence: amino-acid and dipeptide frequencies, residue-group fractions and physico-chemical properties."],
  cutoff: ["Cutoff date", "Training uses entries first published on or before this date; the hold-out uses entries published after."],
  strata: ["Taxonomic sampling", "Each class is sampled proportionally across bacteria, archaea, viruses, fungi, plants, mammals and other animals/eukaryotes, so the model can't cheat by recognising the organism."],
};

function tip(key) {
  if (!GLOSSARY[key]) return "";
  return `<button type="button" class="tip" data-tip="${key}" aria-label="What is ${esc(GLOSSARY[key][0])}?">i</button>`;
}

(function initTooltips() {
  let box = null, pinned = null;
  const show = (el) => {
    const g = GLOSSARY[el.dataset.tip]; if (!g) return;
    if (!box) { box = document.createElement("div"); box.id = "tipbox"; box.setAttribute("role", "tooltip"); document.body.appendChild(box); }
    box.innerHTML = `<b>${esc(g[0])}</b>${esc(g[1])}`;
    box.classList.add("on");
    const r = el.getBoundingClientRect(), bw = box.offsetWidth, bh = box.offsetHeight;
    let x = Math.min(Math.max(8, r.left + r.width / 2 - bw / 2), innerWidth - bw - 8);
    let y = r.bottom + 8; if (y + bh > innerHeight - 8) y = r.top - bh - 8;
    box.style.left = x + "px"; box.style.top = Math.max(8, y) + "px";
  };
  const hide = () => { if (box) box.classList.remove("on"); pinned = null; };
  document.addEventListener("mouseover", e => { const t = e.target.closest?.(".tip"); if (t) show(t); });
  document.addEventListener("mouseout", e => { if (e.target.closest?.(".tip") && !pinned) hide(); });
  document.addEventListener("focusin", e => { const t = e.target.closest?.(".tip"); if (t) show(t); });
  document.addEventListener("focusout", e => { if (e.target.closest?.(".tip")) hide(); });
  document.addEventListener("click", e => {
    const t = e.target.closest?.(".tip");
    if (t) { e.preventDefault(); e.stopPropagation(); if (pinned === t) hide(); else { show(t); pinned = t; } }
    else if (pinned) hide();
  }, true);
  addEventListener("scroll", () => { if (!pinned) hide(); }, { passive: true });
  document.addEventListener("keydown", e => { if (e.key === "Escape") hide(); });
})();

// ---- recently viewed proteins ---------------------------------------------
function getRecent() {
  try { return JSON.parse(localStorage.getItem("pa_recent") || "[]"); } catch (e) { return []; }
}
function addRecent(acc, name) {
  if (!acc) return;
  const list = getRecent().filter(r => r.acc !== acc);
  list.unshift({ acc, name: (name || "").slice(0, 60) });
  try { localStorage.setItem("pa_recent", JSON.stringify(list.slice(0, 12))); } catch (e) {}
}
// Renders chips into `el`; onPick(acc) is called on click (default: open the protein page).
function renderRecent(el, onPick, label = "Recent:") {
  const list = getRecent();
  if (!list.length) { el.innerHTML = ""; return; }
  el.innerHTML = `<span class="muted">${esc(label)}</span>` + list.map(r =>
    `<button type="button" class="pill" data-acc="${esc(r.acc)}" title="${esc(r.name)}">${esc(r.acc)}${r.name ? " · " + esc(r.name.slice(0, 22)) : ""}</button>`).join("")
    + `<button type="button" class="pill" data-clear="1" title="Forget recent proteins">✕</button>`;
  el.onclick = e => {
    const b = e.target.closest("button"); if (!b) return;
    if (b.dataset.clear) { try { localStorage.removeItem("pa_recent"); } catch (x) {} el.innerHTML = ""; return; }
    (onPick || (acc => location.href = "/protein/" + acc))(b.dataset.acc);
  };
}

// ---- charts -------------------------------------------------------------
function barChart(labels, values, { height = 200, unit = "%" } = {}) {
  const W = 640, H = height, L = 36, B = 22, T = 8;
  const max = Math.max(...values, 1e-9), bw = (W - L - 8) / values.length;
  const ticks = niceTicks(0, max, 4);
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img">`;
  for (const t of ticks) {
    const y = T + (H - T - B) * (1 - t / ticks[ticks.length - 1]);
    s += `<line class="axis" x1="${L}" x2="${W}" y1="${y}" y2="${y}"/><text x="${L - 4}" y="${y + 4}" text-anchor="end">${+t.toFixed(1)}</text>`;
  }
  const top = ticks[ticks.length - 1];
  values.forEach((v, i) => {
    const h = (H - T - B) * v / top, x = L + i * bw + 2;
    s += `<rect class="bar" x="${x}" y="${H - B - h}" width="${Math.max(bw - 4, 1)}" height="${h}"><title>${esc(labels[i])}: ${v.toFixed(2)}${unit}</title></rect>`;
    s += `<text x="${x + (bw - 4) / 2}" y="${H - 6}" text-anchor="middle">${esc(labels[i])}</text>`;
  });
  return s + "</svg>";
}

const HYD = { W: 760, H: 230, L: 36, B: 26, T: 8 };
function hydropathyChart(profile, window, segments, threshold = 1.6) {
  if (!profile.length) return `<p class="muted">Sequence shorter than the ${window}-residue window; no profile.</p>`;
  const { W, H, L, B, T } = HYD, half = Math.floor(window / 2);
  const n = profile.length + window - 1;
  const ymin = Math.min(-4, ...profile), ymax = Math.max(4, ...profile);
  const X = r => L + (W - L - 8) * (r - 1) / Math.max(n - 1, 1);
  const Y = v => T + (H - T - B) * (ymax - v) / (ymax - ymin);
  let s = `<svg class="chart hyd" viewBox="0 0 ${W} ${H}" role="img" aria-label="Hydropathy plot" data-n="${n}" data-half="${half}" data-ymin="${ymin}" data-ymax="${ymax}">`;
  segments.forEach(([a, b], i) => s += `<rect class="tm" data-seg="${i}" x="${X(a)}" y="${T}" width="${Math.max(X(b) - X(a), 2)}" height="${H - T - B}"><title>TM segment ${a}–${b} (click to highlight)</title></rect>`);
  for (const t of [-4, -2, 0, 2, 4]) s += `<line class="axis" x1="${L}" x2="${W}" y1="${Y(t)}" y2="${Y(t)}"/><text x="${L - 4}" y="${Y(t) + 4}" text-anchor="end">${t}</text>`;
  s += `<line class="thr" x1="${L}" x2="${W}" y1="${Y(threshold)}" y2="${Y(threshold)}"/>`;
  const pts = profile.map((v, i) => `${X(i + half + 1).toFixed(1)},${Y(v).toFixed(1)}`).join(" ");
  s += `<polyline class="line" points="${pts}"/>`;
  for (const t of niceTicks(1, n, 6)) s += `<text x="${X(t)}" y="${H - 8}" text-anchor="middle">${Math.round(t)}</text>`;
  s += `<line class="cursor" x1="0" x2="0" y1="${T}" y2="${H - B}" visibility="hidden"/><text class="readout" x="${L + 6}" y="${T + 14}"></text>`;
  return s + "</svg>";
}

// Hover read-out on the hydropathy plot; click a TM region -> onSegment(start, end).
function wireHydropathy(svg, r, onSegment) {
  if (!svg) return;
  const { W, L, T, H, B } = HYD, n = +svg.dataset.n, half = +svg.dataset.half;
  const cur = $(".cursor", svg), out = $(".readout", svg);
  const toRes = evt => {
    const pt = svg.createSVGPoint(); pt.x = evt.clientX; pt.y = evt.clientY;
    const p = pt.matrixTransform(svg.getScreenCTM().inverse());
    return Math.round(1 + (p.x - L) / (W - L - 8) * (n - 1));
  };
  svg.addEventListener("pointermove", e => {
    const res = toRes(e);
    if (res < 1 || res > n) { cur.setAttribute("visibility", "hidden"); out.textContent = ""; return; }
    const x = L + (W - L - 8) * (res - 1) / Math.max(n - 1, 1);
    cur.setAttribute("x1", x); cur.setAttribute("x2", x); cur.setAttribute("visibility", "visible");
    const v = r.hydropathy_profile[res - half - 1];
    out.textContent = `Residue ${res} (${r.sequence[res - 1]})` + (v !== undefined ? `: hydropathy ${v.toFixed(2)}` : " (no full window here)");
  });
  svg.addEventListener("pointerleave", () => { cur.setAttribute("visibility", "hidden"); out.textContent = ""; });
  svg.addEventListener("click", e => {
    const rect = e.target.closest(".tm"); if (!rect) return;
    $$(".tm", svg).forEach(x => x.classList.toggle("on", x === rect));
    const [a, b] = r.tm_segments[+rect.dataset.seg];
    onSegment && onSegment(a, b);
  });
}

function niceTicks(lo, hi, n) {
  const span = hi - lo || 1, raw = span / n, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map(m => m * mag).find(s => s >= raw);
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 0.999; t += step) out.push(t);
  return out;
}

// ---- sequence viewer --------------------------------------------------------
// ranges: [{start, end, cls}] (1-based inclusive). 60 residues per line, blocks of 10.
function sequenceHtml(seq, ranges = []) {
  const cls = new Array(seq.length).fill("");
  for (const r of ranges) for (let i = Math.max(1, r.start); i <= Math.min(seq.length, r.end); i++) cls[i - 1] = r.cls;
  const w = String(seq.length).length;
  let out = "";
  for (let i = 0; i < seq.length; i += 60) {
    out += `<span class="n">${String(i + 1).padStart(w)} </span>`;
    for (let j = i; j < Math.min(i + 60, seq.length); ) {
      if ((j - i) % 10 === 0 && j > i) out += " ";
      let k = j; while (k < Math.min(i + 60, seq.length) && cls[k] === cls[j] && (k === j || (k - i) % 10 !== 0)) k++;
      const chunk = seq.slice(j, k);
      out += cls[j] ? `<span class="${cls[j]}" data-pos="${j + 1}">${chunk}</span>` : chunk;
      j = k;
    }
    out += "\n";
  }
  return out;
}
function highlightSequence(view, seq, base, start, end) {
  view.innerHTML = sequenceHtml(seq, [...base, { start, end, cls: "hl" }]);
  const first = $(".hl", view);
  if (first) view.scrollTop = Math.max(0, first.offsetTop - view.offsetTop - 30);
}

// ---- analysis rendering (Analyze and Protein detail pages) -------------------
function renderAnalysis(r, { showSequence = true } = {}) {
  const comp = Object.entries(r.composition);
  const tms = r.tm_segments.length ? r.tm_segments.map(([a, b]) => `${a}–${b}`).join(", ") : "none predicted";
  return `
  <div class="stats">
    <div class="stat"><div class="k">Length${tip("length")}</div><div class="v">${r.length.toLocaleString()} <small>aa</small></div></div>
    <div class="stat"><div class="k">Mol. weight${tip("mw")}</div><div class="v">${fmt(r.molecular_weight / 1000, 2)} <small>kDa</small></div></div>
    <div class="stat"><div class="k">GRAVY${tip("gravy")}</div><div class="v">${fmt(r.gravy, 3)}</div><small class="muted">${esc(r.gravy_class)}</small></div>
    <div class="stat"><div class="k">pI${tip("pi")}</div><div class="v">${fmt(r.isoelectric_point, 2)}</div></div>
    <div class="stat"><div class="k">Charge pH 7${tip("charge")}</div><div class="v">${r.charge_ph7 > 0 ? "+" : ""}${fmt(r.charge_ph7, 2)}</div></div>
    <div class="stat"><div class="k">Aromaticity${tip("aromaticity")}</div><div class="v">${fmt(r.aromaticity, 3)}</div></div>
    <div class="stat"><div class="k">TM segments${tip("tm")}</div><div class="v">${r.n_tm_segments}</div></div>
  </div>
  ${r.other_residues ? `<div class="msg warn">${r.other_residues} non-standard/ambiguous residue(s) (B, Z, X, J, U, O) counted in length/MW only.</div>` : ""}
  <h3>Hydropathy${tip("hydropathy")} <span class="muted small" style="font-weight:400">Kyte–Doolittle, window ${r.hydropathy_window}; shaded = predicted TM: ${esc(tms)}</span></h3>
  ${hydropathyChart(r.hydropathy_profile, r.hydropathy_window, r.tm_segments)}
  ${showSequence ? `<h3>Sequence <span class="muted small" style="font-weight:400">shaded = predicted TM; click a region in the plot to highlight it</span></h3>
  <div class="seqview">${sequenceHtml(r.sequence, r.tm_segments.map(([a, b]) => ({ start: a, end: b, cls: "tmr" })))}</div>` : ""}
  <h3>Amino-acid composition (%)${tip("composition")}</h3>
  ${barChart(comp.map(c => c[0]), comp.map(c => c[1].percent))}
  <details><summary>Composition table</summary><div class="table-wrap"><table>
    <thead><tr><th>Residue</th>${comp.map(c => `<th>${c[0]}</th>`).join("")}</tr></thead>
    <tbody><tr><td>Count</td>${comp.map(c => `<td>${c[1].count}</td>`).join("")}</tr>
    <tr><td>%</td>${comp.map(c => `<td>${c[1].percent.toFixed(1)}</td>`).join("")}</tr></tbody></table></div></details>`;
}

// Render + wire interactivity. Returns {highlight(start, end)}.
function mountAnalysis(el, r, opts = {}) {
  el.innerHTML = renderAnalysis(r, opts);
  const view = $(".seqview", el);
  const base = r.tm_segments.map(([a, b]) => ({ start: a, end: b, cls: "tmr" }));
  const highlight = (a, b) => { if (view) highlightSequence(view, r.sequence, base, a, b); opts.onHighlight && opts.onHighlight(a, b); };
  wireHydropathy($("svg.hyd", el), r, highlight);
  return { highlight };
}

function analysisCsvParts(records) {
  const aas = "ACDEFGHIKLMNPQRSTVWY".split("");
  return [["name", "length", "molecular_weight", "gravy", "gravy_class", "isoelectric_point", "charge_ph7", "aromaticity", "n_tm_segments", "tm_segments", ...aas.map(a => "pct_" + a)],
    records.map(r => [r.name, r.length, r.molecular_weight, r.gravy, r.gravy_class, r.isoelectric_point, r.charge_ph7, r.aromaticity, r.n_tm_segments,
      r.tm_segments.map(s => s.join("-")).join(";"), ...aas.map(a => r.composition[a].percent)])];
}

// ---- feature tracks (Protein detail) ----------------------------------------
const TRACKS = [
  { name: "TM (UniProt)", types: ["Transmembrane", "Intramembrane"], color: "#f59e0b" },
  { name: "TM (this app)", predicted: true, color: "#fbbf24" },
  { name: "Signal peptide", types: ["Signal"], color: "#ef4444" },
  { name: "Topology", types: ["Topological domain"], color: "#94a3b8" },
  { name: "Domains", types: ["Domain", "Zinc finger", "DNA binding"], color: "#3b82f6" },
  { name: "Regions / motifs", types: ["Region", "Motif"], color: "#8b5cf6" },
  { name: "Sites", types: ["Active site", "Binding site"], color: "#10b981" },
];
// Returns {html, items}; items[i] = {start, end, label} for click handlers (data-i).
function featureTracks(length, features, predictedTm) {
  const W = 760, L = 120, R = 8, rowH = 22, X = p => L + (W - L - R) * (p - 1) / Math.max(length - 1, 1);
  const items = []; let rows = "", y = 4;
  for (const t of TRACKS) {
    const feats = t.predicted ? predictedTm.map(([a, b]) => ({ type: "Predicted TM", start: a, end: b, description: "Kyte–Doolittle ≥ 1.6" }))
      : features.filter(f => t.types.includes(f.type));
    if (!feats.length && !t.predicted && !t.types.includes("Transmembrane")) continue;
    rows += `<text x="0" y="${y + 15}">${esc(t.name)}</text><rect class="lane" x="${L}" y="${y + 4}" width="${W - L - R}" height="${rowH - 8}" rx="3"/>`;
    for (const f of feats) {
      const i = items.length;
      const label = `${f.type}${f.description ? ": " + f.description : ""}${f.ligand ? " (" + f.ligand + ")" : ""} ${f.start}–${f.end}`;
      items.push({ start: f.start, end: f.end, label, type: f.type, desc: (f.description || "") + (f.ligand ? ` (${f.ligand})` : "") });
      const x = X(f.start), w = Math.max(X(f.end) - x, 3);
      rows += `<rect class="feat" data-i="${i}" x="${x}" y="${y + 3}" width="${w}" height="${rowH - 6}" rx="2" fill="${t.color}" fill-opacity="${t.types && t.types.includes("Region") && /disorder/i.test(f.description) ? 0.45 : 0.9}"><title>${esc(label)} (click to highlight)</title></rect>`;
    }
    y += rowH;
  }
  for (const t of niceTicks(1, length, 6)) rows += `<text x="${X(t)}" y="${y + 14}" text-anchor="middle">${Math.round(t)}</text>`;
  return { html: `<svg class="tracks" viewBox="0 0 ${W} ${y + 20}" role="img" aria-label="Sequence annotations">${rows}</svg>`, items };
}

// ---- save to project ---------------------------------------------------------
// saveToProject("mutation", "HBB E7V", {...}) opens a small dialog to pick or create a project.
async function saveToProject(kind, title, data) {
  if (window.PA_PUBLIC) { alert("Projects are turned off on this shared server. Run the app on your own computer to save work."); return; }
  let dlg = $("#savedlg");
  if (!dlg) {
    dlg = document.createElement("dialog"); dlg.id = "savedlg"; dlg.className = "card";
    dlg.innerHTML = `<form method="dialog" style="min-width:min(380px,80vw)">
      <h2>Save to project</h2><p class="muted small" id="savewhat"></p>
      <label class="field">Project <select id="saveproj"></select></label>
      <label class="field" id="newwrap" style="margin-top:6px">New project name <input type="text" id="newproj" placeholder="e.g. Haemoglobin coursework"></label>
      <label class="field" style="margin-top:6px">Note (optional) <textarea id="savenote" style="min-height:60px"></textarea></label>
      <div id="savemsg"></div>
      <div class="row" style="margin-top:10px; justify-content:flex-end"><button value="cancel" class="secondary" formnovalidate>Cancel</button><button id="savego" value="ok">Save</button></div></form>`;
    document.body.appendChild(dlg);
    $("#saveproj", dlg).onchange = () => { $("#newwrap", dlg).style.display = $("#saveproj", dlg).value ? "none" : ""; };
  }
  $("#savewhat", dlg).textContent = `${title} (${kind})`;
  $("#savenote", dlg).value = ""; showMsg($("#savemsg", dlg), "");
  try {
    const { projects } = await api("/api/projects");
    let last = null; try { last = localStorage.getItem("pa_last_project"); } catch (e) {}
    $("#saveproj", dlg).innerHTML = projects.map(p => `<option value="${p.id}" ${p.id === last ? "selected" : ""}>${esc(p.name)} (${p.n_items})</option>`).join("") + `<option value="">+ New project…</option>`;
    if (!projects.length) $("#saveproj", dlg).value = "";
    $("#saveproj", dlg).onchange();
  } catch (e) { alert(e.message); return; }
  $("#savego", dlg).onclick = async ev => {
    ev.preventDefault();
    try {
      let pid = $("#saveproj", dlg).value;
      if (!pid) pid = (await api("/api/projects", { json: { name: $("#newproj", dlg).value } })).id;
      await api(`/api/projects/${pid}/items`, { json: { kind, title, data, note: $("#savenote", dlg).value } });
      try { localStorage.setItem("pa_last_project", pid); } catch (e) {}
      dlg.close();
      toast(`Saved to project. <a href="/projects#${pid}">Open project</a>`);
    } catch (e) { showMsg($("#savemsg", dlg), e.message); }
  };
  dlg.showModal();
}
function saveButton(label = "Save to project") {
  return window.PA_PUBLIC ? "" : `<button type="button" class="secondary save-btn">${esc(label)}</button>`;
}
function toast(html) {
  const t = document.createElement("div"); t.className = "toast"; t.innerHTML = html;
  document.body.appendChild(t); setTimeout(() => t.classList.add("on"), 10);
  setTimeout(() => { t.classList.remove("on"); setTimeout(() => t.remove(), 300); }, 4000);
}

// ---- region importance track (explanations) -----------------------------------
function regionTrack(length, regions) {
  const W = 760, H = 90, L = 40, T = 10, B = 22;
  const max = Math.max(0.01, ...regions.map(r => Math.abs(r.drop)));
  const X = p => L + (W - L - 8) * (p - 1) / Math.max(length - 1, 1);
  const Y = v => T + (H - T - B) / 2 * (1 - v / max);
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="Region importance">`;
  s += `<line class="axis" x1="${L}" x2="${W}" y1="${Y(0)}" y2="${Y(0)}"/>`;
  for (const r of regions) {
    const x = X(r.start), w = Math.max(X(r.end) - x, 2), y = Y(Math.max(r.drop, 0)), y0 = Y(0);
    const neg = r.drop < 0;
    s += `<rect x="${x}" y="${neg ? y0 : y}" width="${w}" height="${Math.abs(Y(r.drop) - y0)}" fill="${neg ? "#94a3b8" : "#ef4444"}" fill-opacity=".8"><title>Residues ${r.start}–${r.end}: hiding them changes the probability by ${(-100 * r.drop).toFixed(1)} points</title></rect>`;
  }
  s += `<text x="${L - 4}" y="${T + 8}" text-anchor="end">+${(100 * max).toFixed(0)}</text><text x="${L - 4}" y="${Y(0) + 4}" text-anchor="end">0</text>`;
  for (const t of niceTicks(1, length, 6)) s += `<text x="${X(t)}" y="${H - 6}" text-anchor="middle">${Math.round(t)}</text>`;
  return s + "</svg>";
}

// ---- phylogenetic tree (rectangular phylogram) ----------------------------------
function treeSvg(tree) {
  const leaves = []; let maxDepth = 0;
  (function walk(n, d) { n._d = d; maxDepth = Math.max(maxDepth, d); if (!n.children.length) leaves.push(n); n.children.forEach(c => walk(c, d + c.length)); })(tree, 0);
  const rowH = 24, W = 760, labelW = 230, L = 10, H = leaves.length * rowH + 30;
  const X = d => L + (W - labelW - L - 10) * d / (maxDepth || 1);
  leaves.forEach((n, i) => n._y = 14 + i * rowH);
  (function setY(n) { if (n.children.length) { n.children.forEach(setY); n._y = (n.children[0]._y + n.children[n.children.length - 1]._y) / 2; } })(tree);
  let s = `<svg class="chart tree" viewBox="0 0 ${W} ${H}" role="img" aria-label="Family tree">`;
  (function draw(n) {
    for (const c of n.children) {
      s += `<path d="M${X(n._d)},${c._y} H${X(c._d)}" stroke="var(--text)" fill="none" stroke-width="1.4"><title>branch length ${c.length}</title></path>`;
      draw(c);
    }
    if (n.children.length) s += `<line x1="${X(n._d)}" x2="${X(n._d)}" y1="${n.children[0]._y}" y2="${n.children[n.children.length - 1]._y}" stroke="var(--text)" stroke-width="1.4"/>`;
    else s += `<text x="${X(n._d) + 6}" y="${n._y + 4}" style="fill:var(--text);font-size:12px">${esc(n.name)}</text>`;
  })(tree);
  const bar = niceTicks(0, maxDepth / 3, 1)[1] || maxDepth / 4;
  s += `<line x1="${L}" x2="${X(bar)}" y1="${H - 8}" y2="${H - 8}" stroke="var(--muted)"/><text x="${X(bar) + 4}" y="${H - 4}">${+bar.toFixed(2)} substitutions/site</text>`;
  return s + "</svg>";
}

// ---- quick-jump palette ("/" or Ctrl+K) ---------------------------------------------
const PAGES = [
  ["Home", "/", "start page"], ["Analyze a sequence", "/analyze", "properties, hydropathy"], ["Search UniProt", "/search", "by name or keyword"],
  ["BLAST: find similar sequences", "/blast", "search all of UniProt"], ["Compare two proteins", "/compare", "pairwise alignment"],
  ["Align many proteins + family tree", "/align", "multiple alignment"], ["Compare 3D structures", "/structures", "overlay, TM-score"],
  ["Mutation explorer", "/mutate", "AlphaMissense, variants"], ["Lab tools: DNA → protein", "/lab#translate", "ORFs"],
  ["Lab tools: protease digest", "/lab#digest", "peptide masses"], ["Lab tools: A280 concentration", "/lab#a280", "extinction coefficient"],
  ["Classify protein function", "/classify", "machine learning"], ["Dataset & ML", "/ml", "train, evaluate"], ["Learn: lessons", "/learn", "quizzes"],
  ["Projects (saved work)", "/projects", "notes, reports"], ["Methods & citations", "/methods", "how to cite"],
];
const AA_SEQ_RE = /^[ACDEFGHIKLMNPQRSTVWYXBZUO\s]{20,}$/i;
function paletteItems(q) {
  const s = q.trim(), items = [];
  if (s.startsWith(">") || (AA_SEQ_RE.test(s) && /[EFILPQ]/i.test(s))) items.push(["Analyze this sequence", () => handoff("/analyze", s), "sequence"]);
  if (/^[ACGTUN\s]{30,}$/i.test(s)) items.push(["Translate this DNA", () => handoff("/lab#translate", s, "#dna"), "DNA"]);
  if (ACC_RE.test(s)) {
    items.push([`Open protein ${s.toUpperCase()}`, () => location.href = "/protein/" + s.toUpperCase(), "UniProt entry"]);
    items.push([`Mutations of ${s.toUpperCase()}`, () => location.href = "/mutate?acc=" + s.toUpperCase(), ""]);
  }
  if (s && !s.startsWith(">") && s.length < 200) items.push([`Search UniProt for “${s}”`, () => location.href = "/search?q=" + encodeURIComponent(s), "search"]);
  const low = s.toLowerCase();
  PAGES.filter(([t, , d]) => !low || (t + " " + d).toLowerCase().includes(low)).forEach(([t, u, d]) => items.push([t, () => location.href = u, d]));
  if (!s) getRecent().slice(0, 5).forEach(r => items.push([`${r.acc} · ${r.name}`, () => location.href = "/protein/" + r.acc, "recent"]));
  return items.slice(0, 12);
}
function handoff(url, text, sel) {
  try { sessionStorage.setItem("pa_handoff", text); } catch (e) {}
  location.href = url.includes("#") ? url.replace("#", "?handoff=1#") : url + "?handoff=1";
}
function openPalette() {
  let d = $("#palette");
  if (!d) {
    d = document.createElement("dialog"); d.id = "palette";
    d.innerHTML = `<div class="box"><input type="text" id="pq" placeholder="Type an accession (P69905), a protein name, a sequence, or a tool…" aria-label="Quick search">
      <ul id="pl" role="listbox"></ul><div class="foot">↑ ↓ to choose · Enter to open · Esc to close</div></div>`;
    document.body.appendChild(d);
    d.addEventListener("click", e => { if (e.target === d) d.close(); });
    let sel = 0, items = [];
    const draw = () => {
      items = paletteItems($("#pq").value); sel = Math.min(sel, Math.max(items.length - 1, 0));
      $("#pl").innerHTML = items.map(([t, , h], i) => `<li role="option" data-i="${i}" class="${i === sel ? "on" : ""}">${esc(t)}<small>${esc(h)}</small></li>`).join("");
    };
    $("#pq").addEventListener("input", () => { sel = 0; draw(); });
    $("#pq").addEventListener("keydown", e => {
      if (e.key === "ArrowDown") { sel = Math.min(sel + 1, items.length - 1); draw(); e.preventDefault(); }
      else if (e.key === "ArrowUp") { sel = Math.max(sel - 1, 0); draw(); e.preventDefault(); }
      else if (e.key === "Enter" && items[sel]) { d.close(); items[sel][1](); }
    });
    $("#pl").addEventListener("click", e => { const li = e.target.closest("li"); if (li) { d.close(); items[+li.dataset.i][1](); } });
    d._draw = draw;
  }
  $("#pq").value = ""; d._draw(); d.showModal(); $("#pq").focus();
}

// ---- autosaved drafts: what you type is kept if you leave the page ---------------------
const DRAFT_SELECTOR = "main textarea:not([data-nodraft]):not(#pnotes), main input[type=text]:not([data-nodraft]):not(#newname):not(#pname):not(#email):not(#homeq)";
function initDrafts() {
  const page = location.pathname;
  const els = $$(DRAFT_SELECTOR).filter(el => el.id && !el.closest("#palette, dialog, .item"));
  let restored = 0;
  const params = new URLSearchParams(location.search);
  const fromUrl = params.toString() !== "";  // a link that pre-fills inputs wins over drafts
  for (const el of els) {
    const key = `pa_draft:${page}:${el.id}`;
    let v = null; try { v = localStorage.getItem(key); } catch (e) {}
    if (v && !el.value && !fromUrl) { el.value = v; restored++; }
    el.addEventListener("input", () => { try { el.value ? localStorage.setItem(key, el.value.slice(0, 200000)) : localStorage.removeItem(key); } catch (e) {} });
  }
  const first = els.find(e => e.tagName === "TEXTAREA") || els[0];
  if (restored && first) {
    const note = document.createElement("p"); note.className = "draftnote noprint";
    note.innerHTML = `Restored what you typed last time. <button type="button">Clear</button>`;
    first.insertAdjacentElement("afterend", note);
    $("button", note).onclick = () => { els.forEach(el => { el.value = ""; el.dispatchEvent(new Event("input")); }); note.remove(); };
  }
}

// ---- copy buttons on sequence / alignment blocks --------------------------------------
function addCopyButtons(root = document) {
  $$(".seqview, .aln, .msa, pre", root).forEach(el => {
    if (el.dataset.copy || el.closest(".legend")) return;
    el.dataset.copy = "1";
    const wrap = document.createElement("div"); wrap.className = "copywrap";
    el.parentNode.insertBefore(wrap, el); wrap.appendChild(el);
    const b = document.createElement("button"); b.type = "button"; b.className = "secondary copybtn noprint"; b.textContent = "Copy";
    b.onclick = async () => {
      const text = el.classList.contains("seqview") ? el.innerText.replace(/[\s\d]+/g, "") : el.innerText;
      try { await navigator.clipboard.writeText(text); b.textContent = "Copied ✓"; } catch (e) { b.textContent = "Select & copy"; }
      setTimeout(() => b.textContent = "Copy", 1500);
    };
    wrap.appendChild(b);
  });
}

// ---- "More ▾" menus ---------------------------------------------------------------
document.addEventListener("click", e => {
  const t = e.target.closest(".menu > button");
  $$(".menu.open").forEach(m => { if (!t || m !== t.parentElement) m.classList.remove("open"); });
  if (t) { t.parentElement.classList.toggle("open"); t.setAttribute("aria-expanded", t.parentElement.classList.contains("open")); }
  else if (e.target.closest(".menu-list")) $$(".menu.open").forEach(m => m.classList.remove("open"));
});

document.addEventListener("DOMContentLoaded", () => {
  updateThemeButton();
  $$("[data-theme-set]").forEach(b => b.onclick = () => setTheme(b.dataset.themeSet));
  const mb = $("#menubtn");
  if (mb) mb.onclick = () => { const h = $("#hdr"); h.classList.toggle("open"); mb.setAttribute("aria-expanded", h.classList.contains("open")); };
  const qb = $("#qbtn"); if (qb) qb.onclick = openPalette;
  document.addEventListener("keydown", e => {
    const typing = /^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement?.tagName) || document.activeElement?.isContentEditable;
    if ((e.key === "/" && !typing) || (e.key.toLowerCase() === "k" && (e.ctrlKey || e.metaKey))) { e.preventDefault(); openPalette(); }
    if (e.key === "Escape") { $("#hdr")?.classList.remove("open"); $$(".menu.open").forEach(m => m.classList.remove("open")); }
    // Ctrl+Enter in a text box runs that card's main action
    if (e.key === "Enter" && (e.ctrlKey || e.metaKey) && typing) {
      const card = document.activeElement.closest(".card, section, form");
      const btn = card && $$("button:not(.secondary):not(.tip):not(.copybtn):not(.danger):not(.pill):not(.dismiss)", card)
        .find(b => !b.disabled && b.offsetParent && !b.closest(".draftnote, .recent, .menu-list"));
      if (btn) { e.preventDefault(); btn.click(); }
    }
  });
  const tt = $("#totop");
  if (tt) { tt.onclick = () => scrollTo({ top: 0 }); addEventListener("scroll", () => tt.classList.toggle("on", scrollY > 700), { passive: true }); }
  // Lab tools hand a sequence to Analyze / Classify / BLAST via sessionStorage (?handoff=1).
  if (new URLSearchParams(location.search).get("handoff")) {
    let t = null; try { t = sessionStorage.getItem("pa_handoff"); sessionStorage.removeItem("pa_handoff"); } catch (e) {}
    const box = $("#seq") || $("#t") || $("#q") || $("#dna");
    if (t && box) { box.value = t; box.focus(); }
  }
  initDrafts();
  // copy buttons for blocks rendered later
  new MutationObserver(() => addCopyButtons()).observe(document.querySelector("main"), { childList: true, subtree: true });
  addCopyButtons();
});
