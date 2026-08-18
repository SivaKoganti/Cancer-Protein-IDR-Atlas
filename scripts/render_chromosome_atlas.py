#!/usr/bin/env python3
"""render_chromosome_atlas.py

Builds results/map/chromosome_atlas.html — a self-contained interactive
human karyotype displaying 50 CGC cancer genes at their GRCh38 loci.
Gene markers are color-coded by atlas metrics; clicking opens a full
per-residue IDR / LLPS / conservation detail panel linked to NCBI Entrez.
"""

import argparse
import json
from pathlib import Path

import pandas as pd


# ── GRCh38 reference data ─────────────────────────────────────────────────────

CHROM_SIZES = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16":  90338345,
    "17":  83257441, "18":  80373285, "19":  58617616, "20":  64444167,
    "21":  46709983, "22":  50818468, "X":  156040895, "Y":   57227415,
}

# Centromere midpoints (bp, GRCh38 approximate)
CENTROMERES = {
    "1": 123400000, "2":  93900000, "3":  90900000, "4":  50400000,
    "5":  48400000, "6":  61000000, "7":  59900000, "8":  45200000,
    "9":  43400000, "10": 39800000, "11": 53700000, "12": 35800000,
    "13": 17700000, "14": 17200000, "15": 19000000, "16": 36600000,
    "17": 25100000, "18": 18500000, "19": 26200000, "20": 28100000,
    "21": 12000000, "22": 15000000, "X":  61000000, "Y":  10400000,
}

# GRCh38 gene loci (NCBI Entrez Gene, canonical isoform)
GENE_LOCI = {
    "TP53":    {"chrom": "17", "start":   7661779, "end":   7687538, "band": "17p13.1"},
    "KRAS":    {"chrom": "12", "start":  25205246, "end":  25250929, "band": "12p12.1"},
    "EGFR":    {"chrom":  "7", "start":  55019017, "end":  55207338, "band": "7p11.2"},
    "BRCA1":   {"chrom": "17", "start":  43044295, "end":  43170245, "band": "17q21.31"},
    "BRCA2":   {"chrom": "13", "start":  32315086, "end":  32400268, "band": "13q12.3"},
    "PIK3CA":  {"chrom":  "3", "start": 179148114, "end": 179240093, "band": "3q26.32"},
    "BRAF":    {"chrom":  "7", "start": 140719327, "end": 140924764, "band": "7q34"},
    "APC":     {"chrom":  "5", "start": 112707498, "end": 112846239, "band": "5q22.2"},
    "PTEN":    {"chrom": "10", "start":  89623195, "end":  89728532, "band": "10q23.31"},
    "NRAS":    {"chrom":  "1", "start": 114704469, "end": 114716771, "band": "1p13.2"},
    "CDKN2A":  {"chrom":  "9", "start":  21967752, "end":  22009455, "band": "9p21.3"},
    "RB1":     {"chrom": "13", "start":  48303747, "end":  49058273, "band": "13q14.2"},
    "ERBB2":   {"chrom": "17", "start":  39687914, "end":  39730426, "band": "17q12"},
    "SMAD4":   {"chrom": "18", "start":  51030214, "end":  51085042, "band": "18q21.2"},
    "STK11":   {"chrom": "19", "start":   1205798, "end":   1228434, "band": "19p13.3"},
    "ATM":     {"chrom": "11", "start": 108222484, "end": 108369102, "band": "11q22.3"},
    "CHEK2":   {"chrom": "22", "start":  28687743, "end":  28742311, "band": "22q12.1"},
    "MET":     {"chrom":  "7", "start": 116672196, "end": 116798377, "band": "7q31.2"},
    "VHL":     {"chrom":  "3", "start":  10183318, "end":  10195354, "band": "3p25.3"},
    "GNAS":    {"chrom": "20", "start":  58839718, "end":  58911462, "band": "20q13.32"},
    "IDH1":    {"chrom":  "2", "start": 209113112, "end": 209126214, "band": "2q34"},
    "IDH2":    {"chrom": "15", "start":  90631834, "end":  90645894, "band": "15q26.1"},
    "JAK2":    {"chrom":  "9", "start":   4984098, "end":   5128183, "band": "9p24.1"},
    "KIT":     {"chrom":  "4", "start":  54657928, "end":  54740715, "band": "4q12"},
    "PDGFRA":  {"chrom":  "4", "start":  54229097, "end":  54298245, "band": "4q12"},
    "ALK":     {"chrom":  "2", "start":  29192774, "end":  29921585, "band": "2p23.2"},
    "ROS1":    {"chrom":  "6", "start": 117288070, "end": 117639763, "band": "6q22.1"},
    "RET":     {"chrom": "10", "start":  43572517, "end":  43625797, "band": "10q11.21"},
    "NFE2L2":  {"chrom":  "2", "start": 177230300, "end": 177247972, "band": "2q31.2"},
    "PTCH1":   {"chrom":  "9", "start":  95449086, "end":  95633100, "band": "9q22.32"},
    "FBXW7":   {"chrom":  "4", "start": 153245685, "end": 153505879, "band": "4q31.3"},
    "CTNNB1":  {"chrom":  "3", "start":  41194946, "end":  41279232, "band": "3p22.1"},
    "MAP2K1":  {"chrom": "15", "start":  66436822, "end":  66509310, "band": "15q22.31"},
    "MAP2K2":  {"chrom": "19", "start":   4090931, "end":   4110806, "band": "19p13.3"},
    "NOTCH1":  {"chrom":  "9", "start": 136494433, "end": 136546169, "band": "9q34.3"},
    "NOTCH2":  {"chrom":  "1", "start": 119911553, "end": 120098325, "band": "1p12"},
    "ARID1A":  {"chrom":  "1", "start":  26696848, "end":  26991631, "band": "1p36.11"},
    "ARID2":   {"chrom": "12", "start":  45935503, "end":  46257539, "band": "12q12"},
    "SMARCA4": {"chrom": "19", "start":  11024460, "end":  11127487, "band": "19p13.2"},
    "KMT2A":   {"chrom": "11", "start": 118307205, "end": 118397539, "band": "11q23.3"},
    "KMT2D":   {"chrom": "12", "start":  49012899, "end":  49116520, "band": "12q13.12"},
    "TSC1":    {"chrom":  "9", "start": 132890955, "end": 132970246, "band": "9q34.13"},
    "TSC2":    {"chrom": "16", "start":   2047053, "end":   2143323, "band": "16p13.3"},
    "MLH1":    {"chrom":  "3", "start":  36993332, "end":  37050918, "band": "3p22.2"},
    "MSH2":    {"chrom":  "2", "start":  47403067, "end":  47630948, "band": "2p21"},
    "MSH6":    {"chrom":  "2", "start":  47783088, "end":  47806306, "band": "2p16.3"},
    "PMS2":    {"chrom":  "7", "start":   5969085, "end":   6009766, "band": "7p22.1"},
    "ERCC2":   {"chrom": "19", "start":  45348747, "end":  45382228, "band": "19q13.32"},
    "RNF43":   {"chrom": "17", "start":  56442622, "end":  56479949, "band": "17q22"},
    "NTRK1":   {"chrom":  "1", "start": 156815533, "end": 156881850, "band": "1q23.1"},
}


# ── Data loading ──────────────────────────────────────────────────────────────

def load_data(atlas_dir: Path, summary_path: Path):
    summary = pd.read_csv(summary_path, sep="\t")

    gene_data = {}
    for _, row in summary.iterrows():
        gene = str(row["gene"])
        locus = GENE_LOCI.get(gene, {})
        gene_data[gene] = {
            "chrom":  locus.get("chrom", "?"),
            "band":   locus.get("band", "?"),
            "start":  locus.get("start", 0),
            "end":    locus.get("end", 0),
            "length": int(row["length"]),
            "mean_iupred":             round(float(row["mean_iupred"]), 4),
            "mean_llps":               round(float(row["mean_llps"]), 4),
            "mean_conservation":       round(float(row["mean_conservation"]), 4),
            "mean_structural":         round(float(row["mean_structural"]), 4),
            "low_complexity_fraction": round(float(row["low_complexity_fraction"]), 4),
            "variant_count":           int(row["variant_count"]),
            "virus_interaction_fraction": round(float(row.get("virus_interaction_fraction", 0.0)), 4),
            "mean_vipp_score":         round(float(row.get("mean_vipp_score", 0.0)), 4),
            "max_vipp_score":          round(float(row.get("max_vipp_score", 0.0)), 4),
        }

    residue_data = {}
    for gene in gene_data:
        for suffix in ("_atlas_with_clinvar.tsv", "_atlas.tsv"):
            p = atlas_dir / f"{gene}{suffix}"
            if not p.exists():
                continue
            df = pd.read_csv(p, sep="\t")
            entry = {
                "pos":          df["pos"].tolist(),
                "aa":           df["aa"].tolist(),
                "iupred":       [round(float(v), 3) for v in df["iupred_score"]],
                "is_idr":       [int(v) for v in df["is_idr"]] if "is_idr" in df.columns else [int(float(v) >= 0.5) for v in df["iupred_score"]],
                "idr_class":    df["idr_class"].tolist() if "idr_class" in df.columns else ["disordered" if float(v) >= 0.5 else "structured" for v in df["iupred_score"]],
                "plddt":        [round(float(v), 1) if pd.notna(v) else None for v in df["plddt"]] if "plddt" in df.columns else [None] * len(df),
                "llps":         [round(float(v), 3) for v in df["llps_proxy"]],
                "structural":   [round(float(v), 3) for v in df["structural_proxy"]],
                "conservation": [round(float(v), 3) for v in df["conservation"]],
                "low_cx":       [int(v) for v in df["low_complexity"]],
                "virus_interaction": [int(v) for v in df["virus_interaction"]] if "virus_interaction" in df.columns else [0] * len(df),
                "vipp":         [round(float(v), 3) for v in df["vipp_score"]] if "vipp_score" in df.columns else [0.0] * len(df),
            }
            if "clinvar_count" in df.columns:
                entry["cv_count"] = [int(v) for v in df["clinvar_count"]]
                entry["cv_sig"]   = df["clinvar_significance"].fillna("").tolist()
            residue_data[gene] = entry
            break

    return gene_data, residue_data


# ── HTML generation ───────────────────────────────────────────────────────────

_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Cancer Protein IDR Atlas — Chromosome Map</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;background:#eef1f5;
  color:#222;height:100vh;overflow:hidden;display:flex;flex-direction:column}
header{background:linear-gradient(135deg,#1e3a5f 0%,#2d5a8e 100%);color:#fff;
  padding:9px 18px;flex-shrink:0;display:flex;align-items:center;gap:20px;flex-wrap:wrap}
.title-row h1{font-size:16px;font-weight:800;letter-spacing:.3px}
.title-row p{font-size:10px;opacity:.72;margin-top:1px}
.controls{display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-left:auto}
.controls label{font-size:11px;display:flex;align-items:center;gap:5px}
.controls select,.controls input{font-size:11px;padding:3px 7px;border:none;border-radius:4px;
  background:rgba(255,255,255,.18);color:#fff;cursor:pointer}
.controls select option{color:#222;background:#fff}
.controls input::placeholder{color:rgba(255,255,255,.55)}
#main{display:flex;flex:1;overflow:hidden;min-height:0}
#karyo-panel{background:#fff;border-right:1px solid #dde2ea;overflow:auto;
  padding:10px 6px;width:500px;flex-shrink:0}
#detail-panel{flex:1;overflow-y:auto;padding:16px 20px;background:#fff;min-width:0}
#placeholder{height:100%;display:flex;align-items:center;justify-content:center;
  color:#aab;font-size:13px;text-align:center;flex-direction:column;gap:8px}
.chrom-arm{fill:#cdd3e0;stroke:#9aa3b8;stroke-width:.5}
.centromere{fill:#8891a8}
.chrom-label{font-size:10px;fill:#667;font-weight:500}
.gene-marker{cursor:pointer;transition:opacity .15s}
.gene-marker circle{transition:stroke .1s,stroke-width .1s}
.gene-lbl-sm{font-size:7.5px;fill:#445;pointer-events:none}
.gene-header{display:flex;align-items:baseline;gap:10px;border-bottom:2px solid #1e3a5f;
  padding-bottom:10px;margin-bottom:14px;flex-wrap:wrap}
.gene-name{font-size:24px;font-weight:800;color:#1e3a5f}
.gene-locus{font-size:13px;color:#667}
.entrez-link{font-size:11px;color:#1e6fad;text-decoration:none;border:1px solid #1e6fad;
  padding:2px 8px;border-radius:3px;margin-left:auto;white-space:nowrap}
.entrez-link:hover{background:#1e6fad;color:#fff}
.stats-row{display:flex;gap:9px;margin-bottom:16px;flex-wrap:wrap}
.stat{background:#f0f4f8;border-radius:8px;padding:7px 13px;text-align:center;min-width:74px}
.stat-val{display:block;font-size:18px;font-weight:700;color:#1e3a5f}
.stat-lbl{font-size:9px;color:#789;text-transform:uppercase;letter-spacing:.4px}
.chart-section-title{font-size:11px;font-weight:700;color:#445;margin-bottom:5px;
  letter-spacing:.3px;text-transform:uppercase}
.chart-wrap{background:#fafbfc;border:1px solid #e5eaf0;border-radius:7px;padding:8px 8px 4px}
.cv-legend{display:flex;gap:14px;font-size:10px;color:#667;margin-top:7px;flex-wrap:wrap}
.cv-sw{width:9px;height:9px;border-radius:50%;display:inline-block;margin-right:3px;vertical-align:middle}
#legend-bar{background:#f5f7fa;border-top:1px solid #dde2ea;padding:5px 18px;
  display:flex;align-items:center;gap:8px;font-size:11px;color:#667;flex-shrink:0;flex-wrap:wrap}
#legend-label{font-weight:700;color:#1e3a5f}
#tooltip{position:fixed;background:rgba(15,20,40,.93);color:#fff;padding:8px 11px;
  border-radius:6px;font-size:12px;line-height:1.55;pointer-events:none;display:none;
  z-index:9999;max-width:230px;box-shadow:0 3px 12px rgba(0,0,0,.3)}
.track-bg{fill:#f8f9fb}
</style>
</head>
<body>
<div id="app">
<header>
  <div class="title-row">
    <h1>Cancer Protein IDR Atlas</h1>
    <p>50 CGC cancer genes &bull; GRCh38 loci &bull; IDR / LLPS / Virus / VIPP / ClinVar</p>
  </div>
  <div class="controls">
    <label>Color by:
      <select id="color-by">
        <option value="mean_iupred">IDR disorder</option>
        <option value="mean_llps">LLPS propensity</option>
        <option value="mean_conservation">Conservation</option>
        <option value="virus_interaction_fraction">Virus interaction fraction</option>
        <option value="mean_vipp_score">Mean VIPP score</option>
        <option value="low_complexity_fraction">Low complexity</option>
        <option value="mean_structural">Structural proxy</option>
      </select>
    </label>
    <input id="gene-search" type="search" placeholder="Search gene…" autocomplete="off">
  </div>
</header>
<div id="main">
  <div id="karyo-panel"><svg id="karyo-svg"></svg></div>
  <div id="detail-panel">
    <div id="placeholder">
      <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="#aab" stroke-width="1.5">
        <circle cx="12" cy="12" r="10"/><path d="M12 8v4m0 4h.01"/>
      </svg>
      Click a gene marker on the karyotype to view its full IDR Atlas profile
    </div>
    <div id="gene-detail" style="display:none"></div>
  </div>
</div>
<div id="legend-bar">
  Low&nbsp;
  <svg id="legend-svg" style="vertical-align:middle"></svg>
  &nbsp;High &nbsp;&bull;&nbsp;
  <span id="legend-label"></span>
  &nbsp;&bull;&nbsp;
  <span><span class="cv-sw" style="background:#e34a33"></span>Pathogenic</span>
  <span><span class="cv-sw" style="background:#1a9850"></span>Benign</span>
  <span><span class="cv-sw" style="background:#f4a636"></span>VUS</span>
  <span><span class="cv-sw" style="background:#888"></span>Other</span>
</div>
</div>
<div id="tooltip"></div>

<script>
// ── Embedded atlas data ───────────────────────────────────────────────────────
const CHROM_SIZES  = __CHROM_SIZES__;
const CENTROMERES  = __CENTROMERES__;
const GENE_DATA    = __GENE_DATA__;
const RESIDUE_DATA = __RESIDUE_DATA__;

// ── Layout ────────────────────────────────────────────────────────────────────
const CHROM_ROWS = [
  ['1','2','3','4','5','6'],
  ['7','8','9','10','11','12'],
  ['13','14','15','16','17','18'],
  ['19','20','21','22','X','Y'],
];
const K = {
  maxH: 210, cw: 12, colW: 80, rowH: 268,
  cH: 5,     r: 5.5, spacing: 13,
  mt: 10, ml: 10,
};

// ── Color scale (blue=ordered, red=disordered) ────────────────────────────────
let colorMetric = 'mean_iupred';
const colorScale = d3.scaleSequential()
  .domain([0, 1])
  .interpolator(t => d3.interpolateRdBu(1 - t))
  .clamp(true);

function geneColor(gene) {
  const v = (GENE_DATA[gene] || {})[colorMetric] ?? 0.5;
  return colorScale(v);
}

// ── Tooltip ───────────────────────────────────────────────────────────────────
const tip = d3.select('#tooltip');
function showTip(evt, gene) {
  const d = GENE_DATA[gene] || {};
  tip.style('display','block')
    .html(`<strong style="font-size:13px">${gene}</strong><br>` +
          `Chr ${d.chrom} &bull; ${d.band}<br>` +
          `${d.length} aa &bull; ${d.variant_count} ClinVar variants<br>` +
          `IDR: <b>${(d.mean_iupred||0).toFixed(3)}</b> &nbsp;` +
          `LLPS: <b>${(d.mean_llps||0).toFixed(3)}</b> &nbsp;` +
          `VIPP: <b>${(d.mean_vipp_score||0).toFixed(3)}</b><br>` +
          `Virus frac: <b>${((d.virus_interaction_fraction||0)*100).toFixed(1)}%</b> &nbsp;` +
          `Conservation: <b>${(d.mean_conservation||0).toFixed(3)}</b>`)
    .style('left',(evt.clientX+14)+'px')
    .style('top', (evt.clientY-10)+'px');
}
function moveTip(evt) {
  tip.style('left',(evt.clientX+14)+'px').style('top',(evt.clientY-10)+'px');
}
function hideTip() { tip.style('display','none'); }

// ── Chromosome drawing helpers ────────────────────────────────────────────────
function armPath(w, y0, h, roundTop, roundBot) {
  const r = Math.min(w / 2, h, 6);
  const rt = roundTop ? r : 0, rb = roundBot ? r : 0;
  return `M ${rt},${y0}` +
    (roundTop ? ` Q 0,${y0} 0,${y0+rt}` : ` L 0,${y0}`) +
    ` L 0,${y0+h-rb}` +
    (roundBot ? ` Q 0,${y0+h} ${rb},${y0+h}` : ` L 0,${y0+h}`) +
    ` L ${w-rb},${y0+h}` +
    (roundBot ? ` Q ${w},${y0+h} ${w},${y0+h-rb}` : ` L ${w},${y0+h}`) +
    ` L ${w},${y0+rt}` +
    (roundTop ? ` Q ${w},${y0} ${w-rt},${y0}` : ` L ${w},${y0}`) +
    ' Z';
}

// ── Marker layout (stagger overlapping genes horizontally) ────────────────────
function layoutMarkers(genes, chromLen, chromH) {
  if (!genes.length) return [];
  const minYSep = K.r * 2.4;
  const yPos = genes.map(g => (((g.start + g.end) / 2) / chromLen) * chromH);
  const offsets = new Array(genes.length).fill(0);

  let i = 0;
  while (i < genes.length) {
    let j = i + 1;
    while (j < genes.length && yPos[j] - yPos[i] < minYSep) j++;
    const grpSize = j - i;
    for (let k = 0; k < grpSize; k++) {
      offsets[i + k] = (k - (grpSize - 1) / 2) * K.spacing;
    }
    i = j;
  }

  return genes.map((g, idx) => ({
    gene: g.gene,
    x: K.cw / 2 + offsets[idx],
    y: yPos[idx],
  }));
}

// ── Draw karyotype ────────────────────────────────────────────────────────────
let selectedGene = null;

function drawKaryotype() {
  const maxLen = CHROM_SIZES['1'];
  const svgW = 6 * K.colW + K.ml * 2;
  const svgH = 4 * K.rowH + K.mt * 2;

  const svg = d3.select('#karyo-svg')
    .attr('width', svgW).attr('height', svgH);
  svg.selectAll('*').remove();

  CHROM_ROWS.forEach((row, ri) => {
    row.forEach((chrom, ci) => {
      const cLen = CHROM_SIZES[chrom] || 1;
      const cMid = CENTROMERES[chrom] || cLen * 0.45;
      const cH   = (cLen / maxLen) * K.maxH;
      const cY   = (cMid / cLen) * cH;

      const ox = K.ml + ci * K.colW + (K.colW - K.cw) / 2;
      const oy = K.mt + ri * K.rowH + (K.maxH - cH) / 2 + 14;

      const g = svg.append('g').attr('transform', `translate(${ox},${oy})`);

      // p arm
      const pH = cY - K.cH;
      if (pH > 1)
        g.append('path').attr('d', armPath(K.cw, 0, pH, true, false)).attr('class','chrom-arm');

      // q arm
      const qY = cY + K.cH;
      const qH = cH - qY;
      if (qH > 1)
        g.append('path').attr('d', armPath(K.cw, qY, qH, false, true)).attr('class','chrom-arm');

      // centromere
      g.append('ellipse')
        .attr('cx', K.cw / 2).attr('cy', cY)
        .attr('rx', K.cw / 2).attr('ry', K.cH + 1)
        .attr('class','centromere');

      // chromosome label
      g.append('text')
        .attr('x', K.cw / 2).attr('y', cH + 13)
        .attr('text-anchor','middle').attr('class','chrom-label')
        .text(chrom);

      // gene markers
      const chromGenes = Object.entries(GENE_DATA)
        .filter(([, d]) => d.chrom === chrom)
        .map(([gene, d]) => ({ gene, start: d.start, end: d.end }))
        .sort((a, b) => a.start - b.start);

      const positions = layoutMarkers(chromGenes, cLen, cH);

      positions.forEach(({ gene, x, y }) => {
        const mg = g.append('g')
          .attr('class','gene-marker')
          .attr('data-gene', gene);

        if (Math.abs(x - K.cw / 2) > 2)
          mg.append('line')
            .attr('x1', K.cw / 2).attr('y1', y)
            .attr('x2', x).attr('y2', y)
            .attr('stroke','#aab').attr('stroke-width',.7);

        mg.append('circle')
          .attr('cx', x).attr('cy', y).attr('r', K.r)
          .attr('fill', geneColor(gene))
          .attr('stroke', gene === selectedGene ? '#ffcc00' : '#fff')
          .attr('stroke-width', gene === selectedGene ? 2.5 : 1.2);

        if (chromGenes.length <= 2)
          mg.append('text')
            .attr('x', x + (x >= K.cw / 2 ? K.r + 2 : -(K.r + 2)))
            .attr('y', y + 3)
            .attr('text-anchor', x >= K.cw / 2 ? 'start' : 'end')
            .attr('class','gene-lbl-sm').text(gene);

        mg.on('mouseover', evt => showTip(evt, gene))
          .on('mousemove', moveTip)
          .on('mouseout', hideTip)
          .on('click', () => { selectedGene = gene; refreshMarkers(); showGeneDetail(gene); });
      });
    });
  });
}

function refreshMarkers() {
  d3.selectAll('.gene-marker').each(function() {
    const g = d3.select(this).attr('data-gene');
    d3.select(this).select('circle')
      .attr('stroke',       g === selectedGene ? '#ffcc00' : '#fff')
      .attr('stroke-width', g === selectedGene ? 2.5 : 1.2);
  });
}

function updateColors() {
  d3.selectAll('.gene-marker').each(function() {
    const g = d3.select(this).attr('data-gene');
    d3.select(this).select('circle').attr('fill', geneColor(g));
  });
  updateLegend();
}

// ── Gene detail panel ─────────────────────────────────────────────────────────
function showGeneDetail(gene) {
  const d = GENE_DATA[gene];
  const rdata = RESIDUE_DATA[gene];
  if (!d) return;

  d3.select('#placeholder').style('display','none');
  const panel = d3.select('#gene-detail').style('display','block');
  const ncbiUrl = `https://www.ncbi.nlm.nih.gov/gene/?term=${gene}%5BGene+Name%5D+AND+9606%5BTaxon%5D`;

  panel.html(
    `<div class="gene-header">
       <span class="gene-name">${gene}</span>
       <span class="gene-locus">chr${d.chrom} &bull; ${d.band} &bull;
         ${d.start.toLocaleString()}–${d.end.toLocaleString()} bp</span>
       <a class="entrez-link" href="${ncbiUrl}" target="_blank" rel="noopener noreferrer">
         NCBI Entrez Gene &#8599;
       </a>
     </div>
     <div class="stats-row">
       <div class="stat"><span class="stat-val">${d.length}</span><span class="stat-lbl">Residues</span></div>
       <div class="stat"><span class="stat-val">${d.mean_iupred.toFixed(3)}</span><span class="stat-lbl">Mean IDR</span></div>
       <div class="stat"><span class="stat-val">${d.mean_llps.toFixed(3)}</span><span class="stat-lbl">LLPS</span></div>
       <div class="stat"><span class="stat-val">${d.mean_conservation.toFixed(3)}</span><span class="stat-lbl">Conservation</span></div>
      <div class="stat"><span class="stat-val">${d.mean_vipp_score.toFixed(3)}</span><span class="stat-lbl">Mean VIPP</span></div>
      <div class="stat"><span class="stat-val">${(d.virus_interaction_fraction*100).toFixed(1)}%</span><span class="stat-lbl">Virus-linked</span></div>
       <div class="stat"><span class="stat-val">${(d.low_complexity_fraction*100).toFixed(1)}%</span><span class="stat-lbl">Low Complexity</span></div>
       <div class="stat"><span class="stat-val">${d.variant_count}</span><span class="stat-lbl">ClinVar vars</span></div>
     </div>
     <div class="chart-section-title">Per-residue atlas tracks</div>
     <div class="chart-wrap"><div id="residue-chart"></div></div>
     <div class="cv-legend">
       <span><span class="cv-sw" style="background:#e34a33"></span>Pathogenic</span>
       <span><span class="cv-sw" style="background:#1a9850"></span>Benign</span>
       <span><span class="cv-sw" style="background:#f4a636"></span>VUS / conflicting</span>
       <span><span class="cv-sw" style="background:#888"></span>Other / unknown</span>
     </div>`
  );

  if (rdata) renderResidueChart(rdata);
}

// ── Per-residue multi-track chart ─────────────────────────────────────────────
function sigColor(sig) {
  if (!sig) return '#bbb';
  const s = sig.toLowerCase();
  if (s.includes('pathogenic') && !s.includes('likely_benign') && !s.includes('benign'))
    return '#e34a33';
  if (s.includes('benign'))  return '#1a9850';
  if (s.includes('uncertain') || s.includes('vus') || s.includes('conflicting'))
    return '#f4a636';
  return '#888';
}

function renderResidueChart(rdata) {
  const el = document.getElementById('residue-chart');
  if (!el) return;

  const W  = el.clientWidth || 520;
  const m  = { top: 4, right: 8, bottom: 24, left: 42 };
  const tH = 54, tG = 5;
  const hasCv = !!(rdata.cv_count);
  const cvH   = hasCv ? 24 : 0;

  const IDR_THRESHOLD = 0.5;

  const hasPlddt = rdata.plddt && rdata.plddt.some(v => v != null);

  const tracks = [
    { key: 'iupred',       label: 'IDR score (IUPred2A)', color: '#d73027', fill: true, threshold: IDR_THRESHOLD },
    ...(hasPlddt ? [{ key: 'plddt', label: 'AlphaFold pLDDT (structure confidence)', color: '#2166ac', fill: false, domain: [0,100], threshold70: true }] : []),
    { key: 'llps',         label: 'LLPS propensity',      color: '#7b2d8b', fill: true  },
    { key: 'vipp',         label: 'VIPP composite score', color: '#f39c12', fill: false },
    { key: 'virus_interaction', label: 'Virus interaction flag', color: '#9b59b6', fill: false },
    { key: 'conservation', label: 'Conservation',         color: '#1a9850', fill: false },
    { key: 'structural',   label: 'Structural proxy',     color: '#2166ac', fill: false },
  ];

  const totalH = m.top + tracks.length * (tH + tG) + cvH + m.bottom;
  const n = rdata.pos.length;

  const svg = d3.select('#residue-chart')
    .append('svg').attr('width', W).attr('height', totalH);

  const xRange = [m.left, W - m.right];
  const xScale = d3.scaleLinear().domain([1, n]).range(xRange);

  // shared x-axis
  svg.append('g')
    .attr('transform', `translate(0,${totalH - m.bottom + 2})`)
    .call(d3.axisBottom(xScale).ticks(Math.min(n, 10)).tickFormat(d3.format('d')))
    .call(ag => ag.select('.domain').attr('stroke','#ccc'))
    .call(ag => ag.selectAll('line').attr('stroke','#ccc'))
    .call(ag => ag.selectAll('text').attr('fill','#999').attr('font-size','9px'));

  svg.append('text')
    .attr('x', (xRange[0]+xRange[1])/2).attr('y', totalH-2)
    .attr('text-anchor','middle').attr('fill','#bbb').attr('font-size','9px')
    .text('Amino acid position');

  tracks.forEach((track, ti) => {
    const y0    = m.top + ti * (tH + tG);
    const dom   = track.domain || [0,1];
    const yScale = d3.scaleLinear().domain(dom).range([y0+tH, y0]);
    const vals   = rdata[track.key];
    if (!vals) return;

    const tg = svg.append('g');

    tg.append('rect')
      .attr('x', xRange[0]).attr('y', y0)
      .attr('width', xRange[1]-xRange[0]).attr('height', tH)
      .attr('class','track-bg');

    const fmt = dom[1] > 1 ? d3.format('d') : d3.format('.1f');
    tg.append('g')
      .attr('transform', `translate(${m.left},0)`)
      .call(d3.axisLeft(yScale).ticks(3).tickFormat(fmt))
      .call(ag => ag.select('.domain').attr('stroke','#ddd'))
      .call(ag => ag.selectAll('line').attr('stroke','#eee'))
      .call(ag => ag.selectAll('text').attr('fill','#bbb').attr('font-size','8px'));

    tg.append('text')
      .attr('x', xRange[0]+4).attr('y', y0+11)
      .attr('fill', track.color).attr('font-size','9px').attr('font-weight','700')
      .text(track.label);

    if (track.fill) {
      const thresh = track.threshold || 0;
      const clipped = vals.map(v => Math.max(0, Math.min(1, v) - thresh));
      tg.append('path')
        .datum(clipped)
        .attr('d', d3.area()
          .x((d,i) => xScale(rdata.pos[i]))
          .y0(yScale(thresh))
          .y1((d,i) => yScale(Math.min(1, vals[i])))
          .defined((d) => d > 0)
          .curve(d3.curveLinear))
        .attr('fill', track.color).attr('opacity', .28);
    }

    if (track.threshold != null) {
      tg.append('line')
        .attr('x1', xRange[0]).attr('x2', xRange[1])
        .attr('y1', yScale(track.threshold)).attr('y2', yScale(track.threshold))
        .attr('stroke', '#999').attr('stroke-width', 0.8)
        .attr('stroke-dasharray', '4,3');
      tg.append('text')
        .attr('x', xRange[1] - 2).attr('y', yScale(track.threshold) - 3)
        .attr('text-anchor', 'end').attr('fill', '#999').attr('font-size', '8px')
        .text('IDR threshold (0.5)');
    }

    if (track.threshold70) {
      tg.append('line')
        .attr('x1', xRange[0]).attr('x2', xRange[1])
        .attr('y1', yScale(70)).attr('y2', yScale(70))
        .attr('stroke', '#e67e22').attr('stroke-width', 0.8)
        .attr('stroke-dasharray', '4,3');
      tg.append('text')
        .attr('x', xRange[1] - 2).attr('y', yScale(70) - 3)
        .attr('text-anchor', 'end').attr('fill', '#e67e22').attr('font-size', '8px')
        .text('Confident structure (70)');
    }

    const validVals = vals.map(v => v == null ? dom[0] : v);
    tg.append('path')
      .datum(validVals)
      .attr('d', d3.line()
        .x((d,i) => xScale(rdata.pos[i]))
        .y(d => yScale(Math.max(dom[0], Math.min(dom[1], d))))
        .defined((d,i) => vals[i] != null)
        .curve(d3.curveLinear))
      .attr('fill','none')
      .attr('stroke', track.color)
      .attr('stroke-width', track.fill ? 1.2 : 1.5);

    // Mark conditionally disordered residues on IDR track
    if (track.key === 'iupred' && rdata.idr_class) {
      rdata.idr_class.forEach((cls, i) => {
        if (cls === 'conditionally_disordered') {
          tg.append('circle')
            .attr('cx', xScale(rdata.pos[i]))
            .attr('cy', yScale(Math.min(1, vals[i])))
            .attr('r', 2.5)
            .attr('fill', '#e67e22')
            .attr('stroke', '#fff')
            .attr('stroke-width', 0.5)
            .attr('opacity', 0.85)
            .append('title')
            .text(`${rdata.aa[i]}${rdata.pos[i]}: conditionally disordered (IUPred=${vals[i].toFixed(2)}, pLDDT≥70)`);
        }
      });
    }
  });

  // ClinVar tick track
  if (hasCv) {
    const cvY0   = m.top + tracks.length * (tH + tG);
    const cvYBot = cvY0 + cvH - 4;

    const cvG = svg.append('g');
    cvG.append('text')
      .attr('x', xRange[0]+4).attr('y', cvY0+11)
      .attr('fill','#888').attr('font-size','9px').attr('font-weight','700')
      .text('ClinVar variants');

    rdata.cv_count.forEach((cnt, i) => {
      if (!cnt) return;
      const sig = rdata.cv_sig?.[i] || '';
      const xv  = xScale(rdata.pos[i]);
      cvG.append('line')
        .attr('x1',xv).attr('y1',cvY0+13)
        .attr('x2',xv).attr('y2',cvYBot)
        .attr('stroke', sigColor(sig))
        .attr('stroke-width', Math.min(5, 1+cnt*.6))
        .attr('opacity',.85)
        .append('title')
        .text(`${rdata.aa?.[i]||''}${rdata.pos[i]}: ${sig||'variant'} (n=${cnt})`);
    });
  }
}

// ── Legend ────────────────────────────────────────────────────────────────────
const METRIC_LABELS = {
  mean_iupred:             'IDR disorder score',
  mean_llps:               'LLPS propensity',
  mean_conservation:       'Phylogenetic conservation',
  virus_interaction_fraction: 'Virus interaction fraction',
  mean_vipp_score:         'Mean VIPP score',
  low_complexity_fraction: 'Low-complexity fraction',
  mean_structural:         'Structural proxy',
};

function updateLegend() {
  const W = 180, H = 13;
  const svg = d3.select('#legend-svg').attr('width',W).attr('height',H);
  svg.selectAll('*').remove();
  const gid = 'lg'+Date.now();
  const grad = svg.append('defs').append('linearGradient')
    .attr('id',gid).attr('x1','0%').attr('x2','100%');
  d3.range(11).forEach(i =>
    grad.append('stop').attr('offset',`${i*10}%`).attr('stop-color', colorScale(i/10)));
  svg.append('rect').attr('width',W).attr('height',H).attr('rx',3).attr('fill',`url(#${gid})`);
  d3.select('#legend-label').text(METRIC_LABELS[colorMetric] || colorMetric);
}

// ── Search ────────────────────────────────────────────────────────────────────
document.getElementById('gene-search').addEventListener('input', function() {
  const q = this.value.trim().toUpperCase();
  d3.selectAll('.gene-marker').style('opacity', function() {
    const g = d3.select(this).attr('data-gene') || '';
    return (!q || g.includes(q)) ? 1 : 0.1;
  });
});

// ── Color-by selector ─────────────────────────────────────────────────────────
document.getElementById('color-by').addEventListener('change', function() {
  colorMetric = this.value;
  updateColors();
});

// ── Init ──────────────────────────────────────────────────────────────────────
drawKaryotype();
updateLegend();
</script>
</body>
</html>
"""


def generate_html(gene_data: dict, residue_data: dict) -> str:
    return (
        _TEMPLATE
        .replace("__CHROM_SIZES__",  json.dumps(CHROM_SIZES,  separators=(",", ":")))
        .replace("__CENTROMERES__",  json.dumps(CENTROMERES,  separators=(",", ":")))
        .replace("__GENE_DATA__",    json.dumps(gene_data,    separators=(",", ":")))
        .replace("__RESIDUE_DATA__", json.dumps(residue_data, separators=(",", ":")))
    )


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Render the Cancer Protein IDR Atlas as an interactive chromosome map."
    )
    parser.add_argument("--atlas-dir",      required=True, help="Directory with *_atlas_with_clinvar.tsv files")
    parser.add_argument("--global-summary", required=True, help="global_disorder_phylogeny_atlas.tsv")
    parser.add_argument("--output",         required=True, help="Output HTML path")
    args = parser.parse_args()

    gene_data, residue_data = load_data(Path(args.atlas_dir), Path(args.global_summary))
    html = generate_html(gene_data, residue_data)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    print(f"Wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
