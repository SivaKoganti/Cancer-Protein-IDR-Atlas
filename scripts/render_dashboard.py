#!/usr/bin/env python3
"""render_dashboard.py

Builds results/index.html — the main landing page dashboard for the Cancer Protein
IDR Atlas with navigation buttons, gene search, summary statistics, and quick-access
panels for interactive exploration and data export.
"""

import argparse
import json
from pathlib import Path

import pandas as pd


def load_data(atlas_dir: Path, summary_path: Path):
    """Load and aggregate summary statistics."""
    summary = pd.read_csv(summary_path, sep="\t")
    
    gene_list = []
    for _, row in summary.iterrows():
        gene_list.append({
            "gene": str(row["gene"]),
            "length": int(row["length"]),
            "mean_iupred": round(float(row["mean_iupred"]), 4),
            "mean_llps": round(float(row["mean_llps"]), 4),
            "mean_conservation": round(float(row["mean_conservation"]), 4),
            "mean_structural": round(float(row["mean_structural"]), 4),
            "low_complexity_fraction": round(float(row["low_complexity_fraction"]), 4),
            "variant_count": int(row["variant_count"]),
          "virus_interaction_fraction": round(float(row.get("virus_interaction_fraction", 0.0)), 4),
          "mean_vipp_score": round(float(row.get("mean_vipp_score", 0.0)), 4),
        })
    
    gene_list = sorted(gene_list, key=lambda x: x["gene"])
    
    # Global stats
    stats = {
        "total_genes": len(gene_list),
        "avg_length": round(sum(g["length"] for g in gene_list) / len(gene_list), 0),
        "avg_iupred": round(sum(g["mean_iupred"] for g in gene_list) / len(gene_list), 4),
        "avg_llps": round(sum(g["mean_llps"] for g in gene_list) / len(gene_list), 4),
        "avg_conservation": round(sum(g["mean_conservation"] for g in gene_list) / len(gene_list), 4),
        "avg_vipp": round(sum(g["mean_vipp_score"] for g in gene_list) / len(gene_list), 4),
        "avg_virus_fraction": round(sum(g["virus_interaction_fraction"] for g in gene_list) / len(gene_list), 4),
        "total_variants": sum(g["variant_count"] for g in gene_list),
    }
    
    return gene_list, stats


# ── HTML Template ─────────────────────────────────────────────────────────────

_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Cancer Protein IDR Atlas — Dashboard</title>
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;background:#f5f7fa;
  color:#1a1a1a;line-height:1.6}
.navbar{background:linear-gradient(135deg,#1e3a5f 0%,#2d5a8e 100%);color:#fff;
  padding:14px 20px;display:flex;align-items:center;justify-content:space-between;
  box-shadow:0 2px 8px rgba(0,0,0,.15)}
.navbar h1{font-size:20px;font-weight:800;letter-spacing:.5px}
.navbar-subtitle{font-size:12px;opacity:.8;margin-top:2px}
.navbar-spacer{margin-left:auto}
.navbar-date{font-size:11px;opacity:.72;display:flex;gap:2px}
main{padding:20px 24px;max-width:1400px;margin:0 auto}
.welcome{background:#fff;border-radius:10px;padding:24px;margin-bottom:28px;
  box-shadow:0 2px 6px rgba(0,0,0,.08);border-left:5px solid #1e3a5f}
.welcome h2{font-size:22px;margin-bottom:8px;color:#1e3a5f}
.welcome p{font-size:14px;color:#666;line-height:1.7;max-width:800px}
.stats-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));
  gap:14px;margin-bottom:28px}
.stat-card{background:#fff;border-radius:8px;padding:16px;
  box-shadow:0 2px 6px rgba(0,0,0,.08);border-top:3px solid #2d5a8e;text-align:center}
.stat-card.alt1{border-top-color:#d73027}
.stat-card.alt2{border-top-color:#7b2d8b}
.stat-card.alt3{border-top-color:#1a9850}
.stat-num{font-size:28px;font-weight:800;color:#1e3a5f;display:block;margin-bottom:4px}
.stat-label{font-size:11px;color:#999;text-transform:uppercase;letter-spacing:.3px}
.section-title{font-size:16px;font-weight:800;color:#1e3a5f;margin-bottom:14px;
  margin-top:28px;padding-bottom:8px;border-bottom:2px solid #e5eaf0}
.button-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));
  gap:14px;margin-bottom:28px}
.btn{display:flex;align-items:center;justify-content:center;gap:10px;
  padding:16px 20px;border:none;border-radius:8px;cursor:pointer;
  font-size:14px;font-weight:600;transition:all .2s;text-decoration:none;
  color:#fff;text-align:center;box-shadow:0 2px 8px rgba(0,0,0,.12)}
.btn:hover{transform:translateY(-2px);box-shadow:0 4px 14px rgba(0,0,0,.18)}
.btn:active{transform:translateY(0)}
.btn-primary{background:linear-gradient(135deg,#1e3a5f 0%,#2d5a8e 100%)}
.btn-secondary{background:linear-gradient(135deg,#2d5a8e 0%,#4a80b8 100%)}
.btn-success{background:linear-gradient(135deg,#1a9850 0%,#22b870 100%)}
.btn-info{background:linear-gradient(135deg,#2166ac 0%,#3a8fd7 100%)}
.btn-danger{background:linear-gradient(135deg,#d73027 0%,#f1605e 100%)}
.btn-warning{background:linear-gradient(135deg,#f4a636 0%,#ffc857 100%)}
.btn svg{width:20px;height:20px;stroke:currentColor;fill:none;stroke-width:2}
.gene-table{width:100%;border-collapse:collapse;background:#fff;
  border-radius:8px;overflow:hidden;box-shadow:0 2px 6px rgba(0,0,0,.08);
  margin-bottom:20px}
.gene-table thead{background:#f0f4f8;border-bottom:2px solid #dde2ea}
.gene-table th{padding:12px 14px;text-align:left;font-weight:700;font-size:12px;
  color:#1e3a5f;text-transform:uppercase;letter-spacing:.3px}
.gene-table td{padding:10px 14px;border-bottom:1px solid #e5eaf0;font-size:12px}
.gene-table tr:hover{background:#fafbfc}
.gene-name{font-weight:700;color:#1e3a5f}
.gene-table-search{margin-bottom:14px;display:flex;gap:8px}
.gene-table-search input{padding:10px 12px;border:1px solid #dde2ea;border-radius:6px;
  font-size:13px;flex:1;max-width:300px}
.table-section{max-height:500px;overflow-y:auto;margin-bottom:20px}
.modal{display:none;position:fixed;top:0;left:0;width:100%;height:100%;
  background:rgba(0,0,0,.5);z-index:1000;align-items:center;justify-content:center}
.modal.active{display:flex}
.modal-content{background:#fff;border-radius:12px;padding:28px;
  max-width:600px;width:90%;max-height:85vh;overflow-y:auto;box-shadow:0 10px 40px rgba(0,0,0,.3)}
.modal-header{font-size:18px;font-weight:800;color:#1e3a5f;margin-bottom:16px;
  padding-bottom:12px;border-bottom:2px solid #e5eaf0}
.modal-close{position:absolute;top:12px;right:14px;background:none;border:none;
  font-size:24px;cursor:pointer;color:#999;padding:4px;line-height:1}
.modal-close:hover{color:#1e3a5f}
.info-row{display:flex;justify-content:space-between;padding:8px 0;
  border-bottom:1px solid #e5eaf0;font-size:13px}
.info-label{font-weight:600;color:#1e3a5f}
.info-value{color:#666}
.download-list{list-style:none;margin-top:12px}
.download-list li{padding:8px 0;border-bottom:1px solid #e5eaf0}
.download-list li:last-child{border-bottom:none}
.download-list a{color:#1e6fad;text-decoration:none;font-weight:600;
  display:flex;align-items:center;gap:6px}
.download-list a:hover{text-decoration:underline}
.csv-icon{display:inline-block;width:14px;height:14px;background:#1e3a5f;
  border-radius:2px;color:#fff;font-size:9px;font-weight:700;
  text-align:center;line-height:14px}
.search-help{font-size:11px;color:#999;margin-top:4px}
footer{background:#f0f4f8;border-top:1px solid #dde2ea;padding:14px 20px;
  text-align:center;font-size:11px;color:#999;margin-top:40px}
.feature-list{list-style:none;margin-top:10px}
.feature-list li{padding:4px 0;padding-left:24px;position:relative;font-size:13px;
  color:#666}
.feature-list li::before{content:"✓";position:absolute;left:0;color:#1a9850;
  font-weight:800;font-size:14px}
</style>
</head>
<body>
<div class="navbar">
  <div>
    <h1>🧬 Cancer Protein IDR Atlas</h1>
    <div class="navbar-subtitle">50 CGC genes • IDR • LLPS • Oncovirus • VIPP • ClinVar • GRCh38</div>
  </div>
  <div class="navbar-spacer"></div>
  <div class="navbar-date" id="nav-date"></div>
</div>

<main>
  <div class="welcome">
    <h2>Welcome to the Cancer Protein IDR Atlas</h2>
    <p>
      This integrated resource maps 50 Cancer Gene Census genes to their human chromosome loci
      (GRCh38) and annotates each residue with intrinsic disorder, LLPS propensity, phylogenetic
      conservation, structural properties, and clinical variants. Explore per-residue detail views,
      download annotated atlas tables, and link directly to NCBI Entrez Gene.
    </p>
    <ul class="feature-list">
      <li>Interactive chromosome karyotype with color-coded gene markers</li>
      <li>Per-residue IDR, LLPS, conservation, structural, virus-interaction, and VIPP scoring</li>
      <li>ClinVar integration with pathogenic/benign/VUS classification</li>
      <li>Download per-gene atlas tables and global summary statistics</li>
      <li>Direct links to NCBI Entrez Gene, UniProt, and ENSEMBL</li>
    </ul>
  </div>

  <div class="stats-grid">
    <div class="stat-card">
      <span class="stat-num">__TOTAL_GENES__</span>
      <span class="stat-label">Cancer Genes</span>
    </div>
    <div class="stat-card alt1">
      <span class="stat-num">__AVG_LENGTH__</span>
      <span class="stat-label">Mean Protein Length</span>
    </div>
    <div class="stat-card alt2">
      <span class="stat-num">__AVG_LLPS__</span>
      <span class="stat-label">Mean LLPS Propensity</span>
    </div>
    <div class="stat-card alt3">
      <span class="stat-num">__TOTAL_VARIANTS__</span>
      <span class="stat-label">ClinVar Variants</span>
    </div>
    <div class="stat-card">
      <span class="stat-num">__AVG_VIPP__</span>
      <span class="stat-label">Mean VIPP Score</span>
    </div>
    <div class="stat-card alt2">
      <span class="stat-num">__AVG_VIRUS_PCT__</span>
      <span class="stat-label">Mean Virus-linked Fraction</span>
    </div>
  </div>

  <h3 class="section-title">📊 Main Features</h3>
  <div class="button-grid">
    <a href="map/chromosome_atlas.html" class="btn btn-primary">
      <svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="10"/><path d="M7 12h10M12 7v10"/></svg>
      Interactive Karyotype
    </a>
    <button class="btn btn-secondary" onclick="openModal('genes')">
      <svg viewBox="0 0 24 24"><rect x="3" y="3" width="18" height="18" rx="2"/><path d="M3 9h18M9 3v18"/></svg>
      Gene Table
    </button>
    <button class="btn btn-info" onclick="openModal('download')">
      <svg viewBox="0 0 24 24"><path d="M7 18c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zM17 18c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zM16.75 7l1.35 4H9l-2-4M7 4h10v2H7z"/></svg>
      Download Data
    </button>
    <button class="btn btn-success" onclick="openModal('search')">
      <svg viewBox="0 0 24 24"><circle cx="10" cy="10" r="6"/><path d="m14.5 14.5 5 5"/></svg>
      Search Gene
    </button>
    <button class="btn btn-warning" onclick="openModal('stats')">
      <svg viewBox="0 0 24 24"><rect x="3" y="13" width="3" height="8"/><rect x="10" y="3" width="3" height="18"/><rect x="17" y="8" width="3" height="13"/></svg>
      Statistics
    </button>
    <button class="btn btn-danger" onclick="openModal('help')">
      <svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="10"/><text x="12" y="16" text-anchor="middle" font-size="14" font-weight="bold">?</text></svg>
      Help &amp; Docs
    </button>
  </div>

  <h3 class="section-title">🔍 Gene Summary Table</h3>
  <div class="gene-table-search">
    <input type="search" id="gene-filter" placeholder="Filter by gene name (e.g. TP53, BRCA1)…" autocomplete="off">
    <span class="search-help">Results shown: <span id="gene-count">__TOTAL_GENES__</span></span>
  </div>
  <div class="table-section">
    <table class="gene-table" id="gene-table">
      <thead>
        <tr>
          <th>Gene</th>
          <th>Length (aa)</th>
          <th>IDR Score</th>
          <th>LLPS</th>
          <th>Conservation</th>
          <th>VIPP</th>
          <th>Virus %</th>
          <th>ClinVar Vars</th>
        </tr>
      </thead>
      <tbody id="gene-tbody">
      </tbody>
    </table>
  </div>
</main>

<!-- Modals -->
<div id="modal-genes" class="modal">
  <div class="modal-content">
    <button class="modal-close" onclick="closeModal('genes')">×</button>
    <div class="modal-header">Gene Summary Table</div>
    <p>Sortable table of all 50 cancer genes with their key metrics. Click on a gene to view details.</p>
    <div style="max-height:400px;overflow-y:auto;margin-top:14px">
      <table class="gene-table" style="margin:0">
        <thead>
          <tr>
            <th>Gene</th>
            <th>Length</th>
            <th>IDR</th>
            <th>LLPS</th>
            <th>Conservation</th>
            <th>VIPP</th>
            <th>Virus %</th>
            <th>Variants</th>
          </tr>
        </thead>
        <tbody id="modal-gene-tbody"></tbody>
      </table>
    </div>
  </div>
</div>

<div id="modal-download" class="modal">
  <div class="modal-content">
    <button class="modal-close" onclick="closeModal('download')">×</button>
    <div class="modal-header">📥 Download Data</div>
    <p style="margin-bottom:14px;font-size:13px">
      Download per-gene atlas tables, global summaries, and metadata. Each file includes
      per-residue IDR scores, LLPS propensity, conservation, structural predictions, and
      ClinVar variant annotations.
    </p>
    <div style="background:#f0f4f8;border-radius:6px;padding:12px;margin-bottom:14px">
      <div class="info-row">
        <span class="info-label">Global Summary:</span>
        <a href="atlas/global_disorder_phylogeny_atlas.tsv" download style="color:#1e6fad;text-decoration:none;font-weight:600">
          ⬇ global_disorder_phylogeny_atlas.tsv
        </a>
      </div>
    </div>
    <div style="margin-bottom:14px">
      <strong style="font-size:13px">Per-Gene Atlas Tables (select to download):</strong>
      <ul class="download-list" id="download-list"></ul>
    </div>
    <div style="background:#eef1f5;border-radius:6px;padding:10px;font-size:12px;color:#666">
      <strong>Format:</strong> TSV (tab-separated values) with columns: gene, pos, aa, iupred_score, llps_proxy,
      conservation, virus_interaction, vipp_score, variants, clinvar_count, clinvar_significance
    </div>
  </div>
</div>

<div id="modal-search" class="modal">
  <div class="modal-content">
    <button class="modal-close" onclick="closeModal('search')">×</button>
    <div class="modal-header">🔍 Search Gene</div>
    <input type="search" id="search-input" placeholder="Enter gene symbol (e.g. TP53)…" 
           style="width:100%;padding:10px;border:1px solid #dde2ea;border-radius:6px;font-size:13px;margin-bottom:14px">
    <div id="search-results" style="margin-top:14px"></div>
  </div>
</div>

<div id="modal-stats" class="modal">
  <div class="modal-content">
    <button class="modal-close" onclick="closeModal('stats')">×</button>
    <div class="modal-header">📈 Global Statistics</div>
    <div style="display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:14px">
      <div style="background:#f0f4f8;border-radius:6px;padding:12px">
        <div class="stat-num">__TOTAL_GENES__</div>
        <div class="stat-label">Cancer Genes Analyzed</div>
      </div>
      <div style="background:#f0f4f8;border-radius:6px;padding:12px">
        <div class="stat-num">__TOTAL_VARIANTS__</div>
        <div class="stat-label">ClinVar Variants</div>
      </div>
      <div style="background:#f0f4f8;border-radius:6px;padding:12px">
        <div class="stat-num">__AVG_LENGTH__</div>
        <div class="stat-label">Mean Protein Length</div>
      </div>
      <div style="background:#f0f4f8;border-radius:6px;padding:12px">
        <div class="stat-num">__AVG_IUPRED__</div>
        <div class="stat-label">Mean IDR Score</div>
      </div>
      <div style="background:#f0f4f8;border-radius:6px;padding:12px">
        <div class="stat-num">__AVG_LLPS__</div>
        <div class="stat-label">Mean LLPS</div>
      </div>
      <div style="background:#f0f4f8;border-radius:6px;padding:12px">
        <div class="stat-num">__AVG_CONSERVATION__</div>
        <div class="stat-label">Mean Conservation</div>
      </div>
      <div style="background:#f0f4f8;border-radius:6px;padding:12px">
        <div class="stat-num">__AVG_VIPP__</div>
        <div class="stat-label">Mean VIPP Score</div>
      </div>
      <div style="background:#f0f4f8;border-radius:6px;padding:12px">
        <div class="stat-num">__AVG_VIRUS_PCT__</div>
        <div class="stat-label">Mean Virus-linked Fraction</div>
      </div>
    </div>
  </div>
</div>

<div id="modal-help" class="modal">
  <div class="modal-content">
    <button class="modal-close" onclick="closeModal('help')">×</button>
    <div class="modal-header">❓ Help &amp; Documentation</div>
    <div style="font-size:13px;line-height:1.8;color:#666">
      <h4 style="color:#1e3a5f;margin-top:16px;margin-bottom:8px">What is this resource?</h4>
      <p>
        The Cancer Protein IDR Atlas integrates intrinsic disorder (IDR), liquid-liquid phase separation (LLPS),
        phylogenetic conservation, and ClinVar variants across 50 Cancer Gene Census genes, all mapped to their
        GRCh38 genomic loci.
      </p>
      <h4 style="color:#1e3a5f;margin-top:16px;margin-bottom:8px">How to navigate:</h4>
      <ul style="margin-left:20px;margin-bottom:10px">
        <li><strong>Interactive Karyotype:</strong> Click a gene marker to open its full per-residue detail view with charts.</li>
        <li><strong>Gene Table:</strong> Search or filter the complete list of 50 genes by name.</li>
        <li><strong>Download Data:</strong> Export per-gene atlas TSVs for analysis in R, Python, or Excel.</li>
        <li><strong>Search Gene:</strong> Quick lookup to find a specific gene and jump to its profile.</li>
        <li><strong>Statistics:</strong> View global summary statistics across all genes.</li>
      </ul>
      <h4 style="color:#1e3a5f;margin-top:16px;margin-bottom:8px">Column definitions:</h4>
      <ul style="margin-left:20px">
        <li><strong>IDR Score:</strong> IUPred2A-predicted disorder probability (0–1; higher = more disordered).</li>
        <li><strong>LLPS:</strong> Prion-like and physicochemical LLPS propensity proxy (0–1; higher = more prone to LLPS).</li>
        <li><strong>Conservation:</strong> Phylogenetic conservation across 6 model organisms (0–1; higher = more conserved).</li>
        <li><strong>Virus-linked %:</strong> Fraction of residues mapped to curated oncovirus interaction regions.</li>
        <li><strong>VIPP:</strong> Composite Virus + IDR + LLPS prioritization score (0–1, higher = higher composite risk).</li>
        <li><strong>ClinVar Vars:</strong> Number of clinical variants annotated in ClinVar for this gene.</li>
      </ul>
      <h4 style="color:#1e3a5f;margin-top:16px;margin-bottom:8px">External links:</h4>
      <p>Each gene detail page includes direct links to:</p>
      <ul style="margin-left:20px">
        <li><a href="https://www.ncbi.nlm.nih.gov/gene/" target="_blank" rel="noopener">NCBI Entrez Gene</a></li>
        <li><a href="https://www.uniprot.org/" target="_blank" rel="noopener">UniProt</a></li>
        <li><a href="https://www.ensembl.org/" target="_blank" rel="noopener">Ensembl</a></li>
        <li><a href="https://www.ncbi.nlm.nih.gov/clinvar/" target="_blank" rel="noopener">ClinVar</a></li>
      </ul>
    </div>
  </div>
</div>

<footer>
  Cancer Protein IDR Atlas &bull; GRCh38 &bull; Built __DATE__
</footer>

<script>
const GENES = __GENE_LIST__;

// Populate gene table
function populateGeneTable(genes = GENES) {
  const tbody = document.getElementById('gene-tbody');
  tbody.innerHTML = '';
  genes.forEach(g => {
    const row = document.createElement('tr');
    row.innerHTML = `
      <td class="gene-name">${g.gene}</td>
      <td>${g.length}</td>
      <td>${g.mean_iupred.toFixed(3)}</td>
      <td>${g.mean_llps.toFixed(3)}</td>
      <td>${g.mean_conservation.toFixed(3)}</td>
      <td>${g.mean_vipp_score.toFixed(3)}</td>
      <td>${(g.virus_interaction_fraction * 100).toFixed(1)}%</td>
      <td>${g.variant_count}</td>
    `;
    row.style.cursor = 'pointer';
    row.onclick = () => {
      window.location.href = `map/chromosome_atlas.html`;
    };
    tbody.appendChild(row);
  });
  document.getElementById('gene-count').textContent = genes.length;
}

// Filter genes
document.getElementById('gene-filter').addEventListener('input', e => {
  const q = e.target.value.toUpperCase();
  const filtered = GENES.filter(g => g.gene.includes(q));
  populateGeneTable(filtered);
});

// Modal functions
function openModal(id) {
  const modal = document.getElementById(`modal-${id}`);
  if (modal) modal.classList.add('active');
  if (id === 'genes') populateModalGeneTable();
  if (id === 'download') populateDownloadList();
  if (id === 'search') document.getElementById('search-input').focus();
}

function closeModal(id) {
  const modal = document.getElementById(`modal-${id}`);
  if (modal) modal.classList.remove('active');
}

// Close on escape or outside click
document.addEventListener('keydown', e => {
  if (e.key === 'Escape') {
    document.querySelectorAll('.modal.active').forEach(m => m.classList.remove('active'));
  }
});

document.querySelectorAll('.modal').forEach(m => {
  m.addEventListener('click', e => {
    if (e.target === m) m.classList.remove('active');
  });
});

// Gene table in modal
function populateModalGeneTable() {
  const tbody = document.getElementById('modal-gene-tbody');
  tbody.innerHTML = '';
  GENES.forEach(g => {
    const row = document.createElement('tr');
    row.innerHTML = `
      <td class="gene-name">${g.gene}</td>
      <td>${g.length}</td>
      <td>${g.mean_iupred.toFixed(3)}</td>
      <td>${g.mean_llps.toFixed(3)}</td>
      <td>${g.mean_conservation.toFixed(3)}</td>
      <td>${g.mean_vipp_score.toFixed(3)}</td>
      <td>${(g.virus_interaction_fraction * 100).toFixed(1)}%</td>
      <td>${g.variant_count}</td>
    `;
    tbody.appendChild(row);
  });
}

// Download list
function populateDownloadList() {
  const list = document.getElementById('download-list');
  list.innerHTML = '';
  GENES.forEach(g => {
    const li = document.createElement('li');
    li.innerHTML = `
      <a href="atlas/${g.gene}_atlas_with_clinvar.tsv" download>
        <span class="csv-icon">TSV</span> ${g.gene}_atlas_with_clinvar.tsv
      </a>
    `;
    list.appendChild(li);
  });
}

// Search
document.getElementById('search-input')?.addEventListener('input', e => {
  const q = e.target.value.toUpperCase().trim();
  const resultsDiv = document.getElementById('search-results');
  if (!q) {
    resultsDiv.innerHTML = '<p style="color:#999">Type a gene name to search…</p>';
    return;
  }
  const matches = GENES.filter(g => g.gene.includes(q));
  if (matches.length === 0) {
    resultsDiv.innerHTML = '<p style="color:#d73027">No genes found.</p>';
  } else {
    resultsDiv.innerHTML = matches.map(g => `
      <div style="padding:10px;border:1px solid #dde2ea;border-radius:6px;margin-bottom:8px;cursor:pointer;hover:background:#f0f4f8"
           onclick="location.href='map/chromosome_atlas.html'">
        <strong style="font-size:14px">${g.gene}</strong><br>
        <span style="font-size:12px;color:#666">${g.length} aa • IDR: ${g.mean_iupred.toFixed(3)} • LLPS: ${g.mean_llps.toFixed(3)} • VIPP: ${g.mean_vipp_score.toFixed(3)} • Virus: ${(g.virus_interaction_fraction*100).toFixed(1)}% • ${g.variant_count} variants</span>
      </div>
    `).join('');
  }
});

// Init
populateGeneTable();
document.getElementById('nav-date').textContent = new Date().toLocaleDateString('en-US', {
  weekday: 'short', year: 'numeric', month: 'short', day: 'numeric'
});
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description="Render the dashboard landing page.")
    parser.add_argument("--atlas-dir",      required=True)
    parser.add_argument("--global-summary", required=True)
    parser.add_argument("--output",         required=True)
    args = parser.parse_args()

    gene_list, stats = load_data(Path(args.atlas_dir), Path(args.global_summary))

    html = (
        _TEMPLATE
        .replace("__TOTAL_GENES__",    str(stats["total_genes"]))
        .replace("__AVG_LENGTH__",     f"{stats['avg_length']:.0f}")
        .replace("__AVG_IUPRED__",     f"{stats['avg_iupred']:.3f}")
        .replace("__AVG_LLPS__",       f"{stats['avg_llps']:.3f}")
        .replace("__AVG_CONSERVATION__", f"{stats['avg_conservation']:.3f}")
        .replace("__AVG_VIPP__",       f"{stats['avg_vipp']:.3f}")
        .replace("__AVG_VIRUS_PCT__",  f"{stats['avg_virus_fraction'] * 100:.1f}%")
        .replace("__TOTAL_VARIANTS__", str(stats["total_variants"]))
        .replace("__GENE_LIST__",      json.dumps(gene_list, separators=(",", ":")))
        .replace("__DATE__",           Path(args.global_summary).stat().st_mtime_ns and "August 2026" or "")
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    print(f"Wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
