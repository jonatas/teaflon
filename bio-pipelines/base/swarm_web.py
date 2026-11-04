import os
import json
import time
import urllib.request
from flask import Flask, render_template_string, jsonify, request
import psycopg

app = Flask(__name__)
DB_URI = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MODELS_DIR = "/Users/jonatas/code/ideia.me/assets/models"
STRUCTURES_DIR = os.path.join(REPO_ROOT, "structures")

PIPELINE_PRESETS = {
    "ptfe-to-bioceramic": {
        "title": "PTFE / PFAS -> Fluorapatite Bioceramic (37°C, pH 7.4)",
        "temp_c": 37.0,
        "ph": 7.4,
        "smiles": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
        "ref_smiles": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",
        "sub_label": "Surface C-F Sites (µM)",
        "int_label": "Free Fluoride [F⁻] (µM)",
        "prod_label": "Fluorapatite Ca₅(PO₄)₃F (µM eq)",
        "kcat1": 24.0,
        "km1": 380.0,
        "kcat2": 42.0,
        "km2": 30.0,
        "ki": 22.0,
        "boost": 8.5,
        "f1": 0.50,
        "f2": 0.30,
        "f3": 0.20,
    },
    "plastics-upcycling": {
        "title": "PET & Ester-Urethane -> PCA / Bio-Monomer Sink (68°C, pH 8.0)",
        "temp_c": 68.0,
        "ph": 8.0,
        "smiles": "OC(=O)c1ccc(C(=O)OCCO)cc1",
        "ref_smiles": "OC(=O)c1ccc(C(=O)O)cc1",
        "sub_label": "PET Ester Bonds (µM)",
        "int_label": "MHET / Acidic Monomer (µM)",
        "prod_label": "Upcycled Sink Product (µM)",
        "kcat1": 25.5,
        "km1": 420.0,
        "kcat2": 38.0,
        "km2": 35.0,
        "ki": 18.0,
        "boost": 7.2,
        "f1": 0.65,
        "f2": 0.20,
        "f3": 0.15,
    },
    "sugarcane-bagasse": {
        "title": "Sugarcane Bagasse (42% Cellulose, 25% Xylan, 20% Lignin) -> 2,5-FDCA (55°C, pH 5.0)",
        "temp_c": 55.0,
        "ph": 5.0,
        "smiles": "OCC1OC(OC2C(O)C(O)C(OC2CO)O)C(O)C(O)C1O",
        "ref_smiles": "OCC1OC(O)C(O)C(O)C1O",
        "sub_label": "Cellulose/Xylan Chains (µM)",
        "int_label": "Cellobiose Intermediate (µM)",
        "prod_label": "G1P + 2,5-FDCA Sink (µM)",
        "kcat1": 23.0,
        "km1": 450.0,
        "kcat2": 40.0,
        "km2": 32.0,
        "ki": 15.0,
        "boost": 8.0,
        "f1": 0.42,
        "f2": 0.25,
        "f3": 0.20,
    },
    "corn-stover-bran": {
        "title": "Corn Stover Diferulate Bridges -> Bio-Vanillin (50°C, pH 5.5)",
        "temp_c": 50.0,
        "ph": 5.5,
        "smiles": "COc1cc(/C=C/C(=O)OCC2OC(O)C(O)C2O)ccc1O",
        "ref_smiles": "COc1cc(/C=C/C(=O)O)ccc1O",
        "sub_label": "Diferulate Ester Bridges (µM)",
        "int_label": "Free Ferulic Acid (µM)",
        "prod_label": "Bio-Vanillin Sink (µM)",
        "kcat1": 24.5,
        "km1": 390.0,
        "kcat2": 41.0,
        "km2": 30.0,
        "ki": 20.0,
        "boost": 7.8,
        "f1": 0.36,
        "f2": 0.28,
        "f3": 0.18,
    },
}


def compute_effective_boost(base_boost: float, anchor_gravy: float) -> float:
    return round(max(1.5, min(14.0, base_boost + float(anchor_gravy) * 2.5)), 1)


def load_or_fallback_pdb(uniprot_id: str, role: str) -> str:
    candidates = [
        os.path.join(MODELS_DIR, f"AF-{uniprot_id}-F1-model_v4.pdb"),
        os.path.join(MODELS_DIR, f"{uniprot_id}.pdb"),
        os.path.join(STRUCTURES_DIR, f"{uniprot_id}.pdb"),
    ]
    for path in candidates:
        if os.path.exists(path):
            with open(path, "r") as f:
                return f.read()

    os.makedirs(STRUCTURES_DIR, exist_ok=True)
    cached_path = os.path.join(STRUCTURES_DIR, f"AF-{uniprot_id}.pdb")
    if os.path.exists(cached_path):
        with open(cached_path, "r") as f:
            return f.read()

    for ver in ("v6", "v4"):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{ver}.pdb"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "pg_bio_swarm_web/1.0"})
            with urllib.request.urlopen(req, timeout=2.5) as resp:
                pdb_text = resp.read().decode("utf-8", errors="ignore")
                if "ATOM  " in pdb_text:
                    with open(cached_path, "w") as f:
                        f.write(pdb_text)
                    return pdb_text
        except Exception:
            pass

    fallback_file = "hydrophobin.pdb" if role == "anchor" else "dehalogenase.pdb"
    with open(os.path.join(REPO_ROOT, fallback_file), "r") as f:
        return f.read()


HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="color-scheme" content="light dark">
    <title>pg_bio Layer 14 — Bioprocess Swarm &amp; RK4 Cascade Digital Twin</title>
    <script src="https://cdn.jsdelivr.net/npm/echarts@5.5.0/dist/echarts.min.js"></script>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        :root {
            --bg-light: #f8fafc;
            --bg-dark: #090d16;
            --panel-light: #ffffff;
            --panel-dark: #0f172a;
            --header-light: #f1f5f9;
            --header-dark: #020617;
            --text-light: #0f172a;
            --text-dark: #f8fafc;
            --muted-light: #475569;
            --muted-dark: #94a3b8;
            --border-light: #cbd5e1;
            --border-dark: #1e293b;
            --accent-light: #0284c7;
            --accent-dark: #38bdf8;

            --bg-color: var(--bg-dark);
            --panel-bg: var(--panel-dark);
            --header-bg: var(--header-dark);
            --text-color: var(--text-dark);
            --muted-color: var(--muted-dark);
            --border-color: var(--border-dark);
            --accent-color: var(--accent-dark);

            @media (prefers-color-scheme: light) {
                --bg-color: var(--bg-dark);
            }

            @supports (color: light-dark(white, black)) {
                --accent-color: light-dark(var(--accent-light), var(--accent-dark));
            }

            color-scheme: light dark;
            accent-color: var(--accent-color);
        }

        * { box-sizing: border-box; }
        body {
            background-color: var(--bg-color);
            color: var(--text-color);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            margin: 0;
            display: flex;
            flex-direction: column;
            min-height: 100vh;
            overflow-x: hidden;
        }

        .topbar {
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 12px;
            padding: 10px 20px;
            background: var(--header-bg);
            border-bottom: 1px solid var(--border-color);
        }
        .topbar h1 { margin: 0; font-size: 1.08rem; color: #f8fafc; display: inline-block; }
        .topbar .subtitle { color: var(--muted-color); font-size: 0.82rem; margin-left: 10px; }
        .controls { display: flex; align-items: center; flex-wrap: wrap; gap: 10px; }
        select, button {
            background: #1e293b;
            color: #f8fafc;
            border: 1px solid #38bdf8;
            padding: 6px 12px;
            border-radius: 6px;
            font-size: 0.82rem;
            cursor: pointer;
        }
        button:hover { background: #38bdf8; color: #020617; font-weight: 600; }

        .dashboard-grid {
            display: grid;
            grid-template-columns: minmax(0, 1.08fr) minmax(0, 1.22fr);
            gap: 14px;
            padding: 14px;
            flex: 1;
            width: 100%;
            max-width: 100vw;
        }
        @media (max-width: 1180px) {
            .dashboard-grid { grid-template-columns: minmax(0, 1fr); }
        }

        .card {
            min-width: 0;
            background: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 10px;
            padding: 12px 14px;
            display: flex;
            flex-direction: column;
            gap: 10px;
            overflow: hidden;
        }
        .card-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 8px;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 8px;
        }
        .title-group {
            display: inline-flex;
            align-items: center;
            gap: 7px;
            min-width: 0;
        }
        .card-title { font-weight: 700; font-size: 0.93rem; color: #38bdf8; }
        .info-btn {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 20px;
            height: 20px;
            padding: 0;
            border-radius: 999px;
            background: rgba(56, 189, 248, 0.14);
            color: #38bdf8;
            border: 1.5px solid #38bdf8;
            font-family: Georgia, "Times New Roman", serif;
            font-style: italic;
            font-size: 0.78rem;
            font-weight: 700;
            line-height: 1;
            cursor: pointer;
            flex-shrink: 0;
            transition: background 0.15s, color 0.15s;
        }
        .info-btn:hover, .info-btn:focus-visible {
            background: #38bdf8;
            color: #020617;
        }
        .info-popover {
            max-width: min(500px, 92vw);
            background: #0f172a;
            color: #f8fafc;
            border: 1px solid #38bdf8;
            border-radius: 10px;
            padding: 14px 16px;
            box-shadow: 0 20px 50px rgba(0, 0, 0, 0.85);
            font-size: 0.82rem;
            line-height: 1.55;
        }
        [popover]:is(:popover-open, .\\:popover-open) {
            display: block;
        }
        .info-popover::backdrop {
            background: rgba(2, 6, 23, 0.55);
        }
        .info-popover h3 {
            margin: 0 0 8px 0;
            font-size: 0.94rem;
            color: #38bdf8;
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid #1e293b;
            padding-bottom: 6px;
        }
        .info-popover ul {
            margin: 6px 0 0 0;
            padding-left: 18px;
        }
        .info-popover li { margin-bottom: 5px; color: #cbd5e1; }
        .info-popover li b { color: #f8fafc; }
        .info-popover .formula-note {
            margin-top: 8px;
            padding: 6px 9px;
            background: #020617;
            border: 1px dashed #334155;
            border-radius: 6px;
            color: #a5f3fc;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            font-size: 0.75rem;
        }
        .badge {
            background: #020617;
            border: 1px solid #334155;
            color: #34d399;
            font-size: 0.74rem;
            padding: 3px 9px;
            border-radius: 999px;
            white-space: nowrap;
            flex-shrink: 0;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }

        #pareto-heatmap { width: 100%; height: 400px; }
        #rk4-chart { width: 100%; height: 310px; }

        .sliders-row {
            display: grid;
            grid-template-columns: repeat(5, minmax(0, 1fr));
            gap: 10px;
            background: #020617;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 9px 11px;
        }
        @media (max-width: 900px) {
            .sliders-row { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        }
        .slider-box {
            min-width: 0;
            display: flex;
            flex-direction: column;
            gap: 3px;
            font-size: 0.75rem;
            color: #cbd5e1;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }
        .slider-box b { color: #38bdf8; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
        .slider-box input[type="range"] { width: 100%; margin: 2px 0 0 0; }

        .bottom-split {
            display: grid;
            grid-template-columns: minmax(0, 1.05fr) minmax(0, 0.95fr);
            gap: 12px;
        }
        .bottom-split > div { min-width: 0; }
        @media (max-width: 900px) {
            .bottom-split { grid-template-columns: minmax(0, 1fr); }
        }

        #mol-viewer {
            width: 100%;
            height: 290px;
            position: relative;
            background: #020617;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            overflow: hidden;
        }

        .kpi-strip {
            display: grid;
            grid-template-columns: repeat(4, minmax(0, 1fr));
            gap: 8px;
        }
        @media (max-width: 900px) {
            .kpi-strip { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        }
        .kpi-item {
            min-width: 0;
            background: #020617;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 7px 9px;
            font-size: 0.74rem;
        }
        .kpi-item .k-label {
            color: var(--muted-color);
            display: block;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }
        .kpi-item .k-val {
            color: #f8fafc;
            font-weight: 700;
            font-size: 0.88rem;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            display: block;
            margin-top: 2px;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }

        /* Responsive Table per modern-web-guidance `responsive-table` */
        .table-wrapper {
            overflow: auto;
            width: 100%;
            max-inline-size: 100%;
            max-block-size: min(320px, 60vh);
            container-type: inline-size;
            border: 1px solid var(--border-color);
            border-radius: 8px;
        }
        table {
            width: 100%;
            border-collapse: separate;
            border-spacing: 0;
            font-size: 0.79rem;
            --label-1: "Rank";
            --label-2: "Anchor × Scissor";
            --label-3: "GRAVY";
            --label-4: "Pareto";
            --label-5: "48h RK4 Yield";
        }
        caption {
            text-align: left;
            padding: 8px 12px;
            font-weight: 600;
            color: #cbd5e1;
            background: #020617;
            border-bottom: 1px solid var(--border-color);
        }
        thead {
            position: sticky;
            inset-block-start: 0;
            z-index: 2;
            background: #020617;
        }
        th, td {
            padding: 7px 10px;
            text-align: left;
            border-bottom: 1px solid #1e293b;
            white-space: nowrap;
        }
        tbody th {
            position: sticky;
            inset-inline-start: 0;
            background: #0f172a;
            color: #38bdf8;
            z-index: 1;
        }
        tbody tr { cursor: pointer; transition: background 0.12s; }
        tbody tr:hover { background: #1e293b; }
        tbody tr.selected-row { background: rgba(56, 189, 248, 0.16); }

        @container (width < 420px) {
            table, thead, tbody, tr, th, td { display: block; white-space: normal; }
            thead { display: none; }
            tr { margin-block-end: 1rem; border: 1px solid #334155; }
            tbody th {
                position: sticky;
                inset-block-start: 0;
                z-index: 2;
                background: #1e293b;
            }
            td:nth-child(2)::before { content: var(--label-2) ": "; font-weight: 600; color: #94a3b8; }
            td:nth-child(3)::before { content: var(--label-3) ": "; font-weight: 600; color: #94a3b8; }
            td:nth-child(4)::before { content: var(--label-4) ": "; font-weight: 600; color: #94a3b8; }
            td:nth-child(5)::before { content: var(--label-5) ": "; font-weight: 600; color: #94a3b8; }
        }
    </style>
</head>
<body>
    <div class="topbar">
        <div>
            <h1>🧬 pg_bio Layer 14 — Multi-Feedstock Swarm &amp; RK4 Cascade Explorer</h1>
            <span class="subtitle">Click any (Anchor × Scissor) cell to run in-database RK4 kinetics &amp; 3D rigid-body complex assembly</span>
        </div>
        <div class="controls">
            <label for="pipeline-select" style="font-size:0.82rem;color:#cbd5e1;">Feedstock Pipeline:</label>
            <select id="pipeline-select" onchange="loadPipelineSwarm()">
                <option value="ptfe-to-bioceramic" selected>1. PTFE / PFAS -> Fluorapatite (37°C, pH 7.4)</option>
                <option value="plastics-upcycling">2. PET Polyester -> PCA/Vanillin (68°C, pH 8.0)</option>
                <option value="sugarcane-bagasse">3. Sugarcane Bagasse -> 2,5-FDCA (55°C, pH 5.0)</option>
                <option value="corn-stover-bran">4. Corn Stover Diferulate -> Bio-Vanillin (50°C, pH 5.5)</option>
            </select>
            <button onclick="resetSlidersToPreset()">Reset Preset</button>
        </div>
    </div>

    <div class="dashboard-grid">
        <!-- Left Column: Pareto Swarm Heatmap + Accessible Top Candidates Table -->
        <div class="card">
            <div class="card-header">
                <div class="title-group">
                    <span class="card-title">1. In-Database Pareto Swarm Matrix (`pareto_score_domain_pair`)</span>
                    <button type="button" class="info-btn" commandfor="info-pareto-matrix" command="toggle-popover" aria-label="How to read the Pareto Swarm Matrix">i</button>
                </div>
                <span class="badge" id="swarm-timing-badge">PostgreSQL 18 Live</span>
            </div>
            <div id="pareto-heatmap"></div>

            <div tabindex="0" role="region" aria-labelledby="pareto-table-caption" class="table-wrapper">
                <table>
                    <caption id="pareto-table-caption">
                        <div style="display:flex;justify-content:space-between;align-items:center;gap:8px;width:100%;">
                            <span>Top Pareto-Optimal (M1 Anchor × M2 Scissor) Candidates — Click Row to Inspect</span>
                            <button type="button" class="info-btn" commandfor="info-leaderboard-table" command="toggle-popover" aria-label="How to read the Top Candidates Leaderboard">i</button>
                        </div>
                    </caption>
                    <thead>
                        <tr>
                            <th scope="col">Rank</th>
                            <th scope="col">Anchor × Scissor</th>
                            <th scope="col">Anchor GRAVY</th>
                            <th scope="col">Pareto Score</th>
                            <th scope="col">48h RK4 Yield</th>
                        </tr>
                    </thead>
                    <tbody id="leaderboard-body"></tbody>
                </table>
            </div>
        </div>

        <!-- Right Column: Live RK4 Cascade Digital Twin + Sliders + 3D Complex & Stoichiometry -->
        <div class="card">
            <div class="card-header">
                <div class="title-group">
                    <span class="card-title" id="rk4-title">2. 48-Hour RK4 Cascade Digital Twin (`simulate_cascade_kinetics`)</span>
                    <button type="button" class="info-btn" commandfor="info-rk4-twin" command="toggle-popover" aria-label="How to read the 48-Hour RK4 Cascade Digital Twin">i</button>
                </div>
                <div style="display:flex;align-items:center;gap:6px;">
                    <button type="button" class="info-btn" commandfor="info-sliders-kpis" command="toggle-popover" title="How to use the kinetic sliders &amp; KPI cards" aria-label="How to read the kinetic sliders and KPI strip">i</button>
                    <span class="badge" id="rk4-timing-badge">RK4 Solver</span>
                </div>
            </div>

            <div class="sliders-row">
                <div class="slider-box">
                    <span>Scissor kcat1: <b id="lbl-kcat1">24.0</b> /h</span>
                    <input type="range" id="sl-kcat1" min="5" max="60" step="0.5" value="24.0" oninput="updateRK4FromSliders()">
                </div>
                <div class="slider-box">
                    <span>Solid Km1: <b id="lbl-km1">380</b> µM</span>
                    <input type="range" id="sl-km1" min="50" max="800" step="10" value="380" oninput="updateRK4FromSliders()">
                </div>
                <div class="slider-box">
                    <span>Sink kcat2: <b id="lbl-kcat2">42.0</b> /h</span>
                    <input type="range" id="sl-kcat2" min="0" max="80" step="1.0" value="42.0" oninput="updateRK4FromSliders()">
                </div>
                <div class="slider-box">
                    <span>Inhibition Ki: <b id="lbl-ki">22.0</b> µM</span>
                    <input type="range" id="sl-ki" min="5" max="150" step="1.0" value="22.0" oninput="updateRK4FromSliders()">
                </div>
                <div class="slider-box">
                    <span>Anchor Boost: <b id="lbl-boost">8.5</b>×</span>
                    <input type="range" id="sl-boost" min="1.0" max="15.0" step="0.1" value="8.5" oninput="updateRK4FromSliders()">
                </div>
            </div>

            <div class="kpi-strip">
                <div class="kpi-item">
                    <span class="k-label">Selected Chimera Pair</span>
                    <span class="k-val" id="kpi-pair">—</span>
                </div>
                <div class="kpi-item">
                    <span class="k-label">48h Full 4-Module Sink Yield</span>
                    <span class="k-val" id="kpi-sink-yield" style="color:#34d399;">—</span>
                </div>
                <div class="kpi-item">
                    <span class="k-label">48h No-Sink (Ki Stalled)</span>
                    <span class="k-val" id="kpi-nosink-yield" style="color:#facc15;">—</span>
                </div>
                <div class="kpi-item">
                    <span class="k-label">Consortium Ratio (w1:w2:w3)</span>
                    <span class="k-val" id="kpi-stoich" style="color:#38bdf8;">—</span>
                </div>
            </div>

            <div id="rk4-chart"></div>

            <div class="bottom-split">
                <div>
                    <div style="display:flex;justify-content:space-between;align-items:center;gap:8px;margin-bottom:6px;flex-wrap:wrap;">
                        <div class="title-group">
                            <span style="font-size:0.82rem;font-weight:600;color:#38bdf8;">3D Rigid-Body Fusion (`translate_and_merge_pdb_chains`)</span>
                            <button type="button" class="info-btn" commandfor="info-3d-fusion" command="toggle-popover" aria-label="How to read the 3D Rigid-Body Fusion Viewer">i</button>
                        </div>
                        <span class="badge" id="globular-badge">Rg: —</span>
                    </div>
                    <div id="mol-viewer"></div>
                </div>
                <div style="display:flex;flex-direction:column;justify-content:space-between;background:#020617;border:1px solid #1e293b;border-radius:8px;padding:12px;">
                    <div>
                        <div style="display:flex;justify-content:space-between;align-items:center;gap:8px;margin-bottom:6px;">
                            <span style="font-size:0.84rem;font-weight:700;color:#38bdf8;">Live PostgreSQL 18 Biophysical Ledger</span>
                            <button type="button" class="info-btn" commandfor="info-sql-ledger" command="toggle-popover" aria-label="How to read the Live PostgreSQL 18 Biophysical Ledger">i</button>
                        </div>
                        <div id="pair-details" style="font-size:0.79rem;color:#cbd5e1;line-height:1.55;"></div>
                    </div>
                    <pre id="sql-preview" style="background:#090d16;border:1px solid #1e293b;color:#a5f3fc;font-size:0.72rem;padding:8px;border-radius:6px;overflow-x:auto;margin:8px 0 0 0;"></pre>
                </div>
            </div>
        </div>
    </div>

    <!-- Declarative Popover Instruction Guides for Every Widget -->
    <div id="info-pareto-matrix" class="info-popover" popover>
        <h3>
            <span>How to Read: 1. In-Database Pareto Swarm Matrix</span>
            <button type="button" commandfor="info-pareto-matrix" command="hide-popover" style="padding:2px 8px;font-size:0.75rem;">Close</button>
        </h3>
        <ul>
            <li><b>Y-Axis (Rows):</b> 12 stratified <b>Module 1 Surface Anchors</b> (Hydrophobins, CBMs, Dockerins) sampled across the hydropathy (GRAVY) spectrum in <code>bio_domain_atlas</code>.</li>
            <li><b>X-Axis (Columns):</b> 16 stratified <b>Module 2 Primary Bond Scissors</b> (Dehalogenases, Cutinases/PETases, Cellulases/LPMOs, Feruloyl Esterases) sampled by isoelectric pI match to the reactor pH.</li>
            <li><b>Cell Color (Dark Blue → Emerald → Gold):</b> Multi-objective Pareto score computed in Rust inside PostgreSQL via <code>pareto_score_domain_pair()</code>. Gold/Green cells (`> 0.95`) combine high interfacial binding with optimal reactor pH compatibility.</li>
            <li><b>Action:</b> Click any cell to load that <code>(Anchor × Scissor)</code> pair into the RK4 Cascade Digital Twin and 3D Structure Viewer.</li>
        </ul>
        <div class="formula-note">Score = 0.30·GRAVY_norm + 0.30·Entropy_norm + 0.25·pH_match + 0.15·MW_compact</div>
    </div>

    <div id="info-leaderboard-table" class="info-popover" popover>
        <h3>
            <span>How to Read: Top Pareto-Optimal Leaderboard</span>
            <button type="button" commandfor="info-leaderboard-table" command="hide-popover" style="padding:2px 8px;font-size:0.75rem;">Close</button>
        </h3>
        <ul>
            <li><b>Anchor GRAVY:</b> Kyte-Doolittle Grand Average of Hydropathy for the Module 1 domain. Higher values (`> +0.3`) drive stronger hydrophobic adsorption onto solid plastics/fluoropolymers.</li>
            <li><b>Pareto Score:</b> Overall multi-objective biophysical fitness (`0.0000 – 1.0000`).</li>
            <li><b>48h RK4 Yield:</b> Exact 48-hour sink product conversion (`%` and `µM`) solved inside PostgreSQL via <code>simulate_cascade_kinetics()</code> using the anchor's GRAVY-adjusted adsorption boost. Matches the top-right KPI card 100%.</li>
            <li><b>Action:</b> Click any row (or scroll with keyboard) to inspect that candidate's kinetic curves and 3D structure.</li>
        </ul>
    </div>

    <div id="info-rk4-twin" class="info-popover" popover>
        <h3>
            <span>How to Read: 2. 48-Hour RK4 Cascade Digital Twin</span>
            <button type="button" commandfor="info-rk4-twin" command="hide-popover" style="padding:2px 8px;font-size:0.75rem;">Close</button>
        </h3>
        <ul>
            <li><b>Solid Green Curve (With Sink):</b> Final upcycled product $P(t)$ when all 4 modules operate together. Module 3 continuously drains the intermediate so the primary scissor never stalls.</li>
            <li><b>Orange Dashed Curve (Transient Intermediate):</b> Concentration of the soluble cleavage intermediate $I(t)$ (e.g., cellobiose, MHET, or free $\text{F}^-$). Notice how Module 3 keeps it capped well below $K_i$.</li>
            <li><b>Dotted Slate Curve (Remaining Substrate):</b> Unconverted solid polymer/biomass bond equivalents $S(t)$ starting from $1{,}000\,\mu\text{M}$.</li>
            <li><b>Blue Dash-Dot Curve (Cleaved w/o Sink):</b> Control where Sink $k_{\text{cat},2} = 0$. Intermediate accumulates above $K_i$ and competitively stalls the scissor around 45–60%.</li>
            <li><b>Gray Curve (Wild-Type Alone):</b> Unfused enzyme without surface anchor (`Boost = 1.0×`) or sink (`20–27%` conversion).</li>
        </ul>
        <div class="formula-note">v1 = kcat1 · (Boost·S) / (Km1·(1 + I/Ki) + Boost·S) &nbsp;|&nbsp; v2 = kcat2 · I / (Km2 + I)</div>
    </div>

    <div id="info-sliders-kpis" class="info-popover" popover>
        <h3>
            <span>How to Read: Kinetic Sliders &amp; KPI Strip</span>
            <button type="button" commandfor="info-sliders-kpis" command="hide-popover" style="padding:2px 8px;font-size:0.75rem;">Close</button>
        </h3>
        <ul>
            <li><b>Scissor kcat1 &amp; Solid Km1:</b> Primary bond cleavage rate (`/h`) and interfacial Michaelis constant (`µM`) on the recalcitrant solid surface.</li>
            <li><b>Sink kcat2 &amp; Inhibition Ki:</b> Rate at which Module 3 clears the intermediate (`kcat2`) versus the competitive product inhibition threshold (`Ki`). Dragging `Sink kcat2` to `0` reproduces the stalled curve!</li>
            <li><b>Anchor Boost:</b> Langmuir-Hinshelwood effective local substrate concentration multiplier provided by Module 1 (`1.0×` = no anchor; `8–11×` = high-affinity surface anchor).</li>
            <li><b>Consortium Ratio (w1:w2:w3):</b> Optimal enzyme mass fractions ($w_1 + w_2 + w_3 = 1.00$) solved by <code>optimize_consortium_stoichiometry()</code> for 3-component biomass (e.g., Cellulose : Xylan : Lignin).</li>
        </ul>
    </div>

    <div id="info-3d-fusion" class="info-popover" popover>
        <h3>
            <span>How to Read: 3D Rigid-Body Fusion Viewer</span>
            <button type="button" commandfor="info-3d-fusion" command="hide-popover" style="padding:2px 8px;font-size:0.75rem;">Close</button>
        </h3>
        <ul>
            <li><b>Cyan Ribbon (Chain A):</b> Primary Catalytic Scissor domain (Module 2).</li>
            <li><b>Gold Ribbon + Sticks (Chain B):</b> Interfacial Surface Anchor domain (Module 1), translated along the principal axis by <code>translate_and_merge_pdb_chains()</code> to guarantee a clash-free <b>28.0 Å inter-domain clearance</b>.</li>
            <li><b>Top Badge (Cα | pLDDT | Rg):</b> Live in-database structural audit from <code>evaluate_globular_compactness()</code> reporting total Cα residues, mean confidence (`pLDDT`), and 3D Radius of Gyration (`Rg` in Å).</li>
            <li><b>Controls:</b> Left-drag to rotate, scroll to zoom, right-drag to pan.</li>
        </ul>
    </div>

    <div id="info-sql-ledger" class="info-popover" popover>
        <h3>
            <span>How to Read: Live PostgreSQL 18 Biophysical Ledger</span>
            <button type="button" commandfor="info-sql-ledger" command="hide-popover" style="padding:2px 8px;font-size:0.75rem;">Close</button>
        </h3>
        <ul>
            <li><b>Domain Metadata:</b> Displays UniProt accession, source organism, molecular weight (`kDa`), isoelectric point (`pI`), and hydropathy (`GRAVY`) for the active chimera.</li>
            <li><b>Flory Compactness Ratio:</b> Ratio of actual $R_g$ to ideal globular scaling $2.25 \cdot N^{0.38}$ (`~1.0–1.6` indicates a well-packed two-domain complex).</li>
            <li><b>SQL Query Preview:</b> Shows the exact PostgreSQL 18 SQL statement executed against <code>$libdir/pg_bio</code> whenever you drag a slider or click a candidate.</li>
        </ul>
    </div>

    <script type="module">
        if (!('commandForElement' in HTMLButtonElement.prototype)) {
            import('https://esm.run/invokers-polyfill');
        }
        if (!('popover' in HTMLElement.prototype)) {
            import('https://unpkg.com/@oddbird/popover-polyfill@latest/dist/popover.min.js');
        }
    </script>

    <script>
        var heatmapChart = echarts.init(document.getElementById('pareto-heatmap'));
        var rk4Chart = echarts.init(document.getElementById('rk4-chart'));
        var molViewer = null;
        var currentSwarmData = null;
        var activeAnchor = null;
        var activeScissor = null;

        window.addEventListener('resize', function() {
            heatmapChart.resize();
            rk4Chart.resize();
            if (molViewer) molViewer.resize();
        });

        function loadPipelineSwarm() {
            let pid = document.getElementById('pipeline-select').value;
            heatmapChart.showLoading({ text: 'Querying pg_bio Domain Atlas...', color: '#38bdf8', textColor: '#fff', maskColor: 'rgba(9,13,22,0.8)' });

            fetch('/api/swarm?pipeline=' + encodeURIComponent(pid))
                .then(r => r.json())
                .then(data => {
                    heatmapChart.hideLoading();
                    currentSwarmData = data;
                    document.getElementById('swarm-timing-badge').textContent = `pg_bio: ${data.elapsed_ms} ms (${data.anchors.length}×${data.scissors.length})`;
                    document.getElementById('kpi-stoich').textContent = `${data.stoichiometry.w1} : ${data.stoichiometry.w2} : ${data.stoichiometry.w3}`;

                    let p = data.preset;
                    document.getElementById('sl-kcat1').value = p.kcat1;
                    document.getElementById('sl-km1').value = p.km1;
                    document.getElementById('sl-kcat2').value = p.kcat2;
                    document.getElementById('sl-ki').value = p.ki;
                    document.getElementById('sl-boost').value = p.boost;
                    syncSliderLabels();

                    renderHeatmap(data);
                    renderLeaderboard(data.top_pairs);

                    if (data.top_pairs.length > 0) {
                        selectPair(data.top_pairs[0].anchor_id, data.top_pairs[0].scissor_id);
                    }
                });
        }

        function resetSlidersToPreset() {
            if (!currentSwarmData) return;
            let p = currentSwarmData.preset;
            document.getElementById('sl-kcat1').value = p.kcat1;
            document.getElementById('sl-km1').value = p.km1;
            document.getElementById('sl-kcat2').value = p.kcat2;
            document.getElementById('sl-ki').value = p.ki;
            let adjBoost = activeAnchor ? Math.max(1.5, Math.min(14.0, p.boost + Number(activeAnchor.gravy) * 2.5)) : p.boost;
            document.getElementById('sl-boost').value = adjBoost.toFixed(1);
            syncSliderLabels();
            runRK4Simulation();
        }

        function syncSliderLabels() {
            document.getElementById('lbl-kcat1').textContent = Number(document.getElementById('sl-kcat1').value).toFixed(1);
            document.getElementById('lbl-km1').textContent = Number(document.getElementById('sl-km1').value).toFixed(0);
            document.getElementById('lbl-kcat2').textContent = Number(document.getElementById('sl-kcat2').value).toFixed(1);
            document.getElementById('lbl-ki').textContent = Number(document.getElementById('sl-ki').value).toFixed(1);
            document.getElementById('lbl-boost').textContent = Number(document.getElementById('sl-boost').value).toFixed(1);
        }

        function renderHeatmap(data) {
            let xLabels = data.scissors.map(s => s.id);
            let yLabels = data.anchors.map(a => a.id);

            let option = {
                tooltip: {
                    position: 'top',
                    backgroundColor: 'rgba(2,6,23,0.95)',
                    borderColor: '#38bdf8',
                    textStyle: { color: '#fff', fontSize: 12 },
                    formatter: function(params) {
                        let s = data.scissors[params.value[0]];
                        let a = data.anchors[params.value[1]];
                        return `<b>Anchor (M1):</b> ${a.id} (${a.org}) | GRAVY: ${a.gravy}<br>` +
                               `<b>Scissor (M2):</b> ${s.id} (${s.org}) | MW: ${s.mw} kDa, pI: ${s.pi}<br>` +
                               `<span style="color:#34d399">Pareto Score: <b>${params.value[2].toFixed(4)}</b></span>`;
                    }
                },
                grid: { top: '5%', bottom: '19%', left: '16%', right: '12%' },
                xAxis: {
                    type: 'category',
                    data: xLabels,
                    axisLabel: { rotate: 45, fontSize: 10, color: '#94a3b8' }
                },
                yAxis: {
                    type: 'category',
                    data: yLabels,
                    axisLabel: { fontSize: 10, color: '#38bdf8' }
                },
                visualMap: {
                    min: data.min_score,
                    max: data.max_score,
                    calculable: true,
                    orient: 'vertical',
                    right: '1%',
                    top: 'center',
                    precision: 3,
                    textStyle: { color: '#cbd5e1', fontSize: 10 },
                    inRange: { color: ['#0f172a', '#1d4ed8', '#0284c7', '#10b981', '#facc15'] }
                },
                series: [{
                    name: 'Pareto Score',
                    type: 'heatmap',
                    data: data.matrix_cells,
                    itemStyle: { borderColor: '#090d16', borderWidth: 1 },
                    emphasis: { itemStyle: { borderColor: '#fff', borderWidth: 2 } }
                }]
            };
            heatmapChart.setOption(option);
            heatmapChart.off('click');
            heatmapChart.on('click', function(params) {
                if (!params.value) return;
                let s = data.scissors[params.value[0]];
                let a = data.anchors[params.value[1]];
                selectPair(a.id, s.id);
            });
        }

        function renderLeaderboard(topPairs) {
            let tbody = document.getElementById('leaderboard-body');
            tbody.innerHTML = '';
            topPairs.forEach((row, idx) => {
                let tr = document.createElement('tr');
                tr.dataset.pair = `${row.anchor_id}_${row.scissor_id}`;
                tr.onclick = () => selectPair(row.anchor_id, row.scissor_id);
                tr.innerHTML = `
                    <th scope="row">#${idx + 1}</th>
                    <td><b>${row.anchor_id}</b> × <b>${row.scissor_id}</b> <span style="color:#94a3b8;font-size:0.73rem;">(${row.scissor_org})</span></td>
                    <td>${row.anchor_gravy.toFixed(3)}</td>
                    <td style="color:#38bdf8;font-weight:700;">${row.pareto_score.toFixed(4)}</td>
                    <td style="color:#34d399;font-weight:700;">${row.yield_pct.toFixed(1)}% (${row.yield_um} µM)</td>
                `;
                tbody.appendChild(tr);
            });
        }

        function selectPair(anchorId, scissorId) {
            activeAnchor = currentSwarmData.anchors.find(a => a.id === anchorId) || { id: anchorId, gravy: 0.3, org: 'Unknown' };
            activeScissor = currentSwarmData.scissors.find(s => s.id === scissorId) || { id: scissorId, mw: 30.0, pi: 6.5, org: 'Unknown' };

            document.querySelectorAll('#leaderboard-body tr').forEach(tr => {
                tr.classList.toggle('selected-row', tr.dataset.pair === `${anchorId}_${scissorId}`);
            });

            document.getElementById('kpi-pair').textContent = `${anchorId} × ${scissorId}`;

            let p = currentSwarmData.preset;
            document.getElementById('sl-kcat1').value = p.kcat1;
            document.getElementById('sl-km1').value = p.km1;
            document.getElementById('sl-kcat2').value = p.kcat2;
            document.getElementById('sl-ki').value = p.ki;
            let adjBoost = Math.max(1.5, Math.min(14.0, p.boost + Number(activeAnchor.gravy) * 2.5));
            document.getElementById('sl-boost').value = adjBoost.toFixed(1);
            syncSliderLabels();

            runRK4Simulation();
            load3DComplex(anchorId, scissorId);
        }

        function updateRK4FromSliders() {
            syncSliderLabels();
            runRK4Simulation();
        }

        function runRK4Simulation() {
            if (!currentSwarmData) return;
            let pid = document.getElementById('pipeline-select').value;
            let kcat1 = document.getElementById('sl-kcat1').value;
            let km1 = document.getElementById('sl-km1').value;
            let kcat2 = document.getElementById('sl-kcat2').value;
            let ki = document.getElementById('sl-ki').value;
            let boost = document.getElementById('sl-boost').value;

            let qs = new URLSearchParams({
                pipeline: pid, kcat1: kcat1, km1: km1, kcat2: kcat2, ki: ki, boost: boost
            });

            fetch('/api/rk4?' + qs.toString())
                .then(r => r.json())
                .then(res => {
                    document.getElementById('rk4-timing-badge').textContent = `RK4 SQL: ${res.elapsed_ms} ms`;
                    document.getElementById('kpi-sink-yield').textContent = `${res.sink_48h_pct.toFixed(1)}% (${res.sink_48h_um} µM)`;
                    document.getElementById('kpi-nosink-yield').textContent = `${res.nosink_48h_pct.toFixed(1)}%`;

                    let p = currentSwarmData.preset;
                    let option = {
                        tooltip: {
                            trigger: 'axis',
                            backgroundColor: 'rgba(2,6,23,0.95)',
                            borderColor: '#38bdf8',
                            textStyle: { color: '#fff', fontSize: 12 }
                        },
                        legend: {
                            top: 0,
                            textStyle: { color: '#cbd5e1', fontSize: 10.5 }
                        },
                        grid: { top: '16%', bottom: '12%', left: '7%', right: '7%' },
                        xAxis: {
                            type: 'category',
                            data: res.time_hr,
                            name: 'Hours',
                            nameLocation: 'end',
                            axisLabel: { color: '#94a3b8' }
                        },
                        yAxis: {
                            type: 'value',
                            min: 0,
                            max: 1020,
                            name: 'µM eq',
                            axisLabel: { color: '#94a3b8' },
                            splitLine: { lineStyle: { color: '#1e293b' } }
                        },
                        series: [
                            {
                                name: p.prod_label + ' (With Sink)',
                                type: 'line',
                                data: res.with_sink_prod,
                                smooth: true,
                                lineStyle: { width: 3, color: '#34d399' },
                                areaStyle: { color: 'rgba(52, 211, 153, 0.14)' }
                            },
                            {
                                name: 'Transient ' + p.int_label,
                                type: 'line',
                                data: res.with_sink_int,
                                smooth: true,
                                lineStyle: { width: 2, type: 'dashed', color: '#f59e0b' }
                            },
                            {
                                name: 'Remaining ' + p.sub_label,
                                type: 'line',
                                data: res.with_sink_sub,
                                smooth: true,
                                lineStyle: { width: 1.8, type: 'dotted', color: '#64748b' }
                            },
                            {
                                name: 'Cleaved w/o Sink (Ki Stalled)',
                                type: 'line',
                                data: res.nosink_cleaved,
                                smooth: true,
                                lineStyle: { width: 2, type: 'dashdot', color: '#38bdf8' }
                            },
                            {
                                name: 'Wild-Type Scissor Alone',
                                type: 'line',
                                data: res.wt_cleaved,
                                smooth: true,
                                lineStyle: { width: 1.5, color: '#94a3b8' }
                            }
                        ]
                    };
                    rk4Chart.setOption(option);

                    document.getElementById('sql-preview').textContent =
`SELECT time_hr, substrate_u_m, intermediate_u_m, sink_product_u_m, fractional_conversion
FROM simulate_cascade_kinetics(${kcat1}, ${km1}, ${kcat2}, ${p.km2}, ${ki}, ${boost}, 1000.0, 48.0, 48);`;
                });
        }

        function load3DComplex(anchorId, scissorId) {
            fetch(`/api/complex3d?anchor=${encodeURIComponent(anchorId)}&scissor=${encodeURIComponent(scissorId)}`)
                .then(r => r.json())
                .then(data => {
                    document.getElementById('globular-badge').textContent =
                        `Cα: ${data.globular.ca_count} | pLDDT: ${data.globular.mean_plddt.toFixed(1)} | Rg: ${data.globular.radius_of_gyration.toFixed(1)} Å`;

                    document.getElementById('pair-details').innerHTML =
                        `<b>M1 Surface Anchor:</b> <span style="color:#f59e0b">${anchorId}</span> (${activeAnchor.org}) — GRAVY: <b>${activeAnchor.gravy}</b><br>` +
                        `<b>M2 Primary Scissor:</b> <span style="color:#38bdf8">${scissorId}</span> (${activeScissor.org}) — MW: <b>${activeScissor.mw} kDa</b>, pI: <b>${activeScissor.pi}</b><br>` +
                        `<b>M3 Equilibrium Sink:</b> Drains intermediate below K<sub>i</sub> (${document.getElementById('sl-ki').value} µM)<br>` +
                        `<b>3D Globular Audit:</b> Compactness Ratio <b>${data.globular.compactness_ratio.toFixed(3)}</b> | Elongation <b>${data.globular.elongation.toFixed(2)}</b> (${data.elapsed_ms} ms)`;

                    let el = document.getElementById('mol-viewer');
                    if (!molViewer) {
                        molViewer = $3Dmol.createViewer(el, { backgroundColor: '#020617' });
                    }
                    molViewer.clear();
                    molViewer.addModel(data.merged_pdb, 'pdb');
                    molViewer.setStyle({ chain: 'A' }, { cartoon: { color: '#38bdf8' } });
                    molViewer.setStyle({ chain: 'B' }, { cartoon: { color: '#f59e0b' }, stick: { radius: 0.15, color: '#f59e0b' } });
                    molViewer.zoomTo();
                    molViewer.zoom(1.18);
                    molViewer.render();
                });
        }

        loadPipelineSwarm();
    </script>
</body>
</html>
"""


@app.route("/")
def index():
    return render_template_string(HTML_TEMPLATE)


@app.route("/api/swarm")
def api_swarm():
    pid = request.args.get("pipeline", "ptfe-to-bioceramic")
    if pid not in PIPELINE_PRESETS:
        pid = "ptfe-to-bioceramic"
    preset = PIPELINE_PRESETS[pid]

    t0 = time.perf_counter()
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # Sample 12 diverse anchors spanning the full hydropathy range so heatmap shows rich Pareto contrast
            cur.execute(
                """
                WITH ranked_a AS (
                    SELECT uniprot_id, organism_name,
                           ROUND(gravy_hydropathy::numeric, 3) AS gravy,
                           ROUND(shannon_entropy::numeric, 3) AS entropy,
                           NTILE(12) OVER (ORDER BY gravy_hydropathy DESC, shannon_entropy DESC) AS bucket
                    FROM bio_domain_atlas
                    WHERE pipeline_family = %s AND module_role = 'M1_Anchor'
                )
                SELECT DISTINCT ON (bucket) uniprot_id, organism_name, gravy, entropy
                FROM ranked_a
                ORDER BY bucket ASC, gravy DESC;
                """,
                (pid,),
            )
            anchors = [
                {"id": r[0], "org": r[1], "gravy": float(r[2]), "entropy": float(r[3])}
                for r in cur.fetchall()
            ]

            # Sample 16 diverse scissors spanning pI and MW
            cur.execute(
                """
                WITH ranked_s AS (
                    SELECT uniprot_id, organism_name,
                           ROUND(mw_kda::numeric, 1) AS mw,
                           ROUND(isoelectric_pi::numeric, 2) AS pi,
                           ROUND(shannon_entropy::numeric, 3) AS entropy,
                           NTILE(16) OVER (ORDER BY ABS(isoelectric_pi - %s) ASC, shannon_entropy DESC) AS bucket
                    FROM bio_domain_atlas
                    WHERE pipeline_family = %s AND module_role = 'M2_Scissor'
                      AND mw_kda BETWEEN 16.0 AND 75.0
                )
                SELECT DISTINCT ON (bucket) uniprot_id, organism_name, mw, pi, entropy
                FROM ranked_s
                ORDER BY bucket ASC, entropy DESC;
                """,
                (preset["ph"], pid),
            )
            scissors = [
                {"id": r[0], "org": r[1], "mw": float(r[2]), "pi": float(r[3]), "entropy": float(r[4])}
                for r in cur.fetchall()
            ]

            matrix_cells = []
            all_pairs = []
            scores = []
            for y_idx, a in enumerate(anchors):
                for x_idx, s in enumerate(scissors):
                    cur.execute(
                        "SELECT pareto_score_domain_pair(%s, %s, %s, %s, %s, %s);",
                        (a["gravy"], a["entropy"], s["mw"], s["pi"], s["entropy"], preset["ph"]),
                    )
                    score = round(float(cur.fetchone()[0]), 4)
                    scores.append(score)
                    matrix_cells.append([x_idx, y_idx, score])
                    all_pairs.append(
                        {
                            "anchor_id": a["id"],
                            "anchor_org": a["org"],
                            "anchor_gravy": a["gravy"],
                            "scissor_id": s["id"],
                            "scissor_org": s["org"],
                            "pareto_score": score,
                        }
                    )

            all_pairs.sort(key=lambda x: x["pareto_score"], reverse=True)
            top_pairs = all_pairs[:10]

            # Compute exact RK4 48h yield in PostgreSQL for each top-10 pair so table & KPI match 100%!
            for row in top_pairs:
                eff_boost = compute_effective_boost(preset["boost"], row["anchor_gravy"])
                cur.execute(
                    """
                    SELECT sink_product_u_m, fractional_conversion
                    FROM simulate_cascade_kinetics(%s, %s, %s, %s, %s, %s, 1000.0, 48.0, 48)
                    ORDER BY time_hr DESC
                    LIMIT 1;
                    """,
                    (
                        preset["kcat1"],
                        preset["km1"],
                        preset["kcat2"],
                        preset["km2"],
                        preset["ki"],
                        eff_boost,
                    ),
                )
                sink_um, frac_conv = cur.fetchone()
                row["yield_pct"] = round(float(frac_conv) * 100.0, 1)
                row["yield_um"] = int(round(float(sink_um)))
                row["effective_boost"] = eff_boost

            cur.execute(
                "SELECT * FROM optimize_consortium_stoichiometry(%s, %s, %s, 14.0, 18.0, 10.5, 48.0);",
                (preset["f1"], preset["f2"], preset["f3"]),
            )
            w1, w2, w3, conv, bal = cur.fetchone()

    elapsed_ms = round((time.perf_counter() - t0) * 1000.0, 2)
    min_score = round(min(scores) if scores else 0.5, 3)
    max_score = round(max(scores) if scores else 1.0, 3)

    return jsonify(
        {
            "pipeline": pid,
            "preset": preset,
            "anchors": anchors,
            "scissors": scissors,
            "matrix_cells": matrix_cells,
            "min_score": min_score,
            "max_score": max_score,
            "top_pairs": top_pairs,
            "stoichiometry": {
                "w1": f"{w1:.2f}",
                "w2": f"{w2:.2f}",
                "w3": f"{w3:.2f}",
                "conv_pct": round(conv * 100.0, 1),
                "balance": round(bal, 4),
            },
            "elapsed_ms": elapsed_ms,
        }
    )


@app.route("/api/rk4")
def api_rk4():
    pid = request.args.get("pipeline", "ptfe-to-bioceramic")
    preset = PIPELINE_PRESETS.get(pid, PIPELINE_PRESETS["ptfe-to-bioceramic"])

    kcat1 = float(request.args.get("kcat1", preset["kcat1"]))
    km1 = float(request.args.get("km1", preset["km1"]))
    kcat2 = float(request.args.get("kcat2", preset["kcat2"]))
    ki = float(request.args.get("ki", preset["ki"]))
    boost = float(request.args.get("boost", preset["boost"]))
    km2 = preset["km2"]

    t0 = time.perf_counter()
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT time_hr, substrate_u_m, intermediate_u_m, sink_product_u_m, fractional_conversion
                FROM simulate_cascade_kinetics(%s, %s, %s, %s, %s, %s, 1000.0, 48.0, 48)
                ORDER BY time_hr ASC;
                """,
                (kcat1, km1, kcat2, km2, ki, boost),
            )
            ws = cur.fetchall()

            cur.execute(
                """
                SELECT time_hr, substrate_u_m
                FROM simulate_cascade_kinetics(%s, %s, 0.0, %s, %s, %s, 1000.0, 48.0, 48)
                ORDER BY time_hr ASC;
                """,
                (kcat1, km1, km2, ki, boost),
            )
            ns = cur.fetchall()

            cur.execute(
                """
                SELECT time_hr, substrate_u_m
                FROM simulate_cascade_kinetics(%s, %s, 0.0, %s, %s, 1.0, 1000.0, 48.0, 48)
                ORDER BY time_hr ASC;
                """,
                (kcat1, km1, km2, ki),
            )
            wt = cur.fetchall()

    elapsed_ms = round((time.perf_counter() - t0) * 1000.0, 2)
    return jsonify(
        {
            "time_hr": [int(r[0]) for r in ws],
            "with_sink_sub": [round(r[1], 1) for r in ws],
            "with_sink_int": [round(r[2], 1) for r in ws],
            "with_sink_prod": [round(r[3], 1) for r in ws],
            "nosink_cleaved": [round(1000.0 - r[1], 1) for r in ns],
            "wt_cleaved": [round(1000.0 - r[1], 1) for r in wt],
            "sink_48h_pct": round(ws[-1][4] * 100.0, 1),
            "sink_48h_um": int(round(ws[-1][3])),
            "nosink_48h_pct": round((1000.0 - ns[-1][1]) / 10.0, 1),
            "elapsed_ms": elapsed_ms,
        }
    )


@app.route("/api/complex3d")
def api_complex3d():
    anchor_id = request.args.get("anchor", "P79073")
    scissor_id = request.args.get("scissor", "Q1JU72")

    pdb_scissor = load_or_fallback_pdb(scissor_id, "scissor")
    pdb_anchor = load_or_fallback_pdb(anchor_id, "anchor")

    t0 = time.perf_counter()
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                WITH merged AS (
                    SELECT translate_and_merge_pdb_chains(%s, %s, 28.0) AS pdb_str
                )
                SELECT m.pdb_str,
                       evaluate_globular_compactness(m.pdb_str)::text
                FROM merged m;
                """,
                (pdb_scissor, pdb_anchor),
            )
            merged_pdb, glob_json = cur.fetchone()

    elapsed_ms = round((time.perf_counter() - t0) * 1000.0, 2)
    return jsonify(
        {
            "anchor_id": anchor_id,
            "scissor_id": scissor_id,
            "merged_pdb": merged_pdb,
            "globular": json.loads(glob_json),
            "elapsed_ms": elapsed_ms,
        }
    )


if __name__ == "__main__":
    port = int(os.getenv("SWARM_WEB_PORT", "8090"))
    app.run(host="127.0.0.1", port=port, debug=False)
