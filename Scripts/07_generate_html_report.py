#!/usr/bin/env python3
import os
import sys
import json
import math
import pandas as pd

# ------------------------------------------------------------------
# Usage:
#   python3 07_generate_html_report.py MR_Formatted_Results.csv MR_IVW_OR_Per_SNP.csv style.html MR_CoPe_Report.html
# Expects the following files (if present) in the SAME directory as the output:
#   exposure_manhattan.png, outcome_manhattan.png, exposure_qq.png, outcome_qq.png
#   harmonised_data.csv  (optional; improves SNP table & scatter)
# ------------------------------------------------------------------

summary_csv, snps_csv, html_template, output_html = sys.argv[1:]

# -------------------- Helpers --------------------
def fmt_p(x):
    try:
        x = float(x)
        if math.isnan(x):
            return "—"
        return f"{x:.3g}" if x < 1e-3 else f"{x:.4f}"
    except Exception:
        return "—"

def fmt_num(x, r=3):
    try:
        x = float(x)
        if math.isnan(x):
            return "—"
        # keep integers clean
        if abs(x - round(x)) < 1e-9:
            return str(int(round(x)))
        return f"{x:.{r}f}"
    except Exception:
        return "—"

def CI(l, u):
    if l == "—" or u == "—":
        return "—"
    return f"{l}–{u}"

def js_string(json_obj):
    """
    Prepare a JSON string to be embedded inside single-quoted JS:
       JSON.parse('...here...')
    """
    s = json.dumps(json_obj, separators=(",", ":"))
    # escape backslashes and single quotes so the JS string stays valid
    s = s.replace("\\", "\\\\").replace("'", "\\'")
    return s

def exists_in_cwd(fname):
    return os.path.exists(fname)

# -------------------- Load inputs --------------------
summary = pd.read_csv(summary_csv)
snps_ivw = pd.read_csv(snps_csv)

# Try harmonised_data for richer tables/plots
harm_path = "harmonised_data.csv"
harm = pd.read_csv(harm_path) if os.path.exists(harm_path) else None

# -------------------- Extract KPI/MR metrics --------------------
row = summary.iloc[0] if len(summary) else pd.Series({})

N_SNPs   = int(row.get("N_SNPs", len(snps_ivw))) if not pd.isna(row.get("N_SNPs", float("nan"))) else len(snps_ivw)
QC_pass  = N_SNPs  # adjust if you later track QC separately

IVW_OR   = fmt_num(row.get("IVW_OR"))
IVW_L    = fmt_num(row.get("IVW_CI_Lower"))
IVW_U    = fmt_num(row.get("IVW_CI_Upper"))
IVW_P    = fmt_p(row.get("IVW_Pval"))
IVW_QP   = fmt_p(row.get("IVW_Het_Q_Pval"))
IVW_I2   = fmt_num(row.get("I2_Statistic"), r=2)

WME_OR   = fmt_num(row.get("WM_OR"))
WME_L    = fmt_num(row.get("WM_CI_Lower"))
WME_U    = fmt_num(row.get("WM_CI_Upper"))
WME_P    = fmt_p(row.get("WM_Pval"))

EG_OR    = fmt_num(row.get("Egger_OR"))
EG_L     = fmt_num(row.get("Egger_CI_Lower"))
EG_U     = fmt_num(row.get("Egger_CI_Upper"))
EG_P     = fmt_p(row.get("Egger_Pval"))
EG_INT_P = fmt_p(row.get("Egger_Intercept_Pval"))

# -------------------- SNP table rows --------------------
snp_rows_html = []
if harm is not None and not harm.empty:
    required = {
        "SNP",
        "beta.exposure","se.exposure","pval.exposure",
        "beta.outcome","se.outcome","pval.outcome",
        "effect_allele.exposure","other_allele.exposure"
    }
    if required.issubset(harm.columns):
        for _, r in harm.iterrows():
            snp_rows_html.append(
                "<tr>"
                f"<td>{r['SNP']}</td>"
                f"<td>{fmt_num(r['beta.exposure'],6)}</td>"
                f"<td>{fmt_num(r['se.exposure'],6)}</td>"
                f"<td>{fmt_p(r['pval.exposure'])}</td>"
                f"<td>{fmt_num(r['beta.outcome'],6)}</td>"
                f"<td>{fmt_num(r['se.outcome'],6)}</td>"
                f"<td>{fmt_p(r['pval.outcome'])}</td>"
                f"<td>{r['effect_allele.exposure']}</td>"
                f"<td>{r['other_allele.exposure']}</td>"
                "</tr>"
            )
    else:
        # fallback to the IVW per-SNP OR table if we don't have harmonised fields
        lower = {c.lower(): c for c in snps_ivw.columns}
        c_snp  = lower.get("snp")
        c_or   = lower.get("ivw_or")
        c_l95  = lower.get("ivw_lower_95")
        c_u95  = lower.get("ivw_upper_95")
        if all([c_snp, c_or, c_l95, c_u95]):
            for _, r in snps_ivw.iterrows():
                snp_rows_html.append(
                    "<tr>"
                    f"<td>{r[c_snp]}</td>"
                    f"<td colspan='2'>{fmt_num(r[c_or],4)}</td>"
                    f"<td colspan='6'>{fmt_num(r[c_l95],4)}–{fmt_num(r[c_u95],4)}</td>"
                    "</tr>"
                )
        else:
            snp_rows_html.append(
                "<tr><td colspan='9' style='text-align:center;color:#a00'>"
                "SNP table could not be built: missing expected columns."
                "</td></tr>"
            )
else:
    # no harmonised file — same fallback
    lower = {c.lower(): c for c in snps_ivw.columns}
    c_snp  = lower.get("snp")
    c_or   = lower.get("ivw_or")
    c_l95  = lower.get("ivw_lower_95")
    c_u95  = lower.get("ivw_upper_95")
    if all([c_snp, c_or, c_l95, c_u95]):
        for _, r in snps_ivw.iterrows():
            snp_rows_html.append(
                "<tr>"
                f"<td>{r[c_snp]}</td>"
                f"<td colspan='2'>{fmt_num(r[c_or],4)}</td>"
                f"<td colspan='6'>{fmt_num(r[c_l95],4)}–{fmt_num(r[c_u95],4)}</td>"
                "</tr>"
            )
    else:
        snp_rows_html.append(
            "<tr><td colspan='9' style='text-align:center;color:#a00'>"
            "SNP table could not be built: MR_IVW_OR_Per_SNP.csv missing expected columns."
            "</td></tr>"
        )

# -------------------- Method summary table --------------------
method_rows = [
    f"<tr><td>IVW</td><td>{IVW_OR}</td><td>—</td><td>{IVW_P}</td><td>{IVW_L}</td><td>{IVW_U}</td><td>{N_SNPs}</td></tr>",
    f"<tr><td>Weighted median</td><td>{WME_OR}</td><td>—</td><td>{WME_P}</td><td>{WME_L}</td><td>{WME_U}</td><td>{N_SNPs}</td></tr>",
    f"<tr><td>MR Egger</td><td>{EG_OR}</td><td>—</td><td>{EG_P}</td><td>{EG_L}</td><td>{EG_U}</td><td>{N_SNPs}</td></tr>",
]

# -------------------- Plotly JSON for scatter & forest --------------------
# Scatter (beta.exposure vs beta.outcome)
if harm is not None and {"beta.exposure","beta.outcome","SNP"}.issubset(harm.columns):
    scatter_json = {
        "data": [{
            "x": harm["beta.exposure"].tolist(),
            "y": harm["beta.outcome"].tolist(),
            "text": harm["SNP"].astype(str).tolist(),
            "mode": "markers",
            "type": "scatter",
            "marker": {"size": 8}
        }],
        "layout": {"margin": {"t": 30}, "xaxis": {"title": "β (exposure)"}, "yaxis": {"title": "β (outcome)"}}
    }
else:
    scatter_json = {"data": [], "layout": {"margin": {"t": 30}}}

# Forest from IVW per SNP (OR with error bars to upper)
if {"SNP","IVW_OR","IVW_Lower_95","IVW_Upper_95"}.issubset(snps_ivw.columns):
    x_vals = snps_ivw["IVW_OR"].astype(float).tolist()
    y_vals = snps_ivw["SNP"].astype(str).tolist()
    err_plus = [(float(u) - float(o)) for o, u in zip(snps_ivw["IVW_OR"], snps_ivw["IVW_Upper_95"])]
    forest_json = {
        "data": [{
            "type": "bar",
            "x": x_vals,
            "y": y_vals,
            "orientation": "h",
            "error_x": {"type": "data", "array": err_plus, "visible": True}
        }],
        "layout": {"margin": {"t": 30}, "xaxis": {"title": "IVW OR (95% CI)"}, "yaxis": {"title": "SNP"}}
    }
else:
    forest_json = {"data": [], "layout": {"margin": {"t": 30}}}

# -------------------- Image src placeholders --------------------
# These are *paths only* (the HTML will inject <img> tags via JS).
man_exp_src = "exposure_manhattan.png" if exists_in_cwd("exposure_manhattan.png") else ""
man_out_src = "outcome_manhattan.png"  if exists_in_cwd("outcome_manhattan.png")  else ""
qq_exp_src  = "exposure_qq.png"        if exists_in_cwd("exposure_qq.png")        else ""
qq_out_src  = "outcome_qq.png"         if exists_in_cwd("outcome_qq.png")         else ""

# -------------------- Load HTML template & replace --------------------
with open(html_template, "r", encoding="utf-8") as f:
    html = f.read()

replacements = {
    "<!-- SNP_COUNT -->":            str(N_SNPs),
    "<!-- QC_COUNT -->":             str(QC_pass),
    "<!-- IVW_OR -->":               IVW_OR,
    "<!-- IVW_CI -->":               CI(IVW_L, IVW_U),
    "<!-- IVW_P -->":                IVW_P,
    "<!-- IVW_Q -->":                "—",           # Cochran's Q (if you add it later, wire here)
    "<!-- IVW_QP -->":               IVW_QP,
    "<!-- IVW_I2 -->":               IVW_I2,
    "<!-- WME_OR -->":               WME_OR,
    "<!-- WME_CI -->":               CI(WME_L, WME_U),
    "<!-- WME_P -->":                WME_P,
    "<!-- WME_I2 -->":               "—",
    "<!-- EGGER_OR -->":             EG_OR,
    "<!-- EGGER_CI -->":             CI(EG_L, EG_U),
    "<!-- EGGER_P -->":              EG_P,
    "<!-- EGGER_Q -->":              "—",
    "<!-- EGGER_QP -->":             "—",
    "<!-- EGGER_I2 -->":             "—",
    "<!-- EGGER_INT_P -->":          EG_INT_P,
    "<!-- SNP_TABLE_ROWS -->":       "\n".join(snp_rows_html),
    "<!-- METHOD_TABLE_ROWS -->":    "\n".join(method_rows),
    "<!-- SCATTER_JSON -->":         js_string(scatter_json),
    "<!-- FOREST_JSON -->":          js_string(forest_json),
    "<!-- MANHATTAN_EXPOSURE_SRC -->": man_exp_src or "",
    "<!-- MANHATTAN_OUTCOME_SRC -->":  man_out_src or "",
    "<!-- QQ_EXPOSURE_SRC -->":        qq_exp_src or "",
    "<!-- QQ_OUTCOME_SRC -->":         qq_out_src or "",
}

for k, v in replacements.items():
    html = html.replace(k, v)

with open(output_html, "w", encoding="utf-8") as f:
    f.write(html)

print(f"✅ Report created: {output_html}")
