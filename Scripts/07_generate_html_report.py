#!/usr/bin/env python3
import os
import sys
import math
import pandas as pd

# ------------------------------------------------------------------
# Usage:
#   07_generate_html_report.py <MR_Formatted_Results.csv>
#                               <MR_IVW_OR_Per_SNP.csv>
#                               <template.html>
#                               <out.html>
#
# Assumes: harmonised_data.csv is in CWD (staged by Nextflow).
# ------------------------------------------------------------------

if len(sys.argv) != 5:
    sys.stderr.write(
        "Usage: 07_generate_html_report.py "
        "<MR_Formatted_Results.csv> <MR_IVW_OR_Per_SNP.csv> "
        "<template.html> <out.html>\n"
    )
    sys.exit(1)

summary_csv, snps_csv, html_template, output_html = sys.argv[1:]

# ----------------------------- IO helpers -----------------------------
def read_csv_safe(path):
    try:
        if path and os.path.exists(path) and os.path.getsize(path) > 0:
            return pd.read_csv(path)
    except Exception:
        pass
    return pd.DataFrame()

summary   = read_csv_safe(summary_csv)
snps_ivw  = read_csv_safe(snps_csv)
harm      = read_csv_safe("harmonised_data.csv")

# ----------------------------- formatters -----------------------------
DASH = "—"

def _is_na(x):
    try:
        return pd.isna(x) or (isinstance(x, float) and math.isnan(x))
    except Exception:
        return False

def fmt_float(x, d=3, dash=DASH):
    try:
        if _is_na(x): return dash
        v = float(x)
        if _is_na(v): return dash
        return f"{v:.{d}f}"
    except Exception:
        return dash

def fmt_p(x, dash=DASH):
    try:
        if _is_na(x): return dash
        v = float(x)
        if v < 1e-3:
            return f"{v:.3g}"
        return f"{v:.4f}"
    except Exception:
        return dash

def ci_str(l, u, dash=DASH):
    if l in (None, dash, "") or u in (None, dash, ""):
        return dash
    try:
        if _is_na(l) or _is_na(u): return dash
        return f"{float(l):.3f}–{float(u):.3f}"
    except Exception:
        return f"{l}–{u}"

# ----------------------------- column picking -----------------------------
def pick(cols, *cands):
    """Return first existing column (case-sensitive) from candidates."""
    for c in cands:
        if c in cols:
            return c
    return None

def pick_ci(df, *cands):
    """Case-insensitive column picker -> actual column name or None."""
    if df is None or df.empty:
        return None
    low = {c.lower(): c for c in df.columns}
    for cand in cands:
        if cand.lower() in low:
            return low[cand.lower()]
    return None

# ----------------------------- KPIs / method summary -----------------------------
row = summary.iloc[0] if len(summary) else pd.Series({})

def get_val(df, row, *aliases):
    c = pick_ci(df, *aliases)
    return (row.get(c) if c else None) if not df.empty else None

# N SNPs
N_SNPs = None
c_n = pick_ci(summary, "N_SNPs", "N", "Num_SNPs")
if c_n:
    try: N_SNPs = int(row.get(c_n))
    except Exception: N_SNPs = None
if N_SNPs is None:
    if not snps_ivw.empty and pick_ci(snps_ivw, "SNP"): N_SNPs = snps_ivw.shape[0]
    elif not harm.empty and pick_ci(harm, "SNP", "rsid", "ID"): N_SNPs = harm.shape[0]
    else: N_SNPs = 0
QC_pass = N_SNPs

# IVW fields
IVW_OR   = get_val(summary, row, "IVW_OR", "IVW OR", "IVW_OR_Exp")
IVW_L    = get_val(summary, row, "IVW_CI_Lower", "IVW_CI_L", "IVW_Lower_95")
IVW_U    = get_val(summary, row, "IVW_CI_Upper", "IVW_CI_U", "IVW_Upper_95")
IVW_P    = get_val(summary, row, "IVW_Pval", "IVW_P", "IVW p", "p_IVW")

# Heterogeneity aliases
IVW_Q    = get_val(summary, row, "IVW_Het_Q", "Q", "Q_IVW", "Cochran_Q", "Cochran_Q_IVW")
IVW_QP   = get_val(summary, row, "IVW_Het_Q_Pval", "Q_pval", "Q_Pval", "Q p-value", "Q_P_IVW")
IVW_I2   = get_val(summary, row, "I2_Statistic", "I2", "I^2")

def compute_i2_from_q(Q, k):
    try:
        q = float(Q)
        df = max(1, int(k) - 1)
        if q <= 0: return None
        i2 = max(0.0, (q - df) / q) * 100.0
        return i2
    except Exception:
        return None

if IVW_I2 is None and IVW_Q is not None and N_SNPs and N_SNPs >= 2:
    IVW_I2 = compute_i2_from_q(IVW_Q, N_SNPs)

# Weighted median
WME_OR   = get_val(summary, row, "WM_OR", "WME_OR", "Weighted_median_OR")
WME_L    = get_val(summary, row, "WM_CI_Lower", "WME_CI_Lower")
WME_U    = get_val(summary, row, "WM_CI_Upper", "WME_CI_Upper")
WME_P    = get_val(summary, row, "WM_Pval", "WME_Pval", "WM_P")

# Egger
EG_OR    = get_val(summary, row, "Egger_OR", "MR_Egger_OR")
EG_L     = get_val(summary, row, "Egger_CI_Lower", "MR_Egger_CI_Lower")
EG_U     = get_val(summary, row, "Egger_CI_Upper", "MR_Egger_CI_Upper")
EG_P     = get_val(summary, row, "Egger_Pval", "MR_Egger_Pval")
EG_INT_P = get_val(summary, row, "Egger_Intercept_Pval", "MR_Egger_Intercept_Pval", "Egger_intercept_p")

EG_Q     = get_val(summary, row, "Egger_Q", "Q_Egger", "Q'", "Qprime", "Cochran_Q_Egger")
EG_QP    = get_val(summary, row, "Egger_Q_pval", "Qp_Egger", "Q'_pval", "Qprime_pval", "Egger_QP")
EG_I2    = get_val(summary, row, "Egger_I2", "I2_Egger", "I^2_Egger")
if EG_I2 is None and EG_Q is not None and N_SNPs and N_SNPs >= 2:
    EG_I2 = compute_i2_from_q(EG_Q, N_SNPs)

# Format for display
IVW_OR_f = fmt_float(IVW_OR, 3)
IVW_CI_f = ci_str(IVW_L, IVW_U)
IVW_P_f  = fmt_p(IVW_P)
IVW_Q_f  = fmt_float(IVW_Q, 3) if IVW_Q is not None else DASH
IVW_QP_f = fmt_p(IVW_QP)
IVW_I2_f = fmt_float(IVW_I2, 1) if IVW_I2 is not None else DASH

WME_OR_f = fmt_float(WME_OR, 3)
WME_CI_f = ci_str(WME_L, WME_U)
WME_P_f  = fmt_p(WME_P)

EG_OR_f  = fmt_float(EG_OR, 3)
EG_CI_f  = ci_str(EG_L, EG_U)
EG_P_f   = fmt_p(EG_P)
EG_Q_f   = fmt_float(EG_Q, 3) if EG_Q is not None else DASH
EG_QP_f  = fmt_p(EG_QP)
EG_I2_f  = fmt_float(EG_I2, 1) if EG_I2 is not None else DASH
EG_INT_P_f = fmt_p(EG_INT_P)

# ----------------------------- SNP table -----------------------------
snp_rows_html = []

def snp_rows_from_harmonised(harm_df, ivw_df):
    rows = []
    if harm_df is None or harm_df.empty:
        return rows

    cols = list(harm_df.columns)
    def col(*names): return pick(cols, *names)

    c_id  = col("SNP") or col("rsid","ID")
    c_chr = col("chr.exposure") or col("chr.outcome","CHR","chr","chrom","CHROM")
    c_bp  = col("pos.exposure") or col("pos.outcome","BP","bp","position","POS")

    c_ea_x = col("effect_allele.exposure") or col("EA.exposure","EA_Exposure","effect_allele_exposure")
    c_oa_x = col("other_allele.exposure")  or col("OA.exposure","OA_Exposure","other_allele_exposure")
    c_ea_y = col("effect_allele.outcome")  or col("EA.outcome","EA_Outcome","effect_allele_outcome")
    c_oa_y = col("other_allele.outcome")   or col("OA.outcome","OA_Outcome","other_allele_outcome")

    c_bx = col("beta.exposure") or col("BETA.exposure","beta_exposure","BETA_Exposure","b_exp","b.exposure")
    c_sx = col("se.exposure")   or col("SE.exposure","se_exposure","SE_Exposure","se_exp","se.exposure")
    c_px = col("pval.exposure") or col("P.exposure","p.exposure","pval_exposure","P_Exposure","p_exp","p.exposure")

    c_by = col("beta.outcome")  or col("BETA.outcome","beta_outcome","BETA_Outcome","b_out","b.outcome")
    c_sy = col("se.outcome")    or col("SE.outcome","se_outcome","SE_Outcome","se_out","se.outcome")
    c_py = col("pval.outcome")  or col("P.outcome","p.outcome","pval_outcome","P_Outcome","p_out","p.outcome")

    c_eafx = col("eaf.exposure") or col("EAF.exposure","eaf_exposure","EAF_Exposure")
    c_eafy = col("eaf.outcome")  or col("EAF.outcome","eaf_outcome","EAF_Outcome")

    # per-SNP IVW OR map
    or_map = {}
    if ivw_df is not None and not ivw_df.empty:
        c_snp_ivw = pick_ci(ivw_df, "SNP")
        c_or      = pick_ci(ivw_df, "IVW_OR", "OR", "IVW OR")
        c_l95     = pick_ci(ivw_df, "IVW_Lower_95", "Lower_95", "CI_Lower")
        c_u95     = pick_ci(ivw_df, "IVW_Upper_95", "Upper_95", "CI_Upper")
        if c_snp_ivw and c_or and c_l95 and c_u95:
            for _, rr in ivw_df[[c_snp_ivw, c_or, c_l95, c_u95]].dropna().iterrows():
                or_map[str(rr[c_snp_ivw])] = (
                    fmt_float(rr[c_or], 4),
                    ci_str(rr[c_l95], rr[c_u95])
                )

    if not (c_id and c_bx and c_by):
        rows.append(
            "<tr><td colspan='17' style='text-align:center;color:#a00'>"
            "harmonised_data.csv is present but missing required columns to build the SNP table."
            "</td></tr>"
        )
        return rows

    for _, r in harm_df.iterrows():
        snp = str(r[c_id])
        CHR = r[c_chr] if c_chr else DASH
        BP  = r[c_bp]  if c_bp  else DASH

        EAx = r[c_ea_x] if c_ea_x else DASH
        OAx = r[c_oa_x] if c_oa_x else DASH
        EAy = r[c_ea_y] if c_ea_y else DASH
        OAy = r[c_oa_y] if c_oa_y else DASH

        bx  = fmt_float(r[c_bx], 6) if c_bx else DASH
        sx  = fmt_float(r[c_sx], 6) if c_sx else DASH
        px  = fmt_p(r[c_px])        if c_px else DASH

        by  = fmt_float(r[c_by], 6) if c_by else DASH
        sy  = fmt_float(r[c_sy], 6) if c_sy else DASH
        py  = fmt_p(r[c_py])        if c_py else DASH

        eafx = fmt_float(r[c_eafx], 4) if c_eafx else DASH
        eafy = fmt_float(r[c_eafy], 4) if c_eafy else DASH

        or_val, or_ci = or_map.get(snp, (DASH, DASH))

        rows.append(
            "<tr>"
            f"<td>{snp}</td>"
            f"<td>{CHR}</td><td>{BP}</td>"
            f"<td>{EAx}</td><td>{OAx}</td><td>{EAy}</td><td>{OAy}</td>"
            f"<td>{bx}</td><td>{sx}</td><td>{px}</td>"
            f"<td>{by}</td><td>{sy}</td><td>{py}</td>"
            f"<td>{eafx}</td><td>{eafy}</td>"
            f"<td>{or_val}</td><td>{or_ci}</td>"
            "</tr>"
        )
    return rows

if not harm.empty:
    snp_rows_html = snp_rows_from_harmonised(harm, snps_ivw)
else:
    snp_rows_html = [
        "<tr><td colspan='17' style='text-align:center;color:#a00'>"
        "harmonised_data.csv not found — SNP table not available."
        "</td></tr>"
    ]

# ----------------------------- Method table -----------------------------
method_rows = [
    f"<tr><td>IVW</td><td>{IVW_OR_f}</td><td>{IVW_CI_f}</td><td>{IVW_P_f}</td>"
    f"<td>{IVW_Q_f}</td><td>{IVW_QP_f}</td><td>{IVW_I2_f}</td><td>{DASH}</td></tr>",
    f"<tr><td>Weighted median</td><td>{WME_OR_f}</td><td>{WME_CI_f}</td><td>{WME_P_f}</td>"
    f"<td>{DASH}</td><td>{DASH}</td><td>{DASH}</td><td>{DASH}</td></tr>",
    f"<tr><td>MR Egger</td><td>{EG_OR_f}</td><td>{EG_CI_f}</td><td>{EG_P_f}</td>"
    f"<td>{EG_Q_f}</td><td>{EG_QP_f}</td><td>{EG_I2_f}</td><td>{EG_INT_P_f}</td></tr>",
]

# ----------------------------- Template injection -----------------------------
with open(html_template, "r", encoding="utf-8") as f:
    html = f.read()

repls = {
    "<!-- SNP_COUNT -->":         str(N_SNPs),
    "<!-- QC_COUNT -->":          str(QC_pass),

    "<!-- IVW_OR -->":            IVW_OR_f,
    "<!-- IVW_CI -->":            IVW_CI_f,
    "<!-- IVW_P -->":             IVW_P_f,
    "<!-- IVW_Q -->":             IVW_Q_f,
    "<!-- IVW_QP -->":            IVW_QP_f,
    "<!-- IVW_I2 -->":            IVW_I2_f,

    "<!-- WME_OR -->":            WME_OR_f,
    "<!-- WME_CI -->":            WME_CI_f,
    "<!-- WME_P -->":             WME_P_f,
    "<!-- WME_I2 -->":            "—",

    "<!-- EGGER_OR -->":          EG_OR_f,
    "<!-- EGGER_CI -->":          EG_CI_f,
    "<!-- EGGER_P -->":           EG_P_f,
    "<!-- EGGER_Q -->":           EG_Q_f,
    "<!-- EGGER_QP -->":          EG_QP_f,
    "<!-- EGGER_I2 -->":          EG_I2_f,
    "<!-- EGGER_INT_P -->":       EG_INT_P_f,

    "<!-- SNP_TABLE_ROWS -->":    "\n".join(snp_rows_html),
    "<!-- METHOD_TABLE_ROWS -->": "\n".join(method_rows),

    "<!-- SCATTER_COMBINED_SRC -->":   "MR_Scatter_Combined.png",
    "<!-- LOO_SRC -->":                "MR_LeaveOneOut.png",
    "<!-- SUMMARY_EST_SRC -->":        "mr_summary_estimates.png",
    "<!-- FOREST_STATIC_SRC -->":      "ivw_per_snp_forest_plot.png",
    "<!-- MANHATTAN_EXPOSURE_SRC -->": "exposure_manhattan.png",
    "<!-- MANHATTAN_OUTCOME_SRC -->":  "outcome_manhattan.png",
    "<!-- QQ_EXPOSURE_SRC -->":        "exposure_qq.png",
    "<!-- QQ_OUTCOME_SRC -->":         "outcome_qq.png",
}
for k, v in repls.items():
    html = html.replace(k, v)

with open(output_html, "w", encoding="utf-8") as f:
    f.write(html)

print(f"✅ Report created: {output_html}")
