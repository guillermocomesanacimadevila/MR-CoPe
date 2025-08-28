#!/usr/bin/env python3
import pandas as pd
import sys, os, json
from math import isnan

def num(x, fmt="{:.4g}", na="NA"):
    try:
        if x is None:
            return na
        xv = float(x)
        if xv != xv:  # NaN
            return na
        return fmt.format(xv)
    except Exception:
        return na

def safe_get(df, col, default=None):
    if df is None or df.empty or col not in df.columns:
        return default
    v = df.iloc[0][col]
    return v if pd.notna(v) else default

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 07_generate_html_report.py <MR_Formatted_Results.csv> <MR_IVW_OR_Per_SNP.csv> <style.html> <out.html>")
        sys.exit(1)

    summary_csv, snps_csv, html_template, output_html = sys.argv[1:]

    # Load inputs (robustly)
    summary = pd.read_csv(summary_csv) if os.path.isfile(summary_csv) and os.path.getsize(summary_csv) > 0 else pd.DataFrame()
    snps    = pd.read_csv(snps_csv)    if os.path.isfile(snps_csv)    and os.path.getsize(snps_csv)    > 0 else pd.DataFrame()

    # Counts
    n_instruments = len(snps) if not snps.empty else int(safe_get(summary, "N_SNPs", 0))
    # We don't have per-SNP exposure p-values in current outputs; reuse n_instruments to avoid template breakage
    n_qc_passed = n_instruments

    # Pull summary stats (single row expected)
    ivw_or   = safe_get(summary, "IVW_OR")
    ivw_l    = safe_get(summary, "IVW_CI_Lower")
    ivw_u    = safe_get(summary, "IVW_CI_Upper")
    ivw_p    = safe_get(summary, "IVW_Pval")

    wm_or    = safe_get(summary, "WM_OR")
    wm_l     = safe_get(summary, "WM_CI_Lower")
    wm_u     = safe_get(summary, "WM_CI_Upper")
    wm_p     = safe_get(summary, "WM_Pval")

    eg_or    = safe_get(summary, "Egger_OR")
    eg_l     = safe_get(summary, "Egger_CI_Lower")
    eg_u     = safe_get(summary, "Egger_CI_Upper")
    eg_p     = safe_get(summary, "Egger_Pval")

    # String for IVW OR [CI] used by your template
    ivw_ci_str = f"{num(ivw_l)}–{num(ivw_u)}"
    ivw_or_str = num(ivw_or, "{:.4f}")

    # Build SNP table rows using the columns we actually have
    if not snps.empty:
        # keep only known columns; rename for display
        disp = snps.copy()
        # format numbers
        for col in ("IVW_OR","IVW_Lower_95","IVW_Upper_95"):
            if col in disp.columns:
                disp[col] = disp[col].map(lambda v: num(v, "{:.4f}"))
        snp_rows = "\n".join(
            f"<tr><td>{r.SNP}</td>"
            f"<td>{getattr(r,'IVW_OR','NA')}</td>"
            f"<td>{getattr(r,'IVW_Lower_95','NA')}</td>"
            f"<td>{getattr(r,'IVW_Upper_95','NA')}</td></tr>"
            for r in disp.itertuples(index=False)
        )
    else:
        snp_rows = "<tr><td colspan='4'>No SNP-level results available</td></tr>"

    # Build “method table” rows from the summary row
    method_rows = ""
    if not summary.empty:
        method_rows_list = [
            ("Inverse variance weighted", ivw_or, ivw_p, ivw_l, ivw_u),
            ("Weighted median",          wm_or,  wm_p,  wm_l,  wm_u),
            ("MR Egger",                 eg_or,  eg_p,  eg_l,  eg_u),
        ]
        method_rows = "\n".join(
            f"<tr><td>{meth}</td>"
            f"<td>{num(orv, '{:.4f}')}</td>"
            f"<td>{num(pv,  '{:.4g}')}</td>"
            f"<td>{num(lv,  '{:.4f}')}</td>"
            f"<td>{num(uv,  '{:.4f}')}</td></tr>"
            for (meth, orv, pv, lv, uv) in method_rows_list
        )
    else:
        method_rows = "<tr><td colspan='5'>No summary results available</td></tr>"

    # Minimal JSON payloads so the template JS doesn't crash
    # (we don't have Beta_Exposure/Beta_Outcome here)
    scatter_json = json.dumps([], indent=2)
    forest_json  = json.dumps([], indent=2)

    # Load and populate template
    with open(html_template, "r", encoding="utf-8") as f:
        html = f.read()

    html = html.replace("<!-- SNP_COUNT -->", str(n_instruments))
    html = html.replace("<!-- QC_COUNT -->", str(n_qc_passed))
    html = html.replace("<!-- IVW_OR -->", f"{ivw_or_str} [{ivw_ci_str}]")
    html = html.replace("<!-- IVW_P -->", num(ivw_p, "{:.4g}"))
    html = html.replace("<!-- SNP_TABLE_ROWS -->", snp_rows)
    html = html.replace("<!-- METHOD_TABLE_ROWS -->", method_rows)
    html = html.replace("<!-- SCATTER_JSON -->", scatter_json)
    html = html.replace("<!-- FOREST_JSON -->", forest_json)

    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"✅ Report created: {output_html}")

if __name__ == "__main__":
    main()
