#!/usr/bin/env python3
import pandas as pd
import sys
import json
from datetime import date

summary_csv, snps_csv, html_template, output_html = sys.argv[1:]

summary_df = pd.read_csv(summary_csv)
snps_df = pd.read_csv(snps_csv)

# Dashboard metrics
n_instruments = len(snps_df)
n_qc_passed = snps_df.query("P_Exposure < 5e-8").shape[0]

ivw = summary_df.query("Method == 'IVW'").iloc[0]
ivw_or = round(ivw['OR'], 2)
ivw_ci = f"{round(ivw['CI_lower'], 2)}–{round(ivw['CI_upper'], 2)}"
ivw_p = f"{ivw['P-value']:.4g}"

# SNP table HTML
snp_rows = "\n".join(
    f"<tr><td>{r.SNP}</td><td>{r.Beta_Exposure}</td><td>{r.P_Exposure}</td><td>{r.Beta_Outcome}</td><td>{r.P_Outcome}</td></tr>"
    for r in snps_df.itertuples()
)

# Method table HTML
method_rows = "\n".join(
    f"<tr><td>{r.Method}</td><td>{round(r.Estimate, 3)}</td><td>{r._5:.4g}</td><td>{round(r.CI_lower, 3)}</td><td>{round(r.CI_upper, 3)}</td></tr>"
    for r in summary_df.itertuples()
)

# Plot data (scatter)
scatter_json = json.dumps([{
    "x": snps_df["Beta_Exposure"].tolist(),
    "y": snps_df["Beta_Outcome"].tolist(),
    "text": snps_df["SNP"].tolist(),
    "mode": "markers+text",
    "type": "scatter",
    "textposition": "top center",
    "marker": { "size": 10, "color": "#007aff" }
}], indent=2)

# Plot data (forest)
forest_json = json.dumps([{
    "x": snps_df["Beta_Exposure"].tolist(),
    "y": snps_df["SNP"].tolist(),
    "error_x": {
        "type": "data",
        "array": [0.05] * len(snps_df),  # optional: replace with real SEs
        "visible": True
    },
    "type": "bar",
    "orientation": "h",
    "marker": { "color": "#2ca02c" }
}], indent=2)

# Load and populate template
with open(html_template) as f:
    html = f.read()

html = html.replace("<!-- SNP_COUNT -->", str(n_instruments))
html = html.replace("<!-- QC_COUNT -->", str(n_qc_passed))
html = html.replace("<!-- IVW_OR -->", f"{ivw_or} [{ivw_ci}]")
html = html.replace("<!-- IVW_P -->", ivw_p)
html = html.replace("<!-- SNP_TABLE_ROWS -->", snp_rows)
html = html.replace("<!-- METHOD_TABLE_ROWS -->", method_rows)
html = html.replace("<!-- SCATTER_JSON -->", scatter_json)
html = html.replace("<!-- FOREST_JSON -->", forest_json)

# Write report
with open(output_html, "w") as f:
    f.write(html)

print(f"✅ Report created: {output_html}")
