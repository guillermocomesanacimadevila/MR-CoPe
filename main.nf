#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ========================== PARAMETERS ==========================
params.output_dir      = "./results"
params.script_dir      = "./"
params.log10_flag      = "n"                // Are input p-values already -log10(p)?
params.clump_kb        = 10000              // LD window in kb
params.clump_r2        = 0.001              // LD clumping threshold
params.trait_keywords  = null               // Comma-separated keywords to retain (e.g., "smoking, alcohol")
params.ld_pop          = "EUR"              // LD reference population
params.exposure        = null
params.outcome         = null

// ========================== SCRIPTS ==========================
workflow {

    def script_explore     = file("${params.script_dir}/01_exploratory_analysis.py")
    def script_process     = file("${params.script_dir}/02_gwas_preprocessing.py")
    def script_ld          = file("${params.script_dir}/03_linkage_disequillibrium.R")
    def script_filter_conf = file("${params.script_dir}/04a_filter_confounders.R")
    def script_mr          = file("${params.script_dir}/04_mr_analyses.R")
    def script_vis_sum     = file("${params.script_dir}/05_visualisations.py")
    def script_vis_scatt   = file("${params.script_dir}/06_visualisations.R")
    def script_html_report = file("${params.script_dir}/07_generate_html_report.py")
    def html_template      = file("Frontend/style.html")

    def exposure_file = file(params.exposure)
    def outcome_file  = file(params.outcome)

    // Step 1: Exploratory analysis
    exploratory_analysis(script_explore, exposure_file, outcome_file)

    // Step 2: Preprocessing
    gwas_processing(script_process, exposure_file, outcome_file)

    // Step 3: LD pruning
    ld_filtering(script_ld, gwas_processing.out.filtered)

    // Step 4: Conditional confounder filtering
    def confounder_output

    if (params.trait_keywords == null || params.trait_keywords == ".") {
        skip_confounder_filtering(ld_filtering.out.pruned)
        confounder_output = skip_confounder_filtering.out.filtered
    } else {
        confounder_filtering(script_filter_conf, ld_filtering.out.pruned)
        confounder_output = confounder_filtering.out.filtered
    }

    // Step 5: MR analysis
    mr_analysis(script_mr, confounder_output)

    // Step 6a: Summary visualisation
    visualisation_summary(script_vis_sum, mr_analysis.out.summary, mr_analysis.out.snps)

    // Step 6b: Scatter plots
    visualisation_scatter(script_vis_scatt, mr_analysis.out.harmonised)

    // Step 7: HTML report
    html_report(script_html_report, mr_analysis.out.summary, mr_analysis.out.snps, html_template)
}

// ========================== PROCESSES ==========================

process exploratory_analysis {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path script
    path exposure
    path outcome

    output:
    path("${params.output_dir}/exposure_manhattan.png"), optional: true
    path("${params.output_dir}/outcome_manhattan.png"), optional: true
    path("${params.output_dir}/exposure_qq.png"), optional: true
    path("${params.output_dir}/outcome_qq.png"), optional: true

    script:
    """
    python3 ${script} ${exposure} ${outcome} ${params.output_dir} ${params.log10_flag}
    """
}

process gwas_processing {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path script
    path exposure
    path outcome

    output:
    path("filtered_SNPs_for_MR.csv"), emit: filtered

    script:
    """
    python3 ${script} ${exposure} ${outcome} filtered_SNPs_for_MR.csv
    """
}

process ld_filtering {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path script
    path filtered

    output:
    path("ld_pruned_SNPs.csv"), emit: pruned

    script:
    """
    Rscript ${script} ${filtered} ld_pruned_SNPs.csv ${params.clump_kb} ${params.clump_r2} ${params.ld_pop}
    """
}

process confounder_filtering {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path script
    path pruned

    output:
    path("confounder_filtered_SNPs.csv"), emit: filtered

    script:
    """
    Rscript ${script} ${pruned} confounder_filtered_SNPs.csv "${params.trait_keywords}"
    """
}

process skip_confounder_filtering {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path pruned

    output:
    path("confounder_filtered_SNPs.csv"), emit: filtered

    script:
    """
    cp ${pruned} confounder_filtered_SNPs.csv
    """
}

process mr_analysis {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path script
    path filtered

    output:
    path("MR_Formatted_Results.csv"), emit: summary
    path("MR_IVW_OR_Per_SNP.csv"), emit: snps
    path("harmonised_data.csv"), emit: harmonised

    script:
    """
    Rscript ${script} ${filtered}
    """
}

process visualisation_summary {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path script
    path summary
    path snps

    output:
    path("mr_summary_estimates.png")
    path("ivw_per_snp_forest_plot.png")
    path("ivw_all_snp_ORs.csv")

    script:
    """
    python3 ${script} ${summary} ${snps} ${params.output_dir}
    """
}

process visualisation_scatter {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path script
    path harmonised

    output:
    path("MR_Scatter_IVW.png")
    path("MR_Scatter_Egger.png")
    path("MR_Scatter_WME.png")
    path("MR_LeaveOneOut.png")

    script:
    """
    Rscript ${script} ${harmonised} ${params.output_dir}
    """
}

process html_report {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path script
    path summary
    path snps
    path html_template

    output:
    path("MR_CoPe_Report.html")

    script:
    """
    python3 ${script} ${summary} ${snps} ${html_template} MR_CoPe_Report.html
    """
}
