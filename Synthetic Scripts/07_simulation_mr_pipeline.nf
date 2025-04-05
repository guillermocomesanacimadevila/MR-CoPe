#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --- Parameters ---
params.output_dir  = "./results"
params.script_dir  = "./"  // Folder where your scripts are located

// --- Workflow ---
workflow {

    // Declare script file paths
    def script_sim     = file("${params.script_dir}/01_generate_simulated_data.py")
    def script_explore = file("${params.script_dir}/02_exploratory_analysis.py")
    def script_process = file("${params.script_dir}/03_gwas_processing.py")
    def script_ld      = file("${params.script_dir}/04_linkage_disequillibrium.R")
    def script_mr      = file("${params.script_dir}/05_mr_analyses.R")
    def script_vis_py  = file("${params.script_dir}/06_visualisations.py")
    def script_vis_r   = file("${params.script_dir}/06_visualisations.R")  // New R-based scatter plot

    // Step 1: Simulate GWAS data
    simulate_data(script_sim)

    // Step 2: Manhattan plots
    exploratory_analysis(script_explore, simulate_data.out.exposure, simulate_data.out.outcome)

    // Step 3: Harmonise & filter GWAS
    gwas_processing(script_process, simulate_data.out.exposure, simulate_data.out.outcome)

    // Step 4: LD pruning
    ld_filtering(script_ld, gwas_processing.out.filtered)

    // Step 5: MR analysis (now also emits harmonised_data.csv)
    mr_analysis(script_mr, ld_filtering.out.pruned)

    // Step 6a: Python visual summary (forest plots etc)
    visualisation_summary(script_vis_py, mr_analysis.out.summary, mr_analysis.out.snps)

    // Step 6b: R-based MR scatter plot
    visualisation_scatter(script_vis_r, mr_analysis.out.harmonised)
}

// --- Processes ---

process simulate_data {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path script

    output:
    path("exposure_gwas.csv"), emit: exposure
    path("outcome_gwas.csv"), emit: outcome

    script:
    """
    python3 ${script} exposure_gwas.csv outcome_gwas.csv
    """
}

process exploratory_analysis {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path script
    path exposure
    path outcome

    output:
    path("exposure_manhattan.png")
    path("outcome_manhattan.png")

    script:
    """
    python3 ${script} ${exposure} ${outcome} ${params.output_dir}
    """
}

process gwas_processing {
    publishDir "${params.output_dir}", mode: 'copy'

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
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path script
    path filtered

    output:
    path("ld_pruned_SNPs.csv"), emit: pruned

    script:
    """
    cp ${filtered} filtered_input.csv
    Rscript ${script} filtered_input.csv ld_pruned_SNPs.csv
    """
}

process mr_analysis {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path script
    path pruned

    output:
    path("MR_Formatted_Results.csv"), emit: summary
    path("MR_IVW_OR_Per_SNP.csv"), emit: snps
    path("harmonised_data.csv"), emit: harmonised  // NEW OUTPUT

    script:
    """
    Rscript ${script} ${pruned}
    """
}

process visualisation_summary {
    publishDir "${params.output_dir}", mode: 'copy'

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
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path script
    path harmonised

    output:
    path("MR_Scatter_IVW.png")
    path("MR_Scatter_Egger.png")
    path("MR_Scatter_WME.png")

    script:
    """
    Rscript ${script} ${harmonised}
    """
}