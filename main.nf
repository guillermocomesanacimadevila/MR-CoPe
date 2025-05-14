#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// PARAMETERS
params.output_dir  = "./results"
params.script_dir  = "./"  // Folder where your scripts are located

// SCRIPTS
workflow {

    def script_explore   = file("${params.script_dir}/01_exploratory_analysis.py")
    def script_process   = file("${params.script_dir}/02_gwas_preprocessing.py")
    def script_ld        = file("${params.script_dir}/03_linkage_disequillibrium.R")
    def script_mr        = file("${params.script_dir}/04_mr_analyses.R")
    def script_vis_sum   = file("${params.script_dir}/05_visualisations.py")
    def script_vis_scatt = file("${params.script_dir}/06_visualisations.R")

    def exposure_file = file(params.exposure)
    def outcome_file  = file(params.outcome)

    // Step 1: Exploratory analysis
    exploratory_analysis(script_explore, exposure_file, outcome_file)

    // Step 2: GWAS preprocessing
    gwas_processing(script_process, exposure_file, outcome_file)

    // Step 3: LD pruning
    ld_filtering(script_ld, gwas_processing.out.filtered)

    // Step 4: MR analysis
    mr_analysis(script_mr, ld_filtering.out.pruned)

    // Step 5a: Python visualisation
    visualisation_summary(script_vis_sum, mr_analysis.out.summary, mr_analysis.out.snps)

    // Step 5b: R scatter plots
    visualisation_scatter(script_vis_scatt, mr_analysis.out.harmonised)
}

// =============================== PROCESSES ===============================

process exploratory_analysis {
    publishDir "${params.output_dir}", mode: 'copy'

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
    Rscript ${script} ${filtered} ld_pruned_SNPs.csv
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
    path("harmonised_data.csv"), emit: harmonised

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
    path("MR_LeaveOneOut.png")

    script:
    """
    Rscript ${script} ${harmonised} ${params.output_dir}
    """
}
