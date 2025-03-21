#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --- Parameters ---
params.output_dir = "./results"

// --- Processes ---

process simulate_data {
    publishDir "${params.output_dir}", mode: 'copy'

    output:
    path("exposure_gwas.csv"), emit: exposure
    path("outcome_gwas.csv"), emit: outcome

    script:
    """
    python3 01_generate_simulated_data.py exposure_gwas.csv outcome_gwas.csv
    """
}

process exploratory_analysis {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path exposure
    path outcome

    output:
    path("exposure_manhattan.png")
    path("outcome_manhattan.png")

    script:
    """
    python3 02_exploratory_analysis.py ${exposure} ${outcome} ${params.output_dir}
    """
}

process gwas_processing {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path exposure
    path outcome

    output:
    path("filtered_SNPs_for_MR.csv"), emit: filtered

    script:
    """
    python3 03_gwas_processing.py ${exposure} ${outcome} filtered_SNPs_for_MR.csv
    """
}

process mr_analysis {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path filtered

    output:
    path("MR_Formatted_Results.csv"), emit: summary
    path("MR_IVW_OR_Per_SNP.csv"), emit: snps

    script:
    """
    Rscript 04_mr_analyses.R ${filtered}
    """
}

process visualisation {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path summary
    path snps

    output:
    path("mr_summary_estimates.png")
    path("ivw_per_snp_forest_plot.png")
    path("ivw_all_snp_ORs.csv")

    script:
    """
    python3 05_visualisations.py ${summary} ${snps} ${params.output_dir}
    """
}

// --- Workflow block that defines process order ---
workflow {
    // Step 1: Simulate GWAS data
    simulate_data()

    // Step 2: Exploratory analysis
    exploratory_analysis(simulate_data.out.exposure, simulate_data.out.outcome)

    // Step 3: Harmonize and filter SNPs
    gwas_processing(simulate_data.out.exposure, simulate_data.out.outcome)

    // Step 4: Run MR analysis
    mr_analysis(gwas_processing.out.filtered)

    // Step 5: Plot MR results
    visualisation(mr_analysis.out.summary, mr_analysis.out.snps)
}
