#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define input/output channels
workflow {

    // --- PARAMETERS ---
    params.output_dir = "./results"

    // --- STEP 1: Simulate data ---
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

    // --- STEP 2: Exploratory plots ---
    process exploratory_analysis {
        publishDir "${params.output_dir}", mode: 'copy'

        input:
        path exposure from simulate_data.out.exposure
        path outcome from simulate_data.out.outcome

        output:
        path("exposure_manhattan.png")
        path("outcome_manhattan.png")

        script:
        """
        python3 02_exploratory_analysis.py ${exposure} ${outcome} ${params.output_dir}
        """
    }

    // --- STEP 3: Filter & Harmonize SNPs ---
    process gwas_processing {
        publishDir "${params.output_dir}", mode: 'copy'

        input:
        path exposure from simulate_data.out.exposure
        path outcome from simulate_data.out.outcome

        output:
        path("filtered_SNPs_for_MR.csv"), emit: filtered

        script:
        """
        python3 03_gwas_processing.py ${exposure} ${outcome} filtered_SNPs_for_MR.csv
        """
    }

    // --- STEP 4: MR Analysis (R) ---
    process mr_analysis {
        publishDir "${params.output_dir}", mode: 'copy'

        input:
        path filtered from gwas_processing.out.filtered

        output:
        path("MR_Formatted_Results.csv"), emit: summary
        path("MR_IVW_OR_Per_SNP.csv"), emit: snps

        script:
        """
        Rscript 04_mr_analyses.R ${filtered}
        """
    }

    // --- STEP 5: Visualise MR results ---
    process visualisation {
        publishDir "${params.output_dir}", mode: 'copy'

        input:
        path summary from mr_analysis.out.summary
        path snps from mr_analysis.out.snps

        output:
        path("mr_summary_estimates.png")
        path("ivw_per_snp_forest_plot.png")
        path("ivw_all_snp_ORs.csv")

        script:
        """
        python3 05_visualisations.py ${summary} ${snps} ${params.output_dir}
        """
    }
}
