nextflow.enable.dsl=2

workflow {
    
    // Step 1: Generate GWAS Data
    generate_gwas_data()

    // Step 2: Exploratory Analysis of GWAS Data
    exploratory_analysis()

    // Step 3: Process GWAS Data (QC, filtering, formatting)
    process_gwas_data()

    // Step 4: Perform MR Analyses (TwoSampleMR in R)
    mr_analysis()

    // Step 5: Generate Visualizations
    visualize_results()
}

// Define each process

process generate_gwas_data {
    tag 'Generating GWAS data'
    input:
    output:
    path "ldl_gwas.csv"
    path "ad_gwas.csv"

    script:
    """
    python 01_generate_simulated_data.py
    """
}

process exploratory_analysis {
    tag 'Exploratory Data Analysis'
    input:
    path "ldl_gwas.csv"
    path "ad_gwas.csv"

    script:
    """
    jupyter nbconvert --to notebook --execute 02_exploratory_analysis.ipynb
    """
}

process process_gwas_data {
    tag 'Processing GWAS Data'
    input:
    path "ldl_gwas.csv"
    path "ad_gwas.csv"
    output:
    path "filtered_SNPs_for_MR.csv"

    script:
    """
    python 03_gwas_processing.py
    """
}

process mr_analysis {
    tag 'Performing MR Analysis'
    input:
    path "filtered_SNPs_for_MR.csv"
    output:
    path "MR_Formatted_Results.csv"

    script:
    """
    Rscript 04_mr_analyses.R
    """
}

process visualize_results {
    tag 'Generating MR Visualizations'
    input:
    path "MR_Formatted_Results.csv"

    script:
    """
    python 05_visualisations.py
    """
}

