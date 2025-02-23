-- Load CSVs directly -- 
-- exposure and outcome tables

SELECT * FROM exposure;

SELECT 
    COUNT(*) AS total_variants,
    COUNT(CASE WHEN LENGTH(ALT) > LENGTH(REF) THEN 1 END) AS insertions,
    COUNT(CASE WHEN LENGTH(ALT) < LENGTH(REF) THEN 1 END) AS deletions
FROM exposure;

-- 522 variants - 51 INDELs
-- Clean 51 INDELs in exposure VCF

SELECT *
FROM exposure
WHERE LENGTH(ALT) = LENGTH(REF)
ORDER BY ALT;

-- Let´s count the number of SNPs now  

SELECT COUNT(*)
FROM exposure
WHERE LENGTH(ALT) = LENGTH(REF)
ORDER BY ALT;

-- post-cleaning -> 471 SNPs vs 522 initial (-51 INDELs)

-- Let´s do the same for the outcome VCF -- 
-- outcome.csv

SELECT * FROM outcome;

SELECT COUNT(*) FROM outcome; -- 513 SNPs

SELECT 
    COUNT(*) AS total_variants,
    COUNT(CASE WHEN LENGTH(ALT) > LENGTH(REF) THEN 1 END) AS insertions,
    COUNT(CASE WHEN LENGTH(ALT) < LENGTH(REF) THEN 1 END) AS deletions
FROM outcome; -- 49 INDELs

-- Now let´s clean it! -- 
SELECT *
FROM outcome
WHERE LENGTH(ALT) = LENGTH(REF)
ORDER BY ALT;

-- Check -- 
SELECT COUNT(*)
FROM outcome
WHERE LENGTH(ALT) = LENGTH(REF)
ORDER BY ALT; -- 464 SNPs from initial 513 (-49 SNPs)