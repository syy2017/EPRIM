# EPRIM
ERPIM (epigenetic regulator in immunology) is an computational approach for identifying epigenetic regulators (ERs) associated with immune
function. With gene expression profiles of tumor samples, predefined immune gene sets and interested epigenetic regulator genes, EPRIM
can effectively identify the ERs potentially involved in the immune function.

**Generation of gene rank list**

For each ER, we first generated the ranked coding gene list (without ER genes) based on the expression correlation between those coding genes 
and this ER. We calculated the tumor purity of each sample, and then partial correlation coefficient (PCC) correcting for purity. Combining 
partial correlation coefficient (PCC) and statistical significance (P value), a rank score was calculated to rank coding genes.

**Selection of immune ERs**

Then the ranked gene list was subjected to gene set enrichment analysis (GSEA) on each immune-related pathway, which quantifies the expression 
shift of gene sets in the rank list. Before mapping immune signature genes, we first refined the immune gene sets based on the expression correlation
between single gene and all genes in one gene set. Next, we mapped immune genes in the refined gene sets to the rank list and measured whether these
immune genes significantly enriched in the top or bottom of the ranked gene list. Briefly, up to a given position in the rank list, we respectively
calculated the fraction of genes in one gene set ("hits") with their rank scores as weight and the fraction of genes not in this gene set ("misses").
The enrichment score (ES) was defined as the maximum deviation from zero of  the "hits" fractions subtracting the "misses" fractions. P values were
also calculated for each gene set and were adjusted using the false discovery rate (FDR).

## Usage

    EPRIM (exp.profile, exp.profile.file, interested.ER, signature.list, cor.cutoff = NULL, 

      cor.method = "spearman", platform ="illumina", is.adjusted = FALSE, min.sz = 1, 
       
      perm.times = 100)

**Arguments**

***exp.profile*** A numeric matrix containing the expression of mRNA with
rownames and colnames, row as genes, column as samples.

***exp.profile.file*** A character string representing the pathname of mRNA expression 
file, gene symbol as rownames.

***interested.ER*** A character vector of interested epigenetic regulator genes. 

***signature.list*** A list containing immune gene sets. 

***cor.cutoff*** The cutoff of correlation coefficient used for immune signature
refinement, cor.cutoff = NULL by default. 

***cor.method*** The method used for correlation calculation, either "spearman", "pearson"
or "kendall", "spearman" by default. 

***platform*** Character string indicating platform type, either "affymetrix", "agilent"
or "illumina". Defaults to "illumina".  

***is.adjusted*** Whether or not the 'exp.profile' is preprocessed. If the expression values
were log-transformed and the rows containing many 0 values were removed, is.adjusted = TRUE.

***min.sz*** The least immune signature genes matched in the expression profile. Defaults to 
1. 

***perm.times*** Number of permutation times used for enrichment analysis. 

**Value**

A list including detailed correlation results for epigenetic regulator-signature pairs, and refined 
immune gene sets.

**Example**

``` r
workDir <- "D:/Project/IJob"
exp.profile <- get(load(file=paste0(workDir,"/data/PCG_tpm2.RData")))
exp.profile.file <- paste0(workDir,"/data/PCG_tpm2Name.txt")
interested.ER <- get(load(file=paste0(workDir,"/data/ERgenes.RData"))) 
signature.list <- get(load(file=paste0(workDir,"/data/signature.list.RData"))) 

## Gene expression data
exp.profile[1:3,1:3]
#>                 TCGA.CG.4462 TCGA.CG.4306 TCGA.CG.4436
#> ENSG00000186092      0.00000     0.000000    0.0000000
#> ENSG00000187634     14.32882     1.129263    0.1090498
#> ENSG00000188976     19.17756    28.584328   19.7157980

## Interested epigenetic regulators
head(interested.ER)
#> [1] "ENSG00000133627" "ENSG00000101442" "ENSG00000156802" "ENSG00000140320" "ENSG00000174744" "ENSG00000094804"

## Predetermined immune gene sets
signature.list[1]
#> $Interferons
#> [1] "ENSG00000186803" "ENSG00000233816" "ENSG00000228083" "ENSG00000147885" "ENSG00000234829"
#> [6] "ENSG00000188379" "ENSG00000137080" "ENSG00000236637" "ENSG00000147873" "ENSG00000120235"
#> [11] "ENSG00000214042" "ENSG00000120242" "ENSG00000171855" "ENSG00000184995" "ENSG00000111537"
#> [16] "ENSG00000147896" "ENSG00000177047"

## Running the EPRIM 
test_res <- EPRIM(exp.profile = exp.profile, 
                   exp.profile.file = exp.profile.file,
                   interested.ER = interested.ER,  
                   signature.list = signature.list)

## Showing the GSEA results
test_res$ImmERres[1:3,] 
```



