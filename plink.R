
assoc <- read.table("res_allelic.assoc", header = TRUE)


significant_snps <- nrow(assoc[assoc$P < 0.05, ])


bonferroni_threshold <- 0.05 / nrow(assoc)

bonferroni_snps <- nrow(assoc[assoc$P < bonferroni_threshold, ])


cat("Number of SNPs with p < 0.05:", significant_snps, "\n")
cat("Bonferroni-adjusted p-value threshold:", bonferroni_threshold, "\n")
cat("Number of SNPs with p < Bonferroni threshold:", bonferroni_snps, "\n")

geno <- read.table("res_geno.assoc.logistic", header = TRUE)

geno_significant_snps <- nrow(geno[geno$P < 0.05, ])

geno_bonf_t <- 0.05 / nrow(geno)

geno_bonf_snps <- nrow(geno[geno$P < geno_bonf_t, ])

cat("Number of SNPs with p < 0.05", snps_coch, "\n")
cat("Bonferroni-adjusted p-value threshold:", bonf_t_coch, "\n")
cat("Number of SNPs with p < Bonferroni threshold:", bonf_snps_coch)

for (line in 1:nrow(geno)) {
    if (!is.na(geno$P[line]) && geno$P[line] < geno_bonf_t) {
        print(geno[line, ])
    }
}

for (line in 1:nrow(geno)) {
    if (is.na(geno$P[line])) {
        print(geno[line, ])
    }
}

for (line in 1:nrow(assoc)) {
    if (is.na(assoc$P[line])) {
        print(assoc[line, ])
    }
}

geno$BH_p <- p.adjust(geno$P, method = 'BH')
fdr_thresh <- 0.05
sigma_snps <- geno[geno$BH_p < fdr_thresh, ]

cat("Number of significant SNPs (FDR < 0.05):", nrow(sigma_snps), "\n")
print(sigma_snps)
