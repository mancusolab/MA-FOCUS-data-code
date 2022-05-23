args <- commandArgs(TRUE)
phen <- args[1]
library(tidyverse)

f1 <- paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/EA_", phen, "_GWAMA/EA_", phen,
    "_new_munged.sumstats.gz")

dd1 <- read_table2(f1) %>%
    mutate(chisq = Z**2) %>%
    filter(chisq <= max(0.001 * max(N), 80)) %>%
    select(-chisq)

nr1 <- nrow(dd1)

ct1 <- 0
colnames(dd1)
for (i in 1:22) {
    res1 <- dd1 %>%
        filter(CHR %in% i)
    ct1 <- ct1 + nrow(res1)
    write_tsv(res1,
        file = paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/EA_original/EA_original",
            phen, "_new_munged_chr", i, ".tsv.gz"), quote_escape = FALSE)
}

if(ct1 != nr1) {
    stop("something wrong")
}


f1 <- paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/AA_", phen, "_GWAMA/AA_", phen,
    "_new_munged.sumstats.gz")

dd1 <- read_table2(f1) %>%
    mutate(chisq = Z**2) %>%
    filter(chisq <= max(0.001 * max(N), 80)) %>%
    select(-chisq)

nr1 <- nrow(dd1)

ct1 <- 0
colnames(dd1)
for (i in 1:22) {
    res1 <- dd1 %>%
        filter(CHR %in% i)
    ct1 <- ct1 + nrow(res1)
    write_tsv(res1,
        file = paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/AA_original/AA_original",
            phen, "_new_munged_chr", i, ".tsv.gz"), quote_escape = FALSE)
}

if(ct1 != nr1) {
    stop("something wrong")
}
