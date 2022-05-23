args <- commandArgs(TRUE)
phen <- args[1]
chr <- args[2]

library(tidyverse)

print(paste0("------Running ", phen, "-------------"))
f1 <- paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/EA/EA_",
  phen, "_new_munged_chr", chr, ".tsv.gz")
f2 <- paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/AA/AA_",
  phen, "_new_munged_chr", chr, ".tsv.gz")

dd1 <- read_tsv(f1)
dd2 <- read_tsv(f2)

print("------------Start to combine-----------")
tmp <- dd1 %>%
  rename(Z1 = Z,
    BETA1 = BETA,
    SE1 = SE,
    N1 = N) %>%
  full_join(dd2 %>%
      rename(Z2 = Z,
        BETA2  = BETA,
        SE2 = SE,
        N2 = N),
    by = c("CHR", "SNP", "BP", "A1", "A2")) %>%
  mutate(Z = ifelse(!is.na(BETA1) & !is.na(BETA2) & !is.na(SE1) & !is.na(SE2),
    ((1 / SE1)^2 * BETA1 + (1 / SE2)^2 * BETA2) / sqrt((1 / SE1)^2 + (1 / SE2)^2),
    ifelse(is.na(Z1), Z2,
      ifelse(is.na(Z2), Z1, NA))),
    N = ifelse(!is.na(N1) & !is.na(N2), N1 + N2,
      ifelse(is.na(N1), N2,
        ifelse(is.na(N2), N1, NA)))) %>%
  select(CHR, SNP, BP, A1, A2, Z, N) %>%
  mutate(chisq = Z**2) %>%
  filter(chisq <= max(0.001 * max(N), 80)) %>%
  select(-chisq)

print(tmp[is.na(tmp$Z), ])
print(tmp[is.na(tmp$N), ])

dd <- tmp %>%
  filter(!is.na(Z)) %>%
  distinct()

if(length(unique(dd$SNP)) != nrow(dd)) {
  stop("duplicated SNPs")
}
  
print(paste0("After full join: ", nrow(tmp)))
print(paste0("After filtering NA: ", nrow(dd)))
#print(paste0("After distinct ", nrow(dd1)))

print("-----------finish------------------")
write_tsv(dd, file = paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/meta/meta_", phen,
 "_new_munged_chr", chr, ".tsv.gz"), quote_escape = FALSE)

