args <- commandArgs(TRUE)
phen <- args[1]
chr <- args[2]
library(tidyverse)

dd <- read_table2(paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/EA/EA_",
  phen, "_new_munged_chr", chr, ".tsv.gz")) %>%
  filter(abs(Z) > qnorm(5e-8 / 2, lower.tail = FALSE)) %>%
  mutate(P1 = as.numeric(BP),
    P0 = P1 - 1) %>%
  select(CHR, P0, P1)

write_tsv(dd, file = paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_sig_region/chr/EA_",
  phen ,"_sig_region_chr", chr, ".bed"), quote_escape = FALSE)

dd <- read_table2(paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/AA/AA_",
  phen, "_new_munged_chr", chr, ".tsv.gz")) %>%
  filter(abs(Z) > qnorm(5e-8 / 2, lower.tail = FALSE)) %>%
  mutate(P1 = as.numeric(BP),
    P0 = P1 - 1) %>%
  select(CHR, P0, P1)

write_tsv(dd, file = paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_sig_region/chr/AA_",
  phen ,"_sig_region_chr", chr, ".bed"), quote_escape = FALSE)

dd <- read_table2(paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/meta/meta_",
  phen, "_new_munged_chr", chr, ".tsv.gz")) %>%
  filter(abs(Z) > qnorm(5e-8 / 2, lower.tail = FALSE)) %>%
  mutate(P1 = as.numeric(BP),
    P0 = P1 - 1) %>%
  select(CHR, P0, P1)

write_tsv(dd, file = paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_sig_region/chr/meta_",
  phen ,"_sig_region_chr", chr, ".bed"), quote_escape = FALSE)




