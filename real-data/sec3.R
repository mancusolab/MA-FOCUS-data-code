library(tidyverse)
library(broom)

genoa_ea <- read_tsv("./data/genoa_heritability_ea_all_genes.tsv", col_names = FALSE)
colnames(genoa_ea) <- c("GENE", "GENEID", "POP", "HSQ", "HSQSE", "HSQP", "VG", "VGSE")

genoa_aa <- read_tsv("./data/genoa_heritability_aa_all_genes.tsv", col_names = FALSE)
colnames(genoa_aa) <- c("GENE", "GENEID", "POP", "HSQ", "HSQSE", "HSQP", "VG", "VGSE")


# , across all genes, cis-h_g^2cis- was significantly non-zero with an average of 0.058
# for EA compared to 0.072 for AA (P<1×10^(-100) for both tests).
mean(genoa_ea$HSQ)
tidy(t.test(genoa_ea$HSQ))

mean(genoa_aa$HSQ)
tidy(t.test(genoa_aa$HSQ))

# the 4,646 genes whose expression was significantly heritable in at least one of the cohorts,
# cis-h_g^2cis- estimates were positively correlated across ancestries 
# with r=0.45 (P<1×10^(-100) for both tests against 0 and 1; Figure 4A), 
genoa_all <- read_tsv("./data/genoa_heritability_analyzed.tsv")

tmp <- genoa_all %>%
  select(GENE, HSQ, POP) %>%
  pivot_wider(names_from = POP, values_from = HSQ) %>%
  filter(complete.cases(.)) 
rho <- cor(tmp$ea, tmp$aa)
rho
z_rho <- 0.5*log((1+rho)/(1-rho))
pnorm(z_rho, mean = 0, sd = 1/sqrt(nrow(tmp)-3), lower.tail = FALSE)
pnorm(z_rho, mean = 0.5*log(1.99/0.01), sd = 1/sqrt(nrow(tmp)-3), lower.tail = TRUE)

# We found that CV r^2 was significantly non-zero
# (AA CV r^2=0.11; EA CV r^2=0.10; <1×10^(-100) for both),
# which were strongly correlated with cis-h_g^2 estimates
# (r=0.93 with P<1×10^(-100) for both; Figure S12), 

tot <- read_tsv("./data/total_r2.tsv")
tidy(t.test(tot$ea, na.rm = TRUE))
tidy(t.test(tot$aa, na.rm = TRUE))

genoa_ea <- filter(genoa_all, POP == "ea")
genoa_aa <- filter(genoa_all, POP == "aa")

tidy(cor.test(genoa_ea$HSQ, genoa_ea$MODELCV.R2))
tidy(cor.test(genoa_aa$HSQ, genoa_aa$MODELCV.R2))


# We found  r^2 estimates between measured expression from GEUVADIS individuals
# and expression predicted using GENOA-based weights were significantly correlated
# with estimates of GEUVADIS cis-h_g^2, with r=0.85,0.56 for EUR and YRI (P<1×10^(-100)
# for both tests against 0; P<1×10^(-40)  for testing correlation difference; Figure 4BC).

geuvadis_all <- read_tsv("./data/geuvadis_heritability.tsv")

geuvadis_EUR <- filter(geuvadis_all, POP == "EUR") %>%
  as.data.frame()
geuvadis_YRI <- filter(geuvadis_all, POP == "YRI") %>%
  as.data.frame()

tidy(cor.test(geuvadis_EUR$r2, geuvadis_EUR$HSQ))
tidy(cor.test(geuvadis_YRI$r2, geuvadis_YRI$HSQ))
library(cocor)
data("aptitude")
tmp <- list(EUR = geuvadis_EUR, YRI = geuvadis_YRI)
cocor(~r2 + HSQ | r2 + HSQ, tmp, alternative = "greater")


# Indeed, we found that cis-h_g^2 for AA and YRI were less correlated than EA and EUR
# (r=0.27,0.49 with P<1×10^(-77)  for both tests against 0; P<1×10^(-40)
# for testing correlation difference)
ea <- genoa_ea %>%
  select(GENE, OAHSQ = HSQ) %>%
  inner_join(geuvadis_EUR %>%
      select(GENE, DISHSQ = HSQ),
    by = "GENE")

aa <- genoa_aa %>%
  select(GENE, OAHSQ = HSQ) %>%
  inner_join(geuvadis_YRI %>%
      select(GENE, DISHSQ = HSQ),
    by = "GENE")

tidy(cor.test(ea$OAHSQ, ea$DISHSQ))
tidy(cor.test(aa$OAHSQ, aa$DISHSQ))
tmp <- list(EUR = as.data.frame(ea), YRI = as.data.frame(aa))
cocor(~OAHSQ + DISHSQ | OAHSQ + DISHSQ, tmp, alternative = "greater")

# predicting LCL gene expression levels for GEUVADIS EUR individuals using GENOA AA weights
# (similarly for GEUVADIS YRI and GENOA EA) and estimated an average of r^2=0.040 for EUR
# and 0.033 for YRI (P<1×10^(-100) for both tests; Figure S15).
tidy(t.test(tot$geuv_yri_s))
tidy(t.test(tot$geuv_eur_s))

# Consistent with the previous work59, we found a decrease in accuracy for GEUVADIS YRI
# individuals compared to within-ancestry results (P=1.75×10^(-31)) and similar levels
# of accuracy for GEUVADIS EUR (P=0.09).
tidy(t.test(tot$geuv_eur, tot$geuv_eur_s))
tidy(t.test(tot$geuv_yri, tot$geuv_yri_s))

# In addition, we observed that GENOA data had a higher estimate of LCL cis-h_g^2 and its
# corresponding weight produced higher prediction accuracy for both ancestries than GEUVADIS
# data (P<6.33×10^(-15) for all tests; see Supplementary Note).

# We found lower estimates of cis-h_g^2 in GEUVADIS than in GENOA with estimates of 0.06
# compared to 0.15 for EUR/EA and 0.08 compared to 0.18 for YRI/AA (P<1×10^(-100) for all tests).
# We observed that using GENOA-based weights produced higher prediction accuracy with
# an average r^2 estimate of 0.04 compared to 0.03 for EA/EUR and of 0.05 compared
# to 0.02 for AA/YRI (P=6.33×10^(-15) and P=1.36×10^(-95); Figure S14).

tidy(t.test(ea$OAHSQ, ea$DISHSQ))
tidy(t.test(aa$OAHSQ, aa$DISHSQ))

tidy(t.test(tot$genoa_ea, tot$geuv_eur, alternative = "less"))
tidy(t.test(tot$genoa_aa, tot$geuv_yri, alternative = "less"))

load("./data/twas.RData")

# we conducted multi-ancestry TWAS for each of the 15 blood traits on 4,579 heritable genes
# in 989 unique independent regions (see Methods)
twas_all %>%
 distinct(ID)

read_table2("./data/sig_region/twas_all_ld_region.bed", col_names = FALSE) %>%
  distinct()

# we identified a total of 6,236 (2,009 unique) and 116 (57 unique) genome-wide TWAS
# significant genes in EA and AA, respectively, in 3,032 (622 unique) regions

sig <- twas_all %>%
  filter(POP != "meta") %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/4579)

nrow(filter(sig, POP == "EA"))
length(unique((filter(sig, POP == "EA")$ID)))
nrow(filter(sig, POP == "AA"))
length(unique((filter(sig, POP == "AA")$ID)))

phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")
# # to get TWAS sig gene position
# for (phen in phens) {
#   for (pop in c("EA", "AA")){
#     sig <- twas_all %>%
#       filter(PHEN %in% phen) %>%
#       filter(POP == pop) %>%
#       filter(TWAS.P < 0.05/4579) %>%
#       select(CHR, P0, P1)
#     
#     write_tsv(sig, path = paste0("../data/sig_region/twas_bp_sig/twas_sig_bp_region_",
#       phen, "_", pop, ".bed"), quote_escape = FALSE, col_names = FALSE)
#   }
# }


# after bedtoools, read in the regions
res <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("./data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_AA.bed"), col_names = FALSE) %>%
    distinct() %>%
    mutate(PHEN = phen,
      POO="AA")
  res <- bind_rows(res, tmp)
  tmp <- read_tsv(paste0("./data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_EA.bed"), col_names = FALSE) %>%
    distinct() %>%
    mutate(PHEN = phen,
      POO="EA")
  res <- bind_rows(res, tmp)
}

# twas sig in how many region
# genome-wide TWAS significant genes in EA and AA, respectively, in 3,032 (622 unique) regions
res %>%
  select(-POO) %>%
  distinct()

res %>%
  select(-POO, -PHEN) %>%
  distinct()

# We observed 28 (17 unique) genes significantly associated in both ancestry groups
# across 23 (11 unique) regions.
a <- twas_all %>%
  filter(POP != "meta") %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/4579) %>%
  group_by(ID, PHEN) %>%
  summarize(n = n()) %>%
  filter(n>1)
nrow(a)
length(unique(a$ID))

res %>%
  group_by(X1, X2, X3, PHEN) %>%
  summarize(n = n()) %>%
  filter(n ==2)

res %>%
  group_by(X1, X2, X3, PHEN) %>%
  summarize(n = n()) %>%
  filter(n ==2) %>%
  ungroup() %>%
  select(-PHEN) %>%
  distinct()

# Of the 8,243 trait-matched LD blocks that contained genome-wide significant signals
# (P<5×10^(-8)) in either ancestry or the meta-analysis, 2,940 also exhibited
# transcriptome-wide significant signals in either ancestry.
# Conversely, 115 trait-matched LD blocks contained transcriptome-wide significant
# signals that did not exhibit genome-wide significant signals

# read in regions with GWAS signals 
res2 <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("./data/sig_region/gwas_ld_sig/gwas_sig_ld_region_", phen, ".bed"), col_names = FALSE, col_types = cols()) %>%
    distinct() %>%
    mutate(PHEN = phen)
  res2 <- bind_rows(res2, tmp)
}

res2 %>%
  select(-PHEN) %>%
  distinct()
nrow(distinct(res2))

# preparation
tmp <- res2 %>%
  distinct() %>%
  mutate(GWAS = 1) %>%
  full_join(res %>%
      distinct() %>%
      mutate(TWAS = 1),
    by = c("X1", "X2", "X3", "PHEN")) %>%
  filter(!(is.na(GWAS) & is.na(TWAS)))

# 2940
tmp %>%
  filter(!is.na(TWAS) & !is.na(GWAS))

# how many TWAS is not GWAS
# 115
tmp %>%
  filter(is.na(GWAS) & !is.na(TWAS))

# We observed that these 115 regions exhibited greater GWAS signals on average when
# compared to their trait-matched genomic background (P=0.001,0.02,0.001 for EA, AA,
# and meta-analysis; see Methods), 

gwa <- read_tsv("./data/gwas_all_z2_signals.tsv")
modi_gwa <- read_tsv("./data/sig_region/all_gwas_not_sig_ld.bed", col_names = FALSE) %>%
  rename(CHR=X1,
    P0=X2,
    P1=X3,
    PHEN=X4) %>%
  left_join(gwa, by = c("CHR", "P0", "P1", "PHEN")) %>%
  pivot_wider(names_from = POP, values_from = meanchi:minchi)

ddsig <- read_tsv("./data/115_gwas_signal.tsv") %>%
  pivot_wider(names_from = POP, values_from = meanchi:minchi)

set.seed(123)
n <- 1000
tmp_ea <- NULL
tmp_aa <- NULL
tmp_meta <- NULL

for (i in 1:n) {
  idx <- sample(1:nrow(modi_gwa), 115)
  tmp_ea <- c(tmp_ea, mean(modi_gwa$meanchi_EA[idx]))
  tmp_aa <- c(tmp_aa, mean(modi_gwa$meanchi_AA[idx]))
  tmp_meta <- c(tmp_meta, mean(modi_gwa$meanchi_meta[idx]))
}

(sum(tmp_ea > mean(ddsig$meanchi_EA))+1)/(n+1)
(sum(tmp_aa[!is.na(tmp_aa)] > mean(ddsig$meanchi_AA))+1)/(n+1)
(sum(tmp_meta > mean(ddsig$meanchi_meta))+1)/(n+1)

# Of the 3,032 (622 unique) LD blocks containing TWAS hits, 1,329 (315 unique) contained
# multiple TWAS significant associations (average 3.56 genes per region), thus motivating
# the use of gene fine-mapping. 

# multiple TWAS region
res <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("./data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_AA.bed"), col_names = FALSE) %>%
    mutate(PHEN = phen,
      POO="AA")
  res <- bind_rows(res, tmp)
  tmp <- read_tsv(paste0("./data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_EA.bed"), col_names = FALSE) %>%
    mutate(PHEN = phen,
      POO="EA")
  res <- bind_rows(res, tmp)
}
res %>%
  group_by(X1, X2, X3, PHEN, POO) %>%
  summarize(n = n()) %>%
  filter(n >=2) %>%
  ungroup() %>%
  select(-POO, -n) %>%
  distinct()

res %>%
  group_by(X1, X2, X3, PHEN, POO) %>%
  summarize(n = n()) %>%
  filter(n >=2) %>%
  ungroup() %>%
  distinct(X1, X2, X3)

res %>%
  group_by(X1, X2, X3, PHEN, POO) %>%
  summarize(n = n()) %>%
  filter(n >=2) %>%
  ungroup() %>%
  summarize(nn = mean(n))


# Of the 6,352 significantly associated genes from GENOA across 15 traits and two ancestries,
# 4,315 were assayed in GEUVADIS, and 2,265 exhibited transcriptome-wide significance
# (P<0.05/4,579Z’ at nominal significance). 

load("./data/twas_geuvadis.RData")

twas_all %>%
  filter(POP %in% c("EA", "AA")) %>%
  filter(TWAS.P < 0.05/4579) %>%
  mutate(GENOSQ = TWAS.Z^2,
    GENONSQ = TWAS.Z.NORM^2) %>%
  select(ID, PHEN, POP, GENOSQ, GENONSQ)

twas_all %>%
  filter(POP %in% c("EA", "AA")) %>%
  filter(TWAS.P < 0.05/4579) %>%
  mutate(GENOSQ = TWAS.Z^2,
    GENONSQ = TWAS.Z.NORM^2) %>%
  select(ID, PHEN, POP, GENOSQ, GENONSQ) %>%
  inner_join(twas_geuvadis %>%
      filter(POP %in% c("EA", "AA")) %>%
      # filter(TWAS.P < 0.05/4579) %>%
      mutate(GEUVSQ = TWAS.Z^2,
        GEUVNSQ = TWAS.Z.NORM^2) %>%
      select(ID, PHEN, POP, GEUVSQ, GEUVNSQ),
    by = c("ID", "PHEN", "POP")) 


twas_all %>%
  filter(POP %in% c("EA", "AA")) %>%
  filter(TWAS.P < 0.05/4579) %>%
  mutate(GENOSQ = TWAS.Z^2,
    GENONSQ = TWAS.Z.NORM^2) %>%
  select(ID, PHEN, POP, GENOSQ, GENONSQ) %>%
  inner_join(twas_geuvadis %>%
      filter(POP %in% c("EA", "AA")) %>%
      filter(TWAS.P < 0.05/4579) %>%
      mutate(GEUVSQ = TWAS.Z^2,
        GEUVNSQ = TWAS.Z.NORM^2) %>%
      select(ID, PHEN, POP, GEUVSQ, GEUVNSQ),
    by = c("ID", "PHEN", "POP")) 

twas_total <- twas_all %>%
  filter(POP %in% c("EA", "AA")) %>%
  # filter(TWAS.P < 0.05/4579) %>%
  mutate(GENOSQ = TWAS.Z^2,
    GENONSQ = TWAS.Z.NORM^2) %>%
  select(ID, PHEN, POP, GENOSQ, GENONSQ) %>%
  inner_join(twas_geuvadis %>%
      filter(POP %in% c("EA", "AA")) %>%
      # filter(TWAS.P < 0.05/4579) %>%
      mutate(GEUVSQ = TWAS.Z^2,
        GEUVNSQ = TWAS.Z.NORM^2) %>%
      select(ID, PHEN, POP, GEUVSQ, GEUVNSQ),
    by = c("ID", "PHEN", "POP")) 

aa_twas <- twas_total %>% filter(POP %in% "AA")
ea_twas <- twas_total %>% filter(POP %in% "EA")

tidy(t.test(ea_twas$GENOSQ, ea_twas$GEUVSQ, alternative = "greater"))

tidy(t.test(aa_twas$GENOSQ, aa_twas$GEUVSQ, alternative = "greater"))


#  Interestingly, we found that across-population correlations are higher on average for
# TWAS compared to GWAS (r = 0.061 and 0.052, respectively, P=0.028

tmp <- read_tsv("./data/twas_gwas_corr.tsv") %>%
  mutate(TWAS.W  = (1/(TWAS.SE^2))/sum(1/(TWAS.SE^2)),
    GWAS.W  = (1/(GWAS.SE^2))/sum(1/(GWAS.SE^2)))
mean(tmp$TWAS)
mean(tmp$GWAS)

tmp2 <- tmp %>%
  summarize(TWAS.T = sum(TWAS.W*TWAS),
    TWAS.T.SE = sqrt( 1 / sum( 1 / TWAS.SE^2)),
    GWAS.T = sum(GWAS.W*GWAS),
    GWAS.T.SE = sqrt( 1 / sum( 1 / GWAS.SE^2)))

pnorm((tmp2$TWAS.T - tmp2$GWAS.T)/ sqrt(tmp2$TWAS.T.SE^2 + tmp2$GWAS.T.SE^2),lower.tail = FALSE)

#  In addition, we observed little support that TWAS-based gene effects sizes differ across
# EA and AA (P=0.57)

heter <- twas_all %>%
  mutate(beta = TWAS.Z/(sqrt(EQTL.GWAS.N * HSQ))) %>%
  filter(POP != "meta") %>%
  select(PHEN, ID, POP, beta) %>%
  pivot_wider(values_from = beta, names_from = POP)

t.test(heter$AA, heter$EA)

