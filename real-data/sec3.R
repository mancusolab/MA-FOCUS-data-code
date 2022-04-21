library(tidyverse)
library(broom)

genoa_ea <- read_tsv("../vg/genoa_hr_ea.tsv", col_names = FALSE)
colnames(genoa_ea) <- c("GENE", "GENEID", "POP", "HSQ", "HSQSE", "HSQP", "VG", "VGSE")

genoa_aa <- read_tsv("../vg/genoa_hr_aa.tsv", col_names = FALSE)
colnames(genoa_aa) <- c("GENE", "GENEID", "POP", "HSQ", "HSQSE", "HSQP", "VG", "VGSE")


mean(genoa_ea$HSQ)
tidy(t.test(genoa_ea$HSQ))

mean(genoa_aa$HSQ)
tidy(t.test(genoa_aa$HSQ))

# focusing on the 4,646 genes whose expression was significantly
# heritable in at least one of the cohorts, cis-g2 estimates were positively correlated
# across populations
genoa_all <- read_tsv("../data/genoa_her_total.tsv")

tmp <- genoa_all %>%
  select(GENE, HSQ, POP) %>%
  pivot_wider(names_from = POP, values_from = HSQ) %>%
  filter(complete.cases(.)) 
rho <- cor(tmp$ea, tmp$aa)
rho
z_rho <- 0.5*log((1+rho)/(1-rho))
pnorm(z_rho, mean = 0, sd = 1/sqrt(nrow(tmp)-3), lower.tail = FALSE)
pnorm(z_rho, mean = 0.5*log(1.99/0.01), sd = 1/sqrt(nrow(tmp)-3), lower.tail = TRUE)

tot <- read_tsv("../data/total_r2.tsv")
mean(tot$ea, na.rm = TRUE)
mean(tot$aa, na.rm = TRUE)
tidy(t.test(tot$ea, na.rm = TRUE))
tidy(t.test(tot$aa, na.rm = TRUE))


genoa_ea <- filter(genoa_all, POP == "ea")
genoa_aa <- filter(genoa_all, POP == "aa")

tidy(cor.test(genoa_ea$HSQ, genoa_ea$MODELCV.R2))
tidy(cor.test(genoa_aa$HSQ, genoa_aa$MODELCV.R2))

geuvadis_all <- read_tsv("../data/geuvadis_her_total.tsv")

geuvadis_EUR <- filter(geuvadis_all, POP == "EUR") %>%
  as.data.frame()
geuvadis_YRI <- filter(geuvadis_all, POP == "YRI") %>%
  as.data.frame()

ea <- genoa_ea %>%
  select(GENE, OAVG = HSQ) %>%
  inner_join(geuvadis_EUR %>%
      select(GENE, DISVG = HSQ),
    by = "GENE")
tidy(t.test(ea$OAVG))
tidy(t.test(ea$DISVG))
tidy(t.test(ea$OAVG, ea$DISVG))


aa <- genoa_aa %>%
  select(GENE, OAVG = HSQ) %>%
  inner_join(geuvadis_YRI %>%
      select(GENE, DISVG = HSQ),
    by = "GENE")
tidy(t.test(aa$OAVG))
tidy(t.test(aa$DISVG))
tidy(t.test(aa$OAVG, aa$DISVG))



# We found an average out-of-sample r^2 estimate of 0.04 for EUR and 0.05 for YRI
mean(tot$geuv_eur, na.rm = TRUE)
tidy(t.test(tot$geuv_eur))
mean(tot$geuv_yri, na.rm = TRUE)
tidy(t.test(tot$geuv_yri))

# Estimates were significantly correlated with estimates of GEUVADIS 
tidy(cor.test(geuvadis_EUR$r2, geuvadis_EUR$HSQ))
tidy(cor.test(geuvadis_YRI$r2, geuvadis_YRI$HSQ))
library(cocor)
data("aptitude")
tmp <- list(EUR = geuvadis_EUR, YRI = geuvadis_YRI)
cocor(~r2 + HSQ | r2 + HSQ, tmp, alternative = "greater")


# Indeed, we found that cis-h_g^2 for AA and YRI were less correlated than EA and EUR 
tidy(cor.test(ea$OAVG, ea$DISVG))
tidy(cor.test(aa$OAVG, aa$DISVG))
tmp <- list(EUR = as.data.frame(ea), YRI = as.data.frame(aa))
cocor(~OAVG + DISVG | OAVG + DISVG, tmp, alternative = "greater")


tidy(t.test(tot$geuv_yri_s))
tidy(t.test(tot$geuv_eur_s))
tidy(t.test(tot$geuv_eur, tot$geuv_eur_s))
tidy(t.test(tot$geuv_yri, tot$geuv_yri_s))

load("./data/twas.RData")

twas_all %>%
  filter(POP != "meta") %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/4579) %>%
  group_by(PHEN) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

# TWAS region
# twas_re <- twas_all %>%
#   filter(POP!="meta") %>%
#   select(CHR, P0, P1)
# write_tsv(twas_re, path = "../data/sig_region/twas_all_region.bed", quote_escape = FALSE, col_names = FALSE)

read_table2("../data/sig_region/twas_all_ld_region.bed", col_names = FALSE) %>%
  distinct()

# TWAS sig region
tmp <- twas_all %>%
  group_by(PHEN, POP) %>%
  summarize(n = n())

sig <- twas_all %>%
  filter(POP != "meta") %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/4579)

# report TWAS sig hits
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
  tmp <- read_tsv(paste0("../data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_AA.bed"), col_names = FALSE) %>%
    distinct() %>%
    mutate(PHEN = phen,
      POO="AA")
  res <- bind_rows(res, tmp)
  tmp <- read_tsv(paste0("../data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_EA.bed"), col_names = FALSE) %>%
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

# report shared sig hits
a <- twas_all %>%
  filter(POP != "meta") %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/4579) %>%
  group_by(ID, PHEN) %>%
  summarize(n = n()) %>%
  filter(n>1)
nrow(a)
length(unique(a$ID))

# 28 (17 unique) TWAS significant genes in 23 (11 unique) regions were shared across two ancestries. 
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

# read in regions with GWAS signals 
res2 <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("../data/sig_region/gwas_ld_sig/gwas_sig_ld_region_", phen, ".bed"), col_names = FALSE, col_types = cols()) %>%
    distinct() %>%
    mutate(PHEN = phen)
  res2 <- bind_rows(res2, tmp)
}

# Of the 8,243 (1,052 unique) LD blocks that contained GWAS 
res2 %>%
  select(-PHEN) %>%
  distinct()
nrow(distinct(res2))
nrow(distinct(res2, X1, X2, X3))

# check how many GWAS contains the TWAS
# preparation
tmp <- res2 %>%
  distinct() %>%
  mutate(GWAS = 1) %>%
  full_join(res %>%
      distinct() %>%
      mutate(TWAS = 1),
    by = c("X1", "X2", "X3", "PHEN")) %>%
  filter(!(is.na(GWAS) & is.na(TWAS)))

# how many TWAS is not NA
# 2,940 (598 unique) also exhibited TWAS signals in either ancestry.
tmp %>%
  filter(!is.na(TWAS) & !is.na(GWAS))

tmp %>%
  filter(!is.na(TWAS) & !is.na(GWAS)) %>%
  select(X1, X2, X3) %>%
  distinct()

# how many TWAS is not GWAS
tmp %>%
  filter(is.na(GWAS) & !is.na(TWAS))

tmp %>%
  filter(is.na(GWAS) & !is.na(TWAS)) %>%
  select(X1, X2, X3) %>%
  distinct()


#  these LD blocks contained average GWAS chi-squared statistics of 2.07, 1.05, and 1.61 for EA, AA, and meta-analysis,
# which are larger than other LD blocks that did not show GWAS signals 
phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")
res <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("../data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_AA.bed"), col_names = FALSE) %>%
    distinct() %>%
    mutate(PHEN = phen,
      POO="AA")
  res <- bind_rows(res, tmp)
  tmp <- read_tsv(paste0("../data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_EA.bed"), col_names = FALSE) %>%
    distinct() %>%
    mutate(PHEN = phen,
      POO="EA")
  res <- bind_rows(res, tmp)
}


res2 <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("../data/sig_region/gwas_ld_sig/gwas_sig_ld_region_", phen, ".bed"), col_names = FALSE, col_types = cols()) %>%
    distinct() %>%
    mutate(PHEN = phen)
  res2 <- bind_rows(res2, tmp)
}

tmp <- res2 %>%
  distinct() %>%
  mutate(GWAS = 1) %>%
  full_join(res %>%
      select(-POO) %>%
      distinct() %>%
      mutate(TWAS = 1),
    by = c("X1", "X2", "X3", "PHEN")) %>%
  filter(!(is.na(GWAS) & is.na(TWAS)))

ou <- tmp %>%
  filter(is.na(GWAS) & !is.na(TWAS))
colnames(ou)[1:3] <- c("CHR", "P0","P1")
# write_tsv(ou, "../review/twas_not_gwas_sig_region.tsv", quote_escape = FALSE)

ou %>%
  uniute(P0, P1, sep = )
load("../data/focus38_2.RData")


gwa <- read_tsv("../review/gwas_all_z2_signals.tsv")
modi_gwa <- read_tsv("../review/all_gwas_not_sig_ld.bed", col_names = FALSE) %>%
  rename(CHR=X1,
    P0=X2,
    P1=X3,
    PHEN=X4) %>%
  left_join(gwa, by = c("CHR", "P0", "P1", "PHEN")) %>%
  pivot_wider(names_from = POP, values_from = meanchi:minchi)

ddsig <- read_tsv("../review/gwas_sig_notTWAS.tsv") %>%
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
mean(ddsig$meanchi_EA)
mean(ddsig$meanchi_AA)
mean(ddsig$meanchi_meta)
(sum(tmp_ea > mean(ddsig$meanchi_EA))+1)/(n+1)
(sum(tmp_aa[!is.na(tmp_aa)] > mean(ddsig$meanchi_AA))+1)/(n+1)
(sum(tmp_meta > mean(ddsig$meanchi_meta))+1)/(n+1)

res %>%
  group_by(X1, X2, X3, PHEN) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  distinct(X1, X2, X3, PHEN) %>%
  mutate(idx = 1) %>%
  right_join(ddsig %>%
      select(X1 = CHR, X2 = P0, X3 = P1, PHEN),
    by = c("X1", "X2", "X3", "PHEN")) %>%
  filter(idx == 1)





# how many traits
# 
# tmp %>%
#   filter(is.na(GWAS) & !is.na(TWAS)) %>%
#   distinct(PHEN)
# 
# 
# tmp %>%
#   filter(is.na(GWAS)& !is.na(TWAS)) %>%
#   distinct(X1, X2, X3)

# multiple TWAS region
res <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("../data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_AA.bed"), col_names = FALSE) %>%
    mutate(PHEN = phen,
      POO="AA")
  res <- bind_rows(res, tmp)
  tmp <- read_tsv(paste0("../data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen, "_EA.bed"), col_names = FALSE) %>%
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

# unique
res %>%
  group_by(X1, X2, X3, PHEN, POO) %>%
  summarize(n = n()) %>%
  filter(n >=2) %>%
  ungroup() %>%
  distinct(X1, X2, X3)


#  Interestingly, we found that across-population correlations are higher on average for
# TWAS compared to GWAS (r = 0.061 and 0.052, respectively, P=0.015

tmp <- read_tsv("../data/twas_gwas_corr2.tsv") %>%
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

