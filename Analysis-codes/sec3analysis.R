library(tidyverse)
library(broom)

genoa_ea <- read_tsv("../vg/genoa_hr_ea.tsv", col_names = FALSE)
colnames(genoa_ea) <- c("GENE", "GENEID", "POP", "HSQ", "HSQSE", "HSQP", "VG", "VGSE")
mean(genoa_ea$VG)
tidy(t.test(genoa_ea$VG))

genoa_aa <- read_tsv("../vg/genoa_hr_aa.tsv", col_names = FALSE)
colnames(genoa_aa) <- c("GENE", "GENEID", "POP", "HSQ", "HSQSE", "HSQP", "VG", "VGSE")
mean(genoa_aa$VG)
tidy(t.test(genoa_aa$VG))

# focusing on the 4,646 genes whose expression was significantly
# heritable in at least one of the cohorts, cis-g2 estimates were positively correlated
# across populations
genoa_all <- read_tsv("../data/genoa_her_total.tsv")

tmp <- genoa_all %>%
  select(GENE, VG, POP) %>%
  pivot_wider(names_from = POP, values_from = VG) %>%
  filter(complete.cases(.)) 
rho <- cor(tmp$ea, tmp$aa)
z_rho <- 0.5*log((1+rho)/(1-rho))
pnorm(z_rho, mean = 0, sd = 1/sqrt(nrow(tmp)-3), lower.tail = FALSE)
pnorm(z_rho, mean = 0.5*log(1.99/0.01), sd = 1/sqrt(nrow(tmp)-3), lower.tail = TRUE)

genoa_ea <- filter(genoa_all, POP == "ea")
genoa_aa <- filter(genoa_all, POP == "aa")

tidy(cor.test(genoa_ea$VG, genoa_ea$MODELCV.R2))
tidy(cor.test(genoa_aa$VG, genoa_aa$MODELCV.R2))

geuvadis_all <- read_tsv("../data/geuvadis_her_total.tsv")
tmp <- geuvadis_all %>%
  filter(complete.cases(.)) %>%
  select(GENEID, VG, POP) %>%
  pivot_wider(names_from = POP, values_from = VG) %>%
  filter(complete.cases(.)) 
tidy(t.test(tmp$EUR))
tidy(t.test(tmp$YRI))

geuvadis_EUR <- filter(geuvadis_all, POP == "EUR") %>%
  as.data.frame()
geuvadis_YRI <- filter(geuvadis_all, POP == "YRI") %>%
  as.data.frame()
tidy(cor.test(geuvadis_EUR$r2, geuvadis_EUR$VG))
tidy(cor.test(geuvadis_YRI$r2, geuvadis_YRI$VG))

library(cocor)
data("aptitude")
tmp <- list(EUR = geuvadis_EUR, YRI = geuvadis_YRI)
cocor(~r2 + VG | r2 + VG, tmp, alternative = "greater")

mean(geuvadis_EUR$r2, na.rm = TRUE) / mean(geuvadis_EUR$VG, na.rm = TRUE)
mean(geuvadis_YRI$r2, na.rm = TRUE) / mean(geuvadis_YRI$VG, na.rm = TRUE)

ea <- genoa_ea %>%
  select(GENE, OAVG = VG) %>%
  inner_join(geuvadis_EUR %>%
      select(GENE, DISVG = VG),
    by = "GENE")
tidy(cor.test(ea$OAVG, ea$DISVG))

aa <- genoa_aa %>%
  select(GENE, OAVG = VG) %>%
  inner_join(geuvadis_YRI %>%
      select(GENE, DISVG = VG),
    by = "GENE")
tidy(cor.test(aa$OAVG, aa$DISVG))

tmp <- list(EUR = as.data.frame(ea), YRI = as.data.frame(aa))
cocor(~OAVG + DISVG | OAVG + DISVG, tmp, alternative = "greater")



load("../data/twas.RData")

# TWAS region
twas_re <- twas_all %>%
  select(CHR, P0, P1)
write_tsv(twas_re, path = "twas_all_region.bed", quote_escape = FALSE, col_names = FALSE)
read_table2("../data/indep_region_all.bed", col_names = FALSE) %>%
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

# report shared sig hits
a <- twas_all %>%
  filter(POP != "meta") %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/n()) %>%
  group_by(ID, PHEN) %>%
  summarize(n = n()) %>%
  filter(n>1)
nrow(a)
length(unique(a$ID))

phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")

# to get TWAS sig gene position
for (phen in phens) {
  sig <- twas_all %>%
    filter(PHEN %in% phen) %>%
    filter(POP != "meta") %>%
    group_by(PHEN, POP) %>%
    filter(TWAS.P < 0.05/n()) %>%
    ungroup() %>%
    select(CHR, P0, P1) %>%
    distinct()
  
  write_tsv(sig, path = paste0("../data/sig_region/twas_sig_region_", phen, ".bed"), quote_escape = FALSE, col_names = FALSE)
}

# after bedtoools, read in the regions
res <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("../data/sig_region/indep_region_", phen, ".bed"), col_names = FALSE) %>%
    mutate(PHEN = phen)
  res <- bind_rows(res, tmp)
}

nrow(distinct(res))
nrow(distinct(res, X1, X2, X3))


# read in regions with GWAS signals 
res2 <- tibble()
for (phen in phens) {
  tmp <- read_tsv(paste0("../data/sig_region/indep_region_gw_", phen, ".bed"), col_names = FALSE, col_types = cols()) %>%
    mutate(PHEN = phen)
  res2 <- bind_rows(res2, tmp)
}

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
tmp %>%
  filter(!is.na(TWAS) & !is.na(GWAS))

tmp %>%
  filter(!is.na(TWAS)) %>%
  distinct(X1, X2, X3)

# how many TWAS is not GWAS
tmp %>%
  filter(is.na(GWAS) & !is.na(TWAS))

# how many traits

tmp %>%
  filter(is.na(GWAS) & !is.na(TWAS)) %>%
  distinct(PHEN)


tmp %>%
  filter(is.na(GWAS)& !is.na(TWAS)) %>%
  distinct(X1, X2, X3)

# multiple TWAS region

res %>%
  group_by(X1, X2, X3, PHEN) %>%
  summarize(n = n()) %>%
  filter(n >=2)

# unique
res %>%
  group_by(X1, X2, X3, PHEN) %>%
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



# tot <- re %>%
#   mutate(one = 1) %>%
#   full_join(gw %>%
#       mutate(two=1),
#     by = c("X1", "X2", "X3"))
# 
# load("focus.RData")
# sigt <- twas_all %>%
#   group_by(PHEN, POP) %>%
#   filter(TWAS.P < 0.05/n() & POP != "meta") %>%
#   mutate(TWAS = 1) %>%
#   ungroup() %>%
#   select(ID, PHEN, TWAS, TWAS.P)
# 
# b <- focus_all %>%
#   filter(!grepl("NULL", ID)) %>%
#   select(ID, BLOCK, PHEN) %>%
#   left_join(sigt, by = c("ID", "PHEN")) %>%
#   group_by(BLOCK, PHEN) %>%
#   summarize(nTWAS = sum(TWAS, na.rm = TRUE),
#     n = n()) %>%
#   filter(nTWAS > 0)
# 
# c <- filter(b, nTWAS> 1)
# length(unique(b$BLOCK))
# length(unique(c$BLOCK))
# # testinng big vg
# 
# a <- read_tsv("ea_GRCh38_expr.bed.gz")
# tmp <- unlist(a[a$gene == "ACTG1", 5:377])
# tmp2 <- (tmp - mean(tmp))/sd(tmp)
# names(tmp2) <- NULL
# dd <- data.frame(FID = colnames(a)[5:377], IID = colnames(a)[5:377], GENE = tmp2)
# 
# write_tsv(dd, "/scratch/zeyunlu/tmp_expr.pheno", quote_escape = FALSE, col_names = FALSE)
# 
# 
