library(tidyverse)
library(broom)

load("./data/focus_all.RData")

phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")

focus_analysis <- tibble()
for (phen in phens) {
  # twas_sig_ld_both_LYM.bed
  LD1 <- read_tsv(paste0("./data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen,
    "_EA.bed"), col_names = FALSE)
  LD2 <- read_tsv(paste0("./data/sig_region/twas_ld_sig/twas_sig_ld_region_", phen,
    "_AA.bed"), col_names = FALSE)
  if (nrow(LD1) != 0 & nrow(LD2) != 0) {
    tot <- inner_join(LD1, LD2, by = c("X1", "X2","X3")) %>%
      distinct() %>%
      unite("X0", X1, X2,  sep = ":") %>%
      unite("haha", X0, X3, sep = "..") %>%
      unlist()
    if (length(tot) != 0 ){
      focus_analysis <- focus_analysis2 %>%
        bind_rows(focus_all %>%
            filter(PHEN == phen & BLOCK %in% tot))
    }
  }
}

# save(focus_analysis, file = "./data/focus_analysis.RData")


load("./data/focus_analysis.RData")
load("./data/twas.RData")

focus_analysis <- filter(focus_analysis, !is.na(ID))

# we applied MA-FOCUS to TWAS results for 11 blood traits focusing on 163 genes
# overlapping the 23 (11 unique) regions that contained TWAS signals
# for both EA and AA ancestry for a given trait (see Methods)

length(unique(focus_analysis$PHEN))
length(unique(focus_analysis$ID))

focus_analysis %>%
  distinct(PHEN, BLOCK)
length(unique(focus_analysis$BLOCK))

tmp <- focus_analysis %>%
  pivot_longer(c(PIP.EA, IN.CRED.SET.EA, PIP.AA, IN.CRED.SET.AA,
    PIP.ME, IN.CRED.SET.ME, PIP.meta, IN.CRED.SET.meta)) %>%
  mutate(name = gsub("IN\\.CRED\\.SET\\.", "IN\\.", name)) %>%
  separate(name, into = c("group", "POP"), sep = "\\.") %>%
  pivot_wider(names_from = group, values_from = value)  %>%
  filter(IN == 1) %>%
  ungroup() %>%
  group_by(PHEN, POP, BLOCK) %>%
  arrange(desc(PIP)) %>%
  mutate(rank = row_number())

# ). Across these 23 trait-specific regions, each contained an average of 7.35 TWAS significant
# associations across ancestries and 3.17 genes in the 90%-credible gene set, none of which
# included the null model.

sigt <- twas_all %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/4579 & POP != "meta") %>%
  distinct(ID, PHEN) %>%
  mutate(TWAS = 1) %>%
  ungroup()

focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  select(ID, BLOCK, PHEN) %>%
  left_join(sigt, by = c("ID", "PHEN")) %>%
  group_by(BLOCK, PHEN) %>%
  summarize(n = sum(TWAS, na.rm = TRUE)) %>%
  ungroup() %>%
  summarize(n = mean(n))

# average CGS per locus
tmp %>%
  filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") %>%
  group_by(PHEN, BLOCK) %>%
  summarize(cgs = n()) %>%
  ungroup() %>%
  summarize(n = mean(cgs))

tmp %>%
  filter(grepl("NULL", ID))

# We estimated an average of 2.88 causal genes per region by summing over local PIPs
# in the credible sets, with 19 out of 23 credible sets containing three or fewer genes

tmp %>%
  filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") %>%
  group_by(PHEN, BLOCK) %>%
  summarize(sumPIP = sum(PIP)) %>%
  ungroup() %>%
  summarize(n = mean(sumPIP))

tmp %>%
  filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") %>%
  group_by(PHEN, BLOCK) %>%
  summarize(cgs = n()) %>%
  filter(cgs <= 3)

# The average maximum PIP across credible sets was 0.99 (SD=0.02) and
# retained similar PIPs for the second and the third rank 

tmp %>%
  # filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") %>%
  filter(rank == 1) %>%
  group_by(PHEN, BLOCK) %>%
  ungroup() %>%
  distinct(PIP) %>%
  summarize(maxx = mean(PIP),
    sd = sd(PIP))

tmp %>%
  # filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") %>%
  filter(rank == 2) %>%
  group_by(PHEN, BLOCK) %>%
  ungroup() %>%
  distinct(PIP) %>%
  summarize(maxx = mean(PIP),
    sd = sd(PIP))

tmp %>%
  # filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") %>%
  filter(rank == 3) %>%
  group_by(PHEN, BLOCK) %>%
  ungroup() %>%
  distinct(PIP) %>%
  summarize(maxx = mean(PIP),
    sd = sd(PIP))

# we observed MA-FOCUS output higher means and smaller standard deviations of PIPs
# (P<0.05 for all tests MA-FOCUS obtained a smaller credible gene set on average (3.17)
# compared to the baseline (3.35); however, this result was not significant due to low
# statistical power (P=0.22

tt <- tmp %>%
  filter(IN == 1) %>%
  group_by(PHEN, BLOCK, POP) %>%
  summarise(m = n(),
    vpip = var(PIP),
    mpip = mean(PIP)) %>%
  pivot_wider(names_from = POP, values_from = m:mpip)

tidy(t.test(tt$m_ME, tt$m_meta, alternative = "less"))
# tidy(t.test(tt$m_ME, tt$m_EA, alternative = "less"))
# tidy(t.test(tt$m_ME, tt$m_AA, alternative = "less"))


# var.test(lm(tt$vpip_ME~1), lm(tt$vpip_EA~1), alternative = "less")
var.test(lm(tt$vpip_ME~1), lm(tt$vpip_meta~1), alternative = "less")
# var.test(lm(tt$vpip_ME~1), lm(tt$vpip_AA~1), alternative = "less")

tidy(t.test(tt$mpip_ME, tt$mpip_meta, alternative = "greater"))
# tidy(t.test(tt$mpip_ME, tt$mpip_EA, alternative = "greater"))
# tidy(t.test(tt$mpip_ME, tt$mpip_AA, alternative = "greater"))

# In addition, EA FOCUS did not prioritize 30 out of 73 trait-gene pairs in MA-FOCUS
# credible gene sets and missed 7 out of 23 lead genes

focus_analysis %>% 
  group_by(BLOCK, PHEN) %>%
  arrange(desc(PIP.ME)) %>%
  mutate(rank = row_number()) %>%
  filter(IN.CRED.SET.ME & !IN.CRED.SET.EA)


focus_analysis %>% 
  group_by(BLOCK, PHEN) %>%
  arrange(desc(PIP.ME)) %>%
  mutate(rank = row_number()) %>%
  filter(rank == 1 & !IN.CRED.SET.EA)

# We observed little support for a difference in the percentage of genes co-prioritized
# by AA FOCUS/MA-FOCUS (40.2%) compared with EA FOCUS/MA-FOCUS (49.5%; two-sample
# proportion test P=0.24
b1 <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.EA & IN.CRED.SET.ME)

b2 <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.AA & IN.CRED.SET.ME)

b3 <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.EA | IN.CRED.SET.ME)

b4 <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.AA | IN.CRED.SET.ME)

prop.test(c(nrow(b1), nrow(b2)), n=c(nrow(b3), nrow(b4)))

# We observed an average log-scale BF of 1.44 (SD=3.76), suggesting that credible-set
# genes underlying these blood traits are much more likely to be shared across ancestries
# than ancestry-specific (Figure S20). For instance, NPRL3 in the trait MCV had a
# log-scale BF of 17.1, which we discuss below 

b <- focus_analysis %>%
  filter(IN.CRED.SET.ME == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)),
    logbf = log(bf))

mean(b$logbf)
sd(b$logbf)
t.test(b$logbf, alternative = "greater")
filter(b, PHEN == "MCV" & ID == "NPRL3")

# we re-performed fine-mapping varying the maximum number of causal genes allowed in a
# configuration (see Methods) and found that while inferred PIPs were relatively stable
# (P=7.17×10^(-56)), credible gene sets sizes were sensitive to the upper
# bound on causal genes (see Supplementary Note

load("./data/focus_maxgene5.RData")

load("./data/focus_maxgene1.RData")

compMax <- function(dd1, dd2, method) {
  PIPvar <- sym(paste0("PIP.", method))
  INvar <- sym(paste0("IN.CRED.SET.", method))
  
  tmp1 <- dd1 %>%
    filter(!grepl("NULL", ID)) %>%
    select(PHEN, BLOCK, ID, p1 = !!PIPvar,  i1 = !!INvar) %>%
    inner_join(dd2 %>%
        select(PHEN, BLOCK, ID, p2 = !!PIPvar,  i2 = !!INvar),
      by = c("PHEN", "BLOCK", "ID"))
  print("PIP correlation:")
  print(tidy(cor.test(tmp1$p1, tmp1$p2)))
  
  cgs <- tmp1 %>%
    pivot_longer(cols = c(i1, i2)) %>%
    group_by(PHEN, BLOCK, name) %>%
    summarize(ss = sum(value)) %>%
    pivot_wider(values_from = ss, names_from = name)
  
  print("CGS comparison:")
  print(tidy(t.test(cgs$i1, cgs$i2)))
  
  backgr <- dd1 %>%
    filter(!grepl("NULL", ID)) %>%
    filter(!!INvar) %>%
    group_by(PHEN, BLOCK) %>%
    filter(!!PIPvar == max(!!PIPvar)) %>%
    select(PHEN, BLOCK, ID, !!PIPvar)
  
  tmp2 <- dd1 %>%
    filter(!grepl("NULL", ID)) %>%
    filter(!!INvar) %>%
    group_by(PHEN, BLOCK) %>%
    filter(!!PIPvar == max(!!PIPvar)) %>%
    select(PHEN, BLOCK, ID, !!PIPvar) %>%
    inner_join(dd2 %>%
        filter(!grepl("NULL", ID)) %>%
        filter(!!INvar) %>%
        group_by(PHEN, BLOCK) %>%
        filter(!!PIPvar == max(!!PIPvar)) %>%
        select(PHEN, BLOCK, ID, !!PIPvar),
      by = c("PHEN", "BLOCK", "ID"))
  print(paste0("Lead gene changes from ", nrow(backgr), ":"))
  print(nrow(tmp2))
  
  tmp3 <- dd1 %>%
    filter(!grepl("NULL", ID)) %>%
    filter(!!INvar) %>%
    group_by(PHEN, BLOCK) %>%
    mutate(rank1 = dense_rank(desc(!!PIPvar))) %>%
    filter(!!PIPvar == max(!!PIPvar)) %>%
    select(PHEN, BLOCK, ID, rank1, !!PIPvar) %>%
    right_join(dd2 %>%
        filter(!grepl("NULL", ID)) %>%
        # filter(IN.CRED.SET.ME) %>%
        group_by(PHEN, BLOCK) %>%
        mutate(rank2 = dense_rank(desc(!!PIPvar))) %>%
        select(PHEN, BLOCK, ID, rank2),
      by = c("PHEN", "BLOCK", "ID")) %>%
    filter(!is.na(rank1) & rank2 != 1)
  
  print(paste0("Rank change from lead gene of ", nrow(backgr), ":"))
  print(tmp3$rank2[order(tmp3$rank2)])
  
  tmp4 <- dd1 %>%
    filter(!grepl("NULL", ID)) %>%
    mutate(rank1 = dense_rank(desc(!!PIPvar))) %>%
    select(PHEN, BLOCK, ID, rank1) %>%
    left_join(dd2 %>%
        filter(!grepl("NULL", ID)) %>%
        # filter(IN.CRED.SET.ME) %>%
        mutate(rank2 = dense_rank(desc(!!PIPvar))) %>%
        select(PHEN, BLOCK, ID, rank2),
      by = c("PHEN", "BLOCK", "ID"))
  
  print("Rank corr:")
  print(tidy(cor.test(tmp4$rank1, tmp4$rank2)))
}

compMax(focus_analysis, focus_maxgene5, "ME")
compMax(focus_analysis, focus_maxgene1, "ME")

tmp4 <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.ME) %>%
  group_by(PHEN, BLOCK) %>%
  filter(PIP.ME == max(PIP.ME)) %>%
  mutate(rank1 = dense_rank(desc(PIP.ME))) %>%
  select(PHEN, BLOCK, ID, rank1) %>%
  left_join(focus_maxgene1 %>%
      filter(!grepl("NULL", ID)) %>%
      group_by(PHEN, BLOCK) %>%
      # filter(IN.CRED.SET.ME) %>%
      mutate(rank2 = dense_rank(desc(PIP.ME))) %>%
      select(PHEN, BLOCK, ID, rank2),
    by = c("PHEN", "BLOCK", "ID"))

nrow(filter(tmp4, rank2 == 1))
nrow(filter(tmp4, rank2 == 2)) + nrow(filter(tmp4, rank2 == 3))
nrow(filter(tmp4, rank2 > 10))

compMax(focus_analysis, focus_maxgene1, "meta")
compMax(focus_analysis, focus_maxgene1, "EA")
compMax(focus_analysis, focus_maxgene1, "AA")
compMax(focus_analysis, focus_maxgene5, "meta")
compMax(focus_analysis, focus_maxgene5, "EA")
compMax(focus_analysis, focus_maxgene5, "AA")

# of the 49 genes in GENOA-based credible sets, 17 had GEUVADIS-weight-derived results
# with a PIP correlation estimate of 0.84 (P=2.05×10^(-5)). In addition, 13 out of 17
# genes were the lead genes from GENOA, and among these 13 genes, ten remained lead
# genes from GEUVADIS

load("./data/focus_geuvadis.RData")
focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.ME) %>%
  distinct(ID)

tmp2 <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.ME) %>%
  group_by(PHEN, BLOCK) %>%
  mutate(rank1 = PIP.ME %in% max(PIP.ME)) %>%
  select(PHEN, BLOCK, ID, PIP.ME.N = PIP.ME, rank1) %>%
  left_join(focus_geuvadis %>%
      filter(IN.CRED.SET.ME) %>%
      filter(!grepl("NULL", ID)) %>%
      group_by(PHEN, BLOCK) %>%
      mutate(rank2 = PIP.ME %in% max(PIP.ME)) %>%
      select(PHEN, BLOCK, ID, PIP.ME.GEU = PIP.ME, rank2),
    by = c("PHEN", "BLOCK", "ID")) %>%
  drop_na()

tidy(cor.test(tmp2$PIP.ME.GEU, tmp2$PIP.ME.N))
filter(tmp2, rank1)
filter(tmp2, rank1 & rank2)

# we investigated genes to which MA-FOCUS assigned a high PIP (> 0.75) 
# we found that all 22 baseline-specific genes had low PIPs (< 0.1) from ancestry-specific
# fine-mapping in at least one ancestry, while 12 of these genes had a low PIP in both ancestries.
# On the other hand, only one out of 31 total MA-FOCUS-specific genes had PIPs
# below 0.1 in both AA and EA. We found six out of 31 total MA-FOCUS-specific genes
# achieved a moderate PIP of at least 0.25 in both EA and AA ancestry-specific 
# fine-mapping (ARNT, BAK1, MRPL28, NPRL3, PHTF1, and TARS2

tmp1 <- focus_analysis %>%
  filter(PIP.ME > 0.75 & IN.CRED.SET.ME & !IN.CRED.SET.meta)

tmp2 <- focus_analysis %>%
  filter(PIP.meta > 0.75 & !IN.CRED.SET.ME & IN.CRED.SET.meta)

nrow(tmp2)
all(tmp2$PIP.AA  < 0.1 | tmp2$PIP.EA < 0.1)
sum(tmp2$PIP.AA  < 0.1 & tmp2$PIP.EA < 0.1)

nrow(tmp1)
sum(tmp1$PIP.AA  < 0.1 & tmp1$PIP.EA < 0.1)

tmp1 %>%
  filter(PIP.AA  > 0.25 & PIP.EA > 0.25) %>%
  arrange(ID) %>%
  select(ID)



# enrichment 
# method 1
load("./data/enrich.RData")
mapper <- read_csv("./data/DisGeNET_meta_categories.csv")

enrichment.results.df %>%
  left_join(mapper %>%
      select(Term, `meta_category`),
    by = "Term") %>%
  filter(!is.na(`meta_category`)) %>%
  filter(`meta_category` %in% "Hematological measurement") %>%
  # group_by(TEST, PHENO) %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  group_by(TEST) %>%
  summarize(n = n())

enrichment.results.df %>%
  left_join(mapper %>%
      select(Term, `meta_category`),
    by = "Term") %>%
  filter(!is.na(`meta_category`)) %>%
  filter(`meta_category` %in% "Hematological measurement") %>%
  # group_by(TEST, PHENO) %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  group_by(TEST) %>%
  summarize(X2 = -2*sum(log(Adjusted.P.value)), dof = 2*n(), P = pchisq(X2, dof, lower.tail=FALSE))

# method 2
lk <- c("Blood basophil count (lab test)",
  "Eosinophil count result",
  "Hematocrit procedure",
  "Hemoglobin measurement",
  "Lymphocyte Count measurement",
  "Finding of Mean Corpuscular Hemoglobin",
  "Mean corpuscular hemoglobin concentration determination",
  "Mean Corpuscular Volume (result)",
  "Monocyte count procedure",
  "Platelet mean volume determination (procedure)",
  "Neutrophil count (procedure)",
  "Platelet Count measurement",
  "Red Blood Cell Count measurement",
  "RDW - Red blood cell distribution width result",
  "White Blood Cell Count procedure")

phens <- c("BAS", "EOS","HCT","HGB","LYM","MCH","MCHC","MCV","MON","MPV","NEU", "PLT","RBC","RDW","WBC")
enrichment.results.df %>%
  filter(Term %in% lk) %>%
  mutate(Term = as.character(factor(Term, levels = lk, labels = phens))) %>%
  filter(Term == PHENO) %>%
  group_by(TEST) %>%
  summarize(X2 = -2*sum(log(Adjusted.P.value)), dof = 2*n(), P = pchisq(X2, dof, lower.tail=FALSE))

# silver
dd <- read_tsv("./data/silver.tsv")
mean(dd$PIP.ME)
mean(dd$PIP.meta)


# WBC
tmp %>%
  ungroup() %>%
  distinct(BLOCK, PHEN)

tmp %>%
  ungroup() %>%
  distinct(BLOCK, PHEN) %>%
  group_by(PHEN) %>%
  summarize(n = n())

cgs <- tmp %>%
  filter(PHEN %in% "WBC" & POP %in% c("meta", "ME")) %>%
  group_by(BLOCK, PHEN, POP) %>%
  summarize(cgs = n()) %>%
  pivot_wider(names_from = POP, values_from = cgs)

tidy(t.test(cgs$ME, cgs$meta, alternative = "less"))

pip <- tmp %>%
  filter(PHEN %in% "WBC" & POP %in% c("meta", "ME")) %>%
  pivot_wider(names_from = POP, values_from = PIP)

tidy(t.test(pip$ME, pip$meta, alternative = "greater"))

tmp %>%
  filter(BLOCK %in% "1:151566405..154721871" & PHEN %in% "WBC" & POP %in% "ME") %>%
  mutate(rank = dense_rank(desc(PIP))) %>%
  select(ID, PIP, rank)



