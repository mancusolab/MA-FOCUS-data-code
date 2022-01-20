library(tidyverse)
library(ggpubr)
library(broom)

# table disgenet meta category
cate <- c("Hematological measurement",
  "Abnormality, disorder, or disease of the immune system",
  "Abnormality or disease of blood and blood-forming tissues", 
  "Abnormality or disorder of the cardiovascular system")

cate <- c("Hematological measurement")


tt3 <- read_csv("../enrich/DisGeNET_meta_categories.csv") %>%
  mutate(`Blood-related` = ifelse(`meta_category` %in% cate, TRUE, FALSE))

tt3 <- read_csv("../enrich/DisGeNET_meta_categories.csv") %>%
  filter(meta_category %in% cate)

write_tsv(tt3, "../table/meta_cate.tsv", quote_escape = FALSE)

# get GWAS basic stats 
f <- list.files("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/gwas_stats", full.names = TRUE)
tt5 <- f %>% map_df(read_tsv, col_types = cols())

write_tsv(tt5, file = "/home1/zeyunlu/research/sub_mefocus/table/gwas_stats.tsv", quote_escape = FALSE)

# basic TWAS
load("../data/twas.RData")

tt6 <- twas_all %>%
  filter(!is.na(TWAS.Z) & POP != "meta") %>%
  group_by(PHEN, POP, .drop=FALSE) %>%
  summarize(n = n(),
    meanhsq = mean(HSQ, na.rm = TRUE),
    meanZ = mean(TWAS.Z, na.rm = TRUE),
    meanNZ = mean(TWAS.Z.NORM, na.rm = TRUE)) %>%
  left_join(twas_all %>%
      filter(!is.na(TWAS.Z) & POP != "meta") %>%
      group_by(PHEN, POP) %>%
      filter(TWAS.P < 0.05/n()) %>%
      summarize(sig = n(),
        hsq = mean(HSQ),
        meansigZ = mean(TWAS.Z, na.rm = TRUE),
        meansigNZ = mean(TWAS.Z.NORM, na.rm = TRUE)),
    by = c("POP", "PHEN")) %>%
  mutate(sig = ifelse(is.na(sig), 0, sig))

write_tsv(tt6, path = "../table/twas_basic.tsv", quote_escape = FALSE)


# correlation
dd <- read_tsv("../data/twas_gwas_corr2.tsv") %>%
  mutate(GWAS.t = GWAS/GWAS.SE,
    TWAS.t = TWAS/TWAS.SE,
    GWAS.P = 2*pnorm(abs(GWAS.t), lower.tail = F),
    TWAS.P = 2*pnorm(abs(TWAS.t), lower.tail = F)) %>%
  select(PHEN, GWAS, GWAS.SE, GWAS.P, TWAS, TWAS.SE, TWAS.P)

write_tsv(dd, "../table/twas_gwas_corr_final.tsv", quote_escape = FALSE)

# wd <- "/project/nmancuso_8/zeyunlu/projects/sub_mefocus/corr"
# 
# f1 <- list.files(path = wd, pattern = "raw", full.names = TRUE)
# 
# # add traits
# tmp <- read_tsv(f1[1], col_types = cols()) %>%
#   mutate(phen = gsub("_raw\\.tsv", "",
#     gsub("/project/nmancuso_8/zeyunlu/projects/sub_mefocus/corr/corr_", "", f1[1])))
# 
# for (i in 2:17){
#   tmp <- tmp %>%
#     bind_rows(read_tsv(f1[i], col_types = cols()) %>%
#         mutate(phen = gsub("_raw\\.tsv", "",
#           gsub("/project/nmancuso_8/zeyunlu/projects/sub_mefocus/corr/corr_", "", f1[i]))))
# }
# write_tsv(tmp, file = "/home1/zeyunlu/research/sub_mefocus/data/gwas_raw_corr.tsv", quote_escape = FALSE)
# 
# gw <- read_tsv("../data/gwas_raw_corr.tsv") %>%
#   mutate(SE = estimate/statistic)  %>%
#   filter(type == "norm") %>%
#   select(PHEN = phen, GWAS = estimate, GWAS.SE = SE, GWAS.P = p.value)
# 
# tw <- twas_all %>%
#   group_by(PHEN) %>%
#   filter(POP %in% "EA") %>%
#   mutate(EA = TWAS.Z,
#     EAN = TWAS.Z.NORM) %>%
#   select(PHEN, ID, EA, EAN) %>%
#   inner_join(twas_all %>%
#       group_by(PHEN) %>%
#       filter(POP %in% "AA") %>%
#       mutate(AA = TWAS.Z,
#         AAN = TWAS.Z.NORM) %>%
#       select(PHEN, ID, AA, AAN),
#     by = c("PHEN", "ID")) %>%
#   na.omit %>%
#   group_by(PHEN)
# 
# tw_corr <- tibble()
# for (i in unique(a$PHEN)) {
#   tmp <- tw %>%
#     filter(PHEN == i)
#   tmpC <- cor.test(tmp$EAN, tmp$AAN)
#   tw_corr <- tw_corr %>%
#     bind_rows(tidy(tmpC) %>%
#         mutate(PHEN = i))
# }
# 
# tt7 <- tw_corr %>%
#   mutate(TWAS.SE = estimate/statistic) %>%
#   select(PHEN, TWAS = estimate, TWAS.SE, TWAS.P = p.value) %>%
#   left_join(gw, by = "PHEN") %>%
#   mutate(GWAS.P = ifelse(GWAS.P == 0, "<0.0001", GWAS.P),
#     TWAS.P = ifelse(TWAS.P == 0, "<0.0001", TWAS.P),
#     t = (TWAS - GWAS)/ sqrt(TWAS.SE^2 + GWAS.SE^2),
#     p = pnorm(t, lower.tail = FALSE))
# 
# write_tsv(tt7, path = "../data/twas_gwas_corr.tsv", quote_escape = FALSE)

# basic FOCUS stats

load("../data/focus_aapower.RData")

tmp <- focus_AApower %>%
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


# fine-mapped regions
tmptt8 <- focus_AApower %>%
  distinct(PHEN, BLOCK) %>%
  group_by(PHEN) %>%
  summarize(nRfp = n()) %>%
  # regions that contain true genes
  left_join(tmp %>%
      filter(!grepl("NULL", ID)) %>%
      distinct(PHEN, POP, BLOCK) %>%
      group_by(PHEN, POP) %>%
      summarize(nRgene = n()),
    by = "PHEN") %>%
  select(PHEN, POP, nRfp, nRgene) %>%
  # CGS size without null
  left_join(tmp %>%
      filter(!grepl("NULL", ID)) %>%
      group_by(PHEN, POP) %>%
      summarize(cgs = n()),
    by = c("PHEN", "POP")) %>%
  # lead genes
  left_join(tmp %>%
      group_by(PHEN, POP, BLOCK) %>%
      arrange(-PIP) %>%
      filter(row_number() == 1 & !grepl("NULL", ID)) %>%
      group_by(PHEN, POP) %>%
      summarize(lead = n()),
    by = c("PHEN", "POP")) %>%
  # lead genes, mean PIP
  # lead genes, sd PIP
  left_join(tmp %>%
      group_by(PHEN, POP, BLOCK) %>%
      arrange(-PIP) %>%
      filter(row_number() == 1 & !grepl("NULL", ID)) %>%
      group_by(PHEN, POP) %>%
      summarize(leadmm = mean(PIP),
        leadss = sd(PIP)),
    by = c("PHEN", "POP")) %>%
  # estimated causal
  left_join(tmp %>%
      filter(!grepl("NULL", ID)) %>%
      group_by(POP, PHEN) %>%
      summarize(n = sum(PIP)),
    by = c("PHEN", "POP")) %>%
  # mean PIP without null
  # SD PIP without null
  left_join(tmp %>%
      filter(!grepl("NULL", ID)) %>%
      group_by(PHEN, POP, BLOCK) %>%
      summarize(meanPIP = mean(PIP, na.rm = TRUE),
        sdPIP = sd(PIP, na.rm = TRUE)) %>%
      group_by(PHEN, POP) %>%
      summarize(meanPIP = mean(meanPIP, na.rm = TRUE),
        sdPIP = mean(sdPIP, na.rm = TRUE)),
    by = c("PHEN", "POP")) %>%
  # null model mean PIP of all fine-mapped regions
  # null model SD PIP of all fine-mapped regions
  left_join(tmp %>%
      filter(grepl("NULL", ID)) %>%
      group_by(PHEN, POP) %>%
      summarize(meanPIPnull = mean(PIP, na.rm = TRUE),
        sdPIPnull = sd(PIP, na.rm = TRUE)),
    by = c("PHEN", "POP")) %>%
  mutate(POP = factor(POP, levels = c("ME", "meta", "EA", "AA"),
    labels = c("MA-FOCUS", "Baseline", "EA FOCUS", "AA FOCUS")))

write_tsv(tmptt8, "../table/focus_stats.tsv", quote_escape = FALSE)

tt81 <- tmptt8 %>%
  group_by(POP) %>%
  summarize_at(vars(-PHEN), list(name = mean))

write_tsv(tt81, path = "../table/focus_stats_mean.tsv", quote_escape = FALSE)

tt82 <- tmptt8 %>%
  group_by(POP) %>%
  summarize_at(vars(-PHEN, -leadmm, -leadss, -meanPIP, -sdPIP, -meanPIPnull,
    -sdPIPnull), list(name = sum))

write_tsv(tt82, path = "../table/focus_stats_sum.tsv", quote_escape = FALSE)


# silver

dd <- read_tsv("by_block_together2.tsv")


# Enrichment 
load("../data/enrich_aapower.RData")
mapper <- read_csv("../enrich/DisGeNET_meta_categories.csv")

cate <- c("Hematological measurement",
  "Abnormality, disorder, or disease of the immune system",
  "Abnormality or disease of blood and blood-forming tissues", 
  "Abnormality or disorder of the cardiovascular system")

cate <- c("Hematological measurement")

catelevel <- c("Hematological\nmeasurement",
  "Immune\nsystem",
  "Blood\ntissues", 
  "Cardiovascular\nsystem")

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


ddEn <- enrichment.results.df %>%
  left_join(mapper %>%
      select(Term, `meta_category`),
    by = "Term") %>%
  filter(!is.na(`meta_category`)) %>%
  filter(`meta_category` %in% cate) %>%
  # group_by(TEST, PHENO) %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  select(TEST, PHENO, GENESET, N.test.genes, Term, Adjusted.P.value, Odds.Ratio, Genes,  meta_category)


write_tsv(ddEn, path = "../table/enrich1.tsv", quote_escape = FALSE)

phens <- c("BAS", "EOS","HCT","HGB","LYM","MCH","MCHC","MCV","MON","MPV","NEU", "PLT","RBC","RDW","WBC")
ddEn <- enrichment.results.df %>%
  filter(Term %in% lk) %>%
  mutate(Term2 = as.character(factor(Term, levels = lk, labels = phens))) %>%
  filter(Term2== PHENO) %>%
  select(TEST, PHENO, GENESET, N.test.genes, Term, Adjusted.P.value, Odds.Ratio, Genes)
write_tsv(ddEn, path = "../table/enrich2.tsv", quote_escape = FALSE)


# Morris
tt <- read_csv("../morris/stinggene.csv", col_names = FALSE)
tt <- tt[1:41,]

a <- tmp %>%
  filter(ID %in% tt$X1)
unique(a$ID)

a1 <- tmp %>%
  filter(ID %in% tt$X1) %>%
  filter(POP == "ME")
nrow(a1)
mean(a1$PIP)
sd(a1$PIP)

a2 <- a %>%
  filter(ID %in% tt$X1) %>%
  filter(POP == "meta")

res <- a1 %>%
  ungroup() %>%
  select(BLOCK,  PHEN, ID, PIP) %>%
  mutate(TEST = "MA-FOCUS") %>%
  bind_rows(a2 %>%
      ungroup() %>%
      select(BLOCK,  PHEN, ID, PIP) %>%
      mutate(TEST = "Baseline"))

write_tsv(res, path = "../table/morris.tsv", quote_escape = FALSE)






# END


# PIP correlaiton
phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")
phens <- phens[order(phens)]

tmp2 <- tmp %>%
  filter(!grepl("NULL", ID)) %>%
  pivot_wider(names_from = POP, values_from = PIP)

tt11 <- tidy(cor.test(tmp2$ME, tmp2$meta)) %>%
  mutate(se = estimate/statistic,
    method = "Baseline",
    PHEN = "Total") %>%
  select(PHEN, method, estimate, se, p.value) %>%
  bind_rows(tidy(cor.test(tmp2$ME, tmp2$EA)) %>%
      mutate(se = estimate/statistic,
        method = "EA FOCUS",
        PHEN = "Total") %>%
      select(PHEN, method, estimate, se, p.value)) %>%
  bind_rows(tidy(cor.test(tmp2$ME, tmp2$AA)) %>%
      mutate(se = estimate/statistic,
        method = "AA FOCUS",
        PHEN = "Total") %>%
      select(PHEN, method, estimate, se, p.value))

for (phen in phens) {
  tmptt11 <- tmp2 %>%
    filter(PHEN %in% phen)
  
  tt11 <- tt11 %>%
    bind_rows(tidy(cor.test(tmptt11$ME, tmptt11$meta)) %>%
        mutate(se = estimate/statistic,
          method = "Baseline",
          PHEN = phen) %>%
        select(PHEN, method, estimate, se, p.value)) %>%
    bind_rows(tidy(cor.test(tmptt11$ME, tmptt11$EA)) %>%
        mutate(se = estimate/statistic,
          method = "EA FOCUS",
          PHEN = phen) %>%
        select(PHEN, method, estimate, se, p.value)) %>%
    bind_rows(tidy(cor.test(tmptt11$ME, tmptt11$AA)) %>%
        mutate(se = estimate/statistic,
          method = "AA FOCUS",
          PHEN = phen) %>%
        select(PHEN, method, estimate, se, p.value))
}

write_tsv(tt11, path = "../table/pip_corr.tsv", quote_escape = FALSE)

# enrichment-meta
load("../data/enrich_dis.RData")

mapper <- read_csv("../enrich/DisGeNET_meta_categories.csv")

cate <- c("Hematological measurement",
  "Abnormality, disorder, or disease of the immune system",
  "Abnormality or disease of blood and blood-forming tissues", 
  "Abnormality or disorder of the cardiovascular system")

catelevel <- c("Hematological\nmeasurement",
  "Immune\nsystem",
  "Blood\ntissues", 
  "Cardiovascular\nsystem")

enridd <- enrichment.results.df %>%
  left_join(mapper %>%
      select(Term, `meta_category`),
    by = "Term") %>%
  filter(!is.na(`meta_category`)) %>%
  filter(`meta_category` %in% cate) %>%
  mutate(`meta_category` = factor(`meta_category`, levels = cate, labels = catelevel))

tt12 <- enridd %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  group_by(`meta_category`) %>%
  mutate(total = n()) %>%
  arrange(desc(total)) %>%
  # mutate(`meta_category` = factor(`meta_category`)) %>%
  group_by(`meta_category`, TEST) %>%
  summarize(n = n(),
    total = mean(total)) %>%
  group_by(`meta_category`) %>%
  mutate(ord = desc(n)) %>%
  mutate(TEST = factor(TEST, levels = c("EA", "AA", "ME", "meta"),
    labels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip")))




# silver
dd <- read_tsv("silver/gene_ranking_together.tsv")
tmp <- dd %>%
  # filter(trait %in% filtered_trait$trait) %>%
  # filter(PIP.ME != PIP.meta) %>%
  # mutate(haha = PIP.meta - PIP.ME)
  # filter(!(PIP.ME %in% c(0, 1) | PIP.meta %in% c(0, 1))) %>%
  pivot_longer(c(PIP.ME, PIP.meta))

a <- summary(lm(value ~ name + data + phen, tmp))
head(a$coefficients)

# sting
tt <- read_csv("stinggene.csv", col_names = FALSE)
tt$X2[23:length(tt$X2)] <- 1
tt <- tt %>%
  mutate(IN = ifelse(X2 == 0, TRUE, FALSE)) %>%
  select(ID=X1, IN) %>%
  distinct() %>%
  group_by(ID) %>%
  mutate(n = n()) %>%
  filter(n == 1 | (n >= 2 & IN)) %>%
  select(-n)

res <- tibble()
phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")
for (i in phens) {
  for (j in c("ME", "meta", "EA", "AA")) {
    nn <- sym(paste0("IN.CRED.SET.", j))
    ff <- focus_AApower %>%
      filter(PHEN == i & ID != "NULL") %>%
      select(ID, FM = nn) %>%
      inner_join(tt, by = "ID") %>%
      group_by(FM, IN) %>%
      summarize(n = n())
    gg <- focus_AApower %>%
      filter(PHEN == i & ID != "NULL") %>%
      select(ID, FM = nn) %>%
      inner_join(tt, by = "ID") %>%
      mutate(FM = as.numeric(FM),
        IN = as.numeric(IN))
    tmp <- summary(glm(FM~IN, gg, family = binomial(link = "logit")))
    
    res <- bind_rows(res, tibble(EST = tmp$coefficients[2, 1],
      STD = tmp$coefficients[2, 2], PHEN = i, TEST=j))
    # mm <- c(ifelse(length(ff[ff$FM & ff$IN, ]$n) == 0, 0, ff[ff$FM & ff$IN, ]$n),
    #   ff[!ff$FM & ff$IN, ]$n,
    #   ff[ff$FM & !ff$IN, ]$n,
    #   ff[!ff$FM & !ff$IN, ]$n)
    # # tidy(chisq.test(matrix(mm, nrow = 2)))
    # res <- res %>%
    #   bind_rows(tidy(fisher.test(matrix(mm, nrow = 2))) %>%
    #       mutate(PHEN = i,
    #         POP = j))
  }
}

try <- res %>%
  # mutate(SE = (estimate - conf.low )/1.96) %>%
  filter(EST >= -10) %>%
  group_by(TEST) %>%
  mutate(W  = (1 / (STD^2)) / sum(1 / (STD^2))) %>%
  summarize(W.EST = sum(W*EST),
    W.SE = sqrt( 1 / sum( 1 / STD^2)))
pnorm(try$`W.EST`/ try$`W.SE`, lower.tail = FALSE)

try <- res %>%
  mutate(SE = (estimate - conf.low )/1.96) %>%
  filter(SE!=0) %>%
  group_by(POP) %>%
  mutate(W  = (1 / (SE^2)) / sum(1 / (SE^2))) %>%
  summarize(W.EST = sum(W*estimate),
    W.SE = sqrt( 1 / sum( 1 / SE^2)))

pnorm(try$`W.EST`/ try$`W.SE`, lower.tail = FALSE)


met <- res %>%
  group_by(POP) %>%
  summarize(stat = -2 * sum(log(p.value)))

pchisq(met$stat, 30, lower.tail = FALSE)

met <- res %>%
  filter(EST >-10) %>%
  group_by(TEST) %>%
  mutate(W  = (1 / (STD^2)) / sum(1 / (STD^2))) %>%
  summarize(W.EST = sum(W*EST),
    W.SE = sqrt( 1 / sum( 1 / STD^2)))

pnorm(met$`W.EST`/ met$`W.SE`, lower.tail = FALSE)



res <- tibble()
phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")
for (i in phens) {
  for (j in c("ME", "meta", "EA", "AA")) {
    pips <- sym(paste0("PIP.", j))
    ff <- focus_AApower %>%
      filter(PHEN == i & ID != "NULL") %>%
      select(ID, FM = nn) %>%
      inner_join(tt, by = "ID") %>%
      group_by(FM, IN) %>%
      summarize(n = n())
    gg <- focus_AApower %>%
      filter(PHEN == i & ID != "NULL") %>%
      select(ID, PIP = pips) %>%
      inner_join(tt, by = "ID") %>%
      mutate(
        IN = as.numeric(IN))
    tmp <- glm(IN~PIP, gg, family = binomial(link = "logit"))
    
    res <- bind_rows(res, tidy(tmp) %>%
        mutate(PHEN = i, TEST=j))
    # mm <- c(ifelse(length(ff[ff$FM & ff$IN, ]$n) == 0, 0, ff[ff$FM & ff$IN, ]$n),
    #   ff[!ff$FM & ff$IN, ]$n,
    #   ff[ff$FM & !ff$IN, ]$n,
    #   ff[!ff$FM & !ff$IN, ]$n)
    # # tidy(chisq.test(matrix(mm, nrow = 2)))
    # res <- res %>%
    #   bind_rows(tidy(fisher.test(matrix(mm, nrow = 2))) %>%
    #       mutate(PHEN = i,
    #         POP = j))
  }
}

try <- res %>%
  filter(term == "PIP") %>%
  #  filter(SE!=0) %>%
  group_by(TEST) %>%
  mutate(W  = (1 / (std.error^2)) / sum(1 / (std.error^2))) %>%
  summarize(W.EST = sum(W*estimate),
    W.SE = sqrt( 1 / sum( 1 / std.error^2)))

pnorm(try$`W.EST`/ try$`W.SE`, lower.tail = FALSE)


tt <- read_csv("stinggene.csv", col_names = FALSE)
tt <- tt[1:41,]


tmp1 <- focus_AApower %>%
  pivot_longer(c(PIP.EA, IN.CRED.SET.EA, PIP.AA, IN.CRED.SET.AA,
    PIP.ME, IN.CRED.SET.ME, PIP.meta, IN.CRED.SET.meta)) %>%
  mutate(name = gsub("IN\\.CRED\\.SET\\.", "IN\\.", name)) %>%
  separate(name, into = c("group", "POP"), sep = "\\.") %>%
  pivot_wider(names_from = group, values_from = value)  %>%
  filter(IN == 1) %>%
  ungroup() %>%
  group_by(PHEN, POP, BLOCK) %>%
  arrange(desc(PIP)) %>%
  mutate(rank = row_number()) %>%
  filter(POP %in% c("ME", "meta")) %>%
  ungroup

tmp2 <- focus_AApower %>%
  pivot_longer(c(PIP.EA, IN.CRED.SET.EA, PIP.AA, IN.CRED.SET.AA,
    PIP.ME, IN.CRED.SET.ME, PIP.meta, IN.CRED.SET.meta)) %>%
  mutate(name = gsub("IN\\.CRED\\.SET\\.", "IN\\.", name)) %>%
  separate(name, into = c("group", "POP"), sep = "\\.") %>%
  pivot_wider(names_from = group, values_from = value)  %>%
  # filter(IN == 1) %>%
  ungroup() %>%
  group_by(PHEN, POP, BLOCK) %>%
  arrange(desc(PIP)) %>%
  mutate(rank = row_number()) %>%
  filter(POP %in% c("ME", "meta")) %>%
  ungroup



tmp11 <- tmp1 %>%
  filter(ID %in% tt$X1) %>%
  select(PHEN, POP, PIP, ID)

tmp21 <- tmp2 %>%
  filter(ID %in% tt$X1) %>%
  select(PHEN, POP, PIP, ID)

mean(filter(tmp11, POP == "ME")$PIP)
mean(filter(tmp11, POP == "meta")$PIP)

mean(filter(tmp21, POP == "ME")$PIP)
mean(filter(tmp21, POP == "meta")$PIP)
# %>%
# pivot_wider(names_from = "POP", values_from = "PIP")
summary(lm(PIP ~ POP+PHEN, tmp11))
mean(tmp1$ME, na.rm = TRUE)
mean(tmp1$meta, na.rm = TRUE)

a <- tmp %>%
  filter(ID %in% tt$X1)
unique(a$ID)
a1 <- tmp %>%
  filter(ID %in% tt$X1) %>%
  filter(POP == "ME")
nrow(a1)

a2 <- a %>%
  filter(ID %in% tt$X1) %>%
  filter(POP == "meta")
nrow(a2)

# smith 
gg <- c("RDH13", "ACAD10")

tmp <- focus_AApower %>%
  pivot_longer(c(PIP.EA, IN.CRED.SET.EA, PIP.AA, IN.CRED.SET.AA,
    PIP.ME, IN.CRED.SET.ME, PIP.meta, IN.CRED.SET.meta)) %>%
  mutate(name = gsub("IN\\.CRED\\.SET\\.", "IN\\.", name)) %>%
  separate(name, into = c("group", "POP"), sep = "\\.") %>%
  pivot_wider(names_from = group, values_from = value)  %>%
  filter(IN == 1) %>%
  ungroup() %>%
  group_by(PHEN, POP, BLOCK) %>%
  arrange(desc(PIP)) %>%
  mutate(rank = row_number()) %>%
  filter(POP %in% c("ME", "meta")) %>%
  ungroup()

try <- tmp %>%
  filter(ID %in% gg) %>%
  select(PHEN, POP, PIP, ID) %>%
  pivot_wider(names_from = "POP", values_from = "PIP")

mean(try$ME)
mean(try$meta)

t.test(tmp$PIP.ME, tmp$PIP.meta, alternative = "greater")



