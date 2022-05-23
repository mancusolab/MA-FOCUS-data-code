library(tidyverse)
library(ggpubr)
library(broom)


# Supplementary tables 1, 3 and 5 are either manually typed or directly modified after downloading.

# Supplementary table 2
dd <- read_tsv("data/admixture_sample_size.tsv") %>%
  rename(`Population code` = name) %>%
  left_join(read_tsv("data/igsr_samples.tsv") %>%
      filter(!is.na(`Population code`)) %>%
      filter(!grepl(",", `Superpopulation name`)) %>%
      distinct(`Superpopulation code`, `Superpopulation name`,
        `Population code`, `Population name`),
    by = "Population code") %>%
  select(`Superpopulation code`, `Superpopulation name`,
    `Population code`, `Population name`, value)

# write_tsv(dd, "table-s2.tsv", quote_escape = FALSE)


# Supplementary table 4

tt3 <- read_csv("./data/DisGeNET_meta_categories.csv") %>%
  filter(meta_category %in% "Hematological measurement") %>%
  select(-MedGen_UID, -`HP/EFO`, -DOID)

# write_tsv(tt3, "./table/meta_cate.tsv", quote_escape = FALSE)

# Supplementary table 6 is generated on the high performance computing cluster (HPC)
# since it involves reading and analyzing large GWAS data 

args <- commandArgs(TRUE)
phen <- args[1]

f1 <- list.files("/project/nmancuso_8/zeyunlu/projects/sub_mefocus/filter_chen_gwas",
  pattern = paste0(phen, "_EA"), full.names = TRUE)
f2 <- list.files("/project/nmancuso_8/zeyunlu/projects/sub_mefocus/filter_chen_gwas",
  pattern = paste0(phen, "_AA"), full.names = TRUE)

dd1 <- f1 %>% map_df(read_tsv, col_types = cols(), col_names = FALSE)
dd2 <- f2 %>% map_df(read_tsv, col_types = cols(), col_names = FALSE)
colnames(dd1) <- c("CHR", "BP0", "BP1", "SNP", "A1", "A2", "SE", "BETA", "Z", "N")
colnames(dd2) <- c("CHR", "BP0", "BP1", "SNP", "A1", "A2", "SE", "BETA", "Z", "N")

res1 <- dd1 %>%
  summarize(meanZ = mean(Z),
    sdZ = sd(Z),
    medianZ = median(Z),
    normZ = mean(Z/sqrt(N)),
    meanN = mean(N),
    sdN = sd(N),
    medianN = median(N)) %>%
  mutate(pop = "EA",
    phen = phen) %>%
  select(phen, pop, meanZ, medianZ, sdZ, normZ, meanN, medianN, sdN)

res2 <- dd2 %>%
  summarize(meanZ = mean(Z),
    sdZ = sd(Z),
    medianZ = median(Z),
    normZ = mean(Z/sqrt(N)),
    meanN = mean(N),
    sdN = sd(N),
    medianN = median(N)) %>%
  mutate(pop = "AA",
    phen = phen) %>%
  select(phen, pop, meanZ, medianZ, sdZ, normZ, meanN, medianN, sdN)

res <- bind_rows(res1, res2)

# write_tsv(res,
#   file = paste0("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/gwas_stats/gwas_stats_",
#     phen, ".tsv"), quote_escape = FALSE)


# Supplementary table 7

load("./data/twas.RData")

tt <- twas_all %>%
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

# write_tsv(tt6, path = "./table/twas_basic.tsv", quote_escape = FALSE)

# Supplementary table 8

# GWAS correlation is calculated on the HPC
args <- commandArgs(TRUE)
phen <- args[1]

f1 <- list.files("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/EA",
  pattern = paste0(phen, "_"), full.names = TRUE)
f2 <- list.files("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/AA",
  pattern = paste0(phen, "_"), full.names = TRUE)

dd1 <- f1 %>% map_df(read_tsv, col_types = cols())
dd2 <- f2 %>% map_df(read_tsv, col_types = cols())

dd <- dd1 %>%
  mutate(CHR = as.character(CHR),
    SNP = as.character(SNP),
    BP = as.character(BP),
    EAZ = Z,
    EAnormZ = Z/sqrt(N)) %>%
  select(-Z, -N, -BETA, -SE) %>%
  inner_join(dd2 %>%
      mutate(CHR = as.character(CHR),
        SNP = as.character(SNP),
        BP = as.character(BP),
        AAZ = Z,
        AAnormZ = Z/sqrt(N)) %>%
      select(-Z, -N, -BETA, -SE),
    by = c("CHR", "SNP", "BP", "A1", "A2"))

# Jacknife
res <- tibble()
for (i in 1:floor(nrow(dd)/10000)) {
  idxStart <- 1+(i-1)*10000
  if (i == floor(nrow(dd)/10000)) {
    idxEnd <- nrow(dd)
  } else {
    idxEnd <- idxStart + 9999
  }
  tmp <- dd[-(idxStart:idxEnd),]
  res <- res %>%
    bind_rows(tidy(cor.test(tmp$EAnormZ, tmp$AAnormZ)))
}

ans <- tibble(PHEN = phen, corr = mean(res$estimate),
  SE = sqrt( ((nrow(res) - 1) / nrow(res)) * sum((res$estimate - mean(res$estimate))^2)))

write_tsv(ans,
  file = paste0("/project/nmancuso_8/zeyunlu/projects/sub_mefocus/corr_jack/corr_",
    phen, "_raw_jack.tsv"), quote_escape = FALSE)

# TWAS correlation

phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")

ans <- tibble()
for (phen in phens) {
  dd <- twas_all %>%
    filter(POP == "EA" & PHEN == phen) %>%
    select(ID, CHR, P0, P1, EA = TWAS.Z.NORM) %>%
    inner_join(twas_all %>%
        filter(POP == "AA" & PHEN == phen) %>%
        select(ID, AA = TWAS.Z.NORM),
      by = "ID") %>%
    arrange(CHR, P0, P1)
  # Jacknife
  res <- tibble()
  for (i in 1:floor(nrow(dd)/100)) {
    idxStart <- 1+(i-1)*100
    if (i == floor(nrow(dd)/100)) {
      idxEnd <- nrow(dd)
    } else {
      idxEnd <- idxStart + 99
    }
    tmp <- dd[-(idxStart:idxEnd),]
    res <- res %>%
      bind_rows(tidy(cor.test(tmp$EA, tmp$AA)))
  }
  ans <- ans %>%
    bind_rows(tibble(PHEN = phen, corr = mean(res$estimate),
      SE = sqrt( ((nrow(res) - 1) / nrow(res)) * sum((res$estimate - mean(res$estimate))^2))))
}

corrdd <- gwdd %>%
  select(PHEN, GWAS = corr, GWAS.SE = SE) %>%
  left_join(ans %>%
      select(PHEN, TWAS = corr, TWAS.SE = SE),
    by = "PHEN")

# write_tsv(corrdd, "/home1/zeyunlu/research/sub_mefocus/data/twas_gwas_corr.tsv")

# correlation
dd <- read_tsv("./data/twas_gwas_corr.tsv") %>%
  mutate(GWAS.t = GWAS/GWAS.SE,
    TWAS.t = TWAS/TWAS.SE,
    GWAS.P = 2*pnorm(abs(GWAS.t), lower.tail = F),
    TWAS.P = 2*pnorm(abs(TWAS.t), lower.tail = F)) %>%
  select(PHEN, GWAS, GWAS.SE, GWAS.P, TWAS, TWAS.SE, TWAS.P)

# write_tsv(dd, "../table/twas_gwas_corr_final.tsv", quote_escape = FALSE)

# Supplementary table S9

load("./data/focus_analysis.RData")

tmp1 <- focus_analysis %>%
  filter(PHEN %in% "WBC") %>%
  filter(IN.CRED.SET.meta) %>%
  group_by(BLOCK) %>%
  summarize(Baselinemean = mean(PIP.meta)) %>%
  left_join(focus_analysis %>%
      filter(PHEN %in% "WBC") %>%
      filter(IN.CRED.SET.ME) %>%
      group_by(BLOCK) %>%
      summarize(MEmean = mean(PIP.ME)),
    by = "BLOCK")

tmp2 <- focus_analysis %>%
  filter(PHEN %in% "WBC") %>%
  filter(IN.CRED.SET.ME) %>%
  group_by(BLOCK) %>%
  filter(PIP.ME == max(PIP.ME)) %>%
  select(BLOCK, ID, PIP.ME) %>%
  full_join(focus_analysis %>%
      filter(PHEN %in% "WBC") %>%
      filter(IN.CRED.SET.meta) %>%
      group_by(BLOCK) %>%
      filter(PIP.meta == max(PIP.meta)) %>%
      select(BLOCK, ID, PIP.meta),
    by = c("BLOCK", "ID"))

# Supplementary table s10
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


# fine-mapped regions
fres <- focus_analysis %>%
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

# write_tsv(fres, "./table/focus_stats.tsv", quote_escape = FALSE)


# Supplementary table s11 s12
# Enrichment analysis
options(stringsAsFactors=F)
library(enrichR)


focus_analysis2 <- as.data.frame(focus_analysis)
genesets <- c("DisGeNET")

# define functions ----
# ~ enrichment test function ----
do.enrichR.test <- function(DAT, TEST, PHENO, GENESETS, EXCL.NULL.CSs=T, INCL.TYPE="top.CS.gene", use.wgts=T) {
  dat <- DAT[DAT$PHEN==PHENO,]
  pip.col <- paste0("PIP.",TEST)
  in.CS.col <- paste0("IN.CRED.SET.",TEST)
  if (EXCL.NULL.CSs==F & INCL.TYPE=="all.CS.genes") {
    test.gene.set <- dat %>%
      filter(get(in.CS.col)==TRUE) %>%
      filter(ID!="NULL") %>%
      select(ID, all_of(pip.col))
  } else if (EXCL.NULL.CSs==F & INCL.TYPE=="top.CS.gene") {
    test.gene.set <- dat %>%
      group_by(BLOCK) %>%
      filter(ID!="NULL") %>%
      filter(get(pip.col)==max(get(pip.col))) %>%
      ungroup %>% select(ID, all_of(pip.col))
  } else if (EXCL.NULL.CSs==T & INCL.TYPE=="all.CS.genes") {
    in.blocks <- dat %>%
      group_by(BLOCK) %>%
      filter(get(pip.col)==max(get(pip.col))) %>%
      filter(ID!="NULL") %>%
      select(BLOCK)
    test.gene.set <- dat %>%
      filter(BLOCK %in% in.blocks$BLOCK) %>%
      filter(get(in.CS.col)==TRUE) %>%
      select(ID, all_of(pip.col))
  } else if (EXCL.NULL.CSs==T & INCL.TYPE=="top.CS.gene") {
    test.gene.set <- dat %>%
      group_by(BLOCK) %>%
      filter(get(pip.col)==max(get(pip.col))) %>%
      filter(ID!="NULL") %>%
      ungroup %>% select(ID, all_of(pip.col))
  } else if (INCL.TYPE=="all.genes") {
    EXCL.NULL.CSs <- F
    use.wgts <- T
    warning("'All genes' chosen to define test set - using a weighted list that includes all credible sets")
    test.gene.set <- dat %>%
      filter(ID!="NULL") %>%
      select(ID, all_of(pip.col))
  }
  test.gene.set <- as.data.frame(test.gene.set)
  N.test.genes <- nrow(test.gene.set)
  if (use.wgts==T) {
    result <- enrichr(test.gene.set, databases=GENESETS)
  } else {
    result <- enrichr(test.gene.set$ID, databases=GENESETS)
  }
  for (i in 1:length(result)) {
    if (nrow(result[[i]]) > 0) {
      result[[i]] <- cbind(N.test.genes, result[[i]])
    }
  }
  return(result)
}

# run enrichment analysis ----
PHENOS <- sort(unique(focus_analysis$PHEN))
TESTS <- c("EA", "AA", "ME", "meta")
enrichment.results <- list()
length(enrichment.results) <- length(TESTS)
names(enrichment.results) <- TESTS

for (TEST in TESTS) {
  for (PHENO in PHENOS) {
    cat(paste("Doing enrichment tests for", PHENO, TEST,"\n"))
    enrich.test <- do.enrichR.test(DAT=focus_analysis2,
      TEST=TEST,
      PHENO=PHENO,
      GENESETS=genesets,
      EXCL.NULL.CSs=F,
      INCL.TYPE="all.CS.genes",
      use.wgts=F)
    for (GENESET in genesets) {
      if (nrow(enrich.test[[GENESET]]) > 0) {
        if (is.null(enrichment.results[[TEST]])) {
          enrichment.results[[TEST]] <- cbind(PHENO, GENESET, enrich.test[[GENESET]])
        } else {
          enrichment.results[[TEST]] <- rbind(enrichment.results[[TEST]], cbind(PHENO, GENESET, enrich.test[[GENESET]]))
        }
      }
    }
  }
  enrichment.results[[TEST]] <- enrichment.results[[TEST]][order(enrichment.results[[TEST]]$P.value, decreasing=F),]
}

# reformat as a single data frame
for (TEST in TESTS) {
  if (exists("enrichment.results.df")) {
    enrichment.results.df <- rbind(enrichment.results.df, cbind(TEST, enrichment.results[[TEST]]))
  } else {
    enrichment.results.df <- cbind(TEST, enrichment.results[[TEST]])
  }
}
rm(enrich.test, enrichment.results, TEST)
for (i in 1:nrow(enrichment.results.df)) {
  enrichment.results.df$prop.Overlap[i] <- eval(parse(text=enrichment.results.df$Overlap[i]))
}

# save("enrichment.results.df", file="./data/enrich.RData")


load("./data/enrich.RData")
mapper <- read_csv("./data/DisGeNET_meta_categories.csv")

ddEn <- enrichment.results.df %>%
  left_join(mapper %>%
      select(Term, `meta_category`),
    by = "Term") %>%
  filter(!is.na(`meta_category`)) %>%
  filter(`meta_category` %in% "Hematological measurement") %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  select(TEST, PHENO, GENESET, N.test.genes, Term, Adjusted.P.value, Odds.Ratio, Genes,  meta_category)

# write_tsv(ddEn, path = "./table/enrich1.tsv", quote_escape = FALSE)

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
ddEn <- enrichment.results.df %>%
  filter(Term %in% lk) %>%
  mutate(Term2 = as.character(factor(Term, levels = lk, labels = phens))) %>%
  filter(Term2== PHENO) %>%
  select(TEST, PHENO, GENESET, N.test.genes, Term, Adjusted.P.value, Odds.Ratio, Genes)

# write_tsv(ddEn, path = "./table/enrich2.tsv", quote_escape = FALSE)

# Supplementary table s13 silver analysis
library(SilverStandardPerformance)
library(broom)

omim <- omim_based_silver_standard$table
orphanet <- orphanet_based_silver_standard$table

gene_mapper <- read_tsv("./data/gencode.v26.GRCh38.genes.only.tsv") %>%
  select(gene_name = symbol, gene = encode_id)

load("./data/focus_analysis.RData")
focus_gene <- focus_analysis %>%
  filter(!is.na(ID)) %>%
  filter(!grepl("NULL", ID)) %>%
  inner_join(gene_mapper, by = c("ID" = "gene_name"))

phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")

# Supplementary table s5
filtered_trait <- read_csv("./data/restricted_blood_traits.csv", col_types = cols(phecode = col_character())) %>%
  mutate(trait = ifelse(!is.na(phecode_trait), phecode_trait, hpo_trait)) %>%
  select(HPO, phecode, EFO, trait)

tmp <- focus_gene  %>%
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

res <- tibble()
pGeneMap1 <- omim_based_silver_standard$table %>%
  right_join(filtered_trait, by = c("phecode", "HPO", "EFO"))

pMap1 <- pGeneMap1 %>%
  select(trait, EFO, HPO, phecode) %>%
  distinct()

pGeneMap2 <- orphanet_based_silver_standard$table %>%
  right_join(filtered_trait, by = c("phecode", "HPO", "EFO"))

pMap2 <- pGeneMap2 %>%
  select(trait, EFO, HPO, phecode) %>%
  distinct()

for (phen in phens) {
  block_total1 <- focus_gene %>%
    ungroup() %>%
    filter(PHEN %in% phen) %>%
    mutate(KEEP = gene %in% pGeneMap1$gene) %>%
    group_by(PHEN, BLOCK) %>%
    mutate(ALL = any(KEEP)) %>%
    filter(ALL)
  block_total2 <- focus_gene %>%
    ungroup() %>%
    filter(PHEN %in% phen) %>%
    mutate(KEEP = gene %in% pGeneMap2$gene) %>%
    group_by(PHEN, BLOCK) %>%
    mutate(ALL = any(KEEP)) %>%
    filter(ALL)
  
  # for omim
  for (bb in unique(block_total1$BLOCK)) {
    tmp <- block_total1 %>%
      filter(BLOCK %in% bb) %>%
      ungroup() %>%
      select(gene, PIP.ME, PIP.meta) %>%
      crossing(trait = unique(filtered_trait$trait)) %>%
      as.data.frame()
    if (nrow(tmp) < 3) next
    
    output1 = silver_standard_perf(
      score_table = tmp,
      map_table = pMap1,
      silver_standard = omim_based_silver_standard,
      trait_codes = c("HPO", "phecode", "EFO"))
    
    output1$roc_auc$block <- bb
    output1$roc_auc$phen <- phen
    output1$roc_auc$data <- "omim"
    res <- rbind(res, output1$roc_auc)
  }
  
  # for orphan
  for (bb in unique(block_total2$BLOCK)) {
    tmp <- block_total2 %>%
      filter(BLOCK %in% bb) %>%
      ungroup() %>%
      select(gene, PIP.ME, PIP.meta) %>%
      crossing(trait = unique(filtered_trait$trait)) %>%
      as.data.frame()
    if (nrow(tmp) < 3) next
    
    output1 = silver_standard_perf(
      score_table = tmp,
      map_table = pMap2,
      silver_standard = orphanet_based_silver_standard,
      trait_codes = c("HPO", "phecode", "EFO"))
    
    output1$roc_auc$block <- bb
    output1$roc_auc$phen <- phen
    # output1$roc_auc$trait <- tr
    output1$roc_auc$data <- "orphanet"
    res <- rbind(res, output1$roc_auc)
  }
}

ans <- res %>%
  pivot_wider(names_from = score, values_from = roc_auc)

# write_tsv(ans, "./data/silver.tsv", quote_escape = FALSE)



