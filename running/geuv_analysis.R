library(tidyverse)

# for (pop in c("EUR", "YRI")) {
#   ff <- list.files("/project/nmancuso_8/data/GEUVADIS/GeneExpression",
#     pattern = paste0(pop, "\\..+profile"), full.names = TRUE)
#   
#   gene_name <- gsub(paste0("/project/nmancuso_8/data/GEUVADIS/GeneExpression/", pop, "\\."),
#     "", gsub("\\.predicted\\.gene\\.expression\\.profile", "", ff[1]))
#   tmp <- read_table2(ff[1], col_types = cols()) %>%
#     select(IID, SCORE)
#   colnames(tmp) <- c("IID", gene_name)
#   pred <- tmp
#   for (i in ff[2:length(ff)]) {
#     gene_name <- gsub(paste0("/project/nmancuso_8/data/GEUVADIS/GeneExpression/", pop, "\\."),
#       "", gsub("\\.predicted\\.gene\\.expression\\.profile", "", i))
#     tmp <- read_table2(i, col_types = cols()) %>%
#       select(IID, SCORE)
#     colnames(tmp) <- c("IID", gene_name)
#     pred <- pred %>%
#       left_join(tmp, by = "IID")
#   }
#   write_tsv(pred, file = paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/PRED_", pop, ".tsv"),
#     quote_escape = FALSE)
# }


# Run Switch
# for (pop in c("EUR", "YRI")) {
#   ff <- list.files("/project/nmancuso_8/data/GEUVADIS/GeneExpressionSwitch",
#     pattern = paste0(pop, "\\..+profile"), full.names = TRUE)
# 
#   gene_name <- gsub(paste0("/project/nmancuso_8/data/GEUVADIS/GeneExpressionSwitch/", pop, "\\."),
#     "", gsub("\\.predicted\\.gene\\.expression\\.switch\\.profile", "", ff[1]))
#   tmp <- read_table2(ff[1], col_types = cols()) %>%
#     select(IID, SCORE)
#   colnames(tmp) <- c("IID", gene_name)
#   pred <- tmp
#   for (i in ff[2:length(ff)]) {
#     gene_name <- gsub(paste0("/project/nmancuso_8/data/GEUVADIS/GeneExpressionSwitch/", pop, "\\."),
#       "", gsub("\\.predicted\\.gene\\.expression\\.switch\\.profile", "", i))
#     tmp <- read_table2(i, col_types = cols()) %>%
#       select(IID, SCORE)
#     colnames(tmp) <- c("IID", gene_name)
#     pred <- pred %>%
#       left_join(tmp, by = "IID")
#   }
#   write_tsv(pred, file = paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/PRED_", pop, "_switch.tsv"),
#     quote_escape = FALSE)
# }



for (pop in c("EUR", "YRI")) {
  print(paste0("Running ", pop))
  
  pt <- read_table2(paste0("/project/nmancuso_8/data/GEUVADIS/GenotypeData/GRM/rm_related/GEUVADIS.GRM.all.rm125.",
    pop, ".grm.id"), col_names = FALSE)
  
  pred <- read_tsv(paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/PRED_", pop, ".tsv"))
  pred <- pred %>%
    filter(IID %in% pt$X1)
  gene <- colnames(pred)[colnames(pred) != "IID"]
  
  mapper <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/GENOA.gene.list.tsv",
    col_names = FALSE) %>%
    filter(X1 %in% gene) %>%
    group_by(X1) %>%
    filter(n() == 1)
  gene <- gene[gene %in% mapper$X1]
  real <- read_tsv(paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/PHEN_", pop, "_rm125.tsv"))
  real <- real %>%
    filter(name %in% pt$X1)
  real <- real[c("name", mapper$X2)]
  colnames(real) <- c("badname", mapper$X1)
  real <- real %>%
    pivot_longer(col = !contains("badname"), names_repair = "unique")
  colnames(real) <- c("IID", "gene", "true")
  colnames(pred)[1] <- "badname"
  pred <- pred %>%
    pivot_longer(col = !contains("badname"), names_repair = "unique") %>%
    filter(name %in% gene)
  colnames(pred) <- c("IID", "gene", "pred")
  
  if (nrow(real) == nrow(pred)) {
    print("==========ROW NUMBER MATCHES!============")
  }
  
  res <- pred %>%
    left_join(real, by = c("IID", "gene")) %>%
    group_by(gene) %>%
    summarize(r2 = cor(pred, true)^2)
  
  write_tsv(res, file = paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/", pop, "_r2.tsv"), quote_escape = FALSE)
  rm(pred)
  rm(real)
  ff <- list.files(paste0("/project/nmancuso_8/data/GEUVADIS/Heritability/", pop),
    full.names = TRUE)
  he_ff <- ff[grepl("\\.HE\\.", ff)]
  reml_ff <- ff[grepl("\\.REML\\.", ff)]
  he_ff <- he_ff[order(he_ff)]
  reml_ff <- reml_ff[order(reml_ff)]
  
  if (length(he_ff[!is.na(he_ff)]) == length(reml_ff[!is.na(reml_ff)])){
    print("ID length match!")
  }
  
  heri <- tibble()
  for (i in 1:length(he_ff)) {
    he_id <- gsub(paste0("\\.", pop, "\\.hsq"), "",
      gsub(paste0("/project/nmancuso_8/data/GEUVADIS/Heritability/", pop, "/GEUVADIS\\.HE\\."), "", he_ff[i]))
    reml_id <- gsub(paste0("\\.", pop, "\\.hsq"), "",
      gsub(paste0("/project/nmancuso_8/data/GEUVADIS/Heritability/", pop, "/GEUVADIS\\.REML\\."), "", reml_ff[i]))
    if (he_id != reml_id) {
      stop("id doesn't match!")
    }
    
    tmp_gene <- mapper[mapper$X2 == he_id, ]$X1
    if(length(tmp_gene) == 0) next
    tmp <- read_table2(he_ff[i], col_types = cols()) %>%
      mutate(gene = tmp_gene) %>%
      select(gene, TEST, H2) %>%
      bind_rows(read_table2(reml_ff[i], col_types = cols()) %>%
          mutate(TEST="REML",
            gene = tmp_gene) %>%
          select(gene, TEST, H2)) %>%
      pivot_wider(names_from = TEST, values_from = H2)
    heri <- bind_rows(heri, tmp)
  }
  write_tsv(heri, file = paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/", pop, "_heri.tsv"),
    quote_escape = FALSE)
}


# Switch

for (pop in c("EUR", "YRI")) {
  print(paste0("Running ", pop))
  
  pt <- read_table2(paste0("/project/nmancuso_8/data/GEUVADIS/GenotypeData/GRM/rm_related/GEUVADIS.GRM.all.rm125.",
    pop, ".grm.id"), col_names = FALSE)
  
  pred <- read_tsv(paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/PRED_", pop, "_switch.tsv"))
  pred <- pred %>%
    filter(IID %in% pt$X1)
  gene <- colnames(pred)[colnames(pred) != "IID"]
  
  mapper <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/GENOA.gene.list.tsv",
    col_names = FALSE) %>%
    filter(X1 %in% gene) %>%
    group_by(X1) %>%
    filter(n() == 1)
  gene <- gene[gene %in% mapper$X1]
  real <- read_tsv(paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/PHEN_", pop, "_rm125.tsv"))
  real <- real %>%
    filter(name %in% pt$X1)
  real <- real[c("name", mapper$X2)]
  colnames(real) <- c("badname", mapper$X1)
  real <- real %>%
    pivot_longer(col = !contains("badname"), names_repair = "unique")
  colnames(real) <- c("IID", "gene", "true")
  colnames(pred)[1] <- "badname"
  pred <- pred %>%
    pivot_longer(col = !contains("badname"), names_repair = "unique") %>%
    filter(name %in% gene)
  colnames(pred) <- c("IID", "gene", "pred")
  
  if (nrow(real) == nrow(pred)) {
    print("==========ROW NUMBER MATCHES!============")
  }
  
  res <- pred %>%
    left_join(real, by = c("IID", "gene")) %>%
    group_by(gene) %>%
    summarize(r2 = cor(pred, true)^2) %>%
    filter(!is.na(r2))
  
  write_tsv(res, file = paste0("/project/nmancuso_8/data/GEUVADIS/ProcessedData/", pop, "_r2_switch.tsv"), quote_escape = FALSE)
  
}




# Genoa Switch


# Run switch

# AA
ff <- list.files("/project/nmancuso_8/data/MEFOCUS/pred", pattern = "AA\\..+\\.genoa\\..+profile", full.names = TRUE)

gene_name <- gsub(paste0("/project/nmancuso_8/data/MEFOCUS/pred/AA\\."),
  "", gsub("\\.genoa\\.predicted\\.gene\\.expression\\.switch\\.profile", "", ff[1]))

tmp <- read_table2(ff[1], col_types = cols()) %>%
  select(IID, SCORE)
colnames(tmp) <- c("IID", gene_name)
pred <- tmp
for (i in ff[2:length(ff)]) {
  gene_name <- gsub(paste0("/project/nmancuso_8/data/MEFOCUS/pred/AA\\."),
    "", gsub("\\.genoa\\.predicted\\.gene\\.expression\\.switch\\.profile", "", i))
  tmp <- read_table2(i, col_types = cols()) %>%
    select(IID, SCORE)
  colnames(tmp) <- c("IID", gene_name)
  pred <- pred %>%
    left_join(tmp, by = "IID")
}
write_tsv(pred, file = "/project/nmancuso_8/data/MEFOCUS/genoa_pred_AA_switch.tsv", quote_escape = FALSE)


# EA

ff <- list.files("/project/nmancuso_8/data/MEFOCUS/pred", pattern = "EA\\..+\\.genoa\\..+profile", full.names = TRUE)

gene_name <- gsub(paste0("/project/nmancuso_8/data/MEFOCUS/pred/EA\\."),
  "", gsub("\\.genoa\\.predicted\\.gene\\.expression\\.switch\\.profile", "", ff[1]))

tmp <- read_table2(ff[1], col_types = cols()) %>%
  select(IID, SCORE)
colnames(tmp) <- c("IID", gene_name)
pred <- tmp
for (i in ff[2:length(ff)]) {
  gene_name <- gsub(paste0("/project/nmancuso_8/data/MEFOCUS/pred/EA\\."),
    "", gsub("\\.genoa\\.predicted\\.gene\\.expression\\.switch\\.profile", "", i))
  tmp <- read_table2(i, col_types = cols()) %>%
    select(IID, SCORE)
  colnames(tmp) <- c("IID", gene_name)
  pred <- pred %>%
    left_join(tmp, by = "IID")
}
write_tsv(pred, file = "/project/nmancuso_8/data/MEFOCUS/genoa_pred_EA_switch.tsv", quote_escape = FALSE)

# GENOA
# AA

pt <- read_table2("/project/nmancuso_8/data/GENOA/FUSION/aa_441_pt.id", col_names = FALSE)

pred <- read_tsv("/project/nmancuso_8/data/MEFOCUS/genoa_pred_AA_switch.tsv")
pred <- pred %>%
  filter(IID %in% pt$X1)
gene <- colnames(pred)[colnames(pred) != "IID"]

mapper <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/GENOA.gene.list.tsv",
  col_names = FALSE) %>%
  filter(X1 %in% gene) %>%
  group_by(X1) %>%
  filter(n() == 1)
gene <- gene[gene %in% mapper$X1]
real <- read_tsv("/project/nmancuso_8/data/GENOA/FUSION/expr/aa/aa_GRCh38_expr.bed.gz")
real <- real[c("gene", pt$X1)]
real <- real %>%
  filter(gene %in% mapper$X1)

real <- real %>%
  pivot_longer(cols = -gene)
pred <- pred %>%
  pivot_longer(cols = -IID)
colnames(real) <- c("gene", "iid", "true")
colnames(pred) <- c("iid", "gene", "pred")

res <- pred %>%
  inner_join(real, by = c("iid", "gene")) %>%
  group_by(gene) %>%
  summarize(r2 = cor(pred, true)^2) %>%
  filter(!is.na(r2))

write_tsv(res, file = "/project/nmancuso_8/data/MEFOCUS/genoa_AA_switch_r2.tsv",
  quote_escape = FALSE)

# EA

pt <- read_table2("/project/nmancuso_8/data/GENOA/FUSION/ea_373_pt.id", col_names = FALSE)

pred <- read_tsv("/project/nmancuso_8/data/MEFOCUS/genoa_pred_EA_switch.tsv")
pred <- pred %>%
  filter(IID %in% pt$X1)
gene <- colnames(pred)[colnames(pred) != "IID"]

mapper <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/GENOA.gene.list.tsv",
  col_names = FALSE) %>%
  filter(X1 %in% gene) %>%
  group_by(X1) %>%
  filter(n() == 1)
gene <- gene[gene %in% mapper$X1]
real <- read_tsv("/project/nmancuso_8/data/GENOA/FUSION/expr/ea/ea_GRCh38_expr.bed.gz")
real <- real[c("gene", pt$X1)]
real <- real %>%
  filter(gene %in% mapper$X1)

real <- real %>%
  pivot_longer(cols = -gene)
pred <- pred %>%
  pivot_longer(cols = -IID)
colnames(real) <- c("gene", "iid", "true")
colnames(pred) <- c("iid", "gene", "pred")

res <- pred %>%
  inner_join(real, by = c("iid", "gene")) %>%
  group_by(gene) %>%
  summarize(r2 = cor(pred, true)^2) %>%
  filter(!is.na(r2))

write_tsv(res, file = "/project/nmancuso_8/data/MEFOCUS/genoa_EA_switch_r2.tsv",
  quote_escape = FALSE)


# aggregate

genoa_ea_switch <- read_tsv("/project/nmancuso_8/data/MEFOCUS/genoa_EA_switch_r2.tsv")
genoa_aa_switch <- read_tsv("/project/nmancuso_8/data/MEFOCUS/genoa_AA_switch_r2.tsv")
geuv_yri_switch <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/YRI_r2_switch.tsv")
geuv_eur_switch <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/EUR_r2_switch.tsv")
geuv_yri <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/YRI_r2.tsv")
geuv_eur <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/EUR_r2.tsv")

cvr2 <- read_tsv("/home1/zeyunlu/research/sub_mefocus/data/genoa_her_total.tsv")

tot <- genoa_ea_switch %>%
  select(gene, genoa_ea_s = r2) %>%
  full_join(genoa_aa_switch %>%
      select(gene, genoa_aa_s = r2),
    by = "gene") %>%
  full_join(geuv_yri_switch %>%
      select(gene, geuv_yri_s = r2),
    by = "gene") %>%
  full_join(geuv_eur_switch %>%
      select(gene, geuv_eur_s = r2),
    by = "gene") %>%
  full_join(geuv_eur %>%
      select(gene, geuv_eur = r2),
    by = "gene") %>%
  full_join(geuv_yri %>%
      select(gene, geuv_yri = r2),
    by = "gene") %>%
  full_join(cvr2 %>%
      select(gene = GENE, cvr2 = MODELCV.R2, pop = POP) %>%
      pivot_wider(names_from = pop, values_from = cvr2),
    by = "gene")

write_tsv(tot, "/home1/zeyunlu/research/sub_mefocus/data/total_r2.tsv", quote_escape = FALSE)


# add geuvadis -> genoa 

# AA
ff <- list.files("/project/nmancuso_8/data/MAFOCUS/pred/genoa_aa_geuvadis_yri", pattern = "AA\\..+\\.genoa\\..+profile", full.names = TRUE)

gene_name <- gsub(paste0("/project/nmancuso_8/data/MAFOCUS/pred/genoa_aa_geuvadis_yri/AA\\."),
  "", gsub("\\.genoa\\.predicted\\.gene\\.expression\\.geuvadis\\.yri\\.profile", "", ff[1]))

tmp <- read_table2(ff[1], col_types = cols()) %>%
  select(IID, SCORE)
colnames(tmp) <- c("IID", gene_name)
pred <- tmp
for (i in ff[2:length(ff)]) {
  gene_name <- gsub(paste0("/project/nmancuso_8/data/MAFOCUS/pred/genoa_aa_geuvadis_yri/AA\\."),
    "", gsub("\\.genoa\\.predicted\\.gene\\.expression\\.geuvadis\\.yri\\.profile", "", i))
  tmp <- read_table2(i, col_types = cols()) %>%
    select(IID, SCORE)
  colnames(tmp) <- c("IID", gene_name)
  pred <- pred %>%
    left_join(tmp, by = "IID")
}
write_tsv(pred, file = "/project/nmancuso_8/data/MAFOCUS/genoa_pred_AA_geuvadis_yri.tsv", quote_escape = FALSE)


# EA

ff <- list.files("/project/nmancuso_8/data/MAFOCUS/pred/genoa_ea_geuvadis_eur", pattern = "EA\\..+\\.genoa\\..+profile", full.names = TRUE)

gene_name <- gsub(paste0("/project/nmancuso_8/data/MAFOCUS/pred/genoa_ea_geuvadis_eur/EA\\."),
  "", gsub("\\.genoa\\.predicted\\.gene\\.expression\\.geuvadis\\.eur\\.profile", "", ff[1]))

tmp <- read_table2(ff[1], col_types = cols()) %>%
  select(IID, SCORE)
colnames(tmp) <- c("IID", gene_name)
pred <- tmp
for (i in ff[2:length(ff)]) {
  gene_name <- gsub(paste0("/project/nmancuso_8/data/MAFOCUS/pred/genoa_ea_geuvadis_eur/EA\\."),
    "", gsub("\\.genoa\\.predicted\\.gene\\.expression\\.geuvadis\\.eur\\.profile", "", i))
  tmp <- read_table2(i, col_types = cols()) %>%
    select(IID, SCORE)
  colnames(tmp) <- c("IID", gene_name)
  pred <- pred %>%
    left_join(tmp, by = "IID")
}
write_tsv(pred, file = "/project/nmancuso_8/data/MAFOCUS/genoa_pred_EA_geuvadis_eur.tsv", quote_escape = FALSE)



# AA
pt <- read_table2("/project/nmancuso_8/data/GENOA/FUSION/aa_441_pt.id", col_names = FALSE)

pred <- read_tsv("/project/nmancuso_8/data/MAFOCUS/genoa_pred_AA_geuvadis_yri.tsv")
pred <- pred %>%
  filter(IID %in% pt$X1)
gene <- colnames(pred)[colnames(pred) != "IID"]

mapper <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/GENOA.overlap.gene.list.tsv",
  col_names = FALSE) %>%
  filter(X1 %in% gene) %>%
  group_by(X1) %>%
  filter(n() == 1)
gene <- gene[gene %in% mapper$X1]
real <- read_tsv("/project/nmancuso_8/data/GENOA/FUSION/expr/aa/aa_GRCh38_expr.bed.gz")
real <- real[c("gene", pt$X1)]
real <- real %>%
  filter(gene %in% mapper$X1)

real <- real %>%
  pivot_longer(cols = -gene)
pred <- pred %>%
  pivot_longer(cols = -IID)
colnames(real) <- c("gene", "iid", "true")
colnames(pred) <- c("iid", "gene", "pred")

res <- pred %>%
  inner_join(real, by = c("iid", "gene")) %>%
  group_by(gene) %>%
  summarize(r2 = cor(pred, true)^2) %>%
  filter(!is.na(r2))

write_tsv(res, file = "/project/nmancuso_8/data/MAFOCUS/genoa_AA_geuvadis_yri_r2.tsv",
  quote_escape = FALSE)

# EA

pt <- read_table2("/project/nmancuso_8/data/GENOA/FUSION/ea_373_pt.id", col_names = FALSE)

pred <- read_tsv("/project/nmancuso_8/data/MAFOCUS/genoa_pred_EA_geuvadis_eur.tsv")
pred <- pred %>%
  filter(IID %in% pt$X1)
gene <- colnames(pred)[colnames(pred) != "IID"]

mapper <- read_tsv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/GENOA.overlap.gene.list.tsv",
  col_names = FALSE) %>%
  filter(X1 %in% gene) %>%
  group_by(X1) %>%
  filter(n() == 1)
gene <- gene[gene %in% mapper$X1]
real <- read_tsv("/project/nmancuso_8/data/GENOA/FUSION/expr/ea/ea_GRCh38_expr.bed.gz")
real <- real[c("gene", pt$X1)]
real <- real %>%
  filter(gene %in% mapper$X1)

real <- real %>%
  pivot_longer(cols = -gene)
pred <- pred %>%
  pivot_longer(cols = -IID)
colnames(real) <- c("gene", "iid", "true")
colnames(pred) <- c("iid", "gene", "pred")

res <- pred %>%
  inner_join(real, by = c("iid", "gene")) %>%
  group_by(gene) %>%
  summarize(r2 = cor(pred, true)^2) %>%
  filter(!is.na(r2))

write_tsv(res, file = "/project/nmancuso_8/data/MAFOCUS/genoa_EA_geuvadis_eur_r2.tsv",
  quote_escape = FALSE)

# aggre

dd <- read_tsv("/home1/zeyunlu/research/sub_mefocus/data/total_r2.tsv")
ea <- read_tsv("/project/nmancuso_8/data/MAFOCUS/genoa_EA_geuvadis_eur_r2.tsv")
aa <- read_tsv("/project/nmancuso_8/data/MAFOCUS/genoa_AA_geuvadis_yri_r2.tsv")

tot <- dd %>%
  full_join(ea %>%
      rename(genoa_ea = r2),
    by = "gene") %>%
  full_join(aa %>%
      rename(genoa_aa = r2),
    by = "gene")

write_tsv(tot, "/home1/zeyunlu/research/sub_mefocus/data/total_r2.tsv", quote_escape = FALSE)



theme_corr <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey70"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.title=element_text(size=x_axis_size,face="bold"))
}

eur <- tot %>%
  select(gene, ea, genoa_ea_s,  geuv_eur, geuv_eur_s) %>%
  pivot_longer(cols = c(genoa_ea_s,  geuv_eur, geuv_eur_s)) %>%
  mutate(name = factor(name, levels = c("genoa_ea_s", "geuv_eur", "geuv_eur_s"),
    labels = c("GENOA In-sample with AA Weights",
      "GEUVADIS Out-of-sample with EA Weights",
      "GEUVADIS Out-of-sample with AA Weights")))

eurpp <- ggplot(eur, aes(x = ea, y = value)) +
  geom_point() +
  facet_wrap(~name) +
  geom_smooth(method = "lm") +
  theme_corr() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ylab("R2 between measured and predicted") +
  xlab("GENOA EA CV R2")

yri <- tot %>%
  select(gene, aa, genoa_aa_s,  geuv_yri, geuv_yri_s) %>%
  pivot_longer(cols = c(genoa_aa_s,  geuv_yri, geuv_yri_s)) %>%
  mutate(name = factor(name, levels = c("genoa_aa_s", "geuv_yri", "geuv_yri_s"),
    labels = c("GENOA In-sample with EA Weights",
      "GEUVADIS Out-of-sample with AA Weights",
      "GEUVADIS Out-of-sample with EA Weights")))

yripp <- ggplot(yri, aes(x = aa, y = value)) +
  geom_point() +
  facet_wrap(~name) +
  geom_smooth(method = "lm") +
  theme_corr() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ylab("R2 between measured and predicted") +
  xlab("GENOA AA CV R2")





#####################
# 4x4
ssr2 <- read_tsv("../data/EUR_r2_switch.tsv") %>%
  mutate(POP = "EUR") %>%
  select(GENE = gene, POP, cr2 = r2) %>%
  bind_rows(read_tsv("../data/YRI_r2_switch.tsv") %>%
      mutate(POP = "YRI") %>%
      select(GENE = gene, POP, cr2 = r2))

crossPop <- geuvadisdd %>%
  select(GENE, POP, VG, r2) %>%
  left_join(ssr2, by = c("GENE", "POP"))

p41 <- ggplot(filter(geuvadisdd, POP %in% "EUR"), aes(x = VG, y = r2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  ylab(expression(bold("Prediction")~bolditalic(r)^{2})) +
  xlab(expression(bold("GEUVADIS EUR")~bolditalic(cis)~"-"~bolditalic(sigma)^{2})) +
  ylim(0, 1)

p42 <- ggplot(filter(crossPop, POP %in% "EUR"), aes(x = VG, y = cr2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  ylab(expression(bold("Prediction")~bolditalic(r)^{2}~bold("using YRI weights"))) +
  xlab(expression(bold("GEUVADIS EUR")~bolditalic(cis)~"-"~bolditalic(sigma)^{2})) +
  ylim(0, 1)

p43 <- ggplot(filter(geuvadisdd, POP %in% "YRI"), aes(x = VG, y = r2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  ylab(expression(bold("Prediction")~bolditalic(r)^{2})) +
  xlab(expression(bold("GEUVADIS YRI")~bolditalic(cis)~"-"~bolditalic(sigma)^{2})) +
  ylim(0, 1)

p44 <- ggplot(filter(crossPop, POP %in% "YRI"), aes(x = VG, y = cr2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  ylab(expression(bold("Prediction")~bolditalic(r)^{2}~bold("using EUR weights"))) +
  xlab(expression(bold("GEUVADIS YRI")~bolditalic(cis)~"-"~bolditalic(sigma)^{2})) +
  ylim(0, 1)

p4 <- ggarrange(p41, p42, p43, p44,
  labels = c("A", "B", "C", "D"))

ggsave("44.png", plot = p4, path = "../plot/", height = one_row_height, width = two_col)


p51 <- ggplot(filter(crossPop, POP %in% "EUR"), aes(x = r2, y = cr2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  xlab(expression(bold("EUR Prediction")~bolditalic(r)^{2}~bold("using EUR weights"))) +
  ylab(expression(bold("EUR Prediction")~bolditalic(r)^{2}~bold("using YRI weights"))) +
  ylim(0, 1)

p52 <- ggplot(filter(crossPop, POP %in% "YRI"), aes(x = r2, y = cr2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  xlab(expression(bold("YRI Prediction")~bolditalic(r)^{2}~bold("using YRI weights"))) +
  ylab(expression(bold("YRI Prediction")~bolditalic(r)^{2}~bold("using EUR weights"))) +
  ylim(0, 1)

p5 <- ggarrange(p51, p52,
  labels = c("A", "B"))

ggsave("22.png", plot = p5, path = "../plot/", height = one_row_height, width = two_col)


#####################



