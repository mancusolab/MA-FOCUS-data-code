library(tidyverse)
library(broom)

load("../data/focus_aapower38_2.RData")

load("../data/focus.RData")

load("../data/focus_all.RData")

# load("../data/focus_redu.RData")
load("../data/twas.RData")

focus_AApower38 <- filter(focus_AApower38, !is.na(ID))

genoa_all <- read_tsv("../data/genoa_her_total.tsv")

tt <- focus_AApower38 %>%
  filter(PIP.ME > 0.8) %>%
  select(ID, TWAS.Z.EA, TWAS.Z.AA) %>%
  inner_join(genoa_all %>%
      select(GENE, POP, VG) %>%
      pivot_wider(names_from = POP, values_from = VG),
    by = c("ID" = "GENE")) %>%
  mutate(tmpEA = TWAS.Z.EA/sqrt(ea),
    tmpAA = TWAS.Z.AA/sqrt(aa))

glance(lm(tmpEA ~ tmpAA, tt))
tidy(lm(TWAS.Z.EA ~ TWAS.Z.AA, tt))
glance(lm(TWAS.Z.EA ~ TWAS.Z.AA, tt))



ggplot(tt, aes(x = ))

length(unique(focus_AApower38$PHEN))
length(unique(focus_AApower38$ID))
length(unique(focus_AApower38$BLOCK))

focus_AApower38 %>%
  distinct(PHEN, BLOCK)

tmp <- focus_AApower38 %>%
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


length(unique(tmp$BLOCK))
length(unique(tmp$PHEN))
length(unique(tmp$ID))


# how many sig TWAS in the region
sigt <- twas_all %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/4579 & POP != "meta") %>%
  distinct(ID, PHEN, POP) %>%
  mutate(TWAS = 1) %>%
  ungroup()

focus_AApower38 %>%
  filter(!grepl("NULL", ID)) %>%
  select(ID, BLOCK, PHEN) %>%
  left_join(sigt, by = c("ID", "PHEN")) %>%
  group_by(BLOCK, PHEN) %>%
  summarize(n = sum(TWAS, na.rm = TRUE)) %>%
  ungroup() %>%
  summarize(n = mean(n))



# how manny regions
focus_AApower38 %>%
  filter(IN.CRED.SET.ME) %>%
  filter(!grepl("NULL", ID)) %>%
  group_by(PHEN) %>%
  distinct(BLOCK)

focus_AApower38 %>%
  filter(IN.CRED.SET.ME) %>%
  filter(!grepl("NULL", ID)) %>%
  group_by(PHEN) %>%
  distinct(BLOCK)

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



# average causal genes
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
  filter(cgs == 3)

tmp %>%
  filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") %>%
  group_by(PHEN, BLOCK) %>%
  summarize(cgs = n()) %>%
  filter(cgs <= 3)


# average max rank PIP 
tmp %>%
  # filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") %>%
  filter(rank == 1) %>%
  group_by(PHEN, BLOCK) %>%
  mutate(maxpip = PIP == max(PIP)) %>%
  filter(maxpip) %>%
  ungroup() %>%
  distinct(PIP) %>%
  summarize(maxx = mean(PIP),
    sd = sd(PIP))

tmp %>%
  # filter(!grepl("NULL", ID)) %>%
  filter(POP == "ME") 

focus_AApower38 %>% 
  group_by(BLOCK, PHEN) %>%
  arrange(desc(PIP.ME)) %>%
  mutate(rank = row_number()) %>%
  filter(IN.CRED.SET.ME & !IN.CRED.SET.EA)


focus_AApower38 %>% 
  group_by(BLOCK, PHEN) %>%
  arrange(desc(PIP.ME)) %>%
  mutate(rank = row_number()) %>%
  filter(rank == 1 & !IN.CRED.SET.EA)

focus_AApower38 %>% 
  group_by(BLOCK, PHEN) %>%
  arrange(desc(PIP.ME)) %>%
  mutate(rank = row_number()) %>%
  filter(rank == 1 & !IN.CRED.SET.AA)


# MA-FOCUS cred sets (already done)
b <- focus_AApower38 %>%
  filter(IN.CRED.SET.ME == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)))

mean(log(b$bf))

# Baseline cred sets
b <- focus_AApower38 %>%
  filter(IN.CRED.SET.meta == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)))

mean(log(b$bf))

# EUR-FOCUS cred sets
b <- focus_AApower38 %>%
  filter(IN.CRED.SET.EA == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)))

mean(log(b$bf))

# AFR-FOCUS cred sets
b <- focus_AApower38 %>%
  filter(IN.CRED.SET.AA == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)))

mean(log(b$bf))

# All results MA-FOCUS results [no cred set]
b <- focus_AApower38 %>%
  filter(!grepl("NULL", ID)) %>%
  # filter(IN.CRED.SET.meta == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)))

mean(log(b$bf))

# enrichment 
# method 1
load("../data/enrich_aapower_2.RData")
mapper <- read_csv("../enrich/DisGeNET_meta_categories.csv")

# cate <- c("Hematological measurement",
#   "Abnormality, disorder, or disease of the immune system",
#   "Abnormality or disease of blood and blood-forming tissues", 
#   "Abnormality or disorder of the cardiovascular system")

cate <- c("Hematological measurement")

catelevel <- c("Hematological\nmeasurement",
  "Immune\nsystem",
  "Blood\ntissues", 
  "Cardiovascular\nsystem")

enrichment.results.df %>%
  left_join(mapper %>%
      select(Term, `meta_category`),
    by = "Term") %>%
  filter(!is.na(`meta_category`)) %>%
  filter(`meta_category` %in% cate) %>%
  # group_by(TEST, PHENO) %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  group_by(TEST) %>%
  summarize(n = n())

enrichment.results.df %>%
  left_join(mapper %>%
      select(Term, `meta_category`),
    by = "Term") %>%
  filter(!is.na(`meta_category`)) %>%
  filter(`meta_category` %in% cate) %>%
  # group_by(TEST, PHENO) %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  group_by(TEST) %>%
  summarize(X2 = -2*sum(log(Adjusted.P.value)), dof = 2*n(), P = pchisq(X2, dof, lower.tail=FALSE))

# enrDD <- enrichment.results.df %>%
#   left_join(mapper %>%
#       select(Term, `meta_category`),
#     by = "Term") %>%
#   filter(!is.na(`meta_category`)) %>%
#   filter(`meta_category` %in% cate)#  %>%
#   # group_by(TEST, PHENO) %>%
#   #filter(`Adjusted.P.value` < 0.05/60) 
# 
# enrDD %>% group_by(TEST) %>% summarize(mean(-log10(`Adjusted.P.value`)))
# enrDD %>% group_by(TEST) %>% summarize(N=n(), mean(-log10(`Adjusted.P.value`)))
# library(tidyverse)
# library(ggbeeswarm)
# ggplot(enrDD, aes(x=TEST, y = -log10(`Adjusted.P.value`)))  + geom_beeswarm()

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


dd <- read_tsv("../table/38_by_block_together2.tsv")
mean(dd$PIP.ME)
mean(dd$PIP.meta)



focus_all <- as.data.frame(focus_AApower38)
# ;rm(focus_AApower)

# Negate function
`%ni%` <- Negate(`%in%`)

# PIP locus zoom plot function
plot.locus.PIPs <- function(dat, tests, block, pheno, highlight.gene, colours=c("#e3d691BF", "#5ea9d4BF", "#d68592BF","#34ada1BF"), ...) {
  dat2 <- dat[dat$BLOCK==block & dat$PHEN==pheno,]
  ## Determine x-axis positions
  x.pos <- apply(dat2[,c(grep("P0", colnames(dat2)), grep("P1", colnames(dat2)))], MARGIN=1, mean)
  names(x.pos) <- dat2$ID; x.pos <- sort(x.pos)
  # Plot null model at x=1
  x.pos[1] <-1
  # Start plotting genes at x=2
  x.pos[-1] <- round(x.pos[-1]-x.pos[2]+2)
  # Adjust spacing so total plot is 12 units wide
  x.dist.adj <- as.numeric(max(x.pos) - x.pos[2])/10
  x.pos[-(1:2)] <- signif(x.pos[-(1:2)]/x.dist.adj + x.pos[2],3)
  ## Make plot
  plot(x=x.pos, y=rep(0, length(x.pos)), ylim=c(0,1), type="n", xlab="", ylab="Posterior Inclusion Probability", xaxt="n", yaxt="n", ...)
  abline(v=1.5, lwd=3)
  axis(2, las=2)
  axis(1, las=2, labels=names(x.pos[-match(highlight.gene, names(x.pos))]), at=x.pos[-match(highlight.gene, names(x.pos))], cex.axis=0.8)
  axis(1, las=2, labels=highlight.gene, at=x.pos[highlight.gene], cex.axis=0.8, col.axis ="red", font=2)
  for (test in tests) {
    idx <- match(test, tests)
    points(x.pos, dat2[match(names(x.pos),dat2$ID), paste0("PIP.",test)], pch=16, col=colours[idx], cex=1.4)
  }
}

# P-value locus zoom plot function
plot.locus.Pvals <- function(dat, tests, block, pheno, highlight.gene, colours=c("#e3d691BF", "#5ea9d4BF","#34ada1BF"), ...) {
  dat2 <- dat[dat$BLOCK==block & dat$PHEN==pheno,]
  ## Determine x-axis positions
  x.pos <- apply(dat2[,c(grep("P0", colnames(dat2)), grep("P1", colnames(dat2)))], MARGIN=1, mean)
  names(x.pos) <- dat2$ID; x.pos <- sort(x.pos)
  # Plot null model at x=1
  x.pos[1] <-1
  # Start plotting genes at x=2
  x.pos[-1] <- round(x.pos[-1]-x.pos[2]+2)
  # Adjust spacing so total plot is 12 units wide
  x.dist.adj <- as.numeric(max(x.pos) - x.pos[2])/10
  x.pos[-(1:2)] <- signif(x.pos[-(1:2)]/x.dist.adj + x.pos[2],3)
  ## Make plot
  p.vals <- -log10(na.omit(unlist(dat2[,grep("P.val", colnames(dat2))])))
  y.lim <- c(0, max(p.vals[is.finite(p.vals)]))
  y.lim[2] <- y.lim[2] * 1.01
  plot(x=x.pos, y=rep(0, length(x.pos)), ylim=y.lim, type="n", xlab="", ylab="-log10 P-value", xaxt="n", yaxt="n", ...)
  abline(v=1.5, lwd=3)
  axis(2, las=2)
  axis(1, las=2, labels=names(x.pos[-match(c("NULL",highlight.gene), names(x.pos))]), at=x.pos[-match(c("NULL",highlight.gene), names(x.pos))], cex.axis=0.8)
  axis(1, las=2, labels=highlight.gene, at=x.pos[highlight.gene], cex.axis=0.8, col.axis ="red", font=2)
  for (test in tests) {
    idx <- match(test, tests)
    points(x.pos[-1], -log10(dat2[match(names(x.pos[-1]),dat2$ID), paste0("TWAS.P.val.", test)]), pch=16, col=colours[idx], cex=1.4)
  }
}

### Analysis ----
# Find genes have high PIPs in one method, but are not in the CG set of the other
PIP.thresh <- 0.75
ME.specific.genes <- focus_all[focus_all$IN.CRED.SET.ME==T & focus_all$PIP.ME>PIP.thresh & focus_all$IN.CRED.SET.meta==F & focus_all$ID!="NULL",]
meta.specific.genes <- focus_all[focus_all$IN.CRED.SET.meta==T & focus_all$PIP.meta>PIP.thresh & focus_all$IN.CRED.SET.ME==F & focus_all$ID!="NULL",]

### Figure S19 ----
# plot EA and AA PIP distributions for each set of genes
par(mfrow=c(1,2))
for (test in c("ME","meta")) {
  dat <- get(paste0(test,".specific.genes"))
  if (test=="ME") {name="MA-FOCUS"} else if (test=="meta") {name="baseline"}
  plot(dat$PIP.EA, dat$PIP.AA, xlim=c(0,1), ylim=c(0,1), pch=16, cex=2, col=rgb(0,0,1,0.3), xlab="EA FOCUS PIP", ylab="AA FOCUS PIP", main=name)
}

### Figure S20 ----
# make locus zoom plots
par(mfcol=c(6,2))
for (i in 1:nrow(ME.specific.genes)) {
  block <- ME.specific.genes[i,"BLOCK"]
  pheno <- ME.specific.genes[i,"PHEN"]
  highlight.genes <- ME.specific.genes[ME.specific.genes$BLOCK==block & ME.specific.genes$PHEN==pheno, "ID"]
  plot.locus.PIPs(focus_all, block=block, pheno=pheno, tests=c("EA","AA","ME","meta"), highlight.gene=highlight.genes, main=paste0(pheno,", idx=",i," - MA FOCUS"))
  plot.locus.Pvals(focus_all, block=block, pheno=pheno, tests=c("EA","AA","meta"), highlight.gene=highlight.genes)
}
