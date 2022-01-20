library(tidyverse)

load("../data/focus_aapower.RData")
load("../data/twas.RData")

focus_AApower <- filter(focus_AApower, !is.na(ID))
length(unique(focus_AApower$ID))
length(unique(focus_AApower$BLOCK))

focus_AApower %>%
  distinct(PHEN, BLOCK)

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
  
length(unique(tmp$BLOCK))

# how manny regions
focus_AApower %>%
  filter(IN.CRED.SET.ME) %>%
  filter(!grepl("NULL", ID)) %>%
  group_by(PHEN) %>%
  distinct(BLOCK)

# how many sig TWAS in the region
sigt <- twas_all %>%
  group_by(PHEN, POP) %>%
  filter(TWAS.P < 0.05/4579 & POP != "meta") %>%
  distinct(ID, PHEN, POP) %>%
  mutate(TWAS = 1) %>%
  ungroup()

focus_AApower %>%
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


focus_AApower %>% 
  group_by(BLOCK, PHEN) %>%
  arrange(desc(PIP.ME)) %>%
  mutate(rank = row_number()) %>%
  filter(rank == 1 & !IN.CRED.SET.EA)

focus_AApower %>% 
  group_by(BLOCK, PHEN) %>%
  arrange(desc(PIP.ME)) %>%
  mutate(rank = row_number()) %>%
  filter(rank == 1 & !IN.CRED.SET.AA)


b <- focus_AApower %>%
  filter(IN.CRED.SET.ME == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)))

mean(log(b$bf))



# enrichment 
# method 1
load("../data/enrich_aapower.RData")
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

dd <- read_tsv("by_block_together2.tsv")
mean(dd$PIP.ME)
mean(dd$PIP.meta)

### end

