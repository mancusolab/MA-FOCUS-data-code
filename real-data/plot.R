rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggnewscale)
library(ggthemes)
library(cowplot)
library(ggrepel)

source("figure_style.R")
mean_size = 2
mean_shape = 22
median_size = 2
median_shape = 22
x_axis_size = 8
y_axis_size = 8

point_size = 3
point_shape = 21
errbar_width = 0.5

arrange_label_size = 8
one_row_height = 2.5
two_row_height = 4.5
three_row_height = 6
one_col = 3.35
onehalf_col = 4.49
two_col = 6.85
PIP.threshold <- 0.9

ylabel_PIP <- "Mean of credible set PIP"
ylabel_CGS <- "Credible set size"
ylabel_SD <- "SD of credible set PIP"
ffont <- "sans"
fontsize <- 8

phens <- c("RBC","RDW","WBC","PLT","MPV","LYM","NEU","MON","BAS","EOS","HGB","HCT","MCV","MCH","MCHC")
phens <- phens[order(phens)]

theme_corr <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(size = 0.25, color = "grey70"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    axis.title=element_text(face="bold"),
    text=element_text(size = fontsize,  family = ffont))
}

# shared architecture: correlation between pop
# plot m4 and s12
genoadd <- read_tsv("./data/genoa_heritability.tsv")

genoadd1 <- genoadd %>%
  select(GENE, POP, HSQ) %>%
  pivot_wider(names_from=POP, values_from = HSQ)

p41 <- ggplot(genoadd1, aes(x = ea, y = aa)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  scale_x_continuous(name = expression(bold("GENOA EA")~bolditalic(cis)~"-"~bolditalic(h)[g]^{2}), limits = c(0, 1)) +
  scale_y_continuous(name = expression(bold("GENOA AA")~bolditalic(cis)~"-"~bolditalic(h)[g]^{2}), limits = c(0, 1))

ps121 <- ggplot(filter(genoadd, POP %in% "ea"),
  aes(x = MODELCV.R2, y = HSQ)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  scale_x_continuous(name = expression(bold("GENOA EA CV")~bolditalic(r)^{2}), limits = c(0, 1)) +
  scale_y_continuous(name = expression(bold("GENOA EA")~bolditalic(cis)~"-"~bolditalic(h)[g]^{2}), limits = c(0, 1))

ps122 <- ggplot(filter(genoadd, POP %in% "aa"),
  aes(x = MODELCV.R2, y = HSQ)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  # stat_cor(method = "pearson", label.x = 0.01, label.y = 0.9) +
  scale_x_continuous(name = expression(bold("GENOA AA CV")~bolditalic(r)^{2}), limits = c(0, 1)) +
  scale_y_continuous(name = expression(bold("GENOA AA")~bolditalic(cis)~"-"~bolditalic(h)[g]^{2}), limits = c(0, 1))

geuvadisdd <- read_tsv("./data/geuvadis_heritability.tsv")

p43 <- ggplot(filter(geuvadisdd, POP %in% "EUR"), aes(x = HSQ, y = r2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  ylab(expression(bold("Prediction")~bolditalic(r)^{2})) +
  xlab(expression(bold("GEUVADIS EUR")~bolditalic(cis)~"-"~bolditalic(h)[g]^{2})) +
  ylim(0, 1)

p44 <- ggplot(filter(geuvadisdd, POP %in% "YRI"), aes(x = HSQ, y = r2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_corr() +
  theme(axis.title.x=element_text(size=x_axis_size,face="bold")) +
  ylab(expression(bold("Prediction")~bolditalic(r)^{2})) +
  xlab(expression(bold("GEUVADIS YRI")~bolditalic(cis)~"-"~bolditalic(h)[g]^{2})) +
  ylim(0, 1)

p4 <- ggarrange(p41, p43, p44,
  labels = c("A", "B", "C"), 
  font.label = list(size = arrange_label_size, family = ffont), nrow = 1)

ggsave("figure-m4.pdf", plot = p4, path = "./plot/", height = one_row_height, width = two_col,
  dpi = 300)

ps12 <- ggarrange(ps121, ps122,
  labels = c("A", "B"), nrow = 2)
ggsave("figure-s12.png", plot = ps12, path = "./plot/", height = two_row_height, width = onehalf_col)

# Figure s13 (admixture analysis plot) is drew used pong. See script XXXX

# comparison r2, s14
tot <- read_tsv("./data/total_r2.tsv")

tmptot <- tot %>%
  pivot_longer(-gene) %>%
  mutate(group = factor(ifelse(grepl("(ea)|(eur)", name), "A: EA/EUR participants", "B: AA/YRI participants"),
    levels = c("A: EA/EUR participants", "B: AA/YRI participants")),
    cVar = ifelse(name %in% c("ea", "aa"), "A",
      ifelse(name %in% c("genoa_ea_s", "genoa_aa_s"), "B",
        ifelse(name %in% c("geuv_eur", "geuv_yri"), "C", "D"))))


ps14 <- ggplot(tmptot, aes(x = cVar, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(c(0, 0.3)) +
  facet_wrap(~group, nrow = 2) +
  ylab(expression(bolditalic(r)^{2})) +
  scale_x_discrete(labels = c("A" = expression(bold("Pred. CV")),
    "B" = expression(bold("GENOA Cross-Pop")),
    "C" = expression(bold("GEUVADIS Similar-Pop")),
    "D" = expression(bold("GEUVADIS Cross-Pop")))) +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey70"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title=element_text(size=x_axis_size,face="bold"))

ggsave("figure-s14.png", plot = ps14, path = "./plot/", height = two_row_height, width = two_col)


 load("./data/focus_analysis.RData")
load("./data/twas.RData")

theme_ma <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(size = 0.5, color = "grey70"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    axis.title=element_text(face="bold"),
    text=element_text(size = fontsize,  family = ffont))
}

# TWAS normalized statistics  s15
#  aggregated
a <- twas_all %>%
  filter(!POP %in% "meta") %>%
  select(PHEN, ID, TWAS.Z.NORM, POP) %>%
  pivot_wider(names_from = POP, values_from = TWAS.Z.NORM) %>%
  mutate(PHEN = factor(PHEN, levels = phens))

ps151 <- ggplot(a, aes(x = EA, y = AA)) +
  geom_point() +
  stat_cor(method = "pearson", cor.coef.name = c("r", "rho", "tau")) +
  geom_smooth(method = "lm") +
  theme_ma() +
  xlab("EA Norm. TWAS Z") +
  ylab("AA Norm. TWAS Z")

# by trait

ps152 <- ggplot(a, aes(x = EA, y = AA)) +
  geom_point() +
  stat_cor(method = "pearson", cor.coef.name = c("r", "rho", "tau")) +
  geom_smooth(method = "lm") +
  facet_wrap(~PHEN, nrow = 5) +
  theme_ma() +
  xlab("EA Norm. TWAS Z") +
  ylab("AA Norm. TWAS Z")


ps15 <- ggarrange(ps151, ps152, nrow = 2, labels = c("A", "B"), heights = c(1.3, 3),
  font.label = list(size = arrange_label_size))
ggsave("./figure-s15.png", plot = ps15, path = "./plot/", height = two_row_height*2, width = two_col)


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


# average PIP  s16
a <- tmp %>%
  # filter(POP %in% "ME") %>%
  # filter(!grepl("NULL", ID)) %>%
  select(PHEN, ID, PIP, rank) %>%
  group_by(PHEN, rank, PIP, POP) %>%
  distinct() %>%
  group_by(rank, POP) %>%
  summarize(pip = mean(PIP),
    n = n(),
    se = sd(PIP)/sqrt(n))

ps16 <- ggplot(a, aes(x = rank, y = pip, group = POP, color = POP)) +
  geom_point() +
  scale_color_manual(values = COLS3, labels = LABELS3) +
  geom_errorbar(aes(ymin = pip - se, ymax = pip + se), width = 0.5) +
  ylim(0, 1) +
  xlab("Rank of model PIP in 90% CS") +
  ylab("PIP in 90% CS") +
  scale_x_continuous(breaks = seq(1,9), label = seq(1, 9)) +
  theme_ma()


ggsave("figure-s16.png", plot = ps16, path = "./plot/", height = one_row_height, width = two_col)


# making PIP plots m6
theme_pip <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(size = 0.5, color = "grey70"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title=element_text(face="bold"),
    text=element_text(size = fontsize,  family = ffont))
}

mA <- tmp %>%
  filter(!grepl("NULL", ID)) %>%
  group_by(PHEN, BLOCK, POP) %>%
  summarize(pip = mean(PIP)) %>%
  mutate(PHEN = factor(PHEN, levels = phens),
    POP = factor(POP, levels = c("EA", "AA", "ME", "meta"),
      labels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip")))

p1 <- ggplot(mA, aes(x = POP, y = pip, fill = POP, group = POP)) +
  geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
  scale_fill_manual(values = COLS2, labels = LABELS2) +
  geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.5) +
  theme_pip() +
  new_scale_fill() +
  stat_summary(fun = mean, geom = "point",  aes(fill = "Mean"), size = mean_size, shape = mean_shape,
    position = position_dodge(width = 0.3)) +
  stat_summary(fun = median, geom = "point", aes(fill = "Median"), shape = median_shape, size = median_size,
    position = position_dodge(width = 0.3)) +
  scale_fill_manual(values = c("Mean"="white", "Median"="black")) +
  scale_y_continuous(name = ylabel_PIP)


# sd
mB <- tmp %>%
  filter(!grepl("NULL", ID)) %>%
  group_by(PHEN, BLOCK, POP) %>%
  summarize(sdPIP = sd(PIP)) %>%
  mutate(PHEN = factor(PHEN, levels = phens),
    POP = factor(POP, levels = c("EA", "AA", "ME", "meta"),
      labels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip")))

p2 <- ggplot(mB, aes(x = POP, y = sdPIP, fill = POP, group = POP)) +
  scale_fill_manual(values = COLS2, labels = LABELS2) +
  geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
  theme_pip() +
  new_scale_fill() +
  stat_summary(fun = mean, geom = "point",  aes(fill = "Mean"), size = mean_size, shape = mean_shape,
    position = position_dodge(width = 0.3)) +
  stat_summary(fun = median, geom = "point", aes(fill = "Median"), shape = median_shape, size = median_size,
    position = position_dodge(width = 0.3)) +
  scale_fill_manual(values = c("Mean"="white", "Median"="black")) +
  scale_y_continuous(name = ylabel_SD)

# size
mC1 <- tmp %>%
  filter(!grepl("NULL", ID)) %>%
  group_by(PHEN, BLOCK, POP) %>%
  summarize(cgs = n()) %>%
  group_by(POP, cgs) %>%
  summarize(cc = n()) %>%
  mutate(POP = factor(POP, levels = c("EA", "AA", "ME", "meta"),
    labels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip")))

p3 <- ggplot(mC1, aes(x = cgs, y = cc, fill = POP, group = POP)) +
  scale_fill_manual(values = COLS2, labels = LABELS2) +
  geom_col(alpha = 0.8, position=position_dodge2(preserve = "single"), color = "black") +
  theme_corr() +
  new_scale_fill() +
  scale_y_continuous(name = "Count") +
  scale_x_continuous(name = "90% Credible set size without null models", breaks = 1:10, labels = 1:10)


p6 <- ggarrange(ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, 
  font.label = list(size = arrange_label_size, family = ffont), legend = "none"),
  p3, labels = c("", "C"), font.label = list(size = arrange_label_size, family = ffont),
  nrow=2, common.legend = TRUE, legend = "bottom", legend.grob = ggpubr::get_legend(p1, "bottom"))

ggsave("figure-m6.pdf", plot = p6, path = "./plot/", height = two_row_height, width = two_col,
  dpi = 300)

# PIP correlation s17

ps171 <- ggplot(focus_analysis, aes(x = PIP.ME, y = PIP.meta)) +
  geom_point(size = 0.1) +
  theme_ma() +
  stat_cor(method = "pearson", cor.coef.name = c("r", "rho", "tau"), color = "red", fontface = "bold") +
  geom_smooth(method = "lm") +
  xlab("MA-FOCUS PIP") +
  ylab("Baseline PIP")

ps172 <- ggplot(focus_analysis, aes(x = PIP.ME, y = PIP.EA)) +
  geom_point(size = 0.1) +
  theme_ma() +
  stat_cor(method = "pearson", cor.coef.name = c("r", "rho", "tau"), color = "red") +
  geom_smooth(method = "lm") +
  xlab("MA-FOCUS PIP") +
  ylab("EA PIP")

ps173 <- ggplot(focus_analysis, aes(x = PIP.ME, y = PIP.AA)) +
  geom_point(size = 0.1) +
  theme_ma() +
  stat_cor(method = "pearson", cor.coef.name = c("r", "rho", "tau"), color = "red") +
  geom_smooth(method = "lm") +
  xlab("MA-FOCUS PIP") +
  ylab("AA PIP")


ps17 <- ggarrange(ps171, ps172, ps173, nrow = 1, labels = c("A", "B", "C"),
  font.label = list(size = arrange_label_size))
ggsave("figure-s17.png", plot = ps17, path = "./plot/", height = one_row_height, width = two_col)

# upset plots s18
library(ggvenn)
library(ggpubr)
aa <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.AA) %>%
  mutate(name = paste0(ID, PHEN))

ea <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.EA) %>%
  mutate(name = paste0(ID, PHEN))

me <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  filter(IN.CRED.SET.ME) %>%
  mutate(name = paste0(ID, PHEN))

a1 <- ggvenn(list(AA = aa$name, MA = me$name), text_size = 2)
a2 <- ggvenn(list(EA = ea$name, MA = me$name), text_size = 2)
a3 <- ggvenn(list(AA = aa$name, EA = ea$name), text_size = 2)
vg <- ggarrange(a1, a2, a3, nrow = 1)

library(UpSetR)

ps18dd <- focus_analysis %>%
  filter(!grepl("NULL", ID)) %>%
  mutate(ME = as.numeric(IN.CRED.SET.ME),
    meta = as.numeric(IN.CRED.SET.meta),
    EA = as.numeric(IN.CRED.SET.EA),
    AA = as.numeric(IN.CRED.SET.AA)) %>%
  select(ID, PHEN, ME, meta, EA, AA) %>%
  rename(`MA-FOCUS` = ME,
    Baseline = meta,
    `EA FOCUS` = EA,
    `AA FOCUS` = AA) %>%
  filter(complete.cases(.)) %>%
  as.data.frame()
us <- upset(ps18dd, sets = c("MA-FOCUS", "Baseline", "EA FOCUS", "AA FOCUS"),
  order.by = "freq", keep.order = T, empty.intersections = "on") 

usp <- cowplot::plot_grid(NULL, us$Main_bar, us$Sizes, us$Matrix, nrow=2, align="hv",
  rel_heights = c(1.5,1), rel_widths = c(1,3))
ps18 <- ggarrange(usp, vg, ncol = 1, labels = c("A", "B"), heights = c(3, 2))
ggsave(filename = "figure-s18.png", plot = ps18, path = "./plot/", width = two_col, height = two_row_height)


# TWAS manhattan plots m5
dd <- twas_all %>%
  filter(!is.na(TWAS.P) & TWAS.P != 0) %>%
  mutate(BP = P0+(P1-P0)/2,
    newP = 2*pnorm(abs(TWAS.Z), lower.tail = F),
    TWAS.logP = -1 * log(newP, 10)) %>%
  select(CHR, POP, ID, BP, TWAS.logP)

manP <- function(dd, pop, sig = 0.05/23000) {
  tmp <- dd %>%
    inner_join(dd %>%
        group_by(CHR) %>% 
        summarise(max_bp = max(BP)) %>% 
        mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
        select(CHR, bp_add), by = c("CHR")) %>% 
    mutate(bp_cum = BP + bp_add)
  
  axis_set <- tmp %>% 
    group_by(CHR) %>% 
    summarize(center = mean(bp_cum))
  
  ymaxlim <- max(tmp$TWAS.logP)*1.05
  
  if (pop == "EUR") {
    p <- ggplot(tmp, aes(x = bp_cum, y = TWAS.logP, 
      color = as_factor(CHR))) +
      geom_point(alpha = 0.75) +
      geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
      scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
      ylim(0, ymaxlim) +
      scale_color_manual(values = rep(c("#11A1E0", "#1156E0"), unique(length(axis_set$CHR)))) +
      labs(x = NULL, 
        y = bquote("EA Observed"~~-log[10](p))) + 
      theme_minimal() +
      theme(legend.position = "none",
        panel.border = element_blank(),
        # panel.border = element_rect(fill = NA),
        # axis.text.x = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(size = fontsize,  family = ffont))
  } else {
    p <- ggplot(tmp, aes(x = bp_cum, y = TWAS.logP, 
      color = as_factor(CHR))) +
      geom_point(alpha = 0.75) +
      geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
      scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
      scale_color_manual(values = rep(c("#11A1E0", "#1156E0"), unique(length(axis_set$CHR)))) +
      labs(x = NULL, 
        y = bquote("AA Observed"~~-log[10](p))) + 
      theme_minimal() +
      scale_y_reverse() +
      theme( legend.position = "none",
        panel.border = element_blank(),
        # panel.border = element_rect(fill = NA),
        # axis.text.x = element_text(size = 9),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_blank(),
        text=element_text(size = fontsize,  family = ffont))
  }
  return(p)
}

dd1 <- filter(dd, POP == "EA")
dd2 <- filter(dd, POP == "AA")

p1 <- manP(dd1, "EUR")
p2 <- manP(dd2, "AFR")

p <- ggarrange(p1, p2, ncol = 1, font.label = list(size = arrange_label_size, family = ffont))

dd3 <- read_tsv("./data/twas_gwas_corr.tsv") %>%
  mutate(TWAS.W  = (1/(TWAS.SE^2))/sum(1/(TWAS.SE^2)),
    GWAS.W  = (1/(GWAS.SE^2))/sum(1/(GWAS.SE^2)),
    TWAS.T = sum(TWAS.W*TWAS),
    TWAS.T.SE = sqrt( 1 / sum( 1 / TWAS.SE^2)),
    GWAS.T = sum(GWAS.W*GWAS),
    GWAS.T.SE = sqrt( 1 / sum( 1 / GWAS.SE^2)),
    TWASmin = TWAS - 1.96 * TWAS.T.SE,
    TWASmax = TWAS + 1.96 * TWAS.T.SE,
    GWASmin = GWAS - 1.96 *  GWAS.T.SE,
    GWASmax = GWAS + 1.96 * GWAS.T.SE,
    SIG = ifelse(TWASmin > GWAS | TWASmax < GWAS, PHEN, NA))

p52 <- ggplot(dd3, aes(x = GWAS, y = TWAS, label = SIG)) +
  geom_point(color = "black", size =0.5) +
  # geom_errorbarh(aes(xmin = GWASmin, xmax = GWASmax), width = 0.001) +
  geom_errorbar(aes(ymin = TWASmin, ymax = TWASmax), width = 0.001) +
  geom_abline(slope =1, intercept = 0, color = "red") +
  theme_corr() +
  xlab(expression(bold("GWAS Normalized ")~bolditalic(r))) +
  ylab(expression(bold("TWAS Normalized ")~bolditalic(r))) +
  theme(legend.position = "none") +
  geom_text_repel(size = 2,point.padding = NA,
    box.padding = 0.1) +
  scale_y_continuous(limits = c(0.025, 0.095), breaks = seq(0.025, 0.095, by = 0.01)) +
  scale_x_continuous(limits = c(0.022, 0.0701), breaks = seq(0.02, 0.0701, by = 0.01))

p5 <- ggarrange(p, p52, nrow = 2, labels = c("A", "B"), font.label = list(size = arrange_label_size, family = ffont))

ggsave("figure-m5.pdf", plot = p5, path = "./plot/", height = three_row_height*1.2, width = two_col,
  dpi = 300)


# Enrichment analysis  m7
theme_en <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey70"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    axis.title.y=element_text(face="bold"),
    axis.title.x=element_blank(),
    text=element_text(size = fontsize, face="bold", family = ffont))
}

load("./data/enrich.RData")
mapper <- read_csv("./data/DisGeNET_meta_categories.csv")

enridd <- enrichment.results.df %>%
  left_join(mapper %>%
      select(Term, `meta_category`),
    by = "Term") %>%
  filter(!is.na(`meta_category`)) %>%
  filter(`meta_category` %in% "Hematological measurement")

dd71 <- enridd %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  group_by(`meta_category`) %>%
  mutate(total = n()) %>%
  arrange(desc(total)) %>%
  group_by(`meta_category`, TEST) %>%
  summarize(n = n(),
    total = mean(total)) %>%
  group_by(`meta_category`) %>%
  mutate(ord = desc(n)) %>%
  mutate(TEST = factor(TEST, levels = c("EA", "AA", "ME", "meta"),
    labels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip")))

p71 <- ggplot(data=dd71, aes(fct_reorder(`meta_category`, total, .desc=TRUE), n, fill=TEST)) +
  geom_col(aes(group=ord), alpha = 0.8, position=position_dodge2(preserve = "single"), color = "black") +
  scale_fill_manual(values = COLS2, labels = LABELS2) +
  labs(y = "Number of\nenriched categories") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey70"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    axis.title.y=element_text(face="bold"),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    text=element_text(size = fontsize, face="bold", family = ffont))

# 2nd method
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

dd72 <- enrichment.results.df %>%
  filter(Term %in% lk) %>%
  mutate(Term = as.character(factor(Term, levels = lk, labels = phens))) %>%
  filter(Term == PHENO) %>%
  mutate(p = -1*log(Adjusted.P.value, 10)) %>%
  group_by(TEST) %>%
  arrange(desc(p)) %>%
  mutate(mm = row_number()) %>%
  group_by(Term) %>%
  mutate(kk = ifelse(TEST == "ME", mm, 0),
    kk = max(kk)) %>%
  select(Term, p, TEST, kk) %>%
  mutate(TEST = factor(TEST, levels = c("EA", "AA", "ME", "meta"),
    labels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip")))

p72 <- ggplot(dd72, aes(fct_reorder(Term, kk, .desc=FALSE), y = p, color = TEST)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = COLS2, labels = LABELS2) +
  theme_en() +
  ylab(expression(bold(-log)[10]~bold("(")~bold(p)[enrichment]~bold(")"))) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

p7 <- ggarrange(p71, p72, nrow = 1, labels = c("A", "B"),
  font.label = list(size = arrange_label_size, family = ffont), legend = "bottom", common.legend = TRUE)

ggsave("figure-m7.pdf", plot = p7, path = "./plot/", height = one_row_height, width = two_col,
  dpi = 300)

# s19
b1 <- focus_analysis %>%
  filter(IN.CRED.SET.ME == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)),
    logbf = log(bf))

b <- focus_analysis %>%
  filter(IN.CRED.SET.ME == 1) %>%
  select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
  mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)),
    logbf = log(bf),
    gr = "ME") %>%
  select(gr, logbf) %>%
  bind_rows(focus_analysis %>%
      filter(IN.CRED.SET.meta == 1) %>%
      select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
      mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)),
        logbf = log(bf),
        gr = "Baseline") %>%
      select(gr, logbf)) %>%
  bind_rows(focus_analysis %>%
      filter(IN.CRED.SET.EA == 1) %>%
      select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
      mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)),
        logbf = log(bf),
        gr = "EA") %>%
      select(gr, logbf)) %>%
  bind_rows(focus_analysis %>%
      filter(IN.CRED.SET.AA == 1) %>%
      select(BLOCK, PHEN, ID, PIP.ME, PIP.EA, PIP.AA) %>%
      mutate(bf = PIP.ME/((PIP.EA*(1-PIP.AA)) + PIP.AA*(1-PIP.EA)),
        logbf = log(bf),
        gr = "AA") %>%
      select(gr,logbf)) %>%
  mutate(gr = factor(gr, levels = c("ME", "Baseline", "EA", "AA"),
    labels = c("A: ME", "B: Baseline", "C: EA", "D: AA")))

ps19 <- ggplot(b, aes(x = logbf)) +
  facet_wrap(~gr, scales = "free") +
  geom_histogram() +
  theme_ma() +
  xlab("Log Bayes Factor") +
  ylab("Count") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(size = 0.5, color = "grey70"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    text=element_text(size = fontsize, face="bold", family = ffont))

ggsave("figure-s19.png", plot = ps19, path = "./plot/", height = one_row_height, width = one_col)


