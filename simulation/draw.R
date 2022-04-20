rm(list = ls())
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(ggnewscale)
library(ggthemes)
library(scales)
library(gtable)
library(grid)
library(svglite)
library(ggtext)
source("./figure_style.R")
PIP.threshold <- 0.9
scaleFUN <- function(x) sprintf("%.2f", x)

new_labels1 <- c("shared" = "A", "indep" = "B")

new_labels2 <- c("shared" = "C", "indep" = "D")

new_labels3 <- c("shared" = "E", "indep" = "F")

H2GE <- c(0, 0.00001713863, 0.0001139271, 0.0007573169, 0.005034176)

###figure params
mean_size <- 2
mean_shape <- 22
median_size <- 2
median_shape <- 22
x_axis_size <- 8
y_axis_size <- 8

point_size <- 3
point_shape <- 21
errbar_width <- 0.5

ylabel_PIP <- "Causal gene PIP"
ylabel_num <- "Credible set size"
ylabel_freq <- "Freq. of causal\nin credible set"

arrange_label_size <- 8
one_row_height <- 2.5
two_row_height <- 4.5
three_row_height <- 6.5

one_col <- 3.35
onehalf_col <- 4.49
two_col <- 6.85

ffont <- "sans"
fontsize <- 8
theme_sim <- function() {
    theme(panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.y = element_line(size = 0.5, color = "grey70"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA),
        legend.title = element_blank(),
        axis.title=element_text(face="bold"),
        text=element_text(size = fontsize,  family = ffont))
}

addMeanMedian <- function(p) {
    p <- p +
        new_scale_fill() +
        stat_summary(fun = mean, geom = "point",  aes(fill = "Mean"), size = mean_size, shape = mean_shape,
            position = position_dodge(width = 0.3)) +
        stat_summary(fun = median, geom = "point", aes(fill = "Median"), shape = median_shape, size = median_size,
            position = position_dodge(width = 0.3)) +
        scale_fill_manual(values = c("Mean"="white", "Median"="black"))
    return(p)
}


load("./data/eur_afr.RDat")

dd1 <- dd %>%
    mutate(eqtl.model = factor(eqtl.model, levels = c("shared", "indep"))) %>%
    filter(nge1 == 200 & h2ge == 0.0007573169 & h2g == 0.05)

dd2 <- dd %>%
    mutate(eqtl.model = factor(eqtl.model, levels = c("shared", "indep"))) %>%
    filter(n1 == 100000 & h2ge == 0.0007573169 & h2g == 0.05)

# Plot 2
ddPIPN2 <- dd1 %>%
    filter(true_model %in% 1 & eqtl.model == "indep") %>%
    pivot_longer(cols = c(ME.pip, meta.pip), "PIP") %>%
    select(sim, n1, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(n1, PIP, sep = "_"),
        levels = c("50000_ME.pip", "50000_meta.pip", "1e+05_ME.pip", "1e+05_meta.pip",
            "2e+05_ME.pip", "2e+05_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 0.5,
            ifelse(X %in% c(5, 6), X + 1, X)))

pipN2 <- ggplot(ddPIPN2, aes(x = X, y = value, fill = PIP, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.5) +
    theme_sim() +
    scale_x_continuous(name = "GWAS sample size", breaks= c(1.5, 4, 6.5),
        labels = c("50000", "100000", "200000")) +
    scale_y_continuous(name = ylabel_PIP)
pipN2 <- addMeanMedian(pipN2)

ddGSN2 <- dd1 %>%
    filter(eqtl.model == "indep") %>%
    select(n1, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(n1, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(n1, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(n1, name, sep = "_"),
        levels = c("50000_ME.pip", "50000_meta.pip", "1e+05_ME.pip", "1e+05_meta.pip",
            "2e+05_ME.pip", "2e+05_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(n1, locus, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 0.5,
        ifelse(X %in% c(5, 6), X + 1, X))) %>%
    as.data.frame()

colnames(ddGSN2) <- c("n1", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsN2 <- ggplot(ddGSN2, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_violin(scale = "width", alpha = 0.8) +
    theme_sim() +
    scale_x_continuous(name = "GWAS sample size", breaks= c(1.5, 4, 6.5),
        labels = c("50000", "100000", "200000")) +
    scale_y_continuous(name= ylabel_num)
gsN2 <- addMeanMedian(gsN2)

ddFTN2 <- dd1 %>%
    filter(eqtl.model == "indep") %>%
    select(n1, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(n1, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(n1, name, locus, true_model, eqtl.model) %>%
    group_by(n1, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(n1, name, sep = "_"),
        levels = c("50000_ME.pip", "50000_meta.pip", "1e+05_ME.pip", "1e+05_meta.pip",
            "2e+05_ME.pip", "2e+05_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(n1, name) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTN2) <- c("n1", "test", "eqtl.model", "value", "Label", "X", "se")

ftN2 <- ggplot(ddFTN2, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_point(size = point_size, shape = point_shape) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.5) +
    theme_sim() +
    scale_x_continuous(name = "GWAS sample size", breaks= c(1.5, 3.5, 5.5),
        labels = c("50000", "100000", "200000")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

p2 <- ggarrange(pipN2, gsN2, ftN2, nrow = 1, labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size, family = ffont),
    common.legend = TRUE, legend = "bottom")

# ggsave("figure-m2.pdf", plot = p2, path = "./plot/", height = one_row_height, width = two_col,
#     dpi = 300)


# Plot s2
ddPIPNGES2 <- dd2 %>%
    filter(true_model %in% 1 & eqtl.model == "indep") %>%
    pivot_longer(cols = c(ME.pip, meta.pip), "PIP") %>%
    select(sim, nge1, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(nge1, PIP, sep = "_"),
        levels = c("100_ME.pip", "100_meta.pip", "200_ME.pip", "200_meta.pip",
            "300_ME.pip", "300_meta.pip", "400_ME.pip", "400_meta.pip",
            "500_ME.pip", "500_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 0.5,
            ifelse(X %in% c(5, 6), X + 1,
                ifelse(X %in% c(7, 8), X + 1.5,
                    ifelse(X %in% c(9, 10), X + 2, X)))))

pipNGES2 <- ggplot(ddPIPNGES2, aes(x=X, y=value, fill=PIP, group=Label)) +
    scale_fill_manual(values=COLS4, labels=LABELS4) +
    geom_violin(scale="width", alpha=0.8, position=position_dodge(width=0.3)) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "eQTL sample size",
        breaks= c(1.5, 4, 6.5, 9, 11.5),
        labels = c("100", "200", "300", "400", "500")) +
    scale_y_continuous(name = ylabel_PIP)
pipNGES2 <- addMeanMedian(pipNGES2)

ddGSNGES2 <- dd2 %>%
    filter(eqtl.model == "indep") %>%
    select(nge1, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(nge1, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(nge1, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(nge1, name, sep = "_"),
        levels = c("100_ME.pip", "100_meta.pip", "200_ME.pip", "200_meta.pip",
            "300_ME.pip", "300_meta.pip", "400_ME.pip", "400_meta.pip",
            "500_ME.pip", "500_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(nge1, locus, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 0.5,
        ifelse(X %in% c(5, 6), X + 1,
            ifelse(X %in% c(7, 8), X + 1.5,
                ifelse(X %in% c(9, 10), X + 2, X))))) %>%
    as.data.frame()

colnames(ddGSNGES2) <- c("nge1", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsNGES2 <- ggplot(ddGSNGES2, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS4, labels=LABELS4) +
    geom_violin(scale="width", alpha=0.8) +
    theme_sim() +
    scale_x_continuous(name = "eQTL sample size",
        breaks= c(1.5, 4, 6.5, 9, 11.5),
        labels = c("100", "200", "300", "400", "500")) +
    scale_y_continuous(name=ylabel_num)
gsNGES2 <- addMeanMedian(gsNGES2)

ddFTNGES2 <- dd2 %>%
    filter(eqtl.model == "indep") %>%
    select(nge1, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(nge1, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(nge1, name, locus, true_model, eqtl.model) %>%
    group_by(nge1, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(nge1, name, sep = "_"),
        levels = c("100_ME.pip", "100_meta.pip", "200_ME.pip", "200_meta.pip",
            "300_ME.pip", "300_meta.pip", "400_ME.pip", "400_meta.pip",
            "500_ME.pip", "500_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(nge1, name) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTNGES2) <- c("nge1", "test", "eqtl.model", "value", "Label", "X", "se")

ftNGES2 <- ggplot(ddFTNGES2, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS4, labels=LABELS4) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "eQTL sample size",
        breaks= c(1.5, 3.5, 5.5, 7.5, 9.5),
        labels = c("100", "200", "300", "400", "500")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ps2 <- ggarrange(pipNGES2, gsNGES2, ftNGES2, nrow = 1, labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")
# ggsave("figure-s2.png", plot = ps2, path = "./plot/", height = one_row_height, width = two_col)

# Plot s3
ddPIPNS3 <- dd1 %>%
    filter(true_model %in% 1 & eqtl.model == "shared") %>%
    pivot_longer(cols = c(ME.pip, meta.pip), "PIP") %>%
    select(sim, n1, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(n1, PIP, sep = "_"),
        levels = c("50000_ME.pip", "50000_meta.pip", "1e+05_ME.pip", "1e+05_meta.pip",
            "2e+05_ME.pip", "2e+05_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 0.5,
            ifelse(X %in% c(5, 6), X + 1, X)))

pipNS3 <- ggplot(ddPIPNS3, aes(x = X, y = value, fill = PIP, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "GWAS sample size", breaks= c(1.5, 4, 6.5),
        labels = c("50000", "100000", "200000")) +
    scale_y_continuous(name = ylabel_PIP)
pipNS3 <- addMeanMedian(pipNS3)

ddGSNS3 <- dd1 %>%
    filter(eqtl.model == "shared") %>%
    select(n1, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(n1, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(n1, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(n1, name, sep = "_"),
        levels = c("50000_ME.pip", "50000_meta.pip", "1e+05_ME.pip", "1e+05_meta.pip",
            "2e+05_ME.pip", "2e+05_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(n1, locus, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 0.5,
        ifelse(X %in% c(5, 6), X + 1, X))) %>%
    as.data.frame()

colnames(ddGSNS3) <- c("n1", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsNS3 <- ggplot(ddGSNS3, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_violin(scale = "width", alpha = 0.8) +
    theme_sim() +
    new_scale_fill() +
    scale_x_continuous(name = "GWAS sample size", breaks= c(1.5, 4, 6.5),
        labels = c("50000", "100000", "200000")) +
    scale_y_continuous(name=ylabel_num)
gsNS3 <- addMeanMedian(gsNS3)

ddFTNS3 <- dd1 %>%
    filter(eqtl.model == "shared") %>%
    select(n1, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(n1, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(n1, name, locus, true_model, eqtl.model) %>%
    group_by(n1, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(n1, name, sep = "_"),
        levels = c("50000_ME.pip", "50000_meta.pip", "1e+05_ME.pip", "1e+05_meta.pip",
            "2e+05_ME.pip", "2e+05_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(n1, name) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTNS3) <- c("n1", "test", "eqtl.model", "value", "Label", "X", "se")

ftNS3 <- ggplot(ddFTNS3, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "GWAS sample size", breaks= c(1.5, 3.5, 5.5),
        labels = c("50000", "100000", "200000")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ddPIPNGES3 <- dd2 %>%
    filter(true_model %in% 1 & eqtl.model == "shared") %>%
    pivot_longer(cols = c(ME.pip, meta.pip), "PIP") %>%
    select(sim, nge1, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(nge1, PIP, sep = "_"),
        levels = c("100_ME.pip", "100_meta.pip", "200_ME.pip", "200_meta.pip",
            "300_ME.pip", "300_meta.pip", "400_ME.pip", "400_meta.pip",
            "500_ME.pip", "500_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 0.5,
            ifelse(X %in% c(5, 6), X + 1,
                ifelse(X %in% c(7, 8), X + 1.5,
                    ifelse(X %in% c(9, 10), X + 2, X)))))

pipNGES3 <- ggplot(ddPIPNGES3, aes(x=X, y=value, fill=PIP, group=Label)) +
    scale_fill_manual(values=COLS4, labels=LABELS4) +
    geom_violin(scale="width", alpha=0.8, position=position_dodge(width=0.3)) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "eQTL sample size",
        breaks= c(1.5, 4, 6.5, 9, 11.5),
        labels = c("100", "200", "300", "400", "500")) +
    scale_y_continuous(name = ylabel_PIP)
pipNGES3 <- addMeanMedian(pipNGES3)

ddGSNGES3 <- dd2 %>%
    filter(eqtl.model == "shared") %>%
    select(nge1, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(nge1, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(nge1, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(nge1, name, sep = "_"),
        levels = c("100_ME.pip", "100_meta.pip", "200_ME.pip", "200_meta.pip",
            "300_ME.pip", "300_meta.pip", "400_ME.pip", "400_meta.pip",
            "500_ME.pip", "500_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(nge1, locus, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 0.5,
        ifelse(X %in% c(5, 6), X + 1,
            ifelse(X %in% c(7, 8), X + 1.5,
                ifelse(X %in% c(9, 10), X + 2, X))))) %>%
    as.data.frame()

colnames(ddGSNGES3) <- c("nge1", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsNGES3 <- ggplot(ddGSNGES3, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS4, labels=LABELS4) +
    geom_violin(scale="width", alpha=0.8) +
    theme_sim() +
    scale_x_continuous(name = "eQTL sample size",
        breaks= c(1.5, 4, 6.5, 9, 11.5),
        labels = c("100", "200", "300", "400", "500")) +
    scale_y_continuous(name=ylabel_num)
gsNGES3 <- addMeanMedian(gsNGES3)

ddFTNGES3 <- dd2 %>%
    filter(eqtl.model == "shared") %>%
    select(nge1, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(nge1, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(nge1, name, locus, true_model, eqtl.model) %>%
    group_by(nge1, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(nge1, name, sep = "_"),
        levels = c("100_ME.pip", "100_meta.pip", "200_ME.pip", "200_meta.pip",
            "300_ME.pip", "300_meta.pip", "400_ME.pip", "400_meta.pip",
            "500_ME.pip", "500_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(nge1, name) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTNGES3) <- c("nge1", "test", "eqtl.model", "value", "Label", "X", "se")

ftNGES3 <- ggplot(ddFTNGES3, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS4, labels=LABELS4) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "eQTL sample size",
        breaks= c(1.5, 3.5, 5.5, 7.5, 9.5),
        labels = c("100", "200", "300", "400", "500")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ps3 <- ggarrange(pipNS3, gsNS3, ftNS3, pipNGES3, gsNGES3, ftNGES3,
    labels = c("A", "B", "C", "D", "E","F"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")
# ggsave("figure-s3.png", plot = ps3, path = "./plot/", height = two_row_height, width = two_col)

dd3 <- dd %>%
    filter(n1 == 100000 & nge1 == 200 & eqtl.model == "indep" & h2g == 0.05)

ddPIPH2GE <- dd3 %>%
    filter(true_model %in% 1) %>%
    filter(!h2ge %in% c(4.418771e-05, 2.937327e-04, 1.952554e-03)) %>%
    pivot_longer(cols = c(ME.pip, meta.pip), "PIP") %>%
    select(sim, h2ge, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(h2ge, PIP, sep = "_"),
        levels = c("0_ME.pip", "0_meta.pip", "1.713863e-05_ME.pip",
            "1.713863e-05_meta.pip", "0.0001139271_ME.pip", "0.0001139271_meta.pip",
            "0.0007573169_ME.pip", "0.0007573169_meta.pip", "0.005034176_ME.pip",
            "0.005034176_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 1.5,
            ifelse(X %in% c(5, 6), X + 3,
                ifelse(X %in% c(7, 8), X + 4.5,
                    ifelse(X %in% c(9, 10), X + 6, X)))))

pipH2GE <- ggplot(ddPIPH2GE, aes(x = X, y = value, fill = PIP, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept = PIP.threshold, linetype="dashed", size = 0.25) +
    geom_vline(xintercept = 3.25, linetype = "dotted", size = 0.25) +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait expl. by causal gene expr. (")
        ~bolditalic(h)[GE]^{2}~bold(")")),
        breaks= c(1.5, 5, 8.5, 12, 15.5),
        labels = scientific(H2GE)) +
    scale_y_continuous(name = ylabel_PIP)
pipH2GE <- addMeanMedian(pipH2GE)

ddGSH2GE <- dd3 %>%
    filter(!h2ge %in% c(4.418771e-05, 2.937327e-04, 1.952554e-03)) %>%
    select(h2ge, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(h2ge, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(h2ge, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(h2ge, name, sep = "_"),
        levels = c("0_ME.pip", "0_meta.pip", "1.713863e-05_ME.pip",
            "1.713863e-05_meta.pip", "0.0001139271_ME.pip", "0.0001139271_meta.pip",
            "0.0007573169_ME.pip", "0.0007573169_meta.pip", "0.005034176_ME.pip",
            "0.005034176_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(h2ge, locus, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 1.5,
        ifelse(X %in% c(5, 6), X + 3,
            ifelse(X %in% c(7, 8), X + 4.5,
                ifelse(X %in% c(9, 10), X + 6, X))))) %>%
    as.data.frame()

colnames(ddGSH2GE) <- c("h2ge", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsH2GE <- ggplot(ddGSH2GE, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS4, labels=LABELS4) +
    geom_violin(scale="width", alpha=0.8) +
    geom_vline(xintercept = 3.25, linetype="dotted", size = 0.25) +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait expl. by causal gene expr. (")
        ~bolditalic(h)[GE]^{2}~bold(")")),
        breaks = c(1.5, 5, 8.5, 12, 15.5),
        labels = scientific(H2GE)) +
    scale_y_continuous(name=ylabel_num)
gsH2GE <- addMeanMedian(gsH2GE)

ddFTH2GE <- dd3 %>%
    filter(!h2ge %in% c(4.418771e-05, 2.937327e-04, 1.952554e-03)) %>%
    select(h2ge, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(h2ge, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(h2ge, name, locus, true_model, eqtl.model) %>%
    group_by(h2ge, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(h2ge, name, sep = "_"),
        levels = c("0_ME.pip", "0_meta.pip", "1.713863e-05_ME.pip",
            "1.713863e-05_meta.pip", "0.0001139271_ME.pip", "0.0001139271_meta.pip",
            "0.0007573169_ME.pip", "0.0007573169_meta.pip", "0.005034176_ME.pip",
            "0.005034176_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(h2ge, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 1.5,
        ifelse(X %in% c(5, 6), X + 3,
            ifelse(X %in% c(7, 8), X + 4.5,
                ifelse(X %in% c(9, 10), X + 6, X))))) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTH2GE) <- c("h2g", "test", "eqtl.model", "value", "Label", "X", "se")

ftH2GE <- ggplot(ddFTH2GE, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_vline(xintercept = 3.25, linetype="dotted", size = 0.25) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait expl. by causal gene expr. (")
        ~bolditalic(h)[GE]^{2}~bold(")")),
        breaks = c(1.5, 5, 8.5, 12, 15.5),
        labels = scientific(H2GE)) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ps5 <- ggarrange(pipH2GE, gsH2GE, ftH2GE, nrow = 3, labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")

# ggsave("figure-s5.png", plot = ps5, path = "./plot/", height = three_row_height, width = onehalf_col)

# plot s4
dd4 <- dd %>%
    filter(n1 == 100000 & h2ge == 0.0007573169 & nge1 == 200 & eqtl.model == "indep")

ddPIPH2G <- dd4 %>%
    filter(true_model %in% 1) %>%
    pivot_longer(cols = c(ME.pip, meta.pip), "PIP") %>%
    select(sim, h2g, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(h2g, PIP, sep = "_"),
        levels = c("0.01_ME.pip", "0.01_meta.pip", "0.05_ME.pip", "0.05_meta.pip",
            "0.1_ME.pip", "0.1_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 0.5,
            ifelse(X %in% c(5, 6), X + 1, X)))

pipH2G <- ggplot(ddPIPH2G, aes(x = X, y = value, fill = PIP, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = expression(bolditalic(cis)~"-"~bolditalic(h)[g]^{2}),
        breaks= c(1.5, 4, 6.5),
        labels = c("0.01", "0.05", "0.1")) +
    scale_y_continuous(name = ylabel_PIP)
pipH2G <- addMeanMedian(pipH2G)

ddGSH2G <- dd4 %>%
    select(h2g, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(h2g, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(h2g, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(h2g, name, sep = "_"),
        levels = c("0.01_ME.pip", "0.01_meta.pip", "0.05_ME.pip", "0.05_meta.pip",
            "0.1_ME.pip", "0.1_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(h2g, locus, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 0.5,
        ifelse(X %in% c(5, 6), X + 1, X))) %>%
    as.data.frame()

colnames(ddGSH2G) <- c("h2g", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsH2G <- ggplot(ddGSH2G, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_violin(scale = "width", alpha = 0.8) +
    theme_sim() +
    scale_x_continuous(name = expression(bolditalic(cis)~"-"~bolditalic(h)[g]^{2}),
        breaks= c(1.5, 4, 6.5),
        labels = c("0.01", "0.05", "0.1")) +
    scale_y_continuous(name=ylabel_num)
gsH2G <- addMeanMedian(gsH2G)

ddFTH2G <- dd4 %>%
    select(h2g, locus, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `ME.pip`:`meta.pip`) %>%
    group_by(h2g, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(h2g, name, locus, true_model, eqtl.model) %>%
    group_by(h2g, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(h2g, name, sep = "_"),
        levels = c("0.01_ME.pip", "0.01_meta.pip", "0.05_ME.pip", "0.05_meta.pip",
            "0.1_ME.pip", "0.1_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(h2g, name) %>%
    mutate(name = factor(name, levels = c("ME.pip", "meta.pip"))) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTH2G) <- c("h2g", "test", "eqtl.model", "value", "Label", "X", "se")

ftH2G <- ggplot(ddFTH2G, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS4, labels = LABELS4) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = expression(bolditalic(cis)~"-"~bolditalic(h)[g]^{2}),
        breaks= c(1.5, 3.5, 5.5),
        labels = c("0.01", "0.05", "0.1")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ps4 <- ggarrange(pipH2G, gsH2G, ftH2G, nrow = 1,labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")
# ggsave("figure-s4.png", plot = ps4, path = "./plot/", height = one_row_height, width = two_col)

# Figure S6
dd1 <- dd %>%
    filter(h2ge == 0.0007573169 & nge1 == 200 & h2g == 0.05 & eqtl.model == "indep")

tmp <- dd1 %>%
    filter(true_model %in% 1) %>%
    filter(n1 == 50000) %>%
    mutate(ME = ME.pip) %>%
    select(locus, ME) %>%
    left_join(dd1 %>%
            filter(true_model %in% 1) %>%
            filter(n1 == 100000) %>%
            mutate(SE = pop1.pip) %>%
            select(locus, SE),
        by = c("locus")) %>%
    mutate(group = "100k") %>%
    pivot_longer(c(ME, SE)) %>%
    bind_rows(dd1 %>%
            filter(true_model %in% 1) %>%
            filter(n1 == 100000) %>%
            mutate(ME = ME.pip) %>%
            select(locus, ME) %>%
            left_join(dd1 %>%
                    filter(true_model %in% 1) %>%
                    filter(n1 == 200000) %>%
                    mutate(SE = pop1.pip) %>%
                    select(locus, SE),
                by = c("locus")) %>%
            mutate(group = "200k") %>%
            pivot_longer(c(ME, SE)))

ddPIPCP <- tmp %>%
    mutate(Label = factor(paste(group, name, sep = "_"),
        levels = c("100k_SE", "100k_ME", "200k_SE", "200k_ME")),
        name = factor(name, levels = c("SE", "ME"), labels = c("pop1.pip", "ME.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 0.5, X))

pipCP <- ggplot(ddPIPCP, aes(x = X, y = value, fill = name, group = Label)) +
    scale_fill_manual(values = COLS5, labels = LABELS5) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "Total GWAS sample size",
        breaks= c(1.5, 4),
        labels = c("100000", "200000")) +
    scale_y_continuous(name = ylabel_PIP)
pipCP <- addMeanMedian(pipCP)

tmp <- dd1 %>%
    filter(n1 == 50000) %>%
    mutate(value = ME.pip) %>%
    group_by(sim, locus, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(ME = GS.size[1]) %>%
    ungroup() %>%
    select(locus, ME) %>%
    left_join(dd1 %>%
            filter(n1 == 100000) %>%
            mutate(value = pop1.pip) %>%
            group_by(sim, locus, eqtl.model) %>%
            arrange(desc(value)) %>%
            mutate(value = value/sum(value),
                cumpip = cumsum(value),
                GS.size = which(cumpip >= PIP.threshold)[1],
                rownum = 1:n()) %>%
            filter(GS.size >= rownum) %>%
            summarize(SE = GS.size[1]) %>%
            ungroup() %>%
            select(locus, SE),
        by = c("locus")) %>%
    mutate(group = "100k") %>%
    pivot_longer(c(ME, SE)) %>%
    bind_rows(dd1 %>%
            filter(n1 == 200000) %>%
            mutate(value = ME.pip) %>%
            group_by(sim, locus, eqtl.model) %>%
            arrange(desc(value)) %>%
            mutate(value = value/sum(value),
                cumpip = cumsum(value),
                GS.size = which(cumpip >= PIP.threshold)[1],
                rownum = 1:n()) %>%
            filter(GS.size >= rownum) %>%
            summarize(ME = GS.size[1]) %>%
            ungroup() %>%
            select(locus, ME) %>%
            left_join(dd1 %>%
                    filter(n1 == 100000) %>%
                    mutate(value = pop1.pip) %>%
                    group_by(sim, locus, eqtl.model) %>%
                    arrange(desc(value)) %>%
                    mutate(value = value/sum(value),
                        cumpip = cumsum(value),
                        GS.size = which(cumpip >= PIP.threshold)[1],
                        rownum = 1:n()) %>%
                    filter(GS.size >= rownum) %>%
                    summarize(SE = GS.size[1]) %>%
                    ungroup() %>%
                    select(locus, SE),
                by = c("locus")) %>%
            mutate(group = "200k") %>%
            pivot_longer(c(ME, SE)))

ddGSCP <- tmp %>%
    mutate(Label = factor(paste(group, name, sep = "_"),
        levels = c("100k_SE", "100k_ME", "200k_SE", "200k_ME")),
        name = factor(name, levels = c("SE", "ME"), labels = c("pop1.pip", "ME.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 0.5, X))

gsCP <- ggplot(ddGSCP, aes(x = X, y = value, fill = name, group = Label)) +
    scale_fill_manual(values = COLS5, labels = LABELS5) +
    geom_violin(scale = "width", alpha = 0.8) +
    theme_sim() +
    scale_x_continuous(name= "Total GWAS sample size",
        breaks= c(1.5, 4),
        labels = c("100000", "200000")) +
    scale_y_continuous(name=ylabel_num)
gsCP <- addMeanMedian(gsCP)

tmp <- dd1 %>%
    filter(n1 == 50000) %>%
    mutate(value = ME.pip) %>%
    group_by(sim, locus, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model = sum(true_model)) %>%
    group_by(eqtl.model) %>%
    summarize(ME = sum(true_model) / n()) %>%
    ungroup() %>%
    select(ME) %>%
    bind_cols(dd1 %>%
            filter(n1 == 100000) %>%
            mutate(value = pop1.pip) %>%
            group_by(sim, locus, eqtl.model) %>%
            arrange(desc(value)) %>%
            mutate(value = value/sum(value),
                cumpip = cumsum(value),
                GS.size = which(cumpip >= PIP.threshold)[1],
                rownum = 1:n()) %>%
            filter(GS.size >= rownum) %>%
            summarize(true_model = sum(true_model)) %>%
            group_by(eqtl.model) %>%
            summarize(SE = sum(true_model) / n()) %>%
            ungroup() %>%
            select(SE)) %>%
    mutate(group = "100k") %>%
    pivot_longer(c(ME, SE)) %>%
    bind_rows(dd1 %>%
            filter(n1 == 200000) %>%
            mutate(value = ME.pip) %>%
            group_by(sim, locus, eqtl.model) %>%
            arrange(desc(value)) %>%
            mutate(value = value/sum(value),
                cumpip = cumsum(value),
                GS.size = which(cumpip >= PIP.threshold)[1],
                rownum = 1:n()) %>%
            filter(GS.size >= rownum) %>%
            summarize(true_model = sum(true_model)) %>%
            group_by(eqtl.model) %>%
            summarize(ME = sum(true_model) / n()) %>%
            ungroup() %>%
            select( ME) %>%
            bind_cols(dd1 %>%
                    filter(n1 == 100000) %>%
                    mutate(value = pop1.pip) %>%
                    group_by(sim, locus, eqtl.model) %>%
                    arrange(desc(value)) %>%
                    mutate(value = value/sum(value),
                        cumpip = cumsum(value),
                        GS.size = which(cumpip >= PIP.threshold)[1],
                        rownum = 1:n()) %>%
                    filter(GS.size >= rownum) %>%
                    summarize(true_model = sum(true_model)) %>%
                    group_by(eqtl.model) %>%
                    summarize(SE = sum(true_model) / n()) %>%
                    ungroup() %>%
                    select(SE)) %>%
            mutate(group = "200k") %>%
            pivot_longer(c(ME, SE)))

ddFTCP <- tmp %>%
    mutate(Label = factor(paste(group, name, sep = "_"),
        levels = c("100k_SE", "100k_ME", "200k_SE", "200k_ME")),
        name = factor(name, levels = c("SE", "ME"), labels = c("pop1.pip", "ME.pip")),
        X = as.numeric(Label)) %>%
    mutate(se = sqrt(value*(1-value)/100))

ftCP <- ggplot(ddFTCP, aes(x = X, y = value, fill = name, group = Label)) +
    scale_fill_manual(values = COLS5, labels = LABELS5) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size = 0.25) +
    theme_sim() +
    scale_x_continuous(name= "Total GWAS sample size",
        breaks= c(1.5, 3.5),
        labels = c("100000", "200000")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ps6 <- ggarrange(pipCP, gsCP, ftCP, nrow = 1,labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")

# ggsave("figure-s6.png", plot = ps6, path = "./plot/", height = one_row_height, width = two_col)

# plot s7
load("./data/eur_afr_3pop.RDat")

tmpdd4 <- dd %>%
    filter(h2ge == 0.0007573169 & nge1 == 200 & h2g == 0.05)

load("./data/eur_afr_eas.RDat")

dd4 <- dd %>%
    filter(h2ge == 0.0007573169 & h2g == 0.05 & nge1 == 200 & n1 %in% c(100000, 50000, 25000)) %>%
    mutate(pop = "3-Pop") %>%
    bind_rows(tmpdd4 %>%
            mutate(pop = "2-Pop")) %>%
    filter(eqtl.model %in% "indep") %>%
    mutate(totalGS = ifelse(n1 %in% c(100000, 150000), "300k",
        ifelse(n1 %in% c(75000, 50000), "150k",
            ifelse(n1 %in% c(25000, 37500), "75k", NA)))) 

ddPIP3POPN1 <- dd4 %>%
    filter(true_model %in% 1) %>%
    pivot_longer(cols = ME.pip, "PIP") %>%
    select(sim, totalGS, PIP, value, eqtl.model, pop) %>%
    mutate(Label = factor(paste(totalGS, pop, sep = "_"),
        levels = c("75k_2-Pop", "75k_3-Pop", "150k_2-Pop", "150k_3-Pop",
            "300k_2-Pop", "300k_3-Pop")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(3, 4), X + 0.5,
            ifelse(X %in% c(5, 6), X + 1, X)))

pip3POPN1 <- ggplot(ddPIP3POPN1, aes(x = X, y = value, fill = pop, group = Label)) +
    scale_fill_manual(values=COLS6, labels=LABELS6) +
    geom_violin(scale="width", alpha=0.8, position=position_dodge(width=0.3)) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "Total GWAS sample size",
        breaks= c(1.5, 4, 6.5),
        labels = c("75k", "150k", "300k")) +
    scale_y_continuous(name = ylabel_PIP)
pip3POPN1 <- addMeanMedian(pip3POPN1)

ddGS3POPN1 <- dd4 %>%
    select(totalGS, locus, `ME.pip`, true_model, eqtl.model, pop) %>%
    gather(key = "name", value = "value", `ME.pip`) %>%
    group_by(totalGS, locus, name, eqtl.model, pop) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(totalGS, name, locus, GS.size, true_model, eqtl.model, pop) %>%
    mutate(Label = factor(paste(totalGS, pop, sep = "_"),
        levels = c("75k_2-Pop", "75k_3-Pop", "150k_2-Pop", "150k_3-Pop",
            "300k_2-Pop", "300k_3-Pop")),
        X = as.numeric(Label)) %>%
    arrange(totalGS, locus, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 0.5,
        ifelse(X %in% c(5, 6), X + 1, X))) %>%
    as.data.frame()

colnames(ddGS3POPN1) <- c("totalGS", "test","locus", "value", "contains_true_model", "eqtl.model","pop", "Label", "X")

gs3POPN1 <- ggplot(ddGS3POPN1, aes(x=X, y=value, fill=pop, group=Label)) +
    scale_fill_manual(values=COLS6, labels=LABELS6) +
    geom_violin(scale="width", alpha=0.8) +
    theme_sim() +
    scale_x_continuous(name = "Total GWAS sample size",
        breaks= c(1.5, 4, 6.5),
        labels = c("75k", "150k", "300k")) +
    scale_y_continuous(name=ylabel_num)
pip3POPN1 <- addMeanMedian(pip3POPN1)

ddFT3POPN1 <- dd4 %>%
    select(totalGS, locus, `ME.pip`, true_model, eqtl.model, pop) %>%
    gather(key = "name", value = "value", `ME.pip`) %>%
    group_by(totalGS, locus, name, eqtl.model, pop) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(totalGS, name, locus, true_model, eqtl.model, pop) %>%
    group_by(totalGS, pop, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(totalGS, pop, sep = "_"),
        levels = c("75k_2-Pop", "75k_3-Pop", "150k_2-Pop", "150k_3-Pop",
            "300k_2-Pop", "300k_3-Pop")),
        X = as.numeric(Label)) %>%
    arrange(pop, name) %>%
    mutate(X = ifelse(X %in% c(3, 4), X + 0.5,
        ifelse(X %in% c(5, 6), X + 1,  X))) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFT3POPN1) <- c("totalGS", "pop", "test", "eqtl.model", "value", "Label", "X", "se")

ft3POPN1 <- ggplot(ddFT3POPN1, aes(x = X, y = value, fill = pop, group = Label)) +
    scale_fill_manual(values = COLS6, labels = LABELS6) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "Total GWAS sample size",
        breaks= c(1.5, 4, 6.5),
        labels = c("75k", "150k", "300k")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels = scaleFUN)

ps7 <- ggarrange(pip3POPN1, gs3POPN1, ft3POPN1, nrow = 1,labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")

# ggsave("figure-s7.png", plot = ps7, path = "./plot/", height = one_row_height, width = two_col)


#Plot S8
load("./data/eur_afr_real.RDat")

tmp <- dd %>%
    filter(!h2ge %in% c(4.418771e-05, 2.937327e-04, 1.952554e-03)) %>%
    filter(eqtl.model == "indep")

ddPIPREAL2 <- tmp %>%
    filter(true_model %in% 1) %>%
    pivot_longer(cols = c(pop1.pip, pop2.pip, ME.pip, meta.pip), "PIP") %>%
    select(sim, h2ge, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(h2ge, PIP, sep = "_"),
        levels = c("0_pop1.pip", "0_pop2.pip", "0_ME.pip", "0_meta.pip", "1.713863e-05_pop1.pip",
            "1.713863e-05_pop2.pip","1.713863e-05_ME.pip", "1.713863e-05_meta.pip",
            "0.0001139271_pop1.pip", "0.0001139271_pop2.pip", "0.0001139271_ME.pip",
            "0.0001139271_meta.pip", "0.0007573169_pop1.pip", "0.0007573169_pop2.pip",
            "0.0007573169_ME.pip", "0.0007573169_meta.pip",  "0.005034176_pop1.pip",
            "0.005034176_pop2.pip", "0.005034176_ME.pip", "0.005034176_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(5, 6, 7, 8), X + 1.5,
            ifelse(X %in% c(9, 10, 11, 12), X + 3,
                ifelse(X %in% c(13, 14, 15, 16), X + 4.5,
                    ifelse(X %in% c(17, 18, 19, 20), X + 6, X)))),
        PIP = factor(PIP, levels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip")))

pipREAL2 <- ggplot(ddPIPREAL2, aes(x = X, y = value, fill = PIP, group = Label)) +
    scale_fill_manual(values = COLS2, labels = LABELS2) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    geom_vline(xintercept = 5.25, linetype="dotted", size = 0.25) +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait expl. by causal gene expr. (")
        ~bolditalic(h)[GE]^{2}~bold(")")),
        breaks= c(2.5, 8, 13.5, 19, 24.5),
        labels = scientific(H2GE)) +
    scale_y_continuous(name = ylabel_PIP)
pipREAL2 <- addMeanMedian(pipREAL2)

ddGSREAL2 <- tmp %>%
    select(h2ge, locus, `pop1.pip`, `pop2.pip`, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `pop1.pip`:`meta.pip`) %>%
    group_by(h2ge, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(h2ge, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(h2ge, name, sep = "_"),
        levels = c("0_pop1.pip", "0_pop2.pip", "0_ME.pip", "0_meta.pip", "1.713863e-05_pop1.pip",
            "1.713863e-05_pop2.pip","1.713863e-05_ME.pip", "1.713863e-05_meta.pip",
            "0.0001139271_pop1.pip", "0.0001139271_pop2.pip", "0.0001139271_ME.pip",
            "0.0001139271_meta.pip", "0.0007573169_pop1.pip", "0.0007573169_pop2.pip",
            "0.0007573169_ME.pip", "0.0007573169_meta.pip",  "0.005034176_pop1.pip",
            "0.005034176_pop2.pip", "0.005034176_ME.pip", "0.005034176_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(h2ge, locus, name) %>%
    mutate(X = ifelse(X %in% c(5, 6, 7, 8), X + 1.5,
        ifelse(X %in% c(9, 10, 11, 12), X + 3,
            ifelse(X %in% c(13, 14, 15, 16), X + 4.5,
                ifelse(X %in% c(17, 18, 19, 20), X + 6, X)))),
        name = factor(name, levels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip"))) %>%
    as.data.frame()

colnames(ddGSREAL2) <- c("h2ge", "test", "locus", "value", "contains_true_model", "eqtl.model",
    "Label", "X")

gsREAL2 <- ggplot(ddGSREAL2, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS2, labels=LABELS2) +
    geom_violin(scale="width", alpha=0.8) +
    geom_vline(xintercept = 5.25, linetype="dotted", size = 0.25) +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait expl. by causal gene expr. (")
        ~bolditalic(h)[GE]^{2}~bold(")")),
        breaks= c(2.5, 8, 13.5, 19, 24.5),
        labels = scientific(H2GE)) +
    scale_y_continuous(name=ylabel_num)
gsREAL2 <- addMeanMedian(gsREAL2)

ddFTREAL2 <- tmp %>%
    select(h2ge, locus, `pop1.pip`, `pop2.pip`, `ME.pip`, `meta.pip`, true_model, eqtl.model, pop) %>%
    gather(key = "name", value = "value", `pop1.pip`:`meta.pip`) %>%
    group_by(h2ge, locus, name, eqtl.model, pop) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model = sum(true_model)) %>%
    select(h2ge, name, locus, true_model, eqtl.model, pop) %>%
    group_by(h2ge, pop, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(h2ge, name, sep = "_"),
        levels = c("0_pop1.pip", "0_pop2.pip", "0_ME.pip", "0_meta.pip", "1.713863e-05_pop1.pip",
            "1.713863e-05_pop2.pip","1.713863e-05_ME.pip", "1.713863e-05_meta.pip",
            "0.0001139271_pop1.pip", "0.0001139271_pop2.pip", "0.0001139271_ME.pip",
            "0.0001139271_meta.pip", "0.0007573169_pop1.pip", "0.0007573169_pop2.pip",
            "0.0007573169_ME.pip", "0.0007573169_meta.pip", "0.005034176_pop1.pip",
            "0.005034176_pop2.pip", "0.005034176_ME.pip", "0.005034176_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(pop, name) %>%
    mutate(X = ifelse(X %in% c(5, 6, 7, 8), X + 1.5,
        ifelse(X %in% c(9, 10, 11, 12), X + 3,
            ifelse(X %in% c(13, 14, 15, 16), X + 4.5,
                ifelse(X %in% c(17, 18, 19, 20), X + 6, X)))),
        name = factor(name, levels = c("pop1.pip", "pop2.pip", "ME.pip", "meta.pip"))) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTREAL2) <- c("h2ge", "name", "test", "eqtl.model", "value", "Label", "X", "se")

ftREAL2 <- ggplot(ddFTREAL2, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS2, labels = LABELS2) +
    geom_vline(xintercept = 5.25, linetype="dotted", size = 0.25) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait expl. by causal gene expr. (")
        ~bolditalic(h)[GE]^{2}~bold(")")),
        breaks= c(2.5, 8, 13.5, 19, 24.5),
        labels = scientific(H2GE)) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1.01), labels = scaleFUN)

ps8 <- ggarrange(pipREAL2, gsREAL2, ftREAL2, nrow = 3, labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")

# ggsave("figure-s8.png", plot = ps8, path = "./plot/", height = three_row_height, width = two_col)

# figure 3
load("./data/pop2h2ge.RDat")
dd5 <- dd %>%
    filter(h2ge2 %in% c(0.00001713863, 0.0001139271, 0.0007573169, 0.005034176) &
            eqtl.model %in% "indep")

H2GE <- c(0.00001713863, 0.0001139271, 0.0007573169, 0.005034176)

tmpddPIPH2GE2 <- dd5 %>%
    filter(true_model %in% 1) %>%
    pivot_longer(cols = c(pop2.pip, ME.pip, meta.pip), "PIP") %>%
    select(sim, h2ge2, PIP, value, eqtl.model) %>%
    bind_rows(dd5 %>%
            filter(true_model %in% 1) %>%
            pivot_longer(cols = c(pop1.pip), "PIP") %>%
            filter(h2ge2 == 0.0007573169) %>%
            select(sim, h2ge2, PIP, value, eqtl.model))
eurdd <- tmpddPIPH2GE2 %>%
    filter(PIP == "pop1.pip")
eurddMean <- mean(eurdd$value)
eurddMedian <- median(eurdd$value)

ddPIPH2GE2 <- tmpddPIPH2GE2 %>%
    filter(PIP != "pop1.pip") %>%
    mutate(Label = factor(paste(h2ge2, PIP, sep = "_"),
        levels = c("1.713863e-05_pop2.pip", "1.713863e-05_ME.pip",
            "1.713863e-05_meta.pip", "0.0001139271_pop2.pip", "0.0001139271_ME.pip", "0.0001139271_meta.pip",
            "0.0007573169_pop2.pip", "0.0007573169_ME.pip", "0.0007573169_meta.pip", "0.005034176_pop2.pip",
            "0.005034176_ME.pip", "0.005034176_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(4, 5, 6), X + 1.5,
            ifelse(X %in% c(7, 8, 9), X + 3,
                ifelse(X %in% c(10, 11, 12), X + 4.5, X))),
        PIP = factor(PIP, levels = c("pop2.pip", "ME.pip", "meta.pip")))

pipH2GE2 <- ggplot(ddPIPH2GE2, aes(x = X, y = value, fill = PIP, group = Label)) +
    scale_fill_manual(values = COLS3, labels = LABELS3) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept = PIP.threshold, linetype="dashed", size = 0.5) +
    geom_hline(yintercept = eurddMean, linetype="dotted", size = 1, color = "#F6795E") +
    geom_hline(yintercept = eurddMedian, linetype="dotted", size = 1, color = "#CE15EF") +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait explained by causal gene expression (")
        ~bolditalic(h)[GE]^{2}~bold(") for AFR")),
        breaks= c(2, 6.5, 11, 15.5),
        labels = scientific(H2GE)) +
    scale_y_continuous(name = ylabel_PIP)
pipH2GE2 <- addMeanMedian(pipH2GE2)

p3 <- pipH2GE2
# ggsave("figure-m3.pdf", plot = p3, path = "./plot/", height = one_row_height, width = two_col,
#     dpi = 300)


# figure s9
tmpddGSH2GE2 <- dd5 %>%
    select(h2ge2, locus,`pop1.pip`, `pop2.pip`, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `pop1.pip`:`meta.pip`) %>%
    group_by(h2ge2, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    filter(!(name == "pop1.pip" & h2ge2 != 0.0007573169)) %>%
    select(h2ge2, name, locus, GS.size, true_model, eqtl.model)

eurgsdd <- tmpddGSH2GE2 %>%
    filter(name == "pop1.pip")
eurgsddMean <- mean(eurgsdd$`GS.size`)
eurgsddMedian <- median(eurgsdd$`GS.size`)

ddGSH2GE2 <- tmpddGSH2GE2 %>%
    filter(name != "pop1.pip") %>%
    mutate(Label = factor(paste(h2ge2, name, sep = "_"),
        levels = c("1.713863e-05_pop2.pip", "1.713863e-05_ME.pip",
            "1.713863e-05_meta.pip", "0.0001139271_pop2.pip", "0.0001139271_ME.pip", "0.0001139271_meta.pip",
            "0.0007573169_pop2.pip", "0.0007573169_ME.pip", "0.0007573169_meta.pip", "0.005034176_pop2.pip",
            "0.005034176_ME.pip", "0.005034176_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(h2ge2, locus, name) %>%
    mutate(X = ifelse(X %in% c(4, 5, 6), X + 1.5,
        ifelse(X %in% c(7, 8, 9), X + 3,
            ifelse(X %in% c(10, 11, 12), X + 4.5, X))),
        name = factor(name, levels = c("pop2.pip", "ME.pip", "meta.pip"))) %>%
    as.data.frame()

colnames(ddGSH2GE2) <- c("h2ge2", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsH2GE2 <- ggplot(ddGSH2GE2, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS3, labels=LABELS3) +
    geom_violin(scale="width", alpha=0.8) +
    geom_hline(yintercept = eurgsddMean, linetype="dotted", size = 1, color = "#F6795E") +
    geom_hline(yintercept = eurgsddMedian, linetype="dotted", size = 1, color = "#CE15EF") +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait explained by causal gene expression (")
        ~bolditalic(h)[GE]^{2}~bold(") for AFR")),
        breaks= c(2, 6.5, 11, 15.5),
        labels = scientific(H2GE)) +
    scale_y_continuous(name=ylabel_num)
gsH2GE2 <- addMeanMedian(gsH2GE2)

tmpddFTH2GE2 <- dd5 %>%
    select(h2ge2, locus, `pop1.pip`, `pop2.pip`, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `pop1.pip`:`meta.pip`) %>%
    group_by(h2ge2, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    filter(!(name == "pop1.pip" & h2ge2 != 0.0007573169)) %>%
    select(h2ge2, name, locus, true_model, eqtl.model) %>%
    group_by(h2ge2, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n())

ddeurft <- tmpddFTH2GE2 %>%
    filter(name == "pop1.pip")
ddeurftMean <- ddeurft$freq

ddFTH2GE2 <- tmpddFTH2GE2 %>%
    filter(name != "pop1.pip") %>%
    mutate(Label = factor(paste(h2ge2, name, sep = "_"),
        levels = c("1.713863e-05_pop2.pip", "1.713863e-05_ME.pip",
            "1.713863e-05_meta.pip", "0.0001139271_pop2.pip", "0.0001139271_ME.pip", "0.0001139271_meta.pip",
            "0.0007573169_pop2.pip", "0.0007573169_ME.pip", "0.0007573169_meta.pip", "0.005034176_pop2.pip",
            "0.005034176_ME.pip", "0.005034176_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(h2ge2, name) %>%
    mutate(X = ifelse(X %in% c(4, 5, 6), X + 1.5,
        ifelse(X %in% c(7, 8, 9), X + 3,
            ifelse(X %in% c(10, 11, 12), X + 4.5, X))),
        name = factor(name, levels = c("pop2.pip", "ME.pip", "meta.pip"))) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTH2GE2) <- c("h2ge2", "test", "eqtl.model", "value", "Label", "X", "se")

ftH2GE2 <- ggplot(ddFTH2GE2, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS3, labels = LABELS3) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    geom_hline(yintercept = ddeurftMean, linetype="dotted", size = 1, color = "#F6795E") +
    theme_sim() +
    scale_x_continuous(name = expression(bold("Variation in trait explained by causal gene expression (")
        ~bolditalic(h)[GE]^{2}~bold(") for AFR")),
        breaks= c(2, 6.5, 11, 15.5),
        labels = scientific(H2GE))  +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ps9 <- ggarrange(gsH2GE2, ftH2GE2, nrow = 2, labels = c("A", "B"),
    common.legend = TRUE, legend = "bottom")

# ggsave("figure-s9.png", plot = ps9, path = "./plot/", height = two_row_height, width = two_col)

# plot s10
# proxy
load("./data/eur_afr.RDat")
tmpdd6 <- dd %>%
    filter(nge1 == 200 & h2ge == 0.0007573169 & h2g == 0.05 & n1 == 100000)
load("./data/proxy.RDat")

dd6 <- dd %>%
    bind_rows(tmpdd6) %>%
    filter(eqtl.model == "indep") %>%
    mutate(rho = ifelse(is.na(rho), 1, rho),
        proxy_tissue_pop = "pop2")

ddPIPProxy <- dd6 %>%
    filter(true_model %in% 1) %>%
    pivot_longer(cols = c(pop2.pip, ME.pip, meta.pip), "PIP") %>%
    select(sim, rho, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(rho, PIP, sep = "_"),
        levels = c("0_pop2.pip", "0_ME.pip", "0_meta.pip", "0.3_pop2.pip", "0.3_ME.pip",
            "0.3_meta.pip", "0.6_pop2.pip", "0.6_ME.pip", "0.6_meta.pip", "0.9_pop2.pip",
            "0.9_ME.pip", "0.9_meta.pip", "1_pop2.pip", "1_ME.pip", "1_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(4, 5, 6), X + 1.5,
            ifelse(X %in% c(7, 8, 9), X + 3,
                ifelse(X %in% c(10, 11, 12), X + 4.5,
                    ifelse(X %in% c(13, 14, 15), X + 6, X)))),
        PIP = factor(PIP, levels = c("pop2.pip", "ME.pip", "meta.pip")))

pipProxy <- ggplot(ddPIPProxy, aes(x = X, y = value, fill = PIP, group = Label)) +
    scale_fill_manual(values = COLS3, labels = LABELS3) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept = PIP.threshold, linetype="dashed", size = 0.25) +
    geom_vline(xintercept = 4.25, linetype="dotted", size = 0.25) +
    theme_sim() +
    scale_x_continuous(name = "Correlation between causal and proxy tissue used for AFR",
        breaks= c(2, 6.5, 11, 15.5, 20),
        labels = c("0", "0.3", "0.6", "0.9", "1 (causal tissue used)")) +
    scale_y_continuous(name = ylabel_PIP)
pipProxy <- addMeanMedian(pipProxy)

ddGSProxy <- dd6 %>%
    select(rho, locus, `pop2.pip`, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `pop2.pip`:`meta.pip`) %>%
    group_by(rho, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(rho, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(rho, name, sep = "_"),
        levels = c("0_pop2.pip", "0_ME.pip", "0_meta.pip", "0.3_pop2.pip", "0.3_ME.pip",
            "0.3_meta.pip", "0.6_pop2.pip", "0.6_ME.pip", "0.6_meta.pip", "0.9_pop2.pip",
            "0.9_ME.pip", "0.9_meta.pip", "1_pop2.pip", "1_ME.pip", "1_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(rho, locus, name) %>%
    mutate(X = ifelse(X %in% c(4, 5, 6), X + 1.5,
        ifelse(X %in% c(7, 8, 9), X + 3,
            ifelse(X %in% c(10, 11, 12), X + 4.5,
                ifelse(X %in% c(13, 14, 15), X + 6, X)))),
        name = factor(name, levels = c("pop2.pip", "ME.pip", "meta.pip"))) %>%
    as.data.frame()

colnames(ddGSProxy) <- c("rho", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsProxy <- ggplot(ddGSProxy, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS3, labels=LABELS3) +
    geom_violin(scale="width", alpha=0.8) +
    theme_sim() +
    scale_x_continuous(name = "Correlation between causal and proxy tissue used for AFR",
        breaks= c(2, 6.5, 11, 15.5, 20),
        labels = c("0", "0.3", "0.6", "0.9", "1 (causal tissue used)")) +
    scale_y_continuous(name=ylabel_num)
gsProxy <- addMeanMedian(gsProxy)

ddFTProxy <- dd6 %>%
    select(rho, locus, `pop2.pip`, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `pop2.pip`:`meta.pip`) %>%
    group_by(rho, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(rho, name, locus, true_model, eqtl.model) %>%
    group_by(rho, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(rho, name, sep = "_"),
        levels = c("0_pop2.pip", "0_ME.pip", "0_meta.pip", "0.3_pop2.pip", "0.3_ME.pip",
            "0.3_meta.pip", "0.6_pop2.pip", "0.6_ME.pip", "0.6_meta.pip", "0.9_pop2.pip",
            "0.9_ME.pip", "0.9_meta.pip", "1_pop2.pip", "1_ME.pip", "1_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(rho, name) %>%
    mutate(X = ifelse(X %in% c(4, 5, 6), X + 1.5,
        ifelse(X %in% c(7, 8, 9), X + 3,
            ifelse(X %in% c(10, 11, 12), X + 4.5,
                ifelse(X %in% c(13, 14, 15), X + 6, X)))),
        name = factor(name, levels = c("pop2.pip", "ME.pip", "meta.pip"))) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTProxy) <- c("rho", "test", "eqtl.model", "value", "Label", "X", "se")

ftProxy <- ggplot(ddFTProxy, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS3, labels = LABELS3) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = errbar_width) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "Correlation between causal and proxy tissue used for AFR",
        breaks= c(2, 6.5, 11, 15.5, 20),
        labels = c("0", "0.3", "0.6", "0.9", "1 (causal tissue used)")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ps10 <- ggarrange(pipProxy, gsProxy, ftProxy, nrow = 3,
    labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")

# ggsave("figure-s10.png", plot = ps10, path = "./plot/", height = three_row_height, width = two_col)

# s11
load("./data/eur_afr.RDat")
tmpdd7 <- dd %>%
    filter(nge1 == 200 & h2ge == 0.0007573169 & h2g == 0.05 & n1 == 100000)
load("./data/pop2weights.RDat")

dd7 <- dd %>%
    mutate(pop = "EUR") %>%
    bind_rows(tmpdd7 %>%
            mutate(pop = "AFR")) %>%
    filter(eqtl.model == "indep")

ddPIPPOP2W <- dd7 %>%
    filter(true_model %in% 1) %>%
    pivot_longer(cols = c(pop2.pip, ME.pip, meta.pip), "PIP") %>%
    select(sim, pop, PIP, value, eqtl.model) %>%
    mutate(Label = factor(paste(pop, PIP, sep = "_"),
        levels = c("EUR_pop2.pip", "EUR_ME.pip", "EUR_meta.pip", "AFR_pop2.pip", "AFR_ME.pip",
            "AFR_meta.pip")),
        X = as.numeric(Label),
        X = ifelse(X %in% c(4, 5, 6), X + 1, X),
        PIP = factor(PIP, levels = c("pop2.pip", "ME.pip", "meta.pip")))

pipPOP2W <- ggplot(ddPIPPOP2W, aes(x = X, y = value, fill = PIP, group = Label)) +
    scale_fill_manual(values = COLS3, labels = LABELS3) +
    geom_violin(scale = "width", alpha = 0.8, position=position_dodge(width = 0.3)) +
    geom_hline(yintercept = PIP.threshold, linetype="dashed", size = 0.25) +
    theme_sim() +
    scale_x_continuous(name = "eQTL wgt. used for AFR",
        breaks= c(2, 6),
        labels = c("EUR", "AFR")) +
    scale_y_continuous(name = ylabel_PIP)
pipPOP2W <- addMeanMedian(pipPOP2W)

ddGSPOP2W<- dd7 %>%
    select(pop, locus, `pop2.pip`, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `pop2.pip`:`meta.pip`) %>%
    group_by(pop, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model),
        GS.size = GS.size[1]) %>%
    select(pop, name, locus, GS.size, true_model, eqtl.model) %>%
    mutate(Label = factor(paste(pop, name, sep = "_"),
        levels = c("EUR_pop2.pip", "EUR_ME.pip", "EUR_meta.pip", "AFR_pop2.pip", "AFR_ME.pip",
            "AFR_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(pop, locus, name) %>%
    mutate(X = ifelse(X %in% c(4, 5, 6), X + 1, X),
        name = factor(name, levels = c("pop2.pip", "ME.pip", "meta.pip"))) %>%
    as.data.frame()

colnames(ddGSPOP2W) <- c("pop", "test","locus", "value", "contains_true_model", "eqtl.model", "Label", "X")

gsPOP2W <- ggplot(ddGSPOP2W, aes(x=X, y=value, fill=test, group=Label)) +
    scale_fill_manual(values=COLS3, labels=LABELS3) +
    geom_violin(scale="width", alpha=0.8) +
    theme_sim() +
    scale_x_continuous(name = "eQTL wgt. used for AFR",
        breaks= c(2, 6),
        labels = c("EUR", "AFR")) +
    scale_y_continuous(name=ylabel_num)
gsPOP2W <- addMeanMedian(gsPOP2W)

ddFTPOP2W <- dd7 %>%
    select(pop, locus, `pop2.pip`, `ME.pip`, `meta.pip`, true_model, eqtl.model) %>%
    gather(key = "name", value = "value", `pop2.pip`:`meta.pip`) %>%
    group_by(pop, locus, name, eqtl.model) %>%
    arrange(desc(value)) %>%
    mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
    filter(GS.size >= rownum) %>%
    summarize(true_model=sum(true_model)) %>%
    select(pop, name, locus, true_model, eqtl.model) %>%
    group_by(pop, name, eqtl.model) %>%
    summarize(freq = sum(true_model) / n()) %>%
    mutate(Label = factor(paste(pop, name, sep = "_"),
        levels = c("EUR_pop2.pip", "EUR_ME.pip", "EUR_meta.pip", "AFR_pop2.pip", "AFR_ME.pip",
            "AFR_meta.pip")),
        X = as.numeric(Label)) %>%
    arrange(pop, name) %>%
    mutate(X = ifelse(X %in% c(4, 5, 6), X + 1, X),
        name = factor(name, levels = c("pop2.pip", "ME.pip", "meta.pip"))) %>%
    mutate(se = sqrt(freq*(1-freq)/100)) %>%
    as.data.frame()

colnames(ddFTPOP2W) <- c("pop", "test", "eqtl.model", "value", "Label", "X", "se")

ftPOP2W <- ggplot(ddFTPOP2W, aes(x = X, y = value, fill = test, group = Label)) +
    scale_fill_manual(values = COLS3, labels = LABELS3) +
    geom_point(size = point_size, shape = point_shape) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.2) +
    geom_hline(yintercept=PIP.threshold, linetype="dashed", size=0.25) +
    theme_sim() +
    scale_x_continuous(name = "eQTL wgt. used for AFR",
        breaks= c(2, 6),
        labels = c("EUR", "AFR")) +
    scale_y_continuous(name = ylabel_freq,
        limits = c(0.6, 1), labels=scaleFUN)

ps11 <- ggarrange(pipPOP2W, gsPOP2W, ftPOP2W, nrow = 1, labels = c("A", "B", "C"),
    font.label = list(size = arrange_label_size),
    common.legend = TRUE, legend = "bottom")

# ggsave("figure-s11.png", plot = ps11, path = "./plot/", height = one_row_height, width = two_col)
