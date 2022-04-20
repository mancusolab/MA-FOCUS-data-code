# colour choices
COLS <- c("ME.pip"="#CC6677", "meta.pip"="#009988",
  "pop1.pip"="#DDCC77", "pop2.pip"="#0077BB", "pop3.pip"="#b366cc",
  "2-Pop" = "#24492e", "3-Pop" = "#89689d", "total" = "#F1AC0C")

LABELS <- c("pop1.pip"="EUR FOCUS", "pop2.pip"="AFR FOCUS",
  "ME.pip"="MA-FOCUS", "meta.pip"="Baseline",
  "2-Pop" = "2-Pop MA-FOCUS", "3-Pop" = "3-Pop MA-FOCUS",
  "total" = "Genes in Region")

COLS2 <- c("pop1.pip"="#DDCC77", "pop2.pip"="#0077BB", "ME.pip"="#CC6677", "meta.pip"="#009988")

LABELS2 <- c("pop1.pip"="EA FOCUS", "pop2.pip"="AA FOCUS",
  "ME.pip"="MA-FOCUS", "meta.pip"="Baseline")

COLS3 <- c("EA"="#DDCC77", "AA"="#0077BB", "ME"="#CC6677", "meta"="#009988")

LABELS3 <- c("EA"="EA FOCUS", "AA"="AA FOCUS",
  "ME"="MA-FOCUS", "meta"="Baseline")


scientific <- function(x){
  ifelse(x==0, "0",
    parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}