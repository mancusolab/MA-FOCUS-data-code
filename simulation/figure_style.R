# color choices
COLS <- c("ME.pip"="#CC6677", "meta.pip"="#009988",
  "pop1.pip"="#DDCC77", "pop2.pip"="#0077BB", "pop3.pip"="#b366cc",
  "2-Pop" = "#24492e", "3-Pop" = "#89689d")

LABELS <- c("pop1.pip"="EUR FOCUS", "pop2.pip"="AFR FOCUS",
  "ME.pip"="MA-FOCUS", "meta.pip"="Baseline",
  "2-Pop" = "2-Pop MA-FOCUS", "3-Pop" = "3-Pop MA-FOCUS")

COLS2 <- c("pop1.pip"="#DDCC77", "pop2.pip"="#0077BB", "ME.pip"="#CC6677", "meta.pip"="#009988")

LABELS2 <- c("pop1.pip"="EUR FOCUS", "pop2.pip"="AFR FOCUS",
  "ME.pip"="MA-FOCUS", "meta.pip"="Baseline")

COLS3 <- c( "pop2.pip"="#0077BB", "ME.pip"="#CC6677", "meta.pip"="#009988")

LABELS3 <- c("pop2.pip"="AFR FOCUS",
  "ME.pip"="MA-FOCUS", "meta.pip"="Baseline")

COLS4 <- c("ME.pip"="#CC6677", "meta.pip"="#009988")

LABELS4 <- c("ME.pip"="MA-FOCUS", "meta.pip"="Baseline")

COLS5 <- c("ME.pip"="#CC6677", "pop1.pip"="#DDCC77")

LABELS5 <- c("ME.pip"="MA-FOCUS", "pop1.pip"="EUR FOCUS")

COLS6 <- c("2-Pop" = "#24492e", "3-Pop" = "#89689d")

LABELS6 <- c("2-Pop" = "Two-ancestry MA-FOCUS", "3-Pop" = "Three-ancestry MA-FOCUS")


scientific <- function(x){
  ifelse(x==0, "0",
    parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}
