# colour choices
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


scientific <- function(x){
  ifelse(x==0, "0",
    parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

combinePlots1 <- function(..., position = "right") {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  t1 <- textGrob("Shared eQTLs", gp = gpar(fontface = 4, fontsize = 12))
  t2 <- textGrob("Independent eQTLs", gp = gpar(fontface = 4, fontsize = 12))
  lay <- matrix(c(1, rep(3, 12), rep(4, 12), rep(5, 12), 2, rep(3, 12), rep(4, 12), rep(5, 12)),
    ncol = 2, nrow = 37)
  ga <- grid.arrange(t1, t2, gl[[1]], gl[[2]], gl[[3]], layout_matrix = lay)
  
  combined <- arrangeGrob(ga, legend, ncol = 2,
    widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}


combinePlots2 <- function(..., position = "right") {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  # t1 <- textGrob("eQTL shared across populations", gp = gpar(fontface = 4, fontsize = 12))
  # t2 <- textGrob("eQTLs independent across populations", gp = gpar(fontface = 4, fontsize = 12))
  lay <- matrix(c(rep(1, 12), rep(2, 12), rep(3, 12), rep(4, 12), rep(5, 12), rep(6, 12)),
    ncol = 2, nrow = 36)
  ga <- grid.arrange(gl[[1]], gl[[2]], gl[[3]], gl[[4]], gl[[5]], gl[[6]], layout_matrix = lay)
  
  combined <- arrangeGrob(ga, legend, ncol = 2,
    widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}
combinePlots21 <- function(..., position = "right") {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  # t1 <- textGrob("eQTL shared across populations", gp = gpar(fontface = 4, fontsize = 12))
  # t2 <- textGrob("eQTLs independent across populations", gp = gpar(fontface = 4, fontsize = 12))
  lay <- matrix(c(rep(1, 12), rep(2, 12), rep(3, 12), rep(1, 12), rep(2, 12), rep(3, 12),
    rep(4, 12), rep(5, 12), rep(6, 12)),
    ncol = 3, nrow = 36)
  ga <- grid.arrange(gl[[1]], gl[[2]], gl[[3]], gl[[4]], gl[[5]], gl[[6]], layout_matrix = lay)
  
  combined <- arrangeGrob(ga, legend, ncol = 2,
    widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}



# combinePlots2 <- function(..., position = "right") {
#   plots <- list(...)
#   g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
#   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#   lheight <- sum(legend$height)
#   lwidth <- sum(legend$width)
#   gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
#   ga <- grid.arrange(gl[[1]], gl[[2]], nrow = 2)
#   combined <- arrangeGrob(ga, legend, ncol = 2,
#     widths = unit.c(unit(1, "npc") - lwidth, lwidth))
#   grid.newpage()
#   grid.draw(combined)
#   invisible(combined)
# }

combinePlots3 <- function(..., position = "right") {
  plots <- list(...)
  g <- ggplotGrob(plots[[2]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  ga <- grid.arrange(gl[[1]], gl[[2]], gl[[3]], nrow = 3)
  combined <- arrangeGrob(ga, legend, ncol = 2,
    widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}

combinePlots4 <- function(..., position = "right") {
  plots <- list(...)
  g <- ggplotGrob(plots[[2]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  ga <- grid.arrange(gl[[1]], gl[[2]], gl[[3]], nrow = 3)
  combined <- arrangeGrob(ga, legend, ncol = 2,
    widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}

combinePlots5 <- function(..., position = "right", low = TRUE) {
  plots <- list(...)
  g <- ggplotGrob(plots[[2]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  
  if (low) {
    t <- textGrob("GWAS Sample Size: 400,000 for EUR and 8,000 for AFR",
      gp = gpar(fontface = 4, fontsize = 12))
  } else {
    t <- textGrob("GWAS Sample Size: 500,000 for EUR and 15,000 for AFR",
      gp = gpar(fontface = 4, fontsize = 12))
  }
  lay <- matrix(c(1, rep(2, 12), rep(3, 12), rep(4, 12)), ncol = 1, nrow = 37)
  ga <- grid.arrange(t, gl[[1]], gl[[2]], gl[[3]], layout_matrix = lay)
  
  combined <- arrangeGrob(ga, legend, ncol = 2,
    widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}

combinePlots6 <- function(..., position = "right") {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  t1 <- textGrob("eQTL shared across populations", gp = gpar(fontface = 4, fontsize = 12))
  t2 <- textGrob("eQTLs independent across populations", gp = gpar(fontface = 4, fontsize = 12))
  lay <- matrix(c(1, rep(3, 12), rep(4, 12), rep(5, 12), 2, rep(3, 12), rep(4, 12), rep(5, 12)), ncol = 2, nrow = 37)
  ga <- grid.arrange(t1, t2, gl[[1]], gl[[2]], gl[[3]], layout_matrix = lay)
  
  combined <- arrangeGrob(ga, legend, ncol = 2,
    widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}

combinePlots7 <- function(..., position = "right") {
  plots <- list(...)
  g <- ggplotGrob(plots[[2]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  
  t <- textGrob("GWAS: 511,471 (EUR) and 13,298 (AFR)\neQTL: 373 (EUR) and 441 (AFR)",
    gp = gpar(fontface = 4, fontsize = 12))
  
  lay <- matrix(c(rep(1, 3), rep(2, 12), rep(3, 12), rep(4, 12)), ncol = 1, nrow = 39)
  ga <- grid.arrange(t, gl[[1]], gl[[2]], gl[[3]], layout_matrix = lay)
  
  combined <- arrangeGrob(ga, legend, ncol = 2,
    widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}
