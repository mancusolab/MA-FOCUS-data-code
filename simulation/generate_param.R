library(tidyverse)

# case 1
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- c(50000, 100000, 200000)
NGE1 <- NGE2 <- 200

tt1 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2) %>%
  filter(N1 == N2) 

# case 2
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- 100000
NGE1 <- NGE2 <- c(100, 200, 300, 400, 500)

tt2 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2) %>%
  filter(NGE1 == NGE2) 

# case 3
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- c(0, 0.00001713863, 0.00004418771, 0.0001139271,
  0.0002937327, 0.0007573169, 0.001952554, 0.005034176)
N1 <- N2 <- 100000
NGE1 <- NGE2 <- 200

tt3 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2)

# case 4
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- c(0.01, 0.05, 0.1)
H2GE <- 0.0007573169
N1 <- N2 <- 100000
NGE1 <- NGE2 <- 200

tt4 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2)

tt <- tt1 %>%
  bind_rows(tt2, tt3, tt4) %>%
  distinct() %>%
  rownames_to_column("ROW")

dd <- crossing(ROW = tt$ROW, LOCUS = 1:100) %>%
  left_join(tt, by = "ROW")

write_tsv(dd, path="./codes/param/param_eur_afr.tsv", col_names = FALSE, quote_escape = FALSE)

# compare to 3 pop case1
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- 150000
NGE1 <- NGE2 <- c(100, 200, 300, 400, 500)

tt1 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2) %>%
  filter(NGE1 == NGE2) 

# compare to 3 pop case2
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- c(0, 0.00001713863, 0.00004418771, 0.0001139271,
  0.0002937327, 0.0007573169, 0.001952554, 0.005034176)
N1 <- N2 <- 150000
NGE1 <- NGE2 <- 200

tt2 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2)

# compare to 3 pop case3
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- c(0.01, 0.05, 0.1)
H2GE <- 0.0007573169
N1 <- N2 <- 150000
NGE1 <- NGE2 <- 200

tt3 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2)

# compare to 3 pop case4
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- c(150000, 75000, 37500)
NGE1 <- NGE2 <- 200

tt4 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2) %>%
  filter(N1 == N2)

tt <- tt1 %>%
  bind_rows(tt2, tt3, tt4) %>%
  distinct() %>%
  rownames_to_column("ROW")

dd <- crossing(ROW = tt$ROW, LOCUS = 1:100) %>%
  left_join(tt, by = "ROW")

write_tsv(dd, path="./codes/param/param_eur_afr_3pop.tsv", col_names = FALSE, quote_escape = FALSE)

# real case simulation case 1
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- c(0, 0.00001713863, 0.00004418771, 0.0001139271,
  0.0002937327, 0.0007573169, 0.001952554, 0.005034176)
N1 <- 400000
N2 <- 8000
NGE1 <- 373
NGE2 <- 441

tt1 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2)

# real case simulation case 2
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- c(0, 0.00001713863, 0.00004418771, 0.0001139271,
  0.0002937327, 0.0007573169, 0.001952554, 0.005034176)
N1 <- 500000
N2 <- 15000
NGE1 <- 373
NGE2 <- 441

tt2 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2)

tt <- tt1 %>%
  bind_rows(tt2) %>%
  distinct() %>%
  rownames_to_column("ROW")

dd <- crossing(ROW = tt$ROW, LOCUS = 1:100) %>%
  left_join(tt, by = "ROW")

write_tsv(dd, path="./codes/param/param_eur_afr_real.tsv", col_names = FALSE, quote_escape = FALSE)

# pop2 h2ge

POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
N1 <- N2 <- 100000
NGE1 <- NGE2 <- 200
H2GE1 <- 0.0007573169
H2GE2 <- c(0.00001713863, 0.00004418771, 0.0001139271,
  0.0002937327, 0.0007573169, 0.001952554, 0.005034176)

tt1 <- crossing(POP, MODEL, eQTLMODEL, H2G, N1, NGE1, N2, NGE2, H2GE1, H2GE2)

tt <- tt1 %>%
  distinct() %>%
  rownames_to_column("ROW")

dd <- crossing(ROW = tt$ROW, LOCUS = 1:100) %>%
  left_join(tt, by = "ROW")

write_tsv(dd, path="./codes/param/param_pop2h2ge.tsv", col_names = FALSE, quote_escape = FALSE)

# proxy
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- 100000
NGE1 <- NGE2 <- 200
RHO <- c(0, 0.3, 0.6, 0.9)

tt1 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2, RHO)

tt <- tt1 %>%
  distinct() %>%
  rownames_to_column("ROW")

dd <- crossing(ROW = tt$ROW, LOCUS = 1:100) %>%
  left_join(tt, by = "ROW")

write_tsv(dd, path="./codes/param/param_proxy.tsv", col_names = FALSE, quote_escape = FALSE)

# pop2 weights
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- 100000
NGE1 <- NGE2 <- 200

tt1 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2)

tt <- tt1 %>%
  distinct() %>%
  rownames_to_column("ROW")

dd <- crossing(ROW = tt$ROW, LOCUS = 1:100) %>%
  left_join(tt, by = "ROW")

write_tsv(dd, path="./codes/param/param_pop2weights.tsv", col_names = FALSE, quote_escape = FALSE)


# pop 3 case 1
POP <- "EUR_AFR_EAS"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- N3 <- c(50000, 100000, 200000)
NGE1 <- NGE2 <- NGE3 <- 200

tt1 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2, N3, NGE3) %>%
  filter(N1 == N2 & N1 == N3)

# pop 3 case 2
POP <- "EUR_AFR_EAS"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- N3 <- 100000
NGE1 <- NGE2 <- NGE3 <- c(100, 200, 300, 400, 500)

tt2 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2, N3, NGE3) %>%
  filter(NGE1 == NGE2 & NGE1 == NGE3) 

# pop 3 case 3
POP <- "EUR_AFR_EAS"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- c(0, 0.00001713863, 0.00004418771, 0.0001139271,
  0.0002937327, 0.0007573169, 0.001952554, 0.005034176)
N1 <- N2 <- N3 <- 100000
NGE1 <- NGE2 <- NGE3 <- 200

tt3 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2, N3, NGE3)

# pop 3 case 4
POP <- "EUR_AFR_EAS"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- c(0.01, 0.05, 0.1)
H2GE <- 0.0007573169
N1 <- N2 <- N3 <- 100000
NGE1 <- NGE2 <- NGE3 <- 200

tt4 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2, N3, NGE3)

# pop 3 case 5
POP <- "EUR_AFR_EAS"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- 0.0007573169
N1 <- N2 <- N3 <- c(100000, 50000, 25000)
NGE1 <- NGE2 <- NGE3 <- 200

tt5 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2, N3, NGE3) %>%
  filter(N1 == N2 & N2 == N3)

tt <- tt1 %>%
  bind_rows(tt2, tt3, tt4, tt5) %>%
  distinct() %>%
  rownames_to_column("ROW")

dd <- crossing(ROW = tt$ROW, LOCUS = 1:100) %>%
  left_join(tt, by = "ROW")

write_tsv(dd, path="./codes/param/param_eur_afr_eas.tsv", col_names = FALSE, quote_escape = FALSE)



# real case simulation case 3
POP <- "EUR_AFR"
MODEL <- "pois2"
eQTLMODEL <- c("shared", "indep")
H2G <- 0.05
H2GE <- c(0, 0.00001713863, 0.00004418771, 0.0001139271,
  0.0002937327, 0.0007573169, 0.001952554, 0.005034176)
N1 <- 511471
N2 <- 13298
NGE1 <- 373
NGE2 <- 441

tt1 <- crossing(POP, MODEL, eQTLMODEL, H2G, H2GE, N1, NGE1, N2, NGE2)

tt <- tt1 %>%
  distinct() %>%
  rownames_to_column("ROW")

dd <- crossing(ROW = tt$ROW, LOCUS = 1:100) %>%
  left_join(tt, by = "ROW")

write_tsv(dd, path="./codes/param/param_eur_afr_real2.tsv", col_names = FALSE, quote_escape = FALSE)
