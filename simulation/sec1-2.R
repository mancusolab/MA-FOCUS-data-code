library(tidyverse)
library(broom)
PIP.threshold <- 0.9

load("data/eur_afr.RDat")
# Across all simulation scenarios where causal eQTLs were independent across populations,
# we found MA-FOCUS reported higher PIPs for causal genes than the baseline approach
# (0.62 compared with 0.45; P<110-15)
tmp <- dd %>%
  filter(eqtl.model == "indep") %>%
  filter(true_model %in% 1)
mean(tmp$ME.pip)
mean(tmp$meta.pip)
tmp <- dd %>%
  filter(eqtl.model == "indep") %>%
  filter(true_model %in% 1) %>%
  pivot_longer(cols = c(ME.pip, meta.pip))
tidy(lm(value ~ name + n1 + nge1  + h2g + h2ge, tmp))
1.81e-39/2

# smaller credible sets (4.89 compared to 6.62; P<110-15),
tmp <- dd %>%
  filter(eqtl.model == "indep") %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(GS.size = GS.size[1]) %>%
  mutate(name = paste0(name, "GS")) %>%
  left_join(dd %>%
      select(sim, locus, n1, nge1, h2ge, h2g) %>%
      distinct(),
    by = c("sim", "locus")) 
mean(filter(tmp, name == "ME.pipGS")$`GS.size`)
mean(filter(tmp, name == "meta.pipGS")$`GS.size`)

tidy(lm(GS.size ~ name + n1 + nge1 + h2g + h2ge, tmp))
4.26/2

# higher sensitivity (88.31% compared to 81.31%)
dd %>%
  filter(eqtl.model == "indep") %>%
  # filter(n1 == 100000) %>%
  # filter(h2g == 0.05 & h2ge == 0.0007573169) %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model)) %>%
  # select(n1, name, locus, true_model, eqtl.model) %>%
  group_by(name) %>%
  summarize(freq = sum(true_model) / n()) 

tmp <- dd %>%
  filter(eqtl.model == "indep") %>%
  # filter(n1 == 100000) %>%
  # filter(h2g == 0.05 & h2ge == 0.0007573169) %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model)) %>%
  left_join(dd %>%
      select(sim, locus, n1, nge1, h2ge, h2g) %>%
      distinct(),
    by = c("sim", "locus")) 
tidy(glm(true_model ~ name + n1 + nge1 + h2ge +h2g, tmp, family = binomial("logit")))
1.87/2

# from 200 to 400, improved sensitivity by 6% from 91% to 97%
dd %>%
  filter(eqtl.model == "indep") %>%
  filter(n1 == 100000) %>%
  filter(h2g == 0.05 & h2ge == 0.0007573169) %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name, nge1) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model)) %>%
  # select(n1, name, locus, true_model, eqtl.model) %>%
  group_by(sim, name, nge1) %>%
  summarize(freq = sum(true_model) / n()) %>%
  mutate(name = paste0(name, "freq")) %>%
  pivot_wider(names_from = name, values_from = freq) 

# from 100,000 to 200,000, only increased sensitivity by 2% to 93%
dd %>%
  filter(eqtl.model == "indep") %>%
  filter(nge1 == 200) %>%
  filter(h2g == 0.05 & h2ge == 0.0007573169) %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name, n1) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model)) %>%
  # select(n1, name, locus, true_model, eqtl.model) %>%
  group_by(sim, name, n1) %>%
  summarize(freq = sum(true_model) / n()) %>%
  mutate(name = paste0(name, "freq")) %>%
  pivot_wider(names_from = name, values_from = freq) 

# performance improved as GWAS and eQTL sample sizes increased
# PIP
tmp <- dd %>%
  filter(eqtl.model == "indep") %>%
  filter(true_model %in% 1) %>%
  pivot_longer(cols = c(ME.pip, meta.pip))
tidy(lm(value ~ name + n1 + nge1  + h2g + h2ge, tmp))
# CGS size
tmp <- dd %>%
  filter(eqtl.model == "indep") %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(GS.size = GS.size[1]) %>%
  mutate(name = paste0(name, "GS")) %>%
  # pivot_wider(names_from = name, values_from = GS.size) %>%
  left_join(dd %>%
      select(sim, locus, n1, nge1, h2ge, h2g) %>%
      distinct(),
    by = c("sim", "locus"))
tidy(lm(GS.size ~ name + n1 + nge1 + h2g + h2ge, tmp))
# sensitivity
tmp <- dd %>%
  filter(eqtl.model == "indep") %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model)) %>%
  left_join(dd %>%
      select(sim, locus, n1, nge1, h2ge, h2g) %>%
      distinct(),
    by = c("sim", "locus")) 
tidy(glm(true_model ~ name + n1 + nge1 + h2ge +h2g, tmp, family = binomial("logit")))


# re-performed these simulations assuming that the causal eQTLs are shared across populations
tmp <- dd %>%
  filter(eqtl.model == "shared") %>%
  filter(true_model %in% 1)
mean(tmp$ME.pip)
mean(tmp$meta.pip)
tmp <- dd %>%
  filter(eqtl.model == "shared") %>%
  filter(true_model %in% 1) %>%
  pivot_longer(cols = c(ME.pip, meta.pip))
tidy(lm(value ~ name + n1 + nge1  + h2g + h2ge, tmp))

tmp <- dd %>%
  filter(eqtl.model == "shared") %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(GS.size = GS.size[1]) %>%
  mutate(name = paste0(name, "GS")) %>%
  right_join(dd %>%
      select(sim, locus, n1, nge1, h2ge, h2g) %>%
      distinct(),
    by = c("sim", "locus"))
mean(filter(tmp, name == "ME.pipGS")$`GS.size`)
mean(filter(tmp, name == "meta.pipGS")$`GS.size`)

tidy(lm(GS.size ~ name + n1 + nge1 + h2g + h2ge, tmp))

dd %>%
  filter(eqtl.model == "shared") %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model)) %>%
  group_by(name) %>%
  summarize(freq = sum(true_model) / n()) 

tmp <- dd %>%
  filter(eqtl.model == "shared") %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model)) %>%
  left_join(dd %>%
      select(sim, locus, n1, nge1, h2ge, h2g) %>%
      distinct(),
    by = c("sim", "locus")) 
tidy(glm(true_model ~ name + n1 + nge1 + h2ge +h2g, tmp, family = binomial("logit")))

# this performance advantage was slightly attenuated compared to the independent eQTL setting
tmp <- dd %>%
  mutate(diff = ME.pip - meta.pip) %>%
  filter(true_model %in% 1) 
tidy(lm(diff ~ n1 + nge1 + h2ge + eqtl.model + h2g, tmp))

tmp <- dd %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model),
    GS.size = GS.size[1]) %>%
  pivot_wider(names_from = name, values_from = GS.size) %>%
  mutate(diff = ME.pip - meta.pip) %>%
  left_join(dd %>%
      select(sim, locus, eqtl.model, n1, nge1, h2ge, h2g) %>%
      distinct(),
    by = c("sim", "locus"))
tidy(lm(diff ~ n1 + nge1 + h2ge + eqtl.model + h2g, tmp))

tmp <- dd %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model=sum(true_model)) %>%
  group_by(sim, name) %>%
  summarize(freq = sum(true_model) / n()) %>%
  pivot_wider(names_from = name, values_from = freq) %>%
  mutate(diff = ME.pip - meta.pip) %>%
  left_join(dd %>%
    select(sim, n1, nge1, h2ge, h2g, eqtl.model) %>%
      distinct(),
    by = c("sim"))

mean(filter(tmp, eqtl.model == "indep")$diff)
mean(filter(tmp, eqtl.model == "shared")$diff)

tidy(lm(diff ~ n1 + nge1 + h2ge + eqtl.model + h2g, tmp))

# fixed total GWAS size 
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
      pivot_longer(c(ME, SE))) %>%
  mutate(n = ifelse(group == "100k", 100000, 200000))
tidy(lm(value ~ name + group, tmp))
tmp %>%
  group_by(group, name) %>%
  summarize(n = mean(value)) %>%
  group_by(name) %>%
  summarize(n = mean(n))

tmp %>%
  group_by(name) %>%
  summarize(n = mean(value))
1.22e-2/2

tmp <- dd1 %>%
  filter(n1 == 50000) %>%
  mutate(value = ME.pip) %>%
  group_by(sim, locus) %>%
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
      group_by(sim, locus) %>%
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
      filter(n1 == 100000) %>%
      mutate(value = ME.pip) %>%
      group_by(sim, locus) %>%
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
          filter(n1 == 200000) %>%
          mutate(value = pop1.pip) %>%
          group_by(sim, locus) %>%
          arrange(desc(value)) %>%
          mutate(value = value/sum(value),
            cumpip = cumsum(value),
            GS.size = which(cumpip >= PIP.threshold)[1],
            rownum = 1:n()) %>%
          filter(GS.size >= rownum) %>%
          summarize(SE = GS.size[1]) %>%
          ungroup() %>%
          select(locus,SE),
        by = c("locus")) %>%
      mutate(group = "200k") %>%
      pivot_longer(c(ME, SE))) %>%
  mutate(n = ifelse(group == "100k", 100000, 200000))
tidy(lm(value ~ name + group, tmp))
tmp %>%
  group_by(name) %>%
  summarize(n = mean(value))

tmp <- dd1 %>%
  filter(n1 == 50000) %>%
  mutate(value = ME.pip) %>%
  group_by(locus) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= PIP.threshold)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(ME = sum(true_model)) %>%
  ungroup() %>%
  select(locus, ME) %>%
  left_join(dd1 %>%
      filter(n1 == 100000) %>%
      mutate(value = pop1.pip) %>%
      group_by(locus) %>%
      arrange(desc(value)) %>%
      mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
      filter(GS.size >= rownum) %>%
      summarize(SE = sum(true_model)) %>%
      ungroup() %>%
      select(locus, SE),
    by = c("locus")) %>%
  mutate(group = "100k") %>%
  pivot_longer(c(ME, SE)) %>%
  bind_rows(dd1 %>%
      filter(n1 == 100000) %>%
      mutate(value = ME.pip) %>%
      group_by(locus) %>%
      arrange(desc(value)) %>%
      mutate(value = value/sum(value),
        cumpip = cumsum(value),
        GS.size = which(cumpip >= PIP.threshold)[1],
        rownum = 1:n()) %>%
      filter(GS.size >= rownum) %>%
      summarize(ME = sum(true_model)) %>%
      ungroup() %>%
      select(locus, ME) %>%
      left_join(dd1 %>%
          filter(n1 == 200000) %>%
          mutate(value = pop1.pip) %>%
          group_by(locus) %>%
          arrange(desc(value)) %>%
          mutate(value = value/sum(value),
            cumpip = cumsum(value),
            GS.size = which(cumpip >= PIP.threshold)[1],
            rownum = 1:n()) %>%
          filter(GS.size >= rownum) %>%
          summarize(SE = sum(true_model)) %>%
          ungroup() %>%
          select(locus, SE),
        by = c("locus")) %>%
      mutate(group = "200k") %>%
      pivot_longer(c(ME, SE))) %>%
  mutate(n = ifelse(group == "100k", 100000, 200000))

mean(filter(tmp, name == "ME")$value)
mean(filter(tmp, name == "SE")$value)

tidy(glm(value ~ name + group, tmp, family=binomial("logit")))

# This relative performance advantage held when we compared two- to three-population scenarios
load("data/eur_afr_3pop.RDat")
dd1 <- dd
load("data/eur_afr_eas.RDat")
dd2 <- dd

dd23 <- dd2 %>%
  filter(h2ge == 0.0007573169 & h2g == 0.05 & nge1 == 200 & n1 %in% c(100000, 50000, 25000) & eqtl.model %in% "indep") %>%
  mutate(pop = "3-Pop") %>%
  mutate(sim = as.numeric(factor(sim))) %>%
  bind_rows(dd1 %>%
      filter(h2ge == 0.0007573169 & h2g == 0.05 & nge1 == 200 & eqtl.model %in% "indep") %>%
      mutate(pop = "2-Pop") %>%
      mutate(sim = as.numeric(factor(sim)) *5)) %>%
  mutate(totalGS = ifelse(n1 %in% c(100000, 150000), 300000,
    ifelse(n1 %in% c(75000, 50000), 150000,
      ifelse(n1 %in% c(25000, 37500), 75000, NA)))) %>%
  mutate(sim = as.numeric(factor(sim))) %>%
  select(sim, locus, pop, totalGS, true_model, ME.pip)

tmp <- dd23 %>%
  filter(true_model == 1)

tidy(lm(ME.pip ~ pop+ totalGS, tmp))

tmp <- dd23 %>%
  group_by(sim, locus, pop, totalGS) %>%
  arrange(desc(ME.pip)) %>%
  mutate(value = ME.pip) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(value = GS.size[1])

tidy(lm(value ~ pop+ totalGS, tmp))

tmp <- dd23 %>%
  mutate(value = ME.pip) %>%
  group_by(totalGS, locus, pop) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model = sum(true_model)) 

mean(filter(tmp, pop == "2-Pop")$true_model)
mean(filter(tmp, pop == "3-Pop")$true_model)

tidy(glm(true_model ~ pop + totalGS, tmp, family=binomial("logit")))


# the cis-SNP heritability of gene expression (cis-hg2) and the proportion of
# trait heritability attributable to a causal gene (hGE2). Across architectures,
# MA-FOCUS significantly outperformed the baseline

load("data/eur_afr.RDat")
dd1 <- dd %>%
  filter(eqtl.model == "indep" &
      h2ge %in% c(0, 0.00001713863, 0.0001139271, 0.0007573169, 0.005034176),
    nge1== 200 & n1 == 100000)
tmp <- dd1 %>%
  filter(true_model %in% 1) %>%
  pivot_longer(cols=c(ME.pip, meta.pip))

tidy(lm(value ~ name + h2g + h2ge, tmp))

tmp <- dd1 %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, h2g, h2ge, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(GS.size = GS.size[1]) %>%
  mutate(name = paste0(name, "GS")) %>%
  pivot_wider(names_from = name, values_from = GS.size) %>%
  pivot_longer(cols=c(ME.pipGS, meta.pipGS))

tidy(lm(value ~  name + h2g + h2ge, tmp))

tmp <- dd1 %>% 
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, h2g, h2ge, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(value=sum(true_model))

mean(filter(tmp, name == "ME.pip")$value)
mean(filter(tmp, name == "meta.pip")$value)

tidy(glm(value ~ name + h2g +h2ge, tmp, family=binomial("logit")))

# MA-FOCUS returned larger PIPs for the null model and smaller credible sets
# on average compared with the baseline
dd1 <- dd %>%
  filter(eqtl.model == "indep" &
      h2ge %in% 0, nge1== 200 & n1 == 100000)
tmp <- dd1 %>%
  filter(true_model %in% 1) %>%
  pivot_longer(cols=c(ME.pip, meta.pip))

tidy(lm(value ~ name, tmp))

tmp <- dd1 %>%
  pivot_longer(cols = c(ME.pip, meta.pip)) %>%
  group_by(sim, h2g, h2ge, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(GS.size = GS.size[1]) %>%
  mutate(name = paste0(name, "GS")) %>%
  pivot_wider(names_from = name, values_from = GS.size) %>%
  pivot_longer(cols=c(ME.pipGS, meta.pipGS))

tidy(lm(value ~  name, tmp))


#  we found that MA-FOCUS consistently reported higher PIPs for causal genes and
# smaller 90% credible sets compared with the baseline
load("data/pop2h2ge.RDat")
dd1 <- dd %>%
  filter(h2ge2 %in% c(0.00001713863, 0.0001139271, 0.0007573169, 0.005034176) &
      eqtl.model == "indep")
tmp <- dd1 %>%
  filter(true_model == 1) %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  select(h2ge2, name, value)

tidy(lm(value ~ name + h2ge2, tmp))

tmp <- dd1 %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(value = GS.size[1],
    h2ge2 = first(h2ge2)) %>%
  ungroup() %>%
  select(sim, locus, h2ge2, name, value)

tidy(lm(value ~ name + h2ge2, tmp))

load("./data/eur_afr_real.RDat")
dd1 <- dd %>%
  filter(h2ge %in% c(0, 0.00001713863, 0.0001139271, 0.0007573169, 0.005034176)) %>%
  filter(eqtl.model == "indep")

tmp <- dd1  %>%
  filter(true_model == 1) %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  select(h2ge, name, value)

tidy(lm(value ~ name + h2ge, tmp))

tmp <- dd1 %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(value = GS.size[1],
    h2ge = first(h2ge)) %>%
  ungroup() %>%
  select(sim, locus, name, h2ge, value) %>%
  distinct() 

tidy(lm(value ~ name + h2ge, tmp))

tmp <- dd1 %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, h2ge, name) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model = sum(true_model))

mean(filter(tmp, name == "ME.pip")$true_model)
mean(filter(tmp, name == "meta.pip")$true_model)

tidy(glm(true_model ~ name + h2ge, tmp, family = binomial("logit")))
3.48e-1/2

# We again observed that MA-FOCUS outperformed the baseline approach as well as
# single-pop FOCUS on AFR across all metrics
load("./data/eur_afr.RDat")
tmpdd6 <- dd %>%
  filter(nge1 == 200 & h2ge == 0.0007573169 & h2g == 0.05 & n1 == 100000)
load("./data/proxy.RDat")

dd6 <- dd %>%
  bind_rows(tmpdd6) %>%
  filter(eqtl.model == "indep") %>%
  mutate(rho = ifelse(is.na(rho), 1, rho),
    proxy_tissue_pop = "pop2")

tmp <- dd6 %>%
  filter(true_model == 1) %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  select(rho, name, value) 
tidy(lm(value ~ name + rho, tmp))

tmp <- dd6 %>%
  filter(true_model == 1) %>%
  pivot_longer(c(ME.pip, pop2.pip)) %>%
  select(rho, name, value)
tidy(lm(value ~ name + rho, tmp))

tmp <- dd6 %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name, rho) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(value = GS.size[1],
    rho = first(rho)) %>%
  ungroup() %>%
  select(sim, locus, rho, name, value)
tidy(lm(value ~ name + rho, tmp))

tmp <- dd6 %>%
  pivot_longer(c(ME.pip, pop2.pip)) %>%
  group_by(sim, locus, name, rho) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(value = GS.size[1],
    rho = first(rho)) %>%
  ungroup() %>%
  select(sim, locus, rho, name, value)
tidy(lm(value ~ name + rho, tmp))


tmp <- dd6 %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name, rho) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model = sum(true_model))

tidy(lm(true_model ~ name + rho,tmp))

tmp <- dd6 %>%
  pivot_longer(c(ME.pip, pop2.pip)) %>%
  group_by(sim, locus, name, rho) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model = sum(true_model))

tidy(lm(true_model ~ name + rho,tmp))
1.71e-1/2

load("./data/eur_afr.RDat")
tmpdd7 <- dd %>%
  filter(nge1 == 200 & h2ge == 0.0007573169 & h2g == 0.05 & n1 == 100000)
load("./data/pop2weights.RDat")

dd7 <- dd %>%
  mutate(pop = "EUR") %>%
  bind_rows(tmpdd7 %>%
      mutate(pop = "AFR")) %>%
  filter(eqtl.model == "indep")

tmp <- dd7 %>%
  filter(true_model == 1) %>%
  filter(pop == "EUR") %>%
  pivot_longer(cols=c(ME.pip, meta.pip, pop2.pip))

tidy(lm(value ~ name, tmp))
6.47e-1/2

tmp <- dd7 %>%
  pivot_longer(c(ME.pip, meta.pip)) %>%
  group_by(sim, locus, name, pop) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(value = GS.size[1],
    pop  = first(pop)) %>%
  ungroup() %>%
  select(sim, locus, pop, name, value) %>%
  filter(pop == "EUR")

tidy(lm(value ~ name, tmp))
7.45e-9/2

tmp <- dd7 %>%
  pivot_longer(c(ME.pip, pop2.pip)) %>%
  group_by(sim, locus, name, pop) %>%
  arrange(desc(value)) %>%
  mutate(value = value/sum(value),
    cumpip = cumsum(value),
    GS.size = which(cumpip >= 0.9)[1],
    rownum = 1:n()) %>%
  filter(GS.size >= rownum) %>%
  summarize(true_model = sum(true_model))

tidy(glm(true_model ~ name, tmp, family=binomial("logit")))
7.14e-1/2
