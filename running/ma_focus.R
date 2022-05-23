suppressMessages(library('corpcor')) # for pseudoinverse() function
suppressMessages(library('dplyr')) # for credible set annotation
suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("mvnfast"))
suppressMessages(library('readr')) 

LOG_INFO <- "INFO"
LOG_WARNING <- "WARNING"
LOG_ERROR <- "ERROR"

LOG <- function(x, level=LOG_INFO) {
  str <- paste0("[%D %H:%M:%S - ", level, "]")
  cat(format(Sys.time(), str), x, "\n")
}

log_sum_exp <- function(a, b) {
  x <- c(a, b)
  # Computes log(sum(exp(x))
  # Uses offset trick to avoid numeric overflow: http://jblevins.org/notes/log-sum-exp
  if ( max(abs(x)) > max(x) ) {
    offset <- min(x)
  } else {
    offset <- max(x)
  }
  log(sum(exp(x - offset))) + offset
}

annotate_cred_set <- function(df_s, prb=0.90) {
  # the call to abs  function may seem weird, but there is numerical issues with 1 - cumsum for strange reason.
  # w/o it sometimes small negative numbers appear and throw off CS computation
  df_s <- df_s %>% arrange(-PIP.ME) %>%
    mutate(NPIP = PIP.ME / sum(PIP.ME),
      CSUM = cumsum(NPIP),
      IN.ME.CRED.SET = FALSE)
  for (idx in 1:nrow(df_s)) {
    df_s$IN.ME.CRED.SET[idx] <- TRUE
    if (df_s$CSUM[idx] >= prb) {
      break
    }
  }
  df_s %>% select(-c(NPIP, CSUM))
}

annotate_cred_set_singlePopFOCUS <- function(df_s, pop_number, prb=0.90) {
  # the call to abs  function may seem weird, but there is numerical issues with 1 - cumsum for strange reason.
  # w/o it sometimes small negative numbers appear and throw off CS computation
  df_s <- df_s %>% arrange(-get(paste("PIP",pop_number,sep="."))) %>%
    mutate(NPIP = get(paste("PIP",pop_number,sep=".")) / sum(get(paste("PIP",pop_number,sep="."))),
      CSUM = cumsum(NPIP),
      IN.CRED.SET = FALSE)
  for (idx in 1:nrow(df_s)) {
    df_s$IN.CRED.SET[idx] <- TRUE
    if (df_s$CSUM[idx] >= prb) {
      break
    }
  }
  df_s %>% select(-c(NPIP, CSUM))
}

get_independent <- function(CHR, ID, P0, P1, regions) {
  # Sometimes labels are the same bc of bugs...
  P1 = P0 + 1
  
  # compute overlaps
  g_ranges <- GRanges(seqnames = CHR, ranges= IRanges::IRanges(start = P0, end = P1, names = ID))
  r_ranges <- GRanges(seqnames = regions$CHR, ranges= IRanges::IRanges(start = regions[,2], end = regions[,3]))
  overlaps <- findOverlaps(g_ranges, r_ranges, select="arbitrary", maxgap=1e4)
  
  # some annotations don't overlap exactly
  # just get nearest for now
  if (sum(is.na(overlaps)) > 0) {
    subset <- g_ranges[is.na(overlaps)]
    miss <- nearest(subset, r_ranges, select="arbitrary")
    overlaps[is.na(overlaps)] <- miss
  }
  
  # prettify regions and output them as 'blocks'
  pranges <- paste0(regions$CHR, ":", regions[,2], "..", regions[,3])
  pranges[overlaps]
}

get_local_params <- function(weight.files, IDs, models, genos) {
  # load weights into matrix after QCing...
  Wlist <- lapply(1:length(weight.files), function(i) {
    load(weight.files[i])
    wgt.matrix[is.na(wgt.matrix)] = 0
    
    # Match up the SNPs and weights
    m = match(snps[,2], genos$bim[,2])
    m.keep = !is.na(m)
    snps = snps[m.keep,]
    wgt.matrix = wgt.matrix[m.keep,]
    
    cur.genos = genos$bed[,m[m.keep]]
    cur.bim = genos$bim[m[m.keep],]
    
    # Flip WEIGHTS for mismatching alleles
    qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
    wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
    
    # Predict into reference
    mod = models[i]
    
    wgt <- tibble(SNP=rownames(wgt.matrix), WGT=as.double(wgt.matrix[, mod]))
    colnames(wgt)[2] <- IDs[i]
    wgt
  })
  
  W <- purrr::reduce(Wlist, full_join, by="SNP")
  snps <- W$SNP
  W$SNP <- NULL
  W[is.na(W)] <- 0
  W <- as.matrix(W)
  
  # check if we have single gene
  if (is.null(dim(W)) && length(W) > 0) {
    W <- t(t(W))
  }
  rownames(W) <- snps
  
  # scale weights and compute LD
  m <- match(rownames(W), genos$bim[, 2])
  X <- genos$bed[,m]
  S <- apply(X %*% W, 2, sd)
  if (length(S) == 1) {
    flags <- S != 0
    if (S != 0) {
      S <- t(t(1 / S))
    }
  } else {
    # drop genes with 0 genetic covariance
    flags <- S != 0
    S <- S[flags]
    W <- W[, flags]
    S <- diag(1 / S)
  }
  SW <- W %*% S
  LD <- cor(X)
  
  return (list(SW=SW, LD=LD, FLAGS=flags))
}

allele.qc = function(a1,a2,ref1,ref2) {
  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip
  
  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;
  
  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  
  return(snp)
}

## SSG - MODIFIED TO MATCH PYTHON CODE
# --- Estimate the sample correlation structure for predicted expression ---
# wmat: eQTL weight matrix for a risk region
# ldmat: LD matrix for a risk region
# intercept: bool to return the intercept variable or not
# return: list of length 2 - 1) matrix pred_expr correlation, 2) intercept variable; NULL if intercept=F, matrix (nGene x nSNP dim.) otherwise
estimate_cor <- function(wmat, ldmat, intercept=FALSE) {
  wmat <- as.matrix(wmat) ## SNPs are rows, genes are columns
  ldmat <- as.matrix(ldmat)
  wcov <- ((t(wmat) %*% ldmat) %*% wmat) ## %*% is equivalent to numpy.dot - mdot is the equivalent of >1 numpy.dots chained together (multiplication of > 2 matrices)
  scale <- diag(1 / sqrt(diag(wcov)))
  if (ncol(wmat)==1 & nrow(scale)==0) {
    scale <- wcov;scale[1,1] <- 1
  }
  wcor <- ((scale %*% wcov) %*% scale)
  if (intercept==TRUE) {
    inter = ((scale %*% t(wmat)) %*% ldmat)
    return(list(wcor, inter))
  } else {
    return(list(wcor, NULL))
  }
}

# --- Estimate prior chisq value ---
estimate_prior_chisq <- function(wmat, ldmat, zscores) {
  wmat <- as.matrix(wmat)
  ldmat <- as.matrix(ldmat)
  zscores <- as.matrix(zscores) # z-scores in a single column
  V <- ((t(wmat) %*% ldmat) %*% wmat) # calculate gene covariance matrix from weight and LD matrices
  V_scaled <- diag(1 / sqrt(diag(V))) # scale this matrix
  if (ncol(wmat)==1 & nrow(V_scaled)==0) {
    V_scaled <- V;V_scaled[1,1] <- 1
  }
  cor_mat <- ((V_scaled %*% V) %*% V_scaled)
  biased_prior_chisq <- ((t(zscores) %*% pseudoinverse(cor_mat, tol=1e-15)) %*% zscores) ### need to double-check that pseudoinverse() is equivalent to numpy.linalg.pinv
  prior_chisq <- biased_prior_chisq - nrow(zscores)
  if (prior_chisq < 0) {
    prior_chisq <- biased_prior_chisq
  }
  return(prior_chisq)
}

# --- Compute Bayes factors ---
# zscores: single column matrix of zscores for a pop.
# idx_set: indices representing the gene set being tested
# wcor: predicted expression correlation matrix for a pop. (first output of estimate_cor() function)
# prior_chisq: float prior effect-size variances for a pop. (output of estimate_prior_chisq() function)
# prb: float prior probability for a gene to be causal
# use_log: bool whether to compute the log Bayes Factor
# return: float the Bayes Factor (log Bayes Factor if use_log = TRUE)
bayes_factor <- function(zscores, idx_set, wcor, prior_chisq, prb, use_log=TRUE, tol=2.220446e-14){
  nGenes = length(zscores)
  
  # only need genes in the causal configuration using FINEMAP BF trick
  n_causals = length(idx_set)
  cur_chi2 = prior_chisq / n_causals
  
  cur_wcor = wcor[idx_set, idx_set] # subset correlation matrix to genes being tested
  cur_zscores = zscores[idx_set]
  
  # compute SVD for robust estimation
  if (n_causals > 1) {
    svd_res <- svd(cur_wcor)
    keep <- svd_res$d > tol
    nc <- sum(keep)
    cur_EIG <- svd_res$d[keep]
    cur_U <- svd_res$u[,keep]
    cur_chi2 <- prior_chisq / nc
    scaled_chisq = (t(cur_zscores) %*% cur_U)^2
  } else {
    cur_U <- 1
    cur_EIG <- cur_wcor
    scaled_chisq <- cur_zscores^2
  }
  
  # log BF
  cur_bf <- 0.5 * -sum(log(1 + cur_chi2 * cur_EIG)) +
    0.5 * sum((cur_chi2 / (1 + cur_chi2 * cur_EIG)) * scaled_chisq)
  
  # log prior
  cur_prior = n_causals * log(prb) + (nGenes - n_causals) * log(1 - prb)
  
  if (use_log==TRUE){
    return(c(cur_bf, cur_prior))
  } else {
    return(c(exp(cur_bf), exp(cur_prior)))
  }
}

# --- Perform fine-mapping ---
# zscores = list of Z-scores (1 per pop.) computed from the regression of predicted expression on outcome
# wmats = list of weight matrices (1 per pop.)
# ldmats = list of LD matrices (1 per. pop)
# loc.ID = names of the predicted genes/lncRNAs/etc
# max_causals = maximum number of causal genes to test for
# prb = prior probability of a gene being causal
# tol = tolerance flag for dropping eigenvalues
# verbose = flag for debugging
fine_map <- function(zscores, loc.ID, wmats, ldmats, max_causals=3, prb=1e-3, verbose=F) {
  # How many populations
  n.pops <- length(zscores)
  nGenes <- length(zscores[[1]])
  pips <- rep(0, nGenes)
  
  # Check that all inputs have the right number of populations
  if (length(wmats) != n.pops) {
    stop ("Number of weight matrices does not correspond to number of populations")
  }
  if (length(ldmats) != n.pops) {
    stop ("Number of LD matrices does not correspond to number of populations")
  }
  
  # Calculate prior chi square value for each pop
  prior_chisq <- rep(NA, n.pops)
  for (p in 1:n.pops) {
    prior_chisq[p] <- estimate_prior_chisq(wmats[[p]], ldmats[[p]], zscores[[p]])
  }
  
  # estimate correlation of predicted expression in cohort
  cor_exp <- list();length(cor_exp) <- n.pops # previously [wcor, swld]
  for (p in 1:n.pops) {
    cor_exp[[p]] <- estimate_cor(wmats[[p]], ldmats[[p]])
    rownames(cor_exp[[p]][[1]]) <- loc.ID
    colnames(cor_exp[[p]][[1]]) <- loc.ID
  }  
  
  m <- min(max_causals, nGenes)
  null_res <- m * log(1 - prb)
  marginal <- null_res
  
  # enumerate all subsets of genes
  gene_combos <- unlist(sapply(1:m, function(x) combn(nGenes, x, simplify=F)), recursive=F)
  # Calculate likelihoods...
  for (i in 1:length(gene_combos)) {
    local <- 0
    for (j in 1:n.pops) {
      bf_res = bayes_factor(zscores[[j]], idx_set=gene_combos[[i]], cor_exp[[j]][[1]], prior_chisq=as.numeric(prior_chisq[[j]]), prb=prb)
      local_bf <- bf_res[1]
      local_prior <- bf_res[2]
      if (j == 1) {
        local <- local + local_bf + local_prior
      } else {
        local <- local + local_bf
      }
    }
    # keep track for marginal likelihood
    marginal = log_sum_exp(local, marginal)
    # marginalize the posterior for marginal-posterior on causals
    for (j in gene_combos[[i]]) {
      if (pips[j] == 0) {
        pips[j] <- local
      } else {
        pips[j] <- log_sum_exp(pips[j], local)
      }
    }
  }
  pips = exp(pips - marginal)
  null_res = exp(null_res - marginal)
  result_df <- as.data.frame(c(null_res, pips))
  colnames(result_df) <- "PIP";rownames(result_df) <- c("NULL.MODEL",loc.ID)
  RESULT <- list(result_df, prior_chisq)
  names(RESULT) <- c("PIPs", "prior_chisq")
  return(RESULT)
}

# Perform fine-mapping...
option_list = list(
  make_option("--TWAS", action="store", default=NA, type='character',
    help="Path to TWAS test output files, separated by ':' [required]"),
  make_option("--regions", action="store", default=NA, type='character',
    help="Path to a single .bed file of gene blocks [required]"),
  make_option("--sig_regions", action="store", default=NA, type='character',
    help="Path to regions containing genome-wide significant SNPs, separated by ':' [required]"),
  make_option("--phen", action="store", default=NA, type='character',
    help="Phenotype name"),
  make_option("--out", action="store", default=NA, type='character',
    help="Path to output files [required]"),
  make_option("--ref_ld_chr", action="store", default=NA, type='character',
    help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--minp.input", action="store", default=1.0, type='character',
    help="Minimum p-value to include TWAS Z-score in analysis, separated by ':' [default: %default]"),
  make_option("--use_intercept", action="store_true", default=F, type='logical',
    help="Flag to include intercept for posterior calculation"),
  make_option("--prior_chi2", action="store", default=40, type='double',
    help="Prior chi-sq parameter for causal genes [default: %default]"),
  make_option("--prior_prob", action="store", default=NA, type='character',
    help="Prior probability for a gene to be causal [default: %default]"),
  make_option("--gencode", action="store", default=NA, type='character',
    help="Gene code map"),
  make_option("--maxcausal", action="store", default=3, type='double',
    help="number of max causal genes assumed"),
  make_option("--learn_prior_prob", action="store_true", default=F, type='logical',
    help="Flag to learn the prior causal probability using Empirical Bayes"),
  make_option("--cs_prob", action="store", default=0.9, type='double',
    help="Probability for credible set computation. Annotates genes with in/not-in flag [default: %default]"),
  make_option("--posterior_check", action="store", default=0, type='integer',
    help="Number of posterior simulations for model checking [default: %default]"),
  make_option("--tol", action="store", default=2.220446e-14, type='double',
    help="Error tolerance to determine for non-zero singular values"),
  make_option("--chr", action="store", default=NA, type='integer',
    help="Chromosome to analyze, currently only single chromosome analyses are performed [required]"),
  make_option("--perform_single_pop_FOCUS", action="store_true", default=F, type="logical",
    help="Flag to specify if single-population FOCUS should also be performed [default: %default]")
)

opt = parse_args(OptionParser(option_list=option_list))
options( digits = 5 )
chr <- opt$chr

## Extract regions to fine-map from TWAS results ----
LOG("Loading TWAS summary statistics")
# Load TWAS files, split strings
TWAS_files <- unlist(strsplit(opt$TWAS,":"))
# Load minimum TWAS P-value, split strings
minp.input <- as.numeric(unlist(strsplit(opt$minp.input,":")))
# Determine number of populations
N.pops <- length(TWAS_files)
# Load independent regions (output of LDetect)
blocks <- read.table(opt$regions, header=F);colnames(blocks) <- c("CHR","P0","P1")
blocks <- blocks[blocks$CHR==chr,]

gmap <- read_table2(opt$gencode) %>%
  mutate(chr = as.numeric(chr)) %>%
  filter(chr == as.numeric(opt$chr))

if (nrow(gmap) == 0) stop("gencode map has 0 rows")

# Read in TWAS files (weight lists), annotate with LDetect blocks the genes fall in, do some filtering for missingness
wgtlist <- list();length(wgtlist) <- N.pops
for (i in 1:N.pops) {
  wgtlist[[i]] <- read.table(TWAS_files[i], as.is=T, header=T)
  wgtlist[[i]] <- wgtlist[[i]][which(wgtlist[[i]]$CHR == chr & !is.na(wgtlist[[i]]$TWAS.P)),]
  wgtlist[[i]]$BLOCK <- get_independent(CHR=wgtlist[[i]]$CHR, ID=wgtlist[[i]]$ID, P0=wgtlist[[i]]$P0, P1=wgtlist[[i]]$P1, regions=blocks)
  wgtlist[[i]] <- wgtlist[[i]][order(wgtlist[[i]]$BLOCK, wgtlist[[i]]$P0),]
}
# Filter weight lists down to genes that are common to all datasets/populations...
if (exists("common.genes")) {rm(common.genes)}
if (N.pops > 1) {
  for (i in 2:N.pops) {
    if (!exists("common.genes")) {
      common.genes <- intersect(wgtlist[[i-1]]$ID, wgtlist[[i]]$ID)
    } else {
      common.genes <- intersect(common.genes, wgtlist[[i]]$ID)
    }
  }
} else {
  common.genes <- wgtlist[[1]]$ID
}
# ...and significant in at least one...
sig.genes <- c()
for (i in 1:N.pops) {
  dat <- wgtlist[[i]][match(common.genes, wgtlist[[i]]$ID),]
  sig.genes <- c(sig.genes, dat[dat$TWAS.P<=minp.input[i],"ID"])
}
sig.genes <- unique(sig.genes)
if (length(sig.genes)==0) {
  stop("No significant regions to fine-map, exiting")
}
# ...and have finite Z-scores...
keep.genes <- c()
for (i in 1:N.pops) {
  keep.genes <- c(keep.genes, wgtlist[[i]]$ID[which(is.finite(wgtlist[[i]]$TWAS.Z))])
}
keep.genes <- sort(names(which(table(keep.genes)==N.pops)))
keep.genes <- intersect(keep.genes, sig.genes)
# ...and subset all weight lists
for (i in 1:N.pops) {
  wgtlist[[i]] <- wgtlist[[i]][match(keep.genes, wgtlist[[i]]$ID),]
}

# Load in GWAS significant regions from provided .bed files, and concatenate all
LOG("Loading GWAS significant loci")
GWAS_sig_files <- unlist(strsplit(opt$sig_regions,":"))
sig.regions <- as.data.frame(matrix(ncol=3, nrow=0))
for (i in 1:length(GWAS_sig_files)) {
  if (nrow(sig.regions)==0) {
    sig.regions <- read.table(GWAS_sig_files[i], header=T)
  } else {
    sig.regions <- rbind(sig.regions, read.table(GWAS_sig_files[i], header=T))
  }
}
# Filter to specified chromosome only and sort
sig.regions <- sig.regions[sig.regions$CHR==chr,]
sig.regions <- sig.regions[order(sig.regions$P0),]
sig.regions.GRange <- apply(sig.regions[,2:3], paste, MARGIN=1, collapse="..")
sig.regions.GRange <- GRanges(as.vector(apply(cbind(sig.regions[,1],sig.regions.GRange), paste, MARGIN=1, collapse=":")))

# load in genotype files by chromosome, restrict to matching SNPs and combine
LOG("Loading reference LD panels")
genos <- list();length(genos) <- N.pops
LD_ref_files <- unlist(strsplit(opt$ref_ld_chr,":"))
for (i in 1:N.pops) {
  genos[[i]] = read_plink(LD_ref_files[i],impute="avg")
  genos[[i]]$bed = scale(genos[[i]]$bed)
}
# Determine number of distinct regions for fine-mapping
LOG("Determining unique genomic regions to fine-map")
all.blocks <- c()
for (i in 1:N.pops) {
  all.blocks <- c(all.blocks,wgtlist[[i]]$BLOCK)
}
unique.blocks <- unique(all.blocks)
# Determine number of genes in those regions
all.genes <- c()
for (i in 1:N.pops) {
  all.genes <- c(all.genes,wgtlist[[i]]$ID)
}
unique.genes <- unique(all.genes)

# Re-order blocks
ord <- order(unlist(lapply(strsplit(unique.blocks, ":"), function(x) as.numeric(strsplit(x[2], "\\.\\.")[[1]][1]))))
unique.blocks <- unique.blocks[ord]
unique.blocks.GRange <- GRanges(unique.blocks)

# Run ME FOCUS function
LOG("Performing fine mapping with prior probability for causality")
out.tbl = NULL
# Perform fine mapping over regions with at least one TWAS hit using the prior determined
# either from EmpBayes or whatever the user supplied
for (i in 1:length(unique.blocks)) {
  genes.idx <- wgtlist[[1]]$BLOCK == unique.blocks[i]
  ID <- wgtlist[[1]]$ID[genes.idx]
  CHR <- wgtlist[[1]]$CHR[genes.idx]
  FILE <- list();length(FILE) <- N.pops
  MODEL <- list();length(MODEL) <- N.pops
  TWAS.Z <- list();length(TWAS.Z) <- N.pops
  TWAS.P <- list();length(TWAS.P) <- N.pops
  P0 <- list();length(P0) <- N.pops
  P1 <- list();length(P1) <- N.pops
  for (j in 1:N.pops) {	
    FILE[[j]] <- wgtlist[[j]]$FILE[genes.idx]
    MODEL[[j]] <- wgtlist[[j]]$MODEL[genes.idx]
    TWAS.Z[[j]] <- wgtlist[[j]]$TWAS.Z[genes.idx]
    TWAS.P[[j]] <- wgtlist[[j]]$TWAS.P[genes.idx]
    P0[[j]] <- wgtlist[[j]]$P0[genes.idx]
    P1[[j]] <- wgtlist[[j]]$P1[genes.idx]
  }
  
  # Check that there is GWAS signal at this locus
  GWAS.TWAS.intersect <- intersect(unique.blocks.GRange[i], sig.regions.GRange)
  if (length(GWAS.TWAS.intersect)==0) {
    LOG(paste("Skipping region", unique.blocks[i], "- no GWAS signal"))
    next
  }
  if (length(unique(ID))==0) {
    LOG(paste("Skipping region", unique.blocks[i], "- no genes"))
    next
  }
  LOG(paste("Fine-mapping region", unique.blocks[i], "with", length(unique(ID)), "genes"))
  
  LD <- list();length(LD) <- N.pops
  W <- list();length(W) <- N.pops
  for (j in 1:N.pops) {
    params <- get_local_params(FILE[[j]], ID, MODEL[[j]], genos[[j]])
    LD[[j]] <- params$LD
    W[[j]] <- params$SW
  }
  
  # calculate how many genes in the region based on gencode
  if (!is.na(as.numeric(opt$prior_prob))) {
    prior_prob <- as.numeric(opt$prior_prob)
    num_genes <- 1 / as.numeric(opt$prior_prob)
  } else if (opt$prior_prob == "TWASGENE") {
    prior_prob <- 1 / (length(TWAS.Z[[1]]) + 1e-3)
    num_genes <- length(TWAS.Z[[1]])
  } else if (opt$prior_prob == "BLOCKGENE") {
    tmp_range <- GRanges(seqnames = gmap$chr, ranges= IRanges::IRanges(start = gmap$start, end = gmap$end))
    overlaps <- countOverlaps(tmp_range, unique.blocks.GRange[i])
    prior_prob <- 1 / (sum(overlaps[!is.na(overlaps)]) + 1e-3)
    num_genes <- sum(overlaps[!is.na(overlaps)])
  } else {
    stop("Prior_prob parameter cannot be recognized. Desired input: a numeric probablity, TWASGENE, or ALLGENE.")
  }
  
  LOG(paste0("Prior causal probability is ", prior_prob,
    ", and expected number of genes is ", num_genes))
  
  if (prior_prob >= 1 | prior_prob <= 0) {
    stop("Invalid causal prior. Please specify a number between 0 and 1.")
  }
  
  # Do the multi-ethnic fine-mapping
  result <- fine_map(zscores=TWAS.Z,
    loc.ID=ID,
    wmats=W,
    ldmats=LD,
    max_causals=as.numeric(opt$maxcausal),
    prb=prior_prob)
  
  # Start constructing output table (not being filled yet)
  tbl.FILES <- as.data.frame(matrix(nrow=nrow(result$PIPs), ncol=N.pops));colnames(tbl.FILES) <- paste("FILE",1:N.pops,sep=".")
  tbl.P0 <- tbl.FILES;colnames(tbl.P0) <- paste("P0",1:N.pops,sep=".")
  tbl.P1 <- tbl.FILES;colnames(tbl.P1) <- paste("P1",1:N.pops,sep=".")
  tbl.Zs <- tbl.FILES;colnames(tbl.Zs) <- paste("TWAS.Z",1:N.pops,sep=".")
  tbl.Ps <- tbl.FILES;colnames(tbl.Ps) <- paste("TWAS.P-val",1:N.pops,sep=".")
  tbl.p_chisq <- tbl.FILES;colnames(tbl.p_chisq) <- paste("prior_chi2",1:N.pops,sep=".")
  
  # filling some information here - gene ID, chromosome, LD block
  tbl <- cbind(tbl.FILES,
    ID=c(paste0("NULL.", i), as.character(ID)),
    CHR=rep(chr, nrow(result$PIPs)),
    tbl.P0,
    tbl.P1,
    BLOCK=rep(unique.blocks[i], nrow(result$PIPs)),
    BLOCKGENES=num_genes,
    tbl.Zs,
    tbl.Ps,
    tbl.p_chisq,
    PIP.ME=result$PIPs$PIP)
  
  # filling in file, positional info, TWAS info, and ME fine-mapping results
  for (j in 1:N.pops) {
    tbl[,match(paste("FILE",j,sep="."), colnames(tbl))] <- c("NULL",FILE[[j]])
    tbl[,match(paste("P0",j,sep="."), colnames(tbl))] <- c(0,P0[[j]])
    tbl[,match(paste("P1",j,sep="."), colnames(tbl))] <- c(0,P1[[j]])
    tbl[,match(paste("TWAS.Z",j,sep="."), colnames(tbl))] <- c(0,TWAS.Z[[j]])
    tbl[,match(paste("TWAS.P-val",j,sep="."), colnames(tbl))] <- c(1,TWAS.P[[j]])
    tbl[,match(paste("prior_chi2",j,sep="."), colnames(tbl))] <- rep(result$prior_chisq[[j]], nrow(result$PIPs))
  }
  # annotate with credible set inclusion
  tbl <- annotate_cred_set(tbl, opt$cs_prob)
  # Re-order genes in output table to keep consistent with single-population FOCUS
  gene.order <- c(paste("NULL",i,sep="."),ID)
  tbl <- tbl[match(gene.order, as.character(tbl$ID)),]
  # Perform single-population FOCUS, if user has chosen this
  if (opt$perform_single_pop_FOCUS==T) {
    for (j in 1:N.pops) {
      sp.result <- fine_map(zscores=list(TWAS.Z[[j]]),
        loc.ID=ID,
        wmats=list(W[[j]]),
        ldmats=list(LD[[j]]),
        max_causals=as.numeric(opt$maxcausal),
        prb=prior_prob)
      sp.PIPs <- as.data.frame(sp.result$PIPs$PIP)
      colnames(sp.PIPs) <- paste("PIP",j,sep=".")
      tbl <- cbind(tbl, sp.PIPs)
      tbl <- annotate_cred_set_singlePopFOCUS(tbl, j, opt$cs_prob)
      # re-order again
      tbl <- tbl[match(gene.order, as.character(tbl$ID)),]
      colnames(tbl)[ncol(tbl)] <- paste("IN",j,"CRED.SET",sep=".")
    }
  }
  
  # add results from this LD block to previous results (if any)
  if (is.null(out.tbl)) {
    out.tbl <- tbl
  } else {
    out.tbl <- rbind(out.tbl, tbl)
  }
}

out.tbl$PHEN <- opt$phen

# remove regions w/o finemapping and add credible set flags
if (!is.null(out.tbl)) {
  write.table(out.tbl, quote=F , row.names=F , sep='\t' , file=opt$out)
}

