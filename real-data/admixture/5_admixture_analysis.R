setwd("/project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/admixture_results/")
options(stringsAsFactors=F)

admix.dat <- read.table("K3_rep2.Q")
fam <- read.table("../1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG_LD-LD_200kb-1kbwin-0.3.fam")
key <- read.csv("/project/nmancuso_8/data/GENOA/processed/keys/AA_dbGap_GEO.csv")
study.key <- read.table("/project/nmancuso_8/data/GENOA/FUSION/aa_441_pt.id")
    study.key$dbgap_ID <- key[match(study.key$V1, key$GEO_ID2), "dbgap_PIN2"]
    study.key <- study.key[,2:3];colnames(study.key)[1] <- "GEO_ID"

sum_i=function(i) {
  return(sum(1:i))
}

GRMreader=function(filenm, flag=1) {
  xDt.bin <- paste0(filenm, ".grm.bin")
  xDt.nfl <- paste0(filenm, ".grm.N.bin")
  xDt.gid <- paste0(filenm, ".grm.id")

  xDt.id <- read.table(xDt.gid)
  xDt.n <- nrow(xDt.id)
  xDt.grm <- readBin(file(xDt.bin, "rb"), n=xDt.n*(xDt.n+1)/2, what=numeric(0), size=4)

  sn <- sapply(1:xDt.n, sum_i)
  off <- xDt.grm[-sn]
  diag <- xDt.grm[sn]
  if(flag==1) return(list(diag=diag, off=off, n=xDt.n))
  else {
    xDt.mat <- matrix(data=NA, nrow=xDt.n, ncol=xDt.n)
    xDt.mat[upper.tri(xDt.mat)] <- off
    xDt.mat <- t(xDt.mat)
    xDt.mat[upper.tri(xDt.mat)] <- off
    diag(xDt.mat) <- diag
    xDt.mat <- as.data.frame(xDt.mat)
    colnames(xDt.mat) <- xDt.id$V2
    rownames(xDt.mat) <- xDt.id$V2
    return(list(mat=xDt.mat, n=xDt.n))
  }
}

grm <- GRMreader("/project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/kinship_results/GENOA_AA", flag=0)$mat
study.grm <- grm[study.key$dbgap_ID, study.key$dbgap_ID]

admix.dat <- cbind(fam[,1:2], admix.dat)
colnames(admix.dat)[1:2] <- c("FID", "IID")

for (i in 3:ncol(admix.dat)) {
    pop <- admix.dat[which.max(admix.dat[,i]),"FID"]
    if (pop %in% c("GENOA_AA", "ACB", "GWD", "ESN", "MSL", "YRI", "LWK", "ASW")) {
        anc="AFR"
    } else if (pop %in% c("GENOA_EA", "GBR", "IBS", "CEU", "TSI")) {
        anc="EUR"
    } else if (pop %in% c("PUR", "CLM", "PEL", "MXL")) {
        anc="AMR"
    } else {
        anc=""
    }
    colnames(admix.dat)[i] <- anc
}

high.afr.genoa <- admix.dat[admix.dat$FID=="GENOA_AA" & admix.dat$AFR > 0.75,]
high.afr.genoa.grm <- grm[high.afr.genoa$IID, high.afr.genoa$IID]
related_pairs_0.05 <- as.data.frame(matrix(ncol=2, nrow=0))
for(i in 1:nrow(high.afr.genoa.grm)) {
    relateds <- which(high.afr.genoa.grm[i,] > 0.05)
    related_pairs_0.05[(nrow(related_pairs_0.05)+1):(nrow(related_pairs_0.05)+length(relateds)),] <- as.data.frame(cbind(rep(rownames(high.afr.genoa.grm)[i], length(relateds)), colnames(high.afr.genoa.grm)[relateds]))
}
related_pairs_0.05 <- related_pairs_0.05[which(related_pairs_0.05[,1]!=related_pairs_0.05[,2]),]
relateds <- c()
while (nrow(related_pairs_0.05) > 0) {
    remove.rel <- names(sort(table(related_pairs_0.05[,1]), decreasing=T)[1])
    relateds[length(relateds)+1] <- remove.rel
    related_pairs_0.05 <- related_pairs_0.05[-which(related_pairs_0.05[,1]==remove.rel | related_pairs_0.05[,2]==remove.rel),]
}

high.afr.unrel.genoa <- high.afr.genoa[-which(high.afr.genoa$IID %in% relateds),]
high.afr.unrel.genoa$GEO.id <- key[match(high.afr.unrel.genoa$IID, key$dbgap_PIN2),"GEO_ID2"]
high.afr.unrel.genoa.grm <- grm[high.afr.unrel.genoa$IID, high.afr.unrel.genoa$IID]

write.csv(high.afr.genoa, "AA_GENOA_0.75-AFR-anc.csv", quote=F, row.names=F)

genoa.aa <- admix.dat[admix.dat$IID %in% study.key$dbgap_ID,]
write.csv(genoa.aa, "original_AA_GENOA-AFR-anc_unrelateds.csv", quote=F, row.names=F)

svg("study_AA_GENOA_Afr-anc.svg", height=8, width=12)
    hist(genoa.aa$AFR, breaks=seq(0.45,1,0.05))
dev.off()
svg("highAfr-anc-unrel_AA_GENOA_Afr-anc.svg", height=8, width=12)
    hist(high.afr.unrel.genoa[!is.na(high.afr.unrel.genoa$GEO.id),"AFR"], breaks=seq(0.75,1,0.05))
dev.off()
