options(stringsAsFactors=F)
setwd("/project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/")

fam <- read.table("1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG_LD-LD_200kb-1kbwin-0.3.fam", header=F)
study.EA <- read.table("/project/nmancuso_8/data/GENOA/FUSION/ea_373_pt.id")
study.AA <- read.table("/project/nmancuso_8/data/GENOA/FUSION/aa_441_pt.id")
study.ids <- rbind(study.EA, study.AA)
EA.key <- read.csv("/project/nmancuso_8/data/GENOA/processed/keys/EA_dbGap_GEO.csv", header=T)
AA.key <- read.csv("/project/nmancuso_8/data/GENOA/processed/keys/AA_dbGap_GEO.csv", header=T)
key <- as.data.frame(cbind(c(rep("GENOA_EA", nrow(EA.key)), rep("GENOA_AA", nrow(AA.key))), c(EA.key$dbgap_PIN2, AA.key$dbgap_PIN2), c(EA.key$GEO_ID1, AA.key$GEO_ID2)))
study.ids$V1 <- key[match(study.ids$V2, key$V3), "V2"]

# make admixture files
for (rep in 1:30) {
    Q <- read.table(paste0("admixture_results/K3_rep",rep,".Q"))
    dat <- cbind(fam, Q)
    dat <- dat[dat$V2 %in% study.ids$V1, 7:9]
    write.table(dat, paste0("admixture_results/GENOA-only_K3_rep",rep,".Q"), quote=F, row.names=F, col.names=F)
}
ind2pop <- fam[fam$V2 %in% study.ids$V1,1]
write.table(ind2pop, "admixture_results/GENOA_pong.ind2pop", quote=F, row.names=F, col.names=F)
