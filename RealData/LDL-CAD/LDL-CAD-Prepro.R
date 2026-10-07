library(data.table)
ICD_Pheno <- fread("../../../UKB/phenotypes/ICD_Phenos/UKB_PHENOME_20191022.txt")
colnames(ICD_Pheno)[1] <- "eid"

CAD <- ICD_Pheno$X411.4
CAD[is.na(CAD)] <- FALSE
CAD <- data.frame(ICD_Pheno$eid, CAD)
colnames(CAD) <- c("eid", "CAD")

rm(ICD_Pheno)
gc()

load("../../../ukbb_pheno.Rdata")

DATA <- data.frame(pheno$eid, pheno$used_in_PC, pheno$ethnic, pheno$ldl, pheno$sex, pheno$birth_year, pheno$bmi, pheno$PC1, pheno$PC2, pheno$PC3, pheno$PC4, pheno$PC5)

colnames(DATA) <- c("eid", "used_in_PC", "ethnic", "ldl", "sex_male", "birth_year", "bmi", "PC_1", "PC_2", "PC_3", "PC_4", "PC_5")
DATA <- DATA[!is.na(DATA$used_in_PC), ]
DATA <- DATA[DATA$used_in_PC == "Yes", -2]

CAD_DATA <- merge(CAD, DATA, by = "eid")

# Match genotype variants to the GWAS selection.
GWAS <- fread("GWAS_LDL_AFR.csv")

library(genio)
fam <- read_fam("LDL_geno_AFR.fam")
bim <- read_bim("LDL_geno_AFR.bim")
chp <- paste(bim$chr, ":", bim$pos, sep = "")
G <- read_bed("LDL_geno_AFR.bed",bim$id, fam$id)
G <- G[chp %in% GWAS$chp, ]
gc()
G <- t(G)
G <- cbind(as.numeric(rownames(G)), G)
colnames(G)[1] <- "eid"

CAD_DATA <- merge(CAD_DATA, G, by = "eid")

rm(G)
gc()

ethnic <- list(c("British", "Any other white background", "Irish"),
               c("African", "Caribbean", "Any other Black background"),
               c("Indian", "Pakistani", "Bangladeshi"))

chp <- chp[chp %in% GWAS$chp]

pval <- c()
beta <- c()
for (j in 1:length(chp)) {
  pval[j] <- GWAS$pval[GWAS$chp == chp[j]]
  beta[j] <- GWAS$eff[chp == chp[j]]
}

pval <- data.frame(pval)
beta <- data.frame(beta)
write.csv(pval, "pval_LDL_AFR.csv", row.names = F)
write.csv(beta, "beta_LDL_AFR.csv", row.names = F)

DATA_EUR <- CAD_DATA[CAD_DATA$ethnic %in% ethnic[[1]], -3]
write.csv(DATA_EUR, file = "LDL_CAD_EUR_AFR_SEL.csv", row.names = F)

DATA_SAS <- CAD_DATA[CAD_DATA$ethnic %in% ethnic[[3]], -3]
write.csv(DATA_SAS, file = "LDL_CAD_SAS_AFR_SEL.csv", row.names = F)

DATA_AFR <- CAD_DATA[CAD_DATA$ethnic %in% ethnic[[2]], -3]
write.csv(DATA_AFR, file = "LDL_CAD_AFR_SEL.csv", row.names = F)
