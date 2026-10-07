library(data.table)
load("../../ukbb_pheno.Rdata")

DATA <- data.frame(pheno$eid, pheno$used_in_PC, pheno$sbp, pheno$ethnic, pheno$bmi, pheno$sex, pheno$birth_year, pheno$PC1, pheno$PC2, pheno$PC3, pheno$PC4, pheno$PC5)

colnames(DATA) <- c("eid", "used_in_PC", "sbp", "ethnic", "bmi", "sex_male", "birth_year", "PC_1", "PC_2", "PC_3", "PC_4", "PC_5")
DATA <- DATA[!is.na(DATA$used_in_PC), ]
DATA <- DATA[DATA$used_in_PC == "Yes", -2]

# Match genotype variants to the GWAS selection.
GWAS <- fread("GWAS_BMI_SAS.csv")

library(genio)
fam <- read_fam("BMI_geno_SAS.fam")
bim <- read_bim("BMI_geno_SAS.bim")
chp <- paste(bim$chr, ":", bim$pos, sep = "")
G <- read_bed("BMI_geno_SAS.bed",bim$id, fam$id)
G <- G[chp %in% GWAS$chp, ]
gc()
G <- t(G)
G <- cbind(as.numeric(rownames(G)), G)
colnames(G)[1] <- "eid"

BMI_DATA <- merge(DATA, G, by = "eid")

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
write.csv(pval, "pval_BMI_SAS.csv", row.names = F)
write.csv(beta, "beta_BMI_SAS.csv", row.names = F)

DATA_EUR <- BMI_DATA[BMI_DATA$ethnic %in% ethnic[[1]], -3]
write.csv(DATA_EUR, file = "BMI_SBP_EUR_SAS_SEL.csv", row.names = F)

DATA_SAS <- BMI_DATA[BMI_DATA$ethnic %in% ethnic[[2]], -3]
write.csv(DATA_SAS, file = "BMI_SBP_AFR_SAS_SEL.csv", row.names = F)

DATA_SAS <- BMI_DATA[BMI_DATA$ethnic %in% ethnic[[3]], -3]
write.csv(DATA_SAS, file = "BMI_SBP_SAS_SEL.csv", row.names = F)
