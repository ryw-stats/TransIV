library(genio)

bim_eu <- read_bim("EUR_chr22.bim")
bim_af <- read_bim("AFR_chr22.bim")
bim_sa <- read_bim("SAS_chr22.bim")

GWAS <- read.table("MR_LDL_joint_2013.txt", header = T)
GWAS <- GWAS[GWAS$chr == 22, ]

S_eu <- (bim_eu$pos %in% GWAS$pos[order(GWAS$pval)[1:5000]]) & (bim_eu$pos %in% bim_af$pos) & (bim_eu$pos %in% bim_sa$pos)
S_af <- (bim_af$pos %in% GWAS$pos[order(GWAS$pval)[1:5000]]) & (bim_af$pos %in% bim_eu$pos) & (bim_af$pos %in% bim_sa$pos)
S_sa <- (bim_sa$pos %in% GWAS$pos[order(GWAS$pval)[1:5000]]) & (bim_sa$pos %in% bim_eu$pos) & (bim_sa$pos %in% bim_af$pos)

p_ref <- c()

for (j in 1:sum(S_eu)) {
  p_ref[j] <- GWAS$pval[GWAS$pos == bim_eu$pos[S_eu][j]]
}

G <- read_bed("EUR_chr22.bed", bim_eu$id, 1:120000)
G <- G[S_eu, ]

G_eu <- G

G <- read_bed("AFR_chr22.bed", bim_af$id, 1:120000)
G <- G[S_af, ]

G_af <- G

G <- read_bed("SAS_chr22.bed", bim_sa$id, 1:120000)
G <- G[S_sa, ]

G_sa <- G

S_common <- (rowMeans(G_eu) >= 0.05) & (rowMeans(G_af) >= 0.05) & (rowMeans(G_sa) >= 0.05)

G_eu <- G_eu[S_common, ]
G_af <- G_af[S_common, ]
G_sa <- G_sa[S_common, ]

p_ref <- p_ref[S_common]

gc()

G_tmp <- G_eu - rowMeans(G_eu)

Sig_eu <- matrix(NA, length(p_ref), length(p_ref))

Sig_eu <- G_tmp %*% t(G_tmp)

Sig_eu <- Sig_eu / ncol(G_eu)

gc()

G_tmp <- G_af - rowMeans(G_af)

Sig_af <- G_tmp %*% t(G_tmp)

Sig_af <- Sig_af / ncol(G_af)

gc()

G_tmp <- G_sa - rowMeans(G_sa)

Sig_sa <- G_tmp %*% t(G_tmp)

Sig_sa <- Sig_sa / ncol(G_sa)

a <- diag(Sig_eu)
R_eu <- diag(1 / sqrt(a)) %*% Sig_eu %*% diag(1 / sqrt(a))

a <- diag(Sig_af)
R_af <- diag(1 / sqrt(a)) %*% Sig_af %*% diag(1 / sqrt(a))

a <- diag(Sig_sa)
R_sa <- diag(1 / sqrt(a)) %*% Sig_sa %*% diag(1 / sqrt(a))

ref <- 1 - p_ref

sel <- rep(T, length(ref))

C_eu <- R_eu - diag(diag(R_eu))
C_af <- R_af - diag(diag(R_af))
C_sa <- R_sa - diag(diag(R_sa))

order <- order(ref)

for (i in order) {
  sel[i] <- max(c(C_eu[i, sel], C_af[i, sel], C_sa[i, sel])) < 0.9
}

G_eu <- t(G_eu[sel, ])
G_af <- t(G_af[sel, ])
G_sa <- t(G_sa[sel, ])

save(G_eu, file = "SimData_EU.Rdata")
save(G_af, file = "SimData_AF.Rdata")
save(G_sa, file = "SimData_SA.Rdata")
