library(lme4qtl)
library(data.table)
library(dplyr)
library(lme4)
library(parallel)
library(doSNOW)
library(RLRsim)

G <- readRDS("hap.rds")
hap_cols <- colnames(G)[-1]
gene_id <- sub("_[0-9]+$", "", hap_cols)

hap_per_gene <- table(gene_id)
genes_keep <- names(hap_per_gene[hap_per_gene >= 2])

cols_gene_keep <- hap_cols[gene_id %in% genes_keep]

hap_matrix <- as.matrix(G[, cols_gene_keep])
hap_count  <- colSums(hap_matrix, na.rm = TRUE)
cols_freq_keep <- cols_gene_keep[hap_count >= 2]

G_filtered <- G[, c("Individual", cols_freq_keep)]
hap_matrix <- as.matrix(G_filtered[, -1])

valid_genes <- unique(sub("_[0-9]+$", "", cols_freq_keep))
##############
# 1.Xc
##############
K <- hap_matrix %*% t(hap_matrix)
K <- K * nrow(G_filtered) / sum(diag(K))
rownames(K) <- colnames(K) <- G_filtered$Individual

hap_scaled <- scale(hap_matrix, center = TRUE, scale = FALSE)
pca_res <- prcomp(hap_scaled, center = FALSE)

df_cov <- as.data.frame(pca_res$x[, 1:3])
colnames(df_cov) <- paste0("PC", 1:3)
df_cov$id <- rownames(K)
##############
# 2. EHfixedeffect
##############
pheno <- read.csv("pheno.csv")
pheno <- pheno[match(G_filtered$Individual, pheno$sample), ]

#表型
y_raw <- pheno$EH

y_raw[is.na(y_raw)] <- mean(y_raw, na.rm = TRUE)
df_null <- data.frame(
  y = y_raw,
  df_cov,
  id = G_filtered$Individual
)

EHfixed <- relmatLmer(
  y ~ PC1 + PC2 + PC3 +(1 | id),
  data = df_null,
  relmat = list(id = K),
  REML = FALSE
)

KinEff <- ranef(EHfixed)$id[,1]
y <- y_raw - KinEff
y_vec <- setNames(y, G_filtered$Individual)
##############
# 3. 
##############
cl <- makeCluster(5)      
registerDoSNOW(cl)

results_df <- foreach(
  g = valid_genes,
  .combine = rbind,
  .packages = c("data.table", "lme4")
) %dopar% {
  
  cols <- grep(paste0("^", g, "_"), colnames(hap_matrix))
  if (length(cols) < 2) return(NULL)
  
  dt <- as.data.table(hap_matrix[, cols, drop = FALSE])
  dt$id <- G_filtered$Individual
  
  dt_long <- melt(
    dt,
    id.vars = "id",
    variable.name = "hap_id",
    value.name = "hap_val"
  )[hap_val == 1]
  
  if (nrow(dt_long) < 5) return(NULL)
  
  dt_long <- merge(dt_long, df_cov, by = "id")
  dt_long$y <- y_vec[dt_long$id]
  
  hap_tab <- sort(table(dt_long$hap_id))
  FHlevels <- length(hap_tab)
  FHnumber <- paste(hap_tab, collapse = "/")
  
  fit0 <- lm(y ~ PC1 + PC2 + PC3 , data = dt_long)
  fit1 <- lmer(y ~ PC1 + PC2 + PC3 + (1 | hap_id),
               data = dt_long, REML = FALSE)
  
  vc <- as.data.frame(VarCorr(fit1))
  h2 <- vc$vcov[vc$grp == "hap_id"] / sum(vc$vcov)
  
  lrt <- anova(fit1, fit0)
  
  data.frame(
    pheName = "EH",
    Gene = g,
    h2 = h2,
    FHlevels = FHlevels,
    FHnumber = FHnumber,
    AIC = AIC(fit1),
    BIC = BIC(fit1),
    logLik = logLik(fit1),
    Chisq = lrt$Chisq[2],
    Pvalue = lrt$`Pr(>Chisq)`[2],
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)

write.csv(results_df, 
          file = "GWAS_EHBLUP_PCA3.csv",   
          row.names = FALSE,         
          quote = FALSE)            


bonf_threshold <- 0.05 /  nrow(results_df)
#矫正基因组膨胀
mean.chisq = mean(results_df$Chisq,na.rm=T)
gc.chisq = mean.chisq/0.456
print(gc.chisq)

