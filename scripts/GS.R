
# Load libraries
library(rrBLUP)
library(data.table)

# --------------------------
# Step 1: Load Genotype & Phenotype
# --------------------------
geno_raw <- fread("G.csv")
samples <- geno_raw$Individual
geno <- as.matrix(sapply(geno_raw[, -1, with = FALSE], as.numeric))
#geno <- as.matrix(geno_raw[, 7:ncol(geno_raw), with = FALSE])
rownames(geno) <- samples

phe <- fread("pheno.csv")  

sample_names <- phe[[1]]
phe_mat <- as.matrix(phe[, -1, with = FALSE])
rownames(phe_mat) <- sample_names

# Keep common samples between geno and phe
common_ids <- intersect(rownames(geno), rownames(phe_mat))
geno <- geno[common_ids, , drop = FALSE]
phe_mat <- phe_mat[common_ids, , drop = FALSE]

# --------------------------
# Step 2: Define k-fold CV function
# --------------------------
library(caret)

run_repeated_kfold <- function(pheno_vec, geno_mat, k = 10, n_repeats = 10, seed = 2025) {
  all_acc <- c()
  
  for (r in 1:n_repeats) {
    set.seed(seed + r)
    folds <- createFolds(pheno_vec, k = k, list = TRUE, returnTrain = FALSE)
    
    for (i in seq_along(folds)) {
      test_idx <- folds[[i]]
      train_idx <- setdiff(seq_along(pheno_vec), test_idx)
      
      geno_train <- geno_mat[train_idx, , drop = FALSE]
      geno_test  <- geno_mat[test_idx, , drop = FALSE]
      pheno_train <- pheno_vec[train_idx]
      pheno_test  <- pheno_vec[test_idx]
      
      model <- mixed.solve(y = pheno_train, Z = geno_train)
      u <- as.matrix(model$u)
      pred_test <- as.numeric(geno_test %*% u) + as.numeric(model$beta)
      
      all_acc <- c(all_acc, cor(pred_test, pheno_test, use = "complete.obs"))
    }
  }
  
  list(
    Mean = mean(all_acc, na.rm = TRUE),
    SD = sd(all_acc, na.rm = TRUE),
    Accuracy = all_acc
  )
}

results <- lapply(seq_len(ncol(phe_mat)), function(i) {
  run_repeated_kfold(phe_mat[, i], geno, k = 10, n_repeats = 10, seed = 2025)
})
# --------------------------
# Step 3: save result
# --------------------------
acc_long <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    Trait    = colnames(phe_mat)[i],
    Accuracy = results[[i]]$Accuracy,
    Mean     = results[[i]]$Mean,
    SD       = results[[i]]$SD
  )
}))
write.csv(acc_long,
          file = "GS_accuracy_long_with_mean_sd.csv",
          row.names = FALSE)
