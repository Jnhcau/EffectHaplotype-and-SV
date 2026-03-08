library(vcfR)
library(dplyr)
library(data.table)
ggg <- fread("CDS_transcript_merged.bed", header = FALSE)
setnames(ggg, c("seqid", "start", "end", "gene"))

vcf <- read.vcfR("geno.vcf.gz")  
CHR <- getCHROM(vcf)
POS <- getPOS(vcf)
REF <- getREF(vcf)
ALT <- getALT(vcf)
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

filter_genes <- function(CHR, POS, REF, ALT, gt_matrix, ggg) {
  gene_snp_list <- list()
  sample_names <- colnames(gt_matrix)
  

  if(!is.matrix(gt_matrix)) {
    gt_matrix <- as.matrix(gt_matrix)
  }
  
  for (gene_row in seq_len(nrow(ggg))) {
    gene <- ggg$gene[gene_row]
    seqid <- ggg$seqid[gene_row]
    start_pos <- ggg$start[gene_row]
    end_pos <- ggg$end[gene_row]
    

    snp_mask <- CHR == seqid & POS >= start_pos & POS <= end_pos
    snps_in_gene <- which(snp_mask)
    

    if (length(snps_in_gene) == 0) {
      message("No SNPs found for gene: ", gene)
      next
    }
    

    gene_gt <- tryCatch({
      gt_matrix[snps_in_gene, , drop = FALSE]
    }, error = function(e) {
      message("Error accessing SNPs for ", gene, ": ", e$message)
      matrix(NA, nrow = length(snps_in_gene), ncol = ncol(gt_matrix))
    })

    if (nrow(gene_gt) == 0 || all(is.na(gene_gt))) {
      message("No valid SNPs for gene: ", gene)
      next
    }
    
    gene_ref <- REF[snp_mask]
    gene_alt <- ALT[snp_mask]
    
    haplotypes <- apply(gene_gt, 2, function(geno_col) {
      haplo <- character(length(gene_ref))
      valid <- TRUE
      
      for (i in seq_along(gene_ref)) {

        if (is.na(geno_col[i])) {
          valid <- FALSE
          break
        }
        
        alleles <- unlist(strsplit(geno_col[i], "[/|]"))
        if (length(alleles) != 2) {
          valid <- FALSE
          break
        }
        
        if (alleles[1] != alleles[2]) {
          valid <- FALSE
          break
        }
        
        allele_code <- as.integer(alleles[1])
        if (allele_code == 0) {
          haplo[i] <- gene_ref[i]
        } else {
          alts <- unlist(strsplit(gene_alt[i], ","))
          if (allele_code > length(alts)) {
            valid <- FALSE
            break
          }
          haplo[i] <- alts[allele_code]
        }
      }
      
      if (valid) paste(haplo, collapse = "") else NA
    })
    
    valid_idx <- !is.na(haplotypes)
    if (any(valid_idx)) {
      gene_snp_list[[gene]] <- haplotypes[valid_idx]
    } else {
      message("No valid haplotypes for gene: ", gene)
    }
  }
  
  return(gene_snp_list)
}

gene_snp_list <- filter_genes(CHR, POS, REF, ALT, gt_matrix, ggg)

calculate_haplotype_frequencies <- function(gene_snp_list) {
  haplotype_freq_list <- list()  
  
  for (gene in names(gene_snp_list)) {
    haplotypes <- gene_snp_list[[gene]] 
    
    haplotype_table <- table(haplotypes) 
    haplotype_freq_list[[gene]] <- haplotype_table 
  }
  
  return(haplotype_freq_list)
}

haplotype_frequencies <- calculate_haplotype_frequencies(gene_snp_list)

haplotypes_by_gene <- gene_snp_list

save(haplotypes_by_gene, haplotype_frequencies, 
     file = "hap.RData")
