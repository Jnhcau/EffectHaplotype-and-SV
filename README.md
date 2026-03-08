# Effect Haplotypes and Structural Variants in Modern Maize Breeding
![Main Figure](https://github.com/Jnhcau/EffectHaplotype-and-SV/blob/main/image.jpg)
This repository contains the main scripts used in our paper:  
“Effect haplotypes and structural variants reveal genetic architecture of modern maize breeding”.

---

# Repository Overview

## 1. vcf2EH.R
A script for constructing Effect Haplotypes (EH) from VCF files.  
With this script, you can easily generate EH datasets from your own genotype data.

---

## 2. EHGWAS.R
This script performs GWAS based on EHs, following the method described in our paper.

---

## 3. GS.R
A script for genomic selection (GS) based on SNPs, SVs, or EHs.

---

## 4. Sankey-diagram.R
This is the script I used to draw Sankey diagrams.  
It produces clean and visually appealing figures that can be used to show things like species migration or changes in allele frequencies across populations.  

Feel free to use it in your own studies or papers.

---

## 5. SV-TEanno.sh
This script shows how I analyzed the TE composition of structural variants (SVs).

In short:

- extract SV sequences  
- annotate them using RepeatMasker  
- summarize TE composition using the EDTA pipeline
