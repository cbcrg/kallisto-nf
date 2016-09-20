#!/usr/bin/env Rscript
library("sleuth")

args <- commandArgs(TRUE)

sample_id <- dir(args[1])
sample_id

kal_dirs <- sapply(sample_id, function(id) file.path(args[1], id))
kal_dirs

s2c <- read.table(args[2], header = TRUE, stringsAsFactors=FALSE) 
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c <- s2c[order(s2c$condition), ]

print(s2c)

ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl", version=79)

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = ensembl)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Load the kallisto processed data into the object
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)

# Estimate parameters for the sleuth response error measurement (full) model 
so <- sleuth_fit(so)

# Perform test
so <- sleuth_wt(so,  'conditionscramble')

gene_table <- sleuth_gene_table(so, test = "conditionscramble", test_type = "wt")

write.table(gene_table, paste("gene_table_results.txt"), sep="\t")

save(so, file=paste("sleuth_object.so"))
