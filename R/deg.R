# %%
library(DESeq2)
library(configr)

# %%
HOME = "/home/ubuntu"
source(file.path(HOME, "dev/rnaseqr/R/functions.R"))
deseq_deg_multi

# %%
config = read.config(file = "/home/ubuntu/dev/rnaseqr/config.toml")
config

config
# %%
FIGS <- config$output$figures
project <- config$output$project

counts <- as.matrix(read.table(config$counts$path, 
                               sep = config$counts$delim,
                               header = TRUE,
                               row.names = 1))
counts <- round(counts)

write.table(counts, file = file.path(project, "rna_raw_rounded_counts.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# %%
meta = read_delim(config$metadata$path, 
                  delim = config$metadata$delim)

########################################################
# %% Remove any genes/samples that have 0 counts in all columns/rows.
########################################################
filt_counts <- rm_zero_axes(counts)

########################################################
# %% Remove singleton rows
########################################################
filt_counts <- rm_singleton_genes(filt_counts)

########################################################
# %% Remove genes with duplicate gene symbols
#    And then convert the count matrix to symbols
########################################################
if (!any(grepl("ENSG[0-9]+[.][0-9]+", rownames(filt_counts))))
  stop("Expecting Ensembl gene ids, got something else.")
rownames(filt_counts) <- gsub("(.+)[.].*$", "\\1", rownames(filt_counts))
id2sym = kallisto_t2g_id2sym(t2g_path = "/home/ubuntu/shan_rna_atac/ref/kallisto/homo_sapiens/transcripts_to_genes.txt") 
filt_counts <- remove_genes_with_dup_symbols(filt_counts, id2sym)
rownames(filt_counts) <- id2sym[rownames(filt_counts), ]$symbol

########################################################
# %% Define the minimum count value to keep, based on config
#    And remove genes whose rowSums are that are <that
########################################################
min_count_value = unname(quantile(apply(filt_counts, 1, sum), config$counts$filter$remove_lowest))
keep_genes = names(which(rowSums(filt_counts)>=min_count_value))
filt_counts = filt_counts[keep_genes,]

min_count_value
########################################################
# %% Write out the filtered count matrix
########################################################
write.table(filt_counts, file = file.path(project, "filtered_filt_counts.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

########################################################
# %% Generate DESeq object and filter by group size
########################################################
stopifnot(colnames(filt_counts) == meta$sample)
dds <- DESeqDataSetFromMatrix(countData = filt_counts, colData = meta, design = formula(config$deg$design))
smallestGroupSize <- config$counts$filter$smallest_group_size
keep <- rowSums(counts(dds) >= min_count_value) >= smallestGroupSize
dds <- dds[keep, ]

########################################################
# %% Set up dependant/independant variables
########################################################
indep_vars = strsplit(gsub("^[~]", "", config$deg$design), split = "[+]")[[1]]

########################################################
# %% Run PCA
########################################################
pca_pltdata = DESeq2::plotPCA(vst(dds), 
                              intgroup = c(config$deg$sample_column, indep_vars),
                              ntop = 5000, 
                              returnData = TRUE) 

########################################################
# %% Run DESeq
########################################################
dds = DESeq(dds)

########################################################
# %% Do comparisons and write out results
########################################################
source(file.path(HOME, "dev/rnaseqr/R/functions.R"))
config = read.config(file = "/home/ubuntu/dev/rnaseqr/config.toml")
results = list()
indepvar = indep_vars[1]
meta[,indepvar] = as.factor(pull(meta, !!sym(indepvar)))
combinations = t(combn(levels(meta[[indepvar]]), 2))
results[[indepvar]] = deseq_deg_multi(deseq = dds, 
                          combinations = combinations,
                          indepvar = indepvar,
                          lfccut = config$deg$log2_fc_cutoff,
                          alpha = config$deg$alpha,
                          pval_corr_method = config$deg$pval_corr_method,
                          write_results_to = file.path(project, glue::glue("rna_{indepvar}_deg.csv")))
########################################################
# %% Save RData from analysis
########################################################
save(config, counts, meta, filt_counts, pca_pltdata, dds, results, file = file.path(project, glue::glue("deg.RData")))


# # %%                  w
# combinations = t(combn(levels(meta$treatment),2))
# results = deseq_deg_multi(deseq = dds, combinations = combinations, indepvar = "treatment", write_results_to = "/srv/http/shan/paired/rna_treatment_deg.csv")
#
# # %%
# results = read_delim("/srv/http/shan/paired/rna_treatment_deg.csv", delim = ',') %>%
#   mutate(is_sig = ifelse(padj <= 0.05 & abs(log2FoldChange) > 0.1, TRUE, FALSE)) %>%
#   filter(!is.na(padj))
# deg_gene_list = results %>%
#   filter(is_sig) %>%
#   pull(gene)
#
# deg_gene_list
#
# colnames(results)
# # %%
# which(is.na(results$padj))
# plot_volcano(results,
#              group1 = "BHB",
#              group2 = "control",
#              output_path = file.path(FIGS, "rna_treatment_volcano_lfccut0.1.pdf"),
#              n_label = 20,
#              lfccut = 0.1)
#
# # %%
# length(which(results$is_sig))
# length(which(results$is_sig & results$log2FoldChange > 0))
# length(which(results$is_sig & results$log2FoldChange < 0))
#
# colnames(results)
#
# deg_gene_list
#
# # %%
# list.files(FIGS)
# daTab = as_tibble(read.table(file.path(FIGS, "atac_diffAcc_coding_promoter_4000up_results.tsv"), sep = "\t", header = TRUE)) %>%
#   mutate(fc_numerator = "BHB", fc_denominator = "control") %>%
#   mutate(padj = p.adjust(pvalue, "BH")) %>%
#   mutate(is_sig = ifelse(padj <= 0.05 & abs(log2FoldChange) > lfccut, TRUE, FALSE)) %>%
#   # rename(gene = gene_name) %>%
#   filter(!is.na(padj))
#
# da_promoter_genes = daTab %>%
#   filter(is_sig) %>%
#   pull(gene)
#
# colnames(daTab)
# # %% 
# da_promoter_genes
# deg_gene_list
# which(deg_gene_list %in% da_promoter_genes)
# ########################################################
# # %% Heatmap, DEG
# ########################################################
# degenes = results[results$is_sig,]$gene
# pdf(file.path(FIGS, glue::glue("rna_deg_treatment_lfccut0.1_hm.pdf")), width = 8.5, height = 5.5)
# hm_counts(counts, genes=degenes, show_row_names = TRUE, column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 8))
# dev.off()
#
# # %%

# %%

