# %%
library(DESeq2)

source("/home/ubuntu/dev/ghassemi_bulkrna/scripts/functions.R")

# %%
FIGS <- "/srv/http/shan/paired"
project <- "/home/ubuntu/shan_rna_atac/analysis/rna"

counts <- as.matrix(read.table(file.path(project, "outs", "rna_ensembl_counts_matrix_n21695x8.tsv"), sep = "\t", header = TRUE, row.names = 1))
counts <- round(counts)


colnames(counts)
meta <- cbind(colnames(counts), gsub("[0-9]$", "", colnames(counts)), gsub("[a-zA-Z]+", "", colnames(counts)))
colnames(meta) <- c("sample", "treatment", "biorep")
meta <- as.data.frame(meta)
meta$treatment = as.factor(meta$treatment)

########################################################
# %% Remove any genes/samples that have 0 counts in all columns/rows.
########################################################
counts <- rm_zero_axes(counts)

########################################################
# %% Heatmap raw counts, all
########################################################
pdf(file.path(FIGS, "rna_raw_counts_all.pdf"), width = 8.5, height = 200)
hm_counts(counts, row_names_gp = gpar(fontsize = 0.5), show_row_names = TRUE)
dev.off()

########################################################
# %% Remove singleton rows
########################################################
counts <- rm_singleton_genes(counts)

########################################################
# %% Remove genes with duplicate gene symbols
#    And then convert the count matrix to symbols
########################################################
rownames(counts) <- gsub("(.+)[.].*$", "\\1", rownames(counts))
id2sym = kallisto_t2g_id2sym(t2g_path = "/home/ubuntu/shan_rna_atac/ref/kallisto/homo_sapiens/transcripts_to_genes.txt") 
counts <- remove_genes_with_dup_symbols(counts, id2sym)
rownames(counts) <- id2sym[rownames(counts), ]$symbol

########################################################
# %% Define the count value that defines the lowest 20%
########################################################
q20 = unname(quantile(apply(counts, 1, sum), .2))

########################################################
# %% Heatmap raw counts, minus singletons
########################################################
pdf(file.path(FIGS, glue::glue("rna_raw_counts_rmSingletons.pdf")), width = 8.5, height = 11)
hm_counts(counts)
dev.off()

########################################################
# %% Heatmap, rowSums >q20
########################################################
q20 = unname(quantile(apply(counts, 1, sum), .2))
pdf(file.path(FIGS, glue::glue("rna_raw_counts_rowSumsGtQ20.pdf")), width = 8.5, height = 11)
hm_counts(counts, genes=names(which(rowSums(counts)>=q20)), column_names_gp = gpar(fontsize = 12))
dev.off()

########################################################
# %% Heatmap, rowSums >q20, cpm
########################################################
pdf(file.path(FIGS, glue::glue("rna_raw_counts_rowSumsGtQ20_cpm.pdf")), width = 8.5, height = 11)
sample_norm_counts = apply(counts, 2, function(x) x/sum(counts)*1e6)
hm_counts(sample_norm_counts, genes=names(which(rowSums(counts)>=q20)), column_names_gp = gpar(fontsize = 12))
dev.off()

########################################################
# %% Heatmap top var genes
########################################################
ntopvar = 5000
topvar = names(sort(apply(counts, 1, var), decreasing = TRUE)[1:ntopvar])
pdf(file.path(FIGS, glue::glue("rna_raw_counts_top{ntopvar}var.pdf")), width = 8.5, height = 11)
hm_counts(counts, genes=topvar, column_names_gp = gpar(fontsize = 12))
dev.off()

########################################################
# %% Generate DESeq object and filter
########################################################
stopifnot(colnames(counts) == meta$sample)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~treatment)
smallestGroupSize <- min(table(meta$treatment))
keep <- rowSums(counts(dds) >= q20) >= smallestGroupSize
dds <- dds[keep, ]

########################################################
# %% Make PCA
########################################################
pdf(file.path(FIGS, "raw_counts_pca.pdf"), width = 13.33, height = 7.5)
# DESeq2::plotPCA(vst(dds), intgroup = "biorep", ntop = 2000) + geom_text(aes(label = sample))
pltdata = DESeq2::plotPCA(vst(dds), intgroup = c("sample", "treatment", "biorep"), ntop = 5000, returnData = TRUE) 
ggplot(pltdata, aes(x=PC1, y = PC2, color = treatment)) +
  geom_point(size = 10) +
  geom_text(aes(label = sample), color = "black", size = 6) +
  theme_bw() +
  theme(
        panel.grid = element_blank()
        )
dev.off()



########################################################
# %% Run DESeq
########################################################
dds = DESeq(dds)

# %%                  w
combinations = t(combn(levels(meta$treatment),2))
results = deseq_deg_multi(deseq = dds, combinations = combinations, indepvar = "treatment", write_results_to = "/srv/http/shan/paired/rna_treatment_deg.csv")

# %%
results = read_delim("/srv/http/shan/paired/rna_treatment_deg.csv", delim = ',') %>%
  mutate(is_sig = ifelse(padj <= 0.05 & abs(log2FoldChange) > 0.1, TRUE, FALSE)) %>%
  filter(!is.na(padj))
deg_gene_list = results %>%
  filter(is_sig) %>%
  pull(gene)

deg_gene_list

colnames(results)
# %%
which(is.na(results$padj))
plot_volcano(results,
             group1 = "BHB",
             group2 = "control",
             output_path = file.path(FIGS, "rna_treatment_volcano_lfccut0.1.pdf"),
             n_label = 20,
             lfccut = 0.1)

# %%
length(which(results$is_sig))
length(which(results$is_sig & results$log2FoldChange > 0))
length(which(results$is_sig & results$log2FoldChange < 0))

colnames(results)

deg_gene_list

# %%
list.files(FIGS)
daTab = as_tibble(read.table(file.path(FIGS, "atac_diffAcc_coding_promoter_4000up_results.tsv"), sep = "\t", header = TRUE)) %>%
  mutate(fc_numerator = "BHB", fc_denominator = "control") %>%
  mutate(padj = p.adjust(pvalue, "BH")) %>%
  mutate(is_sig = ifelse(padj <= 0.05 & abs(log2FoldChange) > lfccut, TRUE, FALSE)) %>%
  # rename(gene = gene_name) %>%
  filter(!is.na(padj))

da_promoter_genes = daTab %>%
  filter(is_sig) %>%
  pull(gene)

colnames(daTab)
# %% 
da_promoter_genes
deg_gene_list
which(deg_gene_list %in% da_promoter_genes)
########################################################
# %% Heatmap, DEG
########################################################
degenes = results[results$is_sig,]$gene
pdf(file.path(FIGS, glue::glue("rna_deg_treatment_lfccut0.1_hm.pdf")), width = 8.5, height = 5.5)
hm_counts(counts, genes=degenes, show_row_names = TRUE, column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 8))
dev.off()

# %%

