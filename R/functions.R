# %%
library(DESeq2)
library(patchwork)
library(tidyverse)
library(glue)
library(readxl)
library(tximport)
library(tximeta)
library(gtable)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(tidytable)
library(circlize)

# %%
hm_counts = function(counts, genes = NULL, scale = TRUE, show_row_names = FALSE, show_column_names = TRUE, ...) {
  if (!is.null(genes)) {
    counts = counts[genes,]
  }
  if (scale) {
    counts_scale = t(apply(counts, 1, scale))
    colnames(counts_scale) = colnames(counts)
    rownames(counts_scale) = rownames(counts)
    counts = counts_scale
  }

  hm = ComplexHeatmap::Heatmap(counts,
          show_row_names = show_row_names,
          show_column_names = show_column_names,
          ...
  )
  hm
}

# %%
kallisto_t2g_id2sym = function(t2g_path) {
  id2sym <- unique(read.table(t2g_path, sep = "\t")[, c(2, 3)])
  colnames(id2sym) <- c("ensemblid", "symbol")
  id2sym$ensemblid <- gsub("(.+)[.].*$", "\\1", id2sym$ensemblid)
  rownames(id2sym) <- id2sym$ensemblid
  id2sym
}
# %%
# %%
rm_zero_axes = function(counts) {
  all_zero_genes = which(apply(counts, 1, function(x) sum(x)) == 0)
  all_zero_samples = which(apply(counts, 2, function(x) sum(x)) == 0)
  message("There are ", length(all_zero_genes), " genes that have 0 counts in all samples.")
  message("There are ", length(all_zero_samples), " samples that have 0 counts in all genes")
  if (length(all_zero_genes) > 0) {
    counts = counts[-all_zero_genes,]
  }
  if (length(all_zero_samples) > 0) {
    counts = counts[,-all_zero_samples]
  }
  counts
}
# %%
rm_singleton_genes = function(counts) {
  nsamp = ncol(counts)
  singleton_idx = unname(which(
                        apply(counts,1,function(row) {
                                       nsamp = length(row)
                                       length(which(row==0)) == (nsamp-1)
    })))
  filt_counts = counts[-singleton_idx,]
  message("Removing ", length(singleton_idx), " genes with singleton counts, leaving ", round(sum(filt_counts)/sum(counts)*100, 2), "% of counts.")
  filt_counts
}

# %%
remove_genes_with_dup_symbols = function(counts, id2sym) {
  to_remove = as_tibble(id2sym[rownames(counts),]) %>%
    group_by(symbol) %>%
    summarise(ensemblid = ensemblid, n=n()) %>%
    filter(n>1) %>%
    pull(ensemblid)
  filt_counts = counts[!rownames(counts) %in% to_remove,]
  message("Removing ", length(to_remove), " genes with duplicated gene ids, leaving ", round(sum(filt_counts)/sum(counts)*100, 2), "% of counts.")
  filt_counts
}

# %%
deseq_deg_combinatorial= function(deseq,
                              indepvar,
                              metadata,
                              lfccut = 0.5,
                              alpha = 0.05,
                              pval_corr_method = 'fdr',
                              write_results_to = NULL
                              ) {
  combs = t(combn(unique(metadata[,indepvar]), m=2))
  allres = deseq_deg_multi(
                           deseq,
                           indepvar,
                           combinations = combs,
                           lfccut = lfccut,
                           alpha = alpha,
                           pval_corr_method = pval_corr_method,
                           write_results_to = write_results_to
  )
  allres
}

# %%

deseq_deg_multi = function(
                           deseq,
                           indepvar,
                           combinations,
                           lfccut = 0.5,
                           alpha = 0.05,
                           pval_corr_method = 'fdr',
                           write_results_to = NULL
                           ) {

  allres = data.table()
  for (row in seq_along(combinations[,1])) {
    n = combinations[row,1]
    d = combinations[row,2]
    resdf = results(deseq, contrast=c(indepvar, n, d), alpha=alpha, independentFiltering = FALSE, cooksCutoff = TRUE)
    resdf = as.data.frame(resdf)
    resdf$gene = rownames(resdf)
    resdf$fc_numerator = n
    resdf$fc_denominator = d
    resdf = subset(resdf, select = -c(padj))
    #resdf = resdf[!is.na(resdf$pvalue),]
    allres = bind_rows(allres, resdf)
  }
  allres$padj = p.adjust(allres$pvalue, method = pval_corr_method)
  allres$is_sig = ifelse(allres$padj <= alpha & abs(allres$log2FoldChange) > lfccut, TRUE, FALSE)

  allres = allres %>%
    dplyr::select(gene, fc_numerator, fc_denominator, log2FoldChange, padj, pvalue, everything()) %>%
    arrange(desc(is_sig), log2FoldChange)

  if (!is.null(write_results_to)) {
    if (grepl("[.]csv$", write_results_to)) {
      delim = ','
    } else if (grepl("[.]tsv$", write_results_to)) {
      delim = '\t'
    } else {
      stop("Unknown file extension passed to deseq_deg_multi Use `.csv` or `.tsv`.")
    }
    write_delim(allres, file = write_results_to, delim =  delim, quote = "none")
  }
  allres
}
 
# %%
plot_volcanoes_combinatorial = function(results,
                                        indepvar,
                                        metadata,
                                        output_path= NULL,
                                        n_label = 10,
                                        lfccut = 0.5,
                                        alpha = 0.05,
                                        width = 4.25,
                                        height = 4.25,
                                        dotsize = 0.2,
                                        labelsize = 2,
                                        return_figs = FALSE
                                        ) {
  combs = t(combn(unique(metadata[,indepvar]), m=2))
  plts = plot_volcanoes_multi(results = results,
                       combinations = combs,
                       output_path = output_path,
                       n_label = n_label,
                       lfccut = lfccut,
                       alpha = alpha,
                       width = width,
                       height = height,
                       dotsize = dotsize,
                       labelsize = labelsize,
                       return_figs = return_figs
                       )
  if (return_figs) {
    return(plts)
  } else {
    return (NULL)
  }
}

# %%
plot_volcanoes_multi = function(results,
                                combinations,
                                output_path = NULL,
                                n_label = 10,
                                lfccut = 0.5,
                                alpha = 0.05,
                                width = 4.25,
                                height = 4.25,
                                dotsize = 0.2,
                                labelsize = 2,
                                return_figs = FALSE
                                ) {

  plts = list()
  for (row in seq_len(nrow(combinations))) {
    n = combinations[row,1]
    d = combinations[row,2]
    sub = subset(results, fc_numerator == n & fc_denominator == d)
    plts[[row]] = plot_volcano(reults = sub,
                               output_path = output_path,
                               n_label = n_label,
                               lfccut = lfccut,
                               alpha = alpha,
                               width = width,
                               height = height,
                               dotsize = dotsize,
                               labelsize = labelsize,
                               return_figs = return_figs
    )
  }
  if (return_figs) {
    return(plts)
  } else {
    pdf(output_path, width = width, height = height)
    for (p in plts) {
      print(p)
    }
    dev.off()
  }
}

plot_volcano = function(
                        results,
                        group1,
                        group2,
                        output_path = NULL,
                        n_label = 10,
                        lfccut = 0.5,
                        alpha = 0.05,
                        width = 4.25,
                        height = 4.25,
                        dotsize = 0.2,
                        labelsize = 2,
                        return_figs = FALSE,
                        ylab = "-log10(padj)",
                        xlab = "log2FoldChange",
                        draw_cutoff_lines = TRUE,
                        y_transform_func = function(y) -1*log10(y),
                        ggrepel_max_overlaps = 15
                        ) {
  title = glue("{group1} vs {group2}")
  results = results %>%
    filter(fc_numerator %in% c(group1, group2) & fc_denominator %in% c(group1, group2))
  labels_to_show = results %>% 
    filter(is_sig) %>% 
    slice_min(padj, n=n_label) %>% 
    pull(gene)
  results = results %>% 
    mutate(display_label = ifelse(gene %in% labels_to_show, TRUE, FALSE))
  plt = ggplot(data = results,  aes(x = log2FoldChange, y = y_transform_func(padj), color = is_sig)) +
    geom_point(size = dotsize) +
    geom_text_repel(data = subset(results, display_label), aes(label = gene), color='black', size = labelsize, max.overlaps = ggrepel_max_overlaps) +
    scale_color_manual(values = c("black", "red")) +
    ylab(ylab) +
    xlab(xlab) +
    guides(color = "none") +
    ggtitle(title) +
    theme_classic()
  if (draw_cutoff_lines) {
    plt = plt +
      geom_hline(yintercept = -log10(alpha), linetype = 'dashed') +
      geom_vline(xintercept = -lfccut, linetype = 'dashed') +
      geom_vline(xintercept = lfccut, linetype = 'dashed')
  }
  if (return_figs) {
    return(plt)
  } else {
    pdf(output_path, width = width, height = height)
    print(plt)
    dev.off()
  }
}
