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
hm_counts = function(counts,
                     genes = NULL,
                     scale = TRUE,
                     col = circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red")), 
                     ...) {
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
                               col = col,
                               # width = width,
                               # height = height,
                               ...
  )
  hm
}

# %%

hm_size_and_position = function(figure_width,
                                figure_height,
                                row_dend_width = unit(0.2, "in"),
                                row_names_max_width = unit(1, "in"),
                                column_dend_height = unit(0.2, "in"),
                                column_names_max_height = unit(1, "in"),
                                hm_width_scaling = 0.8
                                ) {
  LEFT_PADDING = unit(0.05, "in")
  values = c()
  values$row_dend_width = row_dend_width
  values$row_names_max_width = row_names_max_width
  values$column_dend_height = column_dend_height
  values$column_names_max_height = column_names_max_height
  values$width = (figure_width - row_dend_width - row_names_max_width)*hm_width_scaling
  values$height = (figure_height - column_dend_height - column_names_max_height)
  values$padding = unit(c(0, LEFT_PADDING, 0, (figure_width-(values$width + row_dend_width + LEFT_PADDING))), "in")
  values

}
# %%
pretty_draw_hm = function(counts,
                          genes = NULL,
                          scale = TRUE,
                          figure_width = unit(6.5, "in"),
                          figure_height = unit(4.5, "in"),
                          col = circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red")),
                          show_row_names = FALSE,
                          ...) {
  legend_title = ifelse(scale, "Scaled expression", "Expression")
  legend = Legend(col_fun = col, title = legend_title, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6))
  hm_plt_params = hm_size_and_position(figure_width = unit(7.5, "in"),
                                       figure_height = unit(5.0, "in"),
                                       row_names_max_width = unit(0, "in"))
  heatmap = hm_counts(counts = counts,
                      genes = genes,
                      scale = scale,
                      col = col,
                      show_row_names = show_row_names,
                      show_heatmap_legend = FALSE,
                      row_dend_width = hm_plt_params$row_dend_width,
                      row_names_max_width = hm_plt_params$row_names_max_width,
                      column_dend_height = hm_plt_params$column_dend_height,
                      column_names_max_height = hm_plt_params$column_names_max_height,
                      width = hm_plt_params$width,
                      height = hm_plt_params$height,
                      ...)
  draw(heatmap,
       padding = hm_plt_params$padding,
  )
  draw(legend, x = unit(0.90, "npc"), y = unit(0.5, "npc"))
}

# %%
top_n_variable_genes = function(counts, ntopvar) {
  names(sort(apply(counts, 1, var), decreasing = TRUE)[1:ntopvar])
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
  message("Removing ", length(singleton_idx), " genes with singleton counts, leaving ", round(sum(filt_counts)/sum(counts)*100, 4), "% of total counts.")
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
  as_tibble(allres)
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
