######## Load Libraries #########
suppressPackageStartupMessages({
  library(BiocParallel)
  library(biomaRt)
  library(CellChat)
  library(cli)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(cowplot)
  library(data.table)
  library(DESeq2)
  library(dplyr)
  library(DropletUtils)
  library(dunn.test)
  library(easyGgplot2)
  library(edgeR)
  library(EnhancedVolcano)
  library(enrichplot)
  library(ggbeeswarm)
  library(ggnewscale)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(gridExtra)
  library(harmony)
  library(here)
  library(ktplots)
  library(limma)
  library(magrittr)
  library(Matrix.utils)
  library(monocle3)
  library(nebula)
  library(Nebulosa)
  library(org.Hs.eg.db)
  library(patchwork)
  library(pheatmap)
  library(purrr)
  library(RColorBrewer)
  library(readxl)
  library(scCustomize)
  library(scDblFinder)
  library(SCP)
  library(Seurat) # should be v 4.4.0
  library(SingleCellExperiment)
  library(slingshot)
  library(speckle)
  library(stringr) 
  library(tibble)
  library(tidymodels)
  library(tidyr)
  library(tidyverse)
  library(tradeSeq)
  library(UCell)
  library(viridis)
})

### Configurations
CFG <- list()
CFG$n_cores <- 128
CFG$random_seed <- 12345   #Random seed for UMAP generation
CFG$mypal<-c("#F28C8F","#ECCF76", '#C3894C',"#6ABEAE", "#F2B8A5", "#147F84", '#C69D9C',  "#E86562", "#B48DAC", "#7CC2C8", "#856889", "#06434B", "#544469", "#3D716A","#78606D", "#AD6279", '#A3B6CB', "#A5C2B2", "#548FA7", "#59707E",   "#6C947E", "#A99688", "#D88E9D","#DD432F", "#D2AE91", "#6A8C35","#910911")
CFG$groupspal<-c("CTRL" = "#A3BF8C", "MAD" = "#D08770", "RES" = "#808EAF", "SAD" = "#813E3C", "FETAL"="#DAAC50")
CFG$fetalpal<-"#DAAC50"
CFG$gradient_white <- c("lightgrey", "#542437")
CFG$gradient_colors<-c("#2C5F71", "#F4D3B5","#542437")

set.seed(CFG$random_seed)

data.table::setDTthreads(CFG$n_cores)
options(future.globals.maxSize = 50 * 1024 ^ 3)

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 6
RNGversion("3.5.0")

### Function to create mito vs velo plot
fun_mito_velo_plot <- function(sObj, sample_id) {
  metadata<-sObj@meta.data
  ggplot(metadata, aes(
    x = qlogis(percent_mito / 100), 
    y = qlogis(velocyto_unsplicedratio)
  )) +
    geom_bin2d() +
    scale_x_continuous(labels = function(x) round(100 * plogis(x), digits = 1)) +
    scale_y_continuous(labels = function(x) round(100 * plogis(x), digits = 1)) +
    geom_vline(xintercept = qlogis(10 / 100)) + 
    geom_hline(yintercept = qlogis(25 / 100)) + 
    ggtitle(paste0("Sample ", sample_id)) + 
    xlab("pct_mito") + 
    ylab("intronic_ratio") +  
    scale_fill_gradientn(colors = CFG$gradient_colors) +
    theme(
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 20),
      title = element_text(size = 20)
    )
}

### Function to create feature violin plot
fun_nFeatureRNA_plot <- function(sObj, sample_id, qc_thresholds) {
  metadata <- sObj@meta.data
  metadata$nFeature_RNA[metadata$nFeature_RNA == 0] <- 0.1
  
  # Calculate medians for each sample
  medians <- metadata %>%
    group_by(sample) %>%
    summarize(median_nFeature = median(nFeature_RNA), .groups = "drop")
  
  ggplot(metadata, aes(x = sample, y = nFeature_RNA, fill = group)) +
    geom_violin(scale = "width", trim = FALSE) + 
    geom_point(data = medians, aes(x = sample, y = median_nFeature), 
               color = "black", size = 1, shape = 21, fill = "black") +  # Median dots
    scale_y_log10(labels = label_scientific()) +  # Apply log scale
    scale_fill_manual(values = CFG$groupspal) +   # Map colors to groups
    geom_hline(yintercept = qc_thresholds$nFeature_RNA_min[which(qc_thresholds$Sample == sample_id)]) +
    geom_hline(yintercept = qc_thresholds$nFeature_RNA_max[which(qc_thresholds$Sample == sample_id)]) +
    theme_minimal() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.text.y = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      plot.title = element_text(size = 24, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + 
    labs(title = "nFeature_RNA", x = "", y = "", fill = "Group") +
    NoLegend()
}


### Function for QC filtering
QC_filtering <- function(seurat_obj, qc_thresholds) {
  
  filtered_list <- list()
  
  # Loop over each sample in the metadata
  for (sample_name in unique(qc_thresholds$Sample)) {
    # Get the filtering parameters for the current sample
    nFeature_RNA_min <- qc_thresholds$nFeature_RNA_min[which(qc_thresholds$Sample == sample_name)]
    nFeature_RNA_max <- qc_thresholds$nFeature_RNA_max[which(qc_thresholds$Sample == sample_name)]
    pct_mito <- qc_thresholds$pct_mito[which(qc_thresholds$Sample == sample_name)]
    intronic_ratio_min <- qc_thresholds$intronic_ratio_min[which(qc_thresholds$Sample == sample_name)]
    doublet_score <- qc_thresholds$scDblFinder_score[which(qc_thresholds$Sample == sample_name)]
    
    # Informing the user about the filtering criteria
    print(paste0("Filtering Sample ", sample_name, " with nFeature_RNA>", nFeature_RNA_min, 
                 ", nFeature_RNA<", nFeature_RNA_max, ", percent_mito<", pct_mito, 
                 ", Velocyto intronic ratio>", intronic_ratio_min,
                 ", Doublet finder score <", doublet_score))
    
    # Filter the merged Seurat object for the current sample
    sample_filtered <- subset(seurat_obj, cells = WhichCells(seurat_obj, expression = (sample == sample_name &
                                                                                         nFeature_RNA > nFeature_RNA_min &
                                                                                         nFeature_RNA < nFeature_RNA_max &
                                                                                         percent_mito < pct_mito &
                                                                                         velocyto_unsplicedratio > intronic_ratio_min &
                                                                                         scDblFinder_score < doublet_score )))
    
    # Append the filtered sample to the list
    filtered_list[[sample_name]] <- sample_filtered
    gc(full=T)
  }
  
  # Merge the filtered samples back into a single Seurat object
  filtered_seurat_obj <- merge(filtered_list[[1]], y = filtered_list[-1])
  
  return(filtered_seurat_obj)
}


### Function for cluster markers with DESeq2 
fun_deseq2_cellType <- function(cluster_name, counts, metadata, 
                                smallestGroupSize = 20, use_vst = TRUE, 
                                padj_threshold = 0.05, log2FC_threshold = 0.5, 
                                minExpression = 10, reference_group = NULL) {
  
  if (!all(rownames(metadata) == colnames(counts))) {
    stop("Mismatch between metadata rownames and counts colnames")
  }
  
  run_DE_for_pair <- function(cluster_name, reference_group, counts, metadata) {
    metadata$comparison <- as.character(metadata$cluster)
    metadata$comparison[metadata$cluster == cluster_name] <- cluster_name
    metadata$comparison[metadata$cluster == reference_group] <- reference_group
    metadata$comparison[!(metadata$comparison %in% c(cluster_name, reference_group))] <- "others"
    
    metadata <- metadata[metadata$comparison %in% c(cluster_name, reference_group), ]
    counts <- counts[, rownames(metadata)]
    metadata$comparison <- factor(gsub("\\/", "_", metadata$comparison))
    
    cat("Running DESeq2 for", cluster_name, "vs", reference_group, "\n")
    
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~comparison)
    dds$comparison <- relevel(dds$comparison, reference_group)
    
    keep <- rowSums(counts(dds) >= minExpression) >= smallestGroupSize
    dds <- dds[keep, ]
    cat("Keeping", sum(keep), "genes after filtering.\n")
    
    if (use_vst) {
      vsd <- vst(dds, blind = FALSE)
      zscores <- t(scale(t(assay(vsd))))
    } else {
      rld <- rlog(dds, blind = FALSE)
      zscores <- t(scale(t(assay(rld))))
    }
    
    dds <- DESeq(dds, test = "Wald")
    
    DE_result <- results(dds, contrast = c("comparison", cluster_name, reference_group)) %>%
      as.data.frame() %>%
      filter(padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold) %>%
      arrange(desc(log2FoldChange)) %>%
      mutate(comparison = paste0(cluster_name, "_vs_", reference_group)) %>%
      rownames_to_column("gene")
    
    # Annotate
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    genes <- getBM(filters = "hgnc_symbol", 
                   attributes = c("hgnc_symbol", "gene_biotype"), 
                   values = DE_result$gene, 
                   mart = mart)
    colnames(genes) <- c("gene", "gene_biotype")
    DE_result <- left_join(DE_result, genes, by = "gene")
    DE_result$gene_biotype[grep("^RP[SL]", DE_result$gene)] <- "ribosomal"
    
    return(DE_result)
  }
  
  if (!is.null(reference_group)) {
    # Run only one comparison — cluster_name vs reference_group
    return(run_DE_for_pair(cluster_name, reference_group, counts, metadata))
  } else {
    # Run for all clusters: cluster vs Rest
    all_clusters <- unique(metadata$cluster)
    result_list <- map(all_clusters, function(cl) {
      meta_copy <- metadata
      meta_copy$comparison <- "Rest"
      meta_copy$comparison[meta_copy$cluster == cl] <- cl
      meta_copy$comparison <- factor(meta_copy$comparison)
      
      dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta_copy, design = ~comparison)
      dds$comparison <- relevel(dds$comparison, "Rest")
      
      keep <- rowSums(counts(dds) >= minExpression) >= smallestGroupSize
      dds <- dds[keep, ]
      cat("Running DESeq2 for", cl, "vs Rest. Keeping", sum(keep), "genes.\n")
      
      if (use_vst) {
        vsd <- vst(dds, blind = FALSE)
      } else {
        rld <- rlog(dds, blind = FALSE)
      }
      
      dds <- DESeq(dds, test = "Wald")
      res <- results(dds, contrast = c("comparison", cl, "Rest")) %>%
        as.data.frame() %>%
        filter(padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold) %>%
        arrange(desc(log2FoldChange)) %>%
        mutate(comparison = paste0(cl, "_vs_Rest")) %>%
        rownames_to_column("gene")
      
      if (nrow(res) > 0) {
        mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
        genes <- getBM(filters = "hgnc_symbol", 
                       attributes = c("hgnc_symbol", "gene_biotype"), 
                       values = res$gene, 
                       mart = mart)
        colnames(genes) <- c("gene", "gene_biotype")
        res <- left_join(res, genes, by = "gene")
        res$gene_biotype[grep("^RP[SL]", res$gene)] <- "ribosomal"
      }
      
      return(res)
    })
    
    names(result_list) <- all_clusters
    return(result_list)
  }
}

### Function for creating FeaturePlot objects
create_feature_plot <- function(feature_name, data, gradient_colors, reduction_method = "harmony.umap", limits = c(0, 5), label=NULL, raster = NULL, order=NULL) {
  FeaturePlot(
    object = data,
    reduction = reduction_method,
    features = feature_name,
    cols = gradient_colors,
    label=label, raster = raster, order=order) + scale_colour_gradientn(limits = limits, colors = gradient_colors)
}

### Function to make UCELL violin plots
plot_ucell_vlnbox <- function(column_name,
                              data,
                              x_axis_column = "cell_type",
                              cell_type_order = NULL,
                              colors = NULL,
                              theme_style = c("bw", "minimal")) {
  theme_style <- match.arg(theme_style)
  
  # Group by sample and x_axis_column to calculate mean per sample within each group
  summary_data <- data %>%
    group_by(sample, .data[[x_axis_column]]) %>%
    summarise(sample_mean = mean(.data[[column_name]]), .groups = "drop")
  
  # Set x-axis order if provided
  if (!is.null(cell_type_order)) {
    summary_data[[x_axis_column]] <- factor(summary_data[[x_axis_column]], levels = cell_type_order)
  }
  
  # Default colors if not provided
  if (is.null(colors)) {
    if (x_axis_column == "group") {
      colors <- CFG$groupspal
    } else {
      colors <- c("#AD6279", "#78606D", "#A3B6CB", "#E4E4E4")
    }
  }
  
  # Create the plot
  p <- ggplot(summary_data, aes(x = .data[[x_axis_column]], y = sample_mean, fill = .data[[x_axis_column]])) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.7) +
    geom_quasirandom(shape = 21, color = "black", fill = "white", size = 2, alpha = 0.8) +
    labs(x = NULL, y = NULL, title = column_name) +
    scale_fill_manual(values = colors) +
    NoLegend()
  
  # Apply theme
  if (theme_style == "bw") {
    p <- p + theme_bw()
  } else if (theme_style == "minimal") {
    p <- p + theme_minimal()
  }
  
  # Common theme settings
  p <- p + theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 22),
    axis.text.y = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
  return(p)
}


### Function to process and cluster seurat object
process_sObj <- function(sObj, n_dims=NULL, n_resolution=NULL, default_assay="RNA", integrate=NULL, lambda=NULL){
  DefaultAssay(sObj)<-default_assay
  sObj <- NormalizeData(object = sObj)
  sObj <- FindVariableFeatures(sObj)
  sObj <- ScaleData(sObj)
  sObj <- RunPCA(sObj)
  
  sObj <- FindNeighbors(sObj, dims = n_dims)
  sObj <- FindClusters(sObj, resolution = n_resolution)
  sObj <- RunUMAP(sObj, dims = n_dims, reduction = "pca", reduction.name = "umap.unintegrated")
  
  sObj <- sObj %>% RunHarmony(integrate, plot_convergence = TRUE,  kmeans_init_nstart=20, kmeans_init_iter_max=150, lambda=lambda)
  sObj <- sObj %>% RunUMAP(reduction = "harmony", dims = n_dims, reduction.name='harmony.umap') %>%
    FindNeighbors(reduction = "harmony", dims = n_dims) %>%
    FindClusters(resolution = n_resolution) %>%
    identity()
  
  return(sObj)
}

### Function to calculate percentages
calculate_cocl_pct <- function(file) {
  meta <- read.table(file, sep = "\t", header = TRUE)
  
  filename <- sub("^.*?_metadata_|\\.txt$", "", file)
  filename <- sub("\\.txt$", "", filename)
  
  summary_per_cluster <- meta %>%
    group_by(seurat_clusters, cls) %>%
    summarise(count = n(), .groups = "drop") %>%
    spread(key = cls, value = count, fill = 0) %>%
    mutate(
      total_cells = `Fetal_Neurolin` + `Adult_ImN` + `Adult_ImN-like` + `Adult_MatN`,
      fetal_pct = `Fetal_Neurolin` / total_cells * 100,
      adult_ImN_pct = `Adult_ImN` / total_cells * 100,
      adult_ImN_like_pct = `Adult_ImN-like` / total_cells * 100,
      adult_MatN_pct = `Adult_MatN` / total_cells * 100
    ) %>%
    select(seurat_clusters, fetal_pct, adult_ImN_pct, adult_ImN_like_pct, adult_MatN_pct) %>%
    mutate(combo = filename)
  
  return(summary_per_cluster)
}

### Function to process co-clustering data
filter_AF_MF <- function(summary_per_cluster, type = "AF") {
  summary_per_cluster <- summary_per_cluster %>%
    rowwise() %>%
    mutate(
      AF_total = fetal_pct + adult_ImN_pct,
      MF_total = fetal_pct + adult_MatN_pct,
      AF_like_total = adult_ImN_like_pct + fetal_pct
    ) %>%
    ungroup()
  
  if (type == "AF") {
    filtered <- summary_per_cluster %>%
      filter(max(AF_total, AF_like_total) > MF_total) %>%
      filter(fetal_pct > 2 & adult_ImN_pct > 2)
  } else if (type == "MF") {
    filtered <- summary_per_cluster %>%
      filter(max(AF_total, AF_like_total) < MF_total) %>%
      filter(fetal_pct > 2 & adult_MatN_pct > 2)
  }
  
  return(filtered %>% select(seurat_clusters, fetal_pct, adult_ImN_pct, adult_MatN_pct, 
                             adult_ImN_like_pct, AF_total, MF_total, AF_like_total, combo))
}

### Function to process UMAP with different parameters
fun_umap_param <- function(SObj, theta, dims, res) {
  
  SObj <- SObj %>% 
    RunHarmony(c("sample", "agegroup"), plot_convergence = TRUE,  
               kmeans_init_nstart=20, kmeans_init_iter_max=150, 
               lambda = c(NULL, NULL), theta=c(theta, theta))
  
  SObj <- SObj %>%
    RunUMAP(reduction = "harmony", dims = dims, reduction.name='harmony.umap') %>%
    FindNeighbors(reduction = "harmony", dims = dims) %>%
    FindClusters(resolution = res) %>%
    identity()
  
  barplots <- dittoSeq::dittoBarPlot(SObj, "sub_cell_type", group.by = "seurat_clusters", 
                                     color.panel = c("#DAAC50", "#548FA7", "#A5C2B2", "lightgray")) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.text.x = element_text(size=14),
                       axis.text.y = element_text(size=14), legend.text = element_text(size=14), 
                       legend.title = element_text(size=14), axis.title.x = element_text(size=14),
                       axis.title.y = element_text(size=14))
  
  meta <- rownames_to_column(SObj@meta.data, var = "cell_id") %>%
    dplyr::select(cell_id, seurat_clusters, sub_cell_type)
  
  meta_filename <- paste0("AF_subset_metadata_res", res, "_dim", max(dims), "_theta", theta, ".txt")
  write.table(meta, here("Data","AF_integration", meta_filename), sep="\t", quote = FALSE, row.names = FALSE)
}

# fix the error with dynamicplot
DynamicPlot2<-function (srt, lineages, features, group.by = NULL, cells = NULL, 
                        slot = "counts", assay = NULL, family = NULL, exp_method = c("log1p", 
                                                                                     "raw", "zscore", "fc", "log2fc"), lib_normalize = identical(slot, 
                                                                                                                                                 "counts"), libsize = NULL, compare_lineages = TRUE, 
                        compare_features = FALSE, add_line = TRUE, add_interval = TRUE, 
                        line.size = 1, line_palette = "Dark2", line_palcolor = NULL, 
                        add_point = TRUE, pt.size = 1, point_palette = "Paired", 
                        point_palcolor = NULL, add_rug = TRUE, flip = FALSE, reverse = FALSE, 
                        x_order = c("value", "rank"), aspect.ratio = NULL, legend.position = "right", 
                        legend.direction = "vertical", theme_use = "theme_scp", 
                        theme_args = list(), combine = TRUE, nrow = NULL, ncol = NULL, 
                        byrow = TRUE, seed = 11) 
{
  set.seed(seed)
  check_R("MatrixGenerics")
  x_order <- match.arg(x_order)
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    stop(group.by, " is not in the meta.data of srt object.")
  }
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", 
                      ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  if (length(exp_method) == 1 && is.function(exp_method)) {
    exp_name <- paste0(as.character(x = formals()$exp_method), 
                       "(", data_nm, ")")
  }
  else {
    exp_method <- match.arg(exp_method)
    exp_name <- paste0(exp_method, "(", data_nm, ")")
  }
  assay <- assay %||% DefaultAssay(srt)
  gene <- features[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  features <- c(gene, meta)
  if (length(features) == 0) {
    stop("No feature found in the srt object.")
  }
  cell_union <- c()
  raw_matrix_list <- list()
  fitted_matrix_list <- list()
  upr_matrix_list <- list()
  lwr_matrix_list <- list()
  for (l in lineages) {
    features_exist <- c()
    raw_matrix <- NULL
    if (paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      raw_matrix <- srt@tools[[paste0("DynamicFeatures_", 
                                      l)]][["raw_matrix"]][, -1]
      fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", 
                                         l)]][["fitted_matrix"]][, -1]
      upr_matrix <- srt@tools[[paste0("DynamicFeatures_", 
                                      l)]][["upr_matrix"]][, -1]
      lwr_matrix <- srt@tools[[paste0("DynamicFeatures_", 
                                      l)]][["lwr_matrix"]][, -1]
      features_exist <- colnames(raw_matrix)
    }
    feature_calcu <- features[!features %in% features_exist]
    if (length(feature_calcu) > 0) {
      srt_tmp <- RunDynamicFeatures(srt, lineages = l, 
                                    features = feature_calcu, assay = assay, slot = slot, 
                                    family = family, libsize = libsize)
      if (is.null(raw_matrix)) {
        raw_matrix <- fitted_matrix <- upr_matrix <- lwr_matrix <- matrix(NA, 
                                                                          nrow = nrow(srt_tmp@tools[[paste0("DynamicFeatures_", 
                                                                                                            l)]][["raw_matrix"]]), ncol = 0)
      }
      raw_matrix <- cbind(raw_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", 
                                                            l)]][["raw_matrix"]][, feature_calcu, drop = FALSE])
      fitted_matrix <- cbind(fitted_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", 
                                                                  l)]][["fitted_matrix"]][, feature_calcu, drop = FALSE])
      upr_matrix <- cbind(upr_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", 
                                                            l)]][["upr_matrix"]][, feature_calcu, drop = FALSE])
      lwr_matrix <- cbind(lwr_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", 
                                                            l)]][["lwr_matrix"]][, feature_calcu, drop = FALSE])
    }
    raw_matrix_list[[l]] <- as_matrix(raw_matrix[, features, 
                                                 drop = FALSE])
    fitted_matrix_list[[l]] <- as_matrix(fitted_matrix[, 
                                                       features, drop = FALSE])
    upr_matrix_list[[l]] <- as_matrix(upr_matrix[, features, 
                                                 drop = FALSE])
    lwr_matrix_list[[l]] <- as_matrix(lwr_matrix[, features, 
                                                 drop = FALSE])
    cell_union <- unique(c(cell_union, rownames(raw_matrix)))
  }
  x_assign <- rowMeans(srt@meta.data[cell_union, lineages, 
                                     drop = FALSE], na.rm = TRUE)
  cell_metadata <- cbind.data.frame(data.frame(row.names = cell_union), 
                                    x_assign = x_assign, srt@meta.data[cell_union, lineages, 
                                                                       drop = FALSE])
  cell_order_list <- list()
  for (l in lineages) {
    cell_metadata_sub <- na.omit(cell_metadata[, l, drop = FALSE])
    cell_metadata_sub <- cell_metadata_sub[order(cell_metadata_sub[[l]], 
                                                 decreasing = FALSE), , drop = FALSE]
    cell_order_list[[l]] <- paste0(rownames(cell_metadata_sub), 
                                   l)
  }
  df_list <- list()
  Y_libsize <- colSums(slot(srt@assays[[assay]], "counts"))
  for (l in lineages) {
    raw_matrix <- raw_matrix_list[[l]]
    fitted_matrix <- fitted_matrix_list[[l]]
    upr_matrix <- upr_matrix_list[[l]]
    lwr_matrix <- lwr_matrix_list[[l]]
    if (isTRUE(lib_normalize) && min(raw_matrix[, gene], 
                                     na.rm = TRUE) >= 0) {
      if (!is.null(libsize)) {
        libsize_use <- libsize
      }
      else {
        libsize_use <- Y_libsize[rownames(raw_matrix)]
        isfloat <- any(libsize_use%%1 != 0, na.rm = TRUE)
        if (isTRUE(isfloat)) {
          libsize_use <- rep(1, length(libsize_use))
          warning("The values in the 'counts' slot are non-integer. Set the library size to 1.", 
                  immediate. = TRUE)
        }
      }
      raw_matrix[, gene] <- raw_matrix[, gene, drop = FALSE]/libsize_use * 
        median(Y_libsize)
    }
    if (is.function(exp_method)) {
      raw_matrix <- t(exp_method(t(raw_matrix)))
      fitted_matrix <- t(exp_method(t(fitted_matrix)))
      upr_matrix <- t(exp_method(t(upr_matrix)))
      lwr_matrix <- t(exp_method(t(lwr_matrix)))
    }
    else if (exp_method == "raw") {
      raw_matrix <- raw_matrix
      fitted_matrix <- fitted_matrix
      upr_matrix <- upr_matrix
      lwr_matrix <- lwr_matrix
    }
    else if (exp_method == "zscore") {
      center <- colMeans(raw_matrix)
      sd <- MatrixGenerics::colSds(raw_matrix)
      raw_matrix <- scale(raw_matrix, center = center, 
                          scale = sd)
      fitted_matrix <- scale(fitted_matrix, center = center, 
                             scale = sd)
      upr_matrix <- scale(upr_matrix, center = center, 
                          scale = sd)
      lwr_matrix <- scale(lwr_matrix, center = center, 
                          scale = sd)
    }
    else if (exp_method == "fc") {
      colm <- colMeans(raw_matrix)
      raw_matrix <- t(t(raw_matrix)/colm)
      fitted_matrix <- t(t(fitted_matrix)/colm)
      upr_matrix <- t(t(upr_matrix)/colm)
      lwr_matrix <- t(t(lwr_matrix)/colm)
    }
    else if (exp_method == "log2fc") {
      colm <- colMeans(raw_matrix)
      raw_matrix <- t(log2(t(raw_matrix)/colm))
      fitted_matrix <- t(log2(t(fitted_matrix)/colm))
      upr_matrix <- t(log2(t(upr_matrix)/colm))
      lwr_matrix <- t(log2(t(lwr_matrix)/colm))
    }
    else if (exp_method == "log1p") {
      raw_matrix <- log1p(raw_matrix)
      fitted_matrix <- log1p(fitted_matrix)
      upr_matrix <- log1p(upr_matrix)
      lwr_matrix <- log1p(lwr_matrix)
    }
    raw_matrix[is.infinite(raw_matrix)] <- max(abs(raw_matrix[!is.infinite(raw_matrix)]), 
                                               na.rm = TRUE) * ifelse(raw_matrix[is.infinite(raw_matrix)] > 
                                                                        0, 1, -1)
    fitted_matrix[is.infinite(fitted_matrix)] <- max(abs(fitted_matrix[!is.infinite(fitted_matrix)])) * 
      ifelse(fitted_matrix[is.infinite(fitted_matrix)] > 
               0, 1, -1)
    upr_matrix[is.infinite(upr_matrix)] <- max(abs(upr_matrix[!is.infinite(upr_matrix)]), 
                                               na.rm = TRUE) * ifelse(upr_matrix[is.infinite(upr_matrix)] > 
                                                                        0, 1, -1)
    lwr_matrix[is.infinite(lwr_matrix)] <- max(abs(lwr_matrix[!is.infinite(lwr_matrix)]), 
                                               na.rm = TRUE) * ifelse(lwr_matrix[is.infinite(lwr_matrix)] > 
                                                                        0, 1, -1)
    raw <- as.data.frame(cbind(cell_metadata[rownames(raw_matrix), 
                                             c(l, "x_assign")], raw_matrix))
    colnames(raw)[1] <- "Pseudotime"
    raw[["Cell"]] <- rownames(raw)
    raw[["Value"]] <- "raw"
    raw <- reshape2::melt(raw, id.vars = c("Cell", "Pseudotime", "x_assign", 
                                 "Value"), value.name = "exp", variable.name = "Features")
    fitted <- as.data.frame(cbind(cell_metadata[rownames(fitted_matrix), 
                                                c(l, "x_assign")], fitted_matrix))
    colnames(fitted)[1] <- "Pseudotime"
    fitted[["Cell"]] <- rownames(fitted)
    fitted[["Value"]] <- "fitted"
    fitted <- reshape2::melt(fitted, id.vars = c("Cell", "Pseudotime", 
                                       "x_assign", "Value"), value.name = "exp", variable.name = "Features")
    upr <- as.data.frame(cbind(cell_metadata[rownames(upr_matrix), 
                                             c(l, "x_assign")], upr_matrix))
    colnames(upr)[1] <- "Pseudotime"
    upr[["Cell"]] <- rownames(upr)
    upr[["Value"]] <- "upr"
    upr <- reshape2::melt(upr, id.vars = c("Cell", "Pseudotime", "x_assign", 
                                 "Value"), value.name = "exp", variable.name = "Features")
    lwr <- as.data.frame(cbind(cell_metadata[rownames(lwr_matrix), 
                                             c(l, "x_assign")], lwr_matrix))
    colnames(lwr)[1] <- "Pseudotime"
    lwr[["Cell"]] <- rownames(lwr)
    lwr[["Value"]] <- "lwr"
    lwr <- reshape2::melt(lwr, id.vars = c("Cell", "Pseudotime", "x_assign", 
                                 "Value"), value.name = "exp", variable.name = "Features")
    raw[["upr"]] <- NA
    raw[["lwr"]] <- NA
    fitted[["upr"]] <- upr[["exp"]]
    fitted[["lwr"]] <- lwr[["exp"]]
    df_tmp <- rbind(raw, fitted)
    df_tmp[["Lineages"]] <- factor(l, levels = lineages)
    df_list[[l]] <- df_tmp
  }
  df_all <- do.call(rbind, df_list)
  rownames(df_all) <- NULL
  if (!is.null(group.by)) {
    cell_group <- srt@meta.data[df_all[["Cell"]], group.by, 
                                drop = FALSE]
    if (!is.factor(cell_group[, group.by])) {
      cell_group[, group.by] <- factor(cell_group[, group.by], 
                                       levels = unique(cell_group[, group.by]))
    }
    df_all <- cbind(df_all, cell_group)
  }
  df_all[["LineagesFeatures"]] <- paste(df_all[["Lineages"]], 
                                        df_all[["Features"]], sep = "-")
  if (!is.null(cells)) {
    df_all <- df_all[df_all[["Cell"]] %in% cells, , drop = FALSE]
  }
  df_all <- df_all[sample(seq_len(nrow(df_all))), , drop = FALSE]
  plist <- list()
  legend <- NULL
  if (isTRUE(compare_lineages)) {
    lineages_use <- list(lineages)
    lineages_formula <- "."
  }
  else {
    lineages_use <- lineages
    lineages_formula <- "Lineages"
  }
  if (isTRUE(compare_features)) {
    features_use <- list(features)
    features_formula <- "."
  }
  else {
    features_use <- features
    features_formula <- "Features"
  }
  formula <- paste(lineages_formula, "~", features_formula)
  fill_by <- "Lineages"
  if (lineages_formula == "." && length(lineages) > 1) {
    lineages_guide <- TRUE
  }
  else {
    lineages_guide <- FALSE
    if (isTRUE(compare_features)) {
      fill_by <- "Features"
    }
  }
  if (features_formula == "." && length(features) > 1) {
    features_guide <- TRUE
  }
  else {
    features_guide <- FALSE
  }
  for (l in lineages_use) {
    for (f in features_use) {
      df <- subset(df_all, df_all[["Lineages"]] %in% l & 
                     df_all[["Features"]] %in% f)
      if (x_order == "rank") {
        df[, "x_assign"] <- rank(df[, "x_assign"])
        df[, "Pseudotime"] <- rank(df[, "Pseudotime"])
      }
      df_point <- unique(df[df[["Value"]] == "raw", c("Cell", 
                                                      "x_assign", "exp", group.by)])
      if (isTRUE(compare_features)) {
        raw_point <- NULL
      }
      else {
        if (isTRUE(add_point)) {
          if (is.null(group.by)) {
            raw_point <- geom_point(data = df_point, 
                                    mapping = aes(x = .data[["x_assign"]], 
                                                  y = .data[["exp"]]), size = pt.size, 
                                    alpha = 0.8)
          }
          else {
            raw_point <- list(geom_point(data = df_point, 
                                         mapping = aes(x = .data[["x_assign"]], 
                                                       y = .data[["exp"]], color = .data[[group.by]]), 
                                         size = pt.size, alpha = 0.8), scale_color_manual(values = palette_scp(df[[group.by]], 
                                                                                                               palette = point_palette, palcolor = point_palcolor)), 
                              scale_fill_manual(values = palette_scp(df[[group.by]], 
                                                                     palette = point_palette, palcolor = point_palcolor), 
                                                guide = guide_legend(override.aes = list(alpha = 1, 
                                                                                         size = 3), order = 1)), new_scale_color(), 
                              new_scale_fill())
          }
        }
        else {
          raw_point <- NULL
        }
      }
      if (isTRUE(add_rug)) {
        if (is.null(group.by)) {
          rug <- list(geom_rug(data = df_point, mapping = aes(x = .data[["x_assign"]]), 
                               alpha = 1, length = unit(0.05, "npc"), show.legend = FALSE))
        }
        else {
          rug <- list(geom_rug(data = df_point, mapping = aes(x = .data[["x_assign"]], 
                                                              color = .data[[group.by]]), alpha = 1, length = unit(0.05, 
                                                                                                                   "npc"), show.legend = isTRUE(compare_features)), 
                      scale_color_manual(values = palette_scp(df[[group.by]], 
                                                              palette = point_palette, palcolor = point_palcolor)), 
                      new_scale_color())
        }
      }
      else {
        rug <- NULL
      }
      if (isTRUE(add_interval)) {
        interval <- list(geom_ribbon(data = subset(df, 
                                                   df[["Value"]] == "fitted"), mapping = aes(x = .data[["Pseudotime"]], 
                                                                                             y = .data[["exp"]], ymin = .data[["lwr"]], 
                                                                                             ymax = .data[["upr"]], fill = .data[[fill_by]], 
                                                                                             group = .data[["LineagesFeatures"]]), alpha = 0.4, 
                                     color = "grey90"), scale_fill_manual(values = palette_scp(df[[fill_by]], 
                                                                                               palette = line_palette, palcolor = line_palcolor), 
                                                                          guide = if (fill_by == "Features" || lineages_guide || 
                                                                                      length(l) == 1) "none" else guide_legend()), 
                         new_scale_fill())
      }
      else {
        interval <- NULL
      }
      if (isTRUE(compare_features)) {
        line <- list(geom_line(data = subset(df, df[["Value"]] == 
                                               "fitted"), mapping = aes(x = .data[["Pseudotime"]], 
                                                                        y = .data[["exp"]], color = .data[["Features"]], 
                                                                        group = .data[["LineagesFeatures"]]), linewidth = line.size, 
                               alpha = 0.8), scale_color_manual(values = palette_scp(df[["Features"]], 
                                                                                     palette = line_palette, palcolor = line_palcolor), 
                                                                guide = if (features_guide) guide_legend(override.aes = list(alpha = 1, 
                                                                                                                             size = 2), order = 2) else "none"), new_scale_color())
      }
      else {
        if (isTRUE(add_line)) {
          line <- list(geom_line(data = subset(df, df[["Value"]] == 
                                                 "fitted"), mapping = aes(x = .data[["Pseudotime"]], 
                                                                          y = .data[["exp"]], color = .data[["Lineages"]], 
                                                                          group = .data[["LineagesFeatures"]]), linewidth = line.size, 
                                 alpha = 0.8), scale_color_manual(values = palette_scp(df[["Lineages"]], 
                                                                                       palette = line_palette, palcolor = line_palcolor), 
                                                                  guide = if (lineages_guide) guide_legend(override.aes = list(alpha = 1, 
                                                                                                                               size = 2), order = 2) else "none"), new_scale_color())
        }
        else {
          line <- NULL
        }
      }
      x_trans <- ifelse(flip, "reverse", "identity")
      x_trans <- ifelse(reverse, setdiff(c("reverse", 
                                           "identity"), x_trans), x_trans)
      p <- ggplot() + scale_x_continuous(trans = x_trans, 
                                         expand = expansion(c(0, 0))) + scale_y_continuous(expand = expansion(c(0.1, 
                                                                                                                0.05))) + raw_point + rug + interval + list(line[[1]], line[[2]]) + 
        labs(x = ifelse(x_order == "rank", "Pseudotime(rank)", 
                        "Pseudotime"), y = exp_name) + facet_grid(formula(formula), 
                                                                  scales = "free") + do.call(theme_use, theme_args) + 
        theme(aspect.ratio = aspect.ratio, legend.position = legend.position, 
              legend.direction = legend.direction)
      if (isTRUE(flip)) {
        p <- p + coord_flip()
      }
      if (is.null(legend)) {
        legend <- get_plot_component(p, 'guide-box-right', return_all = TRUE)
      }
      plist[[paste(paste0(l, collapse = "_"), paste0(f, 
                                                     collapse = "_"), sep = ".")]] <- p + theme(legend.position = "none")
    }
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow,
                         ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    
    if (legend.position != "none") {
      legend <- get_plot_component(p, 'guide-box-right', return_all = TRUE)
      gtable <- as_grob(plot)
      gtable <- gridExtra::grid.arrange(grobs = list(gtable, legend), ncol = 2)
      plot <- wrap_plots(gtable)
      }
    }
    return(plot)
}


unlockBinding("DynamicPlot", asNamespace("SCP"))
assign("DynamicPlot", DynamicPlot2, envir = asNamespace("SCP"))

### function to make barplots of genes from cortical data
genebarplot <- function(gene_symbol) {
  rpkm <- getRPKM(rse_gene)
  gene_data <- rowRanges(rse_gene)
  
  # Find the index for the specified gene
  gene_index <- which(gene_data$Symbol == gene_symbol)
  
  if (length(gene_index) == 0) {
    warning(paste("Gene", gene_symbol, "not found in the dataset."))
    return(NULL)  # Return NULL if the gene is not found
  }
  
  plot_data <- data.frame(
    Day = colData(rse_gene)$DAY,
    rpkm = rpkm[gene_index, ],
    log2rpkm = log2(rpkm[gene_index, ] +1)
  )
  
  plot_data$Day <- paste0("D", plot_data$Day)
  
  plot_data <- plot_data %>%
    mutate(Day_numeric = as.numeric(gsub("D", "", Day))) %>%
    arrange(Day_numeric) %>%
    mutate(Day = factor(Day, levels = unique(Day)))
  
  # Calculate summary statistics: mean and standard error
  summary_stats <- plot_data %>%
    group_by(Day) %>%
    summarise(
      mean_log2rpkm = mean(log2rpkm, na.rm = TRUE),
      se_log2rpkm = sd(log2rpkm, na.rm = TRUE) / sqrt(n())
    )
  
  # Ensure all values are finite before plotting
  plot_data <- plot_data %>% filter(is.finite(log2rpkm))
  summary_stats <- summary_stats %>% filter(is.finite(mean_log2rpkm))
  
  # Create the plot
  p <- ggplot(summary_stats, aes(x = Day, y = mean_log2rpkm)) +
    geom_bar(stat = "identity", fill = "#A0A0A4", color = "black") +
    geom_errorbar(aes(ymin = mean_log2rpkm - se_log2rpkm, ymax = mean_log2rpkm + se_log2rpkm), width = 0.2, color = "black") +
    geom_point(data = plot_data, aes(x = Day, y = log2rpkm), position = position_jitter(width = 0.1), size = 2, alpha = 0.8) +
    labs(title = paste("Expression of", gene_symbol, "over Days in Vitro"),
         x = "Days in Vitro", 
         y = "Expression (log2(RPKM + 1))") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 22), 
          axis.text.y = element_text(size = 22), 
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          plot.title = element_text(size = 24, face = "bold"),
          panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    coord_cartesian(ylim = c(0, 10))
  
  print(p)
}


### Function to calculate DEGs across groups
fun_deseq2_groups <- function(
    counts_list,
    metadata_all,
    metadata_info,
    pairs,
    smallestGroupSize = 4,
    use_rlog = TRUE,
    save_corrected_rlog = FALSE,
    padj_threshold = 0.1,
    log2FC_threshold = 0,
    minExpression = 10
) {
  
  clusters <- unique(metadata_all$cellType)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  DE_result_list <- list()
  
  for (cluster in clusters) {
    for (pair in pairs) {
      target <- pair[1]
      reference <- pair[2]
      
      metadata <- merge(metadata_all, metadata_info, by = c("sample")) %>%
        filter(group %in% c(target, reference)) %>%
        filter(cellType == cluster)
      rownames(metadata) <- metadata$sample
      
      cluster_counts <- aggregated_counts_split[[cluster]][, colnames(aggregated_counts_split[[cluster]]) %in% rownames(metadata)]
      cluster_metadata <- metadata[match(colnames(cluster_counts), rownames(metadata)), ]
      
      cluster_metadata <- cluster_metadata %>%
        group_by(Run) %>%
        mutate(run_count = n()) %>%
        ungroup()
      
      single_run_samples <- cluster_metadata %>% filter(run_count == 1)
      
      if (nrow(single_run_samples) > 1) {
        message("These were single Run samples, they will be grouped together.")
        print(single_run_samples)
        cluster_metadata <- cluster_metadata %>%
          mutate(Run = ifelse(run_count == 1, "pooledRun", Run))
      } else if (nrow(single_run_samples) == 1) {
        message("This is the only sample that has a single run, it will be grouped with next run.")
        print(single_run_samples)
        next_run <- cluster_metadata %>%
          filter(run_count > 1) %>%
          arrange(Run) %>%
          dplyr::slice(1) %>%
          pull(Run)
        cluster_metadata <- cluster_metadata %>%
          mutate(Run = ifelse(run_count == 1, next_run, Run))
      }
      
      cluster_metadata <- cluster_metadata %>%
        dplyr::select(-run_count) %>%
        column_to_rownames("sample")
      
      if (!all(rownames(cluster_metadata) == colnames(cluster_counts))) stop("Sample mismatch!")
      
      message("Running DESeq2 for ", cluster, " ", target, "_vs_", reference, " correcting for Sex and Run")
      cluster_metadata$group <- factor(cluster_metadata$group)
      cluster_metadata$Sex <- factor(cluster_metadata$Sex)
      cluster_metadata$Run <- factor(cluster_metadata$Run)
      
      dds <- DESeqDataSetFromMatrix(cluster_counts,
                                    colData = cluster_metadata,
                                    design = ~0 + Sex + Run + group)
      
      keep <- rowSums(counts(dds) >= minExpression) >= smallestGroupSize
      dds <- dds[keep,]
      
      dds$condition <- factor(cluster_metadata$group)
      
      if (use_rlog) {
        rld <- rlog(dds, blind = FALSE)
      } else {
        vsd <- vst(dds, blind = FALSE)
      }
      
      covariates <- model.matrix(~cluster_metadata$Sex + cluster_metadata$Run)
      covariates <- covariates[,-1]
      mat <- assay(rld)
      mm <- model.matrix(~group, colData(rld))
      mat <- limma::removeBatchEffect(mat, covariates = covariates, design = mm)
      assay(rld) <- mat
      
      if (save_corrected_rlog) {
        saveRDS(rld, here("Data", paste0("rlog_corrected_", cluster, "_", target, "_vs_", reference, ".RDS")))
      }
      
      dds <- DESeq(dds, test = "Wald")
      DE_result <- as.data.frame(results(dds, contrast = c("group", target, reference))) %>%
        filter(padj < padj_threshold) %>%
        arrange(desc(log2FoldChange)) %>%
        mutate(group = paste0(cluster, "-", target, "_vs_", reference)) %>%
        rownames_to_column('gene')
      
      if (nrow(DE_result) > 2) {
        genes <- getBM(filters = 'hgnc_symbol', attributes = c('hgnc_symbol', 'gene_biotype'), values = DE_result$gene, mart = mart)
        colnames(genes) <- c('gene', 'gene_biotype')
        DE_result <- left_join(DE_result, genes, by = "gene") %>% arrange(desc(log2FoldChange))
        DE_result$gene_biotype[grep("^RP[SL]", DE_result$gene)] <- "ribosomal"
        
        subset_data <- DE_result[, !(colnames(DE_result) %in% c("gene_biotype"))]
        duplicated_rows <- duplicated(subset_data) | duplicated(subset_data, fromLast = TRUE)
        DE_result <- DE_result[!duplicated(subset_data), ]
      }
      
      DE_result_list[[paste0(cluster, "-", target, "_vs_", reference)]] <- DE_result
    }
  }
  
  all_columns <- unique(unlist(lapply(DE_result_list, colnames)))
  df_list_fixed <- lapply(DE_result_list, function(df) {
    if (nrow(df) == 0) {
      df <- data.frame(matrix(ncol = length(all_columns), nrow = 0))
      colnames(df) <- all_columns
    } else {
      missing_columns <- setdiff(all_columns, colnames(df))
      df[missing_columns] <- NA
    }
    return(df)
  })
  
  df_list_fixed <- df_list_fixed[sapply(df_list_fixed, nrow) > 0]
  DE_result_all <- do.call(rbind, df_list_fixed)
  
  result_summary <- DE_result_all %>%
    group_by(group) %>%
    summarize(
      num_genes_above_0 = sum(log2FoldChange > 0),
      num_genes_below_0 = sum(log2FoldChange < 0)
    )
  
  return(list(results = DE_result_all, summary = result_summary))
}

list_files_recursive <- function(path) {
  list.files(path, recursive = TRUE, full.names = TRUE)
}

## Functions for cellphoneDB processing
read_and_preprocess <- function(file_paths) {
  data_list <- list()
  for (file_path in file_paths) {
    data <- read.delim(file_path, check.names = FALSE)
    data <- data %>%
      mutate(
        gene_a = na_if(gene_a, ""),
        gene_b = na_if(gene_b, ""),
        gene_pair = paste(gene_a, gene_b, sep = "-")
      )
    data_list[[file_path]] <- data
  }
  return(data_list)
}

# Extract group name from file path
extract_group_name <- function(path) {
  basename(dirname(dirname(path)))
}

# Combine forward and reverse interaction data
combine_pairs <- function(fw, rv) {
  if (!is.null(fw) & !is.null(rv)) {
    return(bind_rows(fw, rv))
  }
  return(NULL)
}

# Create combined name for interaction pair
create_name <- function(fw_name) {
  fw_prefix <- sub("\\|.*", "", fw_name)
  fw_suffix <- sub(".*\\|", "", fw_name)
  return(paste0(fw_prefix, "_and_", fw_suffix))
}

filter_interactions <- function(df) {
  interactions_by_group <- df %>%
    group_by(group) %>%
    summarise(id_cluster = list(id_cluster)) %>%
    deframe()
  
  unique_SAD <- setdiff(interactions_by_group$SAD, c(interactions_by_group$MAD, interactions_by_group$CTRL, interactions_by_group$RES))
  unique_MAD <- setdiff(interactions_by_group$MAD, c(interactions_by_group$SAD, interactions_by_group$CTRL, interactions_by_group$RES))
  unique_CTRL <- setdiff(interactions_by_group$CTRL, c(interactions_by_group$SAD, interactions_by_group$MAD, interactions_by_group$RES))
  unique_RES <- setdiff(interactions_by_group$RES, c(interactions_by_group$SAD, interactions_by_group$MAD, interactions_by_group$CTRL))
  
  shared_SAD_MAD <- intersect(interactions_by_group$SAD, interactions_by_group$MAD)
  shared_RES_MAD <- intersect(interactions_by_group$RES, interactions_by_group$MAD)
  shared_CTRL_RES <- intersect(interactions_by_group$CTRL, interactions_by_group$RES)
  shared_RES_MAD_CTRL <- Reduce(intersect, list(interactions_by_group$RES, interactions_by_group$MAD, interactions_by_group$CTRL))
  
  filtered_df_result <- bind_rows(
    df %>% filter(group == "SAD" & id_cluster %in% unique_SAD) %>% mutate(annotation = "SAD_unique"),
    df %>% filter(group == "MAD" & id_cluster %in% unique_MAD) %>% mutate(annotation = "MAD_unique"),
    df %>% filter(group == "RES" & id_cluster %in% unique_RES) %>% mutate(annotation = "RES_unique"),
    df %>% filter(group == "CTRL" & id_cluster %in% unique_CTRL) %>% mutate(annotation = "CTRL_unique"),
    df %>% filter(group %in% c("SAD", "MAD") & id_cluster %in% shared_SAD_MAD) %>% mutate(annotation = "SAD_MAD_shared"),
    df %>% filter(group %in% c("RES", "MAD") & id_cluster %in% shared_RES_MAD) %>% mutate(annotation = "RES_MAD_shared"),
    df %>% filter(group %in% c("CTRL", "RES") & id_cluster %in% shared_CTRL_RES) %>% mutate(annotation = "CTRL_RES_shared"),
    df %>% filter(group %in% c("CTRL", "RES", "MAD") & id_cluster %in% shared_RES_MAD_CTRL) %>% mutate(annotation = "RES_MAD_CTRL_shared")
  )
  
  return(filtered_df_result)
}