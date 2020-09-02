#' @title Heatmap of Peptide Hits
#' @description takes input of phage results and generates a heatmap of all peptides
#' then returns variable containing the plot
#' @param df phage peptide count matrix
#' @param patient patient of interest
#' @param target.samples vector of target samples of interest from one patient (defined by 'patient' input)
#' @param disease.samples vector of all disease samples in input data frame
#' @param reference.samples vector of all reference samples (healthy and bead controls) in input data frame
#' @param normalization.method method used to normalize data <'log10', 'minmax','meannorm'>
#' @param disease.annot column annotation for disease samples
#' @param reference.annot column annotation for reference samples
#' @param gene.col column name for gene name in input data frame <default: "gene">
#' @param peptide.id.col column name for peptide IDs in input data frame <default: "peptide_id">
#' @return variable containing resulting heatmap.
#' @examples makeHeatmap(df,patient,target.samples,disease.samples,reference.samples)
#' @import pheatmap
#' @import gplots
#' @import data.table
#' @export

makeHeatmap <- function(df,patient,gene.col = "gene",peptide.id.col = "peptide_id",target.samples,
                        disease.samples,reference.samples, normalization.method = "log10",
                        disease.annot = "CoV",reference.annot = "Reference") {
  # Set the margins for heatmaps
  par(mar=c(1,1,1,1))

  row.names(df) <- df[,peptide.id.col]
  df[,peptide.id.col] <- NULL
  df[,gene.col]       <- NULL


  disease.samples.rem   <- disease.samples[! disease.samples %in% target.samples]
  reference.samples.rem <- reference.samples[! reference.samples %in% target.samples]

  df <- df[,c(target.samples,disease.samples.rem,reference.samples.rem)]

  ## Normalize Data

  # Min Max normalize
  min_max_norm <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  # Mean normalize
  mean_norm <- function(x) {
    (x - mean(x)) / (max(x) - min(x))
  }

  # apply Min-Max normalization to first four columns in iris dataset
  df_norm            <- as.data.frame(lapply(df, min_max_norm))
  row.names(df_norm) <- row.names(df)
  df_norm            <- as.matrix(df_norm)

  # apply Mean normalization to first four columns in iris dataset
  df_mnorm            <- as.data.frame(lapply(df, mean_norm))
  row.names(df_mnorm) <- row.names(df)
  df_mnorm            <- as.matrix(df_mnorm)

  # Log10 transform
  df_log10                   <- log10(as.matrix(df))
  df_log10[df_log10 == -Inf] <- 0


  df_fmt <- as.data.frame(df)
  if (normalization.method == "log10")    {df_fmt <- as.data.frame(df_log10)}
  if (normalization.method == "minmax")   {df_fmt <- as.data.frame(df_norm)}
  if (normalization.method == "meannorm") {df_fmt <- as.data.frame(df_mnorm)}



  ## Get Gene order
  heatmap_mtx   <- heatmap.2(as.matrix(df_fmt))
  all_peptides  <- row.names(df_fmt)
  peptide_order <- all_peptides[heatmap_mtx$rowInd]
  gene_list    <- unique(tstrsplit(peptide_order,"_")[[1]])

  df_fmt[,gene.col] <- tstrsplit(row.names(df_fmt),"_")[[1]]

  gene_list_ordered <- c()
  for (gene in gene_list) {
    gene_df      <- df_fmt[df_fmt[,gene.col] %in% gene,]
    gene_df[,gene.col] <- NULL

    gene_peptides_order <- c()
    if (nrow(gene_df) == 1){
      gene_peptides_order <- row.names(gene_df)
    } else {
      if(nrow(gene_df) == 0){next}
      gene_hp_mtx  <- heatmap.2(as.matrix(gene_df))
      gene_peptides       <- row.names(gene_df)
      gene_peptides_order <- gene_peptides[gene_hp_mtx$rowInd]
    }

    gene_list_ordered <- c(gene_list_ordered,gene_peptides_order)
  }

  ## Re-order list of genes in heatmap
  df_fmt <- df_fmt[gene_list_ordered,]


  ## Annotation column Samples
  sample_vec <- names(as.data.frame(df_fmt))
  sample_vec <- sample_vec[! sample_vec %in% c(gene.col,peptide.id.col)]
  diagnosis  <- sample_vec
  diagnosis[diagnosis %in% disease.samples]    <- disease.annot
  diagnosis[diagnosis %in% reference.samples]  <- reference.annot

  diagnosis_annot <- as.data.frame(diagnosis)
  row.names(diagnosis_annot) <- sample_vec

  ## Annotation column Genes
  peptide_vec <- row.names(as.data.frame(df_fmt))
  genes_vec   <- tstrsplit(peptide_vec,"_")[[1]]
  genes_annot <- as.data.frame(genes_vec)
  row.names(genes_annot) <- peptide_vec

  ## Color Scheme

  break_list1 <- seq(-1,1,0.02)  # form minmax norm
  break_list2 <- seq(-2,2,0.04)  # -2 to 2
  break_list3 <- seq(-3,3,0.06)  # -3 to 3
  break_list4 <- seq(-4,4,0.08)  # -4 to 4

  break_list = NULL
  if (normalization.method == "log10")   {break_list = break_list2}
  if (normalization.method == "minmax")  {break_list = break_list1}
  if (normalization.method == "meannorm"){break_list = break_list1}

  color_pal  <- bluered(100)

  ## Heatmap command
  if (gene.col %in% names(df_fmt)){df_fmt[,gene.col] <- NULL}

  # Column annotation bars colors
  colcolor        <- c("purple", "darkgreen")
  names(colcolor) <- c(disease.annot, reference.annot)
  anno_colors <- list(diagnosis = colcolor)

  p <- pheatmap(df_fmt,cluster_cols = F, cluster_rows = F,
                color = color_pal,breaks = break_list,fontsize_row = 4.5,
                annotation_col = diagnosis_annot, annotation_row = genes_annot,
                annotation_legend = F,annotation_names_col = F, annotation_names_row = F,annotation_colors = anno_colors,
                gaps_col = c(2),cellheight = 5)

  return(p)
}
