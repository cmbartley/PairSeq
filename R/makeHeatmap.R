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
#' @param min.rpk For heatmaps set 0 values to an inputted minimum RPK value <default: 0>
#' @param break.list This determines the margins for the resultsing heatmap. the only options are (1,2,3,4) <default: NULL>
#' @param gene.order Vector of gene names to force the order in which peptides are displayed in the heatmap <default: c()>
#' @param avg.non.target.columns Average the replicates for all samples that are not in target.samples <default: FALSE>
#' @param disease.first option for where disease group apears first on the heatmap <default: TRUE>
#' @param cell.height option for cell height in heatmap <default: 5>
#' @param cell.width option for cell width in heatmap <default: 5>
#' @param num.pos.controls option to indicate number of AG Bead controls in the data frame, if you are averaging replicates, keep this number at 1 <default: 1>
#' @param null.ip.col option that indicates a regex to group all AG Bead columns <default: Bead >
#' @param transpose.heatmap option to transpose the heatmaps columns and rows <default: FALSE>
#' @param custom.gene.labels option to change row labels to gene names only labeled once per gene block <default: TRUE>
#' @param chop.length number characters to remove from the end of each sample name to remove the replicate IDs <default: 2>
#' @param split.by insert regex to split sample strings by and keep the first element, which helps simplify complex sample names <default: NULL>
#' @return variable containing resulting heatmap.
#' @examples makeHeatmap(df,patient,target.samples,disease.samples,reference.samples)
#' @import pheatmap
#' @import gplots
#' @import data.table
#' @export

makeHeatmap <- function(df,patient,gene.col = "gene",peptide.id.col = "peptide_id",target.samples,
                        disease.samples,reference.samples, normalization.method = "log10",
                        disease.annot = "COVID-19+",reference.annot = "REFERENCE",min.rpk = 0,
                        break.list = NULL, gene.order = c(),avg.non.target.columns = FALSE,disease.first = TRUE,
                        cell.height=5,cell.width=5, num.pos.controls = 1,null.ip.col = "Bead",transpose.heatmap=FALSE,
                        custom.gene.labels=TRUE,chop.length = 2,split.by = NULL) {
  # Set the margins for heatmaps
  par(mar=c(1,1,1,1))

  row.names(df) <- df[,peptide.id.col]
  df[,peptide.id.col] <- NULL
  df[,gene.col]       <- NULL


  disease.samples.rem   <- disease.samples[! disease.samples %in% target.samples]
  reference.samples.rem <- reference.samples[! reference.samples %in% target.samples]

  # Order by replicate
  orderSamples <- function(sample.vec) {
    sample.vec.order <- c()

    name_vec <- c()
    for (i in 1:length(sample.vec)) {
      sample <- sample.vec[i]
      name  <- ""
      if(grepl(null.ip.col,sample)){
        name <- unlist(strsplit(sample,"_"))[1]
      } else if (grepl("GFAP",sample)){
        name <- unlist(strsplit(sample,"_"))[1]
      } else{
        name_comp <- unlist(strsplit(sample,"_"))
        name <- paste(name_comp[2:(length(name_comp)-1)],collapse = "_")
      }
      name_vec <- c(name_vec,name)
    }
    name_vec <- unique(name_vec)

    # Now order samples
    for (name in name_vec) {
      pt.vec <- sort(sample.vec[grep(name,sample.vec)])
      sample.vec.order <- c(sample.vec.order,pt.vec)
    }

    return(sample.vec.order)
  }

  disease.samples.rem.order   <- c()
  if (length(disease.samples.rem) > 0) {disease.samples.rem.order   <- unique(orderSamples(disease.samples.rem))}
  reference.samples.rem.order <- c()
  if (length(reference.samples.rem) > 0) {reference.samples.rem.order <- unique(orderSamples(reference.samples.rem))}

  df <- df[,c(target.samples,disease.samples.rem.order,reference.samples.rem.order)]

  ## Average non-target columns columns ----
  disease.samples.orig   <- disease.samples
  reference.samples.orig <- reference.samples
  if (avg.non.target.columns == TRUE){
    df <- avgReplicates(df,target.samples,disease.samples,reference.samples, null.ip.col = null.ip.col,chop.length = chop.length,split.by = split.by)
    disease.samples   <- collapseReplicateColumnNames(disease.samples, null.ip.col=null.ip.col,chop.length = chop.length,split.by = split.by)
    reference.samples <- collapseReplicateColumnNames(reference.samples, null.ip.col=null.ip.col,chop.length = chop.length,split.by = split.by)

    if(length(disease.samples.rem.order) > 0){disease.samples.rem.order   <- collapseReplicateColumnNames(disease.samples.rem.order,null.ip.col=null.ip.col,chop.length = chop.length,split.by = split.by)}
    if (length(reference.samples.rem.order) > 0){reference.samples.rem.order <- collapseReplicateColumnNames(reference.samples.rem.order,null.ip.col=null.ip.col,chop.length = chop.length,split.by = split.by)}
  }


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
  mtx           <- as.matrix(df)
  mtx[mtx == 0] <- min.rpk
  df_log10      <- log10(mtx)
  #df_log10[df_log10 == -Inf] <- 0


  df_fmt <- as.data.frame(df)
  if (normalization.method == "log10")    {df_fmt <- as.data.frame(df_log10)}
  if (normalization.method == "minmax")   {df_fmt <- as.data.frame(df_norm)}
  if (normalization.method == "meannorm") {df_fmt <- as.data.frame(df_mnorm)}

  ## Get Gene order
  gene_list     <- c()
  peptide_order <- c()

  # This funtions supports type = c("hclust","sumrpk")
  orderGenes <- function(df_fmt,type = "hclust") {
    gene_list <- c()

    if (type == "hclust"){
      heatmap_mtx   <- heatmap.2(as.matrix(df_fmt))
      all_peptides  <- row.names(df_fmt)
      peptide_order <- all_peptides[heatmap_mtx$rowInd]
      gene_list    <- unique(tstrsplit(peptide_order,"_")[[1]])
    }

    if (type == "sumrpk"){
      genes <- unique(tstrsplit(row.names(df_fmt),"_")[[1]])
      sumrpk_vec <- c()

      for (g in genes) {
        gene_df <- as.data.frame(df_fmt[grep(g,row.names(df_fmt)),target.samples])
        gene_sum_rpk <- sum(rowSums(gene_df))
        sumrpk_vec <- c(sumrpk_vec,gene_sum_rpk)
      }
      names(sumrpk_vec) <- genes
      sumrpk_vec_sort   <- sort(sumrpk_vec,decreasing = T)
      gene_list         <- names(sumrpk_vec_sort)
    }

    return(gene_list)
  }


  if(nrow(df_fmt) == 1){
    peptide_order <- row.names(df_fmt)
    gene_list     <- tstrsplit(peptide_order,"_")[[1]]
  } else{
    gene_list <- orderGenes(df_fmt,type = "sumrpk")

    # If a gene order is provided, use that instead
    if (length(gene.order)> 0){
      gene_list <- gene.order[gene.order %in% gene_list]
    }
  }

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
      par(mar=c(1,1,1,1))

      # hierarchical clustering of peptide IDs
      gene_mtx <- as.matrix(gene_df)
      d        <- dist(gene_mtx)
      hc       <- hclust(d)
      gene_peptides_order <- rownames(gene_mtx)[hc$order]

      # heatmap.2 way
      #gene_hp_mtx  <- heatmap.2(as.matrix(gene_df))
      #gene_peptides       <- row.names(gene_df)
      #gene_peptides_order <- gene_peptides[gene_hp_mtx$rowInd]
    }

    gene_list_ordered <- c(gene_list_ordered,gene_peptides_order)
  }

  ## Re-order list of genes in heatmap
  df_fmt <- df_fmt[gene_list_ordered,]

  # Re-order columns incase disease needs to be towards the right
  if(disease.first != TRUE){
    beads <- reference.samples.rem.order[grep(null.ip.col,reference.samples.rem.order)]
    ref_no_bead <- reference.samples.rem.order[! reference.samples.rem.order %in% beads]
    reordered_cols <- c(target.samples,ref_no_bead,disease.samples.rem.order,beads)

    df_fmt <- df_fmt[,reordered_cols]
  }


  ## Annotation column Samples ----
  sample_vec <- names(as.data.frame(df_fmt))
  sample_vec <- sample_vec[! sample_vec %in% c(gene.col,peptide.id.col)]
  diagnosis  <- sample_vec

  # Target Samples
  diagnosis[diagnosis %in% disease.samples.orig]    <- disease.annot
  diagnosis[diagnosis %in% reference.samples.orig]  <- reference.annot

  # Remaining samples
  diagnosis[diagnosis %in% disease.samples]    <- disease.annot
  diagnosis[diagnosis %in% reference.samples]  <- reference.annot

  diagnosis_annot <- as.data.frame(diagnosis)
  row.names(diagnosis_annot) <- sample_vec

  ## Annotation column Genes ----
  peptide_vec <- row.names(as.data.frame(df_fmt))
  genes_vec   <- tstrsplit(peptide_vec,"_")[[1]]
  genes_annot <- as.data.frame(genes_vec)
  row.names(genes_annot) <- peptide_vec

  ## Color Scheme

  break_list1 <- seq(-1,1,0.02)  # form minmax norm
  break_list2 <- seq(-2,2,0.04)  # -2 to 2
  break_list3 <- seq(-3,3,0.06)  # -3 to 3
  break_list4 <- seq(-4,4,0.08)  # -4 to 4


  # break list -1 to 4
  break_1 <- seq(-1,1,0.04)  # first 51 breaks
  break_2 <- seq(1,4,0.06)  # last 50 breaks
  break_2 <- break_2[3:length(break_2)]
  break_list1_4 <- c(break_1,break_2)

  #break_list1_4 <- break_list1[break_list1 < 0]
  #break_list1_4 <- c(break_list1_4,break_list4[break_list4 >= 0])

  break_list = NULL
  if (normalization.method == "log10")   {break_list = break_list1_4}
  if (normalization.method == "minmax")  {break_list = break_list1}
  if (normalization.method == "meannorm"){break_list = break_list1}

  if(!is.null(break.list)){
    if (break.list == 1){break_list = break_list1}
    if (break.list == 2){break_list = break_list2}
    if (break.list == 3){break_list = break_list3}
    if (break.list == 4){break_list = break_list4}
  }

  color_pal  <- bluered(100)

  ## Heatmap command
  if (gene.col %in% names(df_fmt)){df_fmt[,gene.col] <- NULL}

  # Column annotation bars colors
  colcolor        <- c("purple", "darkgreen")
  names(colcolor) <- c(disease.annot, reference.annot)


  # Row annotation gene block colors
  all_genes <- as.character(unique(genes_annot$genes_vec))
  rowcolor <- c()

  temp_col <- "white"
  for (gene in all_genes) {
    if(temp_col == "white"){
      temp_col = "black"
    } else if (temp_col == "black"){
      temp_col = "white"
    }
    rowcolor <- c(rowcolor,temp_col)
  }
  names(rowcolor) <- all_genes

  # Create annotation color list
  anno_colors <- list(diagnosis = colcolor,genes_vec = rowcolor)

  # Column Gaps
  col_gaps <- c()
  gap1     <- length(target.samples)
  col_gaps <- c(col_gaps,gap1)

  if(disease.first == TRUE){
    gap2 <- length(c(target.samples,disease.samples.rem.order))
    if (length(disease.samples.rem.order)>0){col_gaps <- c(col_gaps,gap2)}
  }
  if(disease.first != TRUE){
    gap2 <- length(c(target.samples,reference.samples.rem.order)) - num.pos.controls
    if (length(reference.samples.rem.order)>0){col_gaps <- c(col_gaps,gap2)}
  }

  gap3 <- length(names(df_fmt)) - num.pos.controls
  col_gaps <- c(col_gaps,gap3)

  # Custom Row names
  cust_labels_row <- row.names(df_fmt)
  cust_labels_row <- tstrsplit(cust_labels_row,"_")[[1]]

  for (gene in unique(cust_labels_row)) {
    gene_sub <- cust_labels_row[cust_labels_row %in% gene]
    if (length(gene_sub)> 1) {
      gene_sub[2:length(gene_sub)] <- ""
    }
    cust_labels_row[cust_labels_row %in% gene] <- gene_sub
  }

  label_row <- NULL
  if(custom.gene.labels == TRUE){label_row <- cust_labels_row}

  p <- pheatmap(df_fmt,cluster_cols = F, cluster_rows = F,
                color = color_pal,breaks = break_list,fontsize_row = 4.5,
                annotation_col = diagnosis_annot, annotation_row = genes_annot,
                annotation_legend = F,annotation_names_col = F, annotation_names_row = F,annotation_colors = anno_colors,
                gaps_col = col_gaps,cellheight = cell.height, cellwidth = cell.width,
                border_color = NA,labels_row = label_row)

  if (transpose.heatmap == TRUE){
    df_fmt_inv <- as.data.frame(t(df_fmt))
    gene_names <- tstrsplit(names(df_fmt_inv),"_")[[1]]
    p <- pheatmap(df_fmt_inv,cluster_cols = F, cluster_rows = F,
                  color = color_pal,breaks = break_list,fontsize_row = 4.5,
                  annotation_col = genes_annot, annotation_row = diagnosis_annot,
                  annotation_legend = F,annotation_names_col = F, annotation_names_row = F,annotation_colors = anno_colors,
                  gaps_row = c(gap1,gap2,gap3),cellheight = cell.height, cellwidth = cell.width,
                  border_color = NA,labels_col = gene_names)
  }

  return(p)
}



#' @title Average Replicates in Peptide data frame
#' @description takes input of phage results used for Heatmaps and averages the replicates on non=target columns
#' @param df phage peptide count matrix
#' @param target.samples vector of target samples of interest from one patient (defined by 'patient' input)
#' @param disease.samples vector of all disease samples in input data frame
#' @param reference.samples vector of all reference samples (healthy and bead controls) in input data frame
#' @param chop.length number characters to remove from the end of each sample name to remove the replicate IDs <default: 2>
#' @param split.by insert regex to split sample strings by and keep the first element, which helps simplify complex sample names <default: NULL>
#' @param null.ip.col option that indicates a regex to group all AG Bead columns <default: Bead >
#' @return data frame with non-target columns averaged
#' @examples avgReplicates(df,target.samples,disease.samples,reference.samples)
#' @import data.table
#' @export
avgReplicates <- function(df,target.samples,disease.samples,reference.samples,chop.length=2, null.ip.col = "Bead",split.by = NULL) {
  df_avg <- df[,target.samples]

  # Format new sample names
  other.samples     <- names(df)[! names(df) %in% target.samples]
  other.samples.fmt <- collapseReplicateColumnNames(other.samples,chop.length, null.ip.col=null.ip.col,split.by=split.by)

  for (sample in other.samples.fmt) {
    sample_reps       <- names(df)[grep(sample,names(df))]
    df_sub            <- as.data.frame(df[,sample_reps])
    names(df_sub)     <- sample_reps
    row.names(df_sub) <- row.names(df)
    if (ncol(df_sub) > 1){
      df_sub[,sample] <- rowMeans(df_sub)
      df_avg <- cbind(df_avg,df_sub[,sample])
    } else{
      df_avg <- cbind(df_avg,df_sub)
    }
  }

  names(df_avg) <- c(target.samples,other.samples.fmt)

  return(df_avg)

}


#' @title Format Column names to remove replicate IDs
#' @description take an input vector of sample names and collapse all replicates
#' @param col.name.vec input column names, sample names only
#' @param chop.length number of characters to chop from the end of each sample name <default: 2>
#' @param split.by insert regex to split sample strings by and keep the first element, which helps simplify complex sample names <default: NULL>
#' @param null.ip.col option that indicates a regex to group all AG Bead columns <default: Bead >
#' @return vector of unique collapsed sample names
#' @examples collapseReplicateColumnNames(col.name.vec,chop.length = 2)
#' @import data.table
collapseReplicateColumnNames <- function(col.name.vec,chop.length = 2,null.ip.col = "Bead",split.by = NULL) {
  for (i in 1:length(col.name.vec)) {
    sample               <- col.name.vec[i]
    if (grepl(null.ip.col,sample)){
      col.name.vec[i] <- tstrsplit(sample,"_")[[1]]
    } else if (grepl("GFAP_",sample)){
      col.name.vec[i] <- tstrsplit(sample,"_")[[1]]
    } else{
      if(!is.null(split.by)){
        sample <- tstrsplit(sample,split.by)[[1]]
      }
      col.name.vec[i] <- substr(sample,1,nchar(sample)-chop.length)
    }

  }
  col.name.vec <- unique(col.name.vec)

  return(col.name.vec)
}

