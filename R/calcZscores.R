#' @title Creates Z-score Columns
#' @description takes input of phage results and calculates Z-scores for every target column based on the input
#' @param df phage peptide count matrix
#' @param target.samples vector of target samples to analyze
#' @param reference.samples vector of reference control samples used for calculating Std. Dev. and average expression
#' @param method method to use for calculating Z-scores ("median","mean") <default: "mean">
#' @param MIN minimum count value for the inputted data frame <default: 1>
#' @return returns phage results with FC columns for each target sample.
#' @examples calcZscores(df, target_columns)
#' @export

## The is code takes a dataframe and calculates the fold change for each column relative to the mean of the control samples
## parameters will need to be adjusted  for other dataframes

calcZscores <- function(df, target.samples, reference.samples, method = "mean", MIN = 1){
  peptide_df <- df[,unique(c(target.samples,reference.samples))]

  # calc Std dev per row
  ref_df <- peptide_df[,unique(c(reference.samples))]
  ref_df <- transform(ref_df, SD=apply(ref_df,1, sd, na.rm = F))
  peptide_df$SD <- ref_df$SD

  # Calc mean or median
  ref_df$SD <- NULL
  ref_df <- transform(ref_df, AVG=apply(ref_df,1, method, na.rm = F))
  peptide_df$AVG <- ref_df$AVG

  zscore_df  <- data.frame(matrix(nrow = nrow(peptide_df), ncol = 0))

  for (sample in target.samples) {
    sample_vec <- peptide_df[,sample]
    zscore_vec <- sample_vec - peptide_df$AVG
    zscore_vec <- zscore_vec / peptide_df$SD
    zscore_vec[is.na(zscore_vec)] <- 0

    zscore_df <- cbind(zscore_df,as.data.frame(zscore_vec))
  }

  names(zscore_df) <- target.samples

  return(zscore_df)
}
