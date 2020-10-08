#' @title Creates Z-score Columns
#' @description takes input of phage results and calculates Z-scores for every target column based on the input
#' @param df phage peptide count matrix
#' @param target.columns list of target samples to analyze
#' @param method method to use for calculating Z-scores ("median","mean") <default: "median">
#' @param MIN minimum count value for the inputted data frame <default: 1>
#' @param onlyZscores option to format the output and return full data frame or only the Z-score columns <default: FALSE>
#' @return returns phage results with FC columns for each target sample.
#' @examples calcZscores(df, target_columns)
#' @export

## The is code takes a dataframe and calculates the fold change for each column relative to the mean of the control samples
## parameters will need to be adjusted  for other dataframes

calcZscores <- function(df, target.samples, method = "median", MIN = 1,peptide.id.col = "peptide_id",onlyZscores = FALSE){
  peptide_df <- df[,c(target.samples, peptide.id.col)]


  for (sample in target.samples) {
    sample_df <- as.data.frame(peptide_df[,c(sample,peptide.id.col)])

    # Filter based on MIN
    sample_df_sub <- sample_df[sample_df[,sample] > MIN,]

    # Z-score calculation
    m <- 0
    if (method == "mean")   {m <- mean(sample_df_sub[,sample])}
    if (method == "median") {m <- median(sample_df_sub[,sample])}
    sd(sample_df_sub[,sample])





  }

  df_out <- cbind(df, FC_df)
  if (onlyFC == TRUE){
    df_out <- FC_df
  }


  #scatter y-pos over plot
  bead_df_log <- log10(peptide_df$CoV_NID1567_Plasma_1)
  bead_df_log[bead_df_log == -Inf] <- 0
  d <- density(bead_df_log)
  plot(d, main="Test")
  polygon(d, col="red", border="blue")

  return(df_out)
}



