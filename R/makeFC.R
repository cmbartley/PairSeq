#' @title Creates Fold Change Columns
#' @description takes input of phage results and calculates Fold Change based on input columns
#' @param df phage peptide count matrix
#' @param target_columns list of target samples to analyze
#' @param mean_column column containing the mean RPK of samples that represent the controls for normalization
#' @return returns phage results with FC columns for each target sample.
#' @examples makeFCdf(df, target_columns, mean_column)
#' @export

## The is code takes a dataframe and calculates the fold change for each column relative to the mean of the control samples
## parameters will need to be adjusted  for other dataframes

makeFCdf <- function(df, target_columns, mean_column){
  FC_colnames <- paste0(target_columns,"_FC")

  FC_df <- df[,target_columns] / df[,mean_column]

  colnames(FC_df) <- FC_colnames

  df_out <<- cbind(df, FC_df)

  return(df_out)
}
