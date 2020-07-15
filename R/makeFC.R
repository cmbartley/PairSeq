## The is code takes a dataframe and calculates the fold change for each column relative to the mean of the control samples
## parameters will need to be adjusted  for other dataframes

makeFCdf <- function(df, target_columns, mean_column){
  FC_colnames <- paste0(target_columns,"_FC")
  
  FC_df <- df[,target_columns] / df[,mean_column]
  
  colnames(FC_df) <- FC_colnames
  
  df_out <<- cbind(df, FC_df)
  
  return(df_out)
}
