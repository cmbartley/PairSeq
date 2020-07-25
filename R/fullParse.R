#' @title Filters phage data
#' @description takes input of phage results and applies RPK and Fold Change filters
#' then returns a subset of phage results
#' @param df phage peptide count matrix
#' @param list_of_samples list of target samples to analyze
#' @param MIN_RPK Threshold for minimum RPK value <default: 0>
#' @param FC_THRESH1 Fold Change threshold 1 <default: 10>
#' @param FC_THRESH2 Fold Change threshold 2 <default: 100>
#' @param SUM_RPK_THRESH Gene level total RPK threshold <default: 50>
#' @return returns a subset of phage results with peptides that pass all filters.
#' @examples fullParse(df,list_of_samples, MIN_RPK = 2,FC_THRESH1 = 10, FC_THRESH2 = 100, SUM_RPK_THRESH = 50)
#' @export

fullParse <- function(df,list_of_samples, MIN_RPK = 0,FC_THRESH1 = 10, FC_THRESH2 = 100, SUM_RPK_THRESH = 50){
        candidate_peptides <- c()

        for(samp in list_of_samples){

                samp_fc <- paste0(samp,"_FC")

                #create sample-specific dataframe of minimum rpK >= 2 and minimum FC >= 10
                df_2rpK      <- subset(df, df[,samp] >= MIN_RPK)
                df_2rpK_FC10 <- subset(df_2rpK, df_2rpK[,samp_fc] >= FC_THRESH1)

                #create separate sample-specific dataframe of minimum rpK >= 2 and minimum FC >= 100
                df_2rpK_FC10_FC100 <- subset(df_2rpK_FC10, df_2rpK_FC10[,samp_fc] >= FC_THRESH2)

                #expand FC100 dataframe to include FC >= 10 phage that mapt to same genes as FC >= 100 phage (goal to identifying overlapping, enriched phage)
                df_extended <- df_2rpK_FC10[which(df_2rpK_FC10$gene %in% df_2rpK_FC10_FC100$gene),]

                #create unique list of genes in the expanded df
                genes <- unique(df_extended$gene)

                #filter out genes whose summed rpK < 50
                for(gene in genes){
                  df_get_50rpK <- df_extended[df_extended$gene %in% gene,]

                  # Verify that sum RPK of genes meets threshold
                  if (sum(df_get_50rpK[,samp]) >= SUM_RPK_THRESH){
                    # Verify at least one peptide has a FC greater than the 2nd threshold
                    if (nrow(df_get_50rpK) > 1){
                      peptide_fc_vec <- df_get_50rpK[,samp_fc]
                      sig_fc         <- peptide_fc_vec[peptide_fc_vec > FC_THRESH2]

                      ## Significant peptide candidates
                      if (length(sig_fc) > 0){
                        # Store peptide IDs
                        candidate_peptides <- c(candidate_peptides,unique(df_get_50rpK$peptide_id))
                      }
                    }

                  }else{
                    next
                  }
                }

        }

        # Subset Significant peptides
        candidate_df <- df[df$peptide_id %in% candidate_peptides,]

        return(candidate_df)
}

