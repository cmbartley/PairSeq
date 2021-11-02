#' @title Filters phage data
#' @description takes input of phage results and applies RPK and Fold Change filters
#' then returns a subset of phage results
#' @param df phage peptide count matrix
#' @param list_of_samples list of target samples to analyze
#' @param MIN_RPK Threshold for minimum RPK value <default: 0>
#' @param FC_THRESH1 Fold Change threshold 1 <default: 10>
#' @param FC_THRESH2 Fold Change threshold 2 <default: 100>
#' @param SUM_RPK_THRESH Gene level total RPK threshold <default: 50>
#' @param ZSCORE_THRESH1 If calculated, applies minimum Z-score threshold for each peptide <default: 0>
#' @param ZSCORE_THRESH2 If calculated, applies additional Z-score threshold similar to the FC_THRESH2 arguement <default: 0>
#' @return returns a subset of phage results with peptides that pass all filters.
#' @examples fullParse(df,list_of_samples, MIN_RPK = 2,FC_THRESH1 = 10, FC_THRESH2 = 100, SUM_RPK_THRESH = 50)
#' @import seqinr
#' @import data.table
#' @export

fullParse <- function(df,list_of_samples, MIN_RPK = 0,FC_THRESH1 = 10, FC_THRESH2 = 100, SUM_RPK_THRESH = 50, ZSCORE_THRESH1 = 0, ZSCORE_THRESH2 = 0){
        candidate_peptides <- c()

        for(samp in list_of_samples){

                samp_fc <- paste0(samp,"_FC")
                samp_z  <- paste0(samp,"_Z")

                #create sample-specific dataframe of minimum rpK >= 2,  minimum FC >= 10 and minium Z-score threshold
                df_2rpK      <- subset(df, df[,samp] >= MIN_RPK)
                if(ZSCORE_THRESH1 != 0){
                  df_2rpk      <- subset(df_2rpK, df_2rpK[,samp_z] >= ZSCORE_THRESH1)
                }
                df_2rpK_FC10 <- subset(df_2rpK, df_2rpK[,samp_fc] >= FC_THRESH1)

                #create separate sample-specific dataframe of minimum rpK >= 2 and minimum FC >= 100
                df_2rpK_FC10_FC100 <- subset(df_2rpK_FC10, df_2rpK_FC10[,samp_fc] >= FC_THRESH2)
                if(ZSCORE_THRESH2 != 0){
                  ddf_2rpK_FC10_FC100 <- subset(df_2rpK_FC10_FC100, df_2rpK_FC10_FC100[,samp_z] >= ZSCORE_THRESH2)
                }

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
                      sig_fc         <- peptide_fc_vec[peptide_fc_vec >= FC_THRESH2]
                      if (ZSCORE_THRESH2 != 0){
                        peptide_fc_vec <- df_get_50rpK[,samp_z]
                        sig_fc         <- peptide_fc_vec[peptide_fc_vec >= ZSCORE_THRESH2]
                      }

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



#' @title Simple Filters for phage data
#' @description takes input of phage results and applies RPK and Fold Change filters
#' then returns a subset of phage results. This is more simple version of 'fullParse'
#' @param peptide_df phage peptide count matrix
#' @param list_of_samples list of all samples. It is assumed the FC columns are sample names with '_FC' at the end
#' @param RPK_THRESH RPK threshold <default: 20>
#' @param FC_THRESH Fold Change threshold <default: 100>
#' @param ZSCORE_THRESH If calculated, applies minimum Z-score threshold for each peptide <default: 0>
#' @param peptide.id.col Name of column with unique Peptide ID <default: "peptide_id">
#' @return returns a subset of phage results with peptides that pass all filters.
#' @examples filterPeptides(peptide_df,list_of_samples, RPK_THRESH = 20,FC_THRESH = 100,peptide.id.col = "peptide_id")
#' @export

filterPeptides <- function(peptide_df, list_of_samples,RPK_THRESH = 20, FC_THRESH = 100, ZSCORE_THRESH = 0,peptide.id.col = "peptide_id"){
  list_of_samples  <- as.list(list_of_samples)

  filterSample <- function(samp, df,RPK_THRESH = 20, FC_THRESH = 100, ZSCORE_THRESH = 0,peptide.id.col = "peptide_id"){
    samp_fc <- paste0(samp,"_FC")
    samp_z  <- paste0(samp,"_Z")
    df_sub  <- subset(df, df[,samp] >= RPK_THRESH)
    df_sub  <- subset(df_sub, df_sub[,samp_fc] >= FC_THRESH)
    if(ZSCORE_THRESH != 0){
      df_sub  <- subset(df_sub, df_sub[,samp_z] >= ZSCORE_THRESH)
    }
    pep_list <- df_sub[,peptide.id.col]

    return(pep_list)
  }

  peptides_to_keep <- unlist(lapply(list_of_samples, filterSample,
                                    df=peptide_df,RPK_THRESH=RPK_THRESH,FC_THRESH=FC_THRESH,ZSCORE_THRESH=ZSCORE_THRESH,
                                    peptide.id.col=peptide.id.col))

  peptide_df_sub <- peptide_df[peptide_df[,peptide.id.col] %in% peptides_to_keep,]

  return(peptide_df_sub)
}
