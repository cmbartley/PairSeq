# Search through a set of peptides for a single gene for Kmer overlaps
#  - If there is no overlap, return an empty data frame
findKmerOverlap <- function(peptide_df, KMER_SIZE = 7){
  # Initialize list of Kmers
  kmer_list         <- list()
  overlapping_kmers <- c()
  peptides_to_keep  <- c()
  
  
  for(i in 1:nrow(peptide_df)){
    row <- peptide_df[i,]
    pep <- row$peptide_id
    seq <- unlist(strsplit(row$sequence,""))
    
    end_pos <- length(seq) - KMER_SIZE + 1
    for(i in 1:end_pos){
      kmer_seq <- seq[i:(i+KMER_SIZE)]
      kmer_seq <- paste(kmer_seq,collapse = "")
      
      # Add kmer into KMER list with peptide ID
      if(kmer_seq %in% names(kmer_list)){
        # Add peptide ID to Kmer list
        current_peptides <- kmer_list[[kmer_seq]]
        all_peptides <- unique(c(current_peptides,pep))
        kmer_list[[kmer_seq]] <- all_peptides
        
        # Store all overlapping kmer sequences and peptide IDs
        overlapping_kmers <- c(overlapping_kmers,kmer_seq)
        peptides_to_keep <- c(peptides_to_keep,all_peptides)
      }else {
        kmer_list[[kmer_seq]] <- pep
      }
      
    }
  }
  peptides_to_keep  <- unique(peptides_to_keep)
  overlapping_kmers <- unique(overlapping_kmers)
  
  # KMER list is more of a sanity check
  kmer_list_overlap <- kmer_list[overlapping_kmers]
  
  peptide_df_clean <- peptide_df[peptide_df$peptide_id %in% peptides_to_keep,]
  
  return(peptide_df_clean)
}



fullParse <- function(df,list_of_samples, FC_THRESH1 = 10, FC_THRESH2 = 100, SUM_RPK_THRESH = 50){
        candidate_peptides <- c()

        for(samp in list_of_samples){
        
                samp_fc <- paste0(samp,"_FC")
                
                #create sample-specific dataframe of minimum rpK >= 2 and minimum FC >= 10
                df_2rpK_FC10 <- subset(df, df[,samp_fc] >= FC_THRESH1)
                
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
                    polygenes <- c(polygenes, gene)
                    
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

