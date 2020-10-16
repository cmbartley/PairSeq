#' @title Identify Kmer overlaps
#' @description takes input of phage results and finds all Kmer overlaps among all peptide sequences
#' then returns a subset of phage results with peptides that contain overlaps
#' @param peptide_df phage peptide count matrix
#' @param KMER_SIZE size of Kmer <default: 7>
#' @param TEMP_FASTA file handle for a temporary FASTA that needs to be generated in order to perform the Kmer analysis <default: "./temp_pep.fasta">
#' @param seq.col column name for peptide AA sequences in input data frame <default: "sequence">
#' @param peptide.col column name for peptide IDs in input data frame <default: "peptide_id">
#' @return returns a subset of phage results with peptides that contain Kmer overlaps.
#' @examples findKmerOverlap(peptide_df,KMER_SIZE = 7,seq.col = "sequence",peptide.col = "peptide_id")
#' @import seqinr
#' @import kmer
#' @import ape
#' @export

findKmerOverlap <- function(peptide_df, KMER_SIZE = 7, TEMP_FASTA = "./temp_pep.fasta",seq.col = "sequence",peptide.col = "peptide_id"){
  seq_list     <- as.list(peptide_df[,seq.col])
  peptide_list <- as.list(peptide_df[,peptide.col])

  write.fasta(sequences = seq_list,names = peptide_list,file.out = TEMP_FASTA)

  # Read in FASTA
  mySequences      <- read.FASTA(TEMP_FASTA,type = "AA")
  file.remove(TEMP_FASTA)

  # Kcount of 7mers
  kmer_count_table <- kcount(mySequences,k = KMER_SIZE)

  # Identify overlaps
  pep_to_keep      <- list()
  for (pep in row.names(kmer_count_table)) {
    kmer_df <- kmer_count_table[row.names(kmer_count_table) %in% pep,]
    target_kmers <- names(kmer_df[kmer_df > 0])

    kmer_df_sub <- kmer_count_table[,colnames(kmer_count_table) %in% target_kmers]
    kmer_df_sub <- as.data.frame(kmer_df_sub)
    # Change the summing somehow
    kmer_df_sub$total <- rowSums(kmer_df_sub)
    kmer_df_sub       <- kmer_df_sub[kmer_df_sub$total > 0,]

    if (nrow(kmer_df_sub) > 1){
      total_overlapping_pep <- length(unique(row.names(kmer_df_sub))) -1
      pep_to_keep[[pep]]    <- total_overlapping_pep
    }
  }

  peptide_df_clean <- peptide_df[peptide_df[,peptide.col] %in% names(pep_to_keep),]

  return(peptide_df_clean)
}



#' @title Identify Kmer overlaps
#' @description takes input of phage results and finds all Kmer overlaps among all peptide sequences
#' then returns a subset of phage results with peptides that contain overlaps
#' @param peptide_df phage peptide count matrix
#' @param KMER_SIZE size of Kmer <default: 7>
#' @return returns phage results with a 'KmerOverlap' column saying whether peptides that contain Kmer overlaps with a 'yes' or 'no'.
#' @examples findKmerOverlap.old(peptide_df,KMER_SIZE = 7)
#' @export

findKmerOverlap.old <- function(peptide_df, KMER_SIZE = 7,peptide.id.col = "peptide_id",sequence.col = "sequence"){
  # Initialize list of Kmers
  kmer_list          <- list()
  overlap_count_list <- list()
  overlapping_kmers  <- c()
  peptides_to_keep   <- c()


  for(i in 1:nrow(peptide_df)){
    row <- peptide_df[i,]
    pep <- row[,peptide.id.col]
    seq <- unlist(strsplit(row[,sequence.col],""))

    end_pos <- length(seq) - KMER_SIZE + 1
    for(i in 1:end_pos){
      kmer_seq <- seq[i:(i+KMER_SIZE -1)]
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

  # Add column for Kmer overlap
  kmer_overlap_col <- peptide_df[,peptide.id.col]

  # Determine number of Kmer overlaps
  for (kmer in names(kmer_list_overlap)){
    peptide_vec <- kmer_list_overlap[[kmer]]
    for (pep in peptide_vec){
      if (pep %in% names(overlap_count_list)){
        current                   <- overlap_count_list[[pep]]
        overlap_count_list[[pep]] <- unique(c(current,kmer))
      } else{
        overlap_count_list[[pep]] <- kmer
      }
    }
  }
  overlap_count_list <- unlist(lapply(overlap_count_list, length))

  kmer_overlap_col[! kmer_overlap_col %in% names(overlap_count_list)] <- 0
  kmer_overlap_col[kmer_overlap_col %in% names(overlap_count_list)] <- overlap_count_list


  peptide_df$KmerOverlap <- as.numeric(kmer_overlap_col)

  #peptide_df_clean <- peptide_df[peptide_df$peptide_id %in% peptides_to_keep,]

  return(peptide_df)
}
