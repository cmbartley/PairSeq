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
