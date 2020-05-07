#Purpose: Create new files, then write them into new folders automatically.
  ##Part 1: New sample names, add columns for calculated SUM, MEAN, FC
      # 1 new file 
  ##Part 2: Merge files by "peptide" (Part 1 w/ peptide_gene_mapping) 
      # 1 new file
  ##Part 3: Filter by rpk and FC thresholds, then create individual sample files
      # 96 new files
  ##Part 4: Pair CSF and plasma by "indicator" in file name
      # 6 new files

## Notes for improvement: create function, generalize columns
# Packages required:
install.packages("rstudioapi") 
install.packages("svDialogs")

# Input for a future function: Choose input folder, output folder, choose file and indicate controls,

############ Part 0: Setup - Assign wd variables, necessary for downstream "dest" variables to create folders and paths.

{
  library(rstudioapi)
  library("svDialogs")
# Select directory for output 
  message_outputdir <- dlg_message("Select/create a folder to save output files")
  phdir1 <-selectDirectory()
  
  message_inputdir <- dlg_message("Select a folder with input dataset and peptide-gene map (Both files must be in the same folder))")
# Select directory with all relevant datasets
  master_dir <-selectDirectory()
  
  message_choosefile <- dlg_message("Select data set")
  phagefile <- file.choose()        # Choose phage dataset for analysis
  numctrl <- 12                     #Indicate number of non-samples (i.e. beads, ctrls)
  posctrl <- 8                      #Indicate number of positive controls
  folder_name <- paste0(sub(".csv", "", (gsub("^.*/","", phagefile)),))
                                          #Create name of folder
  filename <- paste0(folder_name, "_MEAN_SUM_FC.csv", sep = "")
                                          #Create name of part 1's output file
}


#Press RUN here

{
###Start of RUN
{
  ############ Part 1: Create new file, add SUM, MEAN, FC columns
  csv <- read.csv(phagefile)
          
        
#Create an output folder, and create new path to it.
  dest0 <- paste0(phdir1, "/", folder_name)
  dir.create(dest0)
  path_out = dest0
        
#Rearrange bead columns, move BEADS to end for next steps.
  FTD_nobeads <- c(3:(ncol(csv)-numctrl), ((ncol(csv)-(posctrl)+1):ncol(csv)))
      #subset according to numctrl and posctrl
  columns <- c(FTD_nobeads, 88:91)
  input <- csv[ ,columns]
      #completed rearrangement
        
#Assign columns (by number) to variables
  FTD_plasma <- 2:58
  FTD_CSF <-59:85
  FTD_posctrl <- 86:93          # May be redundant?
  FTD_ALL <- c(FTD_plasma ,FTD_CSF, FTD_posctrl)             # all samples and ctrls
  FTD_BEADS <- 94:(ncol(csv)-2)    # beads only
    
#Rename column names: removing last 12 characters
  posctrl_col <- colnames(input)[FTD_posctrl]
  posctrl_name <- names(input)[FTD_posctrl]
  posctrl_name <-substr(posctrl_name, 1, nchar((posctrl_col))-12)
  
  beads_col <- colnames(input)[FTD_BEADS]
  beads_name <- names(input)[FTD_BEADS]
  beads_name <-substr(beads_name, 1, nchar((beads_col))-12)

#Add column: sum of all rows
  input$SUM_ALL <- rowSums(input[ ,(2:ncol(input))])  

#Add column: mean of each group
  input$mean_FTD <- rowMeans(input[,FTD_ALL])       #all FTD samples
  input$mean_BEADS <- rowMeans(input[,FTD_BEADS])   #BEADS only
  input$mean_PLASMA <- rowMeans(input[,FTD_plasma]) #plasma
  input$mean_CSF <- rowMeans(input[,FTD_CSF])       #CSF


#Fold-change dataframe - value divided by mean of other
##Normalize to other group, vice-versa
    FC_FTD <- (input[,FTD_ALL])/ input$mean_BEADS #normalize to beads
    FC_BEADS <- (input[,FTD_BEADS])/ input$mean_FTD #relative to samples

#Fold changes merged with input data frame
        FCbind <- cbind(input, FC_FTD, FC_BEADS)

#Rename FC columns
        colnames(FCbind)[103:198] <- paste("FC_", colnames(FCbind[,c(103:198)]), sep = "")
        #instead of 103:198 -> ncol(input + 1):ncol(FCbind)

#Export as .csv
        write.csv(FCbind, file.path(path_out, filename)) #save as .csv


}

#Clear large csv
  remove(csv,FC_FTD,FC_BEADS,input, FCbind)
  

############ Part 2: Merge Files

{
#Set wd and read newly created file
  phdir2 <- folder_name
  setwd(phdir1)
  setwd(phdir2) #New folder as wd
FTDC <- read.csv(filename) #read output of Part 1
  setwd(master_dir) #Revert wd to master folder
peptide_csv <- "peptide_gene_mapping_TN.csv" 
peptide <- read.csv(peptide_csv) #read peptide file
        #Duplicate file
#Create output file
  subfolder_name <- paste(Sys.Date(), "_Single_Files_", gsub("_.*", "", folder_name), sep = "")
  dest1 <- paste0(dest0, "/",subfolder_name)
  dir.create(dest1)
  path_out2 = dest1

#Rename peptide_gene_mapping seq_id to peptide
  names(peptide)[names(peptide) == "seq_id"] <- "peptide"

#Merge files, create gene&index column, create .csv
  initial_merge <- merge(peptide, FTDC, by = "peptide")
  peptide_id <- paste0(initial_merge$gene, "_", initial_merge$index)
  merged <- cbind(peptide_id, initial_merge)

#mergefilename <- "merge_FTD_Cohort_peptide_2.csv"
  mergefilename <- paste("MERGED_", gsub(".csv", "_", filename),  gsub(".csv", "", peptide_csv), ".csv",sep = "")
  write.csv(merged,file.path(path_out2, mergefilename)) #save as .csv

#Clear large csv
  remove(FTDC, peptide, peptide_id)
}

############ Part 3: Filter by rpk and FC thresholds
{
#Create output folder
  folder_name4 <- paste0("Filtered_Files_", sub("*_._.*.*", "",  mergefilename))
    print(folder_name4)
  dest2 <- paste0(dest1, "/", folder_name4)
  dir.create(dest2)
  path_out3 <- dest2


#Function to create Excel Column format (for output filename, i.e "A", "H", "BO")
letterwrap <- function(n, depth = 1) {
        args <- lapply(1:depth, FUN = function(x) return(LETTERS))
        x <- do.call(expand.grid, args = list(args, stringsAsFactors = F))
        x <- x[, rev(names(x)), drop = F]
        x <- do.call(paste0, x)
        if (n <= length(x)) return(x[1:n])
        return(c(x, letterwrap(n - length(x), depth = depth + 1)))
        }
excel_columns <- letterwrap(26*27) #A to ZZ
}

  
#Loop to run rpk and FC filters for every sample
{
for (i in (seq(9, 104))) {
        y <- 101
        
        which_df <- merged[(merged[,i] >= 2) & (merged[,i + y] >= 10), ]
        FC10_df <- which_df[which(which_df[,i] >=2 & which_df[,i+y] >= 10), ]
        FC100_df <- FC10_df[which(FC10_df[,i +y] >= 100), ]
        
        #Gather peptides from 100-fold enrichment csv for 10-fold threshold
        peptide_check <- FC10_df[which(FC10_df$gene %in% FC100_df$gene), ]
        
        #Filter any gene with total rpk >= 50
        sum_rpk50 <- peptide_check[peptide_check[,i] >= 50, ]
        sum_rpk50_order <- sum_rpk50[order((sum_rpk50)[1], decreasing = F),]
        #Alphabetical order for peptide_id
#Write csv        
write.csv(sum_rpk50_order, file.path(path_out3, 
  paste0(i, "_Col_",excel_columns[i+1],"__",(colnames(sum_rpk50_order[i])),"_", "Filtered", ".csv"))) 
        }
}
  

############Part 4: Pair CSF and plasma
#Create folder for output file. 
{
folder_name5 <- "Paired_CSF_Plasma"
dest3 <- paste0(dest2, "/", folder_name5)
dir.create(dest3)
path_out5 <- dest3
setwd(dest2)
#Function to pair CSF and plasma by "indicator" in name
sample_pair<- function(array_CSF, array_plasma) {
  array_CSF2 <- paste0(array_CSF, "_",sep = "")
  array_plasma2 <- paste0(array_plasma, "_", sep = "")
  
  #Look for files containing indicator.
  col1 <- list.files(pattern=array_CSF2)
  col2 <- list.files(pattern=array_plasma2)
  
  CSF <- read.csv(col1)
  plasma <- read.csv(col2)
  pair <- rbind(CSF,plasma)
  write.csv(pair, file.path(path_out5, 
                            paste0(sub("_CSF", "_paired_samples", col1))))
  
}
#Indicator used = phage plate position
sample_pair("A1", "A2")
sample_pair("A3", "A4")
sample_pair("A5", "A6")
sample_pair("A7", "A8")
sample_pair("A9", "A10")
sample_pair("A11", "A12")
}
####END of RUN
}

