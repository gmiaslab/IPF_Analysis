## Example STAR count-column selection helper.
##
## This code was extracted from the upstream per-study count-assembly scripts
## to show how STAR ReadsPerGene/PerGene count files were converted into
## sample count columns by choosing the unstranded, sense, or antisense column
## based on the relative abundance of sense and antisense counts.
##
## Expected input files are STAR gene-count tables with columns corresponding
## to Ensembl gene ID, Unstranded, Sense, and Antisense counts.
## If a stranded column accounted for more than 80% of stranded gene-assigned
## counts, that column was used. Otherwise, unstranded counts were used.

###### Functions #####
accession_names <- function(fileName) {
  gsub("_GenesReadsPerGene[.]out[.]tab", "", basename(fileName))
}
readFile <- function(fileName) {
  read_tsv(fileName, skip = 4, col_names = c("EnsID", "Unstranded", "Sense", "Antisense"))
}
strandedness <- function(file, threshold) {
  sum_unstranded <- sum(file$Unstranded)
  sum_sense <- sum(file$Sense)
  sum_antisense <- sum(file$Antisense)
  if (sum_sense + sum_antisense == 0) {
    ratio_sense <- 0
    ratio_antisense <- 0
  } else {
    ratio_sense <- sum_sense/(sum_sense + sum_antisense)
    ratio_antisense <- sum_antisense/(sum_sense + sum_antisense)
  }
  if (ratio_sense > threshold && ratio_antisense <= threshold) {
    use_column <- file$Sense
    print("Using sense")
  } else if (ratio_antisense > threshold && ratio_sense <= threshold) {
    use_column <- file$Antisense
    print("Using antisense")
  } else {
    use_column <- file$Unstranded
    print("Using unstranded")
  }
  return(data.frame(counts = use_column, row.names = file$"EnsID"))
}
# function to processes a study, determine strandedness and create dataframe
# for the study
# assuming all files ran the same accessions
 # skip any file with strings in the skip parameter
counts_generator <- function(accession,base_directory = getwd(),skip=c()) {
  print(accession)
  # setwd(paste(base_directory, accession, "/Gene Counts", sep = ""))
  pattern <- paste0(accession, ".*.PerGene.out.tab$")
  local_list <- list.files(path = paste(base_directory, sep = ""), pattern = pattern, full.names = TRUE)
  # Exclude files that match any of the strings in the skip vector
  if (length(skip) > 0) {
    exclude_pattern <- paste(skip, collapse = "|")
    local_list <- local_list[!grepl(exclude_pattern, local_list)]
  }
  print(local_list)
  local_tables <- lapply(local_list, readFile)
  local_expression <- lapply(local_tables, strandedness, threshold = 0.8)
  local_counts <- as.data.frame(local_expression)
  local_names <- unlist(lapply(local_list, accession_names))
  colnames(local_counts) <- c(local_names)
  return(local_counts)
}
merge_rownames <- function(x,y) {
  merged <- merge(x,y,by=0, sort=F)
  row.names(merged) <- merged$`Row.names`
  dplyr::select(merged,-Row.names)
}
