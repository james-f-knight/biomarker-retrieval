args <- commandArgs()
################ IMPORT FILES ############

cat("Importing files...\n")
#Import tilling array expressions and precomputed regulations
tiling_array_expressions <- read.delim("tiling_array_expressions.txt", 
                                       stringsAsFactors = FALSE)
colnames(tiling_array_expressions)[2] <- "Locus"
tiling_array_expressions <- tiling_array_expressions[,-c(3:12)]

tiling_array_regulations <- read.delim("tiling_array_regulations.txt", 
                                       stringsAsFactors = FALSE)
colnames(tiling_array_regulations)[2] <- "Locus"
tiling_array_regulations <- tiling_array_regulations[,-c(3:12)]

# Import RegNet matrix
reg_net_matrix <- read.table("gene_regulations.gb",
                             sep="\t", header = TRUE,
                             colClasses = "character",
                             stringsAsFactors = FALSE)
reg_net_table <- unique(reg_net_matrix[,1, drop = FALSE])
colnames(reg_net_table)[1] <- "Locus" 

#Create vector with positions of stress samples
stress_args <- unlist(strsplit(args[6], '-'))
stress_pos <- NULL
for(arg in stress_args)
{
  if(grepl(':', arg))
  {  
  pos <- unlist(strsplit(arg, ':'))
  stress_pos <- c(stress_pos, pos[1]:pos[2])
  }
  else
  {
    stress_pos <- c(stress_pos, as.numeric(arg))
  }
}

# Use only tiling array data if no RNA-Seq tables are provided
if(!dir.exists("CountTables"))
{
  regulator_data_table <- as.data.frame(merge(reg_net_table, tiling_array_regulations))
  #regulator_data_table <- as.data.frame(tiling_array_regulations)
  stress_row <- c(1,2, lapply(1:(ncol(regulator_data_table) -2), function(n) 
                                                                 {
                                                                   if(n %in% stress_pos) return("STRESS")
                                                                   return("CONTROL")
                                                                 }))
  finished_data_table <- rbind(regulator_data_table, stress_row)
  write.table(finished_data_table, 
              file = "Results/RNAseq_and_tiling_preprocessed",
              quote = FALSE,
              row.names = FALSE,
              sep = "\t")
  
  config <- read.table("Results/configuration.conf", stringsAsFactors = FALSE)
  config["misclassification_cost:",] <- paste(length(unique(stress_pos)), 
                                              (length(tiling_array_regulations) - 2 - length(unique(stress_pos))), sep = ",")
  colnames(config)[1] <- "[parameters]"
  write.table(config, "Results/configuration.conf", row.names = TRUE, quote = FALSE)
  quit()
}


library(stringi)
#Import RNA-Seq count tables

rnaseq_sample_names <- list.files("CountTables/")


rnaseq_split <- sapply(stri_split_fixed(rnaseq_sample_names, '-', n=2), unlist)
rnaseq_conditions <- rnaseq_split[1,]

rnaseq_sample_files <- paste("CountTables/", 
                       rnaseq_sample_names, sep = "")

rnaseq_count_data <- sapply(rnaseq_sample_files, 
                            function (file) 
                            {
                              table<-read.table(file, sep="\t", header = FALSE)
                              return(table[,2])
                            })

colnames(rnaseq_count_data) <- rnaseq_split[2,]

row.names(rnaseq_count_data) <- read.table(rnaseq_sample_files[1], sep="\t", 
                                           header = FALSE,
                                           stringsAsFactors = FALSE)[,1]
rnaseq_count_data <- cbind.data.frame(row.names(rnaseq_count_data), rnaseq_count_data)
colnames(rnaseq_count_data)[1] <- "Locus"
for(name in rnaseq_sample_files) 
{
  cat(name)
  cat('\n')
}

cat("Files loaded\n")
cat("Normalizing RNA-Seq data...\n")

################ MERGE TILING ARRAY EXPRESSIONS AND RNA-SEQ COUNT TABLES ############

merged_data_table <- merge(tiling_array_expressions, rnaseq_count_data, BY="Locus")

################ NORMALIZE RNASEQ COUNTS AND AVERAGE TECHNICAL REPLICATES ############
library(edgeR)

rnaseq_data <- merged_data_table[,(ncol(merged_data_table)-length(rnaseq_sample_files) + 1):ncol(merged_data_table)]
rnaseq_data <- data.frame(cpm(rnaseq_data))
array_data  <- merged_data_table[,3:(ncol(merged_data_table) - length(rnaseq_sample_files))]

duplicates <- which(duplicated(rnaseq_conditions))
replicates <- which(substr(rnaseq_conditions, "4", "4") == "T")
control_pos <- na.omit(match(which(substr(rnaseq_conditions, "3", "3") == "C"), which(!duplicated(rnaseq_conditions))))

result <- NULL
current <- NULL
name <- ""
count <- 1
for(n in 1:length(rnaseq_conditions))
{
  if(n %in% replicates)
  {
    if(n %in% duplicates)
    {
      current <- current + rnaseq_data[,n]
      count <- count + 1
    }
    else
    {
      if(length(current) > 0)
      {
        result <- cbind(result, current/count)
        colnames(result)[ncol(result)] <- name
      }
      current <- rnaseq_data[,n]
      count <- 1
      name <- rnaseq_conditions[n]
    }
  }
  else
  {
    result <- cbind(result, rnaseq_data[,n])
    colnames(result)[ncol(result)] <- rnaseq_conditions[n]
  }
}

if(length(current) > 0)
{
  result <- cbind(result, current/count)
  colnames(result)[ncol(result)] <- name
}

rnaseq_data <- data.frame(result)



################ NORMALIZE TILING ARRAY AND RNA-SEQ DATA ############

library(TDM)

normalized_rnaseq_data <- tdm_transform(ref_data = data.table(cbind(gene=rownames(array_data), array_data)),
                                        target_data = data.table(cbind(gene=rownames(rnaseq_data), rnaseq_data)))

normalized_rnaseq_data <- as.data.frame(lapply(normalized_rnaseq_data[,!1], as.numeric))

normalized_data_table <- cbind(merged_data_table[,1:(ncol(merged_data_table) - length(rnaseq_sample_files))], normalized_rnaseq_data)

write.table(normalized_data_table, 
            file = "Results/RNAseq_and_tiling_normalized",
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

################ CALCULATE REGULATIONS OF RNA-SEQ DATA ############

current <- NULL
rnaseq_regulations <- NULL
for(n in 1:ncol(normalized_rnaseq_data))
{
  if(n %in% control_pos)
  {
    current <- normalized_rnaseq_data[,n]
  }
  else
  {
    rnaseq_regulations <- cbind(rnaseq_regulations, normalized_rnaseq_data[,n] - current)
    colnames(rnaseq_regulations)[ncol(rnaseq_regulations)] <- colnames(normalized_rnaseq_data)[n] 
  }
}

rnaseq_regulations <- cbind(merged_data_table$Locus, rnaseq_regulations)
colnames(rnaseq_regulations)[1] <- "Locus"
finished_data_table <- merge(tiling_array_regulations, rnaseq_regulations, by="Locus")
finished_data_table <- data.frame(lapply(finished_data_table, as.character), stringsAsFactors = FALSE)



stress_row <- c(1,2, lapply(1:(ncol(finished_data_table) -2), function(n) 
{
  if(n %in% stress_pos) return("STRESS")
  return("CONTROL")
}))

finished_data_table <- rbind(as.data.frame(finished_data_table), stress_row)

write.table(finished_data_table, 
            file = "Results/RNAseq_and_tiling_preprocessed",
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

config <- read.table("Results/configuration.conf", stringsAsFactors = FALSE)
config["misclassification_cost:",] <- paste(length(unique(stress_pos)), 
                                            (length(finished_data_table) - 2 - length(unique(stress_pos))), sep = ",")
colnames(config)[1] <- "[parameters]"
write.table(config, "Results/configuration.conf", row.names = TRUE, quote = FALSE)


cat("Normalization complete\n")

