# From tab separated table to .arff
Exp_table <- read.delim(
  "Results/RNAseq_and_tiling_preprocessed",
  stringsAsFactors = FALSE)

# Create the output file
filehand <-file("Results/data.arff")



# Header section
header <- vector ("character", length = nrow(Exp_table))
for (a in 1:nrow(Exp_table)){
  if (a < nrow(Exp_table)){
    header[a]<-paste("@attribute", Exp_table[a,1], "numeric", sep = " ")
  }
  if (a == nrow(Exp_table)){
    header[a]<-paste("@attribute", Exp_table[a,1], "{CONTROL, STRESS}", 
                     sep = " ")
  }
}


# Data section
data <- vector ("character", length = ncol(Exp_table))
for (a in 3:ncol(Exp_table)){
  data[a]<-paste(Exp_table[,a], collapse = ",")
  data[a]<-c(data[a], "\r")
}


# Assembly
writeLines(c(header, "@data", data), filehand)
close(filehand)
print("File converted to .arff")
