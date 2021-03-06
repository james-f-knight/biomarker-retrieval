## Take the biomarkers identified by RGIFE and select the interactions in
#  BacillusRegNet that contain them as target or source node

#  It needs the following input:
#     - The Experimentunion_model.txt file released as output of policies.py
#       contained in ~/Results/Experimentunion_model.txt
#     - The ~/Results/RegNetMatrix_integrated
#       containing the network whose nodes are all contained in the expression
#       data


# Load RegNet
RegNetMatrix_comp_comp <- read.delim("Results/RegNetMatrix_integrated", 
  stringsAsFactors=FALSE)

# Load rgife result
hand<-file("Results/polices_result/Experimentunion_model.txt")
rgife <- readLines(hand)
close(hand)

# Look for the ID in rgife in RegNetMatrix_integrated columns 1 and 8
v<-vector("numeric") #vector for storing matching indexes
for (a in 1:length(rgife)){
  for (z in 1:nrow(RegNetMatrix_comp_comp)){
    pattern1<- paste("^", RegNetMatrix_comp_comp[z,1], "$", sep = "")
    if (grepl(pattern = pattern1, rgife[a])){
      v<-append(v, z)
    }
    pattern2<- paste("^", RegNetMatrix_comp_comp[z,8], "$", sep = "")
    if (grepl(pattern = pattern2, rgife[a])){
      v<-append(v, z)
    }
  }
}


# Keep only the important rows
RegNet_biomarker <- RegNetMatrix_comp_comp[v,]

# Remove interations to empty node
RegNet_biomarker <- RegNet_biomarker[-which(RegNet_biomarker[,8] == ""),]

# Save
write.table(RegNet_biomarker, file = "Results/RGIFE_to_cytoscape", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
