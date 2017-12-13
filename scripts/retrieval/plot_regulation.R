argum <- commandArgs()

table <- read.table(paste(argum[6], "Results/RNAseq_and_tiling_preprocessed", sep= ""), sep="\t", header = TRUE, stringsAsFactors = FALSE)
locus_location <- argum[7]
loci <-  read.delim(paste(argum[6], "Results/", argum[7], sep=""), skip=1, stringsAsFactors = FALSE)[,1]

pos <- which(table[,1] %in% loci) 
name <- table[pos,2]

colour <- sapply(table[nrow(table),3:ncol(table)], function(x)
                                                   {
                                                     if(x=="CONTROL") return("black")
                                                     return("red")
                                                   })
if(!dir.exists(paste(argum[6], "Results/Plots", sep="")))
{
  dir.create(paste(argum[6], "Results/Plots", sep=""))
}
for(n in 1:length(pos))
{
  png(filename = paste(argum[6], "Results/Plots/", loci[n], sep = ""), width = 1920, height = 1080)
  plot(as.numeric(table[pos[n],3:ncol(table)]), col=colour, type = "b",
       xlab = "Sample Index", ylab = paste(name[n],"(",  loci[n], ") regulation", sep=""))
  x <- dev.off()
}

