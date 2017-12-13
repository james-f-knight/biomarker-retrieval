find_positions <- function(data, stress_pos)
{
  below_min <- NULL
  above_min <- NULL
  mins <- NULL
  above_min_bits <- NULL
  for(n in 1:nrow(data))
  {
    n_order <- order(-data[n,], 1:length(data[n,]) %in% stress_pos)

    min_stress_sample <- intersect(n_order, stress_pos)[length(stress_pos)]
    min_stress_rank <- which(n_order == min_stress_sample)
    
    below_min[[length(below_min) + 1]] <- n_order[(min_stress_rank + 1):length(n_order)]
    
    above <- setdiff(n_order[1:min_stress_rank], stress_pos)
    above_min[[length(above_min) + 1]] <- above 
    mins <- c(mins, data[n, min_stress_sample])
    
    vals <- NULL
    for(r in seq(from = 31, by = 31, length.out = ceiling(ncol(data)/31) ))
    {
      vals <- c(vals, strtoi(paste(ifelse(r:(r - 30) %in% above, "1", "0"), collapse = ""), base = 2))
    }
    above_min_bits <- rbind(above_min_bits, vals)
    row.names(above_min_bits)[nrow(above_min_bits)] <- loci[n] 
    if(n %% 100 == 0) {
     # cat(n)
     # cat("\n")
      cat(".")
    }
  }
  return(list(mins, above_min_bits, above_min, below_min)) 
}

n_optimal_thresholds <- function(ps, data, stress_pos, preprocess, draw)
{
  above <- preprocess[[3]][ps]
  best_order <- order(sapply(above, length))
  
  above <- above[best_order]
  mins <- preprocess[[1]][ps][best_order]
  below <- preprocess[[4]][ps][best_order] 
  
  
  res <- optimal_thresholds(ps[best_order], above, below, mins, data, stress_pos)
  if(draw == TRUE) n_draw_graphs(ps[best_order], mins, res, data, stress_pos)
  rev_order <- order(best_order)
  return(c((mins - res)[rev_order], ((mins + res)/2)[rev_order] ))
}


optimal_thresholds <- function(ps, above, below, mins, data, stress_pos)
{
  best_min <- 0
  ls_min = 0
  this_below <- below[[1]]
  this_above <- above[[1]]
  
  if(length(ps) == 1) 
  {
    if(length(this_below) == 0)
    {
      return(-100)
    }
    return(data[ps[1],this_below[1]])
  }
  
  misclass_below <- lapply(below[-1], intersect, y=this_above)
  
  for (n in 0:(length(this_below)))
  {
    this_above <- c(this_above, this_below[n])
    l1 <- data[ps[1],this_below[n+1]]
    if(n == length(this_below)) 
    {
      l1 <- -100  
    }  
    
    if(n > 0) 
    {
      is_above <- sapply(above[-1], is.element, el=this_below[n])
      
      if(all(is_above)) 
        {
          return(best_ls)
        }
      below_features <- which(!is_above)
      misclass_below[below_features] <- lapply(below[-1][below_features], intersect, y=this_above)
    }
    
    ls <- optimal_thresholds(ps[-1], above[-1], misclass_below, mins[-1], data, stress_pos)
    
    test_margins <- mins - c(l1, ls)
    margin_order <- order(test_margins)
    
    test_min <- test_margins[margin_order]
    test_position <- min(which((test_min != best_min)))
    if(!is.infinite(test_position))
    {
      if(test_min[test_position] > best_min[test_position])
      {
        best_ls <- c(l1, ls)
        best_min <- test_min 
      }
    }
    if(margin_order[length(margin_order)] == 1 & ls_min > min(ls))
    {
      return(best_ls)
    }
    ls_min <- min(ls)
    
  }
  return(best_ls)
}

and_classifier_train <- function() {} ### REMOVE SPECIFICS AND PACKAGE ###


n_draw_graphs <- function(ps, mins, ls, data, stress_pos)
{
  par(mfrow = c(1,length(ps)))
  for(n in 1:length(ps))
  {
    col=unlist(sapply(1:ncol(data), n_get_colour, ps=ps, n=n, ls=ls, data=data, stress_pos=stress_pos))
    plot(data[ps[n],], type="p", ylab = loci[ps[n]], 
         col=col)
    lines(c(-10,100), c(ls[n],ls[n]), col="blue")
    lines(c(-10,100), c(mins[n],mins[n]), col="blue")
  }
}

n_get_colour <- function(ps, n, point, ls, data, stress_pos) 
{
  ##cat(ls)
  #cat(" ")
  if(point %in% stress_pos) 
  {
    if(data[ps[n],point] > ls[n])
    {
      return("green") 
    } 
    return("cyan")
  }
  if(data[ps[n],point] <= ls[n])
  {
    if(all(data[ps[-n],point] < ls[-n])) 
    {
      return("black")
    }
    return("orange")
  }
  return("red")
}

test <- function(ps, data, stress_pos, result)
{
  
  boundaries <- result[(length(ps) + 1):length(result)]
  
  if(length(boundaries) == 1)
  {
    true_pos <- sum(data[ps,stress_pos] > boundaries)
    false_pos <- sum(data[ps,-stress_pos] > boundaries)
  }
  else
  {
    if(length(stress_pos) > 0)
    {
      true_pos <- sum(apply(data.frame(data[ps,stress_pos]), 2, function(xs) {all(xs > boundaries)}))
      false_pos <- sum(apply(data.frame(data[ps,-stress_pos]), 2, function(xs) {all(xs > boundaries)}))
    }
    else
    {
      true_pos <- 0
      false_pos <- sum(apply(data[ps,], 2, function(xs) {all(xs > boundaries)}))
    }
  }
  
  
  return(c(true_pos, length(stress_pos) - true_pos, false_pos, ncol(data) - length(stress_pos) - false_pos))
}


library(genefilter)
library(xlsx)
table <- read.delim("Results/RNAseq_and_tiling_preprocessed", stringsAsFactors = FALSE)
data <- table[,3:ncol(table)]
stress_pos <- which(data[nrow(data),] =="STRESS")
data <- sapply(data[1:(nrow(data) - 1),], as.numeric)
loci <- table[,1]
genes <- table[,2]

factor <- factor(sapply(1:ncol(data), function(n) 
{
  if(n %in% stress_pos) return(1)
  return(0)
}))

#################### HIGHEST MARGIN BOUNDARIES ########################

cat("\n")
cat("Preprocessing full data set...")
cat("\n")
#check_pos <- which(apply(data[,stress_pos], 1, min) > 0)
#check_pos_down <- which(apply(-data[,stress_pos], 1, min) > 0)
tests <- rowttests(data[,stress_pos]) 
check_pos <- which(tests$dm > 0 & tests$p.value < 0.05)
check_pos_down <- which(tests$dm < 0 & tests$p.value < 0.05)

biomarkers <- c(check_pos, check_pos_down)
biomarker_loci <- loci[biomarkers]
data[check_pos_down,] <- -data[check_pos_down,]
pairs <- combn(biomarkers, 2)


preprocess <- find_positions(data, stress_pos)
above_min_bits <- preprocess[[2]]

groups <- sample(c(1:(ncol(data) %% 10), rep(1:10, floor(ncol(data) / 10))), ncol(data))
preprocess_cv <- NULL

cat("\n")
cat("Preprocessing cross-validation datasets...")
cat("\n")
for(n in 1:10)
{
  cat(n)
  preprocess_cv[[length(preprocess_cv) + 1]] <- find_positions(data[,which(groups != n)], 
                                                               which(factor[which(groups != n)] == 1))
  cat("\n")
}


{
  cat("\n")
  cat("Computing individual biomarkers...")
  cat("\n")
  separator_individual <- NULL
  for (n in 1:length(biomarkers)) 
  {
    p1 <- biomarkers[n]
    res <- above_min_bits[p1,]
    if(all(res == 0))
    {
      separator_individual <- c(separator_individual, p1)
    }
  }
  cat("Individual biomarkers found: ")
  cat(length(separator_individual))
  cat("\n\n")
  t <- Sys.time()
  cat("Computing biomarker pairs...")
  cat("\n")
  count <- 0
  separator_pairs <- NULL
  for (n in 1:ncol(pairs)) 
  {
    p1 <- pairs[1,n]
    p2 <- pairs[2,n]
    res <- bitwAnd(above_min_bits[p1,], above_min_bits[p2,])
    if(all(res == 0))
    {
      count <- count + 1
      separator_pairs <- cbind(separator_pairs, c(pairs[1,n], pairs[2,n]))
    }
    if(n %% 10000 == 0) cat(".")
  }
  cat("\n")
  cat(paste("Biomarker pairs found:", ncol(separator_pairs), "\nTime taken:", Sys.time() - t, "\n", sep=" "))
}

if(length(separator_individual) > 0)
{
  cat("\n")
  cat("Evaluating individual biomarkers...")
  margins <- NULL
  conf <- matrix(0, length(separator_individual), 4)
  
  for(n in 1:10)
  {
    cat("\n")
    cat("Fold ")
    cat(n)
    data_train <- data[,which(groups != n)]
    stress_pos_train <- which(factor[which(groups != n)] == 1)
    data_test <- data[,which(groups == n)]
    stress_pos_test <- which(factor[which(groups == n)] == 1)
    for(m in 1:length(separator_individual))
    {
      cat(".")
      result <- n_optimal_thresholds(separator_individual[m], data_train, stress_pos_train, preprocess_cv[[n]], 0)
      conf[m,] <- conf[m,] + test(separator_individual[m], data_test, stress_pos_test, result)
    }
  }
  
  cat("\n")
  cat("Computing individual classifiers...")
  cat("\n")
  
  for(p in separator_individual)
  {
    cat(".")
    margins <- rbind(margins, n_optimal_thresholds(p, data, stress_pos, preprocess, 0))
  }
  
  cat("\n\n")
  
  individual_margins <- cbind(separator_individual, margins, conf, sqrt(conf[,1]/(conf[,1] + conf[,2]) * 
                                                                          conf[,4]/(conf[,3] + conf[,4])))
  individual_margins <- individual_margins[order(-individual_margins[,8],-individual_margins[,2]),]
  
  colnames(individual_margins) <- c("Gene", "Margin Size", "Threshold", "TP", "FP", "FN", "TN", "gmean accuracy")
  
  png(filename = "best_individual_plot", width = 1920, height = 1080)
  n_optimal_thresholds(individual_margins[1,1], data, stress_pos, preprocess, 1)
  x <- dev.off()
  
  individual_margins[,1] <- genes[individual_margins[,1]]
  write.csv(individual_margins, file = "individual_biomarkers", row.names = FALSE)
  write.xlsx(write_pm, file="biomarker_pairs.xls", sheetName = "sheet1" , row.names = FALSE)
}




########## PAIRS ############

cat("\n")
cat("Evaluating biomarker pairs...")

margins <- NULL
conf <- matrix(0, ncol(separator_pairs), 4)

for(n in 1:10)
{
  t <- Sys.time()
  cat("\n")
  cat("Fold ")
  cat(n)
  data_train <- data[,which(groups != n)]
  stress_pos_train <- which(factor[which(groups != n)] == 1)
  data_test <- data[,which(groups == n)]
  stress_pos_test <- which(factor[which(groups == n)] == 1)
  count <- 0
  for(m in 1:ncol(separator_pairs))
  {
    result <- n_optimal_thresholds(separator_pairs[,m], data_train, stress_pos_train, preprocess_cv[[n]], 0)
    conf[m,] <- conf[m,] + test(separator_pairs[,m], data_test, stress_pos_test, result)
    if(count %% 1000 == 0) cat(".") 
    count <- count + 1
  }
  cat("\n")
  cat("Time: ")
  cat(Sys.time() - t)
}

cat("\n")
cat("Computing pair classifiers...")
cat("\n")
t <- Sys.time()
count <- 0
for(n in 1:ncol(separator_pairs))
{
  margins <- rbind(margins, n_optimal_thresholds(separator_pairs[,n], data, stress_pos, preprocess, 0))
  if(count %% 1000 == 0) cat(".") 
  count <- count + 1
}
cat("\n")
cat("Time: ")
cat(Sys.time() - t)


pair_margins <- cbind(t(separator_pairs), margins, conf, sqrt(conf[,1]/(conf[,1] + conf[,2]) * 
                                                              conf[,4]/(conf[,3] + conf[,4])))
pair_margins <- pair_margins[order(-pair_margins[,ncol(pair_margins)],-apply(pair_margins[,3:4], 1, min), -apply(pair_margins[,3:4], 1, max)),]

write_pm <- pair_margins

colnames(write_pm) <- c("Gene 1", "Gene 2", "Margin 1", "Margin 2",
                            "Threshold 1", "Threshold 2", "TP", "FP", "FN", "TN", "gmean accuracy")

write_pm[,1] <- genes[pair_margins[,1]]
write_pm[,2] <- genes[pair_margins[,2]]
write.csv(write_pm, file="biomarker_pairs", row.names = FALSE)
write.xlsx(write_pm, file="biomarker_pairs.xls", sheetName = "sheet1" , row.names = FALSE)

png(filename = "best_pair_plot", width = 1920, height = 1080)
n_optimal_thresholds(pair_margins[1,1:2], data, stress_pos, preprocess, 1)
x <- dev.off()


cat("\n")
cat("Finding extensions which improve margin size...")
cat("\n")

previous_margins <- pair_margins[1:100,]
new_margins <- NULL
for(n in 1:nrow(previous_margins))
{
  prev <- previous_margins[n,]
  for(p in biomarkers)
  {
    result <- n_optimal_thresholds(c(prev[1:2], p), data, stress_pos, preprocess, 0)
    if(result[3] < 90)
    {
      new_margins <- rbind(new_margins, c(prev[1:2], p, result))
    }
  }
  cat(".")
}

cat("\n")
cat("Evaluating extensions...")
cat("\n")

conf <- matrix(0, nrow(new_margins), 4)
for(n in 1:10)
{
  t <- Sys.time()
  cat("\n")
  cat("Fold ")
  cat(n)
  data_train <- data[,which(groups != n)]
  stress_pos_train <- which(factor[which(groups != n)] == 1)
  data_test <- data[,which(groups == n)]
  stress_pos_test <- which(factor[which(groups == n)] == 1)
  count <- 0
  for(m in 1:nrow(new_margins))
  {
    result <- n_optimal_thresholds(new_margins[m,1:3], data_train, stress_pos_train, preprocess_cv[[n]], 0)
    conf[m,] <- conf[m,] + test(new_margins[m,1:3], data_test, stress_pos_test, result)
    if(count %% 500 == 0) cat(".")
    count <- count + 1
  }
  cat("\n")
  cat("Time: ")
  cat(Sys.time() - t)
}

new_results <- cbind(new_margins, conf, sqrt(conf[,1]/(conf[,1] + conf[,2]) * 
                                             conf[,4]/(conf[,3] + conf[,4])))

sort_margins <- t(apply(new_results[,4:6], 1, sort)) 

new_results <- new_results[order(-new_results[,ncol(new_results)], -sort_margins[,1], -sort_margins[,2], -sort_margins[,3]),]


write_nm <- new_results

colnames(write_nm) <- c("Gene 1", "Gene 2", "Gene 3", "Margin 1", "Margin 2", "Margin 3",
                        "Threshold 1", "Threshold 2", "Threshold 3", "TP", "FP", "FN", "TN", "gmean accuracy")

write_nm[,1] <- genes[new_results[,1]]
write_nm[,2] <- genes[new_results[,2]]
write_nm[,3] <- genes[new_results[,3]]

write.csv(write_nm, file="biomarker_triples", row.names = FALSE)
write.xlsx(write_nm, file="biomarker_triples.xls", sheetName = "sheet1" , row.names = FALSE)

png(filename = "best_triple_plot", width = 1920, height = 1080)
  n_optimal_thresholds(new_results[1,1:3], data, stress_pos, preprocess, 1)
x <- dev.off()
