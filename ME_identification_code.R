#to get IAP methylation data from bedgraph files:
#bedtools map -a iap_merge_d150.sorted.bed -b "bedgrph file WGBS or WGoxBS" -c 4 -o collapse

meth_data <-c(list.files(pattern="*_modified.csv"))

cpg_16start_average <- function(data_frame) {
  cpg_start <- data.frame(matrix(nrow=nrow(data_frame),ncol=1))
  for (x in 1:16) {
    cpg_start[,x] <- data_frame[,x+2]
  }
  average_vector <- c()
  for (x in 1:nrow(data_frame)) {
    a <- c(as.numeric(cpg_start[x,]))
    average_vector <- c(average_vector,(mean(a)))
  }
  return(average_vector)
}

cpg_16end_average <- function(data_frame) {
  average_vector <- c()
  for (x in 1:nrow(data_frame)) {
    a <- c(data_frame[x,])
    a <- a[3:ncol(data_frame)]
    k <- as.numeric(a[!is.na(a)])
    n <- as.numeric(length(k))
    if (n>15) {
      l <- n-15
      kk <- k[l:n]
      average_vector <- c(average_vector,(mean(kk)))} else
      {average_vector <- c(average_vector, NA)}
  }
  return(average_vector)
}

second_max <- function(data_frame) {
  second_max_vector <- c()
  for (x in 1:nrow(data_frame)) {
    a <- c(as.numeric(data_frame[x,c(1:16)]))
    a <- sort(a, decreasing = TRUE)
    second_max_vector <- c(second_max_vector,a[2])
  }
  return(second_max_vector)
}

second_min <- function(data_frame) {
  second_min_vector <- c()
  for (x in 1:nrow(data_frame)) {
    a <- c(as.numeric(data_frame[x,c(1:16)]))
    a <- sort(a)
    second_min_vector <- c(second_min_vector,a[2])
  }
  return(second_min_vector)  
}


bl6_start_all <- data.frame(row.names = (read.table(meth_data[1], sep="\t", fill= TRUE, stringsAsFactors = FALSE))[,1])
for (x in 1:length(meth_data)) {
  bl6_start_all[,x] <- cpg_16start_average(read.table(meth_data[x], sep="\t", fill= TRUE, stringsAsFactors = FALSE))
}
bl6_start_all[,17] <- (read.table(meth_data[1], sep="\t", fill= TRUE, stringsAsFactors = FALSE))[,2]
colnames(bl6_start_all) <- c(meth_data,"strand")


bl6_end_all <- data.frame(row.names = (read.table(meth_data[1], sep="\t", fill= TRUE, stringsAsFactors = FALSE))[,1])
for (x in 1:length(meth_data)) {
  bl6_end_all[,x] <- cpg_16end_average(read.table(meth_data[x], sep="\t", fill= TRUE, stringsAsFactors = FALSE))
}
bl6_end_all[,17] <- (read.table(meth_data[1], sep="\t", fill= TRUE, stringsAsFactors = FALSE))[,2]
colnames(bl6_end_all) <- c(meth_data,"strand")

#remove rows that have NA - bad sequence
bl6_start_cleaned <- bl6_start_all[complete.cases(bl6_start_all),]
bl6_end_cleaned <- bl6_end_all[complete.cases(bl6_end_all),]

#computational variation range

bl6_start_sec_max <- second_max(bl6_start_cleaned)
bl6_start_sec_min <- second_min(bl6_start_cleaned)
bl6_end_sec_max <- second_max(bl6_end_cleaned)
bl6_end_sec_min <- second_min(bl6_end_cleaned)

delta_start <- c()
for (x in 1:length(bl6_start_sec_max)) {
  a <- bl6_start_sec_max[x] - bl6_start_sec_min[x]
  delta_start <- c(delta_start,a)
}

delta_end <- c()
for (x in 1:length(bl6_end_sec_max)) {
  a <- bl6_end_sec_max[x] - bl6_end_sec_min[x]
  delta_end <- c(delta_end,a)
}

results <- data.frame(row.names = row.names(bl6_end_cleaned))
results[,1] <- bl6_end_cleaned[,17]
results[,2] <- delta_start
results[,3] <- delta_end
colnames(results) <- c("strand","delta_5end","delta_3end")

antisense <- subset(results,results$strand == "-")
sense <- subset(results,results$strand == "+")

final_results <- data.frame(row.names = c(row.names(antisense),row.names(sense)))
final_results <- cbind(final_results,c(antisense$strand,sense$strand),c(antisense$delta_3end,sense$delta_5end),c(antisense$delta_5end,sense$delta_3end))
colnames(final_results) <- c("strand","delta_5LTR","delta_3LTR")

ME_candidates <- final_results[(final_results$delta_5LTR>25)|(final_results$delta_3LTR>25),]

# cell type specific and sex specific DMRs need to be removed from the list of ME candidates

list_ME_LTR5<- row.names(ME_candidates[ME_candidates$delta_5LTR>25,])
list_ME_LTR3<- row.names(ME_candidates[ME_candidates$delta_3LTR>25,])
LTR5 <- rbind(bl6_start_cleaned[bl6_start_cleaned$strand=="+",],bl6_end_cleaned[bl6_end_cleaned$strand=="-",])
LTR3 <- rbind(bl6_start_cleaned[bl6_start_cleaned$strand=="-",],bl6_end_cleaned[bl6_end_cleaned$strand=="+",])
ME_LTR5 <- LTR5[list_ME_LTR5,]
ME_LTR3 <- LTR3[list_ME_LTR3,]

#DMR tests- should identify majority of B/T cell DMRs
#whether all maximum or all minimum values coming from one cell type
#sex specific DMRs can be tested the same way

meth_data_Bcell <- c(list.files(pattern="*_b_[1,2]_d150_modified.csv"))
meth_data_Tcell <- c(list.files(pattern="*_t_[1,2]_d150_modified.csv"))

bl6_start_cleaned <- bl6_start_all[complete.cases(bl6_start_all),]
bl6_end_cleaned <- bl6_end_all[complete.cases(bl6_end_all),]

result_of_DMR_test_max <- c()
for (x in 1:nrow(bl6_start_cleaned)) {
  a <- c(as.numeric(bl6_start_cleaned[x,c(1:16)]))
  a <- sort(a, decreasing = TRUE)
  test_for_DMR <- c()
  for (k in 1:8) {
    if (a[k] %in% bl6_start_cleaned[x, meth_data_Bcell]) {
      test_for_DMR <- c(test_for_DMR, TRUE)}
    else {
      test_for_DMR <- c(test_for_DMR, FALSE)
    }
  }
  FT <- sum(test_for_DMR == FALSE)
  TF <- sum(test_for_DMR == TRUE)
  if ((FT > 0) & (TF > 0)) {
    result_of_DMR_test_max <- c(result_of_DMR_test_max, "not_DMR")
  } else {
    result_of_DMR_test_max <- c(result_of_DMR_test_max, "DMR")
  }
}

bl6_start_cleaned[,18]<- result_of_DMR_test_max


result_of_DMR_test_min <- c()
for (x in 1:nrow(bl6_start_cleaned)) {
  a <- c(as.numeric(bl6_start_cleaned[x,c(1:16)]))
  a <- sort(a)
  test_for_DMR <- c()
  for (k in 1:8) {
    if (a[k] %in% bl6_start_cleaned[x, meth_data_Bcell]) {
      test_for_DMR <- c(test_for_DMR, TRUE)}
    else {
      test_for_DMR <- c(test_for_DMR, FALSE)
    }
  }
  FT <- sum(test_for_DMR == FALSE)
  TF <- sum(test_for_DMR == TRUE)
  if ((FT > 1) & (TF > 1)) {
    result_of_DMR_test_min <- c(result_of_DMR_test_min, "not_DMR")
  } else {
    result_of_DMR_test_min <- c(result_of_DMR_test_min, "DMR")
  }
}

bl6_start_cleaned[,19]<- result_of_DMR_test_min


result_of_DMR_test_max_end <- c()
for (x in 1:nrow(bl6_end_cleaned)) {
  a <- c(as.numeric(bl6_end_cleaned[x,c(1:16)]))
  a <- sort(a, decreasing = TRUE)
  test_for_DMR <- c()
  for (k in 1:8) {
    if (a[k] %in% bl6_end_cleaned[x, meth_data_Bcell]) {
      test_for_DMR <- c(test_for_DMR, TRUE)}
    else {
      test_for_DMR <- c(test_for_DMR, FALSE)
    }
  }
  FT <- sum(test_for_DMR == FALSE)
  TF <- sum(test_for_DMR == TRUE)
  if ((FT > 0) & (TF > 0)) {
    result_of_DMR_test_max_end <- c(result_of_DMR_test_max_end, "not_DMR")
  } else {
    result_of_DMR_test_max_end <- c(result_of_DMR_test_max_end, "DMR")
  }
}

bl6_end_cleaned[,18]<- result_of_DMR_test_max_end


result_of_DMR_test_min_end <- c()
for (x in 1:nrow(bl6_end_cleaned)) {
  a <- c(as.numeric(bl6_end_cleaned[x,c(1:16)]))
  a <- sort(a)
  test_for_DMR <- c()
  for (k in 1:8) {
    if (a[k] %in% bl6_end_cleaned[x, meth_data_Bcell]) {
      test_for_DMR <- c(test_for_DMR, TRUE)}
    else {
      test_for_DMR <- c(test_for_DMR, FALSE)
    }
  }
  FT <- sum(test_for_DMR == FALSE)
  TF <- sum(test_for_DMR == TRUE)
  if ((FT > 0) & (TF > 0)) {
    result_of_DMR_test_min_end <- c(result_of_DMR_test_min_end, "not_DMR")
  } else {
    result_of_DMR_test_min_end <- c(result_of_DMR_test_min_end, "DMR")
  }
}

bl6_end_cleaned[,19]<- result_of_DMR_test_min_end


#f and t tests

meth_data_Bcell <- c(list.files(pattern="*_b_[1,2]_d150_modified.csv"))
meth_data_Tcell <- c(list.files(pattern="*_t_[1,2]_d150_modified.csv"))

f_test_results <- c()
for (x in 1:nrow(bl6_start_cleaned)) {
  b_cell <- as.numeric(c(bl6_start_cleaned[x, meth_data_Bcell]))
  t_cell <- as.numeric(c(bl6_start_cleaned[x, meth_data_Tcell]))
  zz <- var.test(b_cell, t_cell)
  f_test_results <- c(f_test_results, zz$p.value)
}

t_test_result_bvst <- c()
for (x in 1:nrow(bl6_start_cleaned)) {
  b_cell <- as.numeric(c(bl6_start_cleaned[x, meth_data_Bcell]))
  t_cell <- as.numeric(c(bl6_start_cleaned[x, meth_data_Tcell]))
  zz <- t.test(b_cell, t_cell)
  t_test_result_bvst <- c(t_test_result_bvst, zz$p.value)
}

bl6_start_cleaned[,20] <- f_test_results
bl6_start_cleaned[,21] <- t_test_result_bvst

f_test_results <- c()
for (x in 1:nrow(bl6_end_cleaned)) {
  b_cell <- as.numeric(c(bl6_end_cleaned[x, meth_data_Bcell]))
  t_cell <- as.numeric(c(bl6_end_cleaned[x, meth_data_Tcell]))
  zz <- var.test(b_cell, t_cell)
  f_test_results <- c(f_test_results, zz$p.value)
}

t_test_result_bvst <- c()
for (x in 1:nrow(bl6_end_cleaned)) {
  b_cell <- as.numeric(c(bl6_end_cleaned[x, meth_data_Bcell]))
  t_cell <- as.numeric(c(bl6_end_cleaned[x, meth_data_Tcell]))
  zz <- t.test(b_cell, t_cell)
  t_test_result_bvst <- c(t_test_result_bvst, zz$p.value)
}

bl6_end_cleaned[,20] <- f_test_results
bl6_end_cleaned[,21] <- t_test_result_bvst
