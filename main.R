#' Optimize DADA2 filtering and trimming parameters for paired-end reads obtained from the illumina platform.
#' 
#' @param fqF A string, the file path where the forward reads .fastq is located
#'
#' @param fqR A string, the file path where the reverse reads .fastq is located
#'
#' @param inspected_eeMax A vector, the upper limits of expected errors that you want to test
#'
#' @param required_len An integer, overlap length + amplified region length
#' 
#' @param read_len An integer representing the sequencing read length
#' 
#' @return A data.frame that records the optimal filtering parameters of different trimming method, sorted according to the tradeoff of quality and read number (from high to low)
#' 
#' @example  optmize_dada2_qc("./data/sample1_R1.fastq", "./data/sample1_R2.fastq", 1:5, 253 + 20, 175)

library(ShortRead)
library(tidyverse)

# Transform quality score into sequencing error probability
phred_to_errProb <- function (x) {
  return(10^-(x/10))
}
# Transform a FastqQuality object into a matrix, and transform decode quality encoding simultaneously
rfq_to_matrix <- function (rfq) {
  mtx <- as(quality(rfq), "matrix")
  rownames(mtx) <- as.character(rfq@id)
  return(mtx)
}
# Calculate expected errors for each length of reads in a error probability matrix
ee_per_len <- function (errProb, breaks) {
  return(
    lapply(breaks, function (p) { rowSums(errProb[, 1:p]) })
  )
}

# Get post-filtering stats
get_filt_stat <- function (rfqF, rfqR, inspected_eeMax, required_len, read_len) {
  min_len <- required_len - read_len
  # FastqQuality to expected errors 
  errProbF <- phred_to_errProb(rfq_to_matrix(rfqF))
  errProbR<- phred_to_errProb(rfq_to_matrix(rfqR))
  eeF <- ee_per_len(errProbF, min_len:ncol(errProbF))
  eeR <- ee_per_len(errProbR, min_len:ncol(errProbR))
  # Create empty 
  filt <- 
    expand.grid(eeMaxF = inspected_eeMax,
                lenF = min_len:read_len,
                eeMaxR = inspected_eeMax,
                lenR = min_len:read_len)
  filt <- filt[(filt$lenF + filt$lenR) == required_len, ]
  filt$total <- nrow(errProbF)
  filt$pass <- NA
  # Remaining reads are the reads that pass the filtering in both directions
  filt$pass <- 
    unlist(lapply(1:nrow(filt),
                  function (i) {
                    return(sum((eeF[[filt$lenF[i] + 1 - min_len]] <= filt$eeMaxF[i]) * (eeR[[filt$lenR[i] + 1 - min_len]] <=  filt$eeMaxR[i]), na.rm = TRUE))
                  }
    ))
  return(filt)
}

# Add up read number retained after filtering for all samples 
bind_filt_stat <- function(fileF, fileR, inspected_eeMax, required_len, read_len) {
  rfq <- list("F" = readFastq(fileF[1]), "R" = readFastq(fileR[1]))
  filt <- get_filt_stat(rfq[["F"]], rfq[["R"]], inspected_eeMax, required_len, read_len)
  for(i in 2:length(fileF)) {
    rfq <- list("F" = readFastq(fileF[i]), "R" = readFastq(fileR[i]))
    errProb <- lapply(rfq, function (x) { phred_to_errProb(rfq_to_matrix(x)) })
    tmp <- get_filt_stat(rfq[["F"]], rfq[["R"]], inspected_eeMax, required_len, read_len)
    filt[, c("pass", "total")] <- filt[, c("pass", "total")] + tmp[, c("pass", "total")]
    gc()
  }
  filt$pass_ratio <- filt$pass/filt$total
  # Ranking index is used to assess the quality loss of the forward and reverse reads
  filt$ranking_index <- (filt$eeMaxF^2 + filt$eeMaxR^2)^0.5
  return(filt)
}

# Select the best filtering parameters for each length cutting position
optim_len_cutoff <- function (filt) {
  # Exclude parameters that lose more reads under the same filtering criteria
  optim_len <- by(filt, filt[, c("eeMaxR", "eeMaxF")], function(x) rownames(x)[which.max(x$pass)])
  filt <- filt[unlist(optim_len), ]
  # Exclude parameters that lose more reads under the looser filtering criteria
  params <- filt[order(filt$ranking_index, -filt$pass_ratio), ]
  for (i in 2:nrow(params)) {
    if (params[i, "pass_ratio"] <= max(params[1:i-1, "pass_ratio"])) {
      params[i, "ranking_index"] <- NA
    }
  }
  params <- params[!is.na(params$ranking_index), ]
  return(params)
}


# Extract optimal length cutoff and parameter params
sort_params <- function(params) {
  fo_diff <- sapply(1:nrow(params),
                    function(pos){
                      left <- (params[pos, "pass_ratio"] - params[1, "pass_ratio"]) / (params[pos, "ranking_index"])
                      right <- (params[nrow(params), "pass_ratio"] - params[pos, "pass_ratio"]) / (params[nrow(params) , "ranking_index"]- params[pos, "ranking_index"])
                      return(left - right)
                    }
  )
  params$fo_diff <- fo_diff
  params <- params[order(-params$fo_diff), ]
  return(params)
}

optimize_dada2_qc <- function (fqF,fqR,inspected_eeMax, required_len, read_len) {
  filt <- bind_filt_stat(fqF, fqR, inspected_eeMax, required_len, read_len)
  params <- optim_len_cutoff(filt)
  optim_params <- sort_params(params)
  return(optim_params)
}