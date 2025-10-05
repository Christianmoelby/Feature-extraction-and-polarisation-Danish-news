##### Setup
rm(list = ls())
#install.packages("rtools")
library(lubridate)
library(data.table)
library(dplyr)
library(ggplot2)
setwd("WD")
input_matrix_csv <- "DATA.csv"
input_matrix <- read.csv(input_matrix_csv, sep = ";", encoding = "UTF-8")

input_matrix <- as.data.frame(lapply(input_matrix, function(x){
  if (is.character(x)) tolower(x) else x
}))


##### Definer parametre: 
set.seed(1234)                        # for reproducibility
K <- 100                              # number of subsamples
subsample_size_denominator <- 10      # each subsample: approx. one-tenth of full sample
CI_top_limit <- 90                    # top for confidence interval (10%)
CI_bottom_limit <- 11                 # bottom for confidence interval (10%)
min_token_in_j <- 50                   # minimum antal gange token skal optræde i løbet af perioden
min_token_weeks <-10                  # minimum antal forskellige uger et token skal optræde i løbet af et år
maks_pct_out <- 0.5                     # pct. af de mest hypige tokens der filtreres ud
IW <- 3                               # inter week minimum antal gange token optr?der 


##### Udtræk datoer og find uger: 
dr_rows <- input_matrix[grep("^dr", input_matrix[, 1], ignore.case = TRUE), ]
tv2_rows <- input_matrix[grep("^tv2", input_matrix[, 1], ignore.case = TRUE), ]
input_matrix <- rbind(dr_rows, tv2_rows)

extract_date <- function(id) {
  parts <- unlist(strsplit(id, "_")) 
  if (length(parts) >= 2) {
    return(as.Date(parts[2], format = "%Y%m%d")) 
  } else {
    return(NA)  
  }
}


input_matrix$Date <- as.Date(sapply(input_matrix[, 1], extract_date))
input_matrix$Week <- isoweek(input_matrix$Date)
input_matrix$Year <- isoyear(input_matrix$Date)
input_matrix$Week <- ifelse(format(input_matrix$Date, "%m") == "12" & input_matrix$Week == 1, 52, input_matrix$Week)
input_matrix$Week <- ifelse(format(input_matrix$Date, "%m") == "1" & input_matrix$Week == 52, 1, input_matrix$Week)
input_matrix$Week <- ifelse(format(input_matrix$Date, "%m") == "1" & input_matrix$Week == 53, 1, input_matrix$Week)
input_matrix$Week <- ceiling(input_matrix$Week / 2)
input_matrix$Week[input_matrix$Week > 26] <- 26
weeks <- unique(input_matrix$Week)
years <- unique(input_matrix$Year)
input_matrix$BiWeek <- paste0(input_matrix$Year, "-", sprintf("%02d", input_matrix$Week))
BiWeeks <- unique(input_matrix$BiWeek)
print(BiWeeks)


##### start token filtering 
input_table <- as.data.table(input_matrix)
token_columns <- setdiff(colnames(input_table), c("ID", "Theme", "Date", "Week", "Year", "BiWeek"))
if(!all(token_columns %in% colnames(input_table))) { 
  stop("Some of the token columns are missing in the dataset!") 
}
input_table_long <- melt(input_table, id.vars = c("ID", "BiWeek"), 
                         measure.vars = token_columns, 
                         variable.name = "token_col", value.name = "token")
input_table_long <- input_table_long[!is.na(token) & token != "" & token != "NA"]
token_counts <- input_table_long[, .(count = .N, unique_biweeks = uniqueN(BiWeek)), by = token]
filtered_token_counts <- token_counts[count >= min_token_in_j & unique_biweeks >= min_token_weeks]
filtered_total_tokens <- nrow(filtered_token_counts)
top_pct_threshold <- ceiling(filtered_total_tokens * (maks_pct_out / 100))
top_mentioned_tokens <- filtered_token_counts[order(-count)][1:top_pct_threshold, token]
j_filtered_timeline <- filtered_token_counts$token
j_filtered_timeline <- setdiff(j_filtered_timeline, top_mentioned_tokens)



##### Define LO_results og begr?ns horisont
#input_matrix <- input_matrix[input_matrix$Year %in% c(2019,2020,2021,2022,2023,2024), ]                   # T?ND / SLUK for forskellige ?r
#input_matrix <- input_matrix[input_matrix$BiWeek %in% c("2019-07", "2021-07", "2021-21", "2021-23"), ]  # T?ND / SLUK for biweeks ?r
BiWeeks <- unique(input_matrix$BiWeek)
MLE_results <- matrix(NA, nrow = length(BiWeeks), ncol = 7)
colnames(MLE_results) <- c("BiWeek", "X", "CI_lower", "CI_upper", "Xr", "CIr_lower", "CIr_upper")
rownames(MLE_results) <- BiWeeks 



##### Hjælpefunktioner
dot_product <- function(q, p) {
  sum(q * p) 
}
compute_X_matrix <- function(q_matrix, p_matrix) {
  DR_names <- grep("^dr_", colnames(q_matrix), value = TRUE)
  TV2_names <- grep("^tv2_", colnames(q_matrix), value = TRUE)
  
  sum_DR <- sum(vapply(DR_names, function(name) {
    q <- as.numeric(q_matrix[, name])  
    p <- as.numeric(p_matrix[, name])  
    dot_product(q, p)
  }, numeric(1)))
  
  sum_TV2 <- sum(vapply(TV2_names, function(name) {
    q <- as.numeric(q_matrix[, name])  
    p <- as.numeric(p_matrix[, name])  
    dot_product(q, 1 - p)
  }, numeric(1)))
  
  length_DR <- max(1, length(DR_names))
  length_TV2 <- max(1, length(TV2_names))
  
  X <- 0.5 * ((1 / length_DR) * sum_DR + (1 / length_TV2) * sum_TV2)
  return(X)
}

compute_pi <- function(articles_subset, j_filtered) {
  num_tokens <- length(j_filtered)
  
  # Aggregate DR and TV2 tokens
  DR_articles_sub <- articles_subset[grepl("^DR_", names(articles_subset), ignore.case = TRUE)]
  TV2_articles_sub <- articles_subset[grepl("^TV2_", names(articles_subset), ignore.case = TRUE)]
  
  total_DR_tokens <- sum(lengths(DR_articles_sub))
  total_TV2_tokens <- sum(lengths(TV2_articles_sub))
  
  tab_DR <- if (total_DR_tokens > 0) table(unlist(DR_articles_sub)) else NULL
  tab_TV2 <- if (total_TV2_tokens > 0) table(unlist(TV2_articles_sub)) else NULL
  
  counts_DR <- if (!is.null(tab_DR)) as.numeric(tab_DR[j_filtered]) else rep(0, num_tokens)
  counts_TV2 <- if (!is.null(tab_TV2)) as.numeric(tab_TV2[j_filtered]) else rep(0, num_tokens)
  
  counts_DR[is.na(counts_DR)] <- 0
  counts_TV2[is.na(counts_TV2)] <- 0
  
  q_vector_DR <- if (total_DR_tokens > 0) counts_DR / total_DR_tokens else rep(0, num_tokens)
  q_vector_TV2 <- if (total_TV2_tokens > 0) counts_TV2 / total_TV2_tokens else rep(0, num_tokens)
  
  # Create q_matrix_sub with two columns
  q_matrix_sub <- cbind(q_vector_DR, q_vector_TV2)
  colnames(q_matrix_sub) <- c("dr_group", "tv2_group")
  rownames(q_matrix_sub) <- j_filtered
  
  # Create matching p_matrix_sub (shared across groups)
  frac_DR <- if (total_DR_tokens > 0) counts_DR / total_DR_tokens else rep(0, num_tokens)
  frac_TV2 <- if (total_TV2_tokens > 0) counts_TV2 / total_TV2_tokens else rep(0, num_tokens)
  denom <- frac_DR + frac_TV2
  p_vals <- ifelse(denom > 0, frac_DR / denom, 0)
  
  p_matrix_sub <- matrix(p_vals, nrow = num_tokens, ncol = 2)
  rownames(p_matrix_sub) <- j_filtered
  colnames(p_matrix_sub) <- c("dr_group", "tv2_group")
  
  # Compute π
  pi_est <- compute_X_matrix(q_matrix_sub, p_matrix_sub)
  return(pi_est)
}


##### loop over perioder
for (BiWeek in BiWeeks) {
  cat("Processing yearXweek:", BiWeek, "\n")
  BiWeek_data <- input_matrix[input_matrix$BiWeek == BiWeek, ]
  
  # Clean and filter tokens
  articles_raw <- setNames(
    lapply(seq_len(nrow(BiWeek_data)), function(i) {
      tokens <- as.character(BiWeek_data[i, 2:(ncol(BiWeek_data)-2)])
      tokens <- tokens[tokens != "" & tokens != "NA" & !is.na(tokens)]
      tokens
    }),
    BiWeek_data[, 1]
  )
  
  j_filtered <- j_filtered_timeline
  tokens_in_week <- unlist(articles_raw)
  token_week_counts <- table(tokens_in_week)
  valid_tokens_for_week <- names(token_week_counts[token_week_counts >= IW])
  j_filtered <- intersect(j_filtered, valid_tokens_for_week)
  
  # Filter articles to only include j_filtered tokens
  articles <- lapply(articles_raw, function(article) {
    article[article %in% j_filtered]
  })
  
  num_tokens <- length(j_filtered)
  num_articles <- length(articles)
  print(num_tokens)
  
  # Split into DR and TV2 groups
  DR_articles <- articles[grep("^DR_", names(articles), ignore.case = TRUE)]
  TV2_articles <- articles[grep("^TV2_", names(articles), ignore.case = TRUE)]
  
  print("start q/p")
  
  # Compute total token counts
  total_DR_tokens <- sum(lengths(DR_articles))
  total_TV2_tokens <- sum(lengths(TV2_articles))
  
  # Token frequency tables
  tab_DR <- if (total_DR_tokens > 0) table(unlist(DR_articles)) else NULL
  tab_TV2 <- if (total_TV2_tokens > 0) table(unlist(TV2_articles)) else NULL
  
  # Frequency vectors (only j_filtered tokens)
  counts_DR <- if (!is.null(tab_DR)) as.numeric(tab_DR[j_filtered]) else rep(0, num_tokens)
  counts_TV2 <- if (!is.null(tab_TV2)) as.numeric(tab_TV2[j_filtered]) else rep(0, num_tokens)
  counts_DR[is.na(counts_DR)] <- 0
  counts_TV2[is.na(counts_TV2)] <- 0
  
  # Create q_matrix with two columns (DR and TV2 group-level frequencies)
  q_vector_DR <- if (total_DR_tokens > 0) counts_DR / total_DR_tokens else rep(0, num_tokens)
  q_vector_TV2 <- if (total_TV2_tokens > 0) counts_TV2 / total_TV2_tokens else rep(0, num_tokens)
  
  q_matrix <- cbind(q_vector_DR, q_vector_TV2)
  colnames(q_matrix) <- c("dr_group", "tv2_group")
  rownames(q_matrix) <- j_filtered
  
  # Create p_matrix of same shape
  frac_DR <- if (total_DR_tokens > 0) counts_DR / total_DR_tokens else rep(0, num_tokens)
  frac_TV2 <- if (total_TV2_tokens > 0) counts_TV2 / total_TV2_tokens else rep(0, num_tokens)
  
  denom <- frac_DR + frac_TV2
  p_vals <- ifelse(denom > 0, frac_DR / denom, 0)
  
  p_matrix <- matrix(p_vals, nrow = num_tokens, ncol = 2)
  colnames(p_matrix) <- c("dr_group", "tv2_group")
  rownames(p_matrix) <- j_filtered
  
  print("calc X")
  X <- compute_X_matrix(q_matrix, p_matrix)
  
  tau <- length(articles)
  subsample_size <- floor(tau / subsample_size_denominator)
  pi_subsamples <- numeric(K)
  tau_k_vec <- numeric(K)
  
  for (k in 1:K) {
    indices <- sample(seq_along(articles), size = subsample_size, replace = FALSE)
    articles_sub <- articles[indices]
    pi_subsamples[k] <- compute_pi(articles_sub, j_filtered)
    tau_k_vec[k] <- length(articles_sub)
  }
  
  avg_pi_sub <- mean(pi_subsamples)
  Q <- sqrt(tau_k_vec) * (pi_subsamples - avg_pi_sub)
  Q_sorted <- sort(Q)
  Q_90 <- Q_sorted[CI_top_limit]
  Q_11 <- Q_sorted[CI_bottom_limit]
  
  pi_full <- X
  CI_lower <- pi_full - Q_90 / sqrt(tau)
  CI_upper <- pi_full - Q_11 / sqrt(tau)

  ### start random assignment ##########################################################
  
  dr_rows <- BiWeek_data[grep("^dr", BiWeek_data[, 1], ignore.case = TRUE), ]
  tv2_rows <- BiWeek_data[grep("^tv2", BiWeek_data[, 1], ignore.case = TRUE), ]
  dr_rows <- dr_rows[sample(nrow(dr_rows)), ]
  tv2_rows <- tv2_rows[sample(nrow(tv2_rows)), ]
  week_data_random <- rbind(dr_rows, tv2_rows)
  week_data_random[, 1] <- sample(week_data_random[, 1])
  # Clean and filter tokens
  articles_raw <- setNames(
    lapply(seq_len(nrow(week_data_random)), function(i) {
      tokens <- as.character(week_data_random[i, 2:(ncol(week_data_random)-2)])
      tokens <- tokens[tokens != "" & tokens != "NA" & !is.na(tokens)]
      tokens
    }),
    week_data_random[, 1]
  )
  
  j_filtered <- j_filtered_timeline
  tokens_in_week <- unlist(articles_raw)
  token_week_counts <- table(tokens_in_week)
  valid_tokens_for_week <- names(token_week_counts[token_week_counts >= IW])
  j_filtered <- intersect(j_filtered, valid_tokens_for_week)
  
  # Filter articles to only include j_filtered tokens
  articles <- lapply(articles_raw, function(article) {
    article[article %in% j_filtered]
  })
  
  num_tokens <- length(j_filtered)
  num_articles <- length(articles)
  
  # Split into DR and TV2 groups
  DR_articles <- articles[grep("^DR_", names(articles), ignore.case = TRUE)]
  TV2_articles <- articles[grep("^TV2_", names(articles), ignore.case = TRUE)]
  
  # Compute total token counts
  total_DR_tokens <- sum(lengths(DR_articles))
  total_TV2_tokens <- sum(lengths(TV2_articles))
  
  # Token frequency tables
  tab_DR <- if (total_DR_tokens > 0) table(unlist(DR_articles)) else NULL
  tab_TV2 <- if (total_TV2_tokens > 0) table(unlist(TV2_articles)) else NULL
  
  # Frequency vectors (only j_filtered tokens)
  counts_DR <- if (!is.null(tab_DR)) as.numeric(tab_DR[j_filtered]) else rep(0, num_tokens)
  counts_TV2 <- if (!is.null(tab_TV2)) as.numeric(tab_TV2[j_filtered]) else rep(0, num_tokens)
  counts_DR[is.na(counts_DR)] <- 0
  counts_TV2[is.na(counts_TV2)] <- 0
  
  # Create q_matrix with two columns (DR and TV2 group-level frequencies)
  q_vector_DR <- if (total_DR_tokens > 0) counts_DR / total_DR_tokens else rep(0, num_tokens)
  q_vector_TV2 <- if (total_TV2_tokens > 0) counts_TV2 / total_TV2_tokens else rep(0, num_tokens)
  
  q_matrix <- cbind(q_vector_DR, q_vector_TV2)
  colnames(q_matrix) <- c("dr_group", "tv2_group")
  rownames(q_matrix) <- j_filtered
  
  # Create p_matrix of same shape
  frac_DR <- if (total_DR_tokens > 0) counts_DR / total_DR_tokens else rep(0, num_tokens)
  frac_TV2 <- if (total_TV2_tokens > 0) counts_TV2 / total_TV2_tokens else rep(0, num_tokens)
  
  denom <- frac_DR + frac_TV2
  p_vals <- ifelse(denom > 0, frac_DR / denom, 0)
  
  p_matrix <- matrix(p_vals, nrow = num_tokens, ncol = 2)
  colnames(p_matrix) <- c("dr_group", "tv2_group")
  rownames(p_matrix) <- j_filtered
  
  print("calc Xr")
  Xr <- compute_X_matrix(q_matrix, p_matrix)
  
  tau <- length(articles)
  subsample_size <- floor(tau / subsample_size_denominator)
  pi_subsamples <- numeric(K)
  tau_k_vec <- numeric(K)
  
  for (k in 1:K) {
    indices <- sample(seq_along(articles), size = subsample_size, replace = FALSE)
    articles_sub <- articles[indices]
    pi_subsamples[k] <- compute_pi(articles_sub, j_filtered)
    tau_k_vec[k] <- length(articles_sub)
  }
  
  avg_pi_sub <- mean(pi_subsamples)
  Qr <- sqrt(tau_k_vec) * (pi_subsamples - avg_pi_sub)
  Qr_sorted <- sort(Qr)
  Q_90r <- Qr_sorted[CI_top_limit]
  Q_11r <- Qr_sorted[CI_bottom_limit]
  pi_fullr <- Xr
  CIr_lower <- pi_fullr - Q_90r / sqrt(tau)
  CIr_upper <- pi_fullr - Q_11r / sqrt(tau)
  print(X)
  MLE_results[BiWeek,] <- c(BiWeek, X, CI_lower, CI_upper, Xr, CIr_lower, CIr_upper)         # export 
}


##### save: 
save(MLE_results, file="MLE_results_apr18_all.Rdata")
#load("LO_results_samlet_ny_apr.Rdata")

MLE_results <- as.data.frame(MLE_results)
MLE_results$BiWeek <- factor(MLE_results$BiWeek, levels = unique(MLE_results$BiWeek), ordered = TRUE)

ggplot(MLE_results, aes(x = BiWeek)) +  
  # Orange ribbons, lines, and points for the second method (overlay)
  geom_ribbon(aes(ymin = as.numeric(CI_lower), ymax = as.numeric(CI_upper), group = 1), fill = "blue1", alpha = 0.2, data = MLE_results) +
  geom_line(aes(y = as.numeric(X), group = 1), color = "blue1", size = 1, data = MLE_results) +
  geom_point(aes(y = as.numeric(X)), size = 2, colour = "blue3", data = MLE_results) +
  geom_ribbon(aes(ymin = as.numeric(CIr_lower), ymax = as.numeric(CIr_upper), group = 1), fill = "grey50", alpha = 0.2, data = MLE_results) +
  geom_line(aes(y = as.numeric(Xr), group = 1), color = "grey50", size = 1, data = MLE_results) +
  geom_point(aes(y = as.numeric(Xr)), size = 2, colour = "grey55", data = MLE_results) +
  
  
  # Vertical line separating years
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2015-26"), 
             linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2016-26"), 
             linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2017-26"), 
             linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2018-26"), 
             linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2019-25"), ## obs her er 2019-26 før uge 1?
             linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2020-25"), 
             linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2021-26"), 
             linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2022-26"), 
             linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = which(levels(MLE_results$BiWeek) == "2023-26"), 
             linetype = "dashed", color = "black", size = 1) +
  
  # Labels and axis formatting
  labs(x = "Week", y = "Polarisation", title = "LO-estimates with Confidence Intervals. J50x10x05") +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +  # Limit to 2 decimal places
  theme_minimal()  


