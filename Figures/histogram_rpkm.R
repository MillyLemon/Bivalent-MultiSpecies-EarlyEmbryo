# Load required packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Define processing function
process_histone_data <- function(input_file, histone_mark) {
  # Read data
  data <- read.table(input_file, header = TRUE)
  
  # 1. Data preparation and log2 RPKM calculation -------------------------------
  
  calculate_log2_rpkm <- function(count_matrix, peak_lengths) {
    total_reads_m <- colSums(count_matrix) / 1e6
    lengths_kb <- peak_lengths / 1000
    rpkm_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix))
    
    for(i in 1:ncol(count_matrix)) {
      rpkm_matrix[, i] <- count_matrix[, i] / (lengths_kb * total_reads_m[i])
    }
    
    log2_rpkm_matrix <- log2(rpkm_matrix + 1)
    colnames(log2_rpkm_matrix) <- colnames(count_matrix)
    return(log2_rpkm_matrix)
  }
  
  # Extract peak_id and peak lengths
  peak_ids <- data[, 1]
  peak_lengths <- sapply(strsplit(peak_ids, "_"), function(x) {
    as.numeric(x[3]) - as.numeric(x[2])
  })
  
  # Extract sample counts matrix (exclude first column peak_id)
  count_matrix <- data[, -1]
  
  # Calculate log2 RPKM
  log2_rpkm <- calculate_log2_rpkm(count_matrix, peak_lengths)
  
  # Convert log2 RPKM back to data frame and add peak_id
  log2_rpkm_df <- as.data.frame(log2_rpkm)
  log2_rpkm_df$peak_id <- peak_ids
  
  # 2. Sample-level distribution analysis -------------------------------------
  
  # Convert to long format for sample-level analysis
  log2_rpkm_long <- melt(log2_rpkm_df, id.vars = "peak_id", 
                         variable.name = "sample", value.name = "log2_rpkm")
  
  # Create non-zero dataset (sample level)
  log2_rpkm_nonzero <- log2_rpkm_long[log2_rpkm_long$log2_rpkm > 0, ]
  
  # Clean sample names
  log2_rpkm_nonzero$sample_clean <- gsub("^X", "", log2_rpkm_nonzero$sample)
  log2_rpkm_nonzero$stage <- sapply(strsplit(log2_rpkm_nonzero$sample_clean, "_"), function(x) x[1])
  log2_rpkm_nonzero$sample_pretty <- sapply(strsplit(log2_rpkm_nonzero$sample_clean, "_"), function(x) {
    stage <- x[1]
    mark <- x[2]
    rep <- x[3]
    paste0(toupper(stage), " ", mark, " Rep", gsub("rep", "", rep))
  })
  
  stage_order <- c("2cell", "4cell", "8cell", "morula", "ICM", "TE")
  log2_rpkm_nonzero$stage <- factor(log2_rpkm_nonzero$stage, levels = stage_order)
  
  # 3. Calculate distribution statistics for each sample ----------------------
  
  # Calculate number of peaks for each sample at different thresholds
  calculate_sample_threshold_stats <- function(nonzero_data, thresholds) {
    sample_stats <- list()
    
    for(sample_name in unique(nonzero_data$sample_pretty)) {
      sample_data <- nonzero_data[nonzero_data$sample_pretty == sample_name, ]
      total_peaks_in_sample <- nrow(sample_data)
      
      threshold_counts <- sapply(thresholds, function(thresh) {
        sum(sample_data$log2_rpkm > thresh)
      })
      
      sample_stats[[sample_name]] <- list(
        sample = sample_name,
        total_nonzero = total_peaks_in_sample,
        threshold_counts = threshold_counts,
        threshold_percents = round(threshold_counts / total_peaks_in_sample * 100, 1)
      )
    }
    return(sample_stats)
  }
  
  # Add threshold 0.5 as requested
  thresholds <- c(0.5, 1, 1.5, 2, 2.5, 3)
  sample_stats <- calculate_sample_threshold_stats(log2_rpkm_nonzero, thresholds)
  
  # 4. Create improved faceted histograms for samples -------------------------
  
  # Create a function to generate individual sample plots with statistics and letters
  # Create a function to generate individual sample plots with statistics and letters
  # Create a function to generate individual sample plots with statistics and letters
  create_sample_plots_with_stats <- function(nonzero_data, sample_stats, thresholds, total_peaks, free_y = TRUE, start_letter = "A") {
  plots <- list()
  
  # 修复：正确生成字母序列
  all_letters <- LETTERS  # 预定义的字母序列 A-Z
  start_idx <- which(all_letters == start_letter)  # 找到起始字母的索引
  letters <- all_letters[start_idx:(start_idx + length(sample_stats) - 1)]  # 生成连续字母
  
  for(i in 1:length(names(sample_stats))) {
    sample_name <- names(sample_stats)[i]
    current_letter <- letters[i]
    sample_data <- nonzero_data[nonzero_data$sample_pretty == sample_name, ]
    stats <- sample_stats[[sample_name]]
    
    # Calculate percentage of total peaks
    total_nonzero_pct <- round(stats$total_nonzero / total_peaks * 100, 1)
    
    # Create statistics text for this sample
    stats_text <- paste0(
      "Nonzero: ", format(stats$total_nonzero, big.mark = ","), 
      " (", total_nonzero_pct, "%)\n\n",
      "Threshold Counts:\n"
    )
    
    # Add threshold counts
    for(j in 1:length(thresholds)) {
      stats_text <- paste0(stats_text, 
                          ">", thresholds[j], ": ", 
                          format(stats$threshold_counts[j], big.mark = ","), 
                          " (", stats$threshold_percents[j], "%)\n")
    }
    
    # Create plot for this sample
    p <- ggplot(sample_data, aes(x = log2_rpkm)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8, color = "black") +
      # Add threshold lines
      geom_vline(xintercept = thresholds, 
                 linetype = "dashed", 
                 color = c("darkorange", "orange", "yellow", "red", "purple", "brown"),
                 size = 0.8, alpha = 0.7) +
      # Add statistics text
      annotate("text", x = max(sample_data$log2_rpkm) * 0.7, 
               y = Inf, 
               label = stats_text, 
               hjust = 0, vjust = 1.2, size = 2.8, family = "mono", lineheight = 0.8) +
      labs(title = paste0(current_letter, ". ", sample_name),
           x = "log2(RPKM+1)",
           y = "Number of Peaks") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6))
    
    # Set y-axis limits if not free
    if(!free_y) {
      # Calculate global y max across all samples
      all_counts <- sapply(names(sample_stats), function(sname) {
        sdata <- nonzero_data[nonzero_data$sample_pretty == sname, ]
        max(hist(sdata$log2_rpkm, breaks = 50, plot = FALSE)$counts)
      })
      global_y_max <- max(all_counts) * 1.05
      
      p <- p + ylim(0, global_y_max)
    }
    
    plots[[sample_name]] = p
  }
  return(plots)
}
  # Calculate total peaks
  total_peaks <- nrow(log2_rpkm_df)
  
  # Create individual sample plots with statistics - free y-axis (start from A)
  sample_plots_free <- create_sample_plots_with_stats(log2_rpkm_nonzero, sample_stats, thresholds, total_peaks, free_y = TRUE, start_letter = "A")
  
  # Combine all sample plots into a grid - free y-axis
  p_facet_free <- do.call(gridExtra::grid.arrange, c(sample_plots_free, ncol = 3))
  
  # Print and save the faceted plot with free y-axis as PDF and PNG
  print(p_facet_free)
  ggsave(paste0("sample_histograms_free_y_", histone_mark, ".pdf"), p_facet_free, 
         width = 16, height = 12)
  ggsave(paste0("sample_histograms_free_y_", histone_mark, ".png"), p_facet_free, 
         width = 16, height = 12, dpi = 300)
  
  # Create individual sample plots with statistics - fixed y-axis (continue from where free_y left off)
 # 在 create_sample_plots_with_stats 函数调用之后
  n_free_plots <- length(sample_plots_free)
  next_letter <- LETTERS[which(LETTERS == "A") + n_free_plots]  # 继续下一个字母 
  sample_plots_fixed <- create_sample_plots_with_stats(log2_rpkm_nonzero, sample_stats, thresholds, total_peaks, free_y = FALSE, start_letter = next_letter)
  
  # Combine all sample plots into a grid - fixed y-axis
  p_facet_fixed <- do.call(gridExtra::grid.arrange, c(sample_plots_fixed, ncol = 3))
  
  # Print and save the faceted plot with fixed y-axis as PDF and PNG
  print(p_facet_fixed)
  ggsave(paste0("sample_histograms_fixed_y_", histone_mark, ".pdf"), p_facet_fixed, 
         width = 16, height = 12)
  ggsave(paste0("sample_histograms_fixed_y_", histone_mark, ".png"), p_facet_fixed, 
         width = 16, height = 12, dpi = 300)
  
  # 5. Create overall distribution histogram (all samples combined) -----------
  
  # Calculate overall statistics
  total_nonzero_observations <- nrow(log2_rpkm_nonzero)
  
  # Calculate unique peaks (non-redundant count) for each threshold
  calculate_unique_peaks_above_threshold <- function(log2_rpkm_df, thresholds) {
    # Remove peak_id column for calculation
    expr_matrix <- as.matrix(log2_rpkm_df[, -ncol(log2_rpkm_df)])
    
    unique_peaks <- sapply(thresholds, function(thresh) {
      # For each peak, check if it exceeds threshold in any sample
      peak_max <- apply(expr_matrix, 1, max)
      sum(peak_max > thresh)
    })
    
    return(unique_peaks)
  }
  
  # Calculate both total observations and unique peaks for each threshold
  overall_threshold_stats <- sapply(thresholds, function(thresh) {
    sum(log2_rpkm_nonzero$log2_rpkm > thresh)
  })
  
  unique_peaks_stats <- calculate_unique_peaks_above_threshold(log2_rpkm_df, thresholds)
  
  # Create statistics text without black box
  stats_text <- paste0(
    "Dataset Information:\n",
    "Total peaks: ", format(total_peaks, big.mark = ","), "\n",
    "Total non-zero observations: ", format(total_nonzero_observations, big.mark = ","), "\n",
    "Average peaks per sample: ", round(total_nonzero_observations / total_peaks, 1), "\n\n",
    "Threshold Statistics:\n"
  )
  
  # Add threshold statistics with both observation counts and unique peaks
  for(i in 1:length(thresholds)) {
    stats_text <- paste0(stats_text,
                        "log2(RPKM+1) > ", thresholds[i], ":\n",
                        "  Observations: ", format(overall_threshold_stats[i], big.mark = ","), "\n",
                        "  Unique peaks: ", format(unique_peaks_stats[i], big.mark = ","), 
                        " (", round(unique_peaks_stats[i] / total_peaks * 100, 1), "%)\n\n")
  }
  
  # Pre-calculate histogram data to determine y-axis max
  hist_temp <- ggplot(log2_rpkm_nonzero, aes(x = log2_rpkm)) + geom_histogram(bins = 60)
  hist_data <- ggplot_build(hist_temp)
  max_count <- max(hist_data$data[[1]]$count)
  
  # Determine the next letter after fixed_y plots
  total_plots <- length(sample_plots_free) + length(sample_plots_fixed)
  overall_letter <- LETTERS[total_plots + 1]
  
  # Create overall histogram
  p_overall <- ggplot(log2_rpkm_nonzero, aes(x = log2_rpkm)) +
    geom_histogram(bins = 60, fill = "lightblue", alpha = 0.7, color = "black") +
    # Add threshold lines
    geom_vline(xintercept = thresholds, 
               linetype = "dashed", 
               color = c("darkorange", "orange", "yellow", "red", "purple", "brown"),
               size = 1, alpha = 0.7) +
    # Add statistics text without box
    annotate("text", x = max(log2_rpkm_nonzero$log2_rpkm) * 0.65, 
             y = max_count * 0.95, 
             label = stats_text, 
             hjust = 0, vjust = 1, size = 3.2, family = "mono", lineheight = 0.85) +
    labs(title = paste0(overall_letter, ". H3", histone_mark, " - Combined Distribution of Non-zero log2(RPKM+1)"),
         subtitle = "Distribution of log2(RPKM+1) values across all samples for PCA threshold determination",
         x = "log2(RPKM+1)",
         y = "Number of Observations") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6))
  
  print(p_overall)
  ggsave(paste0("overall_distribution_", histone_mark, ".pdf"), p_overall, 
         width = 14, height = 8)
  ggsave(paste0("overall_distribution_", histone_mark, ".png"), p_overall, 
         width = 14, height = 8, dpi = 300)
  
  # 6. Create sample-level threshold statistics table -------------------------
  
  # Create detailed sample statistics table
  sample_summary <- data.frame(
    Sample = names(sample_stats),
    Total_Nonzero_Peaks = sapply(sample_stats, function(x) x$total_nonzero),
    Total_Nonzero_Pct = round(sapply(sample_stats, function(x) x$total_nonzero) / total_peaks * 100, 1),
    stringsAsFactors = FALSE
  )
  
  # Add counts for each threshold
  for(i in 1:length(thresholds)) {
    threshold_col <- paste0("Threshold_", thresholds[i])
    sample_summary[[threshold_col]] <- sapply(sample_stats, function(x) x$threshold_counts[i])
    sample_summary[[paste0(threshold_col, "_Pct")]] <- sapply(sample_stats, function(x) x$threshold_percents[i])
  }
  
  print("Threshold statistics for each sample:")
  print(sample_summary)
  
  # 7. PCA threshold recommendation -------------------------------------------
  
  cat(paste0("\n=== H3", histone_mark, " PCA Threshold Recommendation ===\n"))
  
  # Calculate number of peaks expressed in at least 1, 2, and 3 samples for different thresholds
  calculate_pca_candidates <- function(log2_rpkm_df, thresholds, min_samples = 1) {
    # Remove peak_id column
    expr_matrix <- as.matrix(log2_rpkm_df[, -ncol(log2_rpkm_df)])
    
    pca_stats <- data.frame(
      Threshold = thresholds,
      Peaks_For_PCA = 0,
      Percent_Of_Total = 0
    )
    
    for(i in 1:length(thresholds)) {
      # Calculate how many samples each peak exceeds threshold in
      n_samples_above <- rowSums(expr_matrix > thresholds[i])
      # Count peaks that exceed threshold in at least min_samples
      pca_stats$Peaks_For_PCA[i] <- sum(n_samples_above >= min_samples)
      pca_stats$Percent_Of_Total[i] <- round(pca_stats$Peaks_For_PCA[i] / nrow(log2_rpkm_df) * 100, 1)
    }
    
    return(pca_stats)
  }
  
  # Calculate for different min_samples criteria
  pca_candidates_1sample <- calculate_pca_candidates(log2_rpkm_df, thresholds, min_samples = 1)
  pca_candidates_2samples <- calculate_pca_candidates(log2_rpkm_df, thresholds, min_samples = 2)
  pca_candidates_3samples <- calculate_pca_candidates(log2_rpkm_df, thresholds, min_samples = 3)
  
  # Create a comprehensive summary table
  pca_summary <- data.frame(
    Threshold = thresholds,
    Peaks_1Sample = pca_candidates_1sample$Peaks_For_PCA,
    Percent_1Sample = pca_candidates_1sample$Percent_Of_Total,
    Peaks_2Samples = pca_candidates_2samples$Peaks_For_PCA,
    Percent_2Samples = pca_candidates_2samples$Percent_Of_Total,
    Peaks_3Samples = pca_candidates_3samples$Peaks_For_PCA,
    Percent_3Samples = pca_candidates_3samples$Percent_Of_Total
  )
  
  # Output PCA candidate peaks statistics
  cat("Comprehensive PCA threshold statistics:\n")
  print(pca_summary)
  
  # Multiple recommendation strategies
  cat("\n=== PCA Threshold Recommendation Strategies ===\n")
  
  # Strategy 1: At least 50K peaks in ≥2 samples
  best_threshold_strategy1 <- NULL
  for(thresh in thresholds) {
    candidates <- pca_candidates_2samples[pca_candidates_2samples$Threshold == thresh, "Peaks_For_PCA"]
    if(candidates >= 50000) {
      best_threshold_strategy1 <- thresh
      break
    }
  }
  
  if(is.null(best_threshold_strategy1)) {
    best_threshold_strategy1 <- thresholds[which.max(pca_candidates_2samples$Peaks_For_PCA)]
  }
  
  # Strategy 2: Balance between sensitivity and specificity (≥2 samples)
  # Find threshold that maximizes peaks while maintaining reasonable stringency
  candidate_thresholds <- thresholds[pca_candidates_2samples$Peaks_For_PCA >= 30000]
  if(length(candidate_thresholds) > 0) {
    best_threshold_strategy2 <- min(candidate_thresholds)  # Most sensitive while maintaining 30K peaks
  } else {
    best_threshold_strategy2 <- thresholds[which.max(pca_candidates_2samples$Peaks_For_PCA)]
  }
  
  # Strategy 3: Based on percentage of total peaks (aim for 20-40% retention)
  candidate_pcts <- pca_candidates_2samples$Percent_Of_Total
  ideal_thresholds <- thresholds[candidate_pcts >= 20 & candidate_pcts <= 40]
  if(length(ideal_thresholds) > 0) {
    best_threshold_strategy3 <- min(ideal_thresholds)  # Most sensitive within ideal range
  } else {
    # Find closest to 30%
    closest_idx <- which.min(abs(candidate_pcts - 30))
    best_threshold_strategy3 <- thresholds[closest_idx]
  }
  
  # Final recommendation (prefer strategy 1 if possible, otherwise strategy 2)
  if(!is.null(best_threshold_strategy1)) {
    final_recommendation <- best_threshold_strategy1
    recommendation_reason <- "Retains ≥50K peaks in ≥2 samples"
  } else {
    final_recommendation <- best_threshold_strategy2
    recommendation_reason <- "Balances sensitivity and specificity"
  }
  
  cat(paste0("\nStrategy 1 (≥50K peaks in ≥2 samples): log2(RPKM+1) > ", best_threshold_strategy1, "\n"))
  cat(paste0("  Would retain ", pca_candidates_2samples[pca_candidates_2samples$Threshold == best_threshold_strategy1, "Peaks_For_PCA"], 
             " peaks (", pca_candidates_2samples[pca_candidates_2samples$Threshold == best_threshold_strategy1, "Percent_Of_Total"], "%)\n"))
  
  cat(paste0("Strategy 2 (Balance sensitivity/specificity): log2(RPKM+1) > ", best_threshold_strategy2, "\n"))
  cat(paste0("  Would retain ", pca_candidates_2samples[pca_candidates_2samples$Threshold == best_threshold_strategy2, "Peaks_For_PCA"], 
             " peaks (", pca_candidates_2samples[pca_candidates_2samples$Threshold == best_threshold_strategy2, "Percent_Of_Total"], "%)\n"))
  
  cat(paste0("Strategy 3 (20-40% retention): log2(RPKM+1) > ", best_threshold_strategy3, "\n"))
  cat(paste0("  Would retain ", pca_candidates_2samples[pca_candidates_2samples$Threshold == best_threshold_strategy3, "Peaks_For_PCA"], 
             " peaks (", pca_candidates_2samples[pca_candidates_2samples$Threshold == best_threshold_strategy3, "Percent_Of_Total"], "%)\n"))
  
  cat(paste0("\nFinal Recommendation: log2(RPKM+1) > ", final_recommendation, " (", recommendation_reason, ")\n"))
  
  # Also show the result for at least 1 and 3 samples with final recommendation
  candidates_1sample <- pca_candidates_1sample[pca_candidates_1sample$Threshold == final_recommendation, "Peaks_For_PCA"]
  percent_1sample <- pca_candidates_1sample[pca_candidates_1sample$Threshold == final_recommendation, "Percent_Of_Total"]
  candidates_2samples <- pca_candidates_2samples[pca_candidates_2samples$Threshold == final_recommendation, "Peaks_For_PCA"]
  percent_2samples <- pca_candidates_2samples[pca_candidates_2samples$Threshold == final_recommendation, "Percent_Of_Total"]
  candidates_3samples <- pca_candidates_3samples[pca_candidates_3samples$Threshold == final_recommendation, "Peaks_For_PCA"]
  percent_3samples <- pca_candidates_3samples[pca_candidates_3samples$Threshold == final_recommendation, "Percent_Of_Total"]
  
  cat(paste0("\nWith this threshold:\n"))
  cat(paste0("  - ", candidates_1sample, " peaks (", percent_1sample, "%) are expressed in ≥1 sample\n"))
  cat(paste0("  - ", candidates_2samples, " peaks (", percent_2samples, "%) are expressed in ≥2 samples\n"))
  cat(paste0("  - ", candidates_3samples, " peaks (", percent_3samples, "%) are expressed in ≥3 samples\n"))
  
  # 8. Save filtered PCA matrix -----------------------------------------------
  
  # Apply recommended threshold to filter peaks
  expr_matrix <- as.matrix(log2_rpkm_df[, -ncol(log2_rpkm_df)])
  n_samples_above <- rowSums(expr_matrix > final_recommendation)
  pca_peaks <- log2_rpkm_df[n_samples_above >= 2, ]
  
  # Save PCA matrix
  write.csv(pca_peaks, paste0("pca_matrix_", histone_mark, "_threshold_", final_recommendation, ".csv"), 
            row.names = FALSE)
  
  cat(paste0("\nSaved PCA matrix: pca_matrix_", histone_mark, "_threshold_", final_recommendation, ".csv\n"))
  
  return(list(
    sample_stats = sample_stats,
    pca_summary = pca_summary,
    pca_candidates_1sample = pca_candidates_1sample,
    pca_candidates_2samples = pca_candidates_2samples,
    pca_candidates_3samples = pca_candidates_3samples,
    recommended_threshold = final_recommendation,
    pca_matrix = pca_peaks,
    unique_peaks_stats = data.frame(Threshold = thresholds, Unique_Peaks = unique_peaks_stats)
  ))
}

# Process K4 and K27 data
cat("Processing H3K4me3 data...\n")
k4_results <- process_histone_data("k4.coverage_matrix.txt", "K4me3")

cat("\nProcessing H3K27me3 data...\n")  
k27_results <- process_histone_data("k27.coverage_matrix.txt", "K27me3")

# Output final summary
cat("\n=== Final PCA Analysis Summary ===\n")
summary_df <- data.frame(
  Histone_Mark = c("H3K4me3", "H3K27me3"),
  Recommended_Threshold = c(k4_results$recommended_threshold, k27_results$recommended_threshold),
  Peaks_For_PCA = c(nrow(k4_results$pca_matrix), nrow(k27_results$pca_matrix)),
  Percent_Retained = c(
    k4_results$pca_candidates_2samples[k4_results$pca_candidates_2samples$Threshold == k4_results$recommended_threshold, "Percent_Of_Total"],
    k27_results$pca_candidates_2samples[k27_results$pca_candidates_2samples$Threshold == k27_results$recommended_threshold, "Percent_Of_Total"]
  ),
  Peaks_In_1Sample = c(
    k4_results$pca_candidates_1sample[k4_results$pca_candidates_1sample$Threshold == k4_results$recommended_threshold, "Peaks_For_PCA"],
    k27_results$pca_candidates_1sample[k27_results$pca_candidates_1sample$Threshold == k27_results$recommended_threshold, "Peaks_For_PCA"]
  ),
  Peaks_In_3Samples = c(
    k4_results$pca_candidates_3samples[k4_results$pca_candidates_3samples$Threshold == k4_results$recommended_threshold, "Peaks_For_PCA"],
    k27_results$pca_candidates_3samples[k27_results$pca_candidates_3samples$Threshold == k27_results$recommended_threshold, "Peaks_For_PCA"]
  )
)

print(summary_df)

cat("\nProcessing completed! Generated files:\n")
cat("For H3K4me3:\n")
cat("- sample_histograms_free_y_K4me3.pdf/.png: Faceted histograms with free y-axis\n")
cat("- sample_histograms_fixed_y_K4me3.pdf/.png: Faceted histograms with fixed y-axis\n")
cat("- overall_distribution_K4me3.pdf/.png: Combined distribution\n")
cat("- pca_matrix_K4me3_threshold_X.csv: PCA matrix\n\n")

cat("For H3K27me3:\n")
cat("- sample_histograms_free_y_K27me3.pdf/.png: Faceted histograms with free y-axis\n")
cat("- sample_histograms_fixed_y_K27me3.pdf/.png: Faceted histograms with fixed y-axis\n")
cat("- overall_distribution_K27me3.pdf/.png: Combined distribution\n")
cat("- pca_matrix_K27me3_threshold_X.csv: PCA matrix\n")
