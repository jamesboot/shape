#!/usr/bin/env Rscript

# Script to plot SHAPE data
# H. V. Mears
# Based on clipplotr by A. M. Chakrabarti and C. Capitanchik
# Written with help from chatGPT and tidied by Claude AI

# ==========
# Preamble
# ==========

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(zoo)
  library(patchwork)
  library(scales)
  library(ggforce)
  library(ggbreak)
  library(stringr)
  library(data.table)
  library(ggpointdensity)
  library(viridis)
})

# ==========
# Command Line Options
# ==========

option_list <- list(
  make_option(c("-d", "--dotplot"), type = "character", help = "Input .dp file"),
  make_option(c("-x", "--reactivity"), type = "character", help = "Input reactivity .wig file(s), comma-separated if multiple"),
  make_option(c("-s", "--shannon"), type = "character", help = "Input shannon .wig file(s), comma-separated if multiple"),
  make_option(c("-o", "--output"), type = "character", help = "Output PDF filename"),
  make_option(c("-r", "--region"), type = "character", default = NULL, help = "Region to plot (format: start:end), e.g. '100:500')"),
  make_option(c("--qc"), action = "store_true", default = FALSE, help = "If set, plot additional QC box and whisker plot for per base reactivity"),
  make_option(c("--rc_files"), type = "character", action = "store", help = "Comma-separated list of rc files for SHAPE or DMS-modified samples"),
  make_option(c("--rc_controls"), type = "character", action = "store", help = "Comma-separated list of rc files for controls"),
  make_option(c("--qc_colors"), type = "character", help = "Colors control and treated conditions in QC plots, comma delimited"),
  make_option(c("--qc_max_y"), type = "double", default = 0.15, help = "Fixed maximum y-axis limit for QC plots (default: 0.15)"),
  make_option(c("--qc_auto_y"), action = "store_true", default = FALSE, help = "If set, auto scale y axis for QC plots"),
  make_option(c("--compare_replicates"), action = "store_true", default = FALSE, help = "If set, generate pairwise reactivity scatter plots"),
  make_option(c("--compare_log"), action = "store_true", default = FALSE, help = "If set, log-transform reactivities for comparison plots and linear regression"),
  make_option(c("--break_y"), action = "store_true", default = FALSE, help = "If set, break reactivity y axis at 2.5"),
  make_option(c("--limit_y"), action = "store_true", default = FALSE, help = "If set, limit reactivity y axis to 0-2.5"),
  make_option(c("--smoothing"), action = "store_true", default = FALSE, help = "If set, apply smoothing to the line plot"),
  make_option(c("--window"), type = "integer", default = 25, help = "Window size for smoothing if --smoothing is set [default: %default nt]"),
  make_option(c("--median"), action = "store_true", default = FALSE, help = "If set, normalize data to the median reactivity and Shannon entropy"),
  make_option(c("--colors"), type = "character", help = "Colors for reactivity and shannon entropy graphs, comma delimited"),
  make_option(c("--arc_colors"), type = "character", help = "Colors for arc plots, based on pairing probability, comma delimited (low = 10-30%, mid = 30-80%, high > 80%)"),
  make_option(c("-a", "--alpha"), type = "double", default = 0.8, help = "Transparency of arcs (0 = transparent, 1 = opaque) [default: %default]"),
  make_option(c("--arc_height"), type = "double", default = 0.3, help = "Scaling factor for arc heights [default: %default]"),
  make_option(c("", "--size_x"), type = "integer", default = 210, help = "Plot size in mm (x) [default: %default]"),
  make_option(c("", "--size_y"), type = "integer", default = 160, help = "Plot size in mm (y) [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ==========
# Validation Functions
# ==========

validate_inputs <- function(opt) {
  # Check for required parameters
  if (is.null(opt$region)) {
    stop("ERROR: region of interest must be supplied as e.g. 1000:2400")
  }
  
  if (is.null(opt$output)) {
    stop("ERROR: No output defined")
  }
  
  # Check if any data sources are provided
  if (is.null(opt$dotplot) && is.null(opt$reactivity) && is.null(opt$shannon)) {
    stop("ERROR: No input data provided. Please specify at least one of: --dotplot, --reactivity, or --shannon")
  }
  
  # Check QC requirements
  if (!is.null(opt$qc) && opt$qc) {
    if (is.null(opt$rc_files) || is.null(opt$rc_controls)) {
      stop("ERROR: QC requires --rc_files and --rc_controls to be provided")
    }
  }
  
  # Check if output file exists
  if (file.exists(opt$output)) {
    message("WARNING: Output file '", opt$output, "' exists and will be overwritten!")
  }
  
  # Log input information
  if (!is.null(opt$reactivity)) {
    message("Reactivity file specified: ", opt$reactivity)
  }
  
  if (!is.null(opt$shannon)) {
    message("Shannon entropy file specified: ", opt$shannon)
  }
  
  if (!is.null(opt$dotplot)) {
    message("Dotplot file specified: ", opt$dotplot)
  }
}

# Parse region
parse_region <- function(region_text) {
  coords <- as.integer(strsplit(region_text, ":")[[1]])
  if (length(coords) != 2 || any(is.na(coords))) {
    stop("Region format invalid. Use format like '100:500'")
  }
  list(start = coords[1], end = coords[2])
}

# ==========
# Data Processing Functions
# ==========

# Read .wig files
read_wig <- function(file_path, sample_id = NULL) {
  tryCatch({
    # Read the .wig file
    wig_lines <- readLines(file_path)
    data_lines <- wig_lines[!grepl("^track|^variableStep", wig_lines)]  # Remove headers
    
    # Read data into a data frame
    wig_data <- read.table(text = data_lines, col.names = c("start", "value"))
    
    # Assign sample ID
    wig_data$sample_id <- ifelse(is.null(sample_id), basename(file_path), sample_id)
    
    return(wig_data)
  }, error = function(e) {
    stop("Error reading wig file: ", file_path, " - ", e$message)
  })
}

# Process replicate comparisons:
# Step 1: Process the input data for replicates
process_replicates <- function(data_combined) {
  replicate_data <- data.frame()
  
  for (i in unique(data_combined$sample_id)) {
    reactivity_values <- data_combined[data_combined$sample_id == i, ]
    replicate_df <- data.frame(
      replicate = i,
      start = reactivity_values$start,
      reactivity = reactivity_values$value
    )
    replicate_data <- rbind(replicate_data, replicate_df)
  }
  
  return(replicate_data)
}

# Step 2: Generate all pairwise comparisons for the replicates
generate_comparisons <- function(replicate_data) {
  replicate_labels <- unique(replicate_data$replicate)
  comparisons <- data.frame()
  
  for (i in 1:(length(replicate_labels) - 1)) {
    for (j in (i + 1):length(replicate_labels)) {
      replicate_1_data <- subset(replicate_data, replicate == replicate_labels[i])
      replicate_2_data <- subset(replicate_data, replicate == replicate_labels[j])
      
      comparison_df <- merge(replicate_1_data, replicate_2_data, by = "start", suffixes = c("_1", "_2"))
      comparison_df <- comparison_df[!is.na(comparison_df$reactivity_1) & !is.na(comparison_df$reactivity_2), ]
      comparison_df$pair <- paste(replicate_labels[i], "vs", replicate_labels[j])
      
      comparisons <- rbind(comparisons, comparison_df)
    }
  }
  
  return(comparisons)
}

# Step 3: Prepare data for plotting on a log scale
prepare_plot_data <- function(comparisons, compare_log = TRUE) {
  total_before <- nrow(comparisons)
  
  # Filter out non-positive values only if log is needed
  if (compare_log) {
    comparisons <- comparisons %>%
      filter(reactivity_1 > 0, reactivity_2 > 0)
    
    filtered_out <- total_before - nrow(comparisons)
    if (filtered_out > 0) {
      warning(filtered_out, " data points removed due to zero or negative reactivity values before log transformation.")
    }
    
    comparisons$reactivity_1_transformed <- log10(comparisons$reactivity_1)
    comparisons$reactivity_2_transformed <- log10(comparisons$reactivity_2)
  } else {
    comparisons$reactivity_1_transformed <- comparisons$reactivity_1
    comparisons$reactivity_2_transformed <- comparisons$reactivity_2
  }
  
  return(comparisons)
}

# Step 4: New clean handler for all replicate comparison processing
handle_replicate_comparisons <- function(file_paths, region_range) {
  tryCatch({
    message("Processing replicate comparisons")
    
    # Read files
    file_list <- trimws(unlist(strsplit(file_paths, ",")))
    data_combined <- bind_rows(
      lapply(seq_along(file_list), function(idx) {
        read_wig(file_list[idx], sample_id = paste0("replicate ", idx))
      })
    )
    
    # Subset to region
    data_combined <- data_combined %>%
      filter(start >= region_range$start & start <= region_range$end)
    
    # Process replicate data
    replicate_data <- process_replicates(data_combined)
    
    # Generate comparisons
    comparisons <- generate_comparisons(replicate_data)
    
    return(comparisons)
  }, error = function(e) {
    stop("Error processing replicate comparisons: ", e$message)
  })
}

# Extract p value from linear regression model
overall_p <- function(lm_model) {
  f <- summary(lm_model)$fstatistic
  p_value <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  attributes(p_value) <- NULL
  return(p_value)
}


# Significance threshold bins (for linear regression)
get_significance <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

# Process dotplot data
process_dotplot <- function(file_path, region_range, arc_height_factor) {
  tryCatch({
    message("Processing dotplot data")
    
    dp_data <- read_tsv(
      file_path,
      skip = 2,
      col_names = c("i", "j", "score"),
      col_types = cols(
        i = col_double(),
        j = col_double(),
        score = col_double()
      )
    )
    
    # Clean and transform the data
    dp_data <- dp_data %>%
      filter(!is.na(score)) %>%
      mutate(
        prob = 10^(-score),
        bin = case_when(
          prob > 0.80 ~ "high",
          prob > 0.30 ~ "mid",
          prob > 0.10 ~ "low",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(bin)) %>%
      mutate(bin = factor(bin, levels = c("low", "mid", "high")))
    
    # Subset to region
    dp_data <- dp_data %>% 
      filter(i >= region_range$start & i <= region_range$end &
               j >= region_range$start & j <= region_range$end,
             i != j) # Remove self-arcs
    
    message("Number of arcs after filtering: ", nrow(dp_data))
    
    # Create Bezier control points for arcs
    bezier_data <- dp_data %>%
      rowwise() %>%
      mutate(
        x = list(c(i, (i + j) / 2, j)),
        y = list(c(0, arc_height_factor * abs(j - i), 0))
      ) %>%
      unnest(c(x, y)) %>%
      group_by(i, j, bin) %>%
      mutate(point_id = row_number()) %>%
      ungroup() %>%
      filter(point_id <= 3) %>%  # Ensure only 3 rows per arc
      mutate(arc_group = factor(interaction(i, j))) %>%
      arrange(factor(bin, levels = c("low", "mid", "high")))
    
    return(bezier_data)
  }, error = function(e) {
    stop("Error processing dotplot file: ", e$message)
  })
}

# Function to simplify polygons proportionally based on vertex count
simplify_polygon <- function(polygon, target_vertices = 5000) {
  num_vertices <- nrow(polygon)
  
  if (num_vertices <= target_vertices) return(polygon)
  
  idx <- round(seq(1, num_vertices, length.out = target_vertices))
  return(polygon[idx, ])
}

# Process measurement data (reactivity or shannon)
process_measurement_data <- function(file_paths, region_range, do_smoothing = FALSE, window_size = 25, do_median_normalize = FALSE, data_type = "reactivity") {
  tryCatch({
    message(paste("Processing", data_type, "data"))
    
    # Split comma-separated file paths
    file_list <- trimws(unlist(strsplit(file_paths, ",")))
    
    # Read all files
    data_combined <- bind_rows(
      lapply(seq_along(file_list), function(idx) {
        read_wig(file_list[idx], sample_id = paste0("replicate ", idx))
      })
    )
    
    # Apply median normalization if requested
    if (do_median_normalize && nrow(data_combined) > 0) {
      message("Normalizing to median...")
      median_value <- median(data_combined$value, na.rm = TRUE)
      data_combined$value <- data_combined$value / median_value
    }
    
    # Apply smoothing if requested
    if (do_smoothing) {
      message("Applying smoothing with window size ", window_size)
      
      data_combined <- data_combined %>%
        group_by(sample_id) %>%
        arrange(start) %>%
        mutate(smooth_value = zoo::rollmean(value, k = window_size, 
                                            fill = NA, align = "center", partial = TRUE)) %>%
        ungroup()
      
      # Calculate statistics on smoothed values
      data_summary <- data_combined %>%
        group_by(start) %>%
        summarise(
          mean_value = mean(smooth_value, na.rm = TRUE),
          sd_value   = sd(smooth_value, na.rm = TRUE),
          n          = sum(!is.na(smooth_value)),
          ci_lower   = ifelse(n > 1, 
                              mean_value - 1.96 * sd_value / sqrt(n), 
                              mean_value),
          ci_upper   = ifelse(n > 1, 
                              mean_value + 1.96 * sd_value / sqrt(n), 
                              mean_value),
          .groups    = "drop"
        ) %>%
        mutate(
          ci_lower = ifelse(is.na(ci_lower), mean_value, ci_lower),
          ci_upper = ifelse(is.na(ci_upper), mean_value, ci_upper),
          smooth_value = mean_value  # For consistency in plotting
        )
      
    } else {
      # Calculate statistics directly on raw values
      data_summary <- data_combined %>%
        group_by(start) %>%
        summarise(
          mean_value = mean(value, na.rm = TRUE),
          sd_value   = sd(value, na.rm = TRUE),
          n          = sum(!is.na(value)),
          ci_lower   = ifelse(n > 1, 
                              mean_value - 1.96 * sd_value / sqrt(n), 
                              mean_value),
          ci_upper   = ifelse(n > 1, 
                              mean_value + 1.96 * sd_value / sqrt(n), 
                              mean_value),
          .groups    = "drop"
        ) %>%
        mutate(
          sd_value = ifelse(is.na(sd_value), 0, sd_value),
          ci_lower = ifelse(is.na(ci_lower), mean_value, ci_lower),
          ci_upper = ifelse(is.na(ci_upper), mean_value, ci_upper),
          smooth_value = mean_value  # For consistency in plotting
        )
    }
    
    # Filter CI values to finite values AND inside region
    ribbons <- data_summary %>%
      filter(is.finite(ci_lower), is.finite(ci_upper),
             start >= region_range$start, start <= region_range$end)
    
    # Create ribbon polygon in correct order
    lower <- ribbons %>%
      arrange(start) %>%
      transmute(x = start, y = ci_lower)
    
    upper <- ribbons %>%
      arrange(desc(start)) %>%
      transmute(x = start, y = ci_upper)
    
    ribbon_polygon <- rbind(lower, upper)
    
    # Simplify polygon to reduce file size
    if (nrow(ribbon_polygon) > 0) {
      ribbons <- simplify_polygon(ribbon_polygon, target_vertices = 5000)
    } else {
      ribbons <- NULL
    }
    
    # Subset data to region
    data_summary <- data_summary %>%
      filter(start >= region_range$start & start <= region_range$end)
    
    return(list(
      data_summary = data_summary,
      ribbons = ribbons
    ))
  }, error = function(e) {
    stop(paste("Error processing", data_type, "data:", e$message))
  })
}

# Process QC data
process_qc_data <- function(sample_files, control_files) {
  tryCatch({
    message("Processing QC data")
    
    # Split file paths
    sample_rcs <- trimws(unlist(strsplit(sample_files, split = ",")))
    control_rcs <- trimws(unlist(strsplit(control_files, split = ",")))
    
    # Report counts
    message(sprintf("✓ %d replicate files provided for sample (treated) condition.", length(sample_rcs)))
    message(sprintf("✓ %d replicate files provided for control condition.", length(control_rcs)))
    
    # Check if files exist
    missing_files <- c(sample_rcs, control_rcs)[!file.exists(c(sample_rcs, control_rcs))]
    if (length(missing_files) > 0) {
      stop("Error: The following input files do not exist:\n", paste(missing_files, collapse = "\n"))
    }
    
    # Read RC files
    read_rc <- function(file, condition_label) {
      dt <- fread(file, skip = 1, col.names = c("base", "mutations", "coverage"))
      dt[, reactivity := mutations / coverage]
      dt[, sample := tools::file_path_sans_ext(basename(file))]
      dt[, condition := condition_label]
      return(dt)
    }
    
    # Process files
    sample_dt <- rbindlist(lapply(sample_rcs, read_rc, condition_label = "treated"))
    control_dt <- rbindlist(lapply(control_rcs, read_rc, condition_label = "control"))
    all_dt <- rbindlist(list(sample_dt, control_dt))
    
    return(all_dt)
  }, error = function(e) {
    stop("Error processing QC data: ", e$message)
  })
}

# ==========
# Plotting Functions
# ==========

# Plot arcs
create_arc_plot <- function(bezier_data, region_range, alpha, arc_colors = NULL) {
  message("Creating arc plot")
  
  arc_linewidth <- 200 / (region_range$end - region_range$start)
  
  p <- ggplot(bezier_data, aes(x = x, y = y, group = arc_group, color = bin)) +
    geom_bezier(linewidth = arc_linewidth, alpha = alpha) + 
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    scale_x_continuous(limits = c(region_range$start, region_range$end)) +
    theme_classic() +
    labs(x = "", y = "") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  # Apply custom colors if provided
  if (!is.null(arc_colors)) {
    arc_cols <- strsplit(arc_colors, ",")[[1]]
    p <- p + scale_color_manual(
      name = "Pairing probability",
      values = c("low" = arc_cols[1], "mid" = arc_cols[2], "high" = arc_cols[3]),
      labels = c("low" = "10–30%", "mid" = "30–80%", "high" = ">80%")
    )
  } else {
    p <- p + scale_color_manual(
      name = "Pairing probability",
      values = c("low" = "#d3d3d3", "mid" = "#95CDF0", "high" = "#426AB3"),
      labels = c("low" = "10–30%", "mid" = "30–80%", "high" = ">80%")
    )
  }
  
  return(p)
}

# Plot reactivity/Shannon data
create_measurement_plot <- function(data_list, region_range, data_type = "Reactivity", 
                                    do_smoothing = FALSE, color = NULL, do_median = FALSE,
                                    break_y = FALSE, limit_y = FALSE) {
  
  message(paste("Creating", data_type, "plot"))
  
  # Extract data from the list
  data_summary <- data_list$data_summary
  ribbons <- data_list$ribbons
  
  # Default colors
  default_color <- if(data_type == "Reactivity") "#426AB3" else "#95CDF0"
  plot_color <- ifelse(is.null(color), default_color, color)
  
  # Add a flag for rows where SD is greater than 0
  data_summary$show_error_bar <- data_summary$sd_value > 0  # True if SD is > 0, False if SD is 0
  
  if (do_smoothing) {
    p <- ggplot(data_summary) +
      geom_line(aes(x = start, y = smooth_value), color = plot_color) +
      geom_polygon(data = ribbons,
                   aes(x = x, y = y), 
                   fill = plot_color, alpha = 0.5) +
      labs(x = "", y = paste("Smoothed\n", data_type)) +
      theme_classic() +
      theme(panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5)) +
      theme(panel.grid.minor.x = element_line(color = "grey90", linewidth = 0.4)) +
      scale_y_continuous(labels = comma, expand = c(0, 0)) +
      scale_x_continuous(limits = c(region_range$start, region_range$end))
    
    if (do_median) {
      median_value <- median(data_summary$smooth_value, na.rm = TRUE)
      p <- p + geom_hline(yintercept = median_value, linetype = "dashed", color = "black")
    }
  } else {
    p <- ggplot(data_summary) +
      geom_bar(aes(x = start, y = mean_value), stat = "identity", fill = plot_color) +
      # Only show error bars where SD is greater than 0
      geom_errorbar(data = subset(data_summary, show_error_bar == TRUE),
                    aes(x = start, ymin = mean_value - sd_value, ymax = mean_value + sd_value)) +
      labs(x = "", y = data_type) +
      theme_classic() +
      theme(panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5)) +
      theme(panel.grid.minor.x = element_line(color = "grey90", linewidth = 0.4)) +
      scale_y_continuous(labels = comma, expand = c(0, 0)) +
      scale_x_continuous(limits = c(region_range$start, region_range$end))
    
    if (break_y) {
      message(paste("Splitting", data_type, "y axis at 2.5"))
      p <- p + scale_y_break(breaks=c(2.5), scales = c(1, 0.2), space=.2)
    }
    
    if (limit_y) {
      message(paste("Limiting", data_type, "y axis to 2.5"))
      p <- p + coord_cartesian(ylim = c(0, 2.5))
    }
  }
  
  return(p)
}

# Create pairwise comparison plot
create_comparison_plot <- function(comparison_data, region_range, data_type = "Reactivity") {
  message(paste("Creating pairwise", data_type, "comparison plot"))
  
  comparison_data <- prepare_plot_data(comparison_data, compare_log = opt$compare_log)
  
  label_prefix <- if (opt$compare_log) "Log " else ""
  
  p <- ggplot(data = comparison_data, aes(x = reactivity_1_transformed, y = reactivity_2_transformed)) +
    geom_pointdensity(size = 0.5, adjust = 0.2) +
    geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "solid", linewidth = 0.6) +
    scale_color_viridis(option = "viridis") +
    labs(
      x = paste(label_prefix, data_type),
      y = paste(label_prefix, data_type)
    ) +
    theme_classic(base_size = 12) +
    facet_wrap(~ pair, ncol = 3) +
    theme(
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.position = "none",
      aspect.ratio = 1
    )
  
  return(p)
}

# Save regression stats to a .txt file
save_regression_stats <- function(comparison_data) {
  message("Saving regression stats to a file")
  
  comparison_data <- prepare_plot_data(comparison_data, compare_log = opt$compare_log)
  regression_stats <- list()
  
  for (pair_name in unique(comparison_data$pair)) {
    pair_data <- subset(comparison_data, pair == pair_name)
    
    lm_model <- lm(reactivity_2_transformed ~ reactivity_1_transformed, data = pair_data)
    r_squared <- summary(lm_model)$r.squared
    intercept <- coef(lm_model)[1]
    slope <- coef(lm_model)[2]
    
    slope_deviation_test <- abs(slope - 1) / summary(lm_model)$coefficients[2, 2]
    slope_deviation_p_value <- 2 * pt(-abs(slope_deviation_test), df = lm_model$df.residual)
    
    regression_stats[[pair_name]] <- list(R_squared = r_squared, Intercept = intercept,
                                          Slope = slope, Slope_deviation_p_value = slope_deviation_p_value)
  }
  
  base_output <- ifelse(grepl("\\.pdf$", opt$output, ignore.case = TRUE),
                        sub("\\.pdf$", "", opt$output, ignore.case = TRUE),
                        opt$output)
  regression_file_name <- paste0(base_output, "_comparisons_stats.txt")
  
  write_lines <- c("Pairwise Linear Regression Stats:\n")
  for (pair_name in names(regression_stats)) {
    stats <- regression_stats[[pair_name]]
    write_lines <- c(write_lines, paste(pair_name, ": R-squared =", round(stats$R_squared, 4),
                                        ", Intercept =", round(stats$Intercept, 4),
                                        ", Slope =", round(stats$Slope, 4),
                                        ", Slope deviation p-value =", formatC(stats$Slope_deviation_p_value, format = "e", digits = 2)))
  }
  
  writeLines(write_lines, regression_file_name)
  message(paste("Regression stats saved to", regression_file_name))
}

# Plot QC data

# Calculate the highest whisker (Q3 + 1.5 * IQR) for each group and find the maximum across all groups
qc_highest_whisker <- function(qc_data) {
  # Get the unique groups (e.g., the "base" or "condition" factors)
  unique_groups <- unique(qc_data$base)
  
  # Initialize a vector to store the highest whisker for each group
  highest_whiskers <- numeric(length(unique_groups))
  
  # Loop through each group and calculate the highest whisker
  for (i in seq_along(unique_groups)) {
    # Subset the data for this group
    group_data <- subset(qc_data, base == unique_groups[i])
    
    # Calculate quantiles and IQR for this group
    q1 <- quantile(group_data$reactivity, 0.25, na.rm = TRUE)
    q3 <- quantile(group_data$reactivity, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    
    # Calculate the upper whisker limit as Q3 + 1.5 * IQR
    upper_whisker_limit <- q3 + 1.5 * iqr
    
    # Find the largest value in the group that is <= the upper whisker limit
    highest_whiskers[i] <- max(group_data$reactivity[group_data$reactivity <= upper_whisker_limit], na.rm = TRUE)
  }
  
  # Return the overall highest whisker value across all groups
  return(max(highest_whiskers, na.rm = TRUE))
}


# Helper function to get the next rounded number above a given value
rounded_y <- function(x) {
  # Find the nearest rounded number by rounding up
  if (x == 0) return(1)  # Special case for zero
  magnitude <- 10^floor(log10(x))  # Find the scale of the number (e.g., 1, 10, 100)
  fraction <- x / magnitude  # Scale the number to [1, 10)
  
  # Rounding up to the next rounded number
  if (fraction <= 1.5) {
    return(1.5 * magnitude)
  } else if (fraction <= 2.5) {
    return(2.5 * magnitude)
  } else if (fraction <= 5) {
    return(5 * magnitude)
  } else {
    return(10 * magnitude)
  }
}

create_qc_plot <- function(qc_data, qc_colors = NULL, qc_max_y = 0.15, qc_auto_y = FALSE) {
  message("Creating QC plot")
  
  p <- ggplot(qc_data, aes(x = base, y = reactivity, fill = condition)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.6, outlier.shape = NA) +
    scale_y_continuous(labels = comma) +
    theme_classic(base_size = 24) +
    labs(
      title = "",
      x = "Base",
      y = "Density",
      fill = "Condition"
    )
  
  # Apply custom colors if provided
  if (!is.null(qc_colors)) {
    qc_cols <- strsplit(qc_colors, ",")[[1]]
    p <- p + scale_fill_manual(values = c("control" = qc_cols[1], "treated" = qc_cols[2]))
  } else {
    p <- p + scale_fill_manual(values = c("control" = "#d3d3d3", "treated" = "#426AB3"))
  }
  
  # Calculate highest whisker if qc_auto_y is TRUE
  if (qc_auto_y) {
    highest_whisker <- qc_highest_whisker(qc_data)
    upper_limit <- rounded_y(highest_whisker)
  }
  
  # Decide on y-axis behavior
  if (qc_auto_y) {
    p <- p + coord_cartesian(ylim = c(0, upper_limit))  # Set the dynamic upper limit
  } else {
    p <- p + coord_cartesian(ylim = c(0, qc_max_y))  # Fixed limit
  }
  
  return(p)
}

# Combine plots
combine_plots <- function(plots_list, heights_list, base_output) {
  message("Combining plots...")
  
  # Suppress x-axis for all but the bottom plot
  if (length(plots_list) > 1) {
    for (i in 1:(length(plots_list) - 1)) {
      plots_list[[i]] <- plots_list[[i]] + theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
    }
  }
  
  # Add margins and legend starts
  for (i in 1:length(plots_list)) {
    legend_pos <- if(i == 1) "left" else "none"
    plots_list[[i]] <- plots_list[[i]] + 
      theme(plot.margin = margin(2, 2, 2, 2), legend.start = legend_pos)
  }
  
  # Combine plots
  combined_plot <- wrap_plots(plots_list) +
    plot_layout(ncol = 1, heights = heights_list)
  
  height_per_unit <- 1.2
  total_height <- sum(heights_list) * height_per_unit
  
  # Save combined plot
  ggsave(
    plot = combined_plot,
    filename = paste0(base_output, ".pdf"),
    width = 8,
    height = total_height
  )
  
  message("Plots combined and saved to", paste0(base_output, ".pdf"))
}

# ==========
# Main Execution
# ==========

main <- function() {
  # Validate inputs
  validate_inputs(opt)
  
  # Parse region
  region_range <- parse_region(opt$region)
  message("Using region range: ", paste0(region_range$start, ":", region_range$end))
  
  # Initialize plot and height lists
  plots_list <- list()
  heights_list <- c()
  
  # Process dotplot data and create arc plot if specified
  if (!is.null(opt$dotplot)) {
    bezier_data <- process_dotplot(opt$dotplot, region_range, opt$arc_height)
    p.arcs <- create_arc_plot(bezier_data, region_range, opt$alpha, opt$arc_colors)
    plots_list <- c(plots_list, list(p.arcs))
    heights_list <- c(heights_list, 2)
  }
  
  # Process reactivity data and create plot if specified
  if (!is.null(opt$reactivity)) {
    reactivity_summary <- process_measurement_data(
      opt$reactivity, region_range, 
      do_smoothing = opt$smoothing, 
      window_size = opt$window, 
      do_median_normalize = opt$median,
      data_type = "reactivity"
    )
    
    # Get color if specified
    reactivity_color <- NULL
    if (!is.null(opt$colors)) {
      group_cols <- strsplit(opt$colors, ",")[[1]]
      reactivity_color <- group_cols[1]
    }
    
    p.reactivity <- create_measurement_plot(
      reactivity_summary, region_range,
      data_type = "Reactivity",
      do_smoothing = opt$smoothing,
      color = reactivity_color,
      do_median = opt$median,
      break_y = opt$break_y,
      limit_y = opt$limit_y
    )
    
    plots_list <- c(plots_list, list(p.reactivity))
    heights_list <- c(heights_list, 1)
  }
  
  # Process Shannon entropy data and create plot if specified
  if (!is.null(opt$shannon)) {
    shannon_summary <- process_measurement_data(
      opt$shannon, region_range, 
      do_smoothing = opt$smoothing, 
      window_size = opt$window, 
      do_median_normalize = opt$median,
      data_type = "shannon"
    )
    
    # Get color if specified
    shannon_color <- NULL
    if (!is.null(opt$colors)) {
      group_cols <- strsplit(opt$colors, ",")[[1]]
      shannon_color <- if(length(group_cols) > 1) group_cols[2] else NULL
    }
    
    p.shannon <- create_measurement_plot(
      shannon_summary, region_range,
      data_type = "Shannon Entropy",
      do_smoothing = opt$smoothing,
      color = shannon_color,
      do_median = opt$median
    )
    
    plots_list <- c(plots_list, list(p.shannon))
    heights_list <- c(heights_list, 1)
  }
  
  # Create comparison plot and save stats to separate file if requested
  if (!is.null(opt$compare_replicates) && opt$compare_replicates) {
    comparison_data <- handle_replicate_comparisons(
      file_paths = opt$reactivity,
      region_range = region_range
    )
    
    # Generate comparison plot
    comparison_plot <- create_comparison_plot(comparison_data, region_range, data_type = "Reactivity")
    
    # Save comparison plot separately
    base_output <- ifelse(grepl("\\.pdf$", opt$output, ignore.case = TRUE),
                          sub("\\.pdf$", "", opt$output, ignore.case = TRUE),
                          opt$output)
    
    ggsave(
      filename = paste0(base_output, "_comparisons.pdf"),
      plot = comparison_plot,
      width = 8,
      height = 6
    )
    
    message("Comparison plot saved to", paste0(base_output, "_comparisons.pdf"))
    
    # Save regression stats to a separate file
    save_regression_stats(comparison_data)
  }
  
  # Process QC data and create plot if specified
  if (!is.null(opt$qc) && opt$qc) {
    qc_data <- process_qc_data(opt$rc_files, opt$rc_controls)
    p.qc <- create_qc_plot(qc_data, qc_colors = opt$qc_colors, 
                           qc_max_y = opt$qc_max_y, qc_auto_y = opt$qc_auto_y)
    
    # Save QC plot separately
    base_output <- ifelse(grepl("\\.pdf$", opt$output, ignore.case = TRUE),
                          sub("\\.pdf$", "", opt$output, ignore.case = TRUE),
                          opt$output)
    
    ggsave(
      filename = paste0(base_output, "_qc.pdf"),
      plot = p.qc,
      width = 8,
      height = 6
    )
    
    message("QC plot saved to", paste0(base_output, "_qc.pdf"))
  }
  
  # Combine and save main plots
  if (length(plots_list) > 0) {
    base_output <- ifelse(grepl("\\.pdf$", opt$output, ignore.case = TRUE),
                          sub("\\.pdf$", "", opt$output, ignore.case = TRUE),
                          opt$output)
    
    combine_plots(plots_list, heights_list, base_output)
    message("All plots completed successfully")
  } else {
    message("No main plots were generated. Exiting.")
  }
}

# Run the main function
tryCatch({
  main()
}, error = function(e) {
  message("ERROR: ", e$message)
  quit(status = 1)
})