# Script for plotting delta shape bar plots
# I.e. plot differences in reactivity between conditions 

library(dplyr)
library(ggplot2)
library(XML)

setwd('../trna_shape_v3/')

# Plot difference between Proline and Pro-Dicer within treatment

# Load meta
meta <- read.csv('meta.csv')

# Make a list of samples for each condition
sampleLists <- list(
  Pro_DMS = grep('DMS_', list.dirs('04_rf-norm_Pro_outs'), value = T),
  Pro_5NIA = grep('5NIA_', list.dirs('04_rf-norm_Pro_outs'), value = T),
  ProDic_DMS = grep('DMS_', list.dirs('04_rf-norm_Pro_Dic_outs'), value = T),
  ProDic_5NIA = grep('5NIA_', list.dirs('04_rf-norm_Pro_Dic_outs'), value = T)
) 

# For each list element load the xml files, average
reactivityLists <- list()
for (x in names(sampleLists)) {
  for (y in c(1:3)) {
    # Parse xml file
    input <- XML::xmlParse(file = paste0(sampleLists[[x]][y], '/tRNA.xml'))
    # Extract the reactivity node text
    reactivity_text <- xpathSApply(input, "//reactivity", xmlValue)
    # Clean up: remove whitespace, split by comma
    reactivity_vals <- unlist(strsplit(gsub("\\s+", "", reactivity_text), ","))
    # Drop empty strings (from trailing commas) and convert to numeric
    reactivity_numeric <- as.numeric(reactivity_vals[reactivity_vals != ""])
    # Put into a data.frame
    df <- data.frame(reactivity = reactivity_numeric)
    # Add to list
    reactivityLists[[x]][[y]] <- df
  }
}

# cbind the lists into single dataframe
reactivities <- lapply(reactivityLists, function(x){
  merged <- cbind(x[[1]], x[[2]], x[[3]])
  colnames(merged) <- c('r1', 'r2', 'r3')
  merged <- merged %>%
    rowwise() %>%
    mutate(
      row_mean = mean(c_across(everything()), na.rm = F),
      row_sd   = sd(c_across(everything()), na.rm = F),
      mean_p_sd = row_mean + row_sd,
      mean_m_sd = row_mean - row_sd
    ) %>%
    ungroup() %>%
    mutate(position = row_number()) %>%
    select(position, row_mean, mean_p_sd, mean_m_sd)
  return(merged)
})

# Now calculate delta and plot - for Pro vs Dicer DMS
delta_mean <- reactivities$Pro_DMS$row_mean - reactivities$ProDic_DMS$row_mean
delta_p_sd <- reactivities$Pro_DMS$mean_p_sd - reactivities$ProDic_DMS$mean_p_sd
delta_m_sd <- reactivities$Pro_DMS$mean_m_sd - reactivities$ProDic_DMS$mean_m_sd

ggplot(data.frame(delta_mean), aes(seq_along(delta_mean), delta_mean)) + 
  geom_bar(stat = "identity", fill = 'skyblue') +
  geom_errorbar(aes(ymin = delta_m_sd, ymax = delta_p_sd),
                width = 0.2, colour = "black") +
  theme_bw() +
  xlab('tRNA Position') +
  ylab('Delta Norm SHAPE Reactivity') +
  ggtitle('Pro vs Dicer DMS')

# Now calculate delta and plot - for Pro vs Dicer DMS
delta_mean <- reactivities$Pro_5NIA$row_mean - reactivities$ProDic_5NIA$row_mean
delta_p_sd <- reactivities$Pro_5NIA$mean_p_sd - reactivities$ProDic_5NIA$mean_p_sd
delta_m_sd <- reactivities$Pro_5NIA$mean_m_sd - reactivities$ProDic_5NIA$mean_m_sd

ggplot(data.frame(delta_mean), aes(seq_along(delta_mean), delta_mean)) + 
  geom_bar(stat = "identity", fill = 'skyblue') +
  geom_errorbar(aes(ymin = delta_m_sd, ymax = delta_p_sd),
                width = 0.2, colour = "black") +
  theme_bw() +
  xlab('tRNA Position') +
  ylab('Delta Norm SHAPE Reactivity') +
  ggtitle('Pro vs Dicer DMS')
