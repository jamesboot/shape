# Script for coverage plots of tRNA

# Setwd
setwd('/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v3')

# Packages
library(data.table)
library(ggplot2)
library(ggpubr)

# Locate all coverage files
AlaCovFiles <- list.files(path = '03_rf-count_Ala_outs',
                          pattern = '*_cov.txt',
                          recursive = T,
                          full.names = T) 

ProCovFiles <- list.files(path = '03_rf-count_Pro_outs',
                          pattern = '*_cov.txt',
                          recursive = T,
                          full.names = T)

# Define output folders
AlaOut <- '03_rf-count_Ala_outs'
ProOut <- '03_rf-count_Pro_outs'

# Loop through Ala samples
for (x in AlaCovFiles) {
  
  # Load data
  t <- read.delim(
    x,
    col.names = c('Base', 'n_RT_Stops', 'Coverage'),
    header = F,
    skip = 1
  )
  
  # Add position column
  t$Pos <- as.numeric(rownames(t))
  
  # Define sample
  sample <- gsub('.txt', '', basename(x))
  
  # Plot coverage
  p1 <- ggplot(t, aes(x = Pos, y = Coverage, group = 1)) +
    geom_line(color = 'forestgreen') +
    theme_bw(base_size = 8) +
    xlab('tRNA Position') +
    scale_x_continuous(breaks = seq(0, 75, by = 5))
  
  # Plot mutations
  p2 <- ggplot(t, aes(x = Pos, y = n_RT_Stops, group = 1)) +
    geom_line(color = 'darkblue') +
    theme_bw(base_size = 8) +
    xlab('tRNA Position') +
    ylab('n RT Stops') +
    scale_x_continuous(breaks = seq(0, 75, by = 5))
  
  p3 <- ggarrange(p1, p2, ncol = 1, nrow = 2)
  
  ggsave(
    paste0(AlaOut, '/', sample, '.png'),
    plot = p3,
    height = 3,
    width = 9
  )
}

# Loop through Pro samples
for (x in ProCovFiles) {
  
  # Load data
  t <- read.delim(
    x,
    col.names = c('Base', 'n_RT_Stops', 'Coverage'),
    header = F,
    skip = 1
  )
  
  # Add position column
  t$Pos <- as.numeric(rownames(t))
  
  # Define sample
  sample <- gsub('.txt', '', basename(x))
  
  # Plot coverage
  p1 <- ggplot(t, aes(x = Pos, y = Coverage, group = 1)) +
    geom_line(color = 'forestgreen') +
    theme_bw(base_size = 8) +
    xlab('tRNA Position') +
    scale_x_continuous(breaks = seq(0, 75, by = 5))
  
  # Plot mutations
  p2 <- ggplot(t, aes(x = Pos, y = n_RT_Stops, group = 1)) +
    geom_line(color = 'darkblue') +
    theme_bw(base_size = 8) +
    xlab('tRNA Position') +
    ylab('n RT Stops') +
    scale_x_continuous(breaks = seq(0, 75, by = 5))
  
  p3 <- ggarrange(p1, p2, ncol = 1, nrow = 2)
  
  ggsave(
    paste0(ProOut, '/', sample, '.png'),
    plot = p3,
    height = 3,
    width = 9
  )
}


