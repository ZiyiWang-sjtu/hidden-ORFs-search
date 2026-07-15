# Analysis of the results of hidden ORF search

# 1. Environment setup
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("scales")) install.packages("scales")
library(tidyverse)
library(scales)
library(readr)

# 2. Data loading and preprocessing
output_dir_1 <- "YOUR_PATH"
output_dir_2 <- "YOUR_PATH"

if (!dir.exists(output_dir_1)) dir.create(output_dir_1, recursive = TRUE)
if (!dir.exists(output_dir_2)) dir.create(output_dir_2, recursive = TRUE)

print("Starting analysis workflow...")

df <- read_csv(
  "YOUR_CSV", 
  col_types = cols(
    strand = col_character(),
    is_hidden_ORFs = col_logical(),
    is_putatively_expressible_hidden_ORFs = col_logical()
  )
)

df <- df %>%
  mutate(TGframe = case_when(
    TGframe == "+1" ~ "CF",
    TGframe == "+2" ~ "CF+1",
    TGframe == "+3" ~ "CF+2",
    TGframe == "-1" ~ "RF",
    TGframe == "-2" ~ "RF+1",
    TGframe == "-3" ~ "RF+2",
    TRUE ~ TGframe
  ))

all_plasmids <- df %>%
  select(plasmid_id) %>%
  distinct()

# 3. Set up frame categories
all_possible_frames <- c("CF", "CF+1", "CF+2", "RF", "RF+1", "RF+2")

# 4. Core analysis and visualization pipeline
run_analysis_pipeline <- function(data_subset, category_name, target_dir) {
  
  cat("\nRunning analysis for:", category_name, "...\n")
  cat("Target directory is:", target_dir, "\n")
  
  # 4.1 ORF count per plasmid distribution

  not_applicable <- df %>%
    group_by(plasmid_id) %>%
    summarise(
      not_applicable = all(!has_promoter & !has_terminator),
      .groups = "drop"
    )
  
  hidden_counts <- data_subset %>%
    filter(TGframe != "CF") %>%
    group_by(plasmid_id) %>%
    summarise(
      hidden_count = n(),
      .groups = "drop"
    )
  
  orf_counts <- all_plasmids %>%
    left_join(not_applicable, by = "plasmid_id") %>%
    left_join(hidden_counts, by = "plasmid_id") %>%
    mutate(
      hidden_count = replace_na(hidden_count, 0),
      total_orfs = ifelse(
        not_applicable,
        NA,                 
        hidden_count + 1    
      )
    )
  
  levels_order <- c(
    sort(unique(orf_counts$total_orfs[!is.na(orf_counts$total_orfs)])),
    NA
  )
  
  orf_counts <- orf_counts %>%
    mutate(
      total_orfs_factor = factor(
        ifelse(is.na(total_orfs), "Not Applicable", as.character(total_orfs)),
        levels = c(
          as.character(sort(unique(total_orfs[!is.na(total_orfs)]))),
          "Not Applicable"
        )
      )
    )
  
  p1 <- ggplot(orf_counts, aes(x = total_orfs_factor)) +
    geom_bar(fill = "#4285F4", color = "white", width = 0.6) +
    geom_text(
      stat = "count",
      aes(label = after_stat(count)),
      vjust = -0.5,
      size = 4
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title ="ORF Count per Plasmid Distribution",
      x = "Number of ORFs",
      y = "Frequency (Number of Plasmids)"
    )
  ggsave(file.path(target_dir, "1_ORF_count_distribution.png"), p1, width = 8, height = 6, create.dir = TRUE)
  write_csv(orf_counts, file.path(target_dir, "1_ORF_counts_per_plasmid.csv"))
  
  # 4.2 Distribution of frames 
  df_actual_risk <- data_subset %>% 
    filter(TGframe != "CF")
  
  df_frame_plot <- df_actual_risk %>%
    filter(TGframe %in% all_possible_frames) %>%
    mutate(TGframe = factor(TGframe, levels = all_possible_frames))
  
  p2 <- ggplot(df_frame_plot, aes(x = factor(TGframe, levels = all_possible_frames))) +
    geom_bar(fill = "#4285F4", color = "white", width = 0.7) +
    geom_text(stat = "count",
              aes(label = after_stat(count)),
              vjust = -0.5,
              size = 4) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0)
    ) +
    labs(
      title = paste("Distribution of Frames for", category_name),
      x = "Frame",
      y = "Total Count"
    )
  
  ggsave(file.path(target_dir, "2_frame_bias.png"), p2, width = 8, height = 6, create.dir = TRUE)
  
  # 4.3 Distribution Pie Chart
  risk_comparison_stats <- data_subset %>%
    mutate(risk_status = ifelse(TGframe != "CF", "Hidden ORFs (Frameshift/Antisense)", "Expected ORFs (Target Gene)")) %>%
    group_by(risk_status) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(percentage = count / sum(count) * 100)
  
  p3 <- ggplot(risk_comparison_stats, aes(x = "", y = count, fill = risk_status)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("Hidden ORFs (Frameshift/Antisense)" = "#EA4335", 
                                 "Expected ORFs (Target Gene)" = "#4285F4")) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = paste("Distribution within", category_name)) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 5)
  
  ggsave(file.path(target_dir, "3_distribution_pie.png"), p3, width = 8, height = 8, create.dir = TRUE)
  
  df_export <- data_subset %>% filter(TGframe != "CF")
  write_csv(df_export, file.path(target_dir, "3_summary.csv"))
}


# 5. Run the pipeline for different ORF categories
# Analyze all hidden ORFs
df_hidden_orfs <- df %>% filter(is_hidden_ORFs == TRUE , overlap == "spanning_target_gene")
run_analysis_pipeline(
  data_subset = df_hidden_orfs, 
  category_name = "Hidden ORFs", 
  target_dir = output_dir_1
)

# Analyze all putatively expressible hidden ORFs
df_expressible_orfs <- df %>% filter(is_putatively_expressible_hidden_ORFs == TRUE ,  overlap == "spanning_target_gene")
run_analysis_pipeline(
  data_subset = df_expressible_orfs, 
  category_name = "Putatively Expressible Hidden ORFs", 
  target_dir = output_dir_2
)

print("All tasks completed successfully!")
