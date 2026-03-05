# Analysis of the results of hidden ORF search

# 1. Environment setup
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("scales")) install.packages("scales")
library(tidyverse)
library(scales)
library(readr)

# 2. Data Loading and Preprocessing
output_dir <- "D://data//Mao Lab//hidden ORF search//code//analysis//min_aa 200"
if (!dir.exists(output_dir)) dir.create(output_dir)
print("Starting analysis workflow...")
df <- read_csv("D://data//Mao Lab//hidden ORF search//code//analysis//min_aa 200//results.csv", col_types = cols(strand = col_character()))
df <- df %>%
  mutate(risk_category = case_when(
    expressible == TRUE ~ "High (Expressible)",
    has_promoter == TRUE & expressible == FALSE ~ "Medium (Promoter Only)",
    overlap != "outside_target_gene" ~ "Target Overlap",
    TRUE ~ "Low Risk"
  ))

all_plasmids <- df %>%
  select(plasmid_id) %>%
  distinct()

# 3. Data Analysis and Visualization

# 3.1 Expressible ORF Count per Plasmid Distribution
df_high_risk <- df %>%
  filter(expressible == TRUE)

orf_counts <- df_high_risk %>%
  group_by(plasmid_id) %>%
  summarise(total_orfs = n(), .groups = 'drop') %>%
  full_join(all_plasmids, by = "plasmid_id") %>%
  mutate(total_orfs = replace_na(total_orfs, 0))

unique_vals <- sort(unique(orf_counts$total_orfs))

custom_levels <- c(unique_vals[unique_vals != 0], 0)
orf_counts <- orf_counts %>%
  mutate(total_orfs_factor = factor(total_orfs, levels = custom_levels))

p1 <- ggplot(orf_counts, aes(x = total_orfs_factor)) +
  geom_bar(fill = "#4285F4", color = "white", width = 0.6) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4, color = "black") +
  scale_x_discrete(labels = function(x) ifelse(x == "0", "Not Applicable", x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  labs(title = "Expressible ORF Count per Plasmid Distribution",
       x = "Number of Expressible ORFs", 
       y = "Frequency (Number of Plasmids)")

ggsave(file.path(output_dir, "1_ORF_count_distribution.png"), p1, width = 8, height = 6)
write_csv(orf_counts, file.path(output_dir, "1_ORF_counts_per_plasmid.csv"))

# 3.2 Distribution of Frames for hidden ORFs
df_actual_risk$TGframe <- as.character(df_actual_risk$TGframe)
idx_to_fix <- !grepl("^[+-]", df_actual_risk$TGframe)
df_actual_risk$TGframe[idx_to_fix] <- paste0("+", df_actual_risk$TGframe[idx_to_fix])
all_possible_frames <- c("+1", "+2", "+3", "-1", "-2", "-3")

p3_risk <- ggplot(df_actual_risk, aes(x = factor(TGframe, levels = all_possible_frames))) +
  geom_bar(fill = "#4285F4", color = "white", width = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  labs(
    title = "Distribution of Frames for hidden ORFs",
    x = "Frame",
    y = "Total Count"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0)
  )

ggsave(file.path(output_dir, "2_hidden_ORFs_frame_bias.png"), p3_risk, width = 8, height = 6)

# 3.3 Hidden_ORF_Distribution within Expressible ORFs
risk_comparison_stats <- df %>%
  filter(expressible == TRUE) %>%
  mutate(risk_status = ifelse(is_risk == TRUE, "Hidden ORFs (Frameshift/Antisense)", "Expected ORFs (Target Gene)")) %>%
  group_by(risk_status) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(percentage = count / sum(count) * 100)

p4 <- ggplot(risk_comparison_stats, aes(x = "", y = count, fill = risk_status)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Hidden ORFs (Frameshift/Antisense)" = "#EA4335", 
                               "Expected ORFs (Target Gene)" = "#4285F4")) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "hidden_ORF_Distribution within Expressible ORFs") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 5)

ggsave(file.path(output_dir, "3_hidden_orf_distribution_pie.png"), p4, width = 8, height = 8)
df_export <- df[which(df$is_risk == TRUE), ]
write_csv(df_export, file.path(output_dir, "3_hidden_orf_summary.csv"))

