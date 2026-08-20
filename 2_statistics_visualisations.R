# =============================================================================
# Content: calibration plot, summary stats, Wilcoxon sign tests, 
#and all figures (boxplots, Cleveland plots, ratio plots) built from the outputs of 1_processing.R.

# Steps in this workflow:
#1. Calibration plots
#2. Summary statistics
#3. Wilcoxon signed-rank test
#4. Boxplots
#5. Cleveland plots
#6. Ratio analysis (HumHLr/BaseHLr)
# =============================================================================

library(here)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(writexl)
library(ggpubr)
library(ggrepel)
library(patchwork)

output_data_dir <- here("output_data")
results_dir <- file.path(output_data_dir, "results")

#small buffer for 0 values on log plots
logbuffer <- 1e-6  
#selected threshold
selected_thr <- 0.20

hl_data <- read_xlsx(file.path(output_data_dir, "habitat_losses(R).xlsx"))
calibration_stats <- read_xlsx(file.path(output_data_dir, "calibration_stats.xlsx"))

#----1.Calibration plots---- 
cal_long <- calibration_stats %>%
  select(threshold, MAE, MedAE, RMSE) %>%
  pivot_longer(cols = c(MAE, MedAE, RMSE), names_to = "metric", values_to = "value")

err_colors <- c(
  "MAE"   = "#2166ac",
  "MedAE" = "#b2182b",
  "RMSE"  = "#4d4d4d")
x_breaks <- sort(unique(calibration_stats$threshold))

shared_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.ticks = element_blank(),
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", colour = NA))

#Panel A: MAE, MedAE & RMSE across thresholds
panel_A <- ggplot(cal_long, aes(x = threshold, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = selected_thr, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  scale_color_manual(values = err_colors) +
  scale_x_continuous(breaks = x_breaks) +
  shared_theme +
  labs(x = "HMc threshold", y = expression("Error (km"^2*")"), color = NULL) +
  theme(legend.position = "top")

#Panel B: Bland-Altman at the selected threshold ----
ba_data <- hl_data %>%
  mutate(
    mean_val = (global_humHL + local_humHL) / 2,
    diff_val = global_humHL - local_humHL)

mean_diff <- mean(ba_data$diff_val, na.rm = TRUE)
sd_diff <- sd(ba_data$diff_val,   na.rm = TRUE)
upper <- mean_diff + 1.96 * sd_diff
lower <- mean_diff - 1.96 * sd_diff

top_outliers <- ba_data %>% slice_max(abs(diff_val), n = 2)

panel_B <- ggplot(ba_data, aes(x = mean_val, y = diff_val)) +
  geom_point(size = 2, color = "#4292c6", alpha = 0.7) +
  geom_hline(yintercept = mean_diff, linetype = "solid",  color = "black",  linewidth = 0.6) +
  geom_hline(yintercept = upper,     linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = lower,     linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data  = top_outliers,
    aes(label = island),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey60",
    segment.size  = 0.3) +
  coord_cartesian(clip = "off") +
  shared_theme +
  labs(
    x = expression("Mean HumHL (km"^2*")"),
    y = bquote("Global (" * .(selected_thr) * ") " ~ "−" ~ " Local HumHL (km"^2*")"))

calibration_plot <- panel_A / panel_B + plot_annotation(tag_levels = "A")
print(calibration_plot)

ggsave(file.path(results_dir, "calibration_plot.png"),
       plot = calibration_plot, dpi = 300, height = 12, width = 10)

#----2. Summary stats---- 
hl_summ <- hl_data %>%
  transmute(
    island, archipelago,
    seaHLr, volHLr, local_humHLr, global_humHLr,
    baseHLr) %>%
  mutate(
    diff_local_base   = local_humHLr  - baseHLr,
    diff_global_base  = global_humHLr - baseHLr,
    diff_global_local = global_humHLr - local_humHLr,
    diff_sea_vol      = seaHLr - volHLr)

overall_summ <- hl_summ %>%
  summarise(
    mean_baseHLr = mean(baseHLr, na.rm = TRUE),
    mean_local_humHLr = mean(local_humHLr, na.rm = TRUE),
    mean_global_humHLr = mean(global_humHLr, na.rm = TRUE),
    mean_seaHLr = mean(seaHLr, na.rm = TRUE),
    mean_volHLr = mean(volHLr, na.rm = TRUE),
    
    median_baseHLr = median(baseHLr, na.rm = TRUE),
    median_local_humHLr = median(local_humHLr, na.rm = TRUE),
    median_global_humHLr = median(global_humHLr, na.rm = TRUE),
    median_seaHLr = median(seaHLr, na.rm = TRUE),
    median_volHLr = median(volHLr, na.rm = TRUE),
    
    mean_diff_local_base = mean(diff_local_base, na.rm = TRUE),
    mean_diff_global_base = mean(diff_global_base, na.rm = TRUE),
    mean_diff_global_local = mean(diff_global_local, na.rm = TRUE),
    mean_diff_sea_vol = mean(diff_sea_vol, na.rm = TRUE)
  )

archipelago_summ <- hl_summ %>%
  group_by(archipelago) %>%
  summarise(
    n_islands = n(),
    
    mean_baseHLr = mean(baseHLr, na.rm = TRUE),
    mean_local_humHLr = mean(local_humHLr, na.rm = TRUE),
    mean_global_humHLr = mean(global_humHLr, na.rm = TRUE),
    mean_seaHLr = mean(seaHLr, na.rm = TRUE),
    mean_volHLr = mean(volHLr, na.rm = TRUE),
    
    median_baseHLr = median(baseHLr, na.rm = TRUE),
    median_local_humHLr = median(local_humHLr, na.rm = TRUE),
    median_global_humHLr = median(global_humHLr, na.rm = TRUE),
    median_seaHLr = median(seaHLr, na.rm = TRUE),
    median_volHLr = median(volHLr, na.rm = TRUE),
    
    mean_diff_local_base = mean(diff_local_base, na.rm = TRUE),
    mean_diff_global_base = mean(diff_global_base, na.rm = TRUE),
    mean_diff_global_local = mean(diff_global_local, na.rm = TRUE),
    mean_diff_sea_vol = mean(diff_sea_vol, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(mean_diff_local_base))

island_summ <- hl_summ %>% arrange(archipelago, island)

absdiff_local_base <- hl_summ %>%
  mutate(abs_diff = abs(diff_local_base)) %>%
  arrange(abs_diff) %>%
  select(island, archipelago, baseHLr, local_humHLr, diff_local_base, abs_diff)

absdiff_global_base <- hl_summ %>%
  mutate(abs_diff = abs(diff_global_base)) %>%
  arrange(abs_diff) %>%
  select(island, archipelago, baseHLr, global_humHLr, diff_global_base, abs_diff)

absdiff_global_local <- hl_summ %>%
  mutate(abs_diff = abs(diff_global_local)) %>%
  arrange(desc(abs_diff)) %>%
  select(island, archipelago, global_humHLr, local_humHLr, diff_global_local, abs_diff)

absdiff_sea_vol <- hl_summ %>%
  mutate(abs_diff = abs(diff_sea_vol)) %>%
  arrange(desc(abs_diff)) %>%
  select(island, archipelago, seaHLr, volHLr, baseHLr, diff_sea_vol, abs_diff)

write_xlsx(
  list(
    overall = overall_summ,
    archipelago_summary = archipelago_summ,
    island_level = island_summ,
    absdiff_local_base = absdiff_local_base,
    absdiff_global_base = absdiff_global_base,
    absdiff_global_local = absdiff_global_local,
    absdiff_sea_vol = absdiff_sea_vol),
  file.path(results_dir, "hl_summary_stats.xlsx"))

#----3. Wilcoxon signed-rank tests ---- 
p_to_stars <- function(p) {
  if (is.na(p)) return("n/a")
  if (p < 0.0001) "****"
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else "ns"
}

metrics <- c("global_humHLr", "local_humHLr", "seaHLr", "volHLr")
metric_labels <- c(
  global_humHLr = "Global proxy",
  local_humHLr  = "Local proxy",
  seaHLr        = "SeaHLr",
  volHLr        = "VolHLr")

pairs <- combn(metrics, 2, simplify = FALSE)

wilcox_results <- lapply(pairs, function(p) {
  x <- hl_data[[p[1]]]
  y <- hl_data[[p[2]]]
  
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  
  res <- wilcox.test(x, y, paired = TRUE, exact = FALSE)
  
  data.frame(
    Comparison  = paste(metric_labels[p[1]], "vs", metric_labels[p[2]]),
    V_statistic = res$statistic,
    N_used      = sum(valid),
    P_value     = res$p.value)
}) %>%
  bind_rows()

wilcox_results$P_value_holm <- p.adjust(wilcox_results$P_value, method = "holm")

write_xlsx(wilcox_results, file.path(results_dir, "wilcox_test.xlsx"))

ann <- wilcox_results %>%
  separate(Comparison, into = c("group1", "group2"), sep = " vs ") %>%
  mutate(p = P_value_holm, label = vapply(p, p_to_stars, character(1))) %>%
  select(group1, group2, p, label)

#----4. Boxplots----
hl_long <- hl_data %>%
  select(island, archipelago, global_humHLr, local_humHLr, seaHLr, volHLr) %>%
  pivot_longer(
    cols = c(global_humHLr, local_humHLr, seaHLr, volHLr),
    names_to  = "metric",
    values_to = "rate") %>%
  mutate(
    category = if_else(metric %in% c("seaHLr", "volHLr"), "BaseHLr", "HumHLr"),
    metric_label = recode(metric,
                          seaHLr = "SeaHLr",
                          volHLr = "VolHLr",
                          local_humHLr= "Local proxy",
                          global_humHLr = "Global proxy"),
    metric_label = factor(metric_label,
                          levels = c("SeaHLr", "VolHLr", "Local proxy", "Global proxy")),
    category = factor(category, levels = c("BaseHLr", "HumHLr")),
    rate_logbuff = if_else(metric == "volHLr" & rate == 0, logbuffer, rate))

hl_cols <- c(
  "SeaHLr"                    = "#24bf70",
  "VolHLr"                    = "#167343",
  "BaseHLr (SeaHLr + VolHLr)" = "#2596be",
  "Local proxy"               = "#C2410C",
  "Local HumHLr"              = "#C2410C",
  "Global proxy"              = "#FDBA74",
  "Global HumHLr"             = "#FDBA74")

#linear
ymax_lin <- max(hl_long$rate, na.rm = TRUE)
ann_lin <- ann
ann_lin$y.position <- ymax_lin * 1.1 + (seq_len(nrow(ann_lin)) - 1) * (ymax_lin * 0.06)

hl_boxplot <- ggplot(hl_long, aes(x = metric_label, y = rate, fill = metric_label)) +
  geom_boxplot(width = 0.6, linewidth = 0.35, alpha = 0.65) +
  geom_jitter(width = 0.08, height = 0, size = 1.8, color = "grey30", alpha = 0.5) +
  scale_fill_manual(values = hl_cols) +
  theme_minimal(base_size = 12) +
  labs(
    x = "BaseHLr                                                    HumHLr",
    y = expression("Habitat loss rate (km"^2 * " year"^-1 * ")"),
    subtitle = "Wilcoxon signed-rank test, p.adjust = Holm") +
  theme(
    legend.position  = "none",
    axis.title.x  = element_text(margin = margin(t = 5)),
    axis.title.y  = element_text(margin = margin(r = 5)),
    axis.text.x = element_text(margin = margin(t = 3)),
    axis.text.y = element_text(margin = margin(r = 3)),
    axis.ticks  = element_blank(),
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = NA)) +
  stat_pvalue_manual(
    ann_lin, label = "label", xmin = "group1", xmax = "group2",
    y.position = "y.position", tip.length = 0.01, bracket.size = 0.35,
    size = 4, inherit.aes = FALSE)

print(hl_boxplot)
ggsave(file.path(results_dir, "LINboxplot_wilcoxon.png"),
       plot = hl_boxplot, width = 7, height = 7, dpi = 300)

#log
ymax_log <- max(hl_long$rate_logbuff, na.rm = TRUE)
ann_log <- ann
ann_log$y.position <- ymax_log * 0.12 + (seq_len(nrow(ann_log)) - 1) * (ymax_log * 0.025)

hl_LOGboxplot <- ggplot(hl_long, aes(x = metric_label, y = rate_logbuff, fill = metric_label)) +
  geom_boxplot(width = 0.6, linewidth = 0.35, alpha = 0.65) +
  geom_jitter(width = 0.08, height = 0, size = 1.8, color = "grey30", alpha = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = hl_cols) +
  theme_minimal(base_size = 12) +
  labs(
    x = "BaseHLr                                                    HumHLr",
    y = expression("Habitat loss rate (km"^2 * " year"^-1 * ", log"[10] * " scale)"),
    subtitle = "Wilcoxon signed-rank test, p.adjust = Holm") +
  theme(
    legend.position  = "none",
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.text.x = element_text(margin = margin(t = 3)),
    axis.text.y = element_text(margin = margin(r = 3)),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", colour = NA)) +
  stat_pvalue_manual(
    ann_log, label = "label", xmin = "group1", xmax = "group2",
    y.position = "y.position", tip.length = 0.01, bracket.size = 0.35,
    size = 4, inherit.aes = FALSE)

print(hl_LOGboxplot)
ggsave(file.path(results_dir, "LOGboxplot_wilcoxon.png"),
       plot = hl_LOGboxplot, width = 7, height = 7, dpi = 300)

#----5. Cleveland plots ----
hl_clv <- hl_data %>%
  transmute(island, archipelago, seaHLr, volHLr, baseHLr, local_humHLr, global_humHLr) %>%
  pivot_longer(
    cols = c(seaHLr, volHLr, baseHLr, local_humHLr, global_humHLr),
    names_to = "metric", values_to = "rate") %>%
  mutate(
    metric_label = recode(
      metric,
      seaHLr  = "SeaHLr",
      volHLr  = "VolHLr",
      baseHLr = "BaseHLr (SeaHLr + VolHLr)",
      local_humHLr= "Local HumHLr",
      global_humHLr = "Global HumHLr"),
    metric_label = factor(
      metric_label,
      levels = c("SeaHLr", "VolHLr", "BaseHLr (SeaHLr + VolHLr)", "Local HumHLr", "Global HumHLr")),
    category = case_when(
      metric %in% c("seaHLr", "volHLr", "baseHLr") ~ "BaseHLr",
      TRUE ~ "HumHLr"),
    category = factor(category, levels = c("BaseHLr", "HumHLr")),
    archipelago = factor(archipelago, levels = sort(unique(archipelago))),
    rate_plot = if_else(rate == 0, logbuffer, rate)) %>%
  arrange(archipelago, island) %>%
  group_by(archipelago) %>%
  mutate(island_f = factor(island, levels = rev(unique(island)))) %>%
  ungroup()

hl_shapes <- c(
  "SeaHLr" = 16,
  "VolHLr" = 16,
  "BaseHLr (SeaHLr + VolHLr)" = 16,
  "Local HumHLr" = 17,
  "Global HumHLr"= 17)

cleveland_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position    = "top",
    legend.title       = element_blank(),
    legend.box         = "horizontal",
    panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "#E0E0E0", linewidth = 0.35),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(angle = 0),
    axis.ticks         = element_blank(),
    strip.placement    = "outside",
    strip.background   = element_blank(),
    strip.text.y.right  = element_text(angle = 0, face = "bold", size = 9),
    panel.border       = element_blank(),
    plot.background    = element_rect(fill = "white"),
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.margin        = margin(8, 18, 8, 10))

#linear
hl_cleveland_linear <- ggplot(
  hl_clv, aes(x = rate_plot, y = island_f, color = metric_label, shape = metric_label)) +
  geom_point(size = 2.7, alpha = 0.8, stroke = 0) +
  facet_grid(archipelago ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = hl_cols) +
  scale_shape_manual(values = hl_shapes) +
  labs(x = expression("Habitat loss rate (km"^2*" year"^-1*")"), y = NULL) +
  cleveland_theme

print(hl_cleveland_linear)
ggsave(file.path(results_dir, "LINclv.png"),
       plot = hl_cleveland_linear, width = 10, height = 10, dpi = 300)

#log
hl_cleveland_log <- ggplot(
  hl_clv, aes(x = rate_plot, y = island_f, color = metric_label, shape = metric_label)) +
  geom_point(size = 2.7, alpha = 0.8, stroke = 0) +
  scale_x_log10() +
  facet_grid(archipelago ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = hl_cols) +
  scale_shape_manual(values = hl_shapes) +
  labs(x = expression("Habitat loss rate (km"^2*" year"^-1*", log"[10]*" scale)"), y = NULL) +
  cleveland_theme

print(hl_cleveland_log)
ggsave(file.path(results_dir, "LOGclv.png"),
       plot = hl_cleveland_log, width = 10, height = 10, dpi = 300)

#----6. Ratio analysis (HumHLr/BaseHLr) -----
ratio_clv <- hl_data %>%
  select(island, archipelago, baseHLr, global_ratio, local_ratio) %>%
  pivot_longer(cols = c(global_ratio, local_ratio), names_to = "proxy", values_to = "ratio") %>%
  mutate(
    proxy_label = recode(proxy, global_ratio = "Global proxy", local_ratio = "Local proxy"),
    archipelago = factor(archipelago, levels = sort(unique(archipelago)))) %>%
  arrange(archipelago, island) %>%
  group_by(archipelago) %>%
  mutate(island_f = factor(island, levels = rev(unique(island)))) %>%
  ungroup() %>%
  filter(is.finite(ratio), ratio > 0)

#summary stats
ratio_summary_overall <- ratio_clv %>%
  group_by(proxy_label) %>%
  summarise(
    n_islands = n(),
    mean_ratio = mean(ratio, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE),
    min_ratio = min(ratio, na.rm = TRUE),
    max_ratio = max(ratio, na.rm = TRUE),
    .groups = "drop")

ratio_summary_archipelago <- ratio_clv %>%
  group_by(archipelago, proxy_label) %>%
  summarise(
    n_islands = n(),
    mean_ratio = mean(ratio, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE),
    min_ratio = min(ratio, na.rm = TRUE),
    max_ratio = max(ratio, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(archipelago, proxy_label)

write_xlsx(
  list(
    ratio_values = ratio_clv,
    ratio_summary_overall = ratio_summary_overall,
    ratio_summary_archipelago = ratio_summary_archipelago),
  file.path(results_dir, "humhl_basehl_ratios.xlsx"))

#boxplots

#log
ratio_boxplot_log <- ggplot(ratio_clv, aes(x = proxy_label, y = ratio, fill = proxy_label)) +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey35", linewidth = 0.5) +
  geom_boxplot(width = 0.6, linewidth = 0.35, alpha = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, size = 1.8, color = "grey30", alpha = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = hl_cols) +
  labs(x = NULL, y = expression("HumHLr / BaseHLr (log"[10]*" scale)")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.text.x  = element_text(margin = margin(t = 3)),
    axis.text.y = element_text(margin = margin(r = 3)),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = NA))

print(ratio_boxplot_log)
ggsave(file.path(results_dir, "LOGboxplot_ratio.png"),
       plot = ratio_boxplot_log, width = 7, height = 7, dpi = 300)

#lin
ratio_boxplot_lin <- ratio_boxplot_log +
  scale_y_continuous() +
  labs(y = "HumHLr / BaseHLr")

print(ratio_boxplot_lin)
ggsave(file.path(results_dir, "LIN_boxplot_ratio.png"),
       plot = ratio_boxplot_lin, width = 7, height = 7, dpi = 300)

#Cleveland plots

#log
ratio_cleveland_log <- ggplot(
  ratio_clv, aes(x = ratio, y = island_f, color = proxy_label, fill = proxy_label, shape = proxy_label)) +
  geom_vline(xintercept = 1, linetype = "solid", color = "grey35", linewidth = 0.5) +
  geom_point(size = 2.7, alpha = 0.85, stroke = 0.3) +
  scale_x_log10() +
  scale_shape_manual(values = c("Global proxy" = 17, "Local proxy" = 25)) +
  scale_color_manual(values = hl_cols) +
  scale_fill_manual(values = hl_cols) +
  facet_grid(archipelago ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("HumHLr / BaseHLr  (log"[10]*" scale)"), y = NULL, color = NULL, fill = NULL, shape = NULL) +
  cleveland_theme +
  theme(panel.spacing.y = unit(0, "pt"))

print(ratio_cleveland_log)
ggsave(file.path(results_dir, "LOGclv_ratio.png"),
       plot = ratio_cleveland_log, width = 10, height = 10, dpi = 300)

#lin
ratio_cleveland_lin <- ratio_cleveland_log +
  scale_x_continuous() +
  labs(x = "HumHLr / BaseHLr")

print(ratio_cleveland_lin)
ggsave(file.path(results_dir, "LINclv_ratio.png"),
       plot = ratio_cleveland_lin, width = 10, height = 10, dpi = 300)
