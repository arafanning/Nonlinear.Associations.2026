# Nonlinear Associations Between Cannabis Use and Neurocognition Among People Living With HIV
# Script: 08_tables_figures.R
# Purpose: 08 Tables Figures
#

# ============================================================
# 4. Figure: Distribution of past-year cannabis use by category
# ============================================================

# ------------------------------------------------------------
# 4a. Packages for figure
# ------------------------------------------------------------

packages_needed_fig <- c("dplyr", "ggplot2", "scales", "patchwork")

for (pkg in packages_needed_fig) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)


# ------------------------------------------------------------
# 4b. Create figure-specific data frame
# Does NOT overwrite existing variables in data
# ------------------------------------------------------------

data_fig_cannabis_12m <- data %>%
  mutate(
    cannabis_12m_group = factor(
      du_mar4_12m_aBin,
      levels = c(0, 1, 2),
      labels = c("None", "Low", "High")
    ),
    
    cannabis_12m_log10p1 = log10(du_mar4_12m_a + 1),
    
    cannabis_12m_user_status = ifelse(
      du_mar4_12m_a > 0,
      "Any past-year use",
      "No past-year use"
    )
  )


# ------------------------------------------------------------
# 4c. Summary table for checking n, median, and range
# ------------------------------------------------------------

cannabis_12m_bin_summary <- data_fig_cannabis_12m %>%
  group_by(cannabis_12m_group) %>%
  summarise(
    n = n(),
    min_g = min(du_mar4_12m_a, na.rm = TRUE),
    median_g = median(du_mar4_12m_a, na.rm = TRUE),
    max_g = max(du_mar4_12m_a, na.rm = TRUE),
    .groups = "drop"
  )

cannabis_12m_bin_summary


# ------------------------------------------------------------
# 4d. Create x-axis labels with n and range
# ------------------------------------------------------------

x_labels_12m <- cannabis_12m_bin_summary %>%
  mutate(
    label = paste0(
      cannabis_12m_group,
      "\n",
      "n = ", n,
      "\n",
      min_g, "–", max_g, " g"
    )
  )

x_label_vector_12m <- setNames(
  x_labels_12m$label,
  x_labels_12m$cannabis_12m_group
)


# ------------------------------------------------------------
# 4e. Panel A: Distribution among cannabis users only
# ------------------------------------------------------------

fig_12m_nonzero_distribution <- data_fig_cannabis_12m %>%
  filter(du_mar4_12m_a > 0) %>%
  ggplot(aes(x = du_mar4_12m_a)) +
  geom_histogram(
    bins = 18,
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = median_12m,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_x_log10(
    breaks = c(0.1, 0.5, 180, 288, 10000),
    labels = c("0.1", "0.5", "180", "288", "10,000")
  ) +
  labs(
    x = "Past-Year Cannabis Use (g.)",
    y = "Number of Participants",
    title = "A. Distribution Among Past-Year Cannabis Users"
  ) +
  theme_classic(base_size = 12)

fig_12m_nonzero_distribution


# ------------------------------------------------------------
# 4f. Panel B: Individual values within ordered cannabis-use levels
# Uses log10(x + 1) so zero-use participants can be displayed
# ------------------------------------------------------------

fig_12m_bins_jitter <- ggplot(
  data_fig_cannabis_12m,
  aes(x = cannabis_12m_group, y = cannabis_12m_log10p1)
) +
  geom_boxplot(
    width = 0.35,
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.70,
    size = 2
  ) +
  geom_hline(
    yintercept = log10(median_12m + 1),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_x_discrete(
    labels = x_label_vector_12m
  ) +
  scale_y_continuous(
    breaks = log10(c(0, 0.5, 180, 288, 10000) + 1),
    labels = c("0", "0.5", "180", "288", "10,000")
  ) +
  labs(
    x = "Past-Year Cannabis Use Level",
    y = "Past-Year Cannabis Use (g)",
    title = "B. Past-Year Cannabis Use Within Ordered Categories"
  ) +
  theme_classic(base_size = 12)

fig_12m_bins_jitter


# ------------------------------------------------------------
# 4g. Combine panels
# ------------------------------------------------------------

figure_12m_cannabis_distribution <- 
  fig_12m_nonzero_distribution / fig_12m_bins_jitter +
  plot_layout(heights = c(1, 1.2))

figure_12m_cannabis_distribution


# ------------------------------------------------------------
# 4h. Save figure
# ------------------------------------------------------------

ggsave(
  filename = "fig_cannabis_12m_distribution.png",
  plot = figure_12m_cannabis_distribution,
  width = 7,
  height = 8,
  dpi = 300
)

ggsave(
  filename = "fig_cannabis_12m_distribution.pdf",
  plot = figure_12m_cannabis_distribution,
  width = 7,
  height = 8
)

