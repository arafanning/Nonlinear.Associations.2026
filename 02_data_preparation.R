# Nonlinear Associations Between Cannabis Use and Neurocognition Among People Living With HIV

# ============================================================
# CANDIDACY ANALYSIS SYNTAX
# Neurocognitive outcomes, HIV disease severity, cannabis bins
# ============================================================

# Load data
data <- read_sav("data/HIV+113.sav")
View(data)

# ============================================================
# 0. Packages
# ============================================================

packages_needed <- c("haven", "lavaan", "lmtest", "sandwich")

for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(haven)
library(lavaan)
library(lmtest)
library(sandwich)

# ============================================================
# 1. Create neurocognitive test z-scores and domain composites
# ============================================================

# ------------------------------------------------------------
# 1a. Reverse-code selected variables
# Higher values should reflect better performance across tests.
# ------------------------------------------------------------

data$mst_ct_rev <- -1 * data$mst_s_ct1
data$gpt_time_rev <- -1 * data$gpt_nd_ttime
data$wcpe_rev <- -1 * data$wc_pe_rs


# ------------------------------------------------------------
# 1b. Create z-score variables
# ------------------------------------------------------------

data$bvtr_z   <- as.numeric(scale(data$bv_tr_rs))
data$hvtt_z   <- as.numeric(scale(data$hv_tt_rs))
data$hvdr_z   <- as.numeric(scale(data$hv_dr_rs))
data$bvdr_z   <- as.numeric(scale(data$bv_dr_rs))
data$mstct_z  <- as.numeric(scale(data$mst_ct_rev))
data$gpttim_z <- as.numeric(scale(data$gpt_time_rev))
data$wcpe_z   <- as.numeric(scale(data$wcpe_rev))
data$igtnt_z  <- as.numeric(scale(data$igt_nt_rs))
data$cd_z     <- as.numeric(scale(data$wai_cd_rs))
data$ss_z     <- as.numeric(scale(data$wai_ss_rs))


# ------------------------------------------------------------
# 1c. Check z-score variables
# ------------------------------------------------------------

test_z_vars <- c(
  "bvtr_z", "hvtt_z", "hvdr_z", "bvdr_z",
  "mstct_z", "gpttim_z", "wcpe_z", "igtnt_z",
  "cd_z", "ss_z"
)

head(data[, test_z_vars])

sapply(data[, test_z_vars], function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )
})


# ------------------------------------------------------------
# 1d. Create neurocognitive domain composites
# ------------------------------------------------------------

data$learning_z <- rowMeans(
  data[, c("bvtr_z", "hvtt_z")],
  na.rm = TRUE
)

data$memory_z <- rowMeans(
  data[, c("hvdr_z", "bvdr_z")],
  na.rm = TRUE
)

data$motor_z <- rowMeans(
  data[, c("mstct_z", "gpttim_z")],
  na.rm = TRUE
)

data$executive_z <- rowMeans(
  data[, c("wcpe_z", "igtnt_z")],
  na.rm = TRUE
)

data$processing_speed_z <- rowMeans(
  data[, c("cd_z", "ss_z")],
  na.rm = TRUE
)


# ------------------------------------------------------------
# 1e. Create global neurocognitive composite
# Require at least 9 of 10 tests.
# ------------------------------------------------------------

data$num_tests_available <- rowSums(!is.na(data[, test_z_vars]))

data$global_z <- rowMeans(
  data[, test_z_vars],
  na.rm = TRUE
)

data$global_z[data$num_tests_available < 9] <- NA


# ------------------------------------------------------------
# 1f. Check domain and global composites
# ------------------------------------------------------------

outcome_vars <- c(
  "global_z",
  "learning_z",
  "memory_z",
  "motor_z",
  "executive_z",
  "processing_speed_z"
)

summary(data[, outcome_vars])

colSums(is.na(data[, outcome_vars]))

data$num_domains_missing <- rowSums(is.na(data[, c(
  "learning_z",
  "memory_z",
  "motor_z",
  "executive_z",
  "processing_speed_z"
)]))

table(data$num_domains_missing)


# Check domains calculated from only one available test
sum(rowSums(!is.na(data[, c("bvtr_z", "hvtt_z")])) == 1)       # Learning
sum(rowSums(!is.na(data[, c("hvdr_z", "bvdr_z")])) == 1)       # Memory
sum(rowSums(!is.na(data[, c("mstct_z", "gpttim_z")])) == 1)    # Motor
sum(rowSums(!is.na(data[, c("wcpe_z", "igtnt_z")])) == 1)      # Executive
sum(rowSums(!is.na(data[, c("cd_z", "ss_z")])) == 1)           # Processing speed


# ============================================================
# 2. HIV disease-severity factor
# ============================================================

# ------------------------------------------------------------
# 2a. Prepare HIV disease-severity indicators
# Higher DiseaseSev should reflect worse HIV disease severity.
# ------------------------------------------------------------

data$cd4sqrt <- -sqrt(data$lab_thel)
data$nadirsqrt <- -sqrt(data$mhq_4_cd4_lowest)
data$log_vl2 <- log1p(data$lab_hiv)


# ------------------------------------------------------------
# 2b. Estimate latent disease-severity factor
# ------------------------------------------------------------

DiseaseSev_model <- '
  DiseaseSev =~ nadirsqrt + cd4sqrt + log_vl2
  cd4sqrt ~~ 0.001*cd4sqrt
'

fit <- lavaan::cfa(
  DiseaseSev_model,
  data = data,
  estimator = "MLR",
  missing = "fiml"
)

summary(fit, fit.measures = TRUE, standardized = TRUE)

data$DiseaseSev <- as.numeric(lavaan::lavPredict(fit))

data$DiseaseSev_c <- as.numeric(
  scale(data$DiseaseSev, center = TRUE, scale = FALSE)
)

summary(data$DiseaseSev)
summary(data$DiseaseSev_c)


# ============================================================
# 3. Past-year cannabis-use bins
# Raw variable: du_mar4_12m_a
# 0 = no use
# 1 = low use: >0 and <= median of non-zero past-year grams
# 2 = high use: > median of non-zero past-year grams
# ============================================================

summary(data$du_mar4_12m_a)
hist(data$du_mar4_12m_a)

nonzero_vals_12m <- data$du_mar4_12m_a[data$du_mar4_12m_a > 0]
median_12m <- median(nonzero_vals_12m, na.rm = TRUE)

median_12m

data$du_mar4_12m_aBin <- NA_real_

data$du_mar4_12m_aBin[data$du_mar4_12m_a == 0] <- 0

data$du_mar4_12m_aBin[
  data$du_mar4_12m_a > 0 &
    data$du_mar4_12m_a <= median_12m
] <- 1

data$du_mar4_12m_aBin[
  data$du_mar4_12m_a > median_12m
] <- 2

table(data$du_mar4_12m_aBin, useNA = "ifany")

range_low_12m <- range(
  data$du_mar4_12m_a[data$du_mar4_12m_aBin == 1],
  na.rm = TRUE
)

range_high_12m <- range(
  data$du_mar4_12m_a[data$du_mar4_12m_aBin == 2],
  na.rm = TRUE
)

range_low_12m
range_high_12m

max_low_12m <- max(
  data$du_mar4_12m_a[data$du_mar4_12m_aBin == 1],
  na.rm = TRUE
)

min_high_12m <- min(
  data$du_mar4_12m_a[data$du_mar4_12m_aBin == 2],
  na.rm = TRUE
)

gap_12m <- min_high_12m - max_low_12m

max_low_12m
min_high_12m
gap_12m


# Ordered past-year cannabis variable for polynomial contrasts
data$du_mar4_12m_aBin_ord <- ordered(
  data$du_mar4_12m_aBin,
  levels = c(0, 1, 2),
  labels = c("none", "low", "high")
)

contrasts(data$du_mar4_12m_aBin_ord) <- contr.poly(3)

table(data$du_mar4_12m_aBin_ord, useNA = "ifany")
contrasts(data$du_mar4_12m_aBin_ord)

# ===================================================================
##Added

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



# ============================================================
# 4. Past-30-day cannabis-use bins: median split
# Raw variable: du_mar6_30d_a
# 0 = no use
# 1 = low recent use: >0 and <= median of non-zero 30-day grams
# 2 = high recent use: > median of non-zero 30-day grams
# ============================================================

summary(data$du_mar6_30d_a)
hist(data$du_mar6_30d_a)

nonzero_vals_30d <- data$du_mar6_30d_a[data$du_mar6_30d_a > 0]
median_30d <- median(nonzero_vals_30d, na.rm = TRUE)

median_30d

data$du_mar6_30d_aBin <- NA_real_

data$du_mar6_30d_aBin[data$du_mar6_30d_a == 0] <- 0

data$du_mar6_30d_aBin[
  data$du_mar6_30d_a > 0 &
    data$du_mar6_30d_a <= median_30d
] <- 1

data$du_mar6_30d_aBin[
  data$du_mar6_30d_a > median_30d
] <- 2

table(data$du_mar6_30d_aBin, useNA = "ifany")

range_low_30d <- range(
  data$du_mar6_30d_a[data$du_mar6_30d_aBin == 1],
  na.rm = TRUE
)

range_high_30d <- range(
  data$du_mar6_30d_a[data$du_mar6_30d_aBin == 2],
  na.rm = TRUE
)

range_low_30d
range_high_30d

max_low_30d <- max(
  data$du_mar6_30d_a[data$du_mar6_30d_aBin == 1],
  na.rm = TRUE
)

min_high_30d <- min(
  data$du_mar6_30d_a[data$du_mar6_30d_aBin == 2],
  na.rm = TRUE
)

gap_30d <- min_high_30d - max_low_30d

max_low_30d
min_high_30d
gap_30d


# Ordered 30-day cannabis variable
data$du_mar6_30d_aBin_ord <- ordered(
  data$du_mar6_30d_aBin,
  levels = c(0, 1, 2),
  labels = c("none", "low_recent", "high_recent")
)

contrasts(data$du_mar6_30d_aBin_ord) <- contr.poly(3)

table(data$du_mar6_30d_aBin_ord, useNA = "ifany")
contrasts(data$du_mar6_30d_aBin_ord)


# ============================================================
# 5. Past-30-day cannabis-use bins: 30g threshold sensitivity
# Raw variable: du_mar6_30d_a
# 0 = no use
# 1 = lower recent use: >0 and <=30g
# 2 = heavier recent use: >30g
# ============================================================

thresholds_30d <- c(10, 20, 30, 50, 75, 100, 110)

threshold_check_30d <- sapply(thresholds_30d, function(x) {
  c(
    n_at_or_above = sum(data$du_mar6_30d_a >= x, na.rm = TRUE),
    n_below_but_nonzero = sum(data$du_mar6_30d_a > 0 & data$du_mar6_30d_a < x, na.rm = TRUE),
    n_zero = sum(data$du_mar6_30d_a == 0, na.rm = TRUE)
  )
})

threshold_check_30d

data$du_mar6_30d_aBin30_lowincl <- NA_real_

data$du_mar6_30d_aBin30_lowincl[data$du_mar6_30d_a == 0] <- 0

data$du_mar6_30d_aBin30_lowincl[
  data$du_mar6_30d_a > 0 &
    data$du_mar6_30d_a <= 30
] <- 1

data$du_mar6_30d_aBin30_lowincl[
  data$du_mar6_30d_a > 30
] <- 2

table(data$du_mar6_30d_aBin30_lowincl, useNA = "ifany")

range_low_30d_lowincl <- range(
  data$du_mar6_30d_a[data$du_mar6_30d_aBin30_lowincl == 1],
  na.rm = TRUE
)

range_high_30d_lowincl <- range(
  data$du_mar6_30d_a[data$du_mar6_30d_aBin30_lowincl == 2],
  na.rm = TRUE
)

range_low_30d_lowincl
range_high_30d_lowincl

max_low_30d_lowincl <- max(
  data$du_mar6_30d_a[data$du_mar6_30d_aBin30_lowincl == 1],
  na.rm = TRUE
)

min_high_30d_lowincl <- min(
  data$du_mar6_30d_a[data$du_mar6_30d_aBin30_lowincl == 2],
  na.rm = TRUE
)

gap_30d_lowincl <- min_high_30d_lowincl - max_low_30d_lowincl

max_low_30d_lowincl
min_high_30d_lowincl
gap_30d_lowincl


data$du_mar6_30d_aBin30_lowincl_ord <- ordered(
  data$du_mar6_30d_aBin30_lowincl,
  levels = c(0, 1, 2),
  labels = c("none", "lower_recent_0to30", "heavier_recent_over30")
)

contrasts(data$du_mar6_30d_aBin30_lowincl_ord) <- contr.poly(3)

table(data$du_mar6_30d_aBin30_lowincl_ord, useNA = "ifany")
contrasts(data$du_mar6_30d_aBin30_lowincl_ord)


# ============================================================
# 6. Lifetime cannabis-use bins: median split
# Raw variable: du_mar2_life_a
# 0 = no use
# 1 = low lifetime use: >0 and <= median of non-zero lifetime grams
# 2 = high lifetime use: > median of non-zero lifetime grams
# ============================================================

summary(data$du_mar2_life_a)
hist(data$du_mar2_life_a)

nonzero_vals_life <- data$du_mar2_life_a[data$du_mar2_life_a > 0]
median_life <- median(nonzero_vals_life, na.rm = TRUE)

median_life

data$du_mar2_life_aBin <- NA_real_

data$du_mar2_life_aBin[data$du_mar2_life_a == 0] <- 0

data$du_mar2_life_aBin[
  data$du_mar2_life_a > 0 &
    data$du_mar2_life_a <= median_life
] <- 1

data$du_mar2_life_aBin[
  data$du_mar2_life_a > median_life
] <- 2

table(data$du_mar2_life_aBin, useNA = "ifany")

range_low_life <- range(
  data$du_mar2_life_a[data$du_mar2_life_aBin == 1],
  na.rm = TRUE
)

range_high_life <- range(
  data$du_mar2_life_a[data$du_mar2_life_aBin == 2],
  na.rm = TRUE
)

range_low_life
range_high_life

max_low_life <- max(
  data$du_mar2_life_a[data$du_mar2_life_aBin == 1],
  na.rm = TRUE
)

min_high_life <- min(
  data$du_mar2_life_a[data$du_mar2_life_aBin == 2],
  na.rm = TRUE
)

gap_life <- min_high_life - max_low_life

max_low_life
min_high_life
gap_life


data$du_mar2_life_aBin_ord <- ordered(
  data$du_mar2_life_aBin,
  levels = c(0, 1, 2),
  labels = c("none", "low_lifetime", "high_lifetime")
)

contrasts(data$du_mar2_life_aBin_ord) <- contr.poly(3)

table(data$du_mar2_life_aBin_ord, useNA = "ifany")
contrasts(data$du_mar2_life_aBin_ord)


# ============================================================
# 7. Lifetime cannabis-use bins: q75 sensitivity split
# Raw variable: du_mar2_life_a
# 0 = no use
# 1 = lower/moderate lifetime use: >0 and < q75
# 2 = high lifetime use: >= q75
# ============================================================

quantile(
  data$du_mar2_life_a[data$du_mar2_life_a > 0],
  probs = c(0, .10, .25, .33, .50, .67, .75, .90, .95, 1),
  na.rm = TRUE
)

thresholds_life <- c(100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)

threshold_check_life <- sapply(thresholds_life, function(x) {
  c(
    n_at_or_above = sum(data$du_mar2_life_a >= x, na.rm = TRUE),
    n_below_but_nonzero = sum(data$du_mar2_life_a > 0 & data$du_mar2_life_a < x, na.rm = TRUE),
    n_zero = sum(data$du_mar2_life_a == 0, na.rm = TRUE)
  )
})

threshold_check_life

q75_life <- quantile(
  data$du_mar2_life_a[data$du_mar2_life_a > 0],
  probs = .75,
  na.rm = TRUE
)

q75_life

data$du_mar2_life_aBin_q75 <- NA_real_

data$du_mar2_life_aBin_q75[data$du_mar2_life_a == 0] <- 0

data$du_mar2_life_aBin_q75[
  data$du_mar2_life_a > 0 &
    data$du_mar2_life_a < q75_life
] <- 1

data$du_mar2_life_aBin_q75[
  data$du_mar2_life_a >= q75_life
] <- 2

table(data$du_mar2_life_aBin_q75, useNA = "ifany")

range_lowmod_life_q75 <- range(
  data$du_mar2_life_a[data$du_mar2_life_aBin_q75 == 1],
  na.rm = TRUE
)

range_high_life_q75 <- range(
  data$du_mar2_life_a[data$du_mar2_life_aBin_q75 == 2],
  na.rm = TRUE
)

range_lowmod_life_q75
range_high_life_q75

max_lowmod_life_q75 <- max(
  data$du_mar2_life_a[data$du_mar2_life_aBin_q75 == 1],
  na.rm = TRUE
)

min_high_life_q75 <- min(
  data$du_mar2_life_a[data$du_mar2_life_aBin_q75 == 2],
  na.rm = TRUE
)

gap_life_q75 <- min_high_life_q75 - max_lowmod_life_q75

max_lowmod_life_q75
min_high_life_q75
gap_life_q75


data$du_mar2_life_aBin_q75_ord <- ordered(
  data$du_mar2_life_aBin_q75,
  levels = c(0, 1, 2),
  labels = c("none", "lower_moderate_lifetime", "high_lifetime")
)

contrasts(data$du_mar2_life_aBin_q75_ord) <- contr.poly(3)

table(data$du_mar2_life_aBin_q75_ord, useNA = "ifany")
contrasts(data$du_mar2_life_aBin_q75_ord)


# Continuous log-transformed lifetime exposure
data$du_mar2_life_a_log1p <- log1p(data$du_mar2_life_a)

summary(data$du_mar2_life_a_log1p)
hist(data$du_mar2_life_a_log1p)


# ============================================================
# 8. Decision table comparing cannabis binning options
# ============================================================

cannabis_decision_table <- data.frame(
  exposure = c(
    "Past-year median split",
    "30-day median split",
    "30-day 30g threshold",
    "Lifetime median split",
    "Lifetime q75 split"
  ),
  
  raw_variable = c(
    "du_mar4_12m_a",
    "du_mar6_30d_a",
    "du_mar6_30d_a",
    "du_mar2_life_a",
    "du_mar2_life_a"
  ),
  
  bin_variable = c(
    "du_mar4_12m_aBin",
    "du_mar6_30d_aBin",
    "du_mar6_30d_aBin30_lowincl",
    "du_mar2_life_aBin",
    "du_mar2_life_aBin_q75"
  ),
  
  n_none = c(
    sum(data$du_mar4_12m_aBin == 0, na.rm = TRUE),
    sum(data$du_mar6_30d_aBin == 0, na.rm = TRUE),
    sum(data$du_mar6_30d_aBin30_lowincl == 0, na.rm = TRUE),
    sum(data$du_mar2_life_aBin == 0, na.rm = TRUE),
    sum(data$du_mar2_life_aBin_q75 == 0, na.rm = TRUE)
  ),
  
  n_low_or_lower = c(
    sum(data$du_mar4_12m_aBin == 1, na.rm = TRUE),
    sum(data$du_mar6_30d_aBin == 1, na.rm = TRUE),
    sum(data$du_mar6_30d_aBin30_lowincl == 1, na.rm = TRUE),
    sum(data$du_mar2_life_aBin == 1, na.rm = TRUE),
    sum(data$du_mar2_life_aBin_q75 == 1, na.rm = TRUE)
  ),
  
  n_high_or_heavier = c(
    sum(data$du_mar4_12m_aBin == 2, na.rm = TRUE),
    sum(data$du_mar6_30d_aBin == 2, na.rm = TRUE),
    sum(data$du_mar6_30d_aBin30_lowincl == 2, na.rm = TRUE),
    sum(data$du_mar2_life_aBin == 2, na.rm = TRUE),
    sum(data$du_mar2_life_aBin_q75 == 2, na.rm = TRUE)
  ),
  
  lower_min = c(
    range_low_12m[1],
    range_low_30d[1],
    range_low_30d_lowincl[1],
    range_low_life[1],
    range_lowmod_life_q75[1]
  ),
  
  lower_max = c(
    range_low_12m[2],
    range_low_30d[2],
    range_low_30d_lowincl[2],
    range_low_life[2],
    range_lowmod_life_q75[2]
  ),
  
  higher_min = c(
    range_high_12m[1],
    range_high_30d[1],
    range_high_30d_lowincl[1],
    range_high_life[1],
    range_high_life_q75[1]
  ),
  
  higher_max = c(
    range_high_12m[2],
    range_high_30d[2],
    range_high_30d_lowincl[2],
    range_high_life[2],
    range_high_life_q75[2]
  ),
  
  gap = c(
    gap_12m,
    gap_30d,
    gap_30d_lowincl,
    gap_life,
    gap_life_q75
  )
)

cannabis_decision_table


# ============================================================
# 9. Primary candidacy models
# Primary cannabis predictor: du_mar4_12m_aBin_ord
# ============================================================

# Primary model without additional covariates
# Use this first to make sure the core model runs.

model_global_primary <- lm(
  global_z ~ DiseaseSev_c * du_mar4_12m_aBin_ord,
  data = data
)

summary(model_global_primary)

lmtest::coeftest(
  model_global_primary,
  vcov. = sandwich::vcovHC(model_global_primary, type = "HC3")
)


# Run same model across all neurocognitive outcomes

models_12m_primary <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ DiseaseSev_c * du_mar4_12m_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_12m_primary


# ============================================================
# 10. Optional sensitivity models
# ============================================================

# 30-day median split sensitivity

models_30d_median <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ DiseaseSev_c * du_mar6_30d_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_30d_median


# 30-day 30g threshold sensitivity

models_30d_30g <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ DiseaseSev_c * du_mar6_30d_aBin30_lowincl_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_30d_30g


# Lifetime median split sensitivity

models_life_median <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ DiseaseSev_c * du_mar2_life_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_life_median


# Lifetime q75 sensitivity

models_life_q75 <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ DiseaseSev_c * du_mar2_life_aBin_q75_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_life_q75


# Lifetime log1p continuous sensitivity

models_life_log <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ DiseaseSev_c * du_mar2_life_a_log1p")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_life_log


# ============================================================
# 11. Final decision notes
# ============================================================

# Primary candidacy exposure:
# du_mar4_12m_aBin_ord
#
# Rationale:
# Past-year cannabis grams are retained as the primary cannabis exposure
# because the binning captures none / low / high chronicity of use and is
# most aligned with the candidacy paper's chronic exposure question.
#
# Recent-use sensitivity:
# du_mar6_30d_aBin_ord
# du_mar6_30d_aBin30_lowincl_ord
#
# Lifetime-use sensitivity:
# du_mar2_life_aBin_ord
# du_mar2_life_aBin_q75_ord
# du_mar2_life_a_log1p
#
# Important:
# Raw variables were preserved:
# du_mar4_12m_a
# du_mar6_30d_a
# du_mar2_life_a
#
# Derived variables were newly created.
# ============================================================


# ============================================================
# 12. Individual HIV disease-severity marker sensitivity models
# Purpose: Probe whether current CD4, nadir CD4, or viral load
# show associations that may be obscured by the latent factor.
# ============================================================

# Center individual disease-severity indicators
# These are already transformed so higher values reflect worse severity:
# cd4sqrt   = reversed current CD4
# nadirsqrt = reversed nadir CD4
# log_vl2   = log viral load

data$cd4sqrt_c <- as.numeric(scale(data$cd4sqrt, center = TRUE, scale = FALSE))
data$nadirsqrt_c <- as.numeric(scale(data$nadirsqrt, center = TRUE, scale = FALSE))
data$log_vl2_c <- as.numeric(scale(data$log_vl2, center = TRUE, scale = FALSE))

summary(data[, c("cd4sqrt_c", "nadirsqrt_c", "log_vl2_c")])

# ------------------------------------------------------------
# 12a. Main-effect models for individual HIV markers
# ------------------------------------------------------------

hiv_marker_vars <- c("cd4sqrt_c", "nadirsqrt_c", "log_vl2_c")

models_hiv_markers_main <- list()

for (marker in hiv_marker_vars) {
  
  models_hiv_markers_main[[marker]] <- lapply(outcome_vars, function(y) {
    
    model_formula <- as.formula(
      paste0(y, " ~ ", marker)
    )
    
    model <- lm(model_formula, data = data)
    
    robust_results <- lmtest::coeftest(
      model,
      vcov. = sandwich::vcovHC(model, type = "HC3")
    )
    
    list(
      outcome = y,
      marker = marker,
      model = model,
      robust_results = robust_results
    )
  })
}

models_hiv_markers_main

# ------------------------------------------------------------
# 12b. Individual HIV marker moderation models
# Primary cannabis predictor: past-year cannabis bins
# ------------------------------------------------------------

models_hiv_markers_x_12m <- list()

for (marker in hiv_marker_vars) {
  
  models_hiv_markers_x_12m[[marker]] <- lapply(outcome_vars, function(y) {
    
    model_formula <- as.formula(
      paste0(y, " ~ ", marker, " * du_mar4_12m_aBin_ord")
    )
    
    model <- lm(model_formula, data = data)
    
    robust_results <- lmtest::coeftest(
      model,
      vcov. = sandwich::vcovHC(model, type = "HC3")
    )
    
    list(
      outcome = y,
      marker = marker,
      model = model,
      robust_results = robust_results
    )
  })
}

models_hiv_markers_x_12m

# ------------------------------------------------------------
# 12c. Individual HIV marker moderation models
# Sensitivity cannabis predictor: past-30-day median split
# ------------------------------------------------------------

models_hiv_markers_x_30d <- list()

for (marker in hiv_marker_vars) {
  
  models_hiv_markers_x_30d[[marker]] <- lapply(outcome_vars, function(y) {
    
    model_formula <- as.formula(
      paste0(y, " ~ ", marker, " * du_mar6_30d_aBin_ord")
    )
    
    model <- lm(model_formula, data = data)
    
    robust_results <- lmtest::coeftest(
      model,
      vcov. = sandwich::vcovHC(model, type = "HC3")
    )
    
    list(
      outcome = y,
      marker = marker,
      model = model,
      robust_results = robust_results
    )
  })
}

models_hiv_markers_x_30d


# ============================================================
# 13. Extract robust model results into clean tables
# ============================================================

extract_robust_results <- function(model_list, model_name) {
  
  output <- data.frame()
  
  for (marker_name in names(model_list)) {
    
    for (i in seq_along(model_list[[marker_name]])) {
      
      outcome_i <- model_list[[marker_name]][[i]]$outcome
      marker_i <- model_list[[marker_name]][[i]]$marker
      robust_i <- model_list[[marker_name]][[i]]$robust_results
      
      temp <- data.frame(
        model_set = model_name,
        outcome = outcome_i,
        marker = marker_i,
        term = rownames(robust_i),
        estimate = robust_i[, "Estimate"],
        robust_se = robust_i[, "Std. Error"],
        t_value = robust_i[, "t value"],
        p_value = robust_i[, "Pr(>|t|)"],
        row.names = NULL
      )
      
      output <- rbind(output, temp)
    }
  }
  
  return(output)
}

hiv_marker_main_table <- extract_robust_results(
  model_list = models_hiv_markers_main,
  model_name = "Individual HIV marker main effects"
)

hiv_marker_x_12m_table <- extract_robust_results(
  model_list = models_hiv_markers_x_12m,
  model_name = "Individual HIV marker x past-year cannabis"
)

hiv_marker_x_30d_table <- extract_robust_results(
  model_list = models_hiv_markers_x_30d,
  model_name = "Individual HIV marker x 30-day cannabis"
)

hiv_marker_main_table
hiv_marker_x_12m_table
hiv_marker_x_30d_table

# ------------------------------------------------------------
# Pull only HIV marker main effects and interaction terms
# ------------------------------------------------------------

hiv_marker_x_12m_terms <- hiv_marker_x_12m_table[
  grepl("cd4sqrt_c|nadirsqrt_c|log_vl2_c", hiv_marker_x_12m_table$term),
]

hiv_marker_x_30d_terms <- hiv_marker_x_30d_table[
  grepl("cd4sqrt_c|nadirsqrt_c|log_vl2_c", hiv_marker_x_30d_table$term),
]

hiv_marker_x_12m_terms
hiv_marker_x_30d_terms

# ============================================================
# 14. Correlations among disease-severity indicators
# ============================================================

disease_marker_data <- data[, c(
  "cd4sqrt",
  "nadirsqrt",
  "log_vl2",
  "DiseaseSev"
)]

cor(
  disease_marker_data,
  use = "pairwise.complete.obs"
)

round(
  cor(disease_marker_data, use = "pairwise.complete.obs"),
  3
)

# Missingness in HIV disease markers
colSums(is.na(data[, c(
  "lab_thel",
  "mhq_4_cd4_lowest",
  "lab_hiv",
  "cd4sqrt",
  "nadirsqrt",
  "log_vl2",
  "DiseaseSev"
)]))

# Decision:
# The latent HIV disease-severity factor remains the primary disease-severity
# variable because it represents cumulative HIV burden using current CD4,
# nadir CD4, and viral load. However, because the latent factor was not
# consistently associated with neurocognitive outcomes and did not clearly
# moderate cannabis effects, I probed individual disease markers separately.
# These sensitivity models test whether current immune status, historical immune
# suppression, or current viral burden show distinct associations with cognition
# or cannabis-related cognitive patterns.

# ============================================================
# Check whether the latent HIV Disease Severity factor is supported
# ============================================================

# 1. Full CFA summary with fit indices, standardized loadings, and R-squared
summary(
  fit,
  fit.measures = TRUE,
  standardized = TRUE,
  rsquare = TRUE
)

# 2. Parameter estimates with standardized loadings
lavaan::parameterEstimates(
  fit,
  standardized = TRUE
)

# 3. Key fit measures
lavaan::fitMeasures(
  fit,
  c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
)

# 4. Standardized solution only
lavaan::standardizedSolution(fit)

# 5. Inspect model-implied components
lavaan::inspect(fit, "std")
lavaan::inspect(fit, "theta")   # residual variances
lavaan::inspect(fit, "cov.lv")  # latent variance/covariance

# ============================================================
# Check whether DiseaseSev factor scores are redundant with indicators
# ============================================================

# Correlations among individual markers and factor score
disease_marker_cor <- cor(
  data[, c("cd4sqrt", "nadirsqrt", "log_vl2", "DiseaseSev")],
  use = "pairwise.complete.obs"
)

round(disease_marker_cor, 3)

# Scatterplots: DiseaseSev against each indicator
plot(
  data$cd4sqrt,
  data$DiseaseSev,
  xlab = "Reversed current CD4 sqrt",
  ylab = "DiseaseSev factor score",
  main = "DiseaseSev vs. Current CD4"
)

plot(
  data$nadirsqrt,
  data$DiseaseSev,
  xlab = "Reversed nadir CD4 sqrt",
  ylab = "DiseaseSev factor score",
  main = "DiseaseSev vs. Nadir CD4"
)

plot(
  data$log_vl2,
  data$DiseaseSev,
  xlab = "Log viral load",
  ylab = "DiseaseSev factor score",
  main = "DiseaseSev vs. Viral Load"
)

# ============================================================
# Check whether DiseaseSev factor scores are redundant with indicators
# ============================================================

# Correlations among individual markers and factor score
disease_marker_cor <- cor(
  data[, c("cd4sqrt", "nadirsqrt", "log_vl2", "DiseaseSev")],
  use = "pairwise.complete.obs"
)

round(disease_marker_cor, 3)

# Scatterplots: DiseaseSev against each indicator
plot(
  data$cd4sqrt,
  data$DiseaseSev,
  xlab = "Reversed current CD4 sqrt",
  ylab = "DiseaseSev factor score",
  main = "DiseaseSev vs. Current CD4"
)

plot(
  data$nadirsqrt,
  data$DiseaseSev,
  xlab = "Reversed nadir CD4 sqrt",
  ylab = "DiseaseSev factor score",
  main = "DiseaseSev vs. Nadir CD4"
)

plot(
  data$log_vl2,
  data$DiseaseSev,
  xlab = "Log viral load",
  ylab = "DiseaseSev factor score",
  main = "DiseaseSev vs. Viral Load"
)

# ============================================================
# Decision logic
# ============================================================
#
# The latent DiseaseSev factor is more defensible if:
# 1. All three indicators load meaningfully on the factor.
# 2. Standardized loadings are not dominated by only one marker.
# 3. Residual variances are not strange, negative, or constrained in a way
#    that makes the factor unstable.
# 4. DiseaseSev is not nearly identical to one indicator, especially cd4sqrt.
#
# The latent DiseaseSev factor is questionable if:
# 1. DiseaseSev correlates ~1.00 with cd4sqrt.
# 2. Nadir CD4 and viral load contribute only modestly.
# 3. The factor score appears visually like a direct transformation of current CD4.
# 4. Individual marker models show a different pattern than the latent factor.
#
# If questionable:
# Keep DiseaseSev as the planned primary construct if needed for candidacy continuity,
# but report/probe individual disease markers as sensitivity analyses.
# Based on current results, nadir CD4 appears especially important.
# ============================================================


# ============================================================
# Nadir CD4 sensitivity model
# Tests whether historical immune suppression relates to cognition
# and whether it moderates past-year cannabis effects
# ============================================================

data$nadirsqrt_c <- as.numeric(
  scale(data$nadirsqrt, center = TRUE, scale = FALSE)
)

models_nadir_x_12m <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ nadirsqrt_c * du_mar4_12m_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_nadir_x_12m


# ============================================================
# Nadir CD4 sensitivity models
# Tests whether historical immune suppression relates to cognition
# and whether it moderates past-year cannabis effects
# ============================================================

# Center reversed nadir CD4 variable
# Reminder: nadirsqrt = -sqrt(nadir CD4), so higher values = worse historical immune suppression

data$nadirsqrt_c <- as.numeric(
  scale(data$nadirsqrt, center = TRUE, scale = FALSE)
)

summary(data$nadirsqrt_c)

# ------------------------------------------------------------
# Nadir CD4 main-effect sensitivity models
# ------------------------------------------------------------

models_nadir_main_12m <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ nadirsqrt_c + du_mar4_12m_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    formula = model_formula,
    model = model,
    robust_results = robust_results
  )
})

models_nadir_main_12m


# ------------------------------------------------------------
# Nadir CD4 moderation sensitivity models
# ------------------------------------------------------------

models_nadir_x_12m <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ nadirsqrt_c * du_mar4_12m_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    formula = model_formula,
    model = model,
    robust_results = robust_results
  )
})

models_nadir_x_12m

# ------------------------------------------------------------
# Function to extract robust results from model list
# ------------------------------------------------------------

extract_lmtest_results <- function(model_list, model_name) {
  
  output <- data.frame()
  
  for (i in seq_along(model_list)) {
    
    outcome_i <- model_list[[i]]$outcome
    robust_i <- model_list[[i]]$robust_results
    
    temp <- data.frame(
      model_set = model_name,
      outcome = outcome_i,
      term = rownames(robust_i),
      estimate = robust_i[, "Estimate"],
      robust_se = robust_i[, "Std. Error"],
      t_value = robust_i[, "t value"],
      p_value = robust_i[, "Pr(>|t|)"],
      row.names = NULL
    )
    
    output <- rbind(output, temp)
  }
  
  return(output)
}

nadir_main_table <- extract_lmtest_results(
  model_list = models_nadir_main_12m,
  model_name = "Nadir main effect + past-year cannabis"
)

nadir_interaction_table <- extract_lmtest_results(
  model_list = models_nadir_x_12m,
  model_name = "Nadir x past-year cannabis"
)

nadir_main_table
nadir_interaction_table

# ------------------------------------------------------------
# Key nadir CD4 terms
# ------------------------------------------------------------

nadir_main_key_terms <- nadir_main_table[
  grepl("nadirsqrt_c", nadir_main_table$term),
]

nadir_interaction_key_terms <- nadir_interaction_table[
  grepl("nadirsqrt_c", nadir_interaction_table$term),
]

nadir_main_key_terms
nadir_interaction_key_terms


# ============================================================
# Past-year cannabis + latent HIV disease severity
# Main-effect model
# ============================================================

models_DiseaseSev_12m_main <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ DiseaseSev_c + du_mar4_12m_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_DiseaseSev_12m_main

########################################
########################################
#######tested whether the past-year cannabis association persisted after accounting for latent HIV disease severity.
########################################

# ============================================================
# Past-year cannabis x latent HIV disease severity
# Moderation model
# ============================================================

models_DiseaseSev_x_12m <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ DiseaseSev_c * du_mar4_12m_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_DiseaseSev_x_12m


# ============================================================
# Nadir CD4 sensitivity model
# Main-effect model
# ============================================================

data$nadirsqrt_c <- as.numeric(
  scale(data$nadirsqrt, center = TRUE, scale = FALSE)
)

models_nadir_12m_main <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ nadirsqrt_c + du_mar4_12m_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_nadir_12m_main


# ============================================================
# Nadir CD4 x past-year cannabis
# Moderation sensitivity model
# ============================================================

models_nadir_x_12m <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(y, " ~ nadirsqrt_c * du_mar4_12m_aBin_ord")
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    model = model,
    robust_results = robust_results
  )
})

models_nadir_x_12m

# ============================================================
# Omnibus test: Does past-year cannabis improve model fit?
# Compares disease severity only vs disease severity + cannabis
# ============================================================

omnibus_12m <- lapply(outcome_vars, function(y) {
  
  model_base <- lm(
    as.formula(paste0(y, " ~ DiseaseSev_c")),
    data = data
  )
  
  model_cu <- lm(
    as.formula(paste0(y, " ~ DiseaseSev_c + du_mar4_12m_aBin_ord")),
    data = data
  )
  
  anova_result <- anova(model_base, model_cu)
  
  list(
    outcome = y,
    anova = anova_result
  )
})

omnibus_12m

# ============================================================
# Omnibus test: Does cannabis x disease severity improve model fit?
# ============================================================

omnibus_interaction_12m <- lapply(outcome_vars, function(y) {
  
  model_main <- lm(
    as.formula(paste0(y, " ~ DiseaseSev_c + du_mar4_12m_aBin_ord")),
    data = data
  )
  
  model_int <- lm(
    as.formula(paste0(y, " ~ DiseaseSev_c * du_mar4_12m_aBin_ord")),
    data = data
  )
  
  anova_result <- anova(model_main, model_int)
  
  list(
    outcome = y,
    anova = anova_result
  )
})

omnibus_interaction_12m

# ============================================================
# PREP CHARACTERISTIC VARIABLES FOR COVARIATE SCREENING
# Rule: Do not overwrite raw variables. Create new variables.
# ============================================================

# Install only if needed
if (!requireNamespace("labelled", quietly = TRUE)) {
  install.packages("labelled")
}

library(labelled)

# ------------------------------------------------------------
# 1. Create numeric versions of race/ethnicity source variables
#    Do not overwrite raw labelled variables.
# ------------------------------------------------------------

data$phq_4_latin_num <- as.numeric(labelled::remove_labels(data$phq_4_latin))
data$phq_5_race_num <- as.numeric(labelled::remove_labels(data$phq_5_race))

# ------------------------------------------------------------
# 2. Create race/ethnicity binary factor
#    1 = White non-Hispanic
#    0 = Minoritized
# ------------------------------------------------------------

data$race_eth_binary <- factor(
  ifelse(
    data$phq_5_race_num == 5 & data$phq_4_latin_num == 0,
    1,
    0
  ),
  levels = c(0, 1),
  labels = c("Minoritized", "White non-Hispanic")
)

table(data$race_eth_binary, useNA = "always")

# Numeric version for correlation screening only
data$race_eth_binary_covnum <- as.numeric(data$race_eth_binary) - 1

table(data$race_eth_binary_covnum, useNA = "always")


# ------------------------------------------------------------
# 3. Create sex variable from gender item
#    Based on your coding: phq_6_gender == 2 coded as Female
# ------------------------------------------------------------

table(data$phq_6_gender, useNA = "always")

data$sex <- factor(
  ifelse(data$phq_6_gender == 2, 1, 0),
  levels = c(0, 1),
  labels = c("Male", "Female")
)

table(data$sex, useNA = "always")

# Numeric version for correlation screening only
data$sex_covnum <- as.numeric(data$sex) - 1

table(data$sex_covnum, useNA = "always")


# ------------------------------------------------------------
# 4. Create degree numeric screening variable if needed
# ------------------------------------------------------------

data$phq_7_degree_covnum <- as.numeric(as.factor(data$phq_7_degree))

table(data$phq_7_degree, useNA = "always")
table(data$phq_7_degree_covnum, useNA = "always")

# ============================================================
# COVARIATE SCREENING
# ============================================================

library(psych)
library(dplyr)

# ------------------------------------------------------------
# 1. Define variable groups
# ------------------------------------------------------------

cannabis_vars <- c(
  "du_mar6_30d_a",
  "du_mar4_12m_a",
  "du_mar2_life_a",
  "du_alc4_12m_a",
  "du_alc6_30d_a"
)

# Use numeric versions for correlation screening
covariate_vars_screen <- c(
  "phq_2_age",
  "phq_7_degree_covnum",
  "sex_covnum",
  "race_eth_binary_covnum",
  "bai_total",
  "bdi_total",
  "sc_mde1",
  "sc_mde3",
  "asr_adh_pro_total"
)

# Use the correct processing speed variable name
neurocog_vars <- c(
  "learning_z",
  "memory_z",
  "motor_z",
  "executive_z",
  "processing_speed_z",
  "global_z"
)

# ------------------------------------------------------------
# 2. Check that all variables exist before subsetting
# ------------------------------------------------------------

all_vars_needed <- c(cannabis_vars, covariate_vars_screen, neurocog_vars)

missing_vars <- setdiff(all_vars_needed, names(data))

missing_vars

# ------------------------------------------------------------
# 3. Create screening dataset without altering raw data
# ------------------------------------------------------------

cov_screen_data <- data[, all_vars_needed]

str(cov_screen_data)


# ------------------------------------------------------------
# 4. Cannabis exposure x candidate covariates
# ------------------------------------------------------------

cannabis_covariate_corr <- psych::corr.test(
  cov_screen_data[, c(cannabis_vars, covariate_vars_screen)],
  method = "spearman",
  use = "pairwise.complete.obs",
  adjust = "none"
)

cannabis_covariate_r <- cannabis_covariate_corr$r[
  cannabis_vars,
  covariate_vars_screen
]

cannabis_covariate_p <- cannabis_covariate_corr$p[
  cannabis_vars,
  covariate_vars_screen
]

cannabis_covariate_r
cannabis_covariate_p

# ------------------------------------------------------------
# 5. Candidate covariates x neurocognitive outcomes
# ------------------------------------------------------------

covariate_neurocog_corr <- psych::corr.test(
  cov_screen_data[, c(covariate_vars_screen, neurocog_vars)],
  method = "spearman",
  use = "pairwise.complete.obs",
  adjust = "none"
)

covariate_neurocog_r <- covariate_neurocog_corr$r[
  covariate_vars_screen,
  neurocog_vars
]

covariate_neurocog_p <- covariate_neurocog_corr$p[
  covariate_vars_screen,
  neurocog_vars
]

covariate_neurocog_r
covariate_neurocog_p

# ------------------------------------------------------------
# 6. Create clean long-format correlation tables
# ------------------------------------------------------------

cannabis_covariate_table <- as.data.frame(as.table(cannabis_covariate_r))
names(cannabis_covariate_table) <- c("cannabis_variable", "covariate", "rho")

cannabis_covariate_table$p_value <- as.vector(cannabis_covariate_p)

cannabis_covariate_table <- cannabis_covariate_table %>%
  arrange(p_value)

cannabis_covariate_table


covariate_neurocog_table <- as.data.frame(as.table(covariate_neurocog_r))
names(covariate_neurocog_table) <- c("covariate", "neurocog_outcome", "rho")

covariate_neurocog_table$p_value <- as.vector(covariate_neurocog_p)

covariate_neurocog_table <- covariate_neurocog_table %>%
  arrange(p_value)

covariate_neurocog_table

# ------------------------------------------------------------
# 7. Flag candidate covariates associated with cannabis and cognition
# ------------------------------------------------------------

covariates_related_to_cannabis <- cannabis_covariate_table %>%
  filter(p_value < .10) %>%
  arrange(p_value)

covariates_related_to_neurocog <- covariate_neurocog_table %>%
  filter(p_value < .10) %>%
  arrange(p_value)

covariates_related_to_cannabis
covariates_related_to_neurocog

candidate_covariates_to_consider <- intersect(
  unique(covariates_related_to_cannabis$covariate),
  unique(covariates_related_to_neurocog$covariate)
)

candidate_covariates_to_consider

# ============================================================
# FINAL MODEL COVARIATE VERSIONS
# ============================================================

data$sex_covfac <- as.factor(data$sex)
data$race_eth_binary_covfac <- as.factor(data$race_eth_binary)
data$phq_7_degree_covfac <- as.factor(data$phq_7_degree)



# ============================================================
# COVARIATE MODEL TESTING
# Primary cannabis predictor: du_mar4_12m_aBin_ord
# Primary disease marker: DiseaseSev_c
# Rule: Do not overwrite raw variables
# ============================================================

library(lmtest)
library(sandwich)
library(dplyr)

# Make sure final factor covariates exist
data$sex_covfac <- as.factor(data$sex)
data$race_eth_binary_covfac <- as.factor(data$race_eth_binary)
data$phq_7_degree_covfac <- as.factor(data$phq_7_degree)

# Confirm outcome names
outcome_vars <- c(
  "global_z",
  "learning_z",
  "memory_z",
  "motor_z",
  "executive_z",
  "processing_speed_z"
)

# ------------------------------------------------------------
# Define covariate sets
# ------------------------------------------------------------

covariate_sets <- list(
  
  unadjusted = c(),
  
  age_only = c(
    "phq_2_age"
  ),
  
  primary_demographic = c(
    "phq_2_age",
    "phq_7_degree_covfac"
  ),
  
  expanded_demographic = c(
    "phq_2_age",
    "phq_7_degree_covfac",
    "sex_covfac",
    "race_eth_binary_covfac"
  ),
  
  psych_past_mde = c(
    "phq_2_age",
    "phq_7_degree_covfac",
    "sc_mde3"
  ),
  
  psych_adhd = c(
    "phq_2_age",
    "phq_7_degree_covfac",
    "asr_adh_pro_total"
  )
)

# ------------------------------------------------------------
# Function to run one model
# ------------------------------------------------------------

run_covariate_model <- function(outcome, covariates) {
  
  base_terms <- c(
    "DiseaseSev_c",
    "du_mar4_12m_aBin_ord"
  )
  
  all_terms <- c(base_terms, covariates)
  
  model_formula <- as.formula(
    paste(outcome, "~", paste(all_terms, collapse = " + "))
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = outcome,
    formula = model_formula,
    n = nobs(model),
    model = model,
    robust_results = robust_results
  )
}


# ------------------------------------------------------------
# Run all covariate sets across all outcomes
# ------------------------------------------------------------

covariate_model_results <- list()

for (set_name in names(covariate_sets)) {
  
  covariate_model_results[[set_name]] <- lapply(outcome_vars, function(y) {
    
    run_covariate_model(
      outcome = y,
      covariates = covariate_sets[[set_name]]
    )
  })
}

covariate_model_results
# ------------------------------------------------------------
# Extract robust model results into clean table
# ------------------------------------------------------------

extract_covariate_results <- function(model_results_list) {
  
  output <- data.frame()
  
  for (set_name in names(model_results_list)) {
    
    for (i in seq_along(model_results_list[[set_name]])) {
      
      model_i <- model_results_list[[set_name]][[i]]
      robust_i <- model_i$robust_results
      
      temp <- data.frame(
        covariate_set = set_name,
        outcome = model_i$outcome,
        n = model_i$n,
        term = rownames(robust_i),
        estimate = robust_i[, "Estimate"],
        robust_se = robust_i[, "Std. Error"],
        t_value = robust_i[, "t value"],
        p_value = robust_i[, "Pr(>|t|)"],
        row.names = NULL
      )
      
      output <- rbind(output, temp)
    }
  }
  
  output
}

covariate_results_table <- extract_covariate_results(covariate_model_results)

covariate_results_table

# ------------------------------------------------------------
# Cannabis terms only
# .L = linear cannabis trend
# .Q = quadratic cannabis trend
# ------------------------------------------------------------

cannabis_covariate_terms <- covariate_results_table %>%
  filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  arrange(outcome, covariate_set, term)

cannabis_covariate_terms

# ------------------------------------------------------------
# Key outcomes only
# ------------------------------------------------------------

key_cannabis_covariate_terms <- cannabis_covariate_terms %>%
  filter(outcome %in% c("global_z", "learning_z", "motor_z"))

key_cannabis_covariate_terms

# ------------------------------------------------------------
# Check model N by covariate set and outcome
# ------------------------------------------------------------

model_n_table <- covariate_results_table %>%
  select(covariate_set, outcome, n) %>%
  distinct() %>%
  arrange(outcome, covariate_set)

model_n_table

Outcome ~ DiseaseSev_c + du_mar4_12m_aBin_ord + phq_2_age + phq_7_degree_covfac

#sensitivity for above

Outcome ~ DiseaseSev_c + du_mar4_12m_aBin_ord + phq_2_age + phq_7_degree_covfac + sex_covfac + race_eth_binary_covfac

Outcome ~ DiseaseSev_c + du_mar4_12m_aBin_ord + phq_2_age + phq_7_degree_covfac + sc_mde3

Outcome ~ DiseaseSev_c + du_mar4_12m_aBin_ord + phq_2_age + phq_7_degree_covfac + asr_adh_pro_total


# ============================================================
# FINAL CANDIDACY COVARIATES
# Create new variables; do not overwrite raw variables
# ============================================================

# Age: keep raw age variable, optionally create centered version
data$phq_2_age_c <- as.numeric(
  scale(data$phq_2_age, center = TRUE, scale = FALSE)
)

# Education: use numeric version, NOT factor version
# This assumes phq_7_degree is coded as years/ordered education level.
data$phq_7_degree_num <- as.numeric(data$phq_7_degree)

data$phq_7_degree_c <- as.numeric(
  scale(data$phq_7_degree_num, center = TRUE, scale = FALSE)
)

# Race/ethnicity: use factor version for final model
data$race_eth_binary_covfac <- as.factor(data$race_eth_binary)

# Check distributions
summary(data[, c("phq_2_age", "phq_2_age_c", 
                 "phq_7_degree", "phq_7_degree_num", "phq_7_degree_c")])

table(data$race_eth_binary_covfac, useNA = "ifany")


# ============================================================
# Final candidacy model with selected covariates
# Outcome ~ Disease severity + past-year cannabis + age + education + race/ethnicity
# ============================================================

models_12m_final_covariates <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(
      y,
      " ~ DiseaseSev_c + du_mar4_12m_aBin_ord + ",
      "phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac"
    )
  )
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    outcome = y,
    formula = model_formula,
    n = nobs(model),
    model = model,
    robust_results = robust_results
  )
})

models_12m_final_covariates


# ============================================================
# SUBSTANCE USE COVARIATE CHECK
# Purpose: Check whether alcohol should be included as a covariate
# Rule: Do not overwrite raw variables
# ============================================================

library(psych)
library(dplyr)
library(lmtest)
library(sandwich)

# ------------------------------------------------------------
# 1. Define substance-use variables to check
# ------------------------------------------------------------

substance_covariates_raw <- c(
  "du_alc4_12m_a",   # past-year alcohol
  "du_alc6_30d_a"    # past-30-day alcohol
)

# Check that variables exist
missing_substance_vars <- setdiff(substance_covariates_raw, names(data))
missing_substance_vars

# ------------------------------------------------------------
# 2. Create transformed/centered alcohol variables
#    New variables only; raw variables preserved
# ------------------------------------------------------------

data$du_alc4_12m_a_log1p <- log1p(data$du_alc4_12m_a)
data$du_alc6_30d_a_log1p <- log1p(data$du_alc6_30d_a)

data$du_alc4_12m_a_log1p_c <- as.numeric(
  scale(data$du_alc4_12m_a_log1p, center = TRUE, scale = FALSE)
)

data$du_alc6_30d_a_log1p_c <- as.numeric(
  scale(data$du_alc6_30d_a_log1p, center = TRUE, scale = FALSE)
)

summary(data[, c(
  "du_alc4_12m_a",
  "du_alc4_12m_a_log1p",
  "du_alc6_30d_a",
  "du_alc6_30d_a_log1p"
)])

# ------------------------------------------------------------
# 3. Correlations: alcohol with cannabis and neurocognition
# ------------------------------------------------------------

alcohol_check_vars <- c(
  "du_alc4_12m_a_log1p",
  "du_alc6_30d_a_log1p"
)

primary_cannabis_vars <- c(
  "du_mar4_12m_a",
  "du_mar6_30d_a",
  "du_mar2_life_a"
)

outcome_vars <- c(
  "global_z",
  "learning_z",
  "memory_z",
  "motor_z",
  "executive_z",
  "processing_speed_z"
)

substance_screen_data <- data[, c(
  alcohol_check_vars,
  primary_cannabis_vars,
  outcome_vars
)]

# Alcohol x cannabis
alcohol_cannabis_corr <- psych::corr.test(
  substance_screen_data[, c(alcohol_check_vars, primary_cannabis_vars)],
  method = "spearman",
  use = "pairwise.complete.obs",
  adjust = "none"
)

alcohol_cannabis_r <- alcohol_cannabis_corr$r[
  alcohol_check_vars,
  primary_cannabis_vars
]

alcohol_cannabis_p <- alcohol_cannabis_corr$p[
  alcohol_check_vars,
  primary_cannabis_vars
]

alcohol_cannabis_r
alcohol_cannabis_p


# Alcohol x cognition
alcohol_neurocog_corr <- psych::corr.test(
  substance_screen_data[, c(alcohol_check_vars, outcome_vars)],
  method = "spearman",
  use = "pairwise.complete.obs",
  adjust = "none"
)

alcohol_neurocog_r <- alcohol_neurocog_corr$r[
  alcohol_check_vars,
  outcome_vars
]

alcohol_neurocog_p <- alcohol_neurocog_corr$p[
  alcohol_check_vars,
  outcome_vars
]

alcohol_neurocog_r
alcohol_neurocog_p

# ============================================================
# 4. Model comparison: primary past-year cannabis model
#    without vs. with alcohol covariates
# ============================================================

# Make sure final covariates exist
data$phq_2_age_c <- as.numeric(
  scale(data$phq_2_age, center = TRUE, scale = FALSE)
)

data$phq_7_degree_num <- as.numeric(data$phq_7_degree)

data$phq_7_degree_c <- as.numeric(
  scale(data$phq_7_degree_num, center = TRUE, scale = FALSE)
)

data$race_eth_binary_covfac <- as.factor(data$race_eth_binary)

# ------------------------------------------------------------
# Function to run robust model
# ------------------------------------------------------------

run_robust_lm <- function(model_formula, data) {
  
  model <- lm(model_formula, data = data)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  list(
    formula = model_formula,
    n = nobs(model),
    model = model,
    robust_results = robust_results
  )
}

models_no_alcohol <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(
      y,
      " ~ DiseaseSev_c + du_mar4_12m_aBin_ord + ",
      "phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac"
    )
  )
  
  result <- run_robust_lm(model_formula, data)
  result$outcome <- y
  result
})

models_no_alcohol


models_with_12m_alcohol <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(
      y,
      " ~ DiseaseSev_c + du_mar4_12m_aBin_ord + ",
      "phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac + ",
      "du_alc4_12m_a_log1p_c"
    )
  )
  
  result <- run_robust_lm(model_formula, data)
  result$outcome <- y
  result
})

models_with_12m_alcohol

models_with_30d_alcohol <- lapply(outcome_vars, function(y) {
  
  model_formula <- as.formula(
    paste0(
      y,
      " ~ DiseaseSev_c + du_mar4_12m_aBin_ord + ",
      "phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac + ",
      "du_alc6_30d_a_log1p_c"
    )
  )
  
  result <- run_robust_lm(model_formula, data)
  result$outcome <- y
  result
})

models_with_30d_alcohol

# ============================================================
# 5. Extract cannabis terms across alcohol-adjusted models
# ============================================================

extract_robust_results <- function(model_list, model_name) {
  
  output <- data.frame()
  
  for (i in seq_along(model_list)) {
    
    robust_i <- model_list[[i]]$robust_results
    
    temp <- data.frame(
      model_set = model_name,
      outcome = model_list[[i]]$outcome,
      n = model_list[[i]]$n,
      term = rownames(robust_i),
      estimate = robust_i[, "Estimate"],
      robust_se = robust_i[, "Std. Error"],
      t_value = robust_i[, "t value"],
      p_value = robust_i[, "Pr(>|t|)"],
      row.names = NULL
    )
    
    output <- rbind(output, temp)
  }
  
  output
}

table_no_alcohol <- extract_robust_results(
  models_no_alcohol,
  "No alcohol covariate"
)

table_12m_alcohol <- extract_robust_results(
  models_with_12m_alcohol,
  "Adjusted for past-year alcohol"
)

table_30d_alcohol <- extract_robust_results(
  models_with_30d_alcohol,
  "Adjusted for past-30-day alcohol"
)

alcohol_sensitivity_table <- rbind(
  table_no_alcohol,
  table_12m_alcohol,
  table_30d_alcohol
)

# Pull only cannabis terms
alcohol_sensitivity_cannabis_terms <- alcohol_sensitivity_table %>%
  filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  arrange(outcome, model_set, term)

alcohol_sensitivity_cannabis_terms

# If .Q stays significant/similar:
# The past-year cannabis quadratic effect is robust to alcohol adjustment.
# If .Q weakens substantially:
# Alcohol use may partially account for the cannabis-cognition association.
# If alcohol itself is nonsignificant and .Q is unchanged:
# Alcohol was tested but did not materially affect the primary cannabis finding.


# ============================================================
# Plot significant past-year cannabis quadratic effects
# Outcomes: global_z, learning_z, motor_z
# Predictor: du_mar4_12m_aBin_ord
# Adjusted for: DiseaseSev_c, age, education, race/ethnicity
# ============================================================

# Install packages only if needed
if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(emmeans)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# Confirm ordered cannabis variable
# ------------------------------------------------------------

table(data$du_mar4_12m_aBin_ord, useNA = "ifany")
contrasts(data$du_mar4_12m_aBin_ord)

# ------------------------------------------------------------
# Make sure covariates exist
# ------------------------------------------------------------

data$phq_2_age_c <- as.numeric(
  scale(data$phq_2_age, center = TRUE, scale = FALSE)
)

data$phq_7_degree_num <- as.numeric(data$phq_7_degree)

data$phq_7_degree_c <- as.numeric(
  scale(data$phq_7_degree_num, center = TRUE, scale = FALSE)
)

data$race_eth_binary_covfac <- as.factor(data$race_eth_binary)

# Check
summary(data[, c("phq_2_age_c", "phq_7_degree_c", "DiseaseSev_c")])
table(data$race_eth_binary_covfac, useNA = "ifany")

# ============================================================
# Final adjusted models for significant quadratic outcomes
# ============================================================

model_global_quad <- lm(
  global_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_learning_quad <- lm(
  learning_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_motor_quad <- lm(
  motor_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

# ============================================================
# Estimated marginal means by cannabis group
# ============================================================

emm_global <- as.data.frame(
  emmeans(model_global_quad, ~ du_mar4_12m_aBin_ord)
)

emm_learning <- as.data.frame(
  emmeans(model_learning_quad, ~ du_mar4_12m_aBin_ord)
)

emm_motor <- as.data.frame(
  emmeans(model_motor_quad, ~ du_mar4_12m_aBin_ord)
)

# Add outcome labels
emm_global$outcome <- "Global cognition"
emm_learning$outcome <- "Learning"
emm_motor$outcome <- "Motor"

# Combine
emm_quad_plot_data <- rbind(
  emm_global,
  emm_learning,
  emm_motor
)

emm_quad_plot_data

# ============================================================
# Plot adjusted estimated marginal means
# ============================================================

ggplot(
  emm_quad_plot_data,
  aes(
    x = du_mar4_12m_aBin_ord,
    y = emmean,
    group = outcome
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.10
  ) +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(
    title = "Adjusted Neurocognitive Performance by Past-Year Cannabis Exposure",
    x = "Past-year cannabis exposure",
    y = "Adjusted estimated marginal mean z-score"
  ) +
  theme_classic(base_size = 14)


# ============================================================
# Raw data + adjusted estimated marginal means
# ============================================================

raw_quad_plot_data <- data %>%
  select(
    du_mar4_12m_aBin_ord,
    global_z,
    learning_z,
    motor_z
  ) %>%
  tidyr::pivot_longer(
    cols = c(global_z, learning_z, motor_z),
    names_to = "outcome",
    values_to = "z_score"
  ) %>%
  mutate(
    outcome = case_when(
      outcome == "global_z" ~ "Global cognition",
      outcome == "learning_z" ~ "Learning",
      outcome == "motor_z" ~ "Motor",
      TRUE ~ outcome
    )
  )

ggplot() +
  geom_jitter(
    data = raw_quad_plot_data,
    aes(
      x = du_mar4_12m_aBin_ord,
      y = z_score
    ),
    width = 0.08,
    alpha = 0.35,
    size = 1.8
  ) +
  geom_line(
    data = emm_quad_plot_data,
    aes(
      x = du_mar4_12m_aBin_ord,
      y = emmean,
      group = outcome
    ),
    linewidth = 1
  ) +
  geom_point(
    data = emm_quad_plot_data,
    aes(
      x = du_mar4_12m_aBin_ord,
      y = emmean
    ),
    size = 3
  ) +
  geom_errorbar(
    data = emm_quad_plot_data,
    aes(
      x = du_mar4_12m_aBin_ord,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.10
  ) +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(
    title = "Past-Year Cannabis Exposure and Neurocognitive Performance",
    subtitle = "Points show observed data; lines show covariate-adjusted estimated marginal means",
    x = "Past-year cannabis exposure",
    y = "Neurocognitive z-score"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = "past_year_cannabis_quadratic_effects.png",
  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================
# Pairwise group comparisons for significant quadratic outcomes
# Outcomes: global_z, learning_z, motor_z
# Predictor: du_mar4_12m_aBin_ord
# Groups: none, low, high
# Adjusted for: DiseaseSev_c, age, education, race/ethnicity
# ============================================================

if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
if (!requireNamespace("lmtest", quietly = TRUE)) install.packages("lmtest")
if (!requireNamespace("sandwich", quietly = TRUE)) install.packages("sandwich")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(emmeans)
library(lmtest)
library(sandwich)
library(dplyr)

# Make sure covariates exist
data$phq_2_age_c <- as.numeric(
  scale(data$phq_2_age, center = TRUE, scale = FALSE)
)

data$phq_7_degree_num <- as.numeric(data$phq_7_degree)

data$phq_7_degree_c <- as.numeric(
  scale(data$phq_7_degree_num, center = TRUE, scale = FALSE)
)

data$race_eth_binary_covfac <- as.factor(data$race_eth_binary)

# Fit final models
model_global_quad <- lm(
  global_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_learning_quad <- lm(
  learning_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_motor_quad <- lm(
  motor_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

# Make sure covariates exist
data$phq_2_age_c <- as.numeric(
  scale(data$phq_2_age, center = TRUE, scale = FALSE)
)

data$phq_7_degree_num <- as.numeric(data$phq_7_degree)

data$phq_7_degree_c <- as.numeric(
  scale(data$phq_7_degree_num, center = TRUE, scale = FALSE)
)

data$race_eth_binary_covfac <- as.factor(data$race_eth_binary)

# Fit final models
model_global_quad <- lm(
  global_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_learning_quad <- lm(
  learning_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_motor_quad <- lm(
  motor_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

#test marginal means by cannabis group

emm_global <- emmeans(model_global_quad, ~ du_mar4_12m_aBin_ord)
emm_learning <- emmeans(model_learning_quad, ~ du_mar4_12m_aBin_ord)
emm_motor <- emmeans(model_motor_quad, ~ du_mar4_12m_aBin_ord)

emm_global
emm_learning
emm_motor

#pairwise comparisons w/ tukey adjustment 

# Global cognition
pairs_global_tukey <- pairs(
  emm_global,
  adjust = "tukey"
)

pairs_global_tukey


# Learning
pairs_learning_tukey <- pairs(
  emm_learning,
  adjust = "tukey"
)

pairs_learning_tukey


# Motor
pairs_motor_tukey <- pairs(
  emm_motor,
  adjust = "tukey"
)

pairs_motor_tukey

# Global cognition with HC3 robust SEs
emm_global_robust <- emmeans(
  model_global_quad,
  ~ du_mar4_12m_aBin_ord,
  vcov. = sandwich::vcovHC(model_global_quad, type = "HC3")
)

pairs_global_robust <- pairs(
  emm_global_robust,
  adjust = "tukey"
)

pairs_global_robust


# Learning with HC3 robust SEs
emm_learning_robust <- emmeans(
  model_learning_quad,
  ~ du_mar4_12m_aBin_ord,
  vcov. = sandwich::vcovHC(model_learning_quad, type = "HC3")
)

pairs_learning_robust <- pairs(
  emm_learning_robust,
  adjust = "tukey"
)

pairs_learning_robust


# Motor with HC3 robust SEs
emm_motor_robust <- emmeans(
  model_motor_quad,
  ~ du_mar4_12m_aBin_ord,
  vcov. = sandwich::vcovHC(model_motor_quad, type = "HC3")
)

pairs_motor_robust <- pairs(
  emm_motor_robust,
  adjust = "tukey"
)

pairs_motor_robust

# ============================================================
# Combine robust pairwise comparisons into one table
# ============================================================

global_pairwise_table <- as.data.frame(pairs_global_robust)
global_pairwise_table$outcome <- "Global cognition"

learning_pairwise_table <- as.data.frame(pairs_learning_robust)
learning_pairwise_table$outcome <- "Learning"

motor_pairwise_table <- as.data.frame(pairs_motor_robust)
motor_pairwise_table$outcome <- "Motor"

pairwise_cannabis_group_table <- rbind(
  global_pairwise_table,
  learning_pairwise_table,
  motor_pairwise_table
)

pairwise_cannabis_group_table <- pairwise_cannabis_group_table %>%
  select(
    outcome,
    contrast,
    estimate,
    SE,
    df,
    t.ratio,
    p.value
  )

pairwise_cannabis_group_table

# ============================================================
# Figure: Significant past-year cannabis quadratic effects
# with follow-up pairwise significance annotations
# ============================================================

if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1. Make sure variables are ready
# ------------------------------------------------------------

data$phq_2_age_c <- as.numeric(
  scale(data$phq_2_age, center = TRUE, scale = FALSE)
)

data$phq_7_degree_num <- as.numeric(data$phq_7_degree)

data$phq_7_degree_c <- as.numeric(
  scale(data$phq_7_degree_num, center = TRUE, scale = FALSE)
)

data$race_eth_binary_covfac <- as.factor(data$race_eth_binary)

# Make sure cannabis group is ordered correctly
data$du_mar4_12m_aBin_ord <- ordered(
  data$du_mar4_12m_aBin,
  levels = c(0, 1, 2),
  labels = c("No use", "Low use", "High use")
)

contrasts(data$du_mar4_12m_aBin_ord) <- contr.poly(3)

table(data$du_mar4_12m_aBin_ord, useNA = "ifany")

# ------------------------------------------------------------
# 2. Fit adjusted models for significant quadratic outcomes
# ------------------------------------------------------------

model_global_quad <- lm(
  global_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_learning_quad <- lm(
  learning_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_motor_quad <- lm(
  motor_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

# ------------------------------------------------------------
# 3. Estimated marginal means
# ------------------------------------------------------------

emm_global <- as.data.frame(
  emmeans(model_global_quad, ~ du_mar4_12m_aBin_ord)
)

emm_learning <- as.data.frame(
  emmeans(model_learning_quad, ~ du_mar4_12m_aBin_ord)
)

emm_motor <- as.data.frame(
  emmeans(model_motor_quad, ~ du_mar4_12m_aBin_ord)
)

emm_global$outcome <- "Global cognition"
emm_learning$outcome <- "Learning"
emm_motor$outcome <- "Motor"

emm_quad_plot_data <- rbind(
  emm_global,
  emm_learning,
  emm_motor
)

# Rename x variable for cleaner plotting
emm_quad_plot_data <- emm_quad_plot_data %>%
  rename(cannabis_group = du_mar4_12m_aBin_ord) %>%
  mutate(
    cannabis_group = factor(
      cannabis_group,
      levels = c("No use", "Low use", "High use")
    ),
    outcome = factor(
      outcome,
      levels = c("Global cognition", "Learning", "Motor")
    )
  )

emm_quad_plot_data

# ------------------------------------------------------------
# 4. Create significance annotation data
# Brackets show Low use vs High use follow-up comparisons
# ------------------------------------------------------------

sig_annotations <- data.frame(
  outcome = factor(
    c("Global cognition", "Learning", "Motor"),
    levels = c("Global cognition", "Learning", "Motor")
  ),
  x_start = c(2, 2, 2),  # Low use
  x_end = c(3, 3, 3),    # High use
  label = c("Low > High\np = .032", "Low > High\np = .015", "Low > High\np = .089"),
  stringsAsFactors = FALSE
)

# Place brackets slightly above each outcome's CI range
y_positions <- emm_quad_plot_data %>%
  group_by(outcome) %>%
  summarise(
    y_position = max(upper.CL, na.rm = TRUE) + 0.12,
    .groups = "drop"
  )

sig_annotations <- sig_annotations %>%
  left_join(y_positions, by = "outcome")

sig_annotations

# ------------------------------------------------------------
# 5. Plot adjusted means with significance annotations
# ------------------------------------------------------------

quad_sig_plot <- ggplot(
  emm_quad_plot_data,
  aes(
    x = cannabis_group,
    y = emmean,
    group = 1
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.10,
    linewidth = 0.7
  ) +
  geom_segment(
    data = sig_annotations,
    aes(
      x = x_start,
      xend = x_end,
      y = y_position,
      yend = y_position
    ),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = sig_annotations,
    aes(
      x = x_start,
      xend = x_start,
      y = y_position,
      yend = y_position - 0.04
    ),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = sig_annotations,
    aes(
      x = x_end,
      xend = x_end,
      y = y_position,
      yend = y_position - 0.04
    ),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_text(
    data = sig_annotations,
    aes(
      x = (x_start + x_end) / 2,
      y = y_position + 0.05,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    lineheight = 0.9
  ) +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(
    title = "Adjusted Neurocognitive Performance by Past-Year Cannabis Exposure",
    subtitle = "Brackets show Tukey-adjusted follow-up comparisons for Low vs High use",
    x = "Past-year cannabis exposure",
    y = "Adjusted estimated marginal mean z-score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

quad_sig_plot

ggsave(
  filename = "past_year_cannabis_quadratic_effects_with_significance.png",
  plot = quad_sig_plot,
  width = 11,
  height = 6,
  dpi = 300
)

sig_annotations <- data.frame(
  outcome = factor(
    c("Global cognition", "Learning", "Motor"),
    levels = c("Global cognition", "Learning", "Motor")
  ),
  x_start = c(2, 2, 2),
  x_end = c(3, 3, 3),
  label = c("*", "*", "†"),
  stringsAsFactors = FALSE
)

y_positions <- emm_quad_plot_data %>%
  group_by(outcome) %>%
  summarise(
    y_position = max(upper.CL, na.rm = TRUE) + 0.12,
    .groups = "drop"
  )

sig_annotations <- sig_annotations %>%
  left_join(y_positions, by = "outcome")

# ============================================================
# EFFECT SIZE EXTRACTION FOR SIGNIFICANT QUADRATIC EFFECTS
# ============================================================

if (!requireNamespace("effectsize", quietly = TRUE)) install.packages("effectsize")
if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("lmtest", quietly = TRUE)) install.packages("lmtest")
if (!requireNamespace("sandwich", quietly = TRUE)) install.packages("sandwich")

library(effectsize)
library(emmeans)
library(dplyr)
library(lmtest)
library(sandwich)

# ============================================================
# Final adjusted models for significant quadratic outcomes
# ============================================================

model_global_quad <- lm(
  global_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_learning_quad <- lm(
  learning_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

model_motor_quad <- lm(
  motor_z ~ DiseaseSev_c + du_mar4_12m_aBin_ord +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = data
)

# ============================================================
# Partial eta-squared for model terms
# ============================================================

eta_global <- effectsize::eta_squared(
  model_global_quad,
  partial = TRUE
)

eta_learning <- effectsize::eta_squared(
  model_learning_quad,
  partial = TRUE
)

eta_motor <- effectsize::eta_squared(
  model_motor_quad,
  partial = TRUE
)

eta_global
eta_learning
eta_motor

eta_global$outcome <- "Global cognition"
eta_learning$outcome <- "Learning"
eta_motor$outcome <- "Motor"

eta_quad_table <- rbind(
  as.data.frame(eta_global),
  as.data.frame(eta_learning),
  as.data.frame(eta_motor)
)

eta_quad_table

# ------------------------------------------------------------
# Pull cannabis-related eta-squared rows
# ------------------------------------------------------------

cannabis_eta_table <- eta_quad_table %>%
  filter(grepl("du_mar4_12m_aBin_ord", Parameter))

cannabis_eta_table

# ============================================================
# Specific effect size for quadratic cannabis term using robust HC3 results
# ============================================================

extract_quad_effect_size <- function(model, outcome_name) {
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  term_name <- "du_mar4_12m_aBin_ord.Q"
  
  quad_row <- robust_results[term_name, ]
  
  t_value <- quad_row["t value"]
  p_value <- quad_row["Pr(>|t|)"]
  estimate <- quad_row["Estimate"]
  se <- quad_row["Std. Error"]
  
  df_resid <- df.residual(model)
  
  partial_eta_sq <- (t_value^2) / ((t_value^2) + df_resid)
  partial_r <- sqrt(partial_eta_sq)
  
  data.frame(
    outcome = outcome_name,
    term = term_name,
    estimate = estimate,
    robust_se = se,
    t_value = t_value,
    df = df_resid,
    p_value = p_value,
    partial_eta_sq = partial_eta_sq,
    partial_r = partial_r,
    row.names = NULL
  )
}

quad_effect_sizes <- rbind(
  extract_quad_effect_size(model_global_quad, "Global cognition"),
  extract_quad_effect_size(model_learning_quad, "Learning"),
  extract_quad_effect_size(model_motor_quad, "Motor")
)

quad_effect_sizes

# ============================================================
# Adjusted Cohen's d for pairwise Low vs High comparisons
# ============================================================

emm_global <- emmeans(
  model_global_quad,
  ~ du_mar4_12m_aBin_ord
)

emm_learning <- emmeans(
  model_learning_quad,
  ~ du_mar4_12m_aBin_ord
)

emm_motor <- emmeans(
  model_motor_quad,
  ~ du_mar4_12m_aBin_ord
)

# Cohen's d based on model residual SD
d_global <- eff_size(
  emm_global,
  sigma = sigma(model_global_quad),
  edf = df.residual(model_global_quad)
)

d_learning <- eff_size(
  emm_learning,
  sigma = sigma(model_learning_quad),
  edf = df.residual(model_learning_quad)
)

d_motor <- eff_size(
  emm_motor,
  sigma = sigma(model_motor_quad),
  edf = df.residual(model_motor_quad)
)

d_global
d_learning
d_motor

# ------------------------------------------------------------
# Extract Low vs High Cohen's d
# ------------------------------------------------------------

d_global_table <- as.data.frame(d_global)
d_learning_table <- as.data.frame(d_learning)
d_motor_table <- as.data.frame(d_motor)

d_global_table$outcome <- "Global cognition"
d_learning_table$outcome <- "Learning"
d_motor_table$outcome <- "Motor"

pairwise_d_table <- rbind(
  d_global_table,
  d_learning_table,
  d_motor_table
)

# Look at all pairwise effect sizes first
pairwise_d_table

# Pull low vs high contrast
low_high_d_table <- pairwise_d_table %>%
  filter(grepl("Low use - High use|low - high|Low - High", contrast, ignore.case = TRUE))

low_high_d_table

unique(pairwise_d_table$contrast)

# ============================================================
# Add conventional interpretation labels
# Cohen's d: .20 small, .50 medium, .80 large
# partial eta-squared: .01 small, .06 medium, .14 large
# ============================================================

quad_effect_sizes <- quad_effect_sizes %>%
  mutate(
    eta_interpretation = case_when(
      partial_eta_sq < .01 ~ "very small",
      partial_eta_sq < .06 ~ "small",
      partial_eta_sq < .14 ~ "medium",
      partial_eta_sq >= .14 ~ "large",
      TRUE ~ NA_character_
    )
  )

quad_effect_sizes

low_high_d_table <- low_high_d_table %>%
  mutate(
    d_interpretation = case_when(
      abs(effect.size) < .20 ~ "very small",
      abs(effect.size) < .50 ~ "small",
      abs(effect.size) < .80 ~ "medium",
      abs(effect.size) >= .80 ~ "large",
      TRUE ~ NA_character_
    )
  )

low_high_d_table

##################

# ============================================================
# Effect sizes for past-year cannabis x HIV disease severity models
# DiseaseSev_c is a focal predictor/moderator, not just a covariate
# Demographic covariates: age, education, race/ethnicity
# ============================================================

library(lmtest)
library(sandwich)
library(effectsize)
library(performance)
library(parameters)
library(emmeans)
library(car)
library(dplyr)

# Outcomes with significant cannabis quadratic effects
sig_outcomes <- c("global_z", "learning_z", "motor_z")

# Demographic covariates only
demographic_covariates <- c(
  "phq_2_age",
  "phq_7_degree",
  "race_eth_binary"
)

# Check variables exist
all_needed_vars <- c(
  sig_outcomes,
  "du_mar4_12m_aBin_ord",
  "DiseaseSev_c",
  demographic_covariates
)

setdiff(all_needed_vars, names(data))

# ============================================================
# Run reduced and interaction models
# ============================================================

effect_size_results <- lapply(sig_outcomes, function(y) {
  
  # Reduced/main-effects model:
  # cannabis + DiseaseSev + demographic covariates
  formula_main <- as.formula(
    paste0(
      y,
      " ~ du_mar4_12m_aBin_ord + DiseaseSev_c + ",
      paste(demographic_covariates, collapse = " + ")
    )
  )
  
  # Full moderation model:
  # cannabis x DiseaseSev + demographic covariates
  formula_interaction <- as.formula(
    paste0(
      y,
      " ~ du_mar4_12m_aBin_ord * DiseaseSev_c + ",
      paste(demographic_covariates, collapse = " + ")
    )
  )
  
  model_main <- lm(formula_main, data = data)
  model_interaction <- lm(formula_interaction, data = data)
  
  # Robust coefficient tables
  robust_main <- lmtest::coeftest(
    model_main,
    vcov. = sandwich::vcovHC(model_main, type = "HC3")
  )
  
  robust_interaction <- lmtest::coeftest(
    model_interaction,
    vcov. = sandwich::vcovHC(model_interaction, type = "HC3")
  )
  
  # Robust CIs
  ci_main <- parameters::model_parameters(
    model_main,
    vcov = sandwich::vcovHC(model_main, type = "HC3"),
    ci = 0.95,
    standardize = NULL
  )
  
  ci_interaction <- parameters::model_parameters(
    model_interaction,
    vcov = sandwich::vcovHC(model_interaction, type = "HC3"),
    ci = 0.95,
    standardize = NULL
  )
  
  # R2 for both models
  r2_main <- performance::r2(model_main)
  r2_interaction <- performance::r2(model_interaction)
  
  # Change in R2 from adding cannabis x DiseaseSev interaction
  r2_change <- r2_interaction$R2 - r2_main$R2
  
  # Omnibus comparison: does interaction improve model fit?
  model_compare <- anova(model_main, model_interaction)
  
  # Type III tests
  anova_main <- car::Anova(model_main, type = 3)
  anova_interaction <- car::Anova(model_interaction, type = 3)
  
  eta_main <- effectsize::eta_squared(
    anova_main,
    partial = TRUE,
    ci = 0.95
  )
  
  eta_interaction <- effectsize::eta_squared(
    anova_interaction,
    partial = TRUE,
    ci = 0.95
  )
  
  # Estimated marginal means from reduced model
  # Use this only if interaction is nonsignificant and you are interpreting cannabis main/quadratic pattern
  emm_main <- emmeans::emmeans(
    model_main,
    ~ du_mar4_12m_aBin_ord
  )
  
  pairwise_main <- pairs(
    emm_main,
    adjust = "tukey"
  )
  
  cohens_d_main <- emmeans::eff_size(
    emm_main,
    sigma = sigma(model_main),
    edf = df.residual(model_main)
  )
  
  list(
    outcome = y,
    
    formula_main = formula_main,
    formula_interaction = formula_interaction,
    
    model_main = model_main,
    model_interaction = model_interaction,
    
    robust_main = robust_main,
    robust_interaction = robust_interaction,
    
    ci_main = ci_main,
    ci_interaction = ci_interaction,
    
    r2_main = r2_main,
    r2_interaction = r2_interaction,
    r2_change_interaction = r2_change,
    model_compare_interaction = model_compare,
    
    partial_eta_main = eta_main,
    partial_eta_interaction = eta_interaction,
    
    estimated_marginal_means_main = emm_main,
    pairwise_main = pairwise_main,
    cohens_d_main = cohens_d_main
  )
})

effect_size_results

# Global
effect_size_results[[1]]$formula_interaction
effect_size_results[[1]]$robust_interaction
effect_size_results[[1]]$model_compare_interaction
effect_size_results[[1]]$r2_main
effect_size_results[[1]]$r2_interaction
effect_size_results[[1]]$r2_change_interaction
effect_size_results[[1]]$partial_eta_interaction

# Learning
effect_size_results[[2]]$formula_interaction
effect_size_results[[2]]$robust_interaction
effect_size_results[[2]]$model_compare_interaction
effect_size_results[[2]]$r2_main
effect_size_results[[2]]$r2_interaction
effect_size_results[[2]]$r2_change_interaction
effect_size_results[[2]]$partial_eta_interaction

# Motor
effect_size_results[[3]]$formula_interaction
effect_size_results[[3]]$robust_interaction
effect_size_results[[3]]$model_compare_interaction
effect_size_results[[3]]$r2_main
effect_size_results[[3]]$r2_interaction
effect_size_results[[3]]$r2_change_interaction
effect_size_results[[3]]$partial_eta_interaction
