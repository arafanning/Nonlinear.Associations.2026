# Nonlinear Associations Between Cannabis Use and Neurocognition Among People Living With HIV
# Script: 04_descriptive_statistics.R


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
