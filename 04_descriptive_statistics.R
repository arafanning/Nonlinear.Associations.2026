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
