# Nonlinear Associations Between Cannabis Use and Neurocognition Among People Living With HIV
# Script: 01_setup.R
# Purpose: 01 Setup


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
