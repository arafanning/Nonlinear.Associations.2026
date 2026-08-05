# Nonlinear Associations Between Cannabis Use and Neurocognition Among People Living With HIV
# Script: 03_disease_severity_factor.R
# Purpose: 03 Disease Severity Factor
#

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

