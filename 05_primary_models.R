# Nonlinear Associations Between Cannabis Use and Neurocognition Among People Living With HIV
# Script: 05_primary_models.R

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
