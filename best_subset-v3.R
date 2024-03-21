# Setup to perform a best subset analysis for a linear mixed model
#
# v3 Uses a simplified model w/o the interactions motivated by literature
#  residual analysis during best subset selection
#
# Author: James Henderson
# Updated: December 19, 2023
# 79: -------------------------------------------------------------------------

# libraries: ------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(lme4)

# paths: ----------------------------------------------------------------------
git_path = "~/github/chlorine_review/"

# data: -----------------------------------------------------------------------
foo = load(sprintf("%s/12_14_23_df_encoded.RData", git_path))
df = df_encoded

# formula with all candidate variables: ---------------------------------------
fix0 = log_average_kobs ~
    temp_5int +
    pH_refpka +
    float_year_scaled +
    alpha_0 +
    pH_refpka_x_high_chloride_TRUE +
    pH_refpka_x_temp_5int +
    # group buffers
    buffer_high_organics +
    buffer_natural_or_treated_water +
    buffer_ultrapure +
    # group purification levels
    purification_3 +
    purification_2 +
    high_chloride_TRUE
    
# best subset selection: ------------------------------------------------------
# Approach:
#  1. Create a list with each component a set of columns include/exclude
#  2. Generate all subsets of the indices of that list
#  3. For each subset, create the formula, fit the model, and save the scores.

# this is just a way to generate the list from the formula, you can add
#  items in other ways
iv_list = as.character(fix0)[3] %>%
  str_replace_all("\\n", "") %>%
  str_replace_all(" ", "") %>%
  str_split(., "[+]")
iv_list = as.list(iv_list[[1]])

# group buffers
b_idx = vapply(iv_list, \(x) grepl("^buffer", x), TRUE) |> which()
iv_list[[b_idx[1]]] = unlist(iv_list[b_idx])
iv_list[b_idx[-1]] = NULL

# group purification levels
p_idx = vapply(iv_list, \(x) any(grepl("^pur", x)), TRUE) |> which()
iv_list[[p_idx[1]]] = unlist(iv_list[p_idx])
iv_list[p_idx[-1]] = NULL

# this turns any multiple groups into 1 entry separated by " + "
iv_vec = lapply(iv_list, paste0, collapse = " + ") |> unlist()

n_iv = length(iv_vec)
#iv_subsets = lapply(seq_len(n_iv), combn, x = n_iv)

# w/o grouping largest subsets has ~1e7 and indices take up 1GB of mem
# with grouping down to choose(22, 11) ~1.4e6

# "saturated" model with all terms: -------------------------------------------
f_template =
  "log_average_kobs ~ virus_name_abbrev + %s + (1|paper_ID)"
f_iv = paste(iv_vec, collapse = " + ")
f_sat = sprintf(f_template, f_iv) |> as.formula()
mod_sat = lmer(f_sat, data = df, REML = FALSE)
aic_bic_sat = c(AIC(mod_sat), BIC(mod_sat))

# .11 s
system.time({mod_sat = lmer(f_sat, data = df, REML = FALSE)})

# function to fit a subset: ---------------------------------------------------
fit_subset = function(j, f_temp = f_template) {
  f_iv = iv_vec[iv_subset[, j]] |> paste(collapse = " + ")
  f = sprintf(f_temp, f_iv) |> as.formula()
  mod = update(mod_sat, formula = f)
  c(AIC(mod), BIC(mod))
}

# run through all subsets: ----------------------------------------------------
aic_bic_list = vector("list", length = 1 + n_iv)

for ( k in 0:n_iv ) {

  iv_subset = combn(n_iv, k)

  c# This is where all the computation happens
  aic_bic = vapply(seq_len(ncol(iv_subset)), fit_subset, c(0.0, 0.0))
  aic_bic = t(aic_bic)
  colnames(aic_bic) = c("aic", "bic")

  # format result as a data.table
  aic_bic_dt = data.table(
    row = row,
    k = k,
    i = seq_len(ncol(iv_subset)),
    aic = aic_bic[, "aic"],
    bic = aic_bic[, "bic"]
  )

  aic_bic_list[[k + 1]] = aic_bic_dt
}

# all results: ----------------------------------------------------------------
res = rbindlist(aic_bic_list)

# best model by AIC: ----------------------------------------------------------
aic_idx = res[, which.min(aic)]
best_aic = res[aic_idx, aic]
best_k = res[aic_idx, k]
best_i = res[aic_idx, i]

iv_subset = combn(n_iv, best_k)
iv_idx = iv_subset[, best_i]
iv_aic = paste0(iv_vec[iv_idx], collapse = " + ")
f_best_aic = sprintf(f_template, iv_aic) |> as.formula()
m_best_aic = update(mod_sat, f_best_aic)
stopifnot(near(best_aic, AIC(m_best_aic)))

# best model by BIC: ----------------------------------------------------------
bic_idx = res[, which.min(bic)]
best_bic = res[bic_idx, bic]
best_k = res[bic_idx, k]
best_i = res[bic_idx, i]

iv_subset = combn(n_iv, best_k)
iv_idx = iv_subset[, best_i]
iv_bic = paste0(iv_vec[iv_idx], collapse = " + ")
f_best_bic = sprintf(f_template, iv_bic) |> as.formula()
m_best_bic = update(mod_sat, f_best_bic)
stopifnot(near(best_bic, BIC(m_best_bic)))

save(
  res, best_aic, f_best_aic, m_best_aic, best_bic, f_best_bic, m_best_bic,
  file = sprintf(
    "%s/best_subset/results/best_subset_result_v3.RData",
    git_path
 )
)