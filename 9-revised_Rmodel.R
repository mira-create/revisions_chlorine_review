library(readxl)
library(writexl)
library(lme4)
library(sjPlot)
library(jtools) 
library(lmerTest)
library(performance)
library(MuMIn)
library(dplyr)
library(cowplot)
library(gridExtra)
library(patchwork)
library(tibble)
library(lmerTest)
library(lattice)
library(report)
library(car)
library(ggplot2)


# Read the dataset ----------------------------------------------
df <- read_excel('/Users/mirac/Desktop/Chlorine_Project/August 18 Chlorine Review Paper/Chlorine_Review_Revisions/8-Revisions_merged_dataset_postprocessing.xlsx')

setwd("/Users/mirac/github/chlorine_review")

#View(df)

# Check data type ----------------------------------------------
class(df$purification_level) # "numeric"
class(df$year)               # "character"
class(df$year_float)         # "numeric"
class(df$year_int)           # "numeric"
class(df$float_year_scaled)  # "numeric"
class(df$int_year_scaled)    # "numeric"
class(df$buffer)             # "character"
class(df$genome_length)     # "character"
class(df$diameter)            #numeric
class(df$temp)              #numeric
class(df$pH)
class(df$high_chloride_new)             # 

# Change data type ----------------------------------------------
df$year = as.numeric(df$year)
df$genome_length = as.numeric(df$genome_length)
df$purification_level = as.character(df$purification_level)

# Double check
class(df$purification_level) # "character"
class(df$year)               # "numeric"

#Relevel columns ----------------------------------------------
df$buffer = relevel(factor(df$buffer), "synthetic_buffer")
df$high_chloride_new = relevel(factor(df$high_chloride_new), "FALSE")
df$balt_class = relevel(factor(df$balt_class), "+ssRNA")
df$family = relevel(factor(df$family), "Fiersviridae")
df$genus = relevel(factor(df$genus), "Emesvirus")
df$species = relevel(factor(df$species), "Emesvirus zinderi")
df$virus_name_abbrev = relevel(factor(df$virus_name_abbrev), "MS2 bacteriophage")
df$tail = relevel(factor(df$tail), "no")
df$structure = relevel(factor(df$structure), "non-enveloped")
df$symmetry = relevel(factor(df$symmetry), "icosahedral")
df$purification_level = factor(df$purification_level)

#drop high_chloride from df
df <- df[, !colnames(df) %in% c("high_chloride")]
df$high_chloride

# Set year to a numerical feature
df$year = as.numeric(df$year)

#add new chloride column for sensitivity analysis
df$chlor_sens_anal = df$high_chloride_new
df$chlor_sens_anal[df$phosphate_assumption_needed == 1 & df$chlor_sens_anal == TRUE] <- FALSE
df$chlor_sens_anal[df$phosphate_assumption_needed == 1 & df$chlor_sens_anal == FALSE] <- TRUE


# Set reference pH = 7 and temp  = 20 ----------------------------------------------
ref_pH = function(x) {x - 7.53}
ref_pH_squared = function(x) ref_pH(x)^2
ref_temp = function(x) {(x - 20)/5}
exponentiate = function(x) {10^x}

df$temp_5int = ref_temp(df$temp) #temperature in 5 degree intervals
df$pH_refpka = ref_pH(df$pH)
df$pH_refpka_squared = ref_pH_squared(df$pH)

#Scale features  ----------------------------------------------
df <- transform(df, genome_length_cs=scale(genome_length), diameter_cs=scale(diameter), year_cs = scale(year))

# Create Dummy Variables for each balt_class  ----------------------------------------------
df$dsDNA <- ifelse(df$balt_class == "dsDNA", 1, 0)
df$ssDNA <- ifelse(df$balt_class == "ssDNA", 1, 0)
df$plus_ssRNA <- ifelse(df$balt_class == "+ssRNA", 1, 0)
df$minus_ssRNA <- ifelse(df$balt_class == "-ssRNA", 1, 0)
df$dsRNA <- ifelse(df$balt_class == "dsRNA", 1, 0)

df$double_stranded <- ifelse(df$balt_class == "dsDNA" | df$balt_class == "dsRNA", 1, 0)
df$DNA <- ifelse(df$balt_class == "dsDNA" | df$balt_class == "ssDNA", 1, 0)
df$RNA <- ifelse(df$balt_class == "dsDNA" | df$balt_class == "ssDNA", 0, 1)

# Modify virus_name_abbrev  ----------------------------------------------
#encode variables for additional functions
df_name_mod <- df
df_name_mod$virus_name_abbrev_orig <- df_name_mod$virus_name_abbrev
df_name_mod$virus_name_abbrev <- gsub(" ", "_", df_name_mod$virus_name_abbrev)

#replace special characters that mess things up later on
df_name_mod$virus_name_abbrev = gsub("\\[|\\]", "", df_name_mod$virus_name_abbrev)
df_name_mod$virus_name_abbrev = gsub("\\]", "", df_name_mod$virus_name_abbrev)
df_name_mod$virus_name_abbrev <- gsub("\\.", "", df_name_mod$virus_name_abbrev)
df_name_mod$virus_name_abbrev <- gsub("\\/", "", df_name_mod$virus_name_abbrev)
df_name_mod$virus_name_abbrev <- gsub("\\-", "", df_name_mod$virus_name_abbrev)
df_name_mod$virus_name_abbrev <- gsub("\\'", "", df_name_mod$virus_name_abbrev)
df_name_mod$virus_name_abbrev <- gsub("\\(", "", df_name_mod$virus_name_abbrev)
df_name_mod$virus_name_abbrev <- gsub("\\)", "", df_name_mod$virus_name_abbrev)
df_name_mod$virus_name_abbrev <- tolower(df_name_mod$virus_name_abbrev)

df_name_mod$virus_name_abbrev = relevel(factor(df_name_mod$virus_name_abbrev), "ms2_bacteriophage")

levels(df_name_mod$virus_name_abbrev)

df_name_mod$paper_ID <- as.factor(df_name_mod$paper_ID)

# Create df_encoded for modeling: ----------------------------------------------

add_factor_cols = function(dat, col, stem) {
  to_add = unique(dat[[col]])
  for ( i in seq_along(to_add) ) {
    new_col = paste0(stem, "_", to_add[i])
    dat[[new_col]] = 1L * (dat[[col]] == to_add[i])
  }
  dat
}


df_test = add_factor_cols(df_name_mod, "high_chloride_new", "high_chloride_new")
colnames(df_test)

create_dataset_interactions = function(df, chloride_column){
  #chloride column must be chlor_sens_anal or high_chloride_new
  df_encoded = df
  df_encoded$high_chloride = df_encoded[[chloride_column]]
  df_encoded = add_factor_cols(df_encoded, "virus_name_abbrev", "virus")
  df_encoded = add_factor_cols(df_encoded, "buffer", "buffer")
  df_encoded = add_factor_cols(df_encoded, "purification_level", "purification")
  df_encoded = add_factor_cols(df_encoded, 'high_chloride', 'high_chloride')
  
  df_encoded$pH_refpka_x_high_chloride_TRUE = df_encoded$pH_refpka * df_encoded$high_chloride_TRUE
  df_encoded$pH_refpka_x_temp_5int = df_encoded$pH_refpka * df_encoded$temp_5int
  #add pH interactions
  df_encoded$balt_class_dsDNA_x_pH_refpka = df_encoded$dsDNA * df_encoded$pH_refpka
  df_encoded$balt_class_dsRNA_x_pH_refpka = df_encoded$dsRNA * df_encoded$pH_refpka
  df_encoded$balt_class_ssDNA_x_pH_refpka = df_encoded$ssDNA * df_encoded$pH_refpka
  df_encoded$balt_class_minus_ssRNA_x_pH_refpka = df_encoded$minus_ssRNA * df_encoded$pH_refpka
  df_encoded$DNA_x_pH_refpka = df_encoded$DNA * df_encoded$pH_refpka
  #add temp interactions
  df_encoded$balt_class_dsDNA_x_temp_5int = df_encoded$dsDNA * df_encoded$temp_5int
  df_encoded$balt_class_dsRNA_x_temp_5int = df_encoded$dsRNA * df_encoded$temp_5int
  df_encoded$balt_class_ssDNA_x_temp_5int = df_encoded$ssDNA * df_encoded$temp_5int
  df_encoded$balt_class_minus_ssRNA_x_temp_5int = df_encoded$minus_ssRNA * df_encoded$temp_5int
  df_encoded$DNA_x_temp_5int = df_encoded$DNA * df_encoded$temp_5int
  #add pH2 interactions
  df_encoded$balt_class_dsDNA_x_pH_refpka_squared = df_encoded$dsDNA * df_encoded$pH_refpka_squared
  df_encoded$balt_class_dsRNA_x_pH_refpka_squared = df_encoded$dsRNA * df_encoded$pH_refpka_squared
  df_encoded$balt_class_ssDNA_x_pH_refpka_squared = df_encoded$ssDNA * df_encoded$pH_refpka_squared
  df_encoded$balt_class_minus_ssRNA_x_pH_refpka_squared = df_encoded$minus_ssRNA * df_encoded$pH_refpka_squared
  df_encoded$DNA_x_pH_refpka_squared = df_encoded$DNA * df_encoded$pH_refpka_squared
  
  df_encoded <- df_encoded %>%
    group_by(virus_name_abbrev) %>%
    mutate(unique_paper_count = n_distinct(paper_ID)) %>%
    ungroup()
  
  
  df_virus_chars <- df_encoded[!duplicated(df_encoded[ , c("virus_name_abbrev")]),]
  #from that dataframe, extract only the columns to do with virus characteristics 
  df_virus_chars_subset <- df_virus_chars[c("virus_name_abbrev","dsDNA","dsRNA","ssDNA","plus_ssRNA","minus_ssRNA","DNA","RNA","balt_class", "family", "genus","species", "unique_paper_count")]
  #merge this with a mini data frame showing only the number of data points for each virus
  df_virus_chars_counts <- merge(x = df_virus_chars_subset, y = df_encoded %>% 
                                   count(virus_name_abbrev), by = "virus_name_abbrev", all.x = TRUE)
  
  return(list(df_encoded,df_virus_chars_counts))
}


df_encoded <- create_dataset_interactions(df = df_name_mod, chloride_column = 'high_chloride_new')[[1]]
df_virus_chars_counts <- create_dataset_interactions(df = df_name_mod, chloride_column = 'high_chloride_new')[[2]]
#View(df_virus_chars_counts)
df_virus_chars_counts$n_string <- paste0("n = ", df_virus_chars_counts$n, "; p = ", df_virus_chars_counts$unique_paper_count)

sort(colnames(df_encoded))

#Export df_encoded ----------------------------------------------

#save(df_encoded, file="12_14_23_df_encoded.RData")
#write.csv(df_encoded, "12_14_23_df_encoded.csv", row.names=FALSE)
nrow(df_encoded)
#View(df)

#Initial forward selection AIC  ----------------------------------------------
##Stage 1: check temperature and pH against base model  ----------------------------------------------
stage1_0 <- lm(log_average_kobs ~ virus_name_abbrev , data = df_encoded, REML = FALSE)
stage1_pH <- lm(log_average_kobs ~ virus_name_abbrev + pH_refpka , data = df_encoded, REML = FALSE)
stage1_temp <- lm(log_average_kobs ~ virus_name_abbrev + temp_5int , data = df_encoded, REML = FALSE)
stage1_temp_pH <- lm(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka , data = df_encoded, REML = FALSE)
AIC(stage1_0, stage1_pH, stage1_temp, stage1_temp_pH)
r.squaredGLMM(stage1_0)
r.squaredGLMM(stage1_temp)
r.squaredGLMM(stage1_pH)
r.squaredGLMM(stage1_temp_pH)

#stage1_temp_pH is the model to move forward with

colnames(df)
##Stage 2: check random effects  ----------------------------------------------
stage2_ranef_paperID <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + (1|paper_ID), data = df_encoded, REML = FALSE)
stage2_ranef_corr <-lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + (1|corr_author), data = df_encoded, REML = FALSE)
stage2_ranef_both <-lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + (1|corr_author) + (1|paper_ID), data = df_encoded, REML = FALSE)
AIC(stage2_ranef_paperID, stage2_ranef_corr, stage2_ranef_both)
r.squaredGLMM(stage2_ranef_paperID)
r.squaredGLMM(stage2_ranef_corr)
r.squaredGLMM(stage2_ranef_both)
#stage2_ranef_paperID is the mdodel to move forward with

##Stage 3: check additional effects  ----------------------------------------------
stage3_year <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + float_year_scaled + (1|paper_ID), data = df_encoded, REML = FALSE)
stage3_highchloride <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + high_chloride + (1|paper_ID), data = df_encoded, REML = FALSE)
stage3_buffer <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + buffer + (1|paper_ID), data = df_encoded, REML = FALSE)
stage3_alpha0 <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + alpha_0 + (1|paper_ID), data = df_encoded, REML = FALSE)
stage3_purification <-lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + purification_level + (1|paper_ID), data = df_encoded, REML = FALSE)

AIC(stage2_ranef_paperID, stage3_year, stage3_highchloride, stage3_buffer, stage3_alpha0, stage3_purification)
#stage3_highchloride is the mdodel to move forward with

r.squaredGLMM(stage2_ranef_paperID)
r.squaredGLMM(stage3_year)
r.squaredGLMM(stage3_highchloride)
r.squaredGLMM(stage3_buffer)
r.squaredGLMM(stage3_alpha0)
r.squaredGLMM(stage3_purification)

##Stage 4: interactions with support in the literature  ----------------------------------------------
stage4_inter_cl_pH <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka*high_chloride + (1|paper_ID), data = df_encoded, REML = FALSE)
stage4_inter_temp_pH <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + high_chloride + pH_refpka:temp_5int + (1|paper_ID), data = df_encoded, REML = FALSE)
stage4_bothint <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka + high_chloride + pH_refpka:temp_5int + pH_refpka:high_chloride + (1|paper_ID), data = df_encoded, REML = FALSE)

AIC(stage3_highchloride, stage4_inter_cl_pH, stage4_inter_temp_pH, stage4_bothint)
r.squaredGLMM(stage4_inter_cl_pH)

r.squaredGLMM(stage4_inter_temp_pH)
r.squaredGLMM(stage4_bothint)

r.squaredGLMM()
# move forward with stage4_bothint

tab_model(stage1_temp_pH, digits.re = 10, digits = 10, show.aic = TRUE)

# Primary model   ----------------------------------------------

M_best_subset <-lmer(log_average_kobs ~ 
                       virus_name_abbrev + 
                       temp_5int + 
                       pH_refpka + 
                       pH_refpka_squared + 
                       pH_refpka_x_high_chloride_TRUE + 
                       pH_refpka_x_temp_5int + 
                       balt_class_dsRNA_x_pH_refpka + 
                       DNA_x_pH_refpka + 
                       balt_class_dsDNA_x_temp_5int + 
                       balt_class_ssDNA_x_temp_5int + 
                       DNA_x_pH_refpka_squared + 
                       high_chloride_TRUE + 
                       (1 | paper_ID), data = df_encoded, REML = FALSE)

r.squaredGLMM(M_best_subset)
vif(M_best_subset)
library(performance)
check_collinearity(M_best_subset)

# Sensitivity analyses   ----------------------------------------------

## Forward selection - additional interactions test based on residuals   ----------------------------------------------

## test pH interactions
M1a <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka*high_chloride + pH_refpka:temp_5int + (1|paper_ID), data = df_encoded, REML = FALSE)
AIC(M1a)
r.squaredGLMM(M1a)
M1_allbalt = update(M1a, . ~ . + balt_class:pH_refpka)
AIC(M1a, M1_allbalt)
r.squaredGLMM(M1_allbalt)
#drops AIC from 701.7862 to 630.2217

M1_pH_DNA <- update(M1a, . ~ . + DNA:pH_refpka)
AIC(M1a, M1_pH_DNA)
r.squaredGLMM(M1_pH_DNA)
#drops AIC from 701.7862 to 640.6598

#restart by testing each baltimore class individually
M1_pH_dsDNA <- update(M1a, . ~ . + dsDNA:pH_refpka)
AIC(M1a, M1_pH_dsDNA)
r.squaredGLMM(M1_pH_dsDNA)
#drops AIC from 701.7862 to 654.0918 accepted

M1_pH_dsRNA <- update(M1_pH_dsDNA, . ~ . + dsRNA:pH_refpka)
AIC(M1a, M1_pH_dsDNA, M1_pH_dsRNA)
r.squaredGLMM(M1_pH_dsRNA)
#drops AIC from 654.0918 to 638.6197 accepted

M1_pH_ssDNA <- update(M1_pH_dsRNA, . ~ . + ssDNA:pH_refpka)
AIC(M1a, M1_pH_dsDNA, M1_pH_dsRNA, M1_pH_ssDNA)
r.squaredGLMM(M1_pH_ssDNA)
#drops AIC from 638.6197 to 628.2426 Accepted

M1_pH_minus_ssRNA <- update(M1_pH_ssDNA, . ~ . + minus_ssRNA:pH_refpka)
AIC(M1a, M1_pH_dsDNA, M1_pH_dsRNA, M1_pH_ssDNA, M1_pH_minus_ssRNA)
r.squaredGLMM(M1_pH_minus_ssRNA)
#raises AIC from 628.2426 to 630.2217 Rejected.

## test temperature interactions
M2a <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka*high_chloride +  pH_refpka:temp_5int + dsDNA:pH_refpka + dsRNA:pH_refpka + ssDNA:pH_refpka +(1|paper_ID), data = df_encoded, REML = FALSE)
AIC(M2a)

#test DNA and temp interaction
M2a_DNA <- update(M2a, . ~ . + DNA:temp)
AIC(M2a, M2a_DNA)
r.squaredGLMM(M2a_DNA)

#lowered AIC from 628.2426 to 621.2949

#test baltimore class and temp interactions
M2a_all_balt <- update(M2a, . ~ . + balt_class:temp)
AIC(M2a, M2a_all_balt)
#rank deficient

#test individual baltimore classes and temperature interactions
M2a_temp_dsDNA <- update(M2a, . ~ . + dsDNA:temp_5int)
AIC(M2a, M2a_temp_dsDNA)
r.squaredGLMM(M2a_temp_dsDNA)
#raises AIC from 628.2426 to 629.1038. Rejected

M2a_temp_dsRNA <- update(M2a, . ~ . + dsRNA:temp_5int)
AIC(M2a, M2a_temp_dsRNA)
r.squaredGLMM(M2a_temp_dsRNA)
#raises AIC from 628.2426 to 629.2998. Rejected

M2a_temp_ssDNA <- update(M2a, . ~ . + ssDNA:temp_5int)
AIC(M2a, M2a_temp_ssDNA)
r.squaredGLMM(M2a_temp_ssDNA)
#lowers AIC from 628.2426 to 613.2465

M2a_temp_minus_ssRNA <-update(M2a_temp_ssDNA, . ~ . + minus_ssRNA:temp_5int)
AIC(M2a, M2a_temp_ssDNA, M2a_temp_minus_ssRNA)
#rank deficient

## test pH squared interactions
M3 <-lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka*high_chloride + pH_refpka:temp_5int + dsDNA:pH_refpka + dsRNA:pH_refpka + ssDNA:pH_refpka + ssDNA:temp_5int  +(1|paper_ID), data = df_encoded, REML = FALSE)
M3_ph_squared <- update(M3, .~. + pH_refpka_squared )
AIC(M3, M3_ph_squared)
r.squaredGLMM(M3_ph_squared)
#brings AIC down from 613.2465 to 599.6007

#test pH squared all baltimore class
M3_ph_squared_balt_class <- update(M3_ph_squared, .~. + pH_refpka_squared:balt_class )
AIC(M3_ph_squared, M3_ph_squared_balt_class)
r.squaredGLMM(M3_ph_squared_balt_class)
#brings AIC down from 599.6007 to 590.3595

#test pH squared interactions with DNA 
M3_ph_squared_DNA <- update(M3_ph_squared, .~. +  pH_refpka_squared:DNA)
AIC(M3_ph_squared, M3_ph_squared_DNA)
r.squaredGLMM(M3_ph_squared_DNA)
#brings AIC down from 599.6007 to 586.9065

#test pH squared interactions individual baltimore classes 

M3_ph_squared_dsDNA <-update(M3_ph_squared, .~. + pH_refpka_squared:dsDNA )
AIC(M3, M3_ph_squared, M3_ph_squared_dsDNA)
r.squaredGLMM(M3_ph_squared_dsDNA)
#brings AIC down from 599.6007 to 592.6730 Accepted

M3_ph_squared_dsRNA <- update(M3_ph_squared_dsDNA, .~. + pH_refpka_squared:dsRNA )
AIC(M3, M3_ph_squared, M3_ph_squared_dsDNA, M3_ph_squared_dsRNA)
r.squaredGLMM(M3_ph_squared_dsRNA)
#raises AIC from 592.6730 to 593.2286 Rejected

M3_ph_squared_ssDNA <- update(M3_ph_squared_dsDNA, .~. + pH_refpka_squared:ssDNA )
AIC(M3, M3_ph_squared, M3_ph_squared_dsDNA, M3_ph_squared_ssDNA)
r.squaredGLMM(M3_ph_squared_ssDNA)
#lowers AIC from 592.6730 to 587.5555 Accepted

M3_ph_squared_minus_ssRNA <- update(M3_ph_squared_ssDNA, .~. + pH_refpka_squared:minus_ssRNA )
AIC(M3, M3_ph_squared, M3_ph_squared_ssDNA, M3_ph_squared_minus_ssRNA)
r.squaredGLMM(M3_ph_squared_ssDNA)
#raises aic FROM 587.5555 to 589.5496. Rejected

## forward selection final model    ----------------------------------------------
M_forward <-lmer(log_average_kobs ~ virus_name_abbrev + 
                 temp_5int + 
                 pH_refpka*high_chloride +  
                 pH_refpka:temp_5int +
                 dsDNA:pH_refpka + 
                 dsRNA:pH_refpka + 
                 ssDNA:pH_refpka + 
                 ssDNA:temp_5int + 
                 pH_refpka_squared + 
                 pH_refpka_squared:DNA + 
                 (1|paper_ID), data = df_encoded, REML = FALSE)
AIC(M_forward)
BIC(M_final)
r.squaredGLMM(M_final)
fixef(M_final)
vif(M_final)
#View(df_encoded)

## Redo primary model with virus features     ----------------------------------------------

final_model_novirus <-lmer(log_average_kobs ~ 
                             temp_5int + 
                             pH_refpka + 
                             pH_refpka_squared + 
                             pH_refpka_x_high_chloride_TRUE + 
                             pH_refpka_x_temp_5int + 
                             balt_class_dsRNA_x_pH_refpka + 
                             DNA_x_pH_refpka + 
                             balt_class_dsDNA_x_temp_5int + 
                             balt_class_ssDNA_x_temp_5int + 
                             DNA_x_pH_refpka_squared + 
                             high_chloride_TRUE + 
                             (1 | paper_ID), data = df_encoded, REML = FALSE)

AIC(M_best_subset, final_model_novirus)
r.squaredGLMM(final_model_novirus)
#AIC of model without virus name is 902.2032

M2_vir <- update(final_model_novirus, .~. +  virus_name_abbrev)

M2_spe <- update(final_model_novirus, .~. +  species)
M2_genus <- update(final_model_novirus, .~. +  genus)
M2_fam <- update(final_model_novirus, .~. +  family)

M2_bal <- update(final_model_novirus, .~. +  balt_class)
M2_dia <- update(final_model_novirus, .~. +  diameter)
M2_g_len <- update(final_model_novirus, .~. +  genome_length)

M2_tai <- update(final_model_novirus, .~. +  tail)
M2_sym <- update(final_model_novirus, .~. +  symmetry)
M2_str <- update(final_model_novirus, .~. +  structure)

AIC(final_model_novirus, M2_vir, M2_spe, M2_genus, M2_fam, M2_bal, M2_dia, M2_g_len, M2_tai, M2_sym, M2_str)
r.squaredGLMM(final_model_novirus)

r.squaredGLMM(M2_tai)
r.squaredGLMM(M2_sym)
r.squaredGLMM(M2_str)

r.squaredGLMM(M2_dia)
r.squaredGLMM(M2_g_len)

r.squaredGLMM(M2_bal)

r.squaredGLMM(M2_fam)
r.squaredGLMM(M2_genus)
r.squaredGLMM(M2_spe)

nobs(M2_vir, M2_spe, M2_genus, M2_fam, M2_bal, M2_dia, M2_g_len, M2_tai, M2_sym, M2_str, M2_CG, M2_C, M2_G, M2_A, M2_T, M2_U)

##best subset BIC  ----------------------------------------------

M_best_subset_bic <- lmer(log_average_kobs ~ 
                            virus_name_abbrev + 
                            temp_5int + 
                            pH_refpka + 
                            pH_refpka_x_high_chloride_TRUE + 
                            balt_class_dsRNA_x_pH_refpka + 
                            DNA_x_pH_refpka + 
                            balt_class_ssDNA_x_temp_5int + 
                            DNA_x_pH_refpka_squared + 
                            high_chloride_TRUE +
                            (1 | paper_ID), data = df_encoded, REML = FALSE)


AIC(M_best_subset_bic)
BIC(M_best_subset_bic)

r.squaredGLMM(M_best_subset_bic)

saveRDS(M_best_subset, file = "M_best_subset.rds")

# Multicollinearity      ----------------------------------------------
cor(df_encoded$alpha_0, df_encoded$pH_refpka, method = 'pearson') 


# Function to Predict inactivation rate      ----------------------------------------------
predict_k <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  temp_5int <- rep(ref_temp(temp),82)
  pH_refpka <- rep(ref_pH(pH),82)
  pH_refpka_squared <- rep(ref_pH_squared(pH),82)
  paper_ID <- rep(paper_ID, 82) # for a check
  if (chloride == "TRUE") {
    high_chloride_TRUE <- rep(1,82)
    chloride_level = 'high'
  } else if (chloride == "FALSE") {
    high_chloride_TRUE <- rep(0,82)
    chloride_level = 'low'
  } else {
    print("invalid chloride value")
  }
  
  df_ref <- data.frame(df_virus_counts, temp_5int, pH_refpka, pH_refpka_squared, high_chloride_TRUE, paper_ID)
  
  df_ref$pH_refpka_x_high_chloride_TRUE = df_ref$pH_refpka * df_ref$high_chloride_TRUE
  df_ref$pH_refpka_x_temp_5int = df_ref$pH_refpka * df_ref$temp_5int
  #add pH interactions
  df_ref$balt_class_dsDNA_x_pH_refpka = df_ref$dsDNA * df_ref$pH_refpka
  df_ref$balt_class_dsRNA_x_pH_refpka = df_ref$dsRNA * df_ref$pH_refpka
  df_ref$balt_class_ssDNA_x_pH_refpka = df_ref$ssDNA * df_ref$pH_refpka
  df_ref$balt_class_minus_ssRNA_x_pH_refpka = df_ref$minus_ssRNA * df_ref$pH_refpka
  df_ref$DNA_x_pH_refpka = df_ref$DNA * df_ref$pH_refpka
  #add temp interactions
  df_ref$balt_class_dsDNA_x_temp_5int = df_ref$dsDNA * df_ref$temp_5int
  df_ref$balt_class_dsRNA_x_temp_5int = df_ref$dsRNA * df_ref$temp_5int
  df_ref$balt_class_ssDNA_x_temp_5int = df_ref$ssDNA * df_ref$temp_5int
  df_ref$balt_class_minus_ssRNA_x_temp_5int = df_ref$minus_ssRNA * df_ref$temp_5int
  df_ref$DNA_x_temp_5int = df_ref$DNA * df_ref$temp_5int
  #add pH2 interactions
  df_ref$balt_class_dsDNA_x_pH_refpka_squared = df_ref$dsDNA * df_ref$pH_refpka_squared
  df_ref$balt_class_dsRNA_x_pH_refpka_squared = df_ref$dsRNA * df_ref$pH_refpka_squared
  df_ref$balt_class_ssDNA_x_pH_refpka_squared = df_ref$ssDNA * df_ref$pH_refpka_squared
  df_ref$balt_class_minus_ssRNA_x_pH_refpka_squared = df_ref$minus_ssRNA * df_ref$pH_refpka_squared
  df_ref$DNA_x_pH_refpka_squared = df_ref$DNA * df_ref$pH_refpka_squared
  #add sensitivity analysis interactions
  
  
  #predictions <- data.frame(df_virus_chars_counts$virus_name_abbrev, predict(M_final, newdata = df_ref, re.form = ~0))
  
  form_best <- formula(model)
  f_final <- update(form_best, . ~ . - (1 | paper_ID))
  f_final <- update(f_final, NULL ~ .)
  
  x <- model.matrix(f_final, data = df_ref )
  
  #stopifnot(all(colnames(x) == names(fixef(M_final))))
  
  est_k <- x %*% fixef(model)
  est_SE <- x %*% vcov(model) %*% t(x) |> diag() |> sqrt()
  z <- qnorm(.975)
  est_lower <- est_k - z*est_SE
  est_upper <- est_k + z*est_SE
  
  #make a new dataframe with predictions and confidence intervals
  
  pred_with_CI <- cbind(df_virus_counts, est_k, est_lower, est_upper) 
  return(pred_with_CI)
  
}

predict_k(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)  %>% filter(virus_name_abbrev == 'ms2_bacteriophage')

predict_k(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)  %>% filter(virus_name_abbrev == 'hpai_a2005')

colnames(df_virus_chars_counts)
unique(df_virus_chars_counts$virus_name_abbrev)
10^2.250861

# Find confidence intervals for fold changes---------------------------------------------

#For example, a pH increase from 7.53 to 10 results in a predicted rate constant decreasing 4.03-fold for (+)ssRNA viruses, 1.27-fold for dsRNA viruses, and 61.9-fold for dsDNA viruses. 
#For example, a pH decrease from 10 to 7.53 results in a predicted rate constant increase 4.03-fold for (+)ssRNA viruses, 1.27-fold for dsRNA viruses, and 61.9-fold for dsDNA viruses. 

b <- fixef(M_best_subset)
x <- matrix(0,3,length(b))
colnames(x) <- names(b)
names_pH_cols <- grep("pH",names(b), value=TRUE)
x[,names_pH_cols[1:2]] <- rep(c(10-7.53, ref_pH_squared(10)-ref_pH_squared(7.53)),each=3)
x[2,"balt_class_dsRNA_x_pH_refpka"] <- 10-7.53
x[3,"DNA_x_pH_refpka"] <- 10-7.53
x[3,"DNA_x_pH_refpka_squared"] <-  ref_pH_squared(10)-ref_pH_squared(7.53) 
x <- -1*x
log_fc <- (x %*% b)
#standard errors for fold change
se_log_fc <- x %*% vcov(M_best_subset) %*% t(x) |> diag() |> sqrt()

z <- qnorm(.975)
est_lower_fc <- log_fc - z*se_log_fc
est_upper_fc <- log_fc + z*se_log_fc

10^est_lower_fc
10^est_upper_fc



b <- fixef(M_best_subset)
x <- matrix(0,3,length(b))
colnames(x) <- names(b)
names_pH_cols <- grep("pH",names(b), value=TRUE)
x[,names_pH_cols[1:2]] <- rep(c(9-6, ref_pH_squared(9)-ref_pH_squared(6)),each=3)
x[2,"balt_class_dsRNA_x_pH_refpka"] <- 9-6
x[3,"DNA_x_pH_refpka"] <- 9-6
x[3,"DNA_x_pH_refpka_squared"] <-  ref_pH_squared(9)-ref_pH_squared(6) 
x <- -1*x
log_fc <- (x %*% b)
#standard errors for fold change
se_log_fc <- x %*% vcov(M_best_subset) %*% t(x) |> diag() |> sqrt()

z <- qnorm(.975)
est_lower_fc <- log_fc - z*se_log_fc
est_upper_fc <- log_fc + z*se_log_fc

10^log_fc
10^est_lower_fc
10^est_upper_fc

b[names_pH_cols]

names_pH_cols


# Count whether EPA CT values are exceeded: Table 1      ----------------------------------------------

under_required_k <- function(model, temp_input, pH_input, chloride_input, required_CT){
  prediction <- predict_k(model = model, df_virus_counts = df_virus_chars_counts, temp = temp_input, pH = pH_input, chloride = chloride_input, paper_ID = 32)
  required_k <- -log(1/10000)/required_CT
  num_under_req <- sum(exponentiate(prediction$est_k) < required_k)
  percent_under_req <- num_under_req/82*100
  return(list(num_under_req, percent_under_req))
}

under_required_k(model = M_best_subset, temp_input = 25, pH_input = 10, chloride_input = "FALSE", required_CT = 15)

# plot estimates for all Baltimore classes and (+)ssRNA only. Figure 3      ----------------------------------------------
plot_predictions_all <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  
  plot <- ggplot(subset(pred_with_CI ), aes(x = est_k, y = reorder(virus_name_abbrev, est_k), color = balt_class, shape=balt_class)) + 
    geom_point() +
    scale_x_continuous(breaks=c(-1, 0, 1 ,2, 3),
                       labels=c("0.1", "1", "10", "100","1000")) +
    geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    #xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Inactivation rate constants (L/mg*min)") +
    scale_color_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values=c("#75bbfd", "#fcb001", "#02ab2e",  "#7e1e9c", "#c44240")) +
    scale_shape_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values = c(0,1,2,4,7)) +
    #when the y axis should be blank
    #theme(legend.position = "right", legend.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
    #when the y axis should have virus name
    theme(legend.position = c(0.8, 0.2), legend.title = element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #theme(legend.position = c(0.8, 0.2), legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  #coord_flip()
  #theme(legend.position = "bottom")
  
  #print(plot)
  return(plot)
}

# plot rate constants from review   ----------------------------------------------
plot_initial_points <- function(df, model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  df_test <- merge(x = df, y = pred_with_CI, by = "virus_name_abbrev", all.x = TRUE)
  df_test$balt_class <- df_test$balt_class.x
  plot <- ggplot(df_test, aes(x = log_average_kobs, y = reorder(virus_name_abbrev_orig, est_k), color = balt_class)) + 
    geom_point(alpha=.4) +
    geom_errorbar(data = df_test, aes(xmin = est_lower, xmax = est_upper)) +
    geom_point(data = df_test, aes(x = est_k, y = reorder(virus_name_abbrev_orig, est_k)), color='black', shape = "|", size = 2) +
    scale_x_continuous(breaks=c(-3, -2, -1, 0, 1 ,2, 3),
                       labels=c("0.001","0.01","0.1", "1", "10", "100","1000")) +
    #geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    #xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Inactivation rate constants (L/mg*min)") +
    scale_color_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values=c("#75bbfd", "#fcb001", "#02ab2e",  "#7e1e9c", "#c44240")) +
    scale_shape_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values = c(0,1,2,4,7)) +
    theme_bw()+
    theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
  #theme(legend.position = "")
  #theme(legend.position = c(0.8, 0.2), legend.title = element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())
  return(plot)
}

#plot estimates and predictions on top of each other
jpeg(file="plot1_and_2a.jpeg",width=8,height=11,units="in", res=300)
plot_initial_points(df = df_encoded, model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
dev.off()

rev_plot_initial_just_points <- function(df, model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  df_test <- merge(x = df, y = pred_with_CI, by = "virus_name_abbrev", all.x = TRUE)
  df_test$balt_class <- df_test$balt_class.x
  print
  plot <- ggplot(df_test, aes(x = log_average_kobs, y = reorder(virus_name_abbrev_orig, est_k), color = balt_class)) + 
    geom_point(alpha=.5) +
    #geom_errorbar(data = df_test, aes(xmin = est_lower, xmax = est_upper)) +
    #geom_point(data = df_test, aes(x = est_k, y = reorder(virus_name_abbrev_orig, est_k)), color='black', shape = "|", size = 2) +
    scale_x_continuous(breaks=c(-3, -2, -1, 0, 1 ,2, 3),
                       labels=c("0.001","0.01","0.1", "1", "10", "100","1000")) +
    scale_y_discrete(position = "right") +
    #geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    #xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Observed Rate Constants (L/mg*min)") +
    scale_color_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values=c("#75bbfd", "#fcb001", "#02ab2e",  "#7e1e9c", "#c44240")) +
    scale_shape_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values = c(0,1,2,4,7)) +
    theme_bw()+
    theme(legend.position = c(0.8, 0.2), legend.title = element_blank(),axis.ticks.y=element_blank(), axis.title=element_text(size=16,face="bold"),  plot.margin = unit(c(0, 0, 0, -1), "pt"))
  #theme(legend.position = "")
  #theme(legend.position = c(0.8, 0.2), legend.title = element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())
  return(plot)
}

rev_plot_modeled <- function(df, model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  df_test <- merge(x = pred_with_CI, y = df, by = "virus_name_abbrev", all.x = TRUE)
  #  df_test <- merge(x = df, y = pred_with_CI, by = "virus_name_abbrev", all.x = TRUE)
  #df_test$balt_class <- df_test$balt_class.x
  df_test$balt_class <- df_test$balt_class.y
  df_test <- select(df_test, balt_class, virus_name_abbrev_orig, est_k, est_lower, est_upper, n_string)
  df_test <- distinct(df_test)
  df_test <- df_test[order(df_test$est_k),]
  
  print(nrow(df_test))
  plot <- ggplot(df_test, aes(x = est_k, y = reorder(virus_name_abbrev_orig, est_k), color = balt_class)) + 
    geom_point() +
    geom_errorbar(data = df_test, aes(xmin = est_lower, xmax = est_upper)) +
    #geom_point(data = df_test, aes(x = est_k, y = reorder(virus_name_abbrev_orig, est_k)), color='black', shape = "|", size = 2) +
    scale_x_continuous(breaks=c(-3, -2, -1, 0, 1 ,2, 3),
                       labels=c("0.001","0.01","0.1", "1", "10", "100","1000")) +
    #geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    #xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Estimated Rate Constants (L/mg*min)") +
    scale_color_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values=c("#75bbfd", "#fcb001", "#02ab2e",  "#7e1e9c", "#c44240")) +
    scale_shape_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values = c(0,1,2,4,7)) +
    theme_bw()+
    #theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
    #theme(legend.position = "")
    theme(legend.position = c(0.8, 0.2), legend.title = element_blank(), axis.ticks.y=element_blank(), axis.title=element_text(size=16,face="bold"), plot.margin = unit(c(0, 0, 0, -1), "pt")) + 
    scale_y_discrete(labels = df_test$n_string, position = "right") 
  
  return(plot)
}

just_model<- rev_plot_modeled(df = df_encoded, model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
just_model
just_points <- rev_plot_initial_just_points(df = df_encoded, model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
just_points

#sort(unique(df_encoded$virus_name_abbrev))

jpeg(file="plots_together.jpeg",width=11.5,height=11,units="in", res=300)
just_points +just_model
dev.off()

plot_predictions_ssRNA <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  
    plot <- ggplot(subset(pred_with_CI, balt_class == "+ssRNA" ), aes(x = est_k, y = reorder(virus_name_abbrev, est_k), color = genus, shape=genus)) + 
    geom_point() +
    scale_x_continuous(breaks=c(-1, 0, 1 ,2, 3),
                       labels=c("0.1", "1", "10", "100","1000")) +
    geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    #xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Estimated rate constants (L/mg*min)") +
    scale_color_manual(breaks = c("Enterovirus", "Hepatovirus", "Emesvirus", "Qubevirus", "Norovirus","Vesivirus", "Paslahepevirus" ),
                       values=c("#087804", "#0bf77d", "#9e003a",  "#ff474c", "#95d0fc", "#152eff", "#fac205")) +
    scale_shape_manual(breaks = c("Enterovirus", "Hepatovirus", "Emesvirus", "Qubevirus", "Norovirus","Vesivirus", "Paslahepevirus"),
                       values = c(0,1,2,4,7, 15, 17, 8)) +
    #enterovirus and hepatovirus are picornaviridae and are green
    #emesvirus and qubevirus is fiersviridae and are red
    #norovirus and vesivirus are caliciviridae and are yellow
    #paslahepevirus is hepviridae and is blue
    
    #when the y axis should be blank
    #theme(legend.position = "right", legend.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
    #when the y axis should have virus name
    #theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
    theme(legend.position = c(0.8, 0.2), legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title=element_text(size=16,face="bold"))
  #theme(legend.position = "bottom")
  
  #print(plot)
  return(plot)
}

general_plot_baltclass <- plot_predictions_ssRNA(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
general_plot_baltclass

jpeg(file="plot2_jan8.jpeg",width=5,height=6,units="in", res=300)
general_plot_baltclass
dev.off()

# Function to plot predictions for viruses with more than 5 data points. Figure 4    ----------------------------------------------
plot_predictions_five <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  
  if (chloride == "TRUE") {
    chloride_level = 'high'
  } else if (chloride == "FALSE") {
    chloride_level = 'low'
  } else {
    print("invalid chloride value")
  }
  
  plot <- ggplot(subset(pred_with_CI, n>4 ), aes(x = est_k, y = reorder(virus_name_abbrev, est_k), color = balt_class, shape=balt_class)) + 
    geom_point() +
    scale_x_continuous(breaks=c(-1, 0, 1 ,2, 3),
                       labels=c(".1", "1", "10", "100","1000")) +
    geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    #coord_flip()+
    xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Inactivation rate constant\n (L/mg*min)") +
    scale_color_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values=c("#75bbfd", "#fcb001", "#02ab2e",  "#7e1e9c", "#c44240")) +
    scale_shape_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values = c(0,1,2,4,7)) +
    #when the y axis should be blank
    #theme(legend.position = "right", legend.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
    #when the y axis should have virus name +
    theme(legend.position = "right", legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  #print(plot)
  return(plot)
}

#Figure 4a. Reference. 
cond_1 <- plot_predictions_five(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
cond_1 <- cond_1 + scale_x_continuous(breaks=c(-1, 0, 1 ,2, 3),
                                      labels=c("0.1", "1", "10", "100","1000")) +
  theme(legend.position = "")
#labs(x=("20C, pH 7.53, and low chloride"))

#Figure 4b. Condition 2: Reference and high pH
cond_2 <- plot_predictions_five(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 10, chloride = "FALSE", paper_ID = 32)
cond_2 <- cond_2 + scale_x_continuous(breaks=c( -2, -1 ,0, 1, 2),
                                      labels=c("0.01", "0.1", "1","10", "100"),
                                      limits=c(-2, 4.5)) +
  theme(legend.position = "")
#labs(x=(expression(paste("20 C, ", bold("pH 10"), ", and low chloride"))))

#Figure 4c. Condition 3: Reference and low temperature 
cond_3 <- plot_predictions_five(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 5, pH = 7.53, chloride = "FALSE", paper_ID = 32)
cond_3 <- cond_3 + scale_x_continuous(breaks=c( 0, 1 ,2, 3),
                                      labels=c("1", "10", "100","1000")) +
  theme(legend.position = "")
#labs(x=(expression(paste(bold("5C"), ", pH 7.53, and low chloride"))))

#Figure 4d. Condition 3: Reference and high chloride 
cond_4 <- plot_predictions_five(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "TRUE", paper_ID = 32)
cond_4 <- cond_4 + scale_x_continuous(breaks=c(0, 1,2,3,4),
                                      labels=c("1","10", "100", "1000","10000")) +
  theme(legend.position = "")
#labs(x=(expression(paste("20 C, pH 7.53 and ", bold("high chloride")))))


#for legend extraction
cond_4_with_legend <- plot_predictions_five(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "TRUE", paper_ID = 32)+
  scale_x_continuous(breaks=c(0, 1,2,3,4),
                     labels=c("1","10", "100", "1000","10000"))+
  theme(legend.position="bottom")
library(grid)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
legend <- g_legend(cond_4_with_legend) 
grid.newpage()
grid.draw(legend) 

design <- "
  1122
  1122
  3344
  3344
"

cond_1 + cond_2 + cond_3 + cond_4 + grid.draw(legend) + plot_layout(design = design) 

jpeg(file="plot3_condition.jpeg",width=6,height=6,units="in", res=300)
cond_1 + cond_2 + cond_3 + cond_4
dev.off()

# Functions to return estimates   ----------------------------------------------

return_k <- function(virus_name, model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df_virus_counts, temp, pH, chloride, paper_ID)
  pred <- (pred_with_CI %>% filter(virus_name_abbrev == virus_name))['est_k']
  pred_min_CI <-(pred_with_CI %>% filter(virus_name_abbrev == virus_name))['est_lower']
  pred_max_CI <- (pred_with_CI %>% filter(virus_name_abbrev == virus_name))['est_upper']
  return( list(exponentiate(pred), exponentiate(pred_min_CI), exponentiate(pred_max_CI)))
}

return_k(virus_name = 'ms2_bacteriophage', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
return_k(virus_name = 'ipnv', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 10, pH = 9, chloride = "TRUE", paper_ID = 32)

return_k_ordered <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  #functions to take in experimental conditions and return a dataframe of every virus with inactivation rates predicted under those conditions, ordered from lowest to highest inactivation rate
  pred <- predict_k(model, df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI <- pred[c("virus_name_abbrev","balt_class","n","est_k","est_lower","est_upper")]
  pred_with_CI[c("est_k","est_lower","est_upper")] <- lapply(pred_with_CI[c("est_k","est_lower","est_upper")], exponentiate)
  pred_with_CI <- pred_with_CI[order(pred_with_CI$est_k), ]
  return(pred_with_CI)
  
}

#return an ordered dataframe containing predictions

return_k_df <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  #functions to take in experimental conditions and return a dataframe of every virus with inactivation rates predicted under those conditions, ordered alphabetically by virus name
  pred <- predict_k(model, df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI <- pred[c("virus_name_abbrev","balt_class","n","est_k","est_lower","est_upper")]
  pred_with_CI[c("est_k","est_lower","est_upper")] <- lapply(pred_with_CI[c("est_k","est_lower","est_upper")], exponentiate)
  pred_with_CI <- pred_with_CI[order(pred_with_CI$virus_name_abbrev), ]
  return(pred_with_CI)
  
}

#Export Estimate Dataset   ----------------------------------------------
return_k_ordered(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
ref_cond_pred_ordered <- return_k_ordered(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
virus_names <- df_encoded[, c("virus_name_abbrev", "virus_name_abbrev_orig", "virus_name_full")]

ref_cond_pred_ordered <- unique(merge(x = ref_cond_pred_ordered, y = virus_names, by.x="virus_name_abbrev", by.y="virus_name_abbrev", all.x=TRUE))
ref_cond_pred_ordered <- ref_cond_pred_ordered[, c("virus_name_abbrev_orig", "balt_class","n","est_k","est_lower","est_upper")]
ref_cond_pred_ordered <- ref_cond_pred_ordered[order(ref_cond_pred_ordered$est_k), ]

#names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="virus_name_abbrev"] <- "Virus Abbreviation"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="virus_name_abbrev_orig"] <- "Virus Name"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="balt_class"] <- "Baltimore Class"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="n"] <- "Number of Data Points"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="est_k"] <- "Inactivation Rate Constants (L/mg*min)"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="est_lower"] <- "Lower Bound of 95% CI"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="est_upper"] <- "Upper Bound of 95% CI"

ref_cond_pred_ordered
#View(ref_cond_pred_ordered)
write.csv(ref_cond_pred_ordered, "C:\\Users\\mirac\\Desktop\\ref_cond_pred_ordered_full.csv", row.names=FALSE)
write_xlsx(ref_cond_pred_ordered, "C:\\Users\\mirac\\Desktop\\ref_cond_pred_ordered_full.xlsx")


names(virus_names)[names(virus_names)=="virus_name_abbrev"] <- "Virus Abbreviation for Modeling"
names(virus_names)[names(virus_names)=="virus_name_abbrev_orig"] <- "Virus Abbreviation"
names(virus_names)[names(virus_names)=="virus_name_full"] <- "Virus Name"
write.csv(unique(virus_names), "C:\\Users\\mirac\\Desktop\\virus_abbreviations.csv", row.names=FALSE)
write.csv(unique(virus_names), "C:\\Users\\mirac\\Desktop\\virus_abbreviations.csv", row.names=FALSE)

write_xlsx(unique(virus_names), "C:\\Users\\mirac\\Desktop\\virus_abbreviations.xlsx")



ref_cond_pred_ordered


df$virus_name_abbrev
df_encoded$virus_name_abbrev

# chloride sensitivity analysis   ----------------------------------------------

df_encoded_chloride_sens_anal <- create_dataset_interactions(df = df_name_mod, chloride_column = 'chlor_sens_anal')[[1]]
df_virus_chars_counts_chloride_sens_anal <- create_dataset_interactions(df = df_name_mod, chloride_column = 'chlor_sens_anal')[[2]]

M_best_subset_chloride_sens_anal <-lmer(log_average_kobs ~ 
                                 virus_name_abbrev + 
                                 temp_5int + 
                                 pH_refpka + 
                                 pH_refpka_squared + 
                                 pH_refpka_x_high_chloride_TRUE + 
                                 pH_refpka_x_temp_5int + 
                                 balt_class_dsRNA_x_pH_refpka + 
                                 DNA_x_pH_refpka + 
                                 balt_class_dsDNA_x_temp_5int + 
                                 balt_class_ssDNA_x_temp_5int + 
                                 DNA_x_pH_refpka_squared + 
                                 high_chloride_TRUE + 
                                 (1 | paper_ID), data = df_encoded_chloride_sens_anal, REML = FALSE)

AIC(M_best_subset, M_best_subset_chloride_sens_anal)

#viruses affected by data changes
viruses_changed <- c('cvb5_faulkner', 'cvb2', 'hadv_41_tak','hrv_wa','ms2_bacteriophage' ,'mnv_1')           

predictions_initial <- return_k_df(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)

predictions_initial
#predictions_initial_important_viruses <- predictions_initial[predictions_initial$virus_name_abbrev %in% viruses_changed, ]
colnames(predictions_initial) <- paste0('initial_', colnames(predictions_initial))
predictions_chloride_sens_anal <- return_k_df(model = M_best_subset_chloride_sens_anal, df_virus_counts = df_virus_chars_counts_chloride_sens_anal, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
#predictions_new_important_viruses <- predictions_new[predictions_new$virus_name_abbrev %in% viruses_changed, ]
colnames(predictions_chloride_sens_anal) <- paste0('new_', colnames(predictions_chloride_sens_anal))
merged_predictions <- merge(x = predictions_initial, y = predictions_chloride_sens_anal, by.x = "initial_virus_name_abbrev", by.y = "new_virus_name_abbrev", all = TRUE)
merged_predictions$percent_change <- (merged_predictions$new_est_k - merged_predictions$initial_est_k)/merged_predictions$initial_est_k*100
merged_predictions_important <- merged_predictions[merged_predictions$initial_virus_name_abbrev %in% viruses_changed, ]

chloride_sens_anal_table <- select(merged_predictions, initial_virus_name_abbrev, initial_est_k, new_est_k, percent_change)
chloride_sens_anal_table_greater_than_percent <- chloride_sens_anal_table[abs(chloride_sens_anal_table$percent_change) >3, ]
chloride_sens_anal_table_greater_than_percent<- chloride_sens_anal_table_greater_than_percent[order(chloride_sens_anal_table_greater_than_percent$percent_change), ]
write.csv(chloride_sens_anal_table_greater_than_percent, "C:\\Users\\mirac\\Desktop\\chloride_sens_anal_table_greater_than_percent.csv", row.names=FALSE)


return_k(virus_name = 'hadv_5', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)[[1]]


fixef_init <- fixef(M_best_subset)
fixef_final <- fixef(M_best_subset_chloride_sens_anal)

chlor_sens_anal_change_coef <- data.frame(names(fixef_init), fixef_init, fixef_final)
chlor_sens_anal_change_coef$percent_change <- (chlor_sens_anal_change_coef$fixef_final - chlor_sens_anal_change_coef$fixef_init)/chlor_sens_anal_change_coef$fixef_init*100
chlor_sens_anal_change_coef_greater_than_percent <- chlor_sens_anal_change_coef[abs(chlor_sens_anal_change_coef$percent_change) >3, ]
chlor_sens_anal_change_coef_greater_than_percent
write.csv(chlor_sens_anal_change_coef_greater_than_percent, "C:\\Users\\mirac\\Desktop\\chlor_sens_anal_change_coef_greater_than_percent.csv", row.names=FALSE)

sort(chlor_sens_anal_change_coef)
chlor_sens_anal_change_coef[order(chlor_sens_anal_change_coef$names.fixef_init.), ]

write.csv(chlor_sens_anal_change_coef[order(chlor_sens_anal_change_coef$names.fixef_init.), ] , "C:\\Users\\mirac\\Desktop\\chlor_sens_anal_ordered_coefs.csv", row.names=FALSE) 
#write.csv(ref_cond_pred_ordered, "C:\\Users\\mirac\\Desktop\\ref_cond_pred_ordered.csv", row.names=FALSE)

# Variance inflation factor test #####

vif(M_best_subset)

colnames(df_encoded)
#all terms
M_best_subset
all_terms <- log_average_kobs ~ temp_5int + pH_refpka + pH_refpka_squared +
  float_year_scaled + alpha_0 + pH_refpka_x_high_chloride_TRUE +
  pH_refpka_x_temp_5int + balt_class_dsDNA_x_pH_refpka + balt_class_dsRNA_x_pH_refpka +
  balt_class_ssDNA_x_pH_refpka + balt_class_minus_ssRNA_x_pH_refpka +
  DNA_x_pH_refpka + balt_class_dsDNA_x_temp_5int + balt_class_dsRNA_x_temp_5int +
  balt_class_ssDNA_x_temp_5int + balt_class_minus_ssRNA_x_temp_5int +
  DNA_x_temp_5int + balt_class_dsDNA_x_pH_refpka_squared +
  balt_class_dsRNA_x_pH_refpka_squared + balt_class_ssDNA_x_pH_refpka_squared +
  balt_class_minus_ssRNA_x_pH_refpka_squared + DNA_x_pH_refpka_squared +
  buffer + purification_level + high_chloride_TRUE


# Extract the terms from M_best_subset
existing_terms <- attr(terms(M_best_subset), "term.labels")

additional_terms_list <- setdiff(attr(terms(all_terms), "term.labels"), existing_terms)

# Create an empty vector to store VIF values
#vif_values <- numeric(length(additional_terms_list))
#vif_names <- vector(length = length(additional_terms_list))
#vif_names
gvif_values = c()
a_gsif_values = c()
levels_values = c()
vif_names = c()
# Loop through each additional term
#result_vif_df <- data.frame(term = character(), vif = numeric())

# Loop through each additional term
for (i in seq_along(additional_terms_list)) {
  # Create a new formula with the current additional term
  term <-  additional_terms_list[i]
  formula_with_additional_term <- update(M_best_subset, as.formula(paste(". ~ . + ", term)))
  
  # Fit the model
  model <- lmer(formula_with_additional_term, data = df_encoded)
  
  # Calculate VIF
  vif_matrix <- vif(model)
  vif_matrix <- data.frame(column_name = vif_matrix)
  colnames(vif_matrix) <- c('GVIF','Df','a_GSIF')
  gvif_value <- vif_matrix[[term,"GVIF"]]
  a_gsif_value <- vif_matrix[[term,"a_GSIF"]]
  
  
  # Append results to the data frame
  gvif_values = append(gvif_values, gvif_value)
  a_gsif_values = append(a_gsif_values, a_gsif_value)
  levels_values = append(levels_values, vif_matrix[[term,"Df"]])
  vif_names = append(vif_names, additional_terms_list[i])
  #result_vif_df <- rbind(result_vif_df, data.frame(term = additional_terms_list[i], vif = vif_value))
}

# Print or use result_df as needed

vif_results <- data.frame(vif_names, gvif_values,levels_values, a_gsif_values)
vif_results$squared_a_gsif <- vif_results$a_gsif_values^2
vif_results

return_k(virus_name = 'ms2_bacteriophage', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
vif_best_subset <- vif(M_best_subset)
vif_best_subset_df <- data.frame(column_name = vif_best_subset)
vif_best_subset_df[["high_chloride_TRUE","column_name.GVIF"]]


form_test <- update(M_best_subset, as.formula(paste(". ~ . + ", 'balt_class_minus_ssRNA_x_temp_5int')))

# Fit the model
mod_test <- lmer(form_test, data = df_encoded)
vif(mod_test)


vif_matrix_test <- vif(mod_test)
vif_matrix_test <- data.frame(column_name = vif_matrix_test)
colnames(vif_matrix_test) <- c('GVIF','Df','a_GSIF')
vif_matrix_test




#paper section 3.3 predictions under reference conditions

#View(ref_cond_pred_ordered)

#95% CIs
#ms2
ref_cond_pred_ordered[ref_cond_pred_ordered[['Virus Name']] == 'ms2_bacteriophage', ]

#pv1 mahoney
ref_cond_pred_ordered[ref_cond_pred_ordered[['Virus Name']] == 'pv1_mahoney', ]

#hpai 2005
ref_cond_pred_ordered[ref_cond_pred_ordered[['Virus Name']] == 'hpai_a2005', ]

#rrv g3p3
ref_cond_pred_ordered[ref_cond_pred_ordered[['Virus Name']] == 'rrv_g3p3', ]

#find the number of rows whose confidence intervals don't overlap with ms2
lower_bound_threshold <- 7.712198
upper_bound_threshold <- 18.26051
lower_bound_column <- 'Lower Bound of 95% CI'
upper_bound_column <- 'Upper Bound of 95% CI'


nrow(ref_cond_pred_ordered[
  ref_cond_pred_ordered[[lower_bound_column]] < lower_bound_threshold &
    ref_cond_pred_ordered[[upper_bound_column]] > upper_bound_threshold,
])

#cvb5_l070215
ref_cond_pred_ordered[ref_cond_pred_ordered[['Virus Name']] == 'cvb5_l070215', ]

#ev_70
ref_cond_pred_ordered[ref_cond_pred_ordered[['Virus Name']] == 'ev70', ]

#model quantifies role of changing temp

#at 20C, low chloride, MS2
return_k(virus_name = 'ms2_bacteriophage', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)

#at 25C, low chloride, MS2
return_k(virus_name = 'ms2_bacteriophage', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 25, pH = 7.53, chloride = "FALSE", paper_ID = 32)

percent_increase_temp = 15.25525/11.86712
increase_lower_CI = 9.789094 - 7.712198
increase_upper_CI = 23.77366 - 18.26051

percent_increase_temp
increase_lower_CI
increase_upper_CI

(11.86712 + increase_lower_CI)/11.86712
(11.86712 + increase_upper_CI)/11.86712


#model quantifies role of changing pH

#at 20C, low chloride, MS2
return_k(virus_name = 'ms2_bacteriophage', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)

#at 25C, low chloride, MS2
return_k(virus_name = 'ms2_bacteriophage', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 8.53, chloride = "FALSE", paper_ID = 32)

percent_decrease_pH = (11.86712-7.460159)/11.86712
increase_lower_CI = 7.712198 - 4.817685
increase_upper_CI = 18.26051 - 11.55202

percent_decrease_pH
increase_lower_CI
increase_upper_CI

(11.86712 - (11.86712 - increase_lower_CI))/11.86712
(11.86712 - (11.86712 - increase_upper_CI))/11.86712


#pH increase from 7.53 to 10

fold_change_pH <- function(model, virus_name, df_virus_counts, initial_pH, final_pH, temp, chloride){
  k_pH_init <- return_k(virus_name = virus_name, model = model, df_virus_counts = df_virus_counts, temp = temp, pH = initial_pH, chloride = chloride, paper_ID = 32)[[1]]
  k_pH_final <- return_k(virus_name = virus_name, model = model, df_virus_counts = df_virus_counts, temp = temp, pH = final_pH, chloride = chloride, paper_ID = 32)[[1]]
  k_fold_change <- k_pH_init/k_pH_final
  return(print(paste0("When ", virus_name, " changes from " , initial_pH, " to ", final_pH, " inactivation rate decreases " , k_fold_change , " fold.")))
  
}

fold_change_temp <- function(model, virus_name, df_virus_counts, initial_temp, final_temp, pH, chloride){
  k_temp_init <- return_k(virus_name = virus_name, model = model, df_virus_counts = df_virus_counts, temp = initial_temp, pH = pH, chloride = chloride, paper_ID = 32)[[1]]
  k_temp_final <- return_k(virus_name = virus_name, model = model, df_virus_counts = df_virus_counts, temp = final_temp, pH = pH, chloride = chloride, paper_ID = 32)[[1]]
  k_fold_change <- k_temp_final/k_temp_init
  return(print(paste0("When ", virus_name, " changes from " , initial_temp, " to ", final_temp, " inactivation rate increases " , k_fold_change , " fold.")))
  
}


fold_change_pH(model = M_best_subset, virus_name = 'ms2_bacteriophage', df_virus_counts = df_virus_chars_counts, temp = 20, initial_pH = 7.53, final_pH = 10, chloride = "FALSE")
fold_change_pH(model = M_best_subset, virus_name = 'reo_2_jones', df_virus_counts = df_virus_chars_counts, temp = 20, initial_pH = 7.53, final_pH = 10, chloride = "FALSE")
fold_change_pH(model = M_best_subset, virus_name = 'hadv_7a_s1058', df_virus_counts = df_virus_chars_counts, temp = 20, initial_pH = 7.53, final_pH = 10, chloride = "FALSE")

#which virus more resistant under different conditions
return_k(virus_name = 'e1_farouk', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)[[1]]
return_k(virus_name = 'e1_farouk', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 10, chloride = "FALSE", paper_ID = 32)[[1]]

return_k(virus_name = 'hadv_5', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)[[1]]
return_k(virus_name = 'hadv_5', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 10, chloride = "FALSE", paper_ID = 32)[[1]]

#numbers for table 1
under_required_k(model = M_best_subset, temp_input = 10, pH_input = 6, chloride_input = "FALSE", required_CT = 6)[[1]]
under_required_k(model = M_best_subset, temp_input = 10, pH_input = 6, chloride_input = "FALSE", required_CT = 6)[[2]]


CT_regs_6_9 <- c(12,8,6,4,3,2)
CT_regs_10 <- c(90, 60, 45, 30, 22, 15)
required_CT <- c(rep(CT_regs_6_9, times = 4), CT_regs_10)
required_k <-  -log(1/10000)/(required_CT)
pH_input <- c(rep(6, each= 6), rep(7, each= 6), rep(8, each= 6), rep(9, each= 6), rep(10, each= 6))
temp_input <- c(rep(c(.5, 5, 10, 15, 20, 25), times = 5))

req_df <- data.frame(temp_input, pH_input, required_CT, required_k)
sort(unique(df_encoded$virus_name_abbrev))

result_list <- apply(req_df, 1, function(row) {
  under_required_k(row['temp_input'],
                   row['pH_input'],
                   row['required_CT'],
                   model = M_best_subset,
                   chloride_input = 'FALSE')
})

result_df <- do.call(rbind.data.frame, result_list)
colnames(result_df) <- c("num_under_req", "percent_under_req")

# Combining the result with the original DataFrame
req_df <- cbind(req_df, result_df)

#respond to hav comment
sort(unique(df_encoded$virus_name_abbrev))
hav_hm175_ref_k <- return_k(virus_name = 'hav_hm175', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)[[1]]
prediction_hav_comp <- predict_k(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
sum(exponentiate(prediction_hav_comp$est_k) < hav_hm175_ref_k[[1]])

return_k(virus_name = 'hav_hm175', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 6, chloride = "FALSE", paper_ID = 32)[[1]][[1]]
return_k(virus_name = 'hav_hm175', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7, chloride = "FALSE", paper_ID = 32)[[1]]

required_k <- function(required_CT) {
  return(-log(1/10000)/required_CT)
}


required_k <- function(virus_name, req_df, model, df_virus_counts) {
  k_virus  = c()
  percent_more_recalcitrant = c()
  percent_more_recalcitrant_fos = c()
  below_req_CT  = c()
  for (row in 1:nrow(req_df)) {
    k_temp <- return_k(virus_name = virus_name, model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = req_df$temp_input[row], pH = req_df$pH_input[row], chloride = "FALSE", paper_ID = 32)[[1]][[1]]
    required_k <- -log(1/10000)/(req_df$required_CT[row])
    fos_required_k <- 3*required_k
    k_virus  = c(k_virus, k_temp)
    percent_more_recalcitrant_fos = c(percent_more_recalcitrant_fos, fos_required_k/k_temp)
    if(k_temp > required_k) {
      below_req_CT  = c(below_req_CT, "TRUE")
      percent_more_recalcitrant = c(percent_more_recalcitrant, NA)
    }
    else {
      below_req_CT  = c(below_req_CT, "FALSE")
      percent_temp = required_k/k_temp
      percent_more_recalcitrant = c(percent_more_recalcitrant, percent_temp)
    }
  }
  #df <- data.frame(assign(paste0(virus_name, "_k_virus"), k_virus), assign(paste0(virus_name, "_below_req_CT"), below_req_CT))
  df <- data.frame(k_virus, below_req_CT, percent_more_recalcitrant, percent_more_recalcitrant_fos)
  df <- cbind(req_df, df)
  names(df)[names(df)=="k_virus"] <- paste0(virus_name, "_k_virus")
  names(df)[names(df)=="below_req_CT"] <- paste0(virus_name, "_below_req_CT")
  return(df)
}



df_encoded
df_with_below <- required_k('pv1_mk_500', req_df = req_df, model = M_best_subset, df_virus_counts = df_virus_chars_counts)
View(df_with_below)

write.csv(df_with_below, "/Users/mirac/Desktop/Chlorine_Project/chlorine_review/df_with_below_pv1_mk_500.csv")


#View(df_encoded[df_encoded$virus_name_abbrev == "hav_mbb", ])

sort(unique((df_encoded$virus_name_abbrev)))

return_k(virus_name = 'pv1_mahoney', model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 6, chloride = "FALSE", paper_ID = 32)[[1]]
fold_change_pH(model = M_best_subset, virus_name = 'ms2_bacteriophage', df_virus_counts = df_virus_chars_counts, temp = 20, initial_pH = 6, final_pH = 9, chloride = "FALSE")
fold_change_temp(model = M_best_subset, virus_name = 'ms2_bacteriophage', df_virus_counts = df_virus_chars_counts, pH = 7.53, initial_temp = .5, final_temp = 25, chloride = "FALSE")

####more figures
#View(df[df$balt_class=='ssDNA',])
# Plot Random Effects ####

## for dotplot, qqmath, default ranef plot:
#str(rr1 <- ranef(M_final))
#dotplot(rr1, condVar = T)  ## default

randoms<-ranef(M_best_subset, condVar = TRUE)
qq <- attr(ranef(M_best_subset, condVar = TRUE)[[1]], "postVar")
randoms_df = randoms[[1]]
rand.interc<-randoms_df

randoms_df[[randoms_df<.0005]]
#find intercepts

sd(randoms_df$`(Intercept)`)
# a hypothetical paper drawn from the population has a 32% chance of having a rate constant be biased by 
#at least ±2.8 folds (= 10^.423)

paper_bias_1sd <- 10^sd(randoms_df$`(Intercept)`)
paper_bias_1sd

#find intraclass correlation coefficient (ICC)

ICC_M_best <- 0.4493/(0.4493 + 0.2922)
ICC_M_best

M_best_subset

df_ranef<-data.frame(Intercepts=randoms_df[ ,1],
                     sd.interc=2*sqrt(qq[,,1:length(qq)]),
                     lev.names=rownames(rand.interc))

df_ranef$lev.names<-factor(df_ranef$lev.names,levels=df_ranef$lev.names[order(df_ranef$Intercepts)])
ranef_plot <- ggplot(df_ranef, aes(x = lev.names,y = Intercepts)) + 
  #Added horizontal line at y=0, error bars to points and points with size two
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0,color="black") + 
  geom_point(size = 2, colour = "deepskyblue4") +
  #Removed legends
  theme(legend.position="none")+
  #Changed appearance of plot (black and white theme) and x and y axis labels
  theme_bw() + 
  xlab("") + 
  ylab("Random Intercept (log10 scale)")+
  #Final adjustments of plot
  theme(#axis.text.x=element_text(size=rel(1.2)),
    #axis.title.x=element_text(size=rel(1.3)),
    axis.text.y= element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor=element_blank(),
    # panel.grid.major.x=element_blank(),
    panel.grid.major.x = element_line(size = 0.5))+
  #To put levels on y axis you just need to use coord_flip()
  scale_y_continuous(breaks=seq(-2,1.5,.5), minor_breaks = seq(-2,1.5,.5))+
  coord_flip()

ranef_plot
jpeg(file="random_effects.jpeg",width=6,height=6,units="in", res=300)
ranef_plot
dev.off()



#Violin Plot of Predictions ####
all_preds <- return_k_ordered(model = M_best_subset, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
all_preds

plot_violin <- ggplot(all_preds, aes(x=balt_class, y=est_k)) + 
  geom_violin() +
  xlab("Baltimore class") +
  ylab("Predicted inactivation rate constants (L/mg*min)") +
  scale_y_continuous(trans='log10')

jpeg(file="revisions_violin_plot.jpeg",width=6,height=6,units="in", res=300)
plot_violin
dev.off()

plot_violin


fixef(M_final)

0.20098/(0.20098+0.08611)

sqrt(0.086)


library(AICcmodavg)

#define list of models
models <- list(M_final, M1a)

#specify model names
mod.names <- c('M_final','M1a')

#calculate AIC of each model
aictab(cand.set = models, modnames = mod.names)

AICc(M_final)
extractAIC(M_final)
library(stats)
AIC(M_final)
tab_model(M_final, digits.re = 3, digits = 4)
coef(M_final)

print(VarCorr(M_final),comp="Variance")

# [ PLOT ] Facet plot of each balt_class ####################################################################################################

library('stats')
library('forcats')

residuals_facet <- function(model, df, slope){
  dev.off()
  df_result_new <- data.frame(Predicted = predict(model), Observed = df$log_average_kobs)
  df_result_new$Residuals <-df_result_new$Observed - df_result_new$Predicted
  df_result_new$buffer_type = c(df$buffer_type)
  df_result_new$high_chloride = c(df$high_chloride)
  df_result_new$balt_class = c(df$balt_class)
  df_result_new$pH = c(df$pH_refpka)
  df_result_new$Temperature = c(df$temp_5int)
  
  p_new <- ggplot(df_result_new, 
                  aes(x = Observed,
                      y = Predicted)) +
    geom_point(aes(color = balt_class)) + 
    geom_abline(intercept = 0,
                slope = slope,
                color = "red",
                linewidth = .5) +
    scale_color_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values=c("#75bbfd", "#fcb001", "#02ab2e",  "#7e1e9c", "#c44240")) +
    theme(legend.position="none")+
    #theme(legend.position = c(.75,.75) , legend.text = element_text(size=14), legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    ggtitle("")
  
  p_new + facet_wrap( ~ fct_relevel(balt_class, c("+ssRNA", "-ssRNA", "dsRNA", "dsDNA", "ssDNA")), as.table = FALSE)
}

pre_interactions_residuals <- lmer(log_average_kobs ~ virus_name_abbrev + temp_5int + pH_refpka*high_chloride + pH_refpka:temp_5int + (1|paper_ID), data = df_encoded, REML = FALSE)

residuals_facet(model = pre_interactions_residuals, df = df_encoded, slope = 1)

## testing


residuals_facet(model = M_best_subset, df = df_encoded)

residuals_facet(model = pre_interactions_residuals, df = df_encoded, x_axis = temp_5int)

residuals.merMod(M_final)

data_for_res_qq <- data.frame(Predicted = predict(M_best_subset), Observed = df_encoded$log_average_kobs)
data_for_res_qq$residuals <- data_for_res_qq$Predicted - data_for_res_qq$Observed

#Use the conditional residuals to check the normality of the error term in the model.


ggplot(data_for_res_qq, aes(sample = residuals)) +  # Create QQplot with ggplot2 package
  stat_qq() +
  stat_qq_line(col = "red") +
  ggtitle('qqplot of residuals for normal distribution') +
  theme(plot.title = element_text(hjust = 0.5))


data_for_res_qq_no_ranef <- data.frame(Predicted = predict(M_best_subset, re.form = NA), Observed = df_encoded$log_average_kobs)
data_for_res_qq_no_ranef$residuals <- data_for_res_qq_no_ranef$Predicted - data_for_res_qq_no_ranef$Observed

#Use the conditional residuals to check the normality of the error term in the model.


qq_unc <- ggplot(data_for_res_qq_no_ranef, aes(sample = residuals)) +  # Create QQplot with ggplot2 package
  stat_qq() +
  stat_qq_line(col = "red") +
  ggtitle('qqplot of unconditional residuals for normal distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Theoretical Quantiles") + 
  ylab("Sample Quantiles")



#View(data_for_res_qq)

qq_obs <- ggplot(data_for_res_qq, aes(sample = Observed)) +  # Create QQplot with ggplot2 package
  stat_qq() +
  stat_qq_line(col = "red") +
  ggtitle('qqplot of observed values for normal distribution') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Theoretical Quantiles") + 
  ylab("Sample Quantiles")

qq_obs + qq_unc
qqnorm(residuals(M_final, type = c("pearson")), main='Normal')
qqline(residuals(M_final, type = c("pearson")))


View(df)
write.csv(df, "/Users/mirac/Desktop/Chlorine_Project/chlorine_review/dataset_for_vis.csv")


#### plot model coefficients ################

summary_table = summ(M_best_subset, confint = TRUE, digits = 5)
est_virus = summary_table$coeftable[,c(1,2,3)]
model_coefs_and_preds = data.frame((est_virus))
#View(model_coefs_and_preds)
new_row_names <- gsub("virus_name_abbrev", "", rownames(model_coefs_and_preds))
rownames(model_coefs_and_preds) <- new_row_names

names(model_coefs_and_preds) <- c("estimate", "lower_CI", "upper_CI")
model_coefs_and_preds = model_coefs_and_preds %>% arrange(desc(estimate))
model_coefs_and_preds <- tibble::rownames_to_column(model_coefs_and_preds, "coef_names")
#model_coefs_and_preds$estimate_nonlong <- 10 ^ (model_coefs_and_preds$estimate)
#model_coefs_and_preds$lowerCI_nonlong <- 10 ^ (model_coefs_and_preds$lower_CI)
#model_coefs_and_preds$upperCI_nonlong <- 10 ^ (model_coefs_and_preds$upper_CI)

df_virus_chars_counts[df_virus_chars_counts$n > 4, ]$virus_name_abbrev

#plot_model(M_best_subset)

experimental_terms <- c("temp_5int", "pH_refpka", "pH_refpka_squared" , "pH_refpka_x_high_chloride_TRUE" , "pH_refpka_x_temp_5int", 
                        "balt_class_dsRNA_x_pH_refpka", "DNA_x_pH_refpka", "balt_class_dsDNA_x_temp_5int", 
                        "balt_class_ssDNA_x_temp_5int", "DNA_x_pH_refpka_squared", "high_chloride_TRUE")

#ms2_bacteriophage absent

virus_terms_list <- c("cvb1_49683" , "cvb1_conn5","cvb3_nancy","cvb5","cvb5_ea_80",
                      "cvb5_faulkner","e1_farouk","e11_gregory", "f2_bacteriophage","fcv_f9",
                      "hadv_2_strain_6","hadv_5","hav_hm175","hrv_wa","mnv_1",
                      "parvovirus_h1","pr772_bacteriophage","prd1_bacteriophage","pv1_brunhilde",
                      "pv1_mahoney","pv1_mk_500","pv1_sabin","pv2_mef1","pv2_sabin",
                      "pv3_saukett","qβ_bacteriophage","srv_sa11","φx174_bacteriophage")

virus_terms <- model_coefs_and_preds[model_coefs_and_preds$coef_names %in% virus_terms_list, ]
fixed_terms<- model_coefs_and_preds[model_coefs_and_preds$coef_names %in% experimental_terms, ]

exp_plot <- ggplot(fixed_terms, aes(x = reorder(coef_names, estimate), y = estimate )) + 
  #Added horizontal line at y=0, error bars to points and points with size two
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=0,color="black") + 
  geom_point(size = 2, colour = "deepskyblue4") +
  #Removed legends
  theme(legend.position="none")+
  #Changed appearance of plot (black and white theme) and x and y axis labels
  #theme_bw() + 
  xlab("") + 
  ylab("Coefficient Estimates")+
  scale_y_continuous(limits = c(-1,1))+
  #Final adjustments of plot
  theme(#axis.text.x=element_text(size=rel(1.2)),
    panel.grid.minor=element_blank(),
    panel.grid.major.x = element_line(size = 0.5))+
  coord_flip()+  # Example geom, replace with your actual geoms
  geom_hline(yintercept = 0, color = "black")

vir_plot <- ggplot(virus_terms, aes(x = reorder(coef_names, estimate), y = estimate )) + 
  #Added horizontal line at y=0, error bars to points and points with size two
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=0,color="black") + 
  geom_point(size = 2, colour = "deepskyblue4") +
  #Removed legends
  theme(legend.position="none")+
  #Changed appearance of plot (black and white theme) and x and y axis labels
  #theme_bw() + 
  xlab("") + 
  ylab("Coefficient Estimates")+
  scale_y_continuous(limits = c(-1.5,3))+
  #Final adjustments of plot
  theme(#axis.text.x=element_text(size=rel(1.2)),
    panel.grid.minor=element_blank(),
    panel.grid.major.x = element_line(size = 0.5))+
  coord_flip()+  # Example geom, replace with your actual geoms
  geom_hline(yintercept = 0, color = "black")
#scale_y_continuous(breaks=seq(-2,1.5,.5), minor_breaks = seq(-2,1.5,.5))+



design <- "
  1111
  1111
  2222
"
vir_plot + exp_plot + plot_layout(design = design)


tab_model(M_best_subset)

df_encoded
groupby