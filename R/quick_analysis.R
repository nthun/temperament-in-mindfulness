# cond � condition (control/mindfulness)
# 
# CBs_ch � Corsi backward span 
# HFmix_ch � Hearts and Flowers vegyes r�sz�nek hib�inak k�l�nbs�ge (post � pre)
# HFflowch - Hearts and Flowers vir�gok r�sz�nek hib�inak k�l�nbs�ge (post � pre)
# Shark_ch � Go/No-Go shark error (post � pre)
# TMTe_ch Trail Making Test error (post � pre)
# TMTt_ch � Trail Making Test time (post � pre)
# Loc_e_ch - Location error (post � pre)
# Dir_e_ch � Direction error (post � pre)
# SSTneuch � SST neutral (post � pre)
# Percomch � SST omissions neutral percent (post � pre)
# Cort_ch � cortisol (post � pre)
# 
# 
# 
# EAS � temperamentum k�rdo�v: (ahol van a v�ltoz� v�g�n egy 2-es pl.: EAS_shy2 ott ki van v�ve egy t�tel a jobb reliabilit�s miatt)
# EAS_emot � emocionalit�s sk�la
# EAS_acti � aktivit�s sk�la
# EAS_soc � szociablilit�s 
# EAS_shy � f�l�ns�g 
# SDQ � viselked�sszab�lyoz�s k�rdo�v
# SDQ_tot � �sszprobl�ma pontsz�m
# SDQ_erz � �rzelmi t�netek
# SDQ_vis � viselekd�si probl�m�k
# SDQ_hip � hiperaktivit�s
# SDQ_kapcs � kort�rskapcsolati probl�m�k
# SDQ_szoc � proszoci�lis sk�la ( ez az egy nem probl�m�t jel�l, min�l t�bb, ann�l jobb)

library(tidyverse)
library(haven)
library(lmerTest)
library(broom)
library(broom.mixed)
library(gt)
library(correlation)

theme_set(theme_light())

all_raw <- 
  map(list.files(path = "data/", pattern = ".*.sav",full.names = TRUE), read_spss) |> 
  set_names(1:4)


# Data processing ---------------------------------------------------------
all_proc <-
  all_raw |> 
  map(~mutate(.x, across(where(is.labelled), as_factor))) |> 
  map(~mutate(.x, across(any_of(c("ID", "id")), as.character))) |> 
  map(~rename(.x, cond = any_of(c("Cond", "cond")))) |>
  map(~rename(.x, SDQ_kapcs = any_of(c("SDQ_kapcs", "SDQkapcs")))) |> 
  map(~rename(.x, SDQ_szoc = any_of(c("SDQ_szoc", "SDQszoc"))),
      ~rename(.x, HFflowch = any_of(c("HFflow_ch", "HFflowch")))) |> 
    map(~mutate(.x, cond = case_when(cond %in% c("kontroll", "control") ~ "control",
                                   cond %in% c("mindfulness", "intervention") ~ "mindfulness"
                                   ))) |> 
  map2(.y = names(all_raw), ~mutate(.x, id = paste0(.y, "_", row_number()), .before = 1))
  
# This is a super ugly solution to transform cortisol values in study 3 from �g/dL to ng/mL
all_proc$`3` <- 
  all_proc$`3` |> 
  mutate(across(contains("cort"), ~`*`(., 10)))


map(all_proc, glimpse)

map(all_proc, ~pull(.x, cond))

map(all_proc, ~select(.x, contains("cort")))



# Multilevel models -------------------------------------------------------
# A modelleket �gy �p�tettem fel, hogy a moder�tor az SDQ valamelyik faktora/�sszpontsz�m vagy az EAS valamelyik faktora, a f�ggetletlen v�ltoz� mindig a kond�ci�, �s a f�ggo v�ltoz� a kortizolban, v�grehajt�kban bek�vetkezett v�ltoz�s. Mindegyik kivont �rt�kn�l (change) az a jobb, ha min�l kisebb, kiv�ve a Corsi backward span eset�ben (CBs_ch). 

outcomes <- c("Cort_ch", "CBs_ch", "HFmix_ch", "HFflowch", "Shark_ch", "TMTe_ch", "Loc_e_ch", "Dir_e_ch", "SSTneuch", "Percomch")

moderators <- c("EAS_emot", "EAS_acti", "EAS_soc", "EAS_shy", "SDQ", "SDQ_tot", "SDQ_erz", "SDQ_vis", "SDQ_hip", "SDQ_hip", "SDQ_kapcs", "SDQ_szoc")

# Create aggregated dataset
mod_df <- map(all_proc, ~select(.x, any_of(c("id", "cond", outcomes, moderators))))

# Create models

all_models <-
  mod_df |> 
  bind_rows(.id = "study") |> 
  # Standardize all numeric variables
  mutate(across(where(is.numeric), ~scale(.) |> as.numeric())) |> 
  # Create separate datasets for each model
  pivot_longer(any_of(outcomes), 
               names_to = "outcome", values_to = "outcome_value", values_drop_na = TRUE) |> 
  pivot_longer(any_of(moderators), 
               names_to = "moderator", values_to = "moderator_value", values_drop_na = TRUE) |> 
  group_by(outcome, moderator) |> 
  mutate(study_n = n_distinct(study)) |> 
  group_by(outcome, moderator, study_n) |> 
  nest() |> 
  ungroup() |>
  # Remove all models with less than 2 datasets
  filter(study_n > 1) |> 
  mutate(model = 
           if_else(study_n > 1, 
                   map(data, ~lmer(outcome_value ~ cond * moderator_value + (1|study), data = .x)),
                   map(data, ~lm(outcome_value ~ cond * moderator_value, data = .x))),
         estimates = map(model, tidy, conf.int = TRUE),
         performance = map(model, glance),
         plot = pmap(list(data, outcome, moderator), 
                          ~..1 |> 
                           ggplot() +
                           aes(x = moderator_value, y = outcome_value, color = cond) +
                           geom_point() +
                           geom_smooth(method = "lm") +
                           labs(y = first(..2), x = first(..3),
                           color = "Condition")
         ))


all_models |> 
  unnest(estimates) |> 
  filter(effect == "fixed") |> 
  select(outcome, moderator, study_n, term, estimate, conf.low, conf.high, std.error:p.value) |> 
  # arrange(p.value) |> 
  group_by(outcome, moderator) |>
  gt() |> 
  fmt_number(c("estimate", "std.error", "conf.low", "conf.high", "statistic", "df"), decimals = 2) |> 
  fmt_number(c("p.value"), decimals = 3) |> 
  tab_options(column_labels.font.weight = "bold",
              column_labels.background.color = "lightgrey",
              row_group.background.color = "lightgrey", 
              row_group.font.weight = "bold", 
              row.striping.include_table_body = TRUE,
              row.striping.background_color = "#EEEEEE"
  ) # |> 
  # gtsave("docs/models.docx")
  
plots <- 
  all_models |> 
  filter((outcome == "Cort_ch" & moderator %in% c("EAS_emot", "SDQ_tot")) | 
          (outcome == "HFmix_ch" & moderator == "EAS_emot") |
           (outcome == "Shark_ch" & moderator == "SDQ_vis"))


plt <- pull(plots, plot)

plots |> 
  unnest(performance)
  pull(plot)

# plots$plot[[4]]

# A 4 studyban �sszesen 243 fo 1 IV, 10 DV, 12 moderator. Ez �sszesen 89 kombin�ci�t jelent.
# A modell minden esetben �gy n�z ki: outcome ~ cond * moderator + (1|study)
# Ez azt jelenti, hogy egy modelen bel�l figyelembe vessz�k, hogy m�s studyb�l sz�rmazik az adat.
# �sszesen 4 kombin�ci�ra lehet mind a 4 studyb�l adatot haszn�lni, szint�n 4 studyb�l 3-ra, 29 studyb�l 2 kombin�ci�ra, 52 studyb�l csak 1
# Csak arra lehet t�bbszintu regresszi�t haszn�lni, amihez 1-n�l t�bb study van, teh�t 37 kombin�ci�ra.
  

# Correlations ------------------------------------------------------------
# Megnezned lszi a korrelaciokat a moderatorok (1. temperamentum (EAS) skalak, 2. SDQ total es skalak) es a baseline (pre-test) kortizol es kognitiv teljesitmeny (1. CB, 2. HF mix, 3. HF flower, 4. shark) kozott?

# HFflowch - Cort_ch
# Shark_ch - Cort_ch
# Shark_ch - HFflowch
# Cort_ch - EAS_emot
# HFflowch - EAS_emot
# Shark_ch - EAS_emot
# Cort_ch - SDQ_tot
# HFflowch - SDQ_tot
# Shark_ch - SDQ_tot
# Cort_ch - SDQ_vis
# HFflowch - SDQ_vis
# Shark_ch - SDQ_vis

correlators <- c("EAS_emot", "EAS_acti", "EAS_soc", "EAS_shy", 
                 "SDQ", "SDQ_tot", "SDQ_erz", "SDQ_vis", "SDQ_hip", "SDQ_hip", 
                 "SDQ_kapcs", "SDQ_szoc", "pre_CB_span", 
                 "pre_G_comission_error_shark_error",  "pre_HF_flowers_error",
                 "pre_HF_mix_error", "pre_cort_basal_mean",
                 "Cort_ch", "HFflowch", "Shark_ch")

all_proc |> 
  map(~select(.x, any_of(c("id", correlators)))) |> 
  bind_rows(.id = "study") |> 
  correlation(p_adjust = "none", method = "spearman") |> 
  arrange(p) |> 
  drop_na(rho) |> 
  gt() |> 
  cols_hide(c("CI", "Method", "S")) |> 
  cols_move(n_Obs, CI_high) |> 
  fmt_number(c("rho", "CI_low",	"CI_high", "S"), decimals = 2) |> 
  fmt_number("p", decimals = 3) |> 
  tab_options(column_labels.font.weight = "bold",
              column_labels.background.color = "lightgrey", 
              row.striping.background_color = "#EEEEEE", 
              row.striping.include_table_body = TRUE) |> 
  gtsave("docs/correlations.docx") |>
  force()

all_proc |> 
  map(~select(.x, any_of(c("id", correlators)))) |> 
  bind_rows(.id = "study") |> 
  correlation(p_adjust = "none", method = "spearman", redundant = FALSE) |> 
  summary() |> 
  gt() |> 
  fmt_number(-Parameter, decimals = 2) |> 
  sub_missing(missing_text = "") |> 
  gtsave("docs/cor_matrix.docx")


  

# SANDBOX -------------------------------------------------------------------------------------
