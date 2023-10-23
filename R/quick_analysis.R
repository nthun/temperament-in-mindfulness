# cond – condition (control/mindfulness)
# 
# CBs_ch – Corsi backward span 
# HFmix_ch – Hearts and Flowers vegyes részének hibáinak különbsége (post – pre)
# HFflowch - Hearts and Flowers virágok részének hibáinak különbsége (post – pre)
# Shark_ch – Go/No-Go shark error (post – pre)
# TMTe_ch Trail Making Test error (post – pre)
# TMTt_ch – Trail Making Test time (post – pre)
# Loc_e_ch - Location error (post – pre)
# Dir_e_ch – Direction error (post – pre)
# SSTneuch – SST neutral (post – pre)
# Percomch – SST omissions neutral percent (post – pre)
# Cort_ch – cortisol (post – pre)
# 
# 
# 
# EAS – temperamentum kérdoív: (ahol van a változó végén egy 2-es pl.: EAS_shy2 ott ki van véve egy tétel a jobb reliabilitás miatt)
# EAS_emot – emocionalitás skála
# EAS_acti – aktivitás skála
# EAS_soc – szociablilitás 
# EAS_shy – félénség 
# SDQ – viselkedésszabályozás kérdoív
# SDQ_tot – összprobléma pontszám
# SDQ_erz – érzelmi tünetek
# SDQ_vis – viselekdési problémák
# SDQ_hip – hiperaktivitás
# SDQ_kapcs – kortárskapcsolati problémák
# SDQ_szoc – proszociális skála ( ez az egy nem problémát jelöl, minél több, annál jobb)

library(tidyverse)
library(haven)
library(lmerTest)
library(broom)
library(broom.mixed)
library(gt)
library(correlation)

theme_set(theme_light())

all_raw <- 
  map(list.files("data/", full.names = TRUE), read_spss) |> 
  set_names(1:4)

all_proc <-
  all_raw |> 
  map(~mutate(.x, across(where(is.labelled), as_factor))) |> 
  map(~mutate(.x, across(any_of(c("ID", "id")), as.character))) |> 
  map(~rename(.x, cond = any_of(c("Cond", "cond")))) |> 
  map(~mutate(.x, cond = case_when(cond %in% c("kontroll", "control") ~ "control",
                                   cond %in% c("mindfulness", "intervention") ~ "mindfulness"
                                   ))) |> 
  map2(.y = names(all_raw), ~mutate(.x, id = paste0(.y, "_", row_number()), .before = 1))
  

map(all_proc, glimpse)

map(all_proc, ~pull(.x, cond))

map(all_proc, ~select(.x, contains("cort")))

# A modelleket úgy építettem fel, hogy a moderátor az SDQ valamelyik faktora/összpontszám vagy az EAS valamelyik faktora, a függetletlen változó mindig a kondíció, és a függo változó a kortizolban, végrehajtókban bekövetkezett változás. Mindegyik kivont értéknél (change) az a jobb, ha minél kisebb, kivéve a Corsi backward span esetében (CBs_ch). 

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
  pivot_longer(any_of(outcomes), names_to = "outcome", values_to = "outcome_value", values_drop_na = TRUE) |> 
  pivot_longer(any_of(moderators), names_to = "moderator", values_to = "moderator_value", values_drop_na = TRUE) |> 
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
         plot = map(data, ~ggplot() +
                           aes(x = moderator_value, y = outcome_value, color = cond) +
                           geom_point() +
                           geom_smooth(method = "lm")))

temp <- 
  all_models |> 
  slice(31) |> 
  pull(data)%>%
  .[[1]]

temp |> 
  mutate(mod_group = case_when(moderator_value >= 1 ~ "high",
                               moderator_value <= -1 ~ "low")) |> 
  ggplot() +
  aes(x = moderator_value, y = outcome_value, color = cond) +
  geom_point() +
  geom_smooth(method = "lm")


all_models |> 
  unnest(estimates) |> 
  filter(effect == "fixed") |> 
  select(outcome, moderator, study_n, term, estimate, conf.low, conf.high, std.error:p.value) |> 
  # group_by(outcome, moderator) |> 
  gt() |> 
  fmt_number(c("estimate", "std.error", "conf.low", "conf.high", "statistic", "df"), decimals = 2) |> 
  fmt_number(c("p.value"), decimals = 3) |> 
  tab_options(column_labels.font.weight = "bold",
              column_labels.background.color = "lightgrey",
              row_group.background.color = "lightgrey", 
              row_group.font.weight = "bold", 
              row.striping.include_table_body = TRUE,
              row.striping.background_color = "#EEEEEE"
  )# |> 
  gtsave("docs/models.docx")
  
plots <- 
  all_models |> 
  filter((outcome == "Cort_ch" & moderator %in% c("EAS_emot", "SDQ_tot")) | 
          (outcome == "HFmix_ch" & moderator == "EAS_emot") |
           (outcome == "Shark_ch" & moderator == "SDQ_vis"))
  
plt <- pull(plots, plot)


plots |> 
  slice(2) |> 
  pull(plot) %>%
  .[[1]]
  
xy |>   
  filter(outcome == "TMTe_ch", moderator == "SDQ_tot") |> 
  pull(data) %>% view()



# A 4 studyban összesen 243 fo 1 IV, 10 DV, 12 moderator. Ez összesen 89 kombinációt jelent.
# A modell minden esetben így néz ki: outcome ~ cond * moderator + (1|study)
# Ez azt jelenti, hogy egy modelen belül figyelembe vesszük, hogy más studyból származik az adat.
# Összesen 4 kombinációra lehet mind a 4 studyból adatot használni, szintén 4 studyból 3-ra, 29 studyból 2 kombinációra, 52 studyból csak 1
# Csak arra lehet többszintu regressziót használni, amihez 1-nél több study van, tehát 37 kombinációra.
  

# Correlations ------------------------------------------------------------
# Megnezned lszi a korrelaciokat a moderatorok (1. temperamentum (EAS) skalak, 2. SDQ total es skalak) es a baseline (pre-test) kortizol es kognitiv teljesitmeny (1. CB, 2. HF mix, 3. HF flower, 4. shark) kozott?

correlators <- c("EAS_emot", "EAS_acti", "EAS_soc", "EAS_shy", 
                 "SDQ", "SDQ_tot", "SDQ_erz", "SDQ_vis", "SDQ_hip", "SDQ_hip", "SDQ_kapcs", "SDQ_szoc", 
                 "pre_CB_span", "pre_G_comission_error_shark_error",  "pre_HF_flowers_error", "pre_HF_mix_error", 
                 "pre_cort_basal_mean")

all_proc |> 
  map(all_proc, ~select(.x, any_of(c("id", correlators)))) |> 
  bind_rows(.id = "study") 


  