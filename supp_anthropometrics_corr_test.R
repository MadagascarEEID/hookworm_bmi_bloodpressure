library(corrplot)
library(tidyverse)
long_survey_data <- read.csv("/Users/levkolinski/Library/CloudStorage/Box-Box/NSF Data/Survey Data/long_survey.csv")

long_survey_data$waist_circumference <- long_survey_data$waist_circumference/100
long_survey_data$hip_circumference <- long_survey_data$hip_circumference/100


long_survey_health_data <- long_survey_data |>
  dplyr::select(c(gender, age, height, weight,
                  body_fat, waist_circumference, hip_circumference)) |>
  mutate(across(contains("height"), as.numeric)) |>
  mutate(across(contains("weight"), as.numeric)) |>
  mutate(age = as.numeric(age)) |>
  mutate(body_fat = as.numeric(body_fat)) |> 
  mutate(BMI = weight / (height)^2) |> 
  filter(BMI > 10 & BMI <60) |> 
  mutate(waist_height_ratio = waist_circumference/height,
         waist_hip_ratio = waist_circumference/hip_circumference)

# hist(long_survey_health_data$BMI) |> 
# hist(long_survey_health_data$waist_height_ratio)
# hist(long_survey_health_data$body_fat)


anthropometric_bothsex <- long_survey_health_data |> 
  select(BMI, waist_height_ratio, body_fat, waist_hip_ratio) |> 
  rename(
    "BMI" = BMI,
    "Waist-Height Ratio" =waist_height_ratio,
    "% Body Fat (Bioimpedance)" = body_fat,
    "Waist-Hip Ratio" = waist_hip_ratio) |> 
  drop_na() |> 
  as.matrix()

M_both_sex = cor(anthropometric_bothsex,
                 method = "pearson")

testRes = cor.mtest(anthropometric_bothsex,
                    method = "pearson", 
                    conf.level = 0.95)


cor.test(long_survey_health_data$BMI, long_survey_health_data$body_fat)
M_both_sex
testRes


jpeg("/Users/levkolinski/Desktop/hookworm_bmi_bloodpressure/Figures/R1/anthro_plot_supp.jpeg",
     units = "in", width = 7, height = 7, res = 350, bg = "white")

corrplot(M_both_sex, p.mat = testRes$p,
         order = 'alphabet',
         diag = TRUE,
         type = "upper",
         sig.level = c(0.001, 0.01, 0.05),
         col.lim = c(0,1),
         insig = 'label_sig',
         pch.cex = 0.9,
         tl.col = "black")

dev.off()
