## BEFORE RUNNING THIS SCRIPT, RUN THE "BMI - Model Averaging (Main)" 
## SECTION OF THE RMD (THROUGH LINE 562)

library(dplyr)
library(ggplot2)
library(purrr)
library(marginaleffects)
# 
bmi_model_1 <- lm(
  log_BMI ~ village + age + age2 + gender_f + wealth_index +
    bin_pos_necator + bin_pos_ancylostoma + bin_pos_ancylostoma:bin_pos_necator +
    bin_pos_necator:wealth_index + bin_pos_ancylostoma:wealth_index +
    bin_pos_ancylostoma:bin_pos_necator:wealth_index,
  data = merged_df_for_model, na.action = na.fail
)
# 
dredged_bmi_model_1  <- dredge(bmi_model_1, subset = (!age2 | age) &
                                 (!`bin_pos_ancylostoma:bin_pos_necator` | (bin_pos_necator & bin_pos_ancylostoma)) &
                                 (!`bin_pos_necator:wealth_index` | (bin_pos_necator & wealth_index)) &
                                 (!`bin_pos_ancylostoma:bin_pos_necator:wealth_index` | (bin_pos_necator & bin_pos_ancylostoma & wealth_index)) &
                                 (!`bin_pos_ancylostoma:wealth_index` | (bin_pos_ancylostoma & wealth_index)))

averaged_bmi_model_1 <- model.avg(dredged_bmi_model_1, subset = delta < 2)

# BMI MARGINAL EFFECTS ----

model_list <- MuMIn::get.models(dredged_bmi_model_1, subset = delta < 2)
weights    <- MuMIn::Weights(dredged_bmi_model_1)[seq_along(model_list)]

# ── Wealth grid (use observed range, finely spaced) ─────────────────────────
wealth_range <- range(merged_df_for_model$wealth_index, na.rm = TRUE)
wealth_grid  <- seq(wealth_range[1], wealth_range[2], length.out = 50)

# 4 infection-status combinations
infection_grid <- expand.grid(
  bin_pos_necator     = c(0, 1),
  bin_pos_ancylostoma = c(0, 1)
)

# ── Population-averaged predictions per model, across wealth grid × infection group ──
get_predictions_for_model <- function(model, wealth_grid, infection_grid) {
  map_dfr(seq_len(nrow(infection_grid)), function(i) {
    nec <- infection_grid$bin_pos_necator[i]
    anc <- infection_grid$bin_pos_ancylostoma[i]
    
    preds <- avg_predictions(
      model,
      variables = list(wealth_index = wealth_grid),
      newdata = datagrid(
        model    = model,
        wealth_index = wealth_grid,
        bin_pos_necator     = nec,
        bin_pos_ancylostoma = anc,
        grid_type = "counterfactual",
        newdata   = merged_df_for_model
      ),
      by = "wealth_index"
    )
    
    preds$bin_pos_necator     <- nec
    preds$bin_pos_ancylostoma <- anc
    preds
  })
}

pred_list <- map(model_list, ~ get_predictions_for_model(.x, wealth_grid, infection_grid))

# ── Weighted pooling across models (estimates stay on log_BMI scale) ───────
pooled_preds <- pred_list |>
  imap_dfr(~ mutate(.x, w = weights[.y])) |>
  group_by(wealth_index, bin_pos_necator, bin_pos_ancylostoma) |>
  summarise(
    estimate  = sum(estimate * w) / sum(w),
    std.error = sqrt(sum(w^2 * std.error^2) / sum(w)^2),
    .groups   = "drop"
  ) |>
  mutate(
    conf.low  = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error,
    group = case_when(
      bin_pos_ancylostoma == 1 & bin_pos_necator == 1 ~ "Coinfected",
      bin_pos_ancylostoma == 1                        ~ "A. ceylanicum",
      bin_pos_necator     == 1                        ~ "N. americanus",
      TRUE                                            ~ "Uninfected"
    )
  )

pooled_preds$group <- factor(pooled_preds$group,
                             levels = c("Uninfected", "A. ceylanicum", "N. americanus", "Coinfected"))

# ── Restrict each group's line to its observed wealth range ────────────────
wealth_ranges_by_group <- merged_df_for_model |>
  mutate(
    group = case_when(
      bin_pos_ancylostoma == 1 & bin_pos_necator == 1 ~ "Coinfected",
      bin_pos_ancylostoma == 1                        ~ "A. ceylanicum",
      bin_pos_necator     == 1                        ~ "N. americanus",
      TRUE                                            ~ "Uninfected"
    )
  ) |>
  group_by(group) |>
  summarise(
    wealth_min = min(wealth_index, na.rm = TRUE),
    wealth_max = max(wealth_index, na.rm = TRUE),
    .groups = "drop"
  )

pooled_preds <- pooled_preds |>
  left_join(wealth_ranges_by_group, by = "group") |>
  mutate(
    in_range = wealth_index >= wealth_min & wealth_index <= wealth_max,
    linetype_grp = ifelse(in_range, "observed", "extrapolated")
  )

# ── Observed points (log_BMI scale), grouped by infection status ───────────
observed_df <- merged_df_for_model |>
  mutate(
    group = case_when(
      bin_pos_ancylostoma == 1 & bin_pos_necator == 1 ~ "Coinfected",
      bin_pos_ancylostoma == 1                        ~ "A. ceylanicum",
      bin_pos_necator     == 1                        ~ "N. americanus",
      TRUE                                            ~ "Uninfected"
    )
  ) |>
  mutate(group = factor(group))

# ── Plot, faceted by infection group ────────────────────────────────────────
## necator and ancylostoma colors match SEMs
group_colors <- c(
  "Uninfected"    = "#999999",
  "N. americanus" = "#CC79A7",
  "A. ceylanicum" = "#009E73",
  "Coinfected"    = "#000000"
)

ggplot() +
  geom_point(
    data = observed_df,
    aes(x = wealth_index, y = log_BMI, color = group),
    alpha = 0.25, size = 1.3
  ) +
  geom_ribbon(
    data = pooled_preds,
    aes(x = wealth_index, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15, color = NA
  ) +
  geom_line(
    data = pooled_preds,
    aes(x = wealth_index, y = estimate, color = group, linetype = linetype_grp),
    linewidth = 1
  ) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_linetype_manual(
    values = c(observed = "solid", extrapolated = "dashed"),
    guide = "none"
  ) +
  facet_wrap(~ group, nrow = 2) +
  labs(
    x = "Wealth",
    y = "BMI"
  ) +
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(face = "bold"),
        plot.title  = element_text(size = 16))

## saving----
# ggsave("bmi_marginal_effects.jpeg", path = "/Users/levkolinski/Desktop/hookworm_bmi_bloodpressure/Figures/R1/",
#               units="in", width=10, height=7, dpi=350, device = "jpeg", bg = "white")
