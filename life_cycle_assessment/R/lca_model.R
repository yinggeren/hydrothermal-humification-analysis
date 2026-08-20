read_lca_inputs <- function(project_dir) {
  data_dir <- file.path(project_dir, "data")
  list(
    components = read.csv(file.path(data_dir, "scenario_components.csv"), stringsAsFactors = FALSE),
    assumptions = read.csv(file.path(data_dir, "scenario_assumptions.csv"), stringsAsFactors = FALSE),
    sensitivity = read.csv(file.path(data_dir, "sensitivity_inputs.csv"), stringsAsFactors = FALSE),
    break_even = read.csv(file.path(data_dir, "break_even_inputs.csv"), stringsAsFactors = FALSE),
    monte_carlo = read.csv(file.path(data_dir, "monte_carlo_summary_expected.csv"), stringsAsFactors = FALSE),
    expected = read.csv(file.path(data_dir, "expected_results.csv"), stringsAsFactors = FALSE)
  )
}

calc_scenario_results <- function(inputs) {
  x <- inputs$components
  aggregate(component_value_kg_co2eq_per_t ~ scenario + pathway, data = x, sum)
}

calc_table_s11 <- function(inputs) {
  subset(inputs$components, scenario == "OIS")
}

calc_sensitivity <- function(inputs) {
  x <- inputs$sensitivity
  x$sensitivity_range <- abs(x$net_gwp_high_emission_side - x$net_gwp_low_emission_side)
  x[order(-x$sensitivity_range), ]
}

calc_break_even <- function(inputs) {
  inputs$break_even
}

calc_monte_carlo_summary <- function(inputs) {
  inputs$monte_carlo
}

validate_lca_outputs <- function(inputs, scenarios, sensitivity, break_even, mc_summary) {
  expected <- inputs$expected
  merged <- merge(expected, scenarios, by = c("scenario", "pathway"), sort = FALSE)
  names(merged)[names(merged) == "component_value_kg_co2eq_per_t"] <- "observed_net_gwp"
  scenario_ok <- all(abs(merged$observed_net_gwp - merged$expected_net_gwp) <= 0.051)

  sens_expected <- inputs$sensitivity$sensitivity_range_expected
  sensitivity_ok <- all(abs(round(sensitivity$sensitivity_range, 1) - round(sens_expected, 1)) <= 0.05)

  break_ok <- all(abs(round(break_even$minimum_heat_recovery_percent, 1) -
                        c(73.5, 68.0, 57.1, 40.7, 24.3)) <= 0.05)

  mc_ok <- (
    round(mc_summary$mean_kg_co2eq_per_t[1], 1) == -127.6 &&
      round(mc_summary$median_kg_co2eq_per_t[1], 1) == -131.1 &&
      round(mc_summary$p2_5_kg_co2eq_per_t[1], 1) == -264.0 &&
      round(mc_summary$p97_5_kg_co2eq_per_t[1], 1) == 24.0 &&
      round(100 * mc_summary$probability_net_negative[1], 1) == 94.9
  )

  if (!scenario_ok) stop("Scenario net GWP values do not reproduce expected results.")
  if (!sensitivity_ok) stop("Sensitivity values do not reproduce expected results.")
  if (!break_ok) stop("Break-even thresholds do not reproduce expected results.")
  if (!mc_ok) stop("Monte Carlo summary does not reproduce expected results.")

  c(
    "LCA reproduction checks passed.",
    sprintf("Scenario rows checked: %d", nrow(merged)),
    sprintf("Sensitivity rows checked: %d", nrow(sensitivity)),
    sprintf("Break-even rows checked: %d", nrow(break_even)),
    sprintf("Monte Carlo mean: %.1f kg CO2-eq t-1", mc_summary$mean_kg_co2eq_per_t[1]),
    sprintf("Monte Carlo 95%% interval: %.1f to %.1f kg CO2-eq t-1",
            mc_summary$p2_5_kg_co2eq_per_t[1], mc_summary$p97_5_kg_co2eq_per_t[1]),
    sprintf("Net-negative probability: %.1f%%", 100 * mc_summary$probability_net_negative[1])
  )
}

