# Run all LCA calculations and write reproducible outputs.

args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
if (length(script_flag) > 0) {
  project_dir <- normalizePath(dirname(sub("^--file=", "", script_flag[1])), winslash = "/")
} else {
  project_dir <- normalizePath(getwd(), winslash = "/")
}

source(file.path(project_dir, "R", "lca_model.R"))

output_dir <- file.path(project_dir, "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

inputs <- read_lca_inputs(project_dir)

scenario_results <- calc_scenario_results(inputs)
write.csv(scenario_results, file.path(output_dir, "scenario_results.csv"), row.names = FALSE)

table_s11 <- calc_table_s11(inputs)
write.csv(table_s11, file.path(output_dir, "table_s11_contributions.csv"), row.names = FALSE)

sensitivity <- calc_sensitivity(inputs)
write.csv(sensitivity, file.path(output_dir, "sensitivity_results.csv"), row.names = FALSE)

break_even <- calc_break_even(inputs)
write.csv(break_even, file.path(output_dir, "break_even_results.csv"), row.names = FALSE)

mc_summary <- calc_monte_carlo_summary(inputs)
write.csv(mc_summary, file.path(output_dir, "monte_carlo_summary.csv"), row.names = FALSE)

validation <- validate_lca_outputs(inputs, scenario_results, sensitivity, break_even, mc_summary)
writeLines(validation, file.path(output_dir, "validation_report.txt"))
cat(paste(validation, collapse = "\n"), "\n")

