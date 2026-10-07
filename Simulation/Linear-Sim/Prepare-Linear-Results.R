# Convert the five generator outputs into the inputs expected by the IV analyses.
# Run from this directory or use source(..., chdir = TRUE).
source("../../R/linear_simulation_functions.R", local = TRUE)

result_objects <- c(
  "Results_list_best_new_near",
  "Results_list_best_new_medium0.75",
  "Results_list_best_new_medium0.5",
  "Results_list_best_new_medium0.25",
  "Results_list_best_new_far"
)

for (scenario in seq_along(result_objects)) {
  input <- new.env(parent = emptyenv())
  load(paste0("model_results_", scenario, ".Rdata"), envir = input)
  raw_results <- input[[result_objects[scenario]]]
  if (length(raw_results) != 200L) {
    stop("Scenario ", scenario, " must contain 200 successful replicates.")
  }
  output <- new.env(parent = emptyenv())
  result_name <- paste0("Results", scenario)
  output[[result_name]] <- prepare_linear_results(raw_results)
  save(list = result_name, envir = output,
       file = paste0("linear_", scenario, ".Rdata"), compress = TRUE)
}
