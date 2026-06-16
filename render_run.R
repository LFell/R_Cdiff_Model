# render_run.R
#
# Renders Model1_ODE_Model_Run.qmd with a given run label.
# All outputs (CSVs and plots) are saved to Outputs/<label>/
# The rendered HTML is saved as Model1_ODE_Model_Run_<label>.html
#
# Usage:
#   source("render_run.R")
#   render_run("A")
#   render_run("B")

render_run <- function(label) {

  if (!nzchar(label)) stop("label must be a non-empty string")

  quarto::quarto_render(
    input          = "Model1_ODE_Model_Run.qmd",
    execute_params = list(run_label = label),
    output_file    = paste0("Model1_ODE_Model_Run_", label, ".html")
  )

  message("Run '", label, "' complete.")
  message("  Outputs : Outputs/", label, "/")
  message("  HTML    : Model1_ODE_Model_Run_", label, ".html")
}
