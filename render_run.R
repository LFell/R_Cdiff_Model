# render_run.R
#
# Renders Model1_ODE_Model_Run.qmd with a given run label.
# All outputs (CSVs and plots) are saved to Outputs/<label>/
# The rendered HTML is saved as <label>_C_diff_Model1_Run.html in the project root. Rendering the
# same label again overwrites that file (and Outputs/<label>/), so use a new label, e.g. a date, to
# keep an earlier report for comparison. The HTML is self-contained (embed-resources in the QMD
# YAML), so each copy keeps its own plots even after Outputs/<label>/ is overwritten by a later run.
#
# Usage:
#   source("render_run.R")
#   render_run("A")
#   render_run("B")
#   render_run("YYYYMMDD")

render_run <- function(label) {

  if (!nzchar(label)) stop("label must be a non-empty string")

  output_file <- paste0(label, "_C_diff_Model1_Run", ".html")

  quarto::quarto_render(
    input          = "Model1_ODE_Model_Run.qmd",
    execute_params = list(run_label = label),
    output_file    = output_file
  )

  message("Run '", label, "' complete.")
  message("  Outputs : Outputs/", label, "/")
  message("  HTML    : ", output_file)
}
