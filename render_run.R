# render_run.R
#
# Renders Model1_ODE_Model_Run.qmd with a given run label.
# All outputs (CSVs and plots) are saved to Outputs/<label>/
# The rendered HTML is saved as YYMMDD_Model1_ODE_Model_Run_<label>.html, prefixed with today's date
# (the folder's usual date-prefix convention), so earlier reports are kept for comparison. Rendering
# the same label again on the same day overwrites that day's file. The HTML is self-contained
# (embed-resources in the QMD YAML), so each dated copy keeps its own plots even after Outputs/<label>/
# is overwritten by a later run.
#
# Usage:
#   source("render_run.R")
#   render_run("A")
#   render_run("B")

render_run <- function(label) {

  if (!nzchar(label)) stop("label must be a non-empty string")

  output_file <- paste0(format(Sys.Date(), "%y%m%d"), "_Model1_ODE_Model_Run_", label, ".html")

  quarto::quarto_render(
    input          = "Model1_ODE_Model_Run.qmd",
    execute_params = list(run_label = label),
    output_file    = output_file
  )

  message("Run '", label, "' complete.")
  message("  Outputs : Outputs/", label, "/")
  message("  HTML    : ", output_file)
}
