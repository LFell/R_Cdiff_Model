# render_run.R
#
# Renders the model documents with a given run label:
#   render_run(label)     - Model1_ODE_Model_Run.qmd (main model run)
#                           -> <label>_C_diff_Model1_Run.html, outputs in Outputs/<label>/
#   render_oat(label)     - Model1_ODE_OAT_Run.qmd (one-at-a-time parameter testing)
#                           -> <label>_C_diff_Model1_OAT.html, outputs in Outputs/OAT/<label>/
#   render_run_T(label)   - Model1_ODE_Model_Run_T.qmd (main model run, transfer-state version)
#                           -> <label>_C_diff_Model1_Run_T.html,
#                              outputs in Outputs/Transfer_states/<label>/
#   render_compare(label) - Model1_ODE_Transfer_Compare.qmd (current vs transfer-state version)
#                           -> <label>_Transfer_Compare.html, outputs in Outputs/Transfer_compare/<label>/
# The HTML is saved in the project root. Rendering the same label again overwrites that file (and
# its outputs folder), so use a new label, e.g. a date, to keep an earlier report for comparison.
# Each document has its own outputs folder, so the same label can be used for all of them.
# The HTML is self-contained (embed-resources in the QMD YAML), so each copy keeps its own plots
# even after its outputs folder is overwritten by a later run.
#
# Usage:
#   source("render_run.R")
#   render_run("A")
#   render_run("B")
#   render_run("YYYYMMDD")
#   render_oat("YYYYMMDD")
#   render_run_T("YYYYMMDD")
#   render_compare("YYYYMMDD")

render_labelled <- function(input, label, output_suffix, outputs_path) {

  if (!nzchar(label)) stop("label must be a non-empty string")

  output_file <- paste0(label, output_suffix, ".html")

  quarto::quarto_render(
    input          = input,
    execute_params = list(run_label = label),
    output_file    = output_file
  )

  message("Run '", label, "' complete.")
  message("  Outputs : ", outputs_path, "/")
  message("  HTML    : ", output_file)
}

render_run <- function(label) {
  render_labelled("Model1_ODE_Model_Run.qmd", label, "_C_diff_Model1_Run",
                  file.path("Outputs", label))
}

render_oat <- function(label) {
  render_labelled("Model1_ODE_OAT_Run.qmd", label, "_C_diff_Model1_OAT",
                  file.path("Outputs", "OAT", label))
}

render_run_T <- function(label) {
  render_labelled("Model1_ODE_Model_Run_T.qmd", label, "_C_diff_Model1_Run_T",
                  file.path("Outputs", "Transfer_states", label))
}

render_compare <- function(label) {
  render_labelled("Model1_ODE_Transfer_Compare.qmd", label, "_Transfer_Compare",
                  file.path("Outputs", "Transfer_compare", label))
}
