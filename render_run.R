# render_run.R
#
# Renders the current model (Model 4) documents with a given run label:
#   render_run(label) - Model4_Run.qmd (main model run)
#                       -> <label>_Model4_Run.html, outputs in Outputs/Model4/<label>/
#   render_oat(label) - Model4_OAT_Test.qmd (one-at-a-time parameter test)
#                       -> <label>_Model4_OAT_Test.html, outputs in Outputs/Model4/OAT/<label>/
# The HTML is saved in the project root. Rendering the same label again overwrites that file (and
# its outputs folder), so use a new label, e.g. a date, to keep an earlier report for comparison.
# The HTML is self-contained (embed-resources in the QMD YAML), so each copy keeps its own plots
# even after its outputs folder is overwritten by a later run.
#
# Earlier model versions (and their OAT test and comparison reports) are not rendered from here:
# check out their git tags (model-1 to model-3) instead; see the model versions table in README.md.
#
# Usage:
#   source("render_run.R")
#   render_run("A")
#   render_run("YYYYMMDD")
#   render_oat("YYYYMMDD")

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
  render_labelled("Model4_Run.qmd", label, "_Model4_Run",
                  file.path("Outputs", "Model4", label))
}

render_oat <- function(label) {
  render_labelled("Model4_OAT_Test.qmd", label, "_Model4_OAT_Test",
                  file.path("Outputs", "Model4", "OAT", label))
}
