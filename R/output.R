#' Print model and drug statistics
#'
#' Use this function to pretty print in an R notebook useful statistics for the
#' ensemble model analysis: how many drug combinations were tested by each model,
#' the number of models used and how many nodes each boolean network model had.
#'
#' @param drug.combs integer. Number of drug combinations tested
#' @param models integer. Number of models tested
#' @param nodes integer. Number of network nodes
#' @param html.output logical. If TRUE, the printed output will look nice in an
#' HTML document
#'
#' @importFrom usefun pretty_print_string print_empty_line
#' @export
print_model_and_drug_stats =
  function(drug.combs, models, nodes, html.output) {
    pretty_print_string(paste("Drug combinations tested:", drug.combs))
    print_empty_line(html.output)
    pretty_print_string(paste("Number of models:", models), with.gt = FALSE)
    print_empty_line(html.output)
    pretty_print_string(paste("Number of nodes:", nodes), with.gt = FALSE)
}

#' Print biomarkers for each predicted synergy
#'
#' @param biomarkers.dir string. It specifies the full path name (without the
#' ending character \emph{/}) of the directory which holds the biomarker files
#' for each drug combination in the \code{predicted.synergies}.
#' The biomarker files must be formatted as: \emph{\%drug.comb\%_biomarkers_active} or
#' \emph{\%drug.comb\%_biomarkers_inhibited}, where \%drug.comb\% is an element of
#' the \code{predicted.synergies} vector. If the files are not properly formatted
#' or don't even exist, zero biomarkers are reported.
#' @param predicted.synergies a character vector of the synergies (drug
#' combination names) that were predicted by \strong{at least one} of the models
#' in the dataset.
#' @param html.output logical. If TRUE, it makes the printed output nice for
#' an HTML document. Default value: TRUE.
#'
#' @importFrom usefun pretty_print_string print_empty_line pretty_print_vector_values pretty_print_bold_string
#' @importFrom utils read.table
#' @export
print_biomarkers_per_predicted_synergy =
  function(biomarkers.dir, predicted.synergies, html.output = TRUE) {
    pretty_print_string("")
    for (drug.comb in predicted.synergies) {
      # get the active biomarkers
      active.biomarkers.file =
        paste0(biomarkers.dir, "/", drug.comb, "_biomarkers_active")
      if (file.size(active.biomarkers.file) == 0
          || !file.exists(active.biomarkers.file)) {
        biomarkers.active.names = NULL
      } else {
        biomarkers.active =
          read.table(active.biomarkers.file, stringsAsFactors = FALSE)
        biomarkers.active.names = biomarkers.active[,1]
      }

      # get the inhibited biomarkers
      inhibited.biomarkers.file =
        paste0(biomarkers.dir, "/", drug.comb, "_biomarkers_inhibited")
      if (file.size(inhibited.biomarkers.file) == 0
          || !file.exists(inhibited.biomarkers.file)) {
        biomarkers.inhibited.names = NULL
      } else {
        biomarkers.inhibited =
          read.table(inhibited.biomarkers.file, stringsAsFactors = FALSE)
        biomarkers.inhibited.names = biomarkers.inhibited[,1]
      }

      # print biomarkers
      str = paste("Biomarkers for", drug.comb, "synergy prediction")
      pretty_print_bold_string(str, with.gt = FALSE, html.output = html.output)
      print_empty_line(html.output)
      print_empty_line(html.output)

      pretty_print_bold_string("Active biomarkers", with.gt = FALSE,
                               html.output = html.output)
      print_empty_line(html.output)
      pretty_print_vector_values(biomarkers.active.names, with.gt = FALSE)
      print_empty_line(html.output)
      print_empty_line(html.output)

      pretty_print_bold_string("Inhibited biomarkers", with.gt = FALSE,
                               html.output = html.output)
      print_empty_line(html.output)
      pretty_print_vector_values(biomarkers.inhibited.names, with.gt = FALSE)
      print_empty_line(html.output)
      print_empty_line(html.output)
    }
}
