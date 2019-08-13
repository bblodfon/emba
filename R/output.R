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

#' yeah baby
#'
#' @importFrom usefun pretty_print_string print_empty_line pretty_print_vector_values pretty_print_bold_string
#' @importFrom utils read.table
print_biomarkers_per_predicted_synergy =
  function(biomarkers.dir, drug.comb, predicted.synergies, html.output = TRUE) {
    pretty_print_string("")
    for (drug.comb in predicted.synergies) {
      # get the active biomarkers
      active.biomarkers.file =
        paste0(biomarkers.dir, drug.comb, "_biomarkers_active")
      if (file.size(active.biomarkers.file) == 0) {
        biomarkers.active.names = NULL
      } else {
        biomarkers.active =
          read.table(active.biomarkers.file, stringsAsFactors = FALSE)
        biomarkers.active.names = biomarkers.active[,1]
      }

      # get the inhibited biomarkers
      inhibited.biomarkers.file =
        paste0(biomarkers.dir, drug.comb, "_biomarkers_inhibited")
      if (file.size(inhibited.biomarkers.file) == 0) {
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

#' We get the (previously-found) biomarkers from the respective file.
#' There are 3 policies regarding what to do with the 'new' biomarkers
#' when they share common nodes with the 'old' biomarkers. These are
#' given by the `method`:
#' 1) `replace` (default): if there is at least one new biomarker
#' common with the previously found ones, we discard the old ones and
#' write to the file the new biomarkers
#' 2) `prune.to.common`: keep only the common biomarkers
#' 3) `extend`: we add the non-common biomarkers to the old ones
#'
#' Note: No matter the `method`, if there are no common biomarkers
#' between old and new, we just add the new ones
#'
#' @importFrom usefun save_vector_to_file add_vector_to_df
update_biomarker_files =
  function(biomarkers.dir, drug.comb, biomarkers.active.new,
           biomarkers.inhibited.new, method = "replace") {
    # update the active biomarkers
    active.biomarkers.file =
      paste0(biomarkers.dir, drug.comb, "_biomarkers_active")

    if (file.size(active.biomarkers.file) == 0) {
      save_vector_to_file(vector = biomarkers.active.new,
                          file = active.biomarkers.file,
                          with.row.names = TRUE)
    } else {
      biomarkers.active.prev =
        read.table(active.biomarkers.file, stringsAsFactors = FALSE)
      biomarkers.active.prev.names = biomarkers.active.prev[,1]
      biomarkers.active.new.names = names(biomarkers.active.new)

      biomarkers.active.common = intersect(biomarkers.active.prev.names,
                                           biomarkers.active.new.names)

      if (length(biomarkers.active.common) == 0) {
        biomarkers.active = add_vector_to_df(biomarkers.active.prev,
                                             biomarkers.active.new)
        biomarkers.active = transform(biomarkers.active,
                                      V2 = as.numeric(biomarkers.active$V2))
        save_vector_to_file(vector = biomarkers.active,
                            file = active.biomarkers.file)
      } else {
        if (method == "replace") {
          biomarkers.active = biomarkers.active.new
          save_vector_to_file(vector = biomarkers.active,
                              file = active.biomarkers.file,
                              with.row.names = TRUE)
        } else if (method == "extend") {
          # find the non-common biomarkers and add them to the 'old' ones
          biomarkers.active.to.add = biomarkers.active.new[
            !(biomarkers.active.new.names %in% biomarkers.active.prev.names)
          ]
          biomarkers.active = add_vector_to_df(biomarkers.active.prev,
                                               biomarkers.active.to.add)
          biomarkers.active = transform(biomarkers.active,
                                        V2 = as.numeric(biomarkers.active$V2))
          save_vector_to_file(vector = biomarkers.active,
                              file = active.biomarkers.file)
        } else if (method == "prune.to.common") {
          biomarkers.active = biomarkers.active.new[biomarkers.active.common]
          save_vector_to_file(vector = biomarkers.active,
                              file = active.biomarkers.file,
                              with.row.names = TRUE)
        }
      }
    }

    # update the inhibited biomarkers
    inhibited.biomarkers.file =
      paste0(biomarkers.dir, drug.comb, "_biomarkers_inhibited")

    if (file.size(inhibited.biomarkers.file) == 0) {
      save_vector_to_file(vector = biomarkers.inhibited.new,
                          file = inhibited.biomarkers.file,
                          with.row.names = TRUE)
    } else {
      biomarkers.inhibited.prev =
        read.table(inhibited.biomarkers.file, stringsAsFactors = FALSE)
      biomarkers.inhibited.prev.names = biomarkers.inhibited.prev[,1]
      biomarkers.inhibited.new.names = names(biomarkers.inhibited.new)

      biomarkers.inhibited.common = intersect(biomarkers.inhibited.prev.names,
                                              biomarkers.inhibited.new.names)

      if (length(biomarkers.inhibited.common) == 0) {
        biomarkers.inhibited = add_vector_to_df(biomarkers.inhibited.prev,
                                                biomarkers.inhibited.new)
        biomarkers.inhibited = transform(biomarkers.inhibited,
                                         V2 = as.numeric(biomarkers.inhibited$V2))
        save_vector_to_file(vector = biomarkers.inhibited,
                            file = inhibited.biomarkers.file)
      } else {
        if (method == "replace") {
          biomarkers.inhibited = biomarkers.inhibited.new
          save_vector_to_file(vector = biomarkers.inhibited,
                              file = inhibited.biomarkers.file,
                              with.row.names = TRUE)
        } else if (method == "extend") {
          # find the non-common biomarkers and add them to the 'old' ones
          biomarkers.inhibited.to.add = biomarkers.inhibited.new[
            !(biomarkers.inhibited.new.names %in% biomarkers.inhibited.prev.names)
          ]
          biomarkers.inhibited = add_vector_to_df(biomarkers.inhibited.prev,
                                                  biomarkers.inhibited.to.add)
          biomarkers.inhibited = transform(biomarkers.inhibited,
                                           V2 = as.numeric(biomarkers.inhibited$V2))
          save_vector_to_file(vector = biomarkers.inhibited,
                              file = inhibited.biomarkers.file)
        } else if (method == "prune.to.common") {
          biomarkers.inhibited = biomarkers.inhibited.new[biomarkers.inhibited.common]
          save_vector_to_file(vector = biomarkers.inhibited,
                              file = inhibited.biomarkers.file,
                              with.row.names = TRUE)
        }
      }
    }
}
