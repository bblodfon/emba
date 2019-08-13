#' Load the models predictions data
#'
#' Use this function to read a file that has the model predictions data
#' and output it to a \code{data.frame} object.
#'
#' @param model.predictions.file a tab-delimited file (for the specific format
#' check the example below)
#'
#' @return a \code{data.frame} object with rows the models and columns the
#' drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#'
#' @examples
#'
#' model.predictions.file = system.file("extdata", "model_predictions",
#'   package = "emba", mustWork = TRUE)
#' model.predictions = get_model_predictions(model.predictions.file)
#'
#' @importFrom utils read.table
#' @export
get_model_predictions = function(model.predictions.file) {
  #print(paste("Reading model predictions file:", model.predictions.file))

  lines = readLines(model.predictions.file)
  lines[1] = sub("ModelName\t|#ModelName\t", "", lines[1])
  tmp.file = "model_predictions.tab"
  writeLines(lines, tmp.file)
  model.data = read.table("model_predictions.tab",  check.names = F)

  if (file.exists(tmp.file)) invisible(file.remove(tmp.file))
  for (i in 1:length(colnames(model.data))) {
    colnames(model.data)[i] = gsub("\\[|\\]", "", colnames(model.data)[i])
  }

  return(model.data)
}

#' Load the observed synergies data
#'
#' Use this function to read a file that has the observed synergies data and
#' output it to a character vector. If \code{drug.combinations.tested}
#' is NULL (the default), no data validation is done, otherwise we check that
#' the observed synergies are indeed a subset of the tested drug combinations.
#'
#' @param file string. The name of the file, can be a full path. See example
#' below for the format of an observed synergies file.
#' @param drug.combinations.tested a character vector with drug combinations
#' as elements. Default value: NULL.
#'
#' @return a character vector with elements the names of the drug combinations
#' that were found as synergistic
#'
#' @examples
#' observed.synergies.file = system.file("extdata", "observed_synergies",
#'   package = "emba", mustWork = TRUE)
#' observed.synergies = get_observed_synergies(observed.synergies.file)
#'
#' @export
get_observed_synergies =
  function(file, drug.combinations.tested = NULL) {
    #print(paste("Reading observed synergies file:", file))

    lines = readLines(file)
    observed.synergies = gsub("~", "-", lines)

    if (!is.null(drug.combinations.tested)) {
      validate_observed_synergies_data(observed.synergies, drug.combinations.tested)
    }

    return(observed.synergies)
}

#' Load the models stable state data
#'
#' Use this function to merge the stable states from all models into a single
#' matrix. The models stable states are loaded from \emph{.gitsbe} files that can
#' be found inside the given \code{models.dir} directory.
#'
#' @param models.dir string. A dir with \emph{.gitsbe} files/models
#'
#' @return a matrix (nxm) with n models and m nodes. The row names of the matrix
#' specify the models' names whereas the column names specify the name of the
#' network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
#' element} are either \emph{0} (inactive node) or \emph{1} (active node).
#'
#' @examples
#'
#' models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
#' models.stable.state = get_stable_state_from_models_dir(models.dir)
#'
#' @export
get_stable_state_from_models_dir = function(models.dir) {
  files = list.files(models.dir)
  model.stable.states = character(length(files))

  node.names = get_node_names(models.dir)

  i = 0
  for (file in files) {
    i = i + 1
    lines = readLines(paste0(models.dir, "/", file))
    model.stable.states[i] = gsub("stablestate: ", "", lines[4])
  }

  models.stable.state = data.frame(model.stable.states, row.names = files)
  df = apply(models.stable.state, 1, function(x) {
    as.numeric(strsplit(as.character(x[1]), "")[[1]])
  })
  rownames(df) = node.names

  return(t(df))
}

#' Load the models boolean equation link operator data
#'
#' Use this function to merge the link operator data used in the boolean equations
#' of the models into a single matrix. Every boolean model is defined by a series
#' of boolean equations in the form \eqn{Target *= (Activator OR Activator OR...)
#' AND NOT (Inhibitor OR Inhibitor OR...)"}. The \strong{link operator} can be
#' either \emph{AND NOT}, \emph{OR NOT} or non-existent if the target has only
#' activating regulators or only inhibiting ones. The models are loaded from
#' \emph{.gitsbe} files that can be found inside the given \code{models.dir}
#' directory.
#'
#' @param models.dir string. A dir with \emph{.gitsbe} files/models
#' @param remove.equations.without.link.operator logical. Should we keep the
#' nodes (columns in the returned matrix) which do not have both type of
#' regulators (so no link operator)? Default value: TRUE (remove these nodes).
#'
#' @return a matrix (nxm) with n models and m nodes. The row names of the matrix
#' specify the models' names whereas the column names specify the name of the
#' network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
#' element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
#' (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
#' both activating and inhibiting regulators (no link operator).
#'
#' @examples
#'
#' models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
#' models.link.operator = get_link_operators_from_models_dir(models.dir)
#' models.link.operator.with.extra.nodes =
#'   get_link_operators_from_models_dir(models.dir, FALSE)
#'
#' @export
get_link_operators_from_models_dir =
  function(models.dir, remove.equations.without.link.operator = TRUE) {
    files = list.files(models.dir)
    node.names = get_node_names(models.dir)

    datalist = list(length(files))

    # get the equations
    i = 0
    for (file in files) {
      i = i+1
      lines = readLines(paste0(models.dir, "/", file))
      equations = grep("equation:", lines, value = TRUE)
      values = sapply(equations, function(equation) {
        assign_link_operator_value_to_equation(equation)})
      datalist[[i]] = values
    }

    df = do.call(rbind, datalist)

    rownames(df) = files
    colnames(df) = node.names

    if (remove.equations.without.link.operator) {
      # keep only the equations (columns) that have the 'and not' or 'or not'
      # link operator, i.e. those that can change in the 'link mutations'
      df = df[, colSums(is.na(df)) < nrow(df)]
    } else {
      # keep all equations and put a value of 0.5 for those that don't have a
      # link operator
      df[is.na(df)] = 0.5
    }

    return(df)
}

get_fitness_from_models_dir = function(models.dir) {
  files = list.files(models.dir)
  model.fitness = character(length(files))

  i = 0
  for (file in files) {
    i = i + 1
    lines = readLines(paste0(models.dir, "/", file))
    model.fitness[i] = gsub("fitness: ", "", lines[3])
  }

  model.fitness = as.numeric(model.fitness)
  names(model.fitness) = files

  return(model.fitness)
}

#' Get the node names
#'
#' This function uses the first .gitsbe file that it finds inside the given
#' directory to output a vector of the network node names (which should be the
#' same for every model)
#'
#' @param models.dir string. A dir with \emph{.gitsbe} files/models
#'
#' @return a character vector of the node names (protein and/or gene names)
#'
#' @examples
#'
#' models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
#' nodes = get_node_names(models.dir)
#'
#' @export
get_node_names = function(models.dir) {
  # use the first .gitsbe model file to derive the node names
  file.lines = readLines(paste0(models.dir, "/", list.files(models.dir)[1]))
  node.names = gsub("mapping: (.*) =.*", "\\1",
                    grep("mapping:", file.lines, value = TRUE))
  return(node.names)
}

#' Get the model names
#'
#' @param models.dir string. A dir with \emph{.gitsbe} files/models
#'
#' @return a character vector of the model names, corresponding to the names
#' of the \emph{.gitsbe} files.
#'
#' @examples
#'
#' models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
#' models = get_model_names(models.dir)
#'
#' @export
get_model_names = function(models.dir) {
  return(list.files(models.dir))
}

#' Assign link operator value to boolean equation
#'
#' @param equation string. The boolean equation in the form
#' \eqn{Target *= (Activator OR Activator OR...)
#' AND NOT (Inhibitor OR Inhibitor OR...)"}
#'
#' @return \strong{1} if the \code{equation} has the '\emph{or not}' link operator,
#' \strong{0} if the \code{equation} has the '\emph{and not}' link operator and
#' \strong{NA} if it has neither.
assign_link_operator_value_to_equation = function(equation) {
  if (grepl(".*or not.*", equation)) {
    return(1)
  } else if (grepl(".*and not.*", equation)) {
    return(0)
  } else return(NA)
}

#' A title
#'
#' gettin somethign
#'
#'
#'
#' @importFrom usefun remove_commented_and_empty_lines
get_consensus_steady_state = function(steady.state.file) {
  #print(paste("Reading consensus steady state file:", steady.state.file))

  lines = readLines(steady.state.file)
  lines = remove_commented_and_empty_lines(lines)
  consensus.steady.state = build_consensus_steady_state_vector(lines)

  return(consensus.steady.state)
}

build_consensus_steady_state_vector = function(lines) {
  node.names = character(0)
  activity.states = character(0)
  for (line in lines) {
    values = strsplit(line, "\t")[[1]]
    node.names = c(node.names, values[1])
    activity.states = c(activity.states, values[2])
  }
  activity.states = as.numeric(activity.states)
  stopifnot(length(activity.states) == length(node.names))

  names(activity.states) = node.names
  return(activity.states)
}

#' Is drug combination element of given vector?
#'
#' Use this function to determine if a drug combination is part of a vector of
#' other drug combinations. We take care only of pair-wise drug combinations and
#' an internal check is done for alternative drug names, e.g. we check if
#' \emph{A-B} combination is included, but also for \emph{B-A}.
#'
#' @param drug.comb a string in the form \emph{A-B} (no spaces between the names
#' and the hyphen '-')
#' @param comb.vector a character vector of drug combinations, each one in the
#' form \emph{drugname.1-drugname.2}
#'
#' @return logical, depending if the drug combination is element of the given
#' vector or not.
#'
#' @examples
#' # TRUE
#' is_comb_element_of("A-B", c("E-F", "A-B"))
#' is_comb_element_of("B-A", c("E-F", "A-B"))
#'
#' # FALSE
#' is_comb_element_of("A-B", c("E-F", "A-D"))
#' is_comb_element_of("A-B", c())
#'
#' @export
is_comb_element_of = function(drug.comb, comb.vector) {
  return(is.element(drug.comb, comb.vector) |
           is.element(get_alt_drugname(drug.comb), comb.vector))
}

#' Validate observed synergies data
#'
#' This function checks that the observed synergies are part (a subset) of the
#' tested drug combinations
#'
#' @param observed.synergies a character vector of drug combinations
#' @param drug.combinations.tested a character vector of drug combinations
#'
#' @return NULL if no errors found, otherwise stops execution.
validate_observed_synergies_data =
  function(observed.synergies, drug.combinations.tested) {
    for (drug.comb in observed.synergies) {
      if (!is.element(drug.comb, drug.combinations.tested) &&
          !is.element(get_alt_drugname(drug.comb), drug.combinations.tested)) {
        stop(paste("Drug Combination: ", drug.comb,
                   "is not listed in the observed synergies file"), call. = F)
      }
    }
}

#' Get alternative drug name
#'
#' Use this function on a string \emph{A-B} that represents a drug combination,
#' to get the reverse combination name - \emph{B-A} - for testing/checking data.
#'
#' @param drug.comb a string in the form \emph{drugname.1-drugname.2} (no
#' spaces between the names and the hyphen '-')
#'
#' @return the alternative, yet equivalent drug combination
#'
#' @examples
#' drug.comb = "A-B"
#' alt.drug.comb = get_alt_drugname(drug.comb)
#'
#' @export
get_alt_drugname = function(drug.comb) {
  drug.list = unlist(strsplit(drug.comb,"-"))
  drug.comb.alt = paste0(drug.list[2], "-", drug.list[1])
  return(drug.comb.alt)
}

#' Construct igraph network graph
#'
#' Use this function to create an igraph graph object based on the topology .sif
#' file given. It automatically sets various visualization graph properties and
#' checks if the node names from the topology file are the same as in the models
#' inside the given \code{models.dir} (if not NULL).
#'
#' @param topology.file string. The name of the .sif file (can be a full path
#' name).
#' @param models.dir string. A dir with \emph{.gitsbe} files/models. Default
#' value: NULL. If specified, it is used for the validation of the node names.
#'
#' @return an igraph graph object representing the network as defined in the
#' topology file
#'
#' @seealso \code{\link[igraph]{graph_from_data_frame}},
#' \code{\link[emba]{get_edges_from_topology_file}},
#' \code{\link[emba]{get_node_names}}
#'
#' @importFrom igraph graph_from_data_frame V V<- E E<-
#' @importFrom utils read.table
#' @export
construct_network = function(topology.file, models.dir = NULL) {
  edges = get_edges_from_topology_file(topology.file)

  net = graph_from_data_frame(edges, directed = TRUE)

  # check the vertices/node names if models.dir is not NULL
  if (!is.null(models.dir)) {
    vertices = V(net)$name
    nodes = get_node_names(models.dir)
    stopifnot(all(sort(nodes) == sort(vertices)))
  }

  # set visualization graph properties
  E(net)$width = 1.5
  E(net)$arrow.size = 0.4
  E(net)$curved = 0.4
  V(net)$label.cex = 0.6
  V(net)$size = 10

  return(net)
}

#' Get the edges from a specified topology
#'
#' Use this function to read a topology .sif file (either space or tab-delimited)
#' and get a matrix of network edges specifying the source and target name, the
#' regulation effect (activation or inhibition) and the color (green or red) of
#' each interaction.
#'
#' @param topology.file string. The name of the .sif file (can be a full path
#' name).
#'
#' @return a matrix with as many rows as in the .sif topology file (each row is
#' an edge) and 4 columns defining the source and target node name, the
#' regulation (activation or inhibition) and the color (green or red) of the
#' signed interaction.
#'
#' @export
get_edges_from_topology_file = function(topology.file) {
  #print(paste("Reading topology file:", topology.file))

  edges = read.table(topology.file)

  # reorder&rename columns
  edges = as.matrix(edges[,c(1,3,2)])
  colnames(edges) = c("source", "target", "regulation.effect")

  # change arrow symbols for activation and inhibition to proper name strings
  regulation.effects = edges[,"regulation.effect"]
  regulation.effects = sapply(regulation.effects, function(arrow.symbol) {
    if (arrow.symbol == "->") return("activation")
    else return("inhibition")
  }, USE.NAMES = FALSE)

  edges[,"regulation.effect"] = regulation.effects

  # Set edge.color plotting parameter (igraph) according to regulation effect
  color = sapply(regulation.effects, function(effect) {
    if (effect == "activation") return("green")
    else return("red")
  })
  names(color) = NULL

  edges = cbind(edges, color)

  return(edges)
}

#' returns a data.frame, where the columns represent the network nodes, the rows
#' represent the predicted synergies and the values can either 1 (active biomarker),
#' -1 (inhibited biomarker) or 0 (not a biomarker)
#'
#' @importFrom utils read.table
get_biomarkers_per_synergy =
  function(predicted.synergies, biomarkers.dir, models.dir) {
    # initialize res data.frame
    node.names = get_node_names(models.dir)
    res = as.data.frame(matrix(0, ncol = length(node.names),
                                  nrow = length(predicted.synergies)))
    colnames(res) = node.names
    rownames(res) = predicted.synergies

    for (drug.comb in predicted.synergies) {
      # insert the active biomarkers
      active.biomarkers.file =
        paste0(biomarkers.dir, drug.comb, "_biomarkers_active")

      if (file.size(active.biomarkers.file) != 0) {
        biomarkers.active =
          read.table(active.biomarkers.file, stringsAsFactors = FALSE)
        biomarkers.active.names = biomarkers.active[,1]
        # biomarkers.active.values = biomarkers.active[,2]

        res[drug.comb, biomarkers.active.names] = 1
      }

      # insert the inhibited biomarkers
      inhibited.biomarkers.file =
        paste0(biomarkers.dir, drug.comb, "_biomarkers_inhibited")
      if (file.size(inhibited.biomarkers.file) != 0) {
        biomarkers.inhibited =
          read.table(inhibited.biomarkers.file, stringsAsFactors = FALSE)
        biomarkers.inhibited.names = biomarkers.inhibited[,1]
        # biomarkers.inhibited.values = biomarkers.inhibited[,2]

        res[drug.comb, biomarkers.inhibited.names] = -1
      }
    }

    return(res)
}

#' `biomarkers.dirs` is a vector of the cell lines' biomarker directories
#' and `type` can be either 'active' or 'inhibited'
#'
#' @importFrom utils read.table
get_perf_biomarkers_per_cell_line = function(biomarkers.dirs, type) {
  if (type == "active")
    biomarker.type.extension = "/biomarkers_active"
  else
    biomarker.type.extension = "/biomarkers_inhibited"

  biomarkers.perf = sapply(biomarkers.dirs, function(biomarkers.dir) {
    biomarkers.file = paste0(biomarkers.dir, biomarker.type.extension)
    if (file.size(biomarkers.file) == 0)
      return(list(NULL)) # empty list
    else
      return(
        read.table(biomarkers.file, stringsAsFactors = FALSE)
      )
  })
  # add the cell line name as id
  names(biomarkers.perf) = names(biomarkers.dirs)

  return(biomarkers.perf)
}

# `biomarkers.dirs` is a vector of the cell lines' biomarker directories
# Returns a list of cell-line data frames with rows the true positive
# predicted synergies for each cell line and columns the network nodes (same
# for all). So, cell line (list) => biomarkers of a synergy (row of data.frame)
get_synergy_biomarkers_per_cell_line = function(biomarkers.dirs) {
  biomarkers.per.synergy = list()
  for (i in seq_along(biomarkers.dirs)) {
    biomarkers.file = paste0(biomarkers.dirs[i], "/biomarkers_per_synergy")
    biomarkers.per.synergy[[i]] =
      read.table(file = biomarkers.file, stringsAsFactors = FALSE,
                 check.names = FALSE)
  }
  names(biomarkers.per.synergy) = names(biomarkers.dirs)

  return(biomarkers.per.synergy)
}
