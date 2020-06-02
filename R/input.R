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
#' @importFrom tidyr %>%
#' @importFrom tibble column_to_rownames
#' @export
get_model_predictions = function(model.predictions.file) {
  # Get the column names (drug combination names)
  first.line = readLines(model.predictions.file, n = 1)
  first.line = sub("ModelName\t|#ModelName\t", "", first.line)
  drug.combo.names = unlist(strsplit(first.line, split = '\t'))
  drug.combo.names = sapply(drug.combo.names, function(x) {
    gsub(pattern = "\\[|\\]", replacement = "", x)
    }, USE.NAMES = F)

  model.data = read.table(model.predictions.file, header = FALSE, skip = 1, stringsAsFactors = FALSE)
  model.data = model.data %>% tibble::column_to_rownames(var = "V1")
  colnames(model.data) = drug.combo.names

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
#' \code{data.frame} object. The models stable states are loaded from \emph{.gitsbe}
#' files that can be found inside the given \code{models.dir} directory.
#'
#' @param models.dir string. A dir with \emph{.gitsbe} files/models.
#' \strong{Do not} include the ending path character in the string (\emph{/}).
#'
#' @return a \code{data.frame} (nxm) with n models and m nodes. The row names
#' specify the models' names whereas the column names specify the name of the
#' network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
#' element} are either \emph{0} (inactive node) or \emph{1} (active node). If a
#' \emph{.gitsbe} file/model has zero (0) or more than 1 stable states, a diagnostic
#' message is printed and the corresponding model is discarded, i.e. it will not
#' be included in the returned \code{data.frame} object.
#'
#' @examples
#'
#' models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
#' models.stable.state = get_stable_state_from_models_dir(models.dir)
#'
#' @importFrom stringr str_count
#' @importFrom tidyr %>% separate
#' @importFrom dplyr mutate_if
#'
#' @export
get_stable_state_from_models_dir = function(models.dir) {
  files = list.files(models.dir)
  ss_vec = character(length(files))
  files_to_keep = character(length(files))

  node_names = get_node_names(models.dir)

  index = 1
  for (file in files) {
    lines = readLines(paste0(models.dir, "/", file))
    res = stringr::str_count(string = lines, pattern = "stablestate:")
    if (sum(res) == 1) {
      ss_index = which(res == 1)
      ss_vec[index] = gsub(pattern = "stablestate: ", replacement = "",
        x = lines[ss_index])
      files_to_keep[index] = files[index]
    } else {
      message('Number of stable states different than 1 for model: ', file)
    }
    index = index + 1
  }

  ss_vec = ss_vec[ss_vec != ""]
  files = files_to_keep[files_to_keep != ""]

  df = data.frame(ss_vec)

  df = df %>%
    separate(col = ss_vec, into = node_names, sep = 1:length(node_names)) %>%
    mutate_if(is.character, as.numeric)

  # remove .gitsbe extension from files
  model_names = sapply(files, function(x) {
    sub(pattern = ".gitsbe", replacement = "", x)
  }, USE.NAMES = FALSE)

  rownames(df) = model_names

  return(df)
}

#' Load the models boolean equation link operator data
#'
#' Use this function to merge the link operator data used in the boolean equations
#' of the models into a single \code{data.frame} object. Every boolean model is defined by a series
#' of boolean equations in the form \eqn{Target *= (Activator or Activator or...)
#' and not (Inhibitor or Inhibitor or...)"}. The \strong{link operator} can be
#' either \emph{and not}, \emph{or not} or non-existent if the target has only
#' activating regulators or only inhibiting ones (the \emph{not} remains in the
#' latter case). The models are loaded from
#' \emph{.gitsbe} files that can be found inside the given \code{models.dir}
#' directory.
#'
#' @param models.dir string. A directory path with \emph{.gitsbe} files/models.
#' \strong{Do not} include the ending path character in the string (\emph{/}).
#' @param remove.equations.without.link.operator logical. Should we keep the
#' nodes (columns in the returned matrix) which do not have both type of
#' regulators (so no link operator)? Default value: TRUE (remove these nodes).
#'
#' @return a \code{data.frame} (nxm) with n models and m nodes. The row names
#' specify the models' names whereas the column names specify the
#' network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
#' element} are either \emph{0} (\strong{and not} link operator), \emph{1}
#' (\strong{or not} link operator) or \emph{0.5} if the node is not targeted by
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

    model.names = sapply(files, function(x) {
      sub(pattern = ".gitsbe", replacement = "", x)
    }, USE.NAMES = FALSE)
    node.names = get_node_names(models.dir)

    datalist = list()

    # get the equations
    index = 1
    for (file in files) {
      lines = readLines(paste0(models.dir, "/", file))
      equations = grep("equation:", lines, value = TRUE)
      values = sapply(equations, function(equation) {
        assign_link_operator_value_to_equation(equation)}, USE.NAMES = FALSE)
      datalist[[index]] = values
      index = index + 1
    }

    df = do.call(rbind, datalist)

    rownames(df) = model.names
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

    return(as.data.frame(df))
}

#' Load the models fitness scores
#'
#' Use this function to merge the fitness scores from all models into a single
#' vector (the fitness score is a value between 0 and 1 and denotes how close
#' was the model fitted to one or more training data observations). Each model's
#' fitness value is loaded from the respective \emph{.gitsbe} file that can be
#' found inside the given \code{models.dir} directory.
#'
#' @param models.dir string. A dir with \emph{.gitsbe} files/models
#'
#' @return a numeric vector with elements the fitness scores and the names of the
#' models included in the \emph{names} attribute.
#'
#' @examples
#'
#' models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
#' models.fitness = get_fitness_from_models_dir(models.dir)
#'
#' @export
get_fitness_from_models_dir = function(models.dir) {
  files = list.files(models.dir)

  model.names = sapply(files, function(x) {
    sub(pattern = ".gitsbe", replacement = "", x)
  }, USE.NAMES = FALSE)

  model.fitness = character(length(files))

  i = 0
  for (file in files) {
    i = i + 1
    lines = readLines(paste0(models.dir, "/", file))
    model.fitness[i] = gsub("fitness: ", "", lines[3])
  }

  model.fitness = as.numeric(model.fitness)
  names(model.fitness) = model.names

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
#' of the \emph{.gitsbe} files (the extension is pruned).
#'
#' @examples
#'
#' models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
#' models = get_model_names(models.dir)
#'
#' @export
get_model_names = function(models.dir) {
  model_names = sapply(list.files(models.dir), function(x) {
    sub(pattern = ".gitsbe", replacement = "", x)
  }, USE.NAMES = FALSE)
  return(model_names)
}

#' Assign link operator value to boolean equation
#'
#' @param equation string. The boolean equation in the form
#' \eqn{Target *= (Activator or Activator or...)
#' and not (Inhibitor or Inhibitor or...)"}
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
#' @param observed.synergies a non-empty character vector of drug combinations
#' @param drug.combinations.tested a non-empty character vector of drug combinations
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

#' Get alternative drug combination name
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
#' \code{\link{get_edges_from_topology_file}},
#' \code{\link{get_node_names}}
#'
#' @examples
#' topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
#' net = construct_network(topology.file)
#'
#' @importFrom igraph graph_from_data_frame V V<- E E<-
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
#' @examples
#' topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
#' edges = get_edges_from_topology_file(topology.file)
#'
#' @importFrom utils read.table
#' @export
get_edges_from_topology_file = function(topology.file) {
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

#' Filter the network's vertices
#'
#' Produce an \strong{induced subgraph} of the given \emph{net} igraph object.
#' How many vertices/nodes will be kept in the result graph object is determined
#' by the initial nodes given and the level provided. A level equal to 0 corresponds
#' to a subgraph with only the given nodes, a level equal to 1 to a subgraph with
#' the nodes + their neighbors (the closed neighbourhood set where every node is within
#' 1 edge distnace from the given ones) and a level equal to 2 to a subgraph
#' with the nodes + their neighbors + the nodes neighbor neighbors!
#' (so the neighbourhood of the neighbourhood or every node is within 2 edges
#' distance from the given ones).
#'
#' @param net an igraph object.
#' @param nodes character vector of node names. It must be a subset of the nodes
#' of the \emph{net} object.
#' @param level integer. Can be only 0, 1 or 2 and specifies the neighbourhood
#' depth of the result graph.
#'
#' @return an induced subgraph of the \code{net} igraph object.
#'
#' @importFrom igraph V induced_subgraph
#' @export
filter_network = function(net, nodes, level) {
  stopifnot(class(net) == "igraph")
  stopifnot(level == 0 | level == 1 | level == 2)
  stopifnot(!is_empty(nodes))

  nodes = unique(nodes)
  net.nodes = V(net)$name
  stopifnot(nodes %in% net.nodes)

  if (level == 0) {
    neighbor.nodes = c()
  }
  if (level == 1) {
    neighbor.nodes = get_neighbors(net, nodes)
  } else if (level == 2) {
    neighbor.nodes.lev.1 = get_neighbors(net, nodes)
    neighbor.nodes = c(neighbor.nodes.lev.1, get_neighbors(net, neighbor.nodes.lev.1))
  }

  nodes = unique(c(nodes, neighbor.nodes))
  induced_subgraph(net, nodes)
}

#' Get neighbor nodes
#'
#' Given an igraph network object and vector of node names, this function
#' returns the set of unique neighbor nodes considering both ingoing and outgoing
#' edges (the closed neighbourhood node set).
#'
#' @param net igraph object
#' @param nodes character vector of node names
#'
#' @return a character vector of all the unique neighbors of the given
#' \code{nodes} in the \code{net} graph.
#'
#' @importFrom igraph adjacent_vertices
#' @importFrom purrr flatten_chr
get_neighbors = function(net, nodes) {
  adj.node.list = adjacent_vertices(graph = net, v = nodes, mode = "all")
  return(unique(names(flatten_chr(adj.node.list))))
}

#' Get observed synergies per cell line
#'
#' Use this function to get the observed synergies from the respective
#' files inside the given list of cell line directories.
#'
#' @param cell.line.dirs a character vector of the cell line directories, in the
#' form of \emph{\{path\}/cell_line_name}. The cell line name directory
#' should be different for each element of the vector as we use it to fill in the
#' \code{rownames} of the result \code{data.frame} object. Inside each cell line directory
#' we read the observed synergies from a file called \emph{observed_synergies}
#' (if it exists and is non-empty). This file has the names of the observed
#' drug combinations, one in each line.
#' @param drug.combos a character vector with elements the names of all the drug combinations
#' that were tested in the analysis.
#'
#' @return a data.frame, whose columns represent the drug combinations tested
#' and the rows the cell lines. Possible values for each \emph{cell line-drug combination}
#' element are either \emph{1} (an observed synergy) or \emph{0} (non-observed synergy).
#'
#' @importFrom usefun add_row_to_ternary_df
#' @export
get_observed_synergies_per_cell_line = function(cell.line.dirs, drug.combos) {

  # initialize res data.frame
  res = as.data.frame(matrix(data = NA, nrow = 0, ncol = length(drug.combos)))
  colnames(res) = drug.combos

  for (cell.line.dir in cell.line.dirs) {
    # observed synergies
    observed.synergies.file = paste0(cell.line.dir, "/observed_synergies")

    if (file.exists(observed.synergies.file)
        && file.size(observed.synergies.file) != 0) {
      synergies = get_observed_synergies(observed.synergies.file)
    } else {
      synergies = c()
    }

    cell.line = basename(cell.line.dir)

    res = add_row_to_ternary_df(df = res, values.pos = synergies,
                                values.neg = c(),
                                row.name = cell.line)
  }

  return(res)
}

#' Get synergy scores from file
#'
#' Use this function to read \href{https://druglogics.github.io/druglogics-doc/drabme-install.html#drabme-output}{Drabme's
#' output files} that have synergy scores for a list of tested drug perturbations.
#'
#' Two types of files can be read: the \emph{model-wise} and the \emph{ensemble-wise}
#' synergies:
#' \itemize{
#'   \item The \emph{model-wise} synergies data is structured as a 3-column
#' table, first being the names of the tested drug combinations, second the
#' number of models that predicted that combination as synergistic and the third
#' the number of models that predicted the drug combination as non-synergistic.
#'   \item The \emph{ensemble-wise} synergies data is structured as a 2-column table,
#' the first being the names of drug combinations and the second the ensemble-wise
#' synergy scores for those perturbations.
#'  }
#'
#' Note that no matter the file type, the first line of the file is always skipped
#' and the columns must be \emph{tab-separated}.
#' See example below for the format of such files generated by the
#' \href{https://druglogics.github.io/druglogics-doc/drabme.html}{Drabme} module.
#'
#' @param file_name string. The name of the file, can be a full path.
#' @param file_type string. The type of input file, can be either \emph{ensemblewise}
#' (default option) or \emph{modelwise}.
#'
#' @return a \code{tibble} containing a representation of the data in the file.
#'
#' @examples
#'
#' ensemblewise_synergies_file = system.file("extdata", "ensemblewise_synergies",
#'   package = "emba", mustWork = TRUE)
#' modelwise_synergies_file = system.file("extdata", "modelwise_synergies",
#'   package = "emba", mustWork = TRUE)
#' data_ensemble = get_synergy_scores(ensemblewise_synergies_file) # file_type = "ensemblewise"
#' data_modelwise = get_synergy_scores(modelwise_synergies_file, file_type = "modelwise")
#'
#' @importFrom readr read_tsv
#' @export
get_synergy_scores = function(file_name, file_type = "ensemblewise") {
  stopifnot(file.exists(file_name))
  type = match.arg(arg = file_type, choices = c("ensemblewise", "modelwise"))

  if (type == "ensemblewise") {
    res = read_tsv(file = file_name, skip = 1, col_types = "cd",
      col_names = c('perturbation', 'score'))
  } else { # model-wise
    res = read_tsv(file = file_name, skip = 1, col_types = "cii",
      col_names = c('perturbation', 'synergies', 'non-synergies'))
  }

  res$perturbation = gsub('\\[|\\]', '', res$perturbation)

  return(res)
}
