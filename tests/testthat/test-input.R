library(dplyr)

context("Testing 'get_alt_drugname'")
test_that("it returns proper results", {
  expect_equal(get_alt_drugname("A-B"), "B-A")
  expect_equal(get_alt_drugname("BB-AA"), "AA-BB")
})

context("Testing 'is_comb_element_of'")
test_that("it returns proper results", {
  expect_true(is_comb_element_of("A-B", c("E-F", "A-B")))
  expect_true(is_comb_element_of("B-A", c("E-F", "A-B")))

  expect_false(is_comb_element_of("A-B", c("E-F", "A-D")))
  expect_false(is_comb_element_of("A-B", c()))
})

context("Testing 'get_model_predictions'")
test_that("it returns proper results", {
  model.predictions.file = system.file("extdata", "model_predictions", package = "emba", mustWork = TRUE)
  model.predictions = get_model_predictions(model.predictions.file)

  expect_equal(rownames(model.predictions)[1], "topology.sif_run_961__G19_M5.gitsbe")
  expect_equal(rownames(model.predictions)[2], "topology.sif_run_2402__G19_M9.gitsbe")

  expect_equal(colnames(model.predictions)[1], "5Z-AK")
  expect_equal(colnames(model.predictions)[2], "5Z-BI")

  expect_equal(model.predictions[1,39], 1)
  expect_equal(model.predictions[1,42], 0)
  expect_true(is.na(model.predictions[1,45]))
})

context("Testing 'get_observed_synergies'")
test_that("it returns proper results", {
  observed.synergies.file = system.file("extdata", "observed_synergies",
    package = "emba", mustWork = TRUE)
  obs.synergies = get_observed_synergies(file = observed.synergies.file,
    drug.combinations.tested = c("A-B", "C-A", "D-C"))

  expect_equal(length(obs.synergies), 2)
  expect_equal(obs.synergies, c("A-B", "A-C"))
})

context("Testing 'get_stable_state_from_models_dir'")
test_that("it returns proper results", {
  models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
  models.ss = get_stable_state_from_models_dir(models.dir)

  # `test4` has 2 stable states and is not included
  expect_equal(sort(rownames(models.ss)), c("test", "test2", "test3"))
  expect_equal(ncol(models.ss), 139)
  expect_equal(colnames(models.ss)[1], "MAP3K7")
  expect_equal(colnames(models.ss)[139], "CASP9")
  expect_equal(models.ss["test", "CASP9"], 0)
  expect_equal(models.ss["test2", "CASP9"], 1)
  expect_equal(models.ss["test3", "CASP9"], 1)

  # now `test4` will be included in the returned result
  models.ss.all = get_stable_state_from_models_dir(models.dir, all.ss = TRUE)
  expect_equal(models.ss.all %>% pull(model_name) %>% sort(),
    c("test", "test2", "test3", "test4", "test4"))
  expect_equal(ncol(models.ss.all), 140) # one extra column for the model names
  expect_equal(colnames(models.ss.all)[1], "MAP3K7")
  expect_equal(colnames(models.ss.all)[139], "CASP9")
  expect_equal(models.ss.all %>% filter(model_name == "test") %>% pull(CASP9), 0)
  expect_equal(models.ss.all %>% filter(model_name == "test2") %>% pull(CASP9), 1)
  expect_equal(models.ss.all %>% filter(model_name == "test3") %>% pull(CASP9), 1)
  expect_equal(models.ss.all %>% filter(model_name == "test4") %>% pull(CASP9), c(1 ,0))
})

context("Testing 'get_link_operators_from_models_dir'")
test_that("it returns proper results", {
  models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
  models.link = get_link_operators_from_models_dir(models.dir)
  models.link.extra = get_link_operators_from_models_dir(models.dir, remove.equations.without.link.operator = FALSE)

  expect_equal(sort(rownames(models.link)), c("test", "test2", "test3", "test4"))
  expect_equal(sort(rownames(models.link.extra)), c("test", "test2", "test3", "test4"))
  expect_equal(ncol(models.link), 43)
  expect_equal(ncol(models.link.extra), 139)
  expect_equal(colnames(models.link)[1], "MAP3K4")
  expect_equal(colnames(models.link)[43], "CASP9")
  expect_equal(models.link["test", "CASP9"], 0)
  expect_equal(models.link["test2", "CASP9"], 1)
  expect_equal(models.link["test3", "CASP9"], 1)
  expect_equal(models.link["test4", "CASP9"], 1)
  expect_equal(models.link.extra["test", "MAP3K7"], 0.5)
})

context("Testing 'get_fitness_from_models_dir'")
test_that("it returns proper results", {
  models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
  models.fitness = get_fitness_from_models_dir(models.dir)

  expect_equal(sort(names(models.fitness)), c("test", "test2", "test3", "test4"))
  expect_equal(unname(models.fitness["test"]), 0.640625)
  expect_equal(unname(models.fitness["test2"]), 0.625)
  expect_equal(unname(models.fitness["test3"]), 0.6640625)
})

context("Testing 'get_node_names'")
test_that("it returns proper results", {
  models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
  nodes = get_node_names(models.dir)

  expect_equal(length(nodes), 139)
  expect_equal(nodes[1], "MAP3K7")
  expect_equal(nodes[139], "CASP9")
})

context("Testing 'get_model_names'")
test_that("it returns proper results", {
  models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
  models = get_model_names(models.dir)

  expect_equal(length(models), 4)
  expect_equal(sort(models), c("test", "test2", "test3", "test4"))
})

context("Testing 'assign_link_operator_value_to_equation'")
test_that("it returns proper results", {
  expect_equal(assign_link_operator_value_to_equation("A and not B"), 0)
  expect_equal(assign_link_operator_value_to_equation("C or D or not B"), 1)

  expect_true(is.na(assign_link_operator_value_to_equation("not B")))
  expect_true(is.na(assign_link_operator_value_to_equation("A AND NOT B")))
  expect_true(is.na(assign_link_operator_value_to_equation("C or D OR not B")))
})

context("Testing 'validate_observed_synergies_data'")
test_that("it returns proper results", {
  expect_success(expect_null(validate_observed_synergies_data(observed.synergies = c("A-B", "B-C"),
    drug.combinations.tested = c("B-C", "A-B"))))
  expect_success(expect_null(validate_observed_synergies_data(observed.synergies = c("A-B", "B-C"),
    drug.combinations.tested = c("B-A", "D-E", "C-B"))))
  expect_error(validate_observed_synergies_data(observed.synergies = c("A-B", "B-C"),
    drug.combinations.tested = c("A-B")))
  expect_error(validate_observed_synergies_data(observed.synergies = c(""),
    drug.combinations.tested = c("A-B")))
})

context("Testing 'construct_network', 'filter_network' and 'get_neighbors'")
test_that("it returns proper results", {
  topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
  net = construct_network(topology.file)
  net0 = filter_network(net, nodes = c("A", "D"), level = 0)
  net1 = filter_network(net, nodes = c("A", "D"), level = 1)
  net2 = filter_network(net, nodes = c("A", "D"), level = 2)

  expect_equal(length(V(net)), 10)
  expect_equal(length(E(net)), 11)

  expect_equal(length(V(net0)), 2)
  expect_equal(length(E(net0)), 0)

  expect_equal(length(V(net1)), 5)
  expect_equal(length(E(net1)), 6)

  expect_equal(length(V(net2)), 9)
  expect_equal(length(E(net2)), 10)
})

context("Testing 'get_edges_from_topology_file'")
test_that("it returns proper results", {
  topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
  edges = get_edges_from_topology_file(topology.file)

  expect_equal(dim(edges), c(11, 4))
  expect_equal(colnames(edges), c("source", "target", "regulation.effect", "color"))
})

context("Testing 'get_observed_synergies_per_cell_line'")
test_that("it returns proper results", {
  data.dir = system.file("extdata", package = "emba", mustWork = TRUE)
  res = get_observed_synergies_per_cell_line(cell.line.dirs = c(data.dir, data.dir),
    drug.combos = c("A-C", "D-R", "A-B", "A-R"))

  expect_equal(dim(res), c(2,4))
  expect_equal(rownames(res), c("extdata", "extdata1"))
  expect_equal(colnames(res), c("A-C", "D-R", "A-B", "A-R"))
  expect_equal(sum(res), 4)
})

context("Testing 'get_synergy_scores'")
test_that("it returns proper results", {
  ew_synergies_file = system.file("extdata", "ensemblewise_synergies", package = "emba", mustWork = TRUE)
  mw_synergies_file = system.file("extdata", "modelwise_synergies", package = "emba", mustWork = TRUE)

  ew_syn_scores = get_synergy_scores(ew_synergies_file) # file_type = "ensemblewise"
  mw_syn_scores = get_synergy_scores(mw_synergies_file, file_type = "modelwise")

  ew_scores = ew_syn_scores %>% pull(score)
  mw_synergies = mw_syn_scores %>% pull(synergies)
  mw_non_synergies = mw_syn_scores %>% pull(`non-synergies`)

  expect_equal(dim(ew_syn_scores), c(4,2))
  expect_equal(dim(mw_syn_scores), c(4,3))

  expect_equal(colnames(ew_syn_scores), c("perturbation","score"))
  expect_equal(colnames(mw_syn_scores), c("perturbation","synergies", "non-synergies"))
  expect_equal(ew_syn_scores %>% pull(perturbation), c("A-B", "B-C", "A-C", "A-D"))
  expect_equal(mw_syn_scores %>% pull(perturbation), c("A-B", "B-C", "A-C", "A-D"))
  expect_equal(ew_scores[1], -0.6336996336996337)
  expect_equal(mw_synergies[1], 15)
  expect_equal(mw_non_synergies[1], 78)
})
