context("Testing printing functions")
test_that("they return correct output", {
  expect_output(print_model_and_drug_stats(3, 150, 42, FALSE),
    "> Drug combinations tested: 3\nNumber of models: 150\nNumber of nodes: 42")

  expect_output(print_biomarkers_per_predicted_synergy(biomarkers.dir = tempdir(),
    predicted.synergies = c("A-B", "C-D"), html.output = FALSE),
    "> Biomarkers for A-B synergy prediction\\n\\nActive biomarkers\\n0 nodes: \\n\\nInhibited biomarkers\\n0 nodes: \\n\\nBiomarkers for C-D synergy prediction\\n\\nActive biomarkers\\n0 nodes: \\n\\nInhibited biomarkers\\n0 nodes: \\n")

  biomarkers.dir = system.file("extdata", "biomarkers", package = "emba", mustWork = TRUE)
  expect_output(print_biomarkers_per_predicted_synergy(biomarkers.dir,
    predicted.synergies = c("A-B", "C-D"), html.output = FALSE),
    regexp = "> Biomarkers for A-B synergy prediction\\n\\nActive biomarkers\\n3 nodes: A1, A2, B1\\n\\nInhibited biomarkers\\n3 nodes: A3, A4, B3\\n\\nBiomarkers for C-D synergy prediction\\n\\nActive biomarkers\\n0 nodes: \\n\\nInhibited biomarkers\\n0 nodes: \\n")
})
