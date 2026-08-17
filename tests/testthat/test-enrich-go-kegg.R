test_that("enrich_go_kegg has the documented public interface", {
  expect_true(is.function(enrich_go_kegg))
  expect_identical(
    names(formals(enrich_go_kegg)),
    c(
      "genes", "species", "p_name", "filter_pvalue", "pvalueCutoff",
      "qvalueCutoff", "kegg_analysis_method", "plot_type",
      "scankeyplot_mode", "fill_colors", "go_kegg_topn",
      "description_wrap_width"
    )
  )
})

test_that("invalid inputs fail before enrichment is attempted", {
  expect_error(enrich_go_kegg(character(), "9606"), "non-empty")
  expect_error(enrich_go_kegg("TP53", "rat"), "9606")
  expect_error(enrich_go_kegg("TP53", "9606", p_name = "padj"), "p.value")
  expect_error(enrich_go_kegg("TP53", "9606", plot_type = "treeplot"),
               "plot_type")
})
