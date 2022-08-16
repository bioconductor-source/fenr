ensembl_genes <- function(dataset) {
  mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = dataset)
  biomaRt::getBM(attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id",
    "uniprotswissprot",
    "uniprotsptrembl"
  ), mart = mart) |>
    tibble::as_tibble()
}
