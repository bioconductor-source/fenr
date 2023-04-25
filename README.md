# Fast functional enrichment

Maintainer: Marek Gierlinski (<M.Gierlinski@dundee.ac.uk>)

A simple R package for fast functional enrichment. It is optimised for speed and designed to be used in interactive applications, e.g. Shiny apps.

## Quick example

A Shiny app demonstration `fenr` can be ran directly from GitHub, with the following command in R console:

```
shiny::runGitHub("bartongroup/fenr-shiny-example")
```

## Installation

`fenr` can be installed from GitHub (you need to install `remotes` first).

```
remotes::install_github("bartongroup/fenr", build_vignettes = TRUE)
```

## Usage

The first step is to download functional term data. `fenr` package supports downloads from Gene Ontology, Reactome, KEGG, BioPlanet and WikiPathways. Other ontologies can be used as long as they are converted into a suitable format (see function `prepare_for_enrichment` for details). We will download functional terms and gene mapping from GO.

```
go <- fetch_go(species = "sgd")
```

This is a list with two tibbles containing term information (`term_id` and `term_name`) and gene-term mapping (`term_id` and `gene_symbol`). We convert it into an object suitable for fast functional enrichment. `exmpl_all` is the attached example of gene background - a vector with gene symbols related to all detections in an experiment.

```
data(exmpl_all, exmpl_sel)
go_terms <- prepare_for_enrichment(go$terms, go$mapping, exmpl_all, feature_name = "gene_symbol")
```

`go_terms` is a data structure containing all the mappings in quickly accessible form. From this point on, `go_terms` can be used to do multiple functional enrichments on various gene selections. For example, if `exmpl_all` is a vector with all background gene symbols and `exmpl_sel` is a vector with genes of interest (both attached to the package), functional enrichment can be found using

```
enr <- functional_enrichment(exmpl_all, exmpl_sel, go_terms)
```

The result is a tibble:

```
# A tibble: 51 × 10
  term_id    term_name                                N_with n_with_sel n_expect enrichment odds_ratio ids    p_value p_adjust
   <chr>      <chr>                                     <int>      <int>    <dbl>      <dbl>      <dbl> <chr>    <dbl>    <dbl>
 1 GO:1905356 regulation of snRNA pseudouridine synth…      2          2     0.01      333        Inf   TOR2… 8.61e- 6 5.49e- 5
 2 GO:0031929 TOR signaling                                19         18     0.06      315      41800   TOR2… 0        0       
 3 GO:0031931 TORC1 complex                                11         10     0.03      302       6330   TOR2… 0        0       
 4 GO:0001558 regulation of cell growth                     9          5     0.03      185        544   KOG1… 1.84e-11 2.34e-10
 5 GO:0031932 TORC2 complex                                15          7     0.05      155        435   TOR2… 4.66e-15 7.93e-14
 6 GO:0016242 negative regulation of macroautophagy         9          3     0.03      111        193   KSP1… 1.95e- 6 1.42e- 5
 7 GO:0043666 regulation of phosphoprotein phosphatas…     13          4     0.04      102        182   SAP1… 4.24e- 8 4.33e- 7
 8 GO:0030950 establishment or maintenance of actin c…     15          4     0.05       88.7      149   TOR2… 8.07e- 8 6.86e- 7
 9 GO:0031930 mitochondria-nucleus signaling pathway        9          2     0.03       73.9      105   TOR1… 3.06e- 4 1.3 e- 3
10 GO:0010507 negative regulation of autophagy             12          2     0.04       55.4       73.2 TOR2… 5.58e- 4 2.03e- 3
# ℹ 41 more rows
```

The columns are as follows

 - `N_with` - number of features (genes) with this term in the background of all genes.
 - `n_with_sel` - number of features with this term in the selection.
 - `n_expect` - expected number of features with this term under the null hypothesis (terms are randomly distributed).
 - `enrichment` - ratio of observed to expected.
 - `odds_ratio` - effect size, odds ratio from the contingency table.
 - `ids` - identifiers of features with term in the selection.
 - `p_value` - raw p-value from hypergeometric distribution.
 - `p_adjust` - p-value adjusted for multiple tests using Benjamini-Hochberg approach.
 

