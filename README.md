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

The first step is to download functional term data. `fenr` package supports downloads from Gene Ontology, Reactome, KEGG, BioPlanet and WikiPathways. Other ontologies can be used as long as they are converted into a suitable format (see function `prepare_for_enrichment` for details). We will download functional terms and gene mapping from BioPlanet.

```
bp <- fetch_bp()
```

This is a list with two tibbles containing term information (`term_id` and `term_name`) and gene-term mapping (`term_id` and `gene_symbol`). We convert it into an object suitable for fast functional enrichment. `exmpl_all` is the attached example of gene background - a vector with gene symbols related to all detections in an experiment.

```
data(exmpl_all, exmpl_sel)
bp_terms <- prepare_for_enrichment(bp$terms, bp$mapping, exmpl_all, feature_name = "gene_symbol")
```

`bp_terms` is a data structure containing all the mappings in quickly accessible form. From this point on, `bm_terms` can be used to do multiple functional enrichments on various gene selections. For example, if `exmpl_all` is a vector with all background gene symbols and `exmpl_sel` is a vector with genes of interest (both attached to the package), functional enrichment can be found using

```
enr <- functional_enrichment(exmpl_all, exmpl_sel, bp_terms)
```

The result is a tibble

```
   term_id        term_name                             N_with n_with_sel n_expect enrichment odds_ratio ids    p_value p_adjust
   <chr>          <chr>                                  <int>      <int>    <dbl>      <dbl>      <dbl> <chr>    <dbl>    <dbl>
 1 bioplanet_686  Interferon alpha/beta signaling           40          5     0.11      44.2       66.9  ADAR… 6.5 e- 8 6.5 e- 7
 2 bioplanet_675  Type II interferon signaling (interf…     21          2     0.06      33.7       41.1  IFIT… 1.54e- 3 6.62e- 3
 3 bioplanet_611  Interferon signaling                     128          8     0.36      22.1       38.5  ADAR… 9.68e-10 1.45e- 8
 4 bioplanet_596  Immune system signaling by interfero…    197         10     0.56      18         36.7  ADAR… 3.24e-11 9.71e-10
 5 bioplanet_1663 ERK1/ERK2 MAPK pathway                    33          2     0.09      21.4       25.2  DUSP… 3.8 e- 3 1.34e- 2
 6 bioplanet_1047 Transport of mature mRNAs derived fr…     34          2     0.1       20.8       24.4  NUP2… 4.03e- 3 1.34e- 2
 7 bioplanet_727  Neurotrophic factor-mediated Trk rec…     37          2     0.1       19.1       22.3  KRAS… 4.77e- 3 1.43e- 2
 8 bioplanet_1458 Antiviral mechanism by interferon-st…     63          3     0.18      16.8       20.6  IFIT… 6.89e- 4 4.13e- 3
 9 bioplanet_664  Interferon-gamma signaling pathway        74          3     0.21      14.3       17.4  HLA-… 1.1 e- 3 5.51e- 3
10 bioplanet_1048 Transport of mature transcript to cy…     51          2     0.14      13.9       15.9  NUP2… 8.91e- 3 2.23e- 2
11 bioplanet_905  Immune system                            634         11     1.79       6.14      12.6  ADAR… 2.18e- 7 1.63e- 6
12 bioplanet_1438 EGFR1 pathway                             70          2     0.2       10.1       11.4  DUSP… 1.63e- 2 3.37e- 2
13 bioplanet_517  Signaling by interleukins                 75          2     0.21       9.43      10.6  KRAS… 1.86e- 2 3.37e- 2
14 bioplanet_1675 Interleukin-1 signaling pathway           76          2     0.21       9.31      10.5  SQST… 1.91e- 2 3.37e- 2
15 bioplanet_1340 ERBB1 downstream pathway                  76          2     0.21       9.31      10.5  DUSP… 1.91e- 2 3.37e- 2
16 bioplanet_530  Signaling by NGF                         143          3     0.4        7.42       8.72 DUSP… 7.16e- 3 1.95e- 2
17 bioplanet_544  NGF signaling via TRKA from the plas…     93          2     0.26       7.61       8.5  DUSP… 2.79e- 2 4.64e- 2
18 bioplanet_1138 Messenger RNA processing                 180          3     0.51       5.9        6.86 ADAR… 1.34e- 2 3.1 e- 2
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
 

