#' Create term data for interactive example
#'
#' @param de Differential expression results, use \code{yeast_de} data attached to this package.
#'
#' @return A list of objects containing functional terms for GO, Reactome, KEGG and BioPlanet.
#' @export
#'
#' @examples
#' data(yeast_de)
#' term_data <- fetch_terms_for_example(yeast_de)
fetch_terms_for_example <- function(de) {
  # Binding variables from non-standard evaluation locally
  gene_id <- gene_symbol <- NULL

  # All gene background
  all_genes <- de$gene_id

  # load GO terms
  go <- fetch_go_from_go("sgd")
  go_data <- prepare_for_enrichment(go$terms, go$mapping, all_genes, feature_name = "gene_synonym")

  # load Reactome pathways
  re <- fetch_reactome("Saccharomyces cerevisiae")
  re_data <- prepare_for_enrichment(re$terms, re$mapping, all_genes, feature_name = "gene_id")

  # load KEGG pathways
  kg <- fetch_kegg("sce")
  kg_data <- prepare_for_enrichment(kg$terms, kg$mapping, all_genes, feature_name = "gene_id")

  # load BioPlanet terms
  bp <- fetch_bp()
  # as BP does not use SGD identifiers, we need translate gene symbols to SGD
  # gene ids
  symid <- de |>
    dplyr::select(gene_id, gene_symbol)
  bp$mapping <- bp$mapping |>
    dplyr::left_join(symid, by = "gene_symbol") |>
    tidyr::drop_na()
  bp_data <- prepare_for_enrichment(bp$terms, bp$mapping, all_genes, feature_name = "gene_id")

  # Put all functional term data in one structure; Shiny app will access
  # individual ontologies from this list
  list(
    go = go_data,
    re = re_data,
    kg = kg_data,
    bp = bp_data
  )
}


#' Volcano plot
#'
#' @param d Tibble with x, y and FDR
#' @param fdr_limit FDR limit below which point will be plotted black
#'
#' @return A ggplot object
plot_volcano <- function(d, fdr_limit = 0.05) {
  # Binding variables from non-standard evaluation locally
  x <- y <- FDR <- NULL

  sres <- d |>
    dplyr::filter(FDR <= fdr_limit)
  g <- ggplot(d, aes(x, y)) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 18)
    ) +
    geom_point(size = 0.2) +
    geom_vline(xintercept = 0, colour = "grey70") +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
  if (nrow(sres) > 0) {
    g <- g + geom_hline(yintercept = -log10(max(sres$PValue)), colour = "red",
                        linetype = "dashed", alpha = 0.2)
  }
  g
}

#' MA plot
#'
#' @param d Tibble with x, y and FDR
#' @param fdr_limit FDR limit below which point will be plotted black
#'
#' @return A ggplot object
plot_ma <- function(d, fdr_limit = 0.05) {
  # Binding variables from non-standard evaluation locally
  x <- y <- NULL

  ggplot(d, aes(x, y)) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 18)
    ) +
    geom_point(data = d[d$FDR > fdr_limit,], size = 0.1, colour = "grey50") +
    geom_point(data = d[d$FDR <= fdr_limit,], size = 0.2, colour = "black") +
    geom_hline(yintercept = 0, colour = "grey70") +
    labs(x = expression(log[2]~CPM), y = expression(log[2]~FC))
}



#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param input App input
#'
#' @return A tibble with x, y
get_xy_data <- function(de, input) {
  logFC <- PValue <- logCPM <- NULL
  if (input$plot_type == "volcano") {
    xy_data <- de |>
      dplyr::mutate(x = logFC, y = -log10(PValue))
  } else if (input$plot_type == "ma") {
    xy_data <- de  |>
      dplyr::mutate(x = logCPM, y = logFC)
  }
  xy_data
}

#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param term_data All functional term data
#' @param input App input
#' @param max_points Maximum number of points allowed to select
#'
#' @return A tibble with functional enrichment results
enrichment_table <- function(de, term_data, input, max_points = 10000) {
  # Binding variables from non-standard evaluation locally
  p_adjust <- n_with_sel <- NULL

  xy_data <- get_xy_data(de, input)
  terms <- term_data[[input$ontology]]
  gene2symbol <- rlang::set_names(de$gene_symbol, de$gene_id)
  sel <- NULL
  fe <- NULL
  if (!is.null(input$plot_brush)) {
    brushed <- brushedPoints(xy_data, input$plot_brush)
    sel <- brushed$gene_id
    n <- length(sel)
    if (n > 0 && n <= max_points) {
      fe <- functional_enrichment(de$gene_id, sel, terms, feat2name = gene2symbol)
      if(!is.null(fe))
        fe <- dplyr::filter(fe, p_adjust < 0.05 & n_with_sel > 2)
    } else if (n > 0) {
      fe <- data.frame(Error = paste('only', max_points, 'points can be selected.'))
    }
  }
  fe
}

#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param input App input
#'
#' @return Data for volcano or ma plot
main_plot <- function(de, input) {
  xy_data <- get_xy_data(de, input)
  if (input$plot_type == "volcano") {
    plot_volcano(xy_data)
  } else if (input$plot_type == "ma") {
    plot_ma(xy_data)
  }
}





#' Small Shiny app serving as example for fast enrichment
#'
#' @param de Differential expression results, \code{yeast_de} attached to this
#'   package can be used.
#' @param term_data A list of \code{fenr_terms} objects containing functional
#'   data for various ontologies. \code{fetch_terms_for_example} can be used to
#'   create this object.
#'
#' @return An interactive Shiny app
#' @export
#' @import assertthat
#' @import shiny
#' @import ggplot2
#' @examples
#' \dontrun{
#' data(yeast_de)
#' term_data <- fetch_terms_for_example(yeast_de)
#' enrichment_interactive(yeast_de, term_data)
#' }
enrichment_interactive <- function(de, term_data) {
  ui <- function() {
    shinyUI(fluidPage(
      tags$style("table{font-size: 11px; background-color: #EAF5FF}"),
      titlePanel("fenr example"),
      p("Select a group of genes in the plot to see their functional enrichment."),

      fluidRow(
        column(5,
          radioButtons(
            inputId = "plot_type",
            label = "Plot type:",
            choices = c("Volcano" = "volcano", "MA" = "ma"),
            inline = TRUE
          ),
          plotOutput(
            inputId = "main_plot",
            height = "480px",
            width = "100%",
            brush = "plot_brush",
            hover = "plot_hover"
          )
        ),
        column(7,
          radioButtons(
            inputId = "ontology",
            label = "Ontology:",
            choices = c("GO" = "go", "Reactome" = "re", "KEGG" = "kg", "BioPlanet" = "bp"),
            inline = TRUE
          ),
           div(style = 'height: 480px; overflow-y: scroll', tableOutput("enrichment")),
        )
      )
    ))
  }

  server <- function(input, output, session) {
    # Prevents RStudio from crashing when Shiny window closed manually
    session$onSessionEnded(function() {
      stopApp()
    })

    output$enrichment <- renderTable({
      enrichment_table(de, term_data, input)
    })

    output$main_plot <- renderPlot({
      main_plot(de, input)
    })
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}
