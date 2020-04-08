print('UI Start')
print(Sys.time())

library(shiny)
library(ggplot2)
library(scattermore)
library(pals)
library(shinythemes)
library(cowplot)

shinyUI(
  navbarPage('Ocular scAnthology',
             theme = shinytheme('flatly'),
             selected = 'Overview',
             navbarMenu('Viz',
                        # Expression ---------------
                        tabPanel('UMAP 2D', 
                                 fluidPage(
                                   fluidRow(
                                     column(6,
                                            plotOutput('gene_scatter_plot'),
                                            column(4,
                                                   selectizeInput('Gene', strong('Gene: '),
                                                                  choices=NULL, multiple=FALSE)),
                                            column(3,
                                                   selectizeInput('pt_size', strong('Point Size: '),
                                                                  choices=c(3,5,10,20, 50), multiple=FALSE))),
                                     column(6,
                                            plotOutput('meta_plot'),
                                            shinyWidgets::materialSwitch("label_toggle", label = "Display Cell Type (predicted)", value = TRUE ),
                                            selectizeInput('meta_column', strong('Metadata: '),
                                                           choices= NULL),
                                            radioButtons('meta_column_transform', 
                                                         label = 'Numeric Transform',
                                                         choices = list("None" = "None", "log2" = "log2")))
                                   ),
                                   fluidRow(
                                     column(6, 
                                            actionButton('BUTTON_draw_scatter','(Re)Draw Scatter Plot!', 
                                                         style='background-color: #3399ff; color: #ffffff'),
                                            selectizeInput('grouping_features', strong('Group Table '),
                                                           choices = NULL, 
                                                           multiple = TRUE),
                                            div(DT::dataTableOutput('gene_cluster_stats'), style='font-size:75%')),
                                     column(6,
                                            actionButton('BUTTON_draw_meta','(Re)Draw Meta Plot!', 
                                                         style='background-color: #3399ff; color: #ffffff'))
                                   )
                                 ),
                        ),
                        tabPanel('Heatmap',
                                 column(8,
                                        plotOutput('dotplot', height = 800),
                                        selectizeInput('dotplot_Gene', strong('Genes: '),
                                                       choices=NULL, multiple=TRUE),
                                        actionButton('BUTTON_draw_dotplot','(Re)Draw Dotplot!', 
                                                     style='background-color: #3399ff; color: #ffffff')))
             ),
             tabPanel('Overview',
                      (
                        fluidRow(column(width = 8, offset = 1, h2('Ocular Single Cell Gene Expression Anthology v0.00'))))
             )
  )
)
