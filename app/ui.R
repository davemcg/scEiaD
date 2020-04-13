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
                                            fluidRow(column(5,
                                                            selectizeInput('Gene', strong('Gene: '),
                                                                           choices=NULL, multiple=FALSE)),
                                                     column(5,
                                                            selectizeInput('pt_size', strong('Point Size: '),
                                                                           choices=c(3,5,10,20, 50), 
                                                                           selected = 5, multiple=FALSE)))),
                                     column(6,
                                            plotOutput('meta_plot'),
                                            fluidRow(column(5, 
                                                            selectInput("label_toggle", label = strong("Label: "), 
                                                                        choices = list("None" = 0,
                                                                                       "CellType (published)" = 1,
                                                                                       "CellType (predict)" = 2,
                                                                                       "Cluster" = 3), multiple = FALSE,
                                                                        selected = 2)),
                                                     column(5,selectizeInput('meta_column', strong('Color: '),
                                                                             choices= NULL, selected = 'CellType_predict')),
                                                     column(5, radioButtons('meta_column_transform', 
                                                                            label = 'Numeric Transform', inline = TRUE,
                                                                            choices = list("None" = "None", "log2" = "log2")))))
                                   ),
                                   fluidRow(
                                     column(6, 
                                            fluidRow( 
                                              column(12, actionButton('BUTTON_draw_scatter',' Draw Scatter Plot', icon = icon("arrow-up"),
                                                                      style='background-color: #3399ff; color: #ffffff'),
                                                     actionButton('BUTTON_make_gene_table',' Make Gene Table', icon = icon("arrow-down"),
                                                                  style='background-color: #3399ff; color: #ffffff'))),
                                            br(),
                                            selectizeInput('grouping_features', strong('Gene Table Grouping(s)'),
                                                           choices = NULL, 
                                                           multiple = TRUE),
                                            div(DT::dataTableOutput('gene_cluster_stats'), style='font-size:75%')),
                                     column(6,
                                            fluidRow( 
                                              column(12,
                                                     actionButton('BUTTON_draw_meta',' Draw Meta Plot', icon = icon("arrow-up"),
                                                                  style='background-color: #3399ff; color: #ffffff'),
                                                     actionButton('BUTTON_make_meta_table',' Make Meta Table', icon = icon("arrow-down"),
                                                                  style='background-color: #3399ff; color: #ffffff'))),
                                            br(),
                                            selectizeInput('meta_groupings', strong('Metadata Table Groupings '),
                                                           choices = NULL, 
                                                           multiple = TRUE),
                                            div(DT::dataTableOutput('metadata_stats'), style='font-size:75%')),
                                   )
                                 ),
                        ),
                        tabPanel('Heatmap',
                                 column(8,
                                        plotOutput('dotplot', height = 800),
                                        selectizeInput('dotplot_Gene', strong('Genes: '),
                                                       choices=NULL, multiple=TRUE),
                                        actionButton('BUTTON_draw_dotplot','Draw Dotplot!', 
                                                     style='background-color: #3399ff; color: #ffffff'),
                                        br()))
             ),
             tabPanel('Overview',
                      (
                        fluidRow(column(width = 8, offset = 1, h2('Ocular Single Cell Gene Expression Anthology v0.00'))))
             )
  )
)
