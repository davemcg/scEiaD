# server.R
time <- Sys.time()
cat(file = stderr(), 'Server Go!\n')
#options(shiny.trace=TRUE)
options(shiny.sanitize.errors = FALSE)

library(ggplot2)
library(Cairo)
library(scattermore)
library(pals)
library(ggrepel)
library(patchwork)
library(purrr)
library(pool)
library(RSQLite)
library(dplyr)
library(formattable)
anthology_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-2000-counts-onlyDROPLET-batch-scVI-6-0.1-500-10.sqlite", idleTimeout = 3600000)

# fancy tables
# they come from `tables.Rmd` in analyis/
load('www/formattables.Rdata')
# filter
meta_filter <- left_join(anthology_2020_v01 %>% tbl('metadata_filter'), 
                         anthology_2020_v01 %>% tbl('doublets'), by ='Barcode') %>% 
  as_tibble() %>% 
  mutate(`Doublet Probability` = as.numeric(`Doublet Probability`),
         doublet_score_scran = as.numeric(doublet_score_scran))

# cutdown mf for plotting
mf <- meta_filter %>% sample_frac(0.2)
# get coords for cell labels
celltype_predict_labels <- meta_filter %>% 
  group_by(CellType_predict) %>% 
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
celltype_labels <- meta_filter %>% 
  group_by(CellType) %>% 
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
# get coords for cell labels
cluster_labels <- meta_filter %>% group_by(cluster) %>% summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))


# # attach colors to cell types
# cell_types <- meta_filter %>% 
#   pull(CellType_predict) %>% unique() %>% sort()
# type_val <- setNames(pals::alphabet(n = cell_types %>% length()), cell_types)
# type_col <- scale_colour_manual(values = type_val)
# type_fill <- scale_fill_manual(values = type_val)
cat(file=stderr(), 'Data loaded in ')
cat(file=stderr(), Sys.time() - time)
cat(file=stderr(), ' seconds.\n')


# site begins! ---------
shinyServer(function(input, output, session) {
  observe({
    query <- parseQueryString(session$clientData$url_search)
    ## server help queries ------
    # gene plot updateSelectizeInput -------
    if (is.null(query[['Gene']])){
      updateSelectizeInput(session, 'Gene',
                           choices = anthology_2020_v01 %>% tbl('genes') %>% as_tibble() %>% pull(1),
                           options = list(placeholder = 'Type to search'), 
                           selected = 'CRX',
                           server = TRUE)
    }
    # gene plot category filtering ----
    if (is.null(query[['gene_filter_cat']])){
      updateSelectizeInput(session, 'gene_filter_cat',
                           choices = meta_filter %>% 
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           selected = '',
                           server = TRUE)
    }
    
    observeEvent(input$gene_filter_cat, {
      if (input$gene_filter_cat == ''){
        choice = ''
      } else {
        choice = meta_filter[,input$gene_filter_cat] %>% pull(1) %>% unique() %>% sort()
      }
      updateSelectizeInput(session, 'gene_filter_on',
                           choices = choice,
                           server = TRUE)
    })
    
    
    # meta plot updateSelectizeInput ------
    if (is.null(query[['meta_column']])){
      updateSelectizeInput(session, 'meta_column',
                           choices = meta_filter %>% 
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CellType_predict',
                           server = TRUE)
    }
    # meta plot category filtering ----
    if (is.null(query[['meta_filter_cat']])){
      updateSelectizeInput(session, 'meta_filter_cat',
                           choices = meta_filter %>% 
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           selected = '',
                           server = TRUE)
    }
    
    observeEvent(input$meta_filter_cat, {
      if (input$meta_filter_cat == ''){
        choice = ''
      } else {
        choice = meta_filter[,input$meta_filter_cat] %>% pull(1) %>% unique() %>% sort()
      }
      updateSelectizeInput(session, 'meta_filter_on',
                           choices = choice,
                           server = TRUE)
    })
    
    # dotplot updateSelectizeInput ----
    if (is.null(query[['dotplot_Gene']])){
      updateSelectizeInput(session, 'dotplot_Gene',
                           choices = anthology_2020_v01 %>% tbl('genes') %>% as_tibble() %>% pull(1),
                           options = list(placeholder = 'Type to search'), 
                           selected = c('RHO','WIF1','CABP5', 'AIF1','AQPT4','ARR3','ONECUT1','GRIK1','GAD1','POU4F2'),
                           server = TRUE)
    }
    if (is.null(query[['dotplot_groups']])){
      updateSelectizeInput(session, 'dotplot_groups',
                           choices = meta_filter %>% 
                             select(-Barcode) %>% 
                             select_if(purrr::negate(is.numeric)) %>% 
                             colnames() %>% sort(),
                           options = list(placeholder = 'Type to search',
                                          maxItems = 2), 
                           selected = c('CellType_predict','organism'),
                           server = TRUE)
    }
    
    
    # 
    if (is.null(query[['grouping_features']])){
      updateSelectizeInput(session, 'grouping_features',
                           choices = anthology_2020_v01 %>% tbl('grouped_stats') %>% 
                             select(-Gene, -cell_ct, -cell_exp_ct, -cpm) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'), 
                           selected = c('CellType_predict'),
                           server = TRUE)
    }
    
    if (is.null(query[['meta_groupings']])){
      updateSelectizeInput(session, 'meta_groupings',
                           choices = meta_filter %>% 
                             select(-Barcode, -UMAP_1, -UMAP_2, -nCount_RNA, -nFeature_RNA, -percent_mt) %>% 
                             colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'), 
                           selected = c('CellType_predict', 'organism'),
                           server = TRUE)
    }
    # temporal plot updateSelect -----
    if (is.null(query[['temporal_gene']])){
      updateSelectizeInput(session, 'temporal_gene',
                           choices = anthology_2020_v01 %>% tbl('genes') %>% as_tibble() %>% pull(1),
                           options = list(placeholder = 'Type to search'), 
                           selected = c('PAX6','POU4F2'),
                           server = TRUE)
    }    
    
    
    # facet plot updateSelect --------
    if (is.null(query[['facet']])){
      updateSelectizeInput(session, 'facet',
                           choices = meta_filter %>% 
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = 'organism',
                           server = TRUE)
    }
    
    if (is.null(query[['facet_color']])){
      updateSelectizeInput(session, 'facet_color',
                           choices = meta_filter %>% 
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CellType',
                           server = TRUE)
    }
    # diff table updateSelect ------
    if (is.null(query[['diff_gene']])){
      updateSelectizeInput(session, 'diff_gene',
                           choices = anthology_2020_v01 %>% tbl('genes') %>% as_tibble() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CRX',
                           server = TRUE)
    }
    if (is.null(query[['diff_term']])){
      updateSelectizeInput(session, 'diff_term',
                           choices = anthology_2020_v01 %>% tbl(diff_table()) %>%
                             as_tibble() %>% pull(term) %>% unique() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    
    # BREAK -------
    # gene scatter plot ------------
    gene_scatter_ranges <- reactiveValues(x = c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max()), 
                                          y = c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max()))
    
    gene_scatter_plot <- eventReactive(input$BUTTON_draw_scatter, {
      cat(file=stderr(), paste0(Sys.time(), ' Gene Scatter Plot Call\n'))
      gene <- input$Gene
      pt_size <- input$pt_size_gene %>% as.numeric()
      expression_range <- input$gene_scatter_slider
      p <-  anthology_2020_v01 %>% tbl('cpm') %>% 
        filter(Gene == gene) %>% 
        as_tibble() %>% 
        mutate(cpm = cpm - min(cpm) + 1) %>% 
        filter(cpm > as.numeric(expression_range[1]), 
               cpm < as.numeric(expression_range[2])) %>% 
        left_join(., meta_filter, by = 'Barcode') %>% 
        filter(!is.na(UMAP_1), !is.na(UMAP_2), !is.na(cpm)) 
      if (input$gene_filter_cat != ''){
        p <- p %>% 
          filter(!!as.symbol(input$gene_filter_cat) %in% input$gene_filter_on)
      }
      color_range <- range(p$cpm)
      plot <- p %>% ggplot() + 
        geom_scattermost(cbind(mf$UMAP_1, mf$UMAP_2), color = '#D3D3D333', 
                         pointsize = pt_size * 1.5,
                         pixels=c(750,750)) +
        geom_scattermost(cbind(p$UMAP_1, p$UMAP_2),
                         color = viridis::magma(100, alpha=0.3)
                         [1+99*(p$cpm-color_range[1])/diff(color_range)],
                         pointsize= pt_size,
                         pixels=c(750,750),
                         interpolate=FALSE) + 
        geom_point(data=data.frame(x=double(0)), aes(x,x,color=x)) +
        scale_color_gradientn(  #add the manual guide for the empty aes
          limits=c(min(p$cpm),max(p$cpm)),
          colors=viridis::magma(100),
          name="log2(cpm+1)") +
        theme_cowplot() + 
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank()) +
        annotate("text", -Inf, Inf, label = paste0(gene, ' expression'), hjust = 0, vjust = 1, size = 6)
      
      suppressWarnings(plot)
      
    })
    
    observeEvent(input$gene_scatter_plot_dblclick, {
      brush <- input$gene_scatter_plot_brush
      if (!is.null(brush)) {
        gene_scatter_ranges$x <- c(brush$xmin, brush$xmax)
        gene_scatter_ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        gene_scatter_ranges$x <- c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max())
        gene_scatter_ranges$y <- c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max())
      }
    })
    output$gene_scatter_plot <- renderPlot({
      gene_scatter_plot() + coord_cartesian(xlim = gene_scatter_ranges$x, ylim = gene_scatter_ranges$y)
    })  
    
    
    # metadata plot --------------
    meta_plot <- eventReactive(input$BUTTON_draw_meta, {
      cat(file=stderr(), paste0(Sys.time(), ' Meta Plot Call\n'))
      meta_column <- input$meta_column
      transform <- input$meta_column_transform
      pt_size <- input$pt_size_meta %>% as.numeric() 
      filter_column <- input$meta_column
      if (transform == 'log2' && is.numeric(meta_filter[,meta_column] %>% pull(1))){
        cat('log2 time')
        meta_filter[,meta_column] <- log2(meta_filter[,meta_column] + 1)
      }
      p_data <- meta_filter %>% 
        filter(!is.na(!!as.symbol(meta_column))) 
      # category filtering
      if (input$meta_filter_cat != ''){
        p_data <- p_data %>% 
          filter(!!as.symbol(input$meta_filter_cat) %in% input$meta_filter_on)
      }
      
      # metadata NUMERIC plot --------------
      if (is.numeric(meta_filter[,meta_column] %>% pull(1)) ){
        color_range <- range(p_data[,meta_column] %>% pull(1))
        suppressWarnings(plot <- ggplot() + 
                           geom_scattermost(cbind(mf %>% 
                                                    filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_1),
                                                  mf %>% 
                                                    filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_2)),
                                            pointsize = pt_size, color = '#D3D3D333',
                                            pixels = c(750,750)) +
                           geom_scattermost(cbind(p_data$UMAP_1, p_data$UMAP_2),
                                            color = viridis::viridis(100, alpha=0.3)
                                            [1+99*((p_data[,meta_column] %>% pull(1))-color_range[1])/diff(color_range)],
                                            pointsize= pt_size,
                                            pixels=c(750,750),
                                            interpolate=FALSE) +
                           geom_point(data=data.frame(x=double(0)), aes(x,x,color=x))  + 
                           scale_color_gradientn(  #add the manual guide for the empty aes
                             limits=c(min(p_data[,meta_column] %>% pull(1)),
                                      max(p_data[,meta_column] %>% pull(1))),
                             colors=viridis::viridis(100),
                             name=meta_column) +
                           guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                           theme_cowplot() + 
                           theme(axis.line = element_blank(),
                                 axis.title = element_blank(),
                                 axis.ticks = element_blank(),
                                 axis.text = element_blank()) +
                           annotate("text", -Inf, Inf, label = "Metadata", hjust = 0, vjust = 1, size = 6))
        # metadata CATEGORICAL plot --------------
      } else {
        group <- p_data[,meta_column] %>% pull(1) %>% as.factor()
        p_color = rep(c(pals::alphabet(),pals::alphabet2()), times = 20)[group] %>% paste0(., '33') # alpha 0.33
        color_data <- group %>% levels() %>% tibble::enframe() %>% mutate(x=0) 
        suppressWarnings(plot <- ggplot() +
                           geom_scattermost(cbind(mf %>% 
                                                    filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_1),
                                                  mf %>% 
                                                    filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_2)),
                                            pointsize = pt_size, color = '#D3D3D333',
                                            pixels = c(750,750)) +
                           geom_scattermost(cbind(p_data$UMAP_1, p_data$UMAP_2),
                                            color = p_color ,
                                            pointsize= pt_size,
                                            pixels=c(750,750),
                                            interpolate=FALSE) +
                           geom_point(data=color_data, aes(x,x,color=value), alpha = 0) + 
                           scale_colour_manual(name= meta_column,
                                               values = rep(c(pals::alphabet() %>% unname(),
                                                              pals::alphabet2() %>% unname()), 
                                                            times = 20)) +
                           guides(colour = guide_legend(override.aes = list(alpha = 1, size = 7))) +
                           theme_cowplot() + 
                           theme(axis.line = element_blank(),
                                 axis.title = element_blank(),
                                 axis.ticks = element_blank(),
                                 axis.text = element_blank()) +
                           annotate("text", -Inf, Inf, label = "Metadata", hjust = 0, vjust = 1, size = 6))
      }
      
      more <- NULL
      if ('1' %in% input$label_toggle){
        more <- geom_text_repel(data = celltype_labels, bg.color = 'white',
                                aes(x = UMAP_1, y = UMAP_2, label = CellType)) 
      }
      if ('2' %in% input$label_toggle){
        more <- geom_text_repel(data = celltype_predict_labels, bg.color = 'white',
                                aes(x = UMAP_1, y = UMAP_2, label = CellType_predict)) 
      }
      if ('3' %in% input$label_toggle){
        more <- geom_text_repel(data = cluster_labels, bg.color = 'white',
                                aes(x = UMAP_1, y = UMAP_2, label = cluster),
                                max.iter = 20) 
      } 
      if (meta_column %in% c('cluster','subcluster')){
        suppressWarnings(plot + more + theme(legend.position = 'none'))
      } else {
        suppressWarnings(plot + more)
      }
      
    })
    meta_ranges <- reactiveValues(x = c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max()), 
                                  y = c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max()))
    observeEvent(input$meta_plot_dblclick, {
      brush <- input$meta_plot_brush
      if (!is.null(brush)) {
        meta_ranges$x <- c(brush$xmin, brush$xmax)
        meta_ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        meta_ranges$x <- c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max())
        meta_ranges$y <- c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max())
      }
    })
    
    output$meta_plot <- renderPlot({
      meta_plot() + coord_cartesian(xlim = meta_ranges$x, ylim = meta_ranges$y)
    })  
    
    # gene cluster table  --------
    gene_cluster_stats_maker <- eventReactive(input$BUTTON_make_gene_table, {
      grouping_features <- input$grouping_features
      gene <- input$Gene
      table <- anthology_2020_v01 %>% tbl('grouped_stats') %>% 
        filter(Gene == gene) %>%
        group_by_at(vars(one_of(c('Gene', grouping_features)))) %>% 
        summarise(cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE),
                  cpm = mean(cpm)) %>% 
        as_tibble() %>% 
        tidyr::drop_na() %>% 
        full_join(., meta_filter %>% 
                    group_by_at(vars(one_of(grouping_features))) %>% 
                    summarise(Count = n())) %>% 
        filter(!is.na(Count)) %>% 
        mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
               Expression = round(cpm * (`%` / 100), 2)) %>% 
        select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>% 
        arrange(-Expression) %>% 
        rename(`Cells # Detected` = cell_exp_ct, 
               `Total Cells` = Count,
               `log2(cpm+1)` = Expression)
      
      table %>% DT::datatable(extensions = 'Buttons', rownames = F, 
                              filter = list(position = 'bottom', clear = FALSE),
                              options = list(pageLength = 10, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
    })
    output$gene_cluster_stats <- DT::renderDataTable({ gene_cluster_stats_maker()})
    
    
    # metadata table -----
    metadata_stats <- eventReactive(input$BUTTON_make_meta_table, {
      grouping_features <- input$meta_groupings
      table <- meta_filter %>% 
        group_by_at(vars(one_of(grouping_features))) %>% 
        summarise(Count = n()) %>% 
        select_at(vars(one_of(c(grouping_features, 'Count')))) %>% 
        arrange(-Count) %>% 
        rename(`Total Cells` = Count)
      
      table %>% DT::datatable(extensions = 'Buttons', rownames = F, 
                              filter = list(position = 'bottom', clear = FALSE),
                              options = list(pageLength = 10, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
    })
    output$metadata_stats <- DT::renderDataTable({ metadata_stats()})
    
    # facet plot -----------
    facet_plot <- eventReactive(input$BUTTON_draw_filter, {
      cat(file=stderr(), paste0(Sys.time(), ' Facet Plot Call\n'))
      facet_column <- input$facet
      color_column <- input$facet_color
      #transform <- input$facet_column_transform
      pt_size <- input$pt_size_facet %>% as.numeric() 
      
      gray_data <- meta_filter %>% 
        filter(is.na(!!as.symbol(color_column))) 
      p_data <- meta_filter %>% 
        filter(!is.na(!!as.symbol(facet_column)),
               !is.na(!!as.symbol(color_column)))
      
      suppressWarnings(plot <- ggplot(data = p_data) +
                         geom_scattermore(data = gray_data,
                                          aes(x = UMAP_1, y = UMAP_2),
                                          color = 'gray',
                                          pointsize = pt_size,
                                          pixels = c(750,750),
                                          alpha = 0.4) +
                         geom_scattermore(aes(x = UMAP_1, y = UMAP_2,
                                              color = !!as.symbol(color_column)) ,
                                          pointsize= pt_size,
                                          pixels = c(750,750),
                                          alpha = 0.6) +
                         facet_wrap(vars(!!(as.symbol(facet_column)))) +
                         scale_colour_manual(values = rep(c(pals::alphabet() %>% unname(),
                                                            pals::alphabet2() %>% unname()), 
                                                          times = 20),
                                             na.value = 'gray') +
                         guides(colour = guide_legend(override.aes = list(alpha = 1, size = 7))) +
                         theme_cowplot() + 
                         theme(axis.line = element_blank(),
                               axis.title = element_blank(),
                               axis.ticks = element_blank(),
                               axis.text = element_blank())
      )
      plot
      
    })
    
    output$facet_plot <- renderPlot({
      facet_plot()
    }, height = eventReactive(input$BUTTON_draw_filter, {input$facet_height %>% as.numeric()}))
    
    ## temporal plot -----------
    temporal_plot <- eventReactive(input$BUTTON_draw_temporal, {
      cat(file=stderr(), paste0(Sys.time(), ' Temporal Plot Call\n'))
      gene <- input$temporal_gene
      grouping <- input$temporal_group
      cat(gene)
      cat(grouping)
      if (grouping == 'CellType (predict)'){grouping <- 'CellType_predict'}
      
      temporal_data <- anthology_2020_v01 %>% tbl('cpm') %>%
        filter(Gene %in% gene) %>%
        as_tibble() %>%
        mutate(cpm = cpm - min(cpm) + 1) %>%
        left_join(., meta_filter) %>%
        filter(!is.na(!!as.symbol(grouping)), !grepl('Doub|RPE|Astro', !!as.symbol(grouping))) %>%
        group_by(organism, !!as.symbol(grouping), Age, Gene) %>%
        summarise(cpm = mean(cpm), count = n()) %>%
        ungroup() %>%
        mutate(Age = case_when(organism == 'Mus musculus' & Age == 1000 ~ 20,
                               organism == 'Homo sapiens' & Age == 1000 ~ 65,
                               is.na(Age) ~ 65,
                               Age == 31360 ~ 65,
                               TRUE ~ Age))
      cat(colnames(temporal_data))
      
      temporal_human <-  temporal_data %>% 
        filter(organism == 'Homo sapiens') %>% 
        ggplot(aes(x=Age, y = cpm, color = Gene)) + 
        geom_point(stat = 'identity') +
        ggtitle('Human') +
        geom_line() + ylab('Mean CPM') +
        cowplot::theme_cowplot() + xlab('Age (days from birth)') +
        facet_wrap(vars(!!as.symbol(grouping))) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_colour_manual(values = rep(c(pals::alphabet() %>% unname())))
      temporal_mouse <- temporal_data %>%
        filter(organism == 'Mus musculus') %>%
        ggplot(aes(x=Age, y = cpm, color = Gene)) +
        geom_point(stat = 'identity') +
        ggtitle('Mouse') +
        geom_line() + ylab('Mean CPM') +
        cowplot::theme_cowplot() + xlab('Age (days from birth)') +
        facet_wrap(vars(!!as.symbol(grouping))) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_colour_manual(values = rep(c(pals::alphabet() %>% unname())))
      # draw plot
      temporal_human + temporal_mouse + plot_layout(ncol = 1)
    })
    
    output$temporal_plot <- renderPlot({
      temporal_plot()
    }, height = 700)
    
    ## dotplot ---------
    make_dotplot <- eventReactive(input$BUTTON_draw_dotplot, {
      gene <- input$dotplot_Gene
      grouping_features <- input$dotplot_groups
      dotplot_data <- anthology_2020_v01 %>% tbl('grouped_stats') %>% 
        filter(Gene %in% gene) %>%
        group_by_at(vars(one_of(c('Gene', grouping_features)))) %>% 
        summarise(cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE),
                  cpm = mean(cpm)) %>% 
        as_tibble() %>% 
        tidyr::drop_na() %>% 
        full_join(., meta_filter %>% 
                    group_by_at(vars(one_of(grouping_features))) %>% 
                    summarise(Count = n())) %>% 
        mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
               Expression = cpm * (`%` / 100)) %>% 
        filter(!is.na(Count),
               `%` > 2) %>% 
        select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>% 
        filter(!is.na(Gene)) 
      if (length(grouping_features) == 2){
        colnames(dotplot_data)[c(2,3)] <- c("Group1","Group2")
        dotplot_data <- dotplot_data %>% 
          mutate(Column = paste(Group1,Group2),
                 Column = case_when(grepl('^\\d', Column) ~ paste0('X', Column),
                                    TRUE ~ Column))
      } else {
        colnames(dotplot_data)[c(2)] <- c("Group1")
        dotplot_data <- dotplot_data %>% 
          mutate(Column = Group1,
                 Column = case_when(grepl('^\\d', Column) ~ paste0('X', Column),
                                    TRUE ~ Column)) 
      }
      
      # cluster
      # make data square to calculate euclidean distance
      mat <- dotplot_data %>% 
        select(Gene, Column, Expression) %>%  # drop unused columns to faciliate widening
        tidyr::pivot_wider(names_from = Column, values_from = Expression) %>% 
        data.frame() # make df as tibbles -> matrix annoying
      row.names(mat) <- mat$Gene  # put gene in `row`
      mat <- mat[,-1] #drop gene column as now in rows
      mat[is.na(mat)] <- 0
      h_clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
      v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
      dotplot <- dotplot_data %>% 
        mutate(Gene = factor(Gene, levels = h_clust$labels[h_clust$order]),
               Column = factor(Column, levels = v_clust$labels[v_clust$order] %>% 
                                 gsub('\\.\\.', ' (', .) %>% 
                                 gsub('\\.$', ')', .) %>% 
                                 gsub('\\.', ' ', .))) %>% 
        ggplot(aes(x=Column, y = Gene, size = `%`, color = Expression)) + 
        geom_point() +
        cowplot::theme_cowplot() + 
        scale_color_viridis_c(option = 'magma') + 
        theme(axis.line  = element_blank()) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylab('') + xlab('') +
        theme(axis.ticks = element_blank()) 
      # labels
      order <- left_join(tibble::enframe(v_clust$labels[v_clust$order] %>% 
                                           gsub('\\.\\.', ' (', .) %>% 
                                           gsub('\\.$', ')', .) %>% 
                                           gsub('\\.', ' ', .), value = 'Column'),
                         dotplot_data)
      
      group1_labels <- order %>% 
        ggplot() + 
        geom_tile(aes(x = Column, y = 1, fill = Group1)) + 
        scale_fill_manual(name = grouping_features[1],
                          values = rep(c(pals::alphabet2() %>% unname(),
                                         pals::alphabet() %>% unname()), times = 10)) +
        theme_nothing() +
        aplot::xlim2(dotplot)
      # remove legend if same size as matrix as will be HUGE
      if ((dotplot_data$Group1 %>% unique() %>% length()) == ncol(mat)) { 
        group1_legend <- NULL
      } else {
        group1_legend <- plot_grid(get_legend(group1_labels + theme(legend.position="bottom")))
      }
      
      if (length(grouping_features) == 2){
        group2_labels <- order %>% 
          ggplot() + 
          geom_tile(aes(x = Column, y = 1, fill = Group2)) + 
          scale_fill_manual(name = grouping_features[2],
                            values = rep(c(pals::alphabet() %>% unname(),
                                           pals::alphabet2() %>% unname()), times = 10)) +
          theme_nothing() +
          aplot::xlim2(dotplot)
        if ((dotplot_data$Group2 %>% unique() %>% length()) == ncol(mat)) { 
          group2_legend <- NULL
        } else {
          group2_legend <- plot_grid(get_legend(group2_labels + theme(legend.position="bottom")))
        }
        
        group2_labels +
          group1_labels +
          dotplot + 
          group2_legend +
          group1_legend +
          plot_layout(ncol= 1, heights = c(0.1,0.1, 1, 0.1, 0.2))
      } else {
        group1_labels +
          dotplot + 
          group1_legend +
          plot_layout(ncol= 1, heights = c(0.1, 1, 0.1))
      }
    })
    output$dotplot <- renderPlot({
      make_dotplot()
    }, height = eventReactive(input$BUTTON_draw_dotplot, {input$dotplot_height %>% as.numeric()}))
  })
  
  ## diff testing table ----
  diff_table <- reactive({
    req(input$diff_table_select)
    if (input$diff_table_select == 'Cluster'){
      out <- 'diff_testing_A'
    } else if (input$diff_table_select == 'CellType') {
      out <- 'diff_testing_E'
    } else if (input$diff_table_select == 'SubCluster') {
      out <- 'diff_testing_G'
    } else {out <- 'diff_testing_C'}
    return(out)
  })
  
  output$make_diff_table <- DT::renderDataTable(server = TRUE, {
    gene <- input$diff_gene
    #cat(diff_table)
    if (input$search_by == 'Gene'){
      out <- anthology_2020_v01 %>% tbl(diff_table()) %>% 
        filter(gene_short_name %in% gene) %>% 
        arrange(-abs(normalized_effect)) 
    } else {
      filter_term <- input$diff_term
      out <- anthology_2020_v01 %>% tbl(diff_table()) %>% 
        filter(term %in% filter_term) %>% 
        arrange(-abs(normalized_effect)) 
    }
    out %>% 
      mutate(`AUC Score` = round(count / 55, 1),
             mean_auc = round(mean_auc, 2)) %>% 
      select(-status, -model_component, -count, -med_auc, -estimate, -test_val, -p_value) %>% 
      collect() %>% 
      DT::datatable(extensions = 'Buttons', rownames = F, 
                    filter = list(position = 'bottom', clear = FALSE),
                    options = list(pageLength = 10, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
    
  })
  
  output$formattable01 <- renderFormattable({formattable_01})
  output$formattable02 <- renderFormattable({formattable_02})
})
