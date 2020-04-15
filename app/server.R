# server.R
time <- Sys.time()
cat(file = stderr(), 'Server Go!\n')
#options(shiny.trace=TRUE)
options(shiny.sanitize.errors = FALSE)

library(ggplot2)
library(scattermore)
library(dplyr)
require(pals)
library(pool)
library(RSQLite)
library(cowplot)
library(ggrepel)
library(patchwork)

anthology_2020_v01 <- dbPool(drv = SQLite(), dbname = "./www/anthology_limmaFALSE_nf5000-d50-k7.sqlite", idleTimeout = 3600000)

# filter
meta_filter <- anthology_2020_v01 %>% 
  tbl('metadata') %>% 
  as_tibble() %>% 
  filter(!is.na(CellType_predict), 
         !is.na(study_accession), 
         !CellType_predict %in% c('Astrocytes', 'Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
         !grepl('RPE|Vascul', CellType_predict))  %>% 
  dplyr::select(-UMI) 

# get coords for cell labels
celltype_predict_labels <- meta_filter %>% group_by(CellType_predict) %>% summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
celltype_labels <- meta_filter %>% group_by(CellType) %>% summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
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
    # server gene / tx to UI side
    if (is.null(query[['Gene']])){
      updateSelectizeInput(session, 'Gene',
                           choices = anthology_2020_v01 %>% tbl('genes') %>% as_tibble() %>% pull(1),
                           options = list(placeholder = 'Type to search'), 
                           selected = 'CRX',
                           server = TRUE)
    }
    # server gene / tx to UI side
    if (is.null(query[['meta_column']])){
      updateSelectizeInput(session, 'meta_column',
                           choices = meta_filter %>% 
                             dplyr::select(nCount_RNA:cluster) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    
    # dotplot genes 
    if (is.null(query[['dotplot_Gene']])){
      updateSelectizeInput(session, 'dotplot_Gene',
                           choices = anthology_2020_v01 %>% tbl('genes') %>% as_tibble() %>% pull(1),
                           options = list(placeholder = 'Type to search'), 
                           selected = c('RHO','WIF1','CABP5', 'AIF1','AQPT4','ARR3','ONECUT1','GRIK1','GAD1','POU4F2'),
                           server = TRUE)
    }
    
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
    
    # gene scatter plot ------------
    gene_scatter_plot <- eventReactive(input$BUTTON_draw_scatter, {
      cat(file=stderr(), 'gene scatter plot\n')
      gene <- input$Gene
      pt_size <- as.numeric(input$pt_size) / 10
      query = paste0('select * from cpm where Gene in ("',
                     paste(gene, collapse='","'),'")')
      p <- dbGetQuery(anthology_2020_v01, query) %>% 
        left_join(.,meta_filter, by = 'Barcode') %>% 
        filter(cpm > 1) %>% 
        as_tibble()
      
      plot <- p %>% ggplot() + 
        geom_scattermore(data = meta_filter,
                         aes(x = UMAP_1, y = UMAP_2), 
                         pointsize = 0.3, color = 'gray', alpha = 0.1) +
        geom_scattermore(aes(x = UMAP_1, y = UMAP_2, colour = cpm), 
                         pointsize = pt_size, 
                         alpha = 0.3) +
        #guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
        scale_color_viridis_c(option = 'magma') +
        theme_cowplot() + 
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank()) +
        annotate("text", -Inf, Inf, label = paste0(gene, '\nexpression'), hjust = 0, vjust = 1, size = 6)
      
      plot
      
    })
    output$gene_scatter_plot <- renderPlot({
      gene_scatter_plot()
    })  
    
    
    # metadata plot --------------
    meta_plot <- eventReactive(input$BUTTON_draw_meta, {
      meta_column <- input$meta_column
      transform <- input$meta_column_transform
      if (transform == 'log2' && is.numeric(meta_filter[,meta_column] %>% pull(1))){
        cat('log2 time')
        meta_filter[,meta_column] <- log2(meta_filter[,meta_column] + 1)
      }
      plot <- meta_filter %>% 
        filter(!is.na(!!as.symbol(meta_column))) %>% 
        #sample_frac(0.3) %>% 
        ggplot() + 
        geom_scattermore(aes(x = UMAP_1, 
                             y = UMAP_2),
                         color = 'gray',
                         pointsize = 0.5, 
                         alpha = 0.2) +
        geom_scattermore(aes(x = UMAP_1, 
                             y = UMAP_2, 
                             colour = !!as.symbol(meta_column) ), 
                         pointsize = 0.8, 
                         alpha = 0.2) +
        guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) + 
        theme_cowplot() + 
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank()) +
        annotate("text", -Inf, Inf, label = "Metadata", hjust = 0, vjust = 1, size = 6)
      
      if (is.numeric(meta_filter[,meta_column] %>% pull(1)) ){
        color <- scale_color_viridis_c()
      } else {
        color <- scale_colour_manual(values = rep(c(pals::alphabet() %>% unname(),
                                                    pals::alphabet2() %>% unname()), 
                                                  times = 10))
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
      if (meta_column == 'cluster'){
        plot + more + theme(legend.position = 'none') + color
      } else {
        plot + more + color
      }
      
    })
    output$meta_plot <- renderPlot({
      meta_plot()
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
    
    # NEW SECTION -----------
    ## dotplot ---------
    make_dotplot <- eventReactive(input$BUTTON_draw_dotplot, {
      gene <- input$dotplot_Gene
      grouping_features <- c('CellType_predict', 'organism')
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
        mutate(Column = paste(`CellType_predict`, organism)) %>% 
        filter(!is.na(Gene)) 
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
      org_labels <- order %>% 
        ggplot() + 
        geom_tile(aes(x = Column, y = 1, fill = organism)) + 
        scale_fill_manual(values = rep(c(pals::alphabet() %>% unname(),
                                         pals::alphabet2() %>% unname()))) +
        theme_nothing() +
        ggtree::xlim2(dotplot)
      ct_labels <- order %>% 
        ggplot() + 
        geom_tile(aes(x = Column, y = 1, fill = CellType_predict)) + 
        scale_fill_manual(values = rep(c(pals::alphabet2() %>% unname(),
                                         pals::alphabet() %>% unname()))) +
        theme_nothing() +
        ggtree::xlim2(dotplot)
      
      org_legend <- plot_grid(get_legend(org_labels + theme(legend.position="bottom")))
      ct_legend <- plot_grid(get_legend(ct_labels + theme(legend.position="bottom")))
      org_labels +
        ct_labels +
        dotplot + 
        org_legend +
        ct_legend +
        plot_layout(ncol= 1, heights = c(0.1,0.1, 1, 0.1, 0.2))
      
      
    })
    output$dotplot <- renderPlot({
      make_dotplot()
    })  
    
    
  })
})
