
hm_maker <- function(pseudotime, 
                     num_of_genes = 20, 
                     onlyShowTF = FALSE, 
                     genes = NULL, 
                     output_smooth = FALSE, 
                     round_to = 0.5, 
                     column_title = NULL,
                     max_pseudotime = NULL){
  options(dplyr.summarise.inform=F) 
  if (is.null(genes)){
    genes <- diffPT[[pseudotime]] %>% 
      as_tibble(rownames = 'Gene') %>% 
      arrange(FDR,-abs(logFC)) %>% 
      head(num_of_genes) %>% 
      pull(Gene)
    if (onlyShowTF){
      genes <- diffPT[[pseudotime]] %>% 
        as_tibble(rownames = 'Gene') %>% 
        arrange(FDR,-abs(logFC)) %>% 
        filter(Gene %in% tf$ID) %>% 
        head(num_of_genes) %>% 
        pull(Gene)
    }
  }
  pseudoT <- pseudotime
  bc_data <- colData(sling$sling) %>% 
    as_tibble(rownames = 'Barcode') %>% 
    select(Barcode, one_of(pseudoT)) 
  names(bc_data)[2] <- 'PT'
  bc_time <- bc_data %>% 
    filter(!is.na(PT)) %>% 
    arrange(PT) %>% select(Barcode, PT)
  bc <- bc_time %>% pull(Barcode)
  pt <- bc_time %>% pull(PT)
  ct <- umap %>% right_join(bc %>% enframe(value = 'Barcode'), by = 'Barcode') %>% pull(CellType_predict)
  
  mat <- logcounts(sling$sling)[genes, bc]
  
  long <- mat %>% t() %>% as_tibble(rownames = 'Barcode')  %>% 
    left_join(umap %>% right_join(bc_time, by = 'Barcode'), by = 'Barcode') %>% 
    select(PT, contains('ENSG'), CellType_predict, organism)
  
  if (round_to %% 1 == 0 ){
    multiplier = 1
    digits = round_to
  } else {
    multiplier = 1 / round_to
    digits = 0
  }
  
  long$group <- (round(long$PT * multiplier, digits = digits) / multiplier) 
 #mutate(group = (round(PT * multiplier, digits = round_to)) / multiplier)
  if (!is.null(max_pseudotime)){
    long <- long %>% filter(group <= max_pseudotime)
  }
  long <- long %>% 
    mutate(group = as.character(group)) %>% 
    right_join(seq(min(long$group),max(long$group), by = round_to) %>% 
                 enframe(value = 'group') %>% 
                 mutate(group = as.character(group)), by = 'group') %>% 
    mutate(group = as.numeric(group)) %>% 
    arrange(group)
                               
  new_colname <- long  %>% 
    select(CellType_predict, group) %>% 
    group_by(group, CellType_predict) %>% 
    summarise(count = n()) %>% 
    top_n(1, count) %>% 
    mutate(id = paste0(group, '__', CellType_predict))
  org <- long  %>% 
    select(organism, group) %>% 
    group_by(group, organism) %>% 
    summarise(count = n()) %>% 
    top_n(1, count)
  PT <- long  %>% 
    select(PT, group) %>% 
    group_by(group) %>% 
    summarise(PT = mean(PT)) %>% 
    mutate(PT = round(PT))
  new_colname <- left_join(new_colname, org, by = 'group') %>% 
    mutate(id = paste0(id, '__', organism)) %>% 
    left_join(PT, by = 'group') %>% 
    mutate(id = paste0(id, '__', PT)) 
  
  mat_smooth <- long %>% 
    select(contains("ENSG"), group) %>% 
    group_by(group) %>% 
    summarise_all(mean) %>% 
    left_join(new_colname, by = 'group') %>% 
    select(contains('ENS'), id) 
  mat_smooth_out <- mat_smooth
  id <- mat_smooth$id 
  mat_smooth <- mat_smooth %>% select(-id) %>% t()
  id <- gsub('NA','Unlabelled',id)
  
  
  # set colors
  ct_p_colors <- c(pals::alphabet2(), pals::alphabet())
  ct_p_colors <- ct_p_colors[1:(umap$CellType_predict %>% unique() %>% sort() %>% length())]
  names(ct_p_colors) <- umap$CellType_predict %>% unique() %>% sort()
  ct_p_colors <- c('white', ct_p_colors)
  names(ct_p_colors)[1] <- 'Unlabelled'
  
  anno_mark_at <- seq(1,ncol(mat_smooth), round(ncol(mat_smooth)/long$group %>% max() %>% floor()))
  anno_mark_labels <- (long$group %>% min() %>% floor()):(long$group %>% max() %>% floor())
  if (length(anno_mark_at) != length(anno_mark_labels)){
    anno_mark_labels <- anno_mark_labels[-length(anno_mark_labels)]
  }
  ha_column = HeatmapAnnotation(
    show_annotation_name = FALSE,
    Pseudotime = anno_mark(at = anno_mark_at, labels = anno_mark_labels, labels_rot = 45), 
    df = data.frame(
      #Pseudotime = str_split(id, '__') %>% map(4) %>% unlist() %>% as.numeric(),
      #Pseudotime = anno_mark(at = c(1:10), labels = 1:10),
      CellType = str_split(id, '__') %>% map(2) %>% unlist(),
      Organism = str_split(id, '__') %>% map(3) %>% unlist()),
    
    col = list(CellType = ct_p_colors,
               Organism = c('Homo sapiens' = 'orange',
                            'Mus musculus' = 'green',
                            'Macaca fascicularis' = 'blue',
                            'Unlabelled' = 'white'))
  ) 
  # label as TF or not
  gene_type <- row.names(mat_smooth) %in% tf$ID
  gene_type <- gsub(FALSE, 'No', gene_type)
  gene_type <- gsub(TRUE, 'Yes', gene_type)
  gene_type_col = list('Yes' = 'black', 'No' = 'white')
  # replace ENS with ID
  row.names(mat_smooth) <- 
    row.names(mat_smooth) %>% 
    enframe(value = 'hs_gene_id') %>% 
    left_join(gene_id_converter %>% 
                select(hs_gene_id, hs_gene_name) %>% unique(),
              by = 'hs_gene_id') %>% 
    pull(hs_gene_name)
  # make heatmap
  hm <- Heatmap(mat_smooth,
                cluster_columns = FALSE, 
                cluster_rows = TRUE, 
                na_col = 'white',
                show_column_names = FALSE,
                col=viridis::viridis(100),
                top_annotation = ha_column,
                name = 'log Counts',
                column_title = column_title)
  gene_col <- Heatmap(gene_type, name = "TF", col = gene_type_col)
  if (!output_smooth){
    if (onlyShowTF){
      hm
    } else {
      hm + gene_col
    }
  } else {
    mat_smooth
  }
}

# tf_plots <- list()
# for (i in colData( sling$sling) %>% colnames() %>% grep('slingP',.,value = TRUE)){
#   tf_plots[[i]] <- try({hm_maker(i, 10, onlyShowTF = TRUE)})
# }
# 
# gene_plots <- list()
# for (i in colData( sling$sling) %>% colnames() %>% grep('slingP',.,value = TRUE)){
#   gene_plots[[i]] <- try({hm_maker(i, 40, onlyShowTF = FALSE)})
# }


#hm_maker('slingPseudotime_1', 10, onlyShowTF = TRUE, genes = c('ENSG00000134438','ENSG00000114315', 'ENSG00000148400'))


lm_maker <- function(tidyPT){
  # remove cells in low n groups
  # grouping by GENE, GROUP (pseudotime rounded), and STUDY
  # also remove cells that are ALL ZERO in a grouping
  tidyPT_c <- tidyPT %>% 
    left_join(tidyPT %>% 
                group_by(Gene, group, study_accession) %>% 
                summarise(Count = n(), Zero = sum(CPM)), 
              by = c('Gene','group', 'study_accession')) %>% 
    filter(Count > 20, Zero != 0) 
  
  # split tidyPT_c into two objectds
  # one to run lm with the study_accession covariate
  # and another to run without
  #   (because there's only one study_accession)
  tidyPT_c_low_study_n <- tidyPT_c %>% left_join(tidyPT_c %>% 
                                                   group_by(Gene, group, study_accession) %>% 
                                                   summarise(Count = n()) %>% summarise(Study_Count = n()),
                                                 by = c('Gene', 'group')) %>% 
    filter(Study_Count == 1)
  tidyPT_c_higher_study_n <- tidyPT_c %>% left_join(tidyPT_c %>% 
                                                      group_by(Gene, group, study_accession) %>% 
                                                      summarise(Count = n()) %>% summarise(Study_Count = n()),
                                                    by = c('Gene', 'group')) %>% 
    filter(Study_Count > 1)
  
  # run lm (by gene - group) with study_accession as covariate (where possible)
  print("Running LM with study_accession as covariate")
  tidyPT_lm1 <- tidyPT_c_higher_study_n %>%
    group_by(Gene, group) %>%   
    nest() %>% 
    mutate(lm = map(data, function(df) lm(CPM ~ PT + study_accession, data = df)),
           t_results = map(lm, tidy),
           g_results = map(lm, glance),
           a_results = map(lm, augment, interval = 'confidence'))
  
  tidyPT_lm2 <- tidyPT_c_low_study_n %>%
    group_by(Gene, group) %>%   
    nest() %>% 
    mutate(lm = map(data, function(df) lm(CPM ~ PT, data = df)),
           t_results = map(lm, tidy),
           g_results = map(lm, glance),
           a_results = map(lm, augment, interval = 'confidence'))
  
  print("Running LM for each gene - group")
  tidyPT_lmY <- tidyPT_c %>%
    group_by(Gene, group) %>%   
    nest() %>% 
    mutate(lm = map(data, function(df) lm(CPM ~ PT, data = df)),
           t_results = map(lm, tidy),
           g_results = map(lm, glance),
           a_results = map(lm, augment, interval = 'confidence'))
  
  tidyPT_lm_covariate <- bind_rows(tidyPT_lm1, tidyPT_lm2)
  
  # run separate lm for each gene - group - study_accession set
  print("Running LM for each gene - group - study_accession")
  tidyPT_lmX <- tidyPT_c %>%
    group_by(Gene, group, study_accession) %>%
    nest() %>%
    mutate(lm = map(data, function(df) lm(CPM ~ PT, data = df)),
           t_results = map(lm, tidy),
           g_results = map(lm, glance),
           a_results = map(lm, augment, interval = 'confidence'))
  
  
  out <- list()
  out$lm_covariate <- tidyPT_lm_covariate
  out$lm <- tidyPT_lmY
  out$lm_individual  <- tidyPT_lmX
  out$tidyPT <- tidyPT
  out
  
}

make_tidy_time <- function(pseudoT, genes = NULL, rounding_digits = 0){
  #pseudoT <- 'slingPseudotime_18'
  if (is.null(genes)){
    # if no genes are given, then chose the FDR 0 transcription factors
    pt_tf <- diffPT_tibble %>% 
      filter(Trajectory == pseudoT) %>%
      left_join(gene_id_converter %>% 
                  select(Gene = hs_gene_id, hs_gene_name)) %>% 
      filter(FDR == 0, Gene %in% tf$ID) %>% 
      unique() %>% 
      arrange(-abs(logFC)) %>% pull(Gene)
  } else {
    pt_tf <- genes
  }
  
  bc_data <- colData(sling$sling) %>% 
    as_tibble(rownames = 'Barcode') %>% 
    select(Barcode, one_of(pseudoT)) 
  names(bc_data)[2] <- 'PT'
  bc_time <- bc_data %>% 
    filter(!is.na(PT)) %>% 
    arrange(PT) %>% select(Barcode, PT)
  
  mat <- logcounts(sling$sling)[pt_tf, bc_time$Barcode]
  
  if (rounding_digits %% 1 == 0 ){
    multiplier = 1
  } else {
    multiplier = 1 / rounding_digits
    rounding_digits = 0
  }
  
  mat %>% t() %>% as_tibble(rownames = 'Barcode') %>% 
    left_join(umap %>% right_join(bc_time)) %>% 
    select(PT, contains('ENSG'), Barcode, CellType_predict, organism, study_accession, batch, Platform) %>% 
    pivot_longer(-c(Barcode, PT, CellType_predict, organism, study_accession, batch, Platform), names_to = 'Gene', values_to = 'CPM') %>% 
    mutate(group = (round(PT * multiplier, digits = rounding_digits)) / multiplier)
} 
