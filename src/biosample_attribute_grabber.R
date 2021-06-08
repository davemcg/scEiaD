library(XML)
library(reutils)
library(tidyverse)

attribute_finder <- function(xml_list_obj){
  out <- list()
  for (i in 1:length(xml_list_obj)){
    out[[i]] <- xml_list_obj[[i]]$.attrs['attribute_name'] 
  }
  return(out)
}

value_grabber <- function(attribute, xml_list_obj){
  for (i in 1:length(xml_list_obj)){
    if (attribute %in% (xml_list_obj[[i]]$.attrs)){
      out <- xml_list_obj[[i]]$text
      return(out)
    }
  }
}


attribute_df_maker <- function(id){
  # fetch xml object from NCBI
  eutil_grab <- efetch(uid = id, db = 'biosample', retmode = 'xml') 
  # extract attributes, convert to list
  xml_list_obj <- eutil_grab[["//Attributes"]] %>% XML::xmlToList()
  biosample_title <-  (eutil_grab[["//Description"]]%>% XML::xmlToList())$Title
  # scan through list obj and find all attribures
  attribute_df <-  attribute_finder(xml_list_obj) %>% as.character() %>% data.frame()
  colnames(attribute_df)[1] <- 'attribute'
  # grab the attributes and stick into DF
  attribute_df <- attribute_df %>% rowwise() %>% mutate(value = value_grabber(attribute, xml_list_obj))
  attribute_df$id = id
  attribute_df <- bind_rows(attribute_df, c(attribute = 'biosample_title', value = biosample_title, 'id' = id))
  return(attribute_df)
}

attribute_l <- list()
for (i in unique(bind_rows(new_meta %>% filter(!biosample %in% 
                                               attribute_df$id), 
                           orig_meta %>% mutate(Age = as.character(Age)) %>%  
                           filter(!biosample %in% attribute_df$id)) %>% pull(biosample))){
  print(i)
  attribute_l[[i]] <- try({attribute_df_maker(i)})
  Sys.sleep(1)
}



# remove failed data pulls
failed <- c()
for (i in 1:length(attribute_l)){
  if (class(attribute_l[[i]]) == 'try-error'){
    failed <- c(failed, i)
}}

attribute_c <- attribute_l
attribute_c[failed] <- NULL
attribute_df <- attribute_c %>% bind_rows() %>% as_tibble()
save(attribute_df, attribute_l, file = 'data/2021_06_07_attribute_df.Rdata')

