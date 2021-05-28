library(XML)
library(reutils)
library(tidyverse)

attribute_finder <- function(xml_list_obj){
  out <- list()
  for (i in 1:length(xml_list_obj)){
    out[[i]] <- xml_list_obj[[i]]$.attrs['display_name'] 
  }
  return(out)
}

value_grabber <- function(attribute, xml_list_obj){
  for (i in 1:length(xml_list_obj)){
    if (xml_list_obj[[i]]$.attrs['display_name'] == attribute){
      out <- xml_list_obj[[i]]$text
      return(out)
    }
  }
}

attribute_df_maker <- function(id){
  eutil_grab <- efetch(uid = id, db = 'biosample', retmode = 'xml') 
  xml_list_obj <- eutil_grab[["//Attributes"]] %>% XML::xmlToList()
  
  attribute_df <-  attribute_finder(xml_list_obj) %>% as.character() %>% data.frame()
  colnames(attribute_df)[1] <- 'attribute'
  attribute_df <- attribute_df %>% rowwise() %>% mutate(value = value_grabber(attribute, xml_list_obj))
  attribute_df$id = id
  return(attribute_df)
}

attribute_l <- list()
for (i in unique(new_meta$biosample)){
  print(i)
  attribute_l[[i]] <- try({attribute_df_maker(i)})
  Sys.sleep(2)
  
}
