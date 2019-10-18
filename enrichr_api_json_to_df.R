enrichr_API_output_to_df <- function(enrichr_API_output_jsonlite){
  enrichr_API_output_jsonlite_list <- jsonlite::fromJSON(enrichr_API_output_jsonlite)
  
  database_str <- names(enrichr_API_output_jsonlite_list)
  
  json_file2 <- lapply(enrichr_API_output_jsonlite_list[[1]], function(x) {
    temp_2 <- list()
    
    for (n in seq_along(x)) {
      if (length(x[[n]]) == 1) {
        temp_2[[n]] <- x[[n]]
      }
      else if (length(x[[n]]) > 1) {
        temp_2[[n]] <- paste0(x[[n]], collapse = ', ')
      }
      else {
        stop('error', call. = T)
      }
    }
    
    return(temp_2)
    
  })
  
  unlist_json <-
    lapply(
      X = json_file2,
      FUN = function(x) {
        unlist(x)
      }
    )
  
  unlist_json2 <- as.data.frame(rlist::list.rbind(unlist_json))
  
  colnames(unlist_json2) <-
    c('Database', 'Term', 'P-value', 'Odds_Ratio', 'Combined_Score', 'Genes', 'Adj_P-value', 'no_1', 'no_2')
  
  print(colnames(unlist_json2))
  print(database_str)  

  unlist_json2$no_1 <- NULL
  unlist_json2$no_2 <- NULL
  unlist_json2$Database <- database_str
  
  return(unlist_json2)
}
