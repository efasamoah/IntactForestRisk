

threshold.based.exposure <- function(data, scenario, id, th, baseline){
  
  #Calculate thresholds
  nme <- colnames(data)
  df_xp <- lapply(1:length(scenario),function(i){
    
    df <- cbind(ID=data[,id], data[,nme[grep(scenario[i],nme)]], baseline = data[,baseline])
    m_vars <- colnames(df)

    #Find the marginal baseline velocity
    df <- df %>% 
      group_by(ID) %>% 
      mutate(thresh = round(quantile(.data[[m_vars[grep(baseline,m_vars)]]], probs = th, na.rm = TRUE),1 )) %>% ungroup() %>%
      mutate(dplyr::across(all_of(3:7), ~ ifelse(.x>thresh,1,0)) )
    
    df_n <- cbind(ID=df[,"ID"], df[,m_vars[grep("rcp",m_vars)]])
    df_n <- aggregate(. ~ ID, data = df_n, sum) #sum all pixels for each forest type 
    
    return(df_n) }) %>% purrr::reduce(left_join, by = "ID")
  return(df_xp)
}
