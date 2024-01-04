

#write.csv(file = "All Species Exclusive Peptides.csv", good_pep_excl_tot)
#write.csv(file = "88 Species Exclusive Peptides for RandomForest Model.csv", good_pep_t)
#write.csv(file = "2247 Informative Peptides for Correlation Model.csv", feature_imp_cor_data)

#save(model5, good_pep_excl_tot, good_pep_t, feature_imp_cor_data, file =  "SIPMS_ModelData.RData")
load(SIPMS_ModelData)
data(SIPMS_ModelData)
`%notin%` <- Negate(`%in%`)
# Function to simulate the Prediction process
Species_Correlation <- function(Peptide_list) {
  # Simulating some prediction logic
  Sys.sleep(2)  # Simulating processing time
  load(file = "SIPMS_ModelData.RData",.GlobalEnv)
  ## No. samples
  No_Sample <- ncol(Peptide_list %>% select(-Peptide))
  
  ##Replacing Zeros with NA
  Peptide_list[Peptide_list==0] <- NA
  
  nn_nam <- names(Peptide_list)
  ## Normalization based on total sum
  Peptide_list <- Sum_Normalization(x = Peptide_list %>% select(-Peptide), dff_ = Peptide_list %>% select(Peptide))
  ## Find Imputation value for each sample
  
  imp <- apply(Peptide_list %>% select(-Peptide), 2, function(x){min(x, na.rm = T)})
  ## Change names
  names(Peptide_list)[-1] <- nn_nam[-1]
  
  
  Species_check <- Peptide_list %>% filter(Peptide %in% feature_imp_cor_data$Peptide) %>% distinct()
  Species_check <- left_join(feature_imp_cor_data, Species_check, by = "Peptide")
  cor_m <- cor(x = Species_check %>% select(-Protein.Accession, -Peptide), use = "complete.obs")
  cor_m_res <- as.data.frame(cor_m) %>% select(names(Peptide_list)[-1]) %>% 
    mutate(SampleID = rownames(.)) %>% filter(SampleID %notin% names(Peptide_list)[-1]) %>% 
    rowwise() %>% 
    mutate(Species = strsplit(SampleID, "_")[[1]][1]) %>% ungroup() %>% 
    gather(key = "UnknownSample", value = "Correlation", 1:(ncol(.)-2)) %>% 
    group_by(UnknownSample, Species) %>% summarise(Mean_Cor = mean(Correlation, na.rm = T), nhit = n())
  
  cor_grp <- as.data.frame(cor_m) %>% select(names(Peptide_list)[-1]) %>% 
    mutate(SampleID = rownames(.)) %>% filter(SampleID %notin% names(Peptide_list)[-1]) %>% 
    rowwise() %>% 
    mutate(Species = strsplit(SampleID, "_")[[1]][1]) %>% ungroup() %>% 
    gather(key = "UnknownSample", value = "Correlation", 1:(ncol(.)-2)) %>% 
    group_by(UnknownSample, Species) %>% summarise(Mean_Cor = mean(Correlation, na.rm = T), nhit = n()) %>% 
    ungroup() %>% group_by(UnknownSample) %>% mutate(Rank_Cor = rank(-Mean_Cor))
  
  
  cor_samp <- as.data.frame(cor_m) %>% select(names(Peptide_list)[-1]) %>% 
    mutate(SampleID = rownames(.)) %>% filter(SampleID %notin% names(Peptide_list)[-1]) %>% 
    rowwise() %>% 
    mutate(Species = strsplit(SampleID, "_")[[1]][1]) %>% ungroup() %>% 
    gather(key = "UnknownSample", value = "Correlation", 1:(ncol(.)-2))
  
  cor_p <- left_join(cor_samp, cor_grp, by = c("UnknownSample", "Species")) %>% 
    group_by(UnknownSample) %>% do({
      dd = .
      t = 1
      res <- data.frame()
      while(t < 8){
        dd1 <- dd %>% filter(Rank_Cor == t)
        dd2 <- dd %>% filter(Rank_Cor > t)
        
        rr <- data.frame(Pvalue = t.test(dd1$Correlation, dd2$Correlation)$p.value, 
                         Species = dd1$Species[1], Rank = t)
        res <- bind_rows(res, rr)
        t = t + 1
      }
      res 
    })
  
  
  
  fin_cor <- left_join(cor_m_res, cor_p, by = c("UnknownSample", "Species")) %>% 
    rowwise() %>% mutate(Rank = ifelse(is.na(Rank), 8, Rank), 
                         Pvalue = ifelse(is.na(Pvalue), 1, Pvalue)) %>% 
    ungroup()

  return(fin_cor)
}
Check_Presence <- function(res_cor){
  
  check_result <- res_cor %>% ungroup() %>% group_by(UnknownSample) %>% do({
    dd = .
    sp_check <- dd %>% filter(Mean_Cor>=0.7)
    if (nrow(sp_check>=1)){
      ees <- data.frame(Present = "Yes")
    }
    else{
      ees <- data.frame(Present = "No")
    }
    
  })
  return(check_result)
}
Imputed_value <- function(Peptide_list){
  Peptide_list[Peptide_list==0] <- NA
  
  nn_nam <- names(Peptide_list)
  ## Normalization based on total sum
  Peptide_list <- Sum_Normalization(x = Peptide_list %>% select(-Peptide), dff_ = Peptide_list %>% select(Peptide))
  ## Find Imputation value for each sample
  names(Peptide_list)[-1] <- nn_nam[-1]
  imp <- apply(Peptide_list %>% select(-Peptide), 2, function(x){min(x, na.rm = T)})
  return(imp)
}
Imp_replaced <- function(imp_val, Peptide_list){
  pp_ <- Peptide_list
  for (i in names(Peptide_list[,-1])){
    pp_[is.na(pp_[,i]), i] <- runif(n = length(which(is.na(pp_[,i]))==T),
                                                          min = imp_val[[i]]/500000, 
                                                          max = imp_val[[i]]/5000)
  }
  return(pp_)
}
Species_Prediction <- function(Peptide_list){
  
    
  
  
  
}
Sum_Normalization <- function(x, dff_){
  
  
  for (col in 1:ncol(x)){
    
    dff <- x[,col]/sum(x[,col], na.rm = T)
    dff_ <- bind_cols(dff_, dff)
    
  }
  
  return(dff_)
}

# Install required packages if not already installed
if (!require("shiny")) install.packages("shiny")
if (!require("shinydashboard")) install.packages("shinydashboard")

# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)  # for interactive tables
library(ggplot2)
library(randomForest)
require(dplyr)
require(tidyverse)
require(tidyr)
# Define the UI

ui <- fluidPage(
  data(SIPMS_ModelData, package = "SIPMS"),
  img(src = "https://github.com/hassanakthv/SIPMS/assets/43888767/70437bd0-88f8-4591-8b08-c4f5215e6713",
      alt = "SSE", height = 60, width = 120),
  titlePanel("Species Search Engine - SSE"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Peptide List (CSV)",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      checkboxInput("check_database", "Check against built-in database"),
      actionButton("analyze_btn", "Search")
    ),
    mainPanel(HTML('<img src="https://en.wikipedia.org/wiki/File:Thermo_Fisher_Scientific_logo.svg" width="100%" height="auto">'),
      tabsetPanel(
        tabPanel("Description",verbatimTextOutput("description"), 
                 HTML('<img src="https://github.com/hassanakthv/SIPMS/assets/43888767/b38933a0-56c3-4b79-b6b5-5944f864477b" width="100%" height="auto">')),
        tabPanel("Species Correlation", DTOutput("species_corr"), verbatimTextOutput("presence_species")),
        tabPanel("Search Results - Samples Found in DB", DTOutput("analysis_result"), 
                 plotOutput("prediction_plot")),
        tabPanel("Search Results - Samples Not found in DB", DTOutput("similarity_result"), 
                 plotOutput("similarity_plot"))
      )
    )
  )
)

# Define the server
server <- function(input, output, session) {
  data(SIPMS_ModelData, package = "SIPMS")
  options(shiny.maxRequestSize=30*1024^2)
  # Function to read the uploaded CSV file
  peptides_data <- data.frame()
  peptides_data <- reactive({
    req(input$file)
    if(!is.null(input$file)){
    read.csv(input$file$datapath)[,-1]
    }
    else{
      "Please upload a dataset"
    }
  })
  
  # Function to check the similarity with the built-in database
  sample_correlation <- reactive({
    if(!is.null(input$file)){
    if (input$check_database) {
      # Perform similarity check logic here with the built-in database
      # For demonstration purposes, I'll assume a random number of similar peptides
      Species_Correlation(peptides_data()) 
      
    } 
    } else {
      print("You should see something here")
    }
  })
  
  # Output similar peptides as a table
  output$species_corr <- renderDT({
    if(!is.null(input$file)){
    sample_correlation() %>% select(-Pvalue, -Rank) %>%  datatable(extensions = 'Buttons', options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))) %>%
      formatSignif('Mean_Cor', 2)
  }
  else{
    " "
  }
  })
  
  
  pres_check <- reactive({
    
    r_pres <- Check_Presence(sample_correlation())
    
  })
  output$presence_species <- renderText({
    sampless <- ""
    sample_str_p <- ""
    U_Pres <- unique(pres_check() %>% filter(Present == "Yes"))$UnknownSample
    if (length(U_Pres)!=0){
    for (i in U_Pres){
      sample_str_p <- paste0(sample_str_p, " ", i)
      }
    }
    U_Abs <- unique(pres_check() %>% filter(Present == "No"))$UnknownSample
    sample_str_a <- ""
    if (length(U_Abs)!=0){
      for (j in U_Abs){
        sample_str_a <- paste0(sample_str_a, " ", j)
      }
    }
    if (sample_str_a!="" & sample_str_p!=""){
    sampless <- paste0(sample_str_p,": These samples are present in our current species database.", 
           '\n',sample_str_a, ": These samples are NOT present in our current species database.")
    }
    if (sample_str_a=="" & sample_str_p!=""){
    sampless <-  paste0(sample_str_p,": These samples are present in our current species database.")
    }
    if (sample_str_a!="" & sample_str_p==""){
    sampless <- paste0(sample_str_a, ": These samples are NOT present in our current species database.")
    }
    print(sampless)
    })
    
  
  # Function to perform correlation analysis and prediction
  analysis_result <- reactive({
    if (input$analyze_btn){
    ## Prediction Score for samples in the database
      sp <- unique(pres_check() %>% filter(Present == "Yes"))$UnknownSample 
      cleaned_pep <- peptides_data() %>% select(Peptide, sp)
      cleaned_pep[is.na(cleaned_pep)] <- NA
      nn_nam <- names(cleaned_pep)
      cleaned_pep <- Sum_Normalization(x = cleaned_pep[,-1], dff_ = cleaned_pep[,1])
      names(cleaned_pep) <- nn_nam
      rf_pep <- cleaned_pep %>% filter(Peptide %in% good_pep_t$Peptide) %>% distinct()
      imp_val <- Imputed_value(peptides_data() %>% select(Peptide, sp))
      rf_pep <- Imp_replaced(imp_val = imp_val, Peptide_list = rf_pep)
      rf_pep[,-1] <- log10(rf_pep[,-1])
      
      res_pred <- as.data.frame(predict(model5, t(rf_pep[,-1]), type = "prob"))
      for (i in rownames(res_pred)){
        res_pred[i,] <- res_pred[i,]*(sample_correlation() %>% 
                                        filter(UnknownSample == i))$Mean_Cor
        
        
      }
      res_pred <- round(res_pred/rowSums(res_pred),3) 
      res_pred <- res_pred %>% mutate(UnknownSample = rownames(.)) %>% 
        gather(key = "Species", 
               value = "Probability_Score", 1:8) %>% 
        group_by(UnknownSample) %>% mutate(Rank = rank(-Probability_Score)) %>% ungroup() %>% 
      left_join(.,sample_correlation() %>% select(UnknownSample, Species, Pvalue), 
                by = c("UnknownSample", "Species"))
    }
  })
  
  similarity_result <- reactive({
    if (input$analyze_btn){
      ## Prediction Score for samples in the database
      sp1 <- unique(pres_check() %>% filter(Present == "No"))$UnknownSample 
      sim_res <- sample_correlation() %>% filter(UnknownSample %in% sp1)
      sim_res <- sim_res %>% group_by(UnknownSample) %>% 
        mutate(Similarity_Score = round(Mean_Cor/sum(Mean_Cor, na.rm = T),3)) %>% 
        mutate(Rank = rank(-Similarity_Score)) %>% 
        select(UnknownSample, Species, Similarity_Score, Rank, Pvalue)
      
    }
      
  })
     
  # Output analysis result
  output$analysis_result <- renderDT({
    analysis_result() %>% arrange(Rank, UnknownSample) %>% datatable(extensions = 'Buttons', options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))) #%>%
      #formatSignif('Mean_Cor', 2)
  })
  
  
  output$prediction_plot <- renderPlot({
    gg_data <- analysis_result() 
    gg <- ggplot(data = gg_data, aes(x = reorder(UnknownSample, -Rank), group= Rank, y = 100*Probability_Score, fill = Species, label = Species)) + 
      geom_bar(stat = "identity", position = position_dodge(width = 1)) + theme_bw() + 
      geom_text(position = position_dodge(width = 1,preserve = "total"), size = 4) + 
      labs(x = "Samples", y = "Probablity Score (%)") + coord_flip()
    gg
    
  })
  
  output$similarity_result <- renderDT({
    similarity_result() %>% arrange(Rank, UnknownSample) %>% datatable(extensions = 'Buttons', options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))) #%>%
    #formatSignif('Mean_Cor', 2)
  })
  
  
  output$similarity_plot <- renderPlot({
    gg_data1 <- similarity_result()
    gg1 <- ggplot(data = gg_data1, aes(x = reorder(UnknownSample, -Rank), group= Rank, y = 100*Similarity_Score, fill = Species, label = Species)) + 
      geom_bar(stat = "identity", position = position_dodge(width = 1)) + theme_bw() + 
      geom_text(position = position_dodge(width = 1,preserve = "total"), size = 4) + 
      labs(x = "Samples", y = "Similarity Score (%)") + coord_flip()
    gg1
    
  })
  
  output$description <- renderText({
  paste0("Current Species Search Engine comprise 8 different species: Grey Seal, Harbor Seal
         Southern Elephant Seal, Grey Whale, Cuvier's Beaked Whale,
         Steller's Sea cow, Whooper Swan and Peregrine Falcon")
    })

}
# Run the Shiny app
shinyApp(ui, server)
