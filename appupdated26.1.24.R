
############# DEPENDENCIES #####################

packages <- c("BiocManager", "shinydashboard", "devtools", "tidyverse", "shiny", "ggrepel")

for (package in packages){
  if (!(package %in% installed.packages()[,"Package"])) {
    # If not installed, install the package
    install.packages(package)
  }
  }

if (!("BiocManager" %in% installed.packages()[,"Package"])) {
  # If not installed, install the package
  install.packages("BiocManager")
}

if (!("shinydashboard" %in% installed.packages()[,"Package"])) {
  # If not installed, install the package
  install.packages("shinydashboard")
}

if (!("devtools" %in% installed.packages()[,"Package"])) {
  # If not installed, install the package
  install.packages("devtools")
}

if (!("edgeR" %in% installed.packages()[,"Package"])) {
  # If not installed, install the package from Bioconductor
  BiocManager::install("edgeR")
}

if (!("annotables" %in% installed.packages()[,"Package"])){
  devtools::install_github("stephenturner/annotables")
}

library(shinydashboard)
library(shiny)
library(tidyverse)
library(edgeR)
library(annotables)
library(ggrepel)


######################################



pcnames <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

########## U I #######

ui <- dashboardPage(
  dashboardHeader(title = "RNA-seq analyser"),
  
  
  dashboardSidebar(
            sidebarMenu( 
                        menuItem("Data upload", 
                                 tabName = "dashboard", 
                                 icon = icon("dashboard")),
                       
                        menuItem("DGE analysis", 
                                 tabName = "dge_analysis", 
                                 icon = icon("th"))
                        )
                  ),
  
  
  dashboardBody(  
    tabItems(
      tabItem(tabName = "dashboard",
              h1("Data upload"),
          fluidRow(
    
                  box(
                    title = "Inputs", status = "warning",
                    "Please upload your file with count matrix. The first column should contain Ensmbl gene IDs. The column names should be in following format: group1, group2 etc, where group corresponds to group name and numer corresponds to individual sample. For example control1, control2, lps1, lps2 etc.", 
                    fileInput("file123", "Choose a file"),
                    br(), 
                    "Please provide names of genes that you want to view.",
                    textInput("input_names", "Enter names (comma-separated):")
                      ), 
                  
                  box(title = "Density plot of counts", 
                      status = "primary", 
                      plotOutput("plot2", height = 250)
                      )
                  
                    ),
                 
          
          fluidRow( box(title = "Counts of genes", 
                        "Table with gene counts per sample",
                        plotOutput("filtered_plot", height = 325)
                        ),
                    tabBox(height = "250px",
                           tabPanel("Principal components", "Variance explained by each component",
                                    plotOutput("varPCplot", height = 300)),
                           
                           tabPanel("PCA", "PCA analysis",
                                    selectInput("pc1", "Choose first PC", pcnames),
                                    selectInput("pc2", "Choose second PC", pcnames),
                                    plotOutput("scoreplot", height = 200)),
                           
                           tabPanel("Correlation", "Correlation plot of samples.",
                                    plotOutput("corrplot", height = 350)))
                    )
                ),
      
       
      
        tabItem(tabName = "dge_analysis",
                fluidRow(
                  box(
                    title = "Initialize DGE analysis", status = "warning",
                    "Please select the groups of your interest and click on the button below to initialise the DGE analysis.",
                    br(),
                    br(),
                    selectInput("selected_group", "Select a group of interest:", choices = NULL),
                    selectInput("selected_group2", "Select a group to compare against:", choices = NULL),
                    actionButton("initialise_button", "Initialize")
                  ),
                  box(title = "Dispersion plot",
                      "Graph for pairwise dispersion will be generated below",
                      br(),
                      plotOutput("dispersion_plot", height = 250)
                  )
                ),
                fluidRow(
                  box(
                    title = "Volcano plot of DEG", status = "warning",
                    "Please click on the button below to initialise the DGE analysis.",
                    br(),
                    plotOutput("volcano_plot", height = 350)
                  ),
                  box(
                    title = "Heatmap of top 20 DEGs", status = "warning",
                    "Please click on the button below to initialise the DGE analysis.",
                    br(),
                    plotOutput("heatmap", height = 350)
                  )
                )
                )
        )
    )
)


########## SERVER ############################

server <- function(input, output, session) {
  
  dane_loaded <- reactive ({
    req(input$file123, message = "Please upload a file first.")
    a <- read.csv(input$file123$datapath, header = TRUE)
    colnames(a)[1] = "ensgene"
    return(a)
  })
  
  ##############  CPM plots of individual genes - Page 1##################
  plotting <- function(dataa, symbole) {
    
    if (length(symbole) <= 10) {
      
      dataa %>% 
        filter(symbol %in% symbole) %>% 
        ggplot(aes(x = symbol, y = log2(count), color = group))+
        geom_boxplot()+
        labs(y = "log2 CPM")
      
    } else {
      print("Too many names")
    }
    
  }
  
  ################ Reactive data for PCA ######################
  summary1 <- reactive({
    file_content <- dane_loaded()
    
    symbolee <- rnor6 %>% 
      filter(ensgene %in% file_content$ensgene) %>% 
      select("ensgene", "symbol")
    
    elo2 <- merge(file_content, symbolee, by = "ensgene")
    elo2 <- elo2[!duplicated(elo2$symbol), ] 
    row.names(elo2) <- elo2$symbol
    elo2 <- elo2[, colnames(elo2) != "ensgene"]
    
    cpmed <- data.frame(cpm(elo2[, !(colnames(elo2) == "symbol")]))
    
    cpmed$symbol <- elo2$symbol
    return(cpmed)
  })
    
  
  ############### Menu items icons ######
  output$menuitem <- renderMenu({
    menuItem("Menu item", icon = icon("calendar"))
  })
  
  ########## PAGE 1 ###################
  
  output$filtered_plot <- renderPlot({
    req(input$file123, input$input_names, message = "Please upload the file.")
    
    
    final_data <- pivot_longer(summary1(), cols = -"symbol" , names_to = "samples", values_to = "count") %>% 
      mutate(group = sub("\\d+$", "", samples))
    
      names_to_filter <- as.vector(strsplit(input$input_names, ", ")[[1]])
      
      plotting(final_data, names_to_filter)
    })
  
  output$plot2 <- renderPlot({
    req(input$file123, message = "Please upload the file.")
    
    file_content <- dane_loaded()
    
    
    file1 <- pivot_longer(data = file_content, 
                          cols = 2:length(file_content), 
                          names_to = "sample", 
                          values_to = "counts")
  
    ggplot(file1, aes(x = counts, color = sample))+
    geom_density()+
    scale_x_continuous(limits = c(0,1000))
  })

  ####################### PCA_SUMMARY ####################
  
  PCA_data <- reactive({
    dane <- dane_loaded()
    dane <- dane[,-1]
    dane_norm <- dane[filterByExpr(dane),]
    sumy <- colSums(dane_norm)
    dane_norm1 <- apply(dane_norm, 1, FUN = function(CX) {log((CX+1)/(sumy+1))} )
    dane_norm1 <- scale(dane_norm1)
    pca_r<- prcomp(dane_norm1)
    
    X <- pca_r$x
    
    PCs <- data.frame(pca_r$x,
                      samples = sub("([A-Za-z]+)\\d+", "\\1", row.names(data.frame(pca_r$x))))
    
    
    variability <- (pca_r$sdev ^ 2) / sum(pca_r$sdev ^ 2) *100
    
    vari <- data.frame(PC = names(PCs[,1:12]),
                       variability)
    return(list(PC = PCs, var = vari, varia = variability, x = X))
  })
  
  
  ########## Score PC plotting function #####
  pc_plot <- function(data, x1, x2, vari) {
    x1_var <- as.integer(vari |> 
                           filter(PC == x1) |> 
                           pull(variability))
    
    x2_var <- as.integer(vari |> 
                           filter(PC == x2) |> 
                           pull(variability))
    
    ggplot(data, aes(x = data[[x1]], y = data[[x2]], color = samples))+
      geom_point(size = 4)+
      labs(x = paste0(x1,": ",x1_var, "%"), 
           y= paste0(x2,": ",x2_var, "%"))+
      theme_bw()
  }
   
  ###### Score plot and Var plot - Page 1 ##########
    output$scoreplot <- renderPlot({
    req(input$file123, message = "Please upload the file.")
  
  PCA_results <- PCA_data()
  
  pc_plot(data = PCA_results$PC, input$pc1, input$pc2, vari = PCA_results$var)
  
  })
  
  
  output$varPCplot <- renderPlot({
  req(input$file123, message = "Please upload the file.")
  
  PCA_results <- PCA_data()
  
  
  barplot(PCA_results$varia, xlab = "PC", ylab = "Variability explained %", 
          names.arg = c(1:12) )
  })
  
  ############# Correlation Plot - Page 1 ###########
  output$corrplot <- renderPlot({
    req(input$file123, message = "Please upload the file.")
   
    PCA_results <- PCA_data()
    
    correlation <- cor(t(data.frame(PCA_results$x)))
    heatmap(correlation)
  })
  
  
  ######### PAGE 2 #######
  disp_calc <- reactive({
    jazda <- dane_loaded()
    rownames(jazda) <- jazda[, 1]
    jazda <- jazda[, -1]
    group <- (gsub("([A-Za-z]+)([0-9]+)", "\\1", colnames(jazda)))
    design_matrix <- data.frame("Condition" = group)
    condition <- factor(design_matrix$Condition)
    model_matrix <- model.matrix(~condition)
    lecone  <- DGEList(counts = jazda, group = group)
    zatrzymaj <- filterByExpr(lecone)
    filtered_counts <- lecone[zatrzymaj, ]
    normalised <- normLibSizes(filtered_counts) 
    normalised$design <- model_matrix
    return(normalised)
  })
  
  
  
  observe({
    updateSelectInput(session, "selected_group", choices = unique(disp_calc()$samples$group))
    updateSelectInput(session, "selected_group2", choices = unique(disp_calc()$samples$group))
  })
  
  
  DGE1 <- reactiveVal(NULL)
  stat_temp <- reactiveVal(NULL)
  
  
  observeEvent(input$initialise_button, {
    dge <- estimateDisp(disp_calc(), design = disp_calc()$design)
    data_frame_DGE <- data.frame("Mean_counts" = rowMeans(dge$counts), 
                                 "Disp" = dge$tagwise.dispersion, 
                                 "Trended" = dge$trended.dispersion)
    
  DGE1(data_frame_DGE)
  stat_temp(dge)
  
              })
  
  DGE_data <- reactive({
    # Access the stored DGE data
    data_frame_DGE <- DGE1()
    
    # Check if data_frame_DGE is not NULL before returning
    if (!is.null(data_frame_DGE)) {
      return(data.frame(data_frame_DGE))
    } else {
      return(NULL)
    }
  })
  
  
      
  
  buttonClicks <- reactiveVal(0)
  
  # Event handler for the button click
  observeEvent(input$initialise_button, {
    # Increment the button clicks counter
    buttonClicks(buttonClicks() + 1)
  })
  
  
  output$dispersion_plot <- renderPlot({
    lol <- DGE_data()
    req(input$file123)
    if (buttonClicks() > 0){
    ggplot(lol, aes(x = log10(Mean_counts), y = log10(Disp)))+
      geom_point(aes(alpha = 0.01))+
      geom_line(aes(x = log10(Mean_counts), y = log10(Trended), color = "red"), 
                linewidth = 1.2)+
      labs(x = "log10 mean gene counts", y = "log10 tagwise dispersion")
    } else {
      print("ehh")
    }
  })
  
  
  
  stat_test <- reactive({
    # Access the stored DGE data
    stat_123 <- stat_temp()
    
    # Check if data_frame_DGE is not NULL before returning
    if (!is.null(stat_123)) {
      return(stat_123)
    } else {
      return(NULL)
    }
  })
  
  
  selected_values <- reactiveValues(value1 = NULL,
                                    value2 = NULL)
  
  observeEvent(input$selected_group, {
    selected_values$value1 <- input$selected_group
  })
  
  observeEvent(input$selected_group2, {
    selected_values$value2 <- input$selected_group2
  })
  
  
  
  labelsa <- reactiveVal(NULL)
  
  
  
  output$volcano_plot <- renderPlot({
    req(input$file123)
    if (!is.null(stat_test())) {
      fit <- glmQLFit(stat_test(), stat_test()$design)
      
      lata <- data.frame(groups = unique(fit$samples$group))
      
      vec = rep(0, length(lata$groups) )
      vec[which(lata$groups == selected_values$value1)] = 1
      vec[which(lata$groups == selected_values$value2)] = -1
      
      
      if (which(lata$groups == selected_values$value2) == 1) {
        xvsc <- data.frame(topTags(glmQLFTest(fit, 
                                              coef = which(lata$groups == selected_values$value1)), 
                                   n = 50000))
      } else {
        xvsc <- data.frame(topTags(glmQLFTest(fit, contrast = vec), n = 50000))
      }
      
      
      xvsc$ensgene <- row.names(xvsc)
      
      names <- rnor6 %>% 
        filter(ensgene %in% row.names(xvsc)) %>% 
        select(1,3)
      xvsc <- merge(xvsc, names, by = "ensgene")
      
      xvsc$diffexpressed <- "Not regulated"
      xvsc$diffexpressed[xvsc$logFC > 0.1 & xvsc$FDR < 0.05] <- "Upregulated"
      xvsc$diffexpressed[xvsc$logFC < -0.1 & xvsc$FDR < 0.05] <- "Downregulated"
      
      
      labels <- xvsc %>% 
        arrange(desc(logFC)) %>% 
        select("symbol")
      
      labels <- rbind(
        head(labels, 10),
        tail(labels, 10)
      )
      
      labelsa(labels)
      
      xvsc$delabel <- ifelse(xvsc$symbol %in% labels$symbol, xvsc$symbol, NA)
      
      ggplot(xvsc, aes(x = logFC, y = abs(log10(FDR)), col = diffexpressed))+
        geom_point(alpha = 0.5, size = 2)+
        labs(y ="-log10 adjusted pvalue", legend = FALSE)+
        geom_vline(xintercept= 0, col="grey", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="grey", linetype = "longdash") + 
        scale_color_manual(values=c("#21908CFF", "grey", "#FFA500"))+
        geom_text_repel(aes(label = delabel), color = "black")+
        theme_bw() +
        theme(axis.title.y = element_text(face = "bold"),
              axis.title.x = element_text(face = "bold"))+
        guides(color = guide_legend(title = NULL))
      } else {
      return(NULL)
    }
    
  })
  
  
  labels_react <- reactive({
    # Access the stored DGE data
    labe <- labelsa()
    
    # Check if data_frame_DGE is not NULL before returning
    if (!is.null(labe)) {
      print(labe)
      return(data.frame(labe))
    } else {
      return(NULL)
    }
  })
  
  
  
  output$heatmap <- renderPlot({
    file_content <- dane_loaded()
    labels <- labels_react()
    
    if (!is.null(stat_test())) {
    
    symbolee <- rnor6 %>% 
      filter(ensgene %in% file_content$ensgene) %>% 
      select("ensgene", "symbol")
    
    elo2 <- merge(file_content, symbolee, by = "ensgene")
    elo2 <- elo2[!duplicated(elo2$symbol), ] 
    row.names(elo2) <- elo2$symbol
    elo2 <- elo2[, colnames(elo2) != "ensgene"]
    
    cpmed <- data.frame(cpm(elo2[, !(colnames(elo2) == "symbol")]))
    
    head(cpmed)
    zatrzymane <- t(scale(t(cpmed)))
    zatrzymane <- zatrzymane[row.names(zatrzymane)  %in% labels$symbol, ]
    heatmap(zatrzymane)
    }
  })
 
  }

############## UI SERVER APP ############
shinyApp(ui = ui, server = server)
