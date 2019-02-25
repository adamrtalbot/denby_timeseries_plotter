#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(shiny)
library(scales)
library(RSQLite)

########################################################
# Define UI for application that plots gene expression #
########################################################

ui <- fluidPage(
  
   # Application title
  titlePanel("Denby Lab Timeseries Data Plot"),
   
   # Sidebar with input fields 
   verticalLayout(
      #inputPanel(
         textAreaInput(inputId = "genes", 
                       label = "Genes:", 
                       value = "Lsat_1_v5_gn_4_1201, Bcin02g05780", 
                       resize = "both"),
         br(),
         selectInput('pathosystem', 
                     "Timeseries data:", 
                     c("Lettuce-Botrytis" = "letbot", 
                       "Lettuce-Sclerotinia" = "letsclero", 
                       "Arabidopsis-Botrytis" = "athalbot"), 
                     multiple = FALSE),
         selectInput('dataformat', 'Data format:',
                     c('Counts' = 'counts',
                       'Log(2) + 1' = 'log2',
                       'VST' = 'vst'),
                     multiple = FALSE),
         checkboxInput("pointButton", "Points", TRUE),
         checkboxInput("smoothButton", "Loess Regression", FALSE),
         checkboxInput("meanButton", "Mean and SE", TRUE),
         checkboxInput("fixyaxis", "Fix y-Axis?", FALSE),
     # ),
      
      # Show a plot of the generated distribution
     
     plotOutput("genePlot"),
     downloadButton(outputId = 'downloadPlot', label = "Download Plot"),
     downloadButton(outputId = 'downloadData', label = "Download CSV")
   )
)


#############################################
# Define server logic required to draw plot #
#############################################

server <- function(input, output) {
  
  timeseries_db <- dbConnect(drv = RSQLite::SQLite(), "../data/timeseries-db.sqlite")
  
  output$genePlot <- renderPlot({
    
    # Make gene_list for further use:
    gene_list <- input$genes %>%
      gsub(pattern = "\n", replacement = " ", x = ., perl = TRUE) %>%
      gsub(pattern = ",", replacement = " ", x = ., perl = TRUE) %>%
      gsub(pattern = "[ ]+", replacement = " ", x = ., perl = TRUE)
    gene_list <- strsplit(x = gene_list, split = " ")[[1]]
    
    # Retrieve data from SQL database using parameterised query:
    
    rs <- dbSendQuery(timeseries_db, paste('SELECT * FROM', input$pathosystem, 'WHERE gene in (:genes)'))
    dbBind(rs, param = list(genes = gene_list))
    select.data <- dbFetch(rs)
    
    # Select the expression data (since the array doesn't have counts or VST)
    
    if(input$pathosystem != "athalbot"){
      select.data <- select.data %>% mutate_("expression" = input$dataformat)
    }
    
    # Sort data by input
    
    select.data <- select.data %>% mutate(gene = factor(gene, levels = gene_list))
    
    # Create plot for data:
    
    gene.plot <-  ggplot(data = select.data, aes(x = timepoint, 
                                                 y = expression, 
                                                 colour = treatment, 
                                                 group = treatment, 
                                                 fill = treatment)) +
      ylab("Expression") +
      xlab("Hours Post Infection (hpi)") +
      scale_y_continuous(breaks = pretty_breaks(n = 8), 
                         expand = c(0.001, 1)) +
      scale_x_continuous(limits = c(min(select.data$timepoint), max(select.data$timepoint)), 
                         breaks = seq(9, 54, by = 3)) + 
      labs(colour = "Treatment", 
           group = "Treatment", 
           fill = "Treatment") +
      facet_wrap(~gene, scales = "free") +
      theme(strip.background = element_blank())

    if(input$pointButton == TRUE){
      gene.plot <- gene.plot + geom_point()
    }
    
    if(input$smoothButton == TRUE){
      gene.plot <- gene.plot + geom_smooth(method = "loess", 
                                           mapping = aes(fill = treatment))
    }
    
    if(input$meanButton == TRUE){
      gene.plot <- gene.plot + 
        stat_summary(fun.y = "mean", 
                     fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                     fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                     geom = "ribbon", alpha = 0.3) + 
        stat_summary(fun.y = "mean", geom = "line")
    }
    
    if(input$dataformat == "log2"){
      gene.plot <- gene.plot + ylab(bquote(~Log[2]~" (Count + 1)"))
    } else if (input$dataformat == "count") {
      gene.plot <- gene.plot + ylab("Count")
    } else if (input$dataformat == "vst") {
      gene.plot <- gene.plot + ylab("VST counts")
    }
    
    if(input$pathosystem == "athalbot"){
      gene.plot <- gene.plot + 
        scale_x_continuous(limits = c(min(select.data$timepoint), max(select.data$timepoint)),
                           breaks = seq(0, 48, by = 2)) +
        ylab("Relative Expression")
    }
    
    if(input$fixyaxis == TRUE){
      gene.plot <- gene.plot +
        scale_y_continuous(limits = c(floor(min(select.data$expression)), ceiling(max(select.data$expression))))
    }
    
    out.data <- out.data <- select.data %>% 
      arrange(timepoint) %>%
      mutate(sample = paste0(treatment, "_tp", timepoint, "_rep", replicate)) %>% 
      select(gene, sample, expression) %>% 
      spread(key = sample, value = expression, fill = NA)              
            
    output$downloadPlot <- downloadHandler(filename = "plot.png", 
                                           content = function(file){
                                             ggsave(file, plot = gene.plot, device = "png", width = 16)
                                           }
    )
    
    output$downloadData <- downloadHandler(
      filename = "gene_exp.csv", 
      content = function(file){
        write.csv(x = out.data, file, row.names = FALSE)
      }
    )
    gene.plot

  }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

#####################################
# To do:
# - Conversion of ATG to Lsat
# - Use plotly for plotting
