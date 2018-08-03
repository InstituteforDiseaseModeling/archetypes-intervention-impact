
## prototype shiny app for intervention impact

library(shiny)
library(data.table)


main_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables")
data <- fread(file.path(main_dir, "lookup_ind_ints_het_biting.csv"))
data <- data[(is.na(Hates_Nets) | Hates_Nets==0) & (is.na(ACT_HS_Rate) | ACT_HS_Rate==0.15),
             list(Site_Name, Run_Number, Intervention, Coverage, x_Temporary_Larval_Habitat, initial_prev, final_prev, mean_initial, mean_final)]


ui <- fluidPage(
  
  titlePanel("Intervention Impact"),
  
  sidebarLayout(
    sidebarPanel(
      
      checkboxGroupInput("ITN",
                         h3("ITN"),
                         choices=list("0"=0,
                                      "20"=0.2,
                                      "40"=0.4,
                                      "60"=0.6,
                                      "80"=0.8),
                         selected=0.2),
      
      checkboxGroupInput("IRS",
                         h3("IRS"),
                         choices=list("0"=0,
                                      "20"=0.2,
                                      "40"=0.4,
                                      "60"=0.6,
                                      "80"=0.8),
                         selected=0.2),
      
      checkboxGroupInput("ACT",
                         h3("ACT"),
                         choices=list("0"=0,
                                      "20"=0.2,
                                      "40"=0.4,
                                      "60"=0.6,
                                      "80"=0.8),
                         selected=0.2)
      
    ),
    mainPanel(plotOutput("curves")
              )
  )
  
)

server <- function(input, output){
  
  output$curves <- renderPlot({
    ggplot(data[(Intervention=="ITN" & Coverage==input$ITN) |
                (Intervention=="IRS" & Coverage==input$IRS) |
                (Intervention=="ACT" & Coverage==input$ACT)], aes(x=mean_initial, y=mean_final, color=Coverage)) +
      geom_line(aes(linetype=Intervention), size=1) +
      facet_wrap(~Site_Name)
  })
  
}

shinyApp(ui, server)