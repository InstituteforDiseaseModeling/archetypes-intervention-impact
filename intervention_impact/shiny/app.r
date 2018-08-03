
## prototype shiny app for intervention impact

library(shiny)
library(data.table)
library(ggplot2)

theme_set(theme_minimal(base_size = 18))

# main_dir <- file.path(Sys.getenv("USERPROFILE"), 
#                       "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables")
data <- fread(file.path("data", "lookup_ind_ints_het_biting.csv"))
data <- data[(is.na(Hates_Nets) | Hates_Nets==0) & (is.na(ACT_HS_Rate) | ACT_HS_Rate==0.15),
             list(Site_Name, Run_Number, Intervention, Coverage, x_Temporary_Larval_Habitat, initial_prev, final_prev, mean_initial, mean_final)]


ui <- fluidPage(
  
  plotOutput("curves", height="600px"),
  
  hr(),
  
  fluidRow(
    column(3,
           checkboxGroupInput("site",
                              h3("Sites"),
                              choices=list("Aba (Nigeria)"="aba",
                                           "Bajonapo (Ecuador)"="bajonapo",
                                           "Djibo (Burkina Faso)"="djibo",
                                           "Gode (Ethiopia)"="gode",
                                           "Kananga (DRC)"="kananga",
                                           "Karen (Myanmar)"="karen",
                                           "Martae (Cameroon)" = "martae",
                                           "Moine (Mozambique)"="moine"),
                              selected=c("aba", "kananga", "moine"))
           ),
    column(2, 
           checkboxGroupInput("ITN",
                              h3("ITN (%)"),
                              choices=list("0"=0,
                                           "20"=0.2,
                                           "40"=0.4,
                                           "60"=0.6,
                                           "80"=0.8),
                              selected=0.2)
           ),
    
    column(2,
           checkboxGroupInput("IRS",
                              h3("IRS (%)"),
                              choices=list("0"=0,
                                           "20"=0.2,
                                           "40"=0.4,
                                           "60"=0.6,
                                           "80"=0.8),
                              selected=0.4)
    ),
    
    column(2,
           checkboxGroupInput("ACT",
                              h3("ACT (%)"),
                              choices=list("0"=0,
                                           "20"=0.2,
                                           "40"=0.4,
                                           "60"=0.6,
                                           "80"=0.8),
                              selected=0.6)
    ),
    
    column(3,  
           radioButtons("interactions", 
                        h3("Include Interactions"),
                        choices=list("Yes"=T,
                                     "No"=F),
                        selected=F),
           radioButtons("random_seeds",
                        h3("Show All Runs"),
                        choices=list("Yes"=T,
                                     "No"=F),
                        selected=F)
           )
    
  )
  
  
)

server <- function(input, output){
  
  plotData <- reactive({
    subset <- data[((Intervention=="ITN" & Coverage %in% input$ITN) |
                    (Intervention=="IRS" & Coverage %in% input$IRS) |
                    (Intervention=="ACT" & Coverage %in% input$ACT)) & 
                   (Site_Name %in% input$site)]
    subset
  })
  
  
  
  output$curves <- renderPlot({
    out_plot <- ggplot(plotData(), aes(x=mean_initial, y=mean_final, color=Coverage)) +
                  geom_line(aes(linetype=Intervention), size=1) +
                  geom_abline() + 
                  facet_wrap(~Site_Name) +
                  labs(x="Initial Prevalence",
                       y="Final Prevalence")
    if (input$random_seeds){
      out_plot <- out_plot + geom_line(aes(x=initial_prev, y=final_prev, group=interaction(Intervention, Coverage, Run_Number)), alpha=0.25)
    }
    print(out_plot)
  })
  
}

shinyApp(ui, server)