
## prototype shiny app for intervention impact

library(shiny)
library(data.table)
library(rsconnect)
library(ggplot2)

theme_set(theme_minimal(base_size = 18))

data <- fread(file.path("data", "lookup_full_interactions_v3.csv"))

data[, names:=gsub("[0-9]", "", Intervention)]

int_choices = list("20"=0.2,
                   "40"=0.4,
                   "60"=0.6,
                   "80"=0.8)

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
                              selected=c("djibo", "kananga", "karen"))
           ),
    column(2, 
           checkboxGroupInput("ITN",
                              h3("ITN (%)"),
                              choices=int_choices,
                              selected=0.2)
           ),
    
    column(2,
           checkboxGroupInput("IRS",
                              h3("IRS (%)"),
                              choices=int_choices,
                              selected=0.4)
    ),
    
    column(2,
           checkboxGroupInput("ACT",
                              h3("ACT (%)"),
                              choices=int_choices,
                              selected=0.6)
    ),
    
    column(3,  
           radioButtons("interactions", 
                        h3("Include Interactions"),
                        choices=list("Yes"=T,
                                     "No"=F),
                        selected=T),
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
    
    if (input$interactions){
      
      combos <- data.table(expand.grid(ITN=c(0, input$ITN), IRS=c(0, input$IRS), ACT=c(0, input$ACT)))
      combos[, label:=""]
      
      for (int in c("ITN", "IRS", "ACT")){
        combos[get(int)!=0, label:= paste0(label, int, " ", get(int), "; ") ]
      }
      combos[label=="", label:="None"]
      combos[, label:=gsub(" $", "", label)]
      
      subset <- data[(Intervention %in% combos$label) &
                      (Site_Name %in% input$site)
                     ]
      
    }else{
      subset <- data[((ITN_Coverage %in% input$ITN & IRS_Coverage==0 & ACT_Coverage==0) |
                        (IRS_Coverage %in% input$IRS & ITN_Coverage==0 & ACT_Coverage==0) |
                        (ACT_Coverage %in% input$ACT & IRS_Coverage==0 & ITN_Coverage==0)) &
                       (Site_Name %in% input$site)]
    }
    
    # set color order
    for_colors <- subset[x_Temporary_Larval_Habitat== unique(subset$x_Temporary_Larval_Habitat)[40] &
                           Run_Number==0 &
                           Site_Name==input$site[1] ]
    for_colors <- for_colors[order(mean_final, decreasing=T)]
    
    subset[, Intervention:= factor(Intervention, levels=for_colors$Intervention)]
    
    subset
    
  })
  
  
  
  output$curves <- renderPlot({
    out_plot <- ggplot(plotData(), aes(x=mean_initial, y=mean_final, color=Intervention)) +
                  geom_line(aes(linetype=names), size=1.25) +
                  scale_linetype(guide="none") + 
                  geom_abline() + 
                  facet_wrap(~Site_Name) +
                  labs(x="Initial Prevalence",
                       y="Final Prevalence")
    if (input$random_seeds){
      out_plot <- out_plot + geom_line(aes(x=initial_prev, y=final_prev, group=interaction(Intervention, Run_Number)), alpha=0.25)
    }
    print(out_plot)
  })
  
}

shinyApp(ui, server)