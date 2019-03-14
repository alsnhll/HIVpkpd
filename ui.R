library(shiny)
library(ggplot2)
library(plotly)
library(xtable)
require(visNetwork)
shinyUI(fluidPage(
  titlePanel("HIV Antiretroviral Therapy Visualization"),
  sidebarPanel(
    h3("Parameter Selection:",align="center"),
      helpText("*This tool is for research purposes only and is not intended for patient care."),
             selectInput("drugs",
                         label = "Select drug to visualize:",
                         choices = c("3TC", "ABC", "ATV", "ATV/r", "AZT", "CAB", "d4T", "ddI", "DRV/r", "DTG", 
                                     "EFV", "ENF", "ETV", "EVG", "FTC", "IDV", "IDV/r","LPV/r", "NFV", "NVP",
                                     "RAL", "RPV", "SQV", "SQV/r", "TDF", "TPV/r"), 
                         selected = c("3TC")
             ),
    conditionalPanel(condition = "input.conditionedPanel == 7|| input.conditionedPanel != 3 && input.conditionedPanel !=6 && input.conditionedPanel !=7",
             radioButtons("units", "Concentration units:",
                          choices = list("µM" = "µM","µg/ml" = "µg/ml","ng/ml"="ng/ml"),selected = "µM"
                          )),
    
    conditionalPanel(condition = "input.conditionedPanel == 1| input.conditionedPanel == 3",      
          h4("Pharmacokinetic visualization:",align="center"),
                 selectInput("model", 
                             label = "Initial condition:",
                             choices = c("No Drug", "Steady State"),
                             selected = "No Drug"
                             ),
    numericInput("doses", 
                label = "Number of doses:",
                value = 1
                )),
    conditionalPanel(condition = "input.conditionedPanel ==7 || input.conditionedPanel != 1 && input.conditionedPanel !=5 && input.conditionedPanel !=7&& input.conditionedPanel !=6",
                     h4("Pharmacodynamic visualization:",align="center"),
                     

                 sliderInput("R0",
                             HTML(paste("Baseline fitness of wild-type virus (R", tags$sub(0),'):', sep = "")),
                             min = 1.25,
                             max = 25,
                             value = 10
                             )),
    conditionalPanel(condition = "input.conditionedPanel == 2 | input.conditionedPanel == 3",
                     radioButtons("fitness_disp","Display viral fitness as:",
                                  c("Intantaneous viral fitness" = 'r0',
                                    'Instantaneous inhibitory potential (IIP)' = 'iip')
                     ) ),
    
    conditionalPanel(condition = "input.conditionedPanel == 2",
                     
                     sliderInput("conc_range",
                                 HTML(paste("Log",tags$sub(10),' concentration range:',sep="")),
                                 min = -5, max = 5, value = c(-3, 3)
                     )),
    conditionalPanel(condition = "input.conditionedPanel !=1 && input.conditionedPanel !=5 && input.conditionedPanel != 6 && input.conditionedPanel !=7",
    h4("Resistant strains to visualize:",align="center"),
                 selectInput("show_muts",
                             label = "Select mutant strains to show:",
                             choices = c("None","All Known Mutants","Single Known Mutant","Test Mutant"),
                             selected = "All Known Mutants"
                             )),
  
   conditionalPanel(
      condition="input.show_muts !== 'None'&& input.show_muts !== 'Test Mutant' && input.conditionedPanel !=1 && input.conditionedPanel !=5 && input.conditionedPanel != 6 && input.conditionedPanel !=7",
      checkboxInput("disp_nas",label = "Include mutants with incomplete data?",
                    value= TRUE)
      ),
  
  conditionalPanel(
    condition = "input.conditionedPanel !=1 && input.show_muts == 'Single Known Mutant'",
    selectInput("mut_select",
                label = "Choose a mutant strain to visualize:",
                "")),
  conditionalPanel(condition = "input.conditionedPanel != 1 && input.conditionedPanel !=5 && input.conditionedPanel != 6 && input.conditionedPanel !=7",
                   sliderInput("maxfit","Maximum allowed relative fitness of mutant strain (if greater than WT):",
                               min=0,max=0.99,value = 0.95)),

  conditionalPanel(
     h4("Parameters to assume for missing data:",align="center"),
     
     condition = "input.conditionedPanel ==7|| input.conditionedPanel !=1 && input.show_muts !== 'None'&& input.disp_nas==true ||input.conditionedPanel !=1 && input.show_muts == 'Test Mutant'",
     
     numericInput("sigma",
                  "Fractional change in slope:",
                  value = 0)
     ),

  conditionalPanel(
     condition = "input.conditionedPanel ==7||input.conditionedPanel !=1 && input.show_muts !== 'None'&& input.disp_nas==true||input.conditionedPanel !=1 && input.show_muts == 'Test Mutant'",
     numericInput("s",
                  "Relative fitness (0 < x < 1):",
                  value = 0.9)
     ),

  conditionalPanel(
     condition = "input.conditionedPanel ==7||input.conditionedPanel !=1  && input.show_muts !== 'None'&& input.disp_nas==true||input.conditionedPanel !=1 && input.show_muts == 'Test Mutant'",
     numericInput("rho",
                  HTML(paste("Fold change in IC", tags$sub(50),':', sep = "")), #""
                  value = 2)
     ),
  conditionalPanel(condition = "input.conditionedPanel ==7",
                   h4("Parameters for mutation matrix construction / visualization:",align="center"),
                   checkboxInput("selfmut",label = "Display self-transitions?",
                                   value=FALSE),
                   checkboxInput("backmut",label = "Allow back mutations?",
                                 value= FALSE),
                   checkboxInput("directmulti",label = "Allow instantaneous multi-hit mutations?",
                                 value= FALSE),
                   checkboxInput("directmultmult",label = "Allow mutants to directly mutate to one another?",
                                 value= FALSE)
  ),

      
      width = 3
    ),
  
  
  
  mainPanel(
    navbarPage("Select output:",
               navbarMenu("Figures:",
                tabPanel("Pharmacokinetic Profile",br(),plotlyOutput("plotPK"),value = 1),
                
                tabPanel("Pharmacodynamic Profile",plotlyOutput("drPlot"),br(),br(),br(),br(),br(),br(),br(),br(),verbatimTextOutput("warning"),value = 2), 
                
                tabPanel("PK/PD Profile",plotlyOutput("plotpd"),br(),br(),br(),br(),br(),br(),br(),br(),verbatimTextOutput("warning2"),value = 3),
                
               tabPanel("Mutation Network Visualization",br(),
                        fluidRow(column(12,align = "center",
                                        verbatimTextOutput("visText"),
                                        verbatimTextOutput('param_warning'),
                                        verbatimTextOutput("viswarning"),
                                        
                                        visNetworkOutput('network'))),br(),br(),value =7),
               
                tabPanel("Selection Windows", plotlyOutput("mswPlot"),br(),br(),br(),br(),br(),br(),br(),br(),verbatimTextOutput("warning3"),br(),
                         radioButtons("window_units", "Axis Units:",
                                      choices = list("Time since last dose" = "days", "Drug concentration" = "concs"),selected = "days"),
                         value = 4)),
              
    
    navbarMenu("Tables:",
               tabPanel("PK/PD Values", br(),
                                  fluidRow(column(12,align = "center",
                                                  h3("Pharmacokinetic and pharmacodynamic parameters:"),br(),
                                                  tableOutput("drugvalues"),br(),br(),
                                                  
                                                  h3("Parameters calculated from one-compartment pharmacokinetic model:",align="center"),br(),
                                                  tableOutput("fittedparams"),br()
                                  )),value = 5),
               
               tabPanel("Resistant Strain Data", br(),
                        fluidRow(column(12,align = "center",
                                        h3("Resistant strain data:"),
                                        tableOutput("mutvalues"))),value = 6),
               tabPanel("Mutation Matrix", br(),
                        fluidRow(column(12,align = "center",
                                        h3("Mutation Matrix (Q):"),
                                        tableOutput('Qmatrix'))),br(),br(),verbatimTextOutput("matwarning"),value = 7)
               # ,id = "conditionedPanel")
               )
    ,id = "conditionedPanel")
  )
  
  
  
  
  
  
  
  
  
  
  
  
  )
)