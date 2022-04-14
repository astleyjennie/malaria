source("R/packages.R")
source("R/model.R")

# Define UI for application 
ui <- fluidPage(
  
  # Title
  titlePanel("Malaria in Thailand"),
  
  # Sidebar with slider inputs 
  sidebarLayout(
    sidebarPanel(
      actionButton("go", "Go"),
      sliderInput(inputId="t1", label = "Migration", value = 1, min=0, max=1,step=1),
      sliderInput(inputId="t6", label = "Border testing", value = 0, min=0, max=1,step=1),
      sliderInput(inputId="t7", label = "G6PD testing", value = 0, min=0, max=1,step=1),
      sliderInput(inputId="t8", label = "Treat imported cases", value = 0, min=0, max=1,step=1),
      
      bsCollapse(
        # Mosquito parameters
        bsCollapsePanel("Transmission parameters",
                        sliderInput(inputId="a", label = "human feeding rate per mosquito", value = 0.3, min=0, max=10,step=0.01),
                        sliderInput(inputId="b", label = "transmission efficiency mosquito to human", value = 0.5, min=0, max=1,step=0.01),
                        sliderInput(inputId="c", label = "transmission efficiency human to mosquito", value = 0.5, min=0, max=1,step=0.01),
                        sliderInput(inputId="m", label = "Mosquitoes per human", value = 0.8, min=0.1, max=3,step=0.1),
                        sliderInput(inputId="relv", label = "Relative infectiousness of vivax compared to falciparum", value = 0.29, min=0.1, max=1,step=0.01),
                        sliderInput(inputId="nu2", label = "Migration rate: leaving Thailand", value = 0.001, min=0, max=0.002,step=0.0005),
                        sliderInput(inputId="nu3", label = "Emigration rate: entering Thailand", value = 0.000285, min=0.0001, max=0.0004,step=0.0001)
                        
        )
      )
    ),
    
    # Tabs for incidence and costs
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Cases", plotOutput("modelPlot") %>% withSpinner()),
                  tabPanel("Incidence", plotOutput("incPlot") %>% withSpinner()),
                  tabPanel("RDT Cost", plotOutput("rdtPlot") %>% withSpinner()),
                  tabPanel("G6PD Screening Cost", plotOutput("g6pdPlot") %>% withSpinner()),
                  tabPanel("ACT Cost", plotOutput("actPlot") %>% withSpinner()),
                  tabPanel("Primaquine cost", plotOutput("primPlot") %>% withSpinner())
      )
    )
  )
)

# Server
server <- function(input, output) {
  modelOut <- eventReactive(input$go, {
    # Overwrite the parms with inputs
    parameters["t1"] <- input$t1
    parameters["t6"] <- input$t6
    parameters["t7"] <- input$t7
    parameters["t8"] <- input$t8
    parameters["a"] <- input$a
    parameters["b"] <- input$b
    parameters["c"] <- input$c
    parameters["m"] <- input$m
    parameters["relv"] <- input$relv
    parameters["nu2"] <- input$nu2
    parameters["nu3"] <- input$nu3
    #parms["gamma_m"] <- input$gamma_m
    #parms["mu_m"] <- input$mu_m
    
    # Run the model
    result <- ode(times=times, y=istate, func=thailand_model, parms=parameters) %>% 
    as.data.frame() %>%
    as_tibble() %>%
    mutate(PThaif = (Sf+Ef+Af+Cf+Tf+Rf+E2f),
           PThaiv = (Sv+Ev+Av+Cv+TPv+TAv+Lv+Rv+E2v),
           Pf=(Outf+Testf+PThaif),
           Pv=(Outv+Testv+PThaiv),
           If = (Ef+Af+Cf+Tf),
           Iv = (Ev+Av+Cv+TPv+TAv),
           P = PThaif+PThaiv,
           CTrtv=CACT+CPrim,
           Incf = c(0, diff(CIncf)),
           Incv = c(0, diff(CIncv)),
           Trtf = c(0, diff(CTrtf)),
           Trtv = c(0, diff(CTrtv)),
           CostPrim = c(0, diff(CPrim)),
           CostACT = c(0,diff(CACT)),
           CostG6PD = c(0, diff(CCv)),
           CostRDTf = c(0, diff(CTestf)),
           CostRDTv = c(0, diff(CTestv))) %>% 
      pivot_longer(names_to = "variable", cols = !1) %>%
      mutate(SP = ifelse(str_ends(variable, "f"), "Pf", "Pv")
      )
  })
  
  #Define Plots
  
  predicted_cases <- as.data.frame(seq(from=2010,to=2037,by=1))
  
  output$modelPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("If", "Iv")) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Cases", y =("population"), colour="Compartment")+
      facet_wrap(~SP)
  })
  
  output$incPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("Incf", "Incv")) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Incidence", y =("population"), colour="Compartment")+
      facet_wrap(~SP)
  })

  output$primPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("CostPrim"), time>100) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Daily Primaquine Cose", y =("USD"), colour="Compartment")+
      facet_wrap(~SP)
  })
  
  output$actPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("CostACT"), time>100) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Daily ACT Cost", y =("USD"), colour="Compartment")+
      facet_wrap(~SP)
  })
  
  output$rdtPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("CostRDTf"), time>100) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Daily RDT Cost", y =("USD"), colour="Compartment")+
      facet_wrap(~SP)
  })
  
  output$g6pdPlot <- renderPlot({
    modelOut() %>% 
      filter(variable %in% c("CostG6PD"), time>100) %>% 
      group_by(variable) %>%
      ggplot()+
      geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
      theme_minimal() +
      labs(title = "Daily G6PD Screening Cost", y =("USD"), colour="Compartment")+
      facet_wrap(~SP)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
