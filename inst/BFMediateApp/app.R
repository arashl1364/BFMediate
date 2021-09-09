#' @include Mediate.R
NULL

# mediation image
img <- "www/mediation.png"

# simulation function (partial mediation)
simPartialMed = function(beta_1,beta_2, sigma_M, sigma_Y,N,X) {
  eps_M = rnorm(N)*sigma_M
  eps_Y = rnorm(N)*sigma_Y
  M = beta_1[1] + beta_1[2] * X + eps_M # generate mediator M
  Y = beta_2[1] + beta_2[2] * M + beta_2[3] * X + eps_Y # generate dependent variable
  list(X = X, M = M, Y = Y)
}

shiny::shinyApp(

  # Define UI
  ui = shiny::fluidPage(title="BFMediate",
    withMathJax(),
    # Application title
    shiny::h4("Mediation Analysis using Bayes Factors"),
    # Sidebar
    shiny::sidebarLayout(
    #shiny::fluidRow(
      shiny::sidebarPanel(
      #shiny::column(4,
        shiny::numericInput(inputId= 'seed', label='Seed (Reproducibility)',
                            value=42),
        shiny::selectInput(inputId='model', label='Choose model type',
                           choices=c('Simple mediation model'='Simple',
                                     'Latent variable model'='Latent')
                          ),
        shiny::helpText('$$M=\\beta_{0M}+\\beta_{1}X+\\varepsilon_{M}$$'),
        shiny::helpText('$$Y=\\beta_{0Y}+\\beta_{2}M+\\beta_{3}X+\\varepsilon_{X}$$'),
        shiny::hr(),
        # Simulation parameters ------------------------------------------------
        shiny::h4('Simulation Parameters (DGP)'),
        column(12,
          shiny::numericInput(inputId='N', label='Number of observations',
                              value=1000)
        ), column(6,
          shiny::numericInput(inputId="sigma_M", label = "eps_M",
                              value=1^.5),
        ), column(6,
          shiny::numericInput(inputId="sigma_Y", label = "eps_Y",
                              value=1^.5),
        ), column(6,
        shiny::numericInput(inputId="beta_0M", label = "beta_0M",
                            value=1),
        ), column(6,
        shiny::numericInput(inputId="beta_1", label = "beta_1",
                            value=.5),
        ), column(6,
        shiny::numericInput(inputId="beta_0Y", label = "beta_0Y",
                            value=1),
        ), column(6,
        shiny::numericInput(inputId="beta_2", label = "beta_2",
                            value=1.5),
        ), column(12,
        shiny::numericInput(inputId="beta_3", label = "beta_3 (direct effect)",
                            value=0),
        ),
        shiny::hr(),
        shiny::actionButton(inputId = "simulate", label = "Simulate data",
                            # database
                            # angle-right
                            icon("database"), class = "btn-info"),
        shiny::hr(),
        # Priors ---------------------------------------------------------------
        shiny::h4('Priors'),
        column(6,
        shiny::numericInput(inputId='prior_0M', label='Prior beta_0M',
                            value=100),
        ), column(6,
        shiny::numericInput(inputId='prior_1', label='Prior beta_1',
                            value=100),
        ), column(6,
        shiny::numericInput(inputId='prior_0Y', label='Prior beta_0Y',
                            value=100),
        ), column(6,
        shiny::numericInput(inputId='prior_2', label='Prior beta_2',
                            value=100),
        ), column(12,
        shiny::numericInput(inputId='prior_3', label='Prior beta_3 (reference prior)',
                            value=1),
        ),
        shiny::hr(),
        # Estimation parameters ------------------------------------------------
        shiny::h4('Estimation Parameters'),
        shiny::numericInput(inputId='R', label='MCMC Iterations',
                            value=10000),
        shiny::numericInput(inputId='burnin', label='Burn-In (Warm-Up)',
                            value=2000),
        shiny::hr(),
        shiny::actionButton(inputId = "estimatemodel", label = "Estimate model",
                            # calculator
                            # angle-double-right
                            # chevron-circle-right
                            # forward
                            icon("calculator"), class = "btn-warning"),
        shiny::br(),
        shiny::tags$small(shiny::HTML("&nbsp;&nbsp; (might take a while)"))
      ),

      # main panel
      shiny::mainPanel(
      #shiny::column(8,
        shiny::tags$style(type="text/css", "img{display: block; margin-left: auto; margin-right: auto;}"),
        # show mediation image
        shiny::imageOutput(outputId = "img", height = "100%"),
        shiny::br(),
        # shiny::h4('Data summary'),
        shiny::htmlOutput('summary_heading'),
        # Verbatim text for data summary
        shiny::verbatimTextOutput("summary"),
        shiny::br(),
        # data table (interactive)
        DT::dataTableOutput("df"),
        # shiny::dataTableOutput("df")
        shiny::br(),
        shiny::htmlOutput('results_heading'),
        shiny::verbatimTextOutput("results")
      )
    ),
  shiny::hr(),
  shiny::div(style="text-align:center; color:gray", shiny::tags$small("(c) 2021 Stefan Mayer (University of Tübingen)"))
  ),

  server = function(input, output, session) {

    # display image
    output$img <- shiny::renderImage({list(src = img, height = 100)},
                                     deleteFile = FALSE)

    # data simulation

    sim <- shiny::eventReactive(input$simulate,{
      if(input$model == 'Simple'){
        set.seed(input$seed)
        N <- input$N # number of observations
        sigma_M <- input$sigma_M # error std M
        sigma_Y <- input$sigma_Y # error std Y
        beta_1 <- c(input$beta_0M, input$beta_1) # beta_0M and beta_1
        beta_2 <- c(input$beta_0Y, input$beta_2, input$beta_3) # beta_0Y, beta_2, beta_3
        X <- stats::rnorm(N, mean = 1,sd = 1) # generate random X
        # generate data based on parameters
        simPartialMed(beta_1, beta_2, sigma_M, sigma_Y, N, X)
      }
    })

    out <- shiny::eventReactive(input$estimatemodel,{
      shiny::req(sim())
      if(input$model == 'Simple'){
        # Estimation and Bayes factor computation
        #
        # Choosing the reference prior for the direct effect (beta_3)
        A_M <- c(input$prior_0M, input$prior_1); # Prior variance for beta_0M, beta_1
        A_Y_ref <- c(input$prior_0Y, input$prior_2, input$prior_3) # Prior variance for beta_0Y, beta_2, beta_3(carefully chosen reference prior, please see xxx & xxx (2020) before changing)
        prior <- list(A_M = A_M, A_Y = A_Y_ref)
        shiny::showModal(shiny::modalDialog("Estimating ...", footer=NULL))
        res <- Mediate(Data = sim(), Model = input$model, Prior = prior, R = input$R, burnin = input$burnin)
        shiny::removeModal()
        res
      }
    })

    output$summary_heading <- shiny::renderText({
      shiny::req(sim())
      paste0("<h4>Data summary</h4>")
    })

    output$summary <- shiny::renderPrint({
      if(input$model != 'Simple') return('So far only the simple mediation model is implemented.')
      shiny::req(sim())
      df <- data.frame(sim())
      summary(df)
    })

    output$df <- DT::renderDataTable({
      if(input$model != 'Simple') return(NULL)
      DT::datatable(data.frame(sim()), options = list(
        searching = FALSE, # no searching / filtering the data
        pageLength=5))     # show 5 entries as default
    })

    output$results_heading <- shiny::renderText({
      shiny::req(sim())
      shiny::req(out())
      paste0("<h4>Results</h4>")
    })

    output$results <- shiny::renderPrint({
      shiny::req(sim())
      shiny::req(out())
      if(input$model != 'Simple') return('So far only the simple mediation model is implemented.')
      cat(
        "Indirect effect 95% posterior credible interval:",
        out()$Simple$Indirect_CI,
        "\nDirect effect 95% posterior credible interval:",
        out()$Simple$Direct_CI,
        "\nBayes factor:",
        out()$Simple$BF,
        "\nVerdict:",
        out()$Simple$evidence
      )
    })
  }
)

# Run the application
# shinyApp(ui = ui, server = server)

