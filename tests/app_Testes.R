library(shiny)
library(plotly)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(Metrics)
library(RcppFaddeeva)
library(rlang)

#Define function
f1 <- function(colnames,rownumber) {
  tbl_colnames <- colnames
  tbl_colnames %>% purrr::map_dfc(~tibble::tibble(!!.x := numeric())) %>%
    bind_rows(as_tibble(matrix(data = rep(0,length(colnames)*rownumber),
                               nrow = rownumber,
                               ncol = length(colnames),
                               dimnames = list(NULL,tbl_colnames))))
}
finite.differences <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }

  n <- length(x)

  # Initialize a vector of length n to enter the derivative approximations
  fdx <- vector(length = n)

  # Iterate through the values using the forward differencing method
  for (i in 2:n) {
    fdx[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i])
  }

  # For the last value, since we are unable to perform the forward differencing method
  # as only the first n values are known, we use the backward differencing approach
  # instead. Note this will essentially give the same value as the last iteration
  # in the forward differencing method, but it is used as an approximation as we
  # don't have any more information
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])

  return(fdx)
}
norm01 <- function(x){
  stp1 <- x - min(x)
  stp1/max(stp1)

}
# Define UI for application
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "darkly"),
  # Application title
  titlePanel("Preview Deconvolutions"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      tags$hr(),
      sliderInput("alpha","Experimental transparency",
                  min = 0.05,
                  max = 1,
                  value = 1,
                  step = 0.05),
      numericInput(inputId = "npeaks",
                   label = "Number of Peaks",
                   value = 1,min = 0,max = 15),
      uiOutput("secondSelection"),
      uiOutput("x0"),
      numericInput("A",
                   label = "Amplitude",
                   value = 1,min = 0,max = 100,step = .5),
      numericInput("sigma",
                   label = "Sigma",
                   value = 1,min = 0,max = 100,step = .25),
      numericInput("gamma",
                   label = "Gamma",
                   value = 1,min = 0,max = 100,step = .25),
      actionButton(inputId = "add",label = "Add band",
                   icon = icon("plus-square"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotlyOutput("plot")),
        tabPanel("Derivatives", plotlyOutput("derivatives")),
        tabPanel("Peaks stats", tableOutput("Peaks_tab"))
      ),
      downloadButton('download',"Download the peaks data")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$secondSelection <- renderUI({
    if(is.null(data())) return(NULL)
    selectInput(label = "Curve",
                paste("band",1:input$npeaks,sep = ""),
                inputId = "curve")
  })

  #Start x0 in the minimium val
  output$x0 <- renderUI({
    req(data())
    data_fix <- data()
    x_min <- min(data_fix[,1])
    x_max <- max(data_fix[,1])

    numericInput("x0",
                 label = "Center",
                 value = x_min,
                 min = x_min,
                 max = x_max)
  })

  #Reactive file
  inFile <- reactive({
    input$file1
  })

  #Reactive experimental data
  data <- reactive({
    req(input$file1)
    read.csv(inFile()$datapath,
             header = T,
             sep = ";")
  })

  #Define base tibble
  tbl <- reactive({
    req(data(), input$npeaks)
    cols <- c("Wavenumber",paste("band",1:input$npeaks,sep = ""))
    f1(colnames = cols,rownumber = nrow(data())) %>%
      mutate("Wavenumber" = data()[,1])
  })

  #### - Peaks update - ################################################
  dt <- tibble("Peak" = character(),
               "Center" = numeric(),
               "Amplitude" = numeric(),
               "Sigma" = numeric(),
               "Gamma" = numeric())

  DF1 <- reactiveValues(data=dt,
                        peaks = NULL) #Empty tibble of peak stats

  observeEvent(input$add, {
    req(input$curve, input$x0, input$A, input$sigma, input$gamma)

    dt <- DF1$data

    if (input$curve %in% dt$Peak) {
      dt <- dt %>%
        as_tibble() %>%
        rows_update(tibble("Peak" = (input$curve),
                           "Center" = input$x0,
                           "Amplitude" = input$A,
                           "Sigma" = input$sigma,
                           "Gamma" = input$gamma)) %>%
        arrange(Peak)
    } else {
      dt <- dt %>%
        as_tibble() %>%
        add_row("Peak" = (input$curve),
                "Center" = input$x0,
                "Amplitude" = input$A,
                "Sigma" = input$sigma,
                "Gamma" = input$gamma) %>%
        arrange(Peak)
    }

    DF1$data <- dt

    # Update peaks
    tmp1 <- if(is.null(DF1$peaks)) tbl() else DF1$peaks
    cv <- paste0(input$curve)

    DF1$peaks <- tmp1 %>%
      mutate(!!cv := input$A*Voigt(x = data()[,1],
                                   x0 = input$x0,
                                   sigma = input$sigma,
                                   gamma = input$gamma))
  })

  output$plot <- renderPlotly({
    req(data())

    data_fix <- data()
    x_min <- min(data_fix[,1])
    x_max <- max(data_fix[,1])
    by = (x_max-x_min)/10

    # Initialize fit1 based on whether peaks exist
    if(input$npeaks == 0 || is.null(DF1$peaks)) {
      fit1 <- rep(0, nrow(data_fix))
    } else {
      fit1 <- DF1$peaks[,-1] %>% rowSums()
    }

    # Create basic plot
    p <- ggplot(data = tibble(Wavenumber = data_fix[,1], y = data_fix[,2]),
                aes(x = Wavenumber, y = y)) +
      geom_point(color = "red", alpha = input$alpha) +
      theme_bw() +
      scale_x_reverse(breaks = seq(from = x_min, to = x_max, by = by)) +
      xlab("Wavenumber cm-1") +
      ylab("Normalized absorbance")

    # Add components if they exist
    if(input$npeaks > 0 && !is.null(DF1$peaks)) {
      MAE_value <- mae(actual = data_fix[,2], predicted = fit1)
      Prev <- input$A*Voigt(x = data_fix[,1],
                            x0 = input$x0,
                            sigma = input$sigma,
                            gamma = input$gamma)

      p <- p +
        geom_line(aes(y = fit1), color = "black", size = 1) +
        geom_area(data = pivot_longer(DF1$peaks,
                                      cols = starts_with("band"),
                                      names_to = "Band",
                                      values_to = "Values"),
                  aes(x = Wavenumber, y = Values, fill = Band),
                  alpha = .5, position = "identity") +
        viridis::scale_fill_viridis(discrete = TRUE) +
        geom_area(data = tibble(x = data_fix[,1], y = Prev),
                  aes(x = x, y = y),
                  fill = "navy", alpha = .5) +
        geom_text(data = tibble(x = data_fix[,1], y = Prev),
                  aes(x = input$x0, y = max(y), label = input$curve)) +
        annotate(geom = "text",
                 x = x_min + 10,
                 y = .05,
                 label = paste("MAE =", round(MAE_value, 5), sep = ""))
    }

    ggplotly(p)
  })

  output$derivatives <- renderPlotly({
    req(data())
    df1 <- data()
    colnames(df1) <- c("Wavenumber","Experimental")

    x_min <- min(df1[,1])
    x_max <- max(df1[,1])
    by = (x_max-x_min)/10

    df1 %>%
      mutate(First = finite.differences(x = Wavenumber, y = Experimental)) %>%
      mutate(Second = finite.differences(x = Wavenumber, y = First)) %>%
      select(-First) %>%
      mutate(Experimental = norm01(Experimental),
             Second = norm01(Second)) %>%
      pivot_longer(cols = -Wavenumber,
                   names_to = "Spec",
                   values_to = "Values") %>%
      ggplot(aes(x = Wavenumber, y = Values, col = Spec)) +
      geom_line() +
      facet_wrap(~Spec, nrow = 2, scales = "free") +
      theme_bw() +
      scale_x_reverse(breaks = seq(from = x_min,
                                   to = x_max,
                                   by = by)) +
      xlab("Wavenumber cm-1") +
      ylab("Normalized absorbance") +
      scale_color_manual(values = c("red","navy"))
  })

  output$Peaks_tab <- renderTable({
    req(DF1$data)
    DF1$data
  })

  output$download <- downloadHandler(
    filename = function() {
      paste('PeakStats', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(DF1$data, con, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
