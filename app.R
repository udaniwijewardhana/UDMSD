# Loading libraries
library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(breakpoint)
library(cumSeg)
library(changepoint)
library(bcp)

################################################################################################################
# Shiny App for Single Species Univariate Changepoint Analysis
################################################################################################################

### Shiny user interface ###

ui <- fluidPage(
  
titlePanel(strong("UDMCA - Shiny App for Single Species Univariate Changepoint Analysis", style = "color:#3474A7")),
hr(),
  
div(style="display: inline-block;vertical-align:top; width: 300px;", fileInput("file", "Choose data CSV File", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
div(style="display: inline-block;vertical-align:top; width: 300px;", selectInput("datanorm", "Count data normalization:", choices=c("rnorm", "stand", "none"), selected = "none")),
div(style="display: inline-block;vertical-align:top; width: 300px;", selectInput("prednorm", "Predictors normalization:", choices=c("rnorm", "stand", "none"), selected = "none")),

tabsetPanel(
    
    tabPanel("Auunual Data",
             fluidRow(column(12, DT::dataTableOutput("contents")))
    ),
    
    tabPanel("Changepoint Analysis",
             sidebarLayout(
             sidebarPanel(div(style='height:900px; overflow: scroll',
             selectInput("changepoint", "changepoint method", choices=c("changepoint", "breakpoint", "cumSeg", "bcp"), selected = "changepoint"),
             
             # changepoint method                   
             conditionalPanel(
             condition = "input.changepoint == 'changepoint'",
             tags$h3("Algorithms of changepoint package:"),
             tags$ol(
             tags$li("cpt.mean(data,penalty='MBIC',pen.value=0,method='AMOC', Q=5,test.stat='Normal',class=TRUE,param.estimates=TRUE, minseglen=1)"),
             tags$li("cpt.var(data,penalty='MBIC',pen.value=0,know.mean=FALSE, mu=NA,method='AMOC',Q=5,test.stat='Normal', class=TRUE,param.estimates=TRUE,minseglen=2)"),
             tags$li("cpt.meanvar(data,penalty='MBIC',pen.value=0,method='AMOC', Q=5,test.stat='Normal',class=TRUE,param.estimates=TRUE, shape=1,minseglen=2)")),
             selectInput("changes", "type of changes to identify", choices=c("mean", "variance", "mean and variance"), selected = "mean"),
             selectInput("penalty", "penalty", choices=c("BIC", "MBIC", "AIC"), selected = "AIC"),
             selectInput("method", "method", choices=c("AMOC", "PELT", "SegNeigh", "BinSeg"), selected = "PELT"),
             numericInput("Q", "maximum number of changepoints to search", 5, min = 0, max = NA, step = 1),
             helpText("The maximum number of changepoints to search for using the BinSeg method.
             The maximum number of segments (number of changepoints + 1) to search for using the SegNeigh method."),
             selectInput("test.stat", "test statistic", choices=c("Normal", "Gamma", "Exponential", "Poisson"), selected = "Normal"),
             helpText("Gamma, Exponential and Poisson distribution only applicable for mean and variance function."),
             selectInput("class", "class", choices=c("TRUE", "FALSE"), selected = "TRUE"),
             selectInput("param.estimates", "param.estimates", choices=c("TRUE", "FALSE"), selected = "TRUE"),
             numericInput("shape", "shape parameter for Gamma distribution", 1, min = -100, max = NA, step = 0.001),
             numericInput("minseglen", "minseglen", 1, min = 1, max = NA, step = 1),
             helpText("Positive integer giving the minimum segment length (no. of observations between
                                           changes), default is the minimum allowed by theory. For cpt.var and cpt.meanvar,
                                           the minimum minseglen is 2 and for cpt.mean, the minimum minseglen is 1.")),
                                
             # breakpoint method
             conditionalPanel(
             condition = "input.changepoint == 'breakpoint'",
             tags$h3("Algorithms of breakpoint package:"),
             tags$ol(
             tags$li("CE.NB(data, Nmax = 10, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8, distyp = 1, penalty = 'BIC', parallel = FALSE)"),
             tags$li("CE.NB.Init(data, init.locs, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8, distyp = 1, penalty = 'BIC', var.init = 1e+05, parallel = FALSE)"),
             tags$li("CE.Normal.Init.Mean(data, init.locs, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8, distyp = 1, penalty = 'BIC', var.init = 1e+05, parallel = FALSE)"),
             tags$li("CE.Normal.Init.MeanVar(data, init.locs, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8, distyp = 1, penalty = 'BIC', var.init = 1e+05, parallel = FALSE)"),
             tags$li("CE.Normal.Mean(data, Nmax = 10, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8, distyp = 1, penalty = 'BIC', parallel = FALSE)"),
             tags$li("CE.Normal.MeanVar(data, Nmax = 10, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8, distyp = 1, penalty = 'BIC', parallel = FALSE)"),
             tags$li("CE.ZINB(data, Nmax = 10, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8, distyp = 1, penalty = 'BIC', parallel = FALSE)"),
             tags$li("CE.ZINB.Init(data, init.locs, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8, distyp = 1, penalty = 'BIC', var.init = 1e+05, parallel = FALSE)")),
             selectInput("algorithm", "breakpoint function", choices=c("CE.NB", "CE.NB.Init", "CE.Normal.Init.Mean", "CE.Normal.Init.MeanVar", "CE.ZINB", "CE.ZINB.Init"), selected = "CE.NB"),
             numericInput("Nmax", "maximum number of breakpoints to search", 10, min = 0, max = NA, step = 1),
             numericInput("eps", "the cut-off value for the stopping criterion in the CE method", 0.01, min = 0, max = 100, step = 0.00001),
             numericInput("rho", "the fraction which is used to obtain the best performing set of sample solutions", 0.05, min = 0, max = NA, step = 0.00001),
             numericInput("M", "sample size to be used in simulating the locations of break-points", 200, min = 0, max = NA, step = 1),
             numericInput("h", "minimum aberration width", 5, min = 0, max = NA, step = 1),
             numericInput("a", "a smoothing parameter value", 0.8, min = 0, max = NA, step = 0.00001),
             helpText("It is used in the four parameter beta distribution to smooth both shape parameters.
                      When simulating from the truncated normal distribution,
                      this value is used to smooth the estimates of the mean values."),
             numericInput("b", "a smoothing parameter value", 0.8, min = 0, max = NA, step = 0.00001),
             helpText("It is used in the truncated normal distribution to smooth the estimates of the standard deviation."),
             selectInput("distyp", "distribution to simulate break-point locations", choices=c("1", "2"), selected = "1"),
             helpText("Options: 1 = four parameter beta distribution, 2 = truncated normal distribution"),
             selectInput("penalty", "penalty", choices=c("BIC", "AIC"), selected = "BIC"),
             numericInput("var.init", "Initial variance value to facilitate the search process", 100000, min = 0, max = NA, step = 0.00001),
             selectInput("parallel", "parallel", choices=c("TRUE", "FALSE"), selected = "FALSE"),
             textInput("init.locs", "Initial break-point locations - enter a vector (comma delimited) - e.g. '0,1,2'", NULL)),
                                
             # cumSeg method                 
             conditionalPanel(
             condition = "input.changepoint == 'cumSeg'",
             tags$h3("Algorithms of cumSeg package:"),
             tags$li("fit.control(toll = 0.001, it.max = 5, display = FALSE, last = TRUE, maxit.glm = 25, h = 1, stop.if.error = FALSE)"),
             tags$li("sel.control(display = FALSE, type = c('bic', 'mdl', 'rss'), S = 1, Cn = log(log(n)), alg = c('stepwise', 'lasso'), edf.psi = TRUE)"),
             tags$li("jumpoints(y, x, k = min(30, round(length(y)/10)), output = '2', psi = NULL, round = TRUE, control = fit.control(), selection = sel.control())"),
             hr(),
             numericInput("k", "k = the starting number of changepoints", 2, min = 0, max = NA, step = 1),
             helpText("It should be quite larger than the supposed number of (true) changepoints.
                      This argument is ignored if starting values of the
                      changepoints are specified via psi."),
             selectInput("output", "output", choices=c(1,2,3), selected = 1),
             textInput("psi", "psi = starting values for the changepoints - enter a vector (comma delimited) - e.g. '0,1,2'", NULL),
             helpText("When psi=NULL (default), k quantiles are assumed."),
             selectInput("round", "round", choices=c("TRUE", "FALSE"), selected = "TRUE"),
             numericInput("toll", "toll = positive convergence tolerance", 0.001, min = 0, max = NA, step = 0.0001),
             numericInput("it.max", "it.max = integer giving the maximal number of iterations", 5, min = 0, max = NA, step = 1),
             selectInput("display", "display", choices=c("TRUE", "FALSE"), selected = "FALSE"),
             selectInput("last", "last", choices=c("TRUE", "FALSE"), selected = "TRUE"),
             numericInput("maxit.glm", "maxit.glm", 25, min = 0, max = NA, step = 1),
             numericInput("h", "h", 1, min = 0, max = NA, step = 1),
             selectInput("stop.if.error", "stop.if.error", choices=c("TRUE", "FALSE"), selected = "FALSE"),
             helpText("logical indicating if the algorithm should stop when one or more estimated
                      changepoints do not assume admissible values."),
             selectInput("type", "type", choices=c("bic", "mdl", "rss"), selected = "bic"),
             numericInput("S", "S", 1, min = 0, max = NA, step = 1),
             helpText("If type = 'rss' the optimal model is selected when the residual sum of squares
                      decreases by the threshold S."),
             numericInput("Cn", "Cn", 1, min = 0, max = NA, step = 0.0000001),
             helpText("If type= 'bic' a character string (as a function of 'n') to specify to generalized
                      BIC. If Cn=1 the standard BIC is used."),
             selectInput("alg", "alg", choices=c("stepwise", "lasso"), selected = "stepwise"),
             selectInput("edf.psi", "edf.psi", choices=c("TRUE", "FALSE"), selected = "TRUE")),
                                
             # bcp method
             conditionalPanel(
             condition = "input.changepoint == 'bcp'",
             tags$h3("Algorithm of bcp package:"),
             tags$li("bcp(y, x = NULL, id = NULL, adj = NULL, w0 = NULL, p0 = 0.2, d = 10, burnin = 50, mcmc = 500, return.mcmc = FALSE,
             boundaryType = 'node', p1 = 1, freqAPP = 20, nreg = -1)"),
             #numericInput("adjk", "k = the number of neighbors assumed for a typical vertex in adj list", 4, min = 4, max = 8, step = 4),
             #helpText("An adjacency list. Indexing the observations from 1 to n, the i-th
             #         element of the list is a vector of indices (offset by 1) for nodes that share an
             #         edge with node i. Generates an adjacency list for a n node by m node grid,
             #         assuming a maximum of k neighbors. k must be always 4 or 8. n and m are always >= 2.
             #         Adjacency list algorithm = makeAdjGrid(n,m,k)"),
             selectInput("pred", "predictors", choices=c("Yes", "No"), selected = "No"),
             helpText("w0 is applicable when predictors are available."),
             textInput("w0", "w0 - enter a vector (comma delimited) - e.g. '0.2,0.5,0.7'", NULL),
             numericInput("p0", "p0", 0.2, min = 0, max = 1, step = 0.01),
             helpText("w0 and p0 are optional values which are between 0 and 1
                      for Barry and Haritgan's hyperparameters; these
                      default to the value 0.2, which has been found to work well."),
             numericInput("p1", "p1 = The proportion of Active Pixel Passes run that are the actual
                          Active Pixel Passes", 1, min = 0, max = 1, step = 0.01),
             numericInput("d", "d", 10, min = 0, max = NA, step = 1),
             helpText("a positive number only used for linear regression change point models.
                      Lower d means higher chance of fitting the full linear model
                      (instead of the intercept-only model)."),
             numericInput("burnin", "burnin", 50, min = 0, max = NA, step = 1),
             numericInput("mcmc", "mcmc", 500, min = 0, max = NA, step = 1),
             selectInput("return.mcmc", "return.mcmc", choices=c("TRUE", "FALSE"), selected = "FALSE"),
             numericInput("freqAPP", "freqAPP", 20, min = 0, max = NA, step = 1),
             selectInput("boundaryType", "boundaryType", choices=c("node", "edge"), selected = "node"),
             numericInput("nreg", "nreg", -1, min = 2, max = NA, step = 1),
             helpText("only applicable for regression; related to parameter d describing the
                      minimum number of observations needed in a block to allow for fitting a regression
                      model. Defaults to 2*number of predictors.")
))),
               
mainPanel(
          tags$h3("Summary results of changepoint method:"),
          fluidRow(div(style='height:400px; overflow: scroll', verbatimTextOutput("changepoint"))),
          tags$h3("Changepoint Plot:"),
          fluidRow(plotOutput("changepointPlot")))
))))

### Shiny server ###

server <- function(input, output) {
  
# Input data csv file
  
filedata <- reactive({
    inFile <- input$file
    if (is.null(inFile)){return(NULL)}
    
    x <- data.frame(read.csv(inFile$datapath, fileEncoding="UTF-8-BOM"))

    x$Count0 <- x$Count
    x$Count0[is.na(x$Count0)] <- 0
    
    if(ncol(x)>3){
    p = subset(x, select = -c(Year,Count,Count0))
    }else {p = NULL}
    
    if(!is.null(p)){
    for(i in 1:ncol(p)){
    if(input$prednorm == "rnorm"){p[,i] <- round(rnorm(p[,i]), digits = 4)
    }else if(input$prednorm == "stand") {p[,i] <- round(scale(p[,i]), digits = 4)
    }else {p[,i] <- round(p[,i], digits = 4)}}
    }
    
    if(input$datanorm == "rnorm"){x$Count0 <- round(rnorm(x$Count0), digits = 4)
    } else if (input$datanorm == "stand"){x$Count0 <- round(scale(x$Count0), digits = 4) 
    } else {x$Count0 <- round(x$Count0, digits = 4)}
    
    z = subset(x, select = c(Year,Count,Count0))
    
    if(!is.null(p)){
      Final = cbind(z, p)
    }else {
      Final = x
    }
    return(Final)
})

predictors <- reactive({
  x <- filedata()
  p <- subset(x, select = -c(Year,Count,Count0))
  if(is.null(p)){pred = NULL}else{pred = p}
  return(p)
})
  
fit <- reactive({
       fit <- fit.control(toll = input$toll, it.max = input$it.max, display = input$display, last = input$last,
       maxit.glm = input$maxit.glm, h = input$h, stop.if.error = input$stop.if.error)
return(fit)
})
  
sel <- reactive({
       sel <- sel.control(display = input$display, type = input$type, S = input$S,
       Cn = input$Cn, alg = input$alg, edf.psi = input$edf.psi)
return(sel)
})
  
psi <- reactive({
    
    data <- filedata()
    x <- 1:nrow(data)
    
    if(input$psi == ""){
      psi = quantile(x, prob= seq(0,1,l=input$k+2)[-c(1,input$k+2)], names=FALSE)
      return(psi)
    } else {
      psi = as.numeric(unlist(strsplit(input$psi,",")))
      return(psi)
    }
})
  
w0 <- reactive({
    
    inFile <- input$file
    if (is.null(inFile)){return(NULL)}
    
    x <- filedata()
    predictors <- predictors()
    
    w = as.numeric(unlist(strsplit(input$w0,",")))
    if (any(w > 1) | any(w < 0)){stop("Each element in w0 must be between 0 and 1.")}
    
    if(input$pred == "No"){return(NULL)}
    if(input$pred == "Yes"){
    if(length(w) == 0){
        w <- NULL
    } else {
    if(length(w) == ncol(predictors)){
        w <- w
    } else if(length(w) > ncol(predictors)){
        stop("Length of w0 is greater than the number of predictors.")
    } else if(length(w) == 1){
        w <- rep(w, ncol(predictors))
        print("Assume you wanted each error-to-signal ratio to be iid from U(0, w0).")
    } else {
        stop("Length of w0 is less than the number of predictors.")
    }
    }}
    return(w)
})
  
adjx <- reactive({
  
  x <- filedata()
  predictors <- predictors()
  n = nrow(x)
  
  if (is.null(ncol(x)>3)){
    m = nrow(x)
    adj =  makeAdjGrid(n, m, input$adjk)
  }else {
    m <- ncol(predictors)
    adj = makeAdjGrid(n, m, input$adjk)
  }
  return(adj)
})
  
# Changepoint analysis
  
changepoint <- reactive({
    
    data <- filedata()
    predictors = predictors()
    w0 <- w0()
    #adj <- adjx()
    
    if(input$changepoint == "changepoint"){
      
    if(input$method == "Gamma"){
        
    if(input$changes == "mean"){
          cpt = cpt.mean(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                         test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                         shape=input$shape, minseglen=input$minseglen)
    } else if(input$changes == "variance"){
          cpt = cpt.var(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                        test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                        shape=input$shape, minseglen=input$minseglen)
    } else{
          cpt = cpt.meanvar(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                            test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                            shape=input$shape, minseglen=input$minseglen)
    }} else{
          
    if(input$changes == "mean"){
          cpt = cpt.mean(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                         test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                         minseglen=input$minseglen)
    } else if(input$changes == "variance"){
          cpt = cpt.var(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                        test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                        minseglen=input$minseglen)
    } else{
          cpt = cpt.meanvar(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                            test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                            minseglen=input$minseglen)
          }}
    return(cpt)
      
    } else if(input$changepoint == "breakpoint"){
      
    if(input$algorithm == "CE.NB"){
        bp = CE.NB(as.data.frame(data$Count0), Nmax = input$Nmax, eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                   a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty, parallel = input$parallel)
    } else if(input$algorithm == "CE.NB.Init"){
        bp = CE.NB.Init(as.data.frame(data$Count0), as.numeric(unlist(strsplit(input$init.locs,","))),
                        eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                        a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty,
                        var.init = input$var.init, parallel = input$parallel)
    } else if(input$algorithm == "CE.NB.Init.Mean"){
        bp = CE.Normal.Init.Mean(as.data.frame(data$Count0), as.numeric(unlist(strsplit(input$init.locs,","))),
                                 eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                                 a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty,
                                 var.init = input$var.init, parallel = input$parallel)
    } else if(input$algorithm == "CE.NB.Init.MeanVar"){
        bp = CE.Normal.Init.MeanVar(as.data.frame(data$Count0), as.numeric(unlist(strsplit(input$init.locs,","))),
                                    eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                                    a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty,
                                    var.init = input$var.init, parallel = input$parallel)
    } else if(input$algorithm == "CE.Normal.Mean"){
        bp = CE.Normal.Mean(as.data.frame(data$Count0), Nmax = input$Nmax, eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                            a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty, parallel = input$parallel)
    } else if(input$algorithm == "CE.Normal.MeanVar"){
        bp = CE.Normal.MeanVar(as.data.frame(data$Count0), Nmax = input$Nmax, eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                               a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty, parallel = input$parallel)
    } else if(input$algorithm == "CE.ZINB"){
        bp = CE.ZINB(as.data.frame(data$Count0), Nmax = input$Nmax, eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                     a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty, parallel = input$parallel)
    } else {
        bp = CE.ZINB.Init(as.data.frame(data$Count0), as.numeric(unlist(strsplit(input$init.locs,","))),
                          eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                          a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty,
                          var.init = input$var.init, parallel = input$parallel)
    }
    return(bp)
      
    } else if(input$changepoint == "cumSeg"){
      
    x <- 1:nrow(data)
      
    cumSeg = jumpoints(as.matrix(data$Count0), x, output = input$output,
                       k = input$k, round = input$round, psi = psi(), control = fit(), selection = sel())
    return(cumSeg)
    } else {
      
    if(input$pred == "Yes" && !is.null(predictors)){pred <- as.matrix(predictors)
    } else {pred <- NULL}
      
    bcp = bcp::bcp(y = as.matrix(data$Count0), x = pred, id = NULL,
                   adj = NULL, w0 = w0(),
                   p0 = input$p0, d = input$d, burnin = input$burnin, mcmc = input$mcmc,
                   return.mcmc = input$return.mcmc, boundaryType = input$boundaryType, p1 = input$p1,
                   freqAPP = input$freqApp, nreg = input$nreg)
    return(bcp)
    }
})
  
# Output of data table
  
output$contents <- DT::renderDataTable({
    if (is.null(input$file)) {return(NULL)}
    df <- DT::datatable(filedata())
    df
})
  
# Output of changepoint method summary results
  output$changepoint <- renderPrint({
  req(input$file)
  changepoint()
})
  
# Output of changepoint method plot
  
output$changepointPlot <- renderPlot({
  req(input$file)
  if(input$changepoint == "breakpoint"){
  data <- filedata()
  profilePlot(changepoint(), as.data.frame(data$Count0))
  } else {
  plot(changepoint())}
})
  
}

shinyApp(ui, server)
