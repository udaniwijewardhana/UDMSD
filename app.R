library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(lattice)
library(gridExtra)
library(leaflet)
library(INLA)

### Shiny user interface ###

ui <- fluidPage(
  
titlePanel(strong("UDMSD - Single species annual spatial and temporal joint model to geostatistical zero-inflated data", titleWidth = 350)),
hr(),
  
div(style="display: inline-block;vertical-align:top; width: 300px;", fileInput("file", "Choose data CSV File", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
div(style="display: inline-block;vertical-align:top; width: 300px;", selectInput("prednorm", "Predictors normalization:", choices=c("rnorm", "stand", "none"), selected = "none")),
  
tabsetPanel(
    
tabPanel("Auunual Data",
         fluidRow(style = "margin-top: 25px;",
         column(7, p(tags$b('Annual Count Data', style = "font-size: 150%; font-family:Helvetica; color:#4c4c4c; text-align:left;"))),
         column(5, p(tags$b('Annual Predictor Data', style = "font-size: 150%; font-family:Helvetica; color:#4c4c4c; text-align:left;")))),
         fluidRow(
         column(7, DT::dataTableOutput("contents")),
         column(5, verbatimTextOutput("predictors")))
),
    
tabPanel("Single-species Joint Model",
         sidebarLayout(
         sidebarPanel(div(style='height:840px; overflow: scroll',
         sliderInput(inputId = "offset", label = "offsets for automatic boundaries", value = c(0.1, 0.3), min = 0.05, max = 2),
         sliderInput(inputId = "max.edge", label = "max.edge", value = c(0.05, 0.5), min = 0.01, max = 2),
         sliderInput(inputId = "cutoff", label = "cutoff", value = 0.01, min = 0, max = 0.2),
         sliderInput(inputId = "min.angle", label = "min.angle", value = c(21, 30), min = 1, max = 35),
         numericInput("convex", "convex for boundary", 0.5, min = NA, max = NA, step = 0.0001),
         actionButton("mesh", "mesh"),
         selectInput("distribution", "Distribution of count model:",
                     choices=c("Negative Binomial", "Zeroinflated Negative Binomial",
                               "Negative Binomial Hurdle"), selected = "Negative Binomial Hurdle"),
         helpText("Posterior plots are only applicable for spatial or spatio-temporal models."),
         selectInput("speffect", "spatial random effect model:",
                     choices=c("spde", "iid"), selected = "iid"),
         selectInput("tempeffect", "temporal random effect model:",
                     choices=c("ar1", "iid", "rw1", "rw2"), selected = "ar1"),
         selectInput("int", "interaction between predictor variables:", choices = c("0","1"), selected = "0"),
         selectInput("term", "variables which have interactions:", choices=c("1 and 2","1 and 3", "1 and 4", "1 and 5",
                     "2 and 3", "2 and 4", "2 and 5", "3 and 4", "3 and 5", "4 and 5", "all"), selected = ""),
         actionButton("summary", "Summary")
)),
         
mainPanel(fluidRow(column(4, tags$h3("Mesh:"), uiOutput("tab")),
                   column(8, tags$h3("Posterior Plots:"))),
          fluidRow(column(4, div(style='height:360px', plotOutput("mesh", "mesh", width = "80%", height = "360px"))),
                   column(4, div(style='height:360px', plotOutput("posteriormPlot"))),
                   column(4, div(style='height:360px', plotOutput("posteriorsdPlot")))),
          fluidRow(column(12, tags$h3("Summary results of species joint model:"), div(style='height:400px; overflow: scroll',
                   verbatimTextOutput("summary")))))
))))

### Shiny server ###

server <- function(input, output) {
  
# Input data csv file
  
filedata1 <- reactive({
    inFile <- input$file
    if (is.null(inFile)){return(NULL)}
    
    x <- as.data.frame(read.csv(inFile$datapath, fileEncoding="UTF-8-BOM"))
    
    y = dplyr::select_if(x, is.numeric)
    if(ncol(x)>5){p = subset(y, select = -c(Year, Count, Latitude, Longitude))
    }else {p = NULL}
    
    if(!is.null(p)){
    for(i in 1:ncol(p)){
    if(input$prednorm == "rnorm"){p[,i] <- round(rnorm(p[,i]), digits = 4)
    }else if(input$prednorm == "stand") {p[,i] <- round(scale(p[,i]), digits = 4)
    }else {p[,i] <- round(p[,i], digits = 4)}}}
    
    if(is.null(p)){
    Final = x[ , (names(x) %in% c("Year", "Count", "Locality", "Latitude", "Longitude"))]
    }else {
    Final = cbind(x[ , (names(x) %in% c("Year", "Count", "Locality", "Latitude", "Longitude"))], p)
    }
    return(Final)
})
  
# Subset possible numeric predictor variables

filedata2 <- reactive({
  req(input$file)
  x <- filedata1()
  
  y = dplyr::select_if(x, is.numeric)
  if(ncol(x)>5){
  p = subset(y, select = -c(Year, Count, Latitude, Longitude))
  }else {p = NULL}
  
  if(!is.null(p)){
  for(i in 1:ncol(p)){
  if(input$prednorm == "rnorm"){p[,i] <- round(rnorm(p[,i]), digits = 4)
  }else if(input$prednorm == "stand"){p[,i] <- round(scale(p[,i]), digits = 4)
  }else {p[,i] <- round(p[,i], digits = 4)}}
  }
  return(p)
})

# Output of the data table

output$contents <- DT::renderDataTable({
  req(input$file)
  df <- filedata1()
  return(DT::datatable(df, options = list(scrollX = TRUE)))
})

# Output of the numeric predictors summary table

output$predictors <- renderPrint({
  req(input$file)
  df <- filedata2()
  if (is.null(df)){return(NULL)}
  return(summary(df))
})
  
# Create INLA mesh
  
mesh <- reactive({
    m <- filedata1()
    coords <- cbind(m$Latitude, m$Longitude)
    bnd = inla.nonconvex.hull(coords, convex = input$convex)
    
    out <- INLA::inla.mesh.2d(
      loc = coords,
      boundary = bnd,
      max.edge = input$max.edge,
      min.angle = rev(input$min.angle),
      cutoff = input$cutoff,
      offset = input$offset)
    return(out)
})
  
mm <- eventReactive(input$mesh, {
  m <- mesh()
  df <- filedata1()
  coords <- cbind(df$Latitude, df$Longitude)
  plot(m)
  points(coords, col = "red")
})
  
# Output of mesh
  
output$mesh <- renderPlot({
  return(mm())
})
  
filedata3 <- reactive({
  if (is.null(input$file)){return(NULL)}
  df <- filedata1()
  pred <- filedata2()
  n = nrow(pred)
  
  pred.z = lapply(seq_along(1:ncol(pred)), function(x) c(pred[,x], rep(NA, n)))
  pred.y = lapply(seq_along(1:ncol(pred)), function(x) c(rep(NA, n), pred[,x]))
  pred.z <- data.frame(matrix(unlist(pred.z), nrow=n*2),stringsAsFactors=FALSE)
  pred.y <- data.frame(matrix(unlist(pred.y), nrow=n*2),stringsAsFactors=FALSE)

  colnames(pred.z) = gsub(" ", "", paste("p", ".z", 1:ncol(pred.z)))
  colnames(pred.y) = gsub(" ", "", paste("p", ".y", 1:ncol(pred.y)))
  return(cbind(pred.z, pred.y))
})
  
filedata4 <- reactive({
  
    if (is.null(filedata2())){return(NULL)}
    df <- filedata1()
    count = df$Count
    n.y = n.z = n = nrow(df)
    z <- rep(0,n)
    y = count
    
    if(input$distribution == "Negative Binomial Hurdle"){
      y[y == 0] <- NA
    }else {y = count}
    
    if(input$distribution == "Negative Binomial Hurdle"){
      z[which(count > 0)]<-1 } else {z[which(count > 0 | count == 0)]<-1}
    
    z = c(z, rep(NA,n))
    y = c(rep(NA,n), y)
    Y = data.frame(cbind(z,y))
    return(Y)
})
  
# Create joint zero inflation model
  
fitsummary <- reactive({
  
  df1 <- filedata1()
  df2 <- filedata3()
  Y <- filedata4()
  Year = df1$Year
  coords <- cbind(df1$Latitude, df1$Longitude)
  n = nrow(df1)
  
  mu.z = rep(1:0, c(n, n))
  mu.y = rep(0:1, c(n, n))
  
  coords$hhid = mesh()$idx$loc
  S = rep(coords$hhid, 2)
  year = rep(Year, 2)
  
  spde = INLA::inla.spde2.matern(mesh(), constr = TRUE)
  
  sp <- if(input$speffect == "spde"){spde}else{"iid"}
  tempeffect <- as.character(if(input$tempeffect == "ar1"){"ar1"
  }else if(input$tempeffect == "rw1"){"rw1"
  }else if(input$tempeffect == "rw1"){"rw2"
  }else{"iid"})
       
  distribution <- as.character(if(input$distribution == "Negative Binomial"){"nbinomial"
  } else if(input$distribution == "Zeroinflated Negative Binomial") {"zeroinflatednbinomial1"
  } else {"zeroinflatednbinomial0"})
  
if (is.null(filedata2())){

  data = list(Y = Y, mu.z = mu.z, mu.y = mu.y)
  formula = Y ~ 0 + mu.z + mu.y + f(S, model = sp) + f(year, model = tempeffect)
      
  model <- INLA::inla(formula, data = data, family = c("binomial", distribution),
           control.family = list(list(link = "logit"), list(link = "log")),
           control.compute = list(dic = TRUE,cpo = TRUE, po = TRUE), verbose = TRUE)
  return(model)
      
} else {
      
      pred.z <- df2[,c(1:(ncol(df2)/2))]
      pred.y <- df2[,c(((ncol(df2)/2)+1):ncol(df2))]
      p.z = matrix(NA, ncol = 5, nrow = n*2)
      p.y = matrix(NA, ncol = 5, nrow = n*2)
      for( i in 1:ncol(pred.z)){p.z[,i] = pred.z[,i]}
      for( i in 1:ncol(pred.y)){p.y[,i] = pred.y[,i]}
      
      p.z1 = p.z[,1]; p.z2 = p.z[,2]; p.z3 = p.z[,3]; p.z4 = p.z[,4]; p.z5 = p.z[,5]
      p.y1 = p.y[,1]; p.y2 = p.y[,2]; p.y3 = p.y[,3]; p.y4 = p.y[,4]; p.y5 = p.y[,5]
      
      npred = (ncol(df2)/2)
      
      int <- as.integer(if(input$int == "1"){1}else {0})
      
      if(npred == 2 & int == 1){if(input$term == "1 and 2"){term <- as.integer(1)}else{return("error")}
      }else if(npred == 3 & int == 1){
      if(input$term == "1 and 2"){term <- as.integer(1)
        }else if(input$term == "1 and 3"){term <- as.integer(2)
        }else if(input$term == "2 and 3"){term <- as.integer(5)
        }else if(input$term == "all"){term <- as.character("all")
        }else{return("error")}
      }else if(npred == 4 & int == 1){
        if(input$term == "1 and 2"){term <- as.integer(1)
        }else if(input$term == "1 and 3"){term <- as.integer(2)
        }else if(input$term == "1 and 4"){term <- as.integer(3)
        }else if(input$term == "2 and 3"){term <- as.integer(5)
        }else if(input$term == "2 and 4"){term <- as.integer(6)
        }else if(input$term == "3 and 4"){term <- as.integer(8)
        }else if(input$term == "all"){term <- as.character("all")
        }else{return("error")}
      }else if(npred == 5 & int == 1){
        if(input$term == "1 and 2"){term <- as.integer(1)
        }else if(input$term == "1 and 3"){term <- as.integer(2)
        }else if(input$term == "1 and 4"){term <- as.integer(3)
        }else if(input$term == "1 and 5"){term <- as.integer(4)
        }else if(input$term == "2 and 3"){term <- as.integer(5)
        }else if(input$term == "2 and 4"){term <- as.integer(6)
        }else if(input$term == "2 and 5"){term <- as.integer(7)
        }else if(input$term == "3 and 4"){term <- as.integer(8)
        }else if(input$term == "3 and 5"){term <- as.integer(9)
        }else if(input$term == "4 and 5"){term <- as.integer(10)
        }else if(input$term == "all"){term <- as.character("all")
        }else{return("error")}
      }
      
      data = list(Y = Y, mu.z = mu.z, mu.y = mu.y,
                  p.z1 = p.z1, p.y1 = p.y1, p.z2 = p.z2, p.y2 = p.y2,
                  p.z3 = p.z3, p.y3 = p.y3, p.z4 = p.z4, p.y4 = p.y4,
                  p.z5 = p.z5, p.y5 = p.y5)
      
      if(npred == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.y1 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 2 & int == 0){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.y1 + p.y2 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 2 & int == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 + p.y1 * p.y2 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 3 & int == 0){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 + p.y1 + p.y2 + p.y3 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 3 & int == 1 & term == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 + p.z3 + p.y1 * p.y2 + p.y3 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 3 & int == 1 & term == 2){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z3 + p.z2 + p.y1 * p.y3 + p.y2 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 3 & int == 1 & term == 5){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z3 + p.y1 + p.y2 * p.y3 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 3 & int == 1 & term == "all"){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 * p.z3 + p.y1 * p.y2 * p.y3 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 4 & int == 0){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 + p.z4 + p.y1 + p.y2 + p.y3 + p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 4 & int == 1 & term == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 + p.z3 + p.z4 + p.y1 * p.y2 + p.y3 + p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 4 & int == 0 & term == 2){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z3 + p.z2 + p.z4 + p.y1 * p.y3 + p.y2 + p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 4 & int == 0 & term == 3){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z4 + p.z2 + p.z3 + p.y1 * p.y4 + p.y2 + p.y3 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 4 & int == 0 & term == 5){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z3 + p.z4 + p.y1 + p.y2 * p.y3 + p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 4 & int == 0 & term == 6){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z4 + p.z3 + p.y1 + p.y2 * p.y4 + p.y3 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 4 & int == 0 & term == 8){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 * p.z4 + p.y1 + p.y2 + p.y3 * p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 4 & int == 0 & term == "all"){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 * p.z3 * p.z4 + p.y1 * p.y2 * p.y3 * p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 + p.z4 + p.z5 + p.y1 + p.y2 + p.y3 + p.y4 + p.y5 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 1 & term == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 + p.z3 + p.z4 + p.z5 + p.y1 * p.y2 + p.y3 + p.y4 + p.y5 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 2){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z3 + p.z2 + p.z4 + p.z5 + p.y1 * p.y3 + p.y2 + p.y4 + p.y5 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 3){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z4 + p.z2 + p.z3 + p.z5 + p.y1 * p.y4 + p.y2 + p.y3 + p.y5 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 4){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z5 + p.z2 + p.z3 + p.z4 + p.y1 * p.y5 + p.y2 + p.y3 + p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 5){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z3 + p.z4 + p.z5 + p.y1 + p.y2 * p.y3 + p.y4 + p.y5 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 6){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z4 + p.z3 + p.z5 + p.y1 + p.y2 * p.y4 + p.y3 + p.y5 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 7){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z5 + p.z3 + p.z4 + p.y1 + p.y2 * p.y5 + p.y3 + p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 8){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 * p.z4 + p.z5 + p.y1 + p.y2 + p.y3 * p.y4 + p.y5 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 9){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 * p.z5 + p.z4 + p.y1 + p.y2 + p.y3 * p.y5 + p.y4 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 0 & term == 10){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 + p.z4 * p.z5 + p.y1 + p.y2 + p.y3 + p.y4 * p.y5 + f(S, model = sp) + f(year, model = tempeffect)
      } else if(npred == 5 & int == 1 & term == "all"){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 * p.z3 * p.z4 * p.z5 + p.y1 * p.y2 * p.y3 * p.y4 * p.z5 + f(S, model = sp) + f(year, model = tempeffect)
      } else {formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 * p.z3 * p.z4 * p.z5 + p.y1 * p.y2 * p.y3 * p.y4 * p.z5 + f(S, model = sp) + f(year, model = tempeffect)
      }

      model <- INLA::inla(formula, data = data, family = c("binomial", distribution),
                          control.family = list(list(link = "logit"), list(link = "log")),
                          control.compute = list(dic = TRUE,cpo = TRUE, po = TRUE), verbose = TRUE)
      return(model)
}
})
  
proj <- reactive({
  if (is.null(fitsummary())){(NULL)}
  proj <- INLA::inla.mesh.projector(mesh(), projection = "longlat", dims = c(150, 100))
  return(proj)
})
  
# Summary output of joint zero inflation model
  
fitsum <- eventReactive(input$summary, {
  fitsummary()
})
  
output$summary <- renderPrint({
  return(fitsum()$summary.fixed)
})
  
# Create posterior mean plot
  
pmPlot <- reactive({
    if(input$speffect == "iid"){return(NULL)
    } else {
    mp <- levelplot(row.values=proj()$x, column.values=proj()$y,
                      inla.mesh.project(proj(), fitsummary()$summary.random$S$mean),
                      xla='Latitude', yla='Longitude',
                      main='posterior mean plot', contour=TRUE,
                      xlim=range(proj()$x), ylim=range(proj()$y))
    return(mp)}
})
  
# Create posterior standard deviation plot
  
psdPlot <- reactive({
    if(input$speffect == "iid"){return(NULL)
    } else {
    sdp <- levelplot(row.values=proj()$x, column.values=proj()$y,
                       inla.mesh.project(proj(), fitsummary()$summary.random$S$sd),
                       xla='Latitude', yla='Longitude',
                       main='posterior standard deviation plot', contour=TRUE,
                       xlim=range(proj()$x), ylim=range(proj()$y))
    return(sdp)}
})
  
# Output of posterior standard deviation plot
  
output$posteriormPlot <- renderPlot({pmPlot()})
  
# Output of posterior standard deviation plot
  
output$posteriorsdPlot <- renderPlot({psdPlot()})
  
url <- a("Definition", href="https://rdrr.io/github/andrewzm/INLA/man/inla.mesh.2d.html")
  output$tab <- renderUI({
  tagList("URL link:", url)
})
  
}

shinyApp(ui, server)
