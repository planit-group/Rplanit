# APPS -------------------------------------------------------------------------
# (GUI ed altro)

#' Values eploration (Web App)
#' 
#' Start an interactive GUI (web app) to explore the evaluated values ditribution.
#' It can use directly  a plan object as argument (deriving from it the corresponding values, ct and contours objects to visualize).
#' It can accept also a list of
#' plan objects, to explore simultaneusly different plans.
#' @param plan The plan object
#' @param values The values object (optional if plan is given)
#' @param ct The CT object (optional)
#' @param contours The contour object (optional)
#' @family Apps
#' @export
#' @import shiny
values.app <- function(plan=NULL, values=NULL, ct=NULL, contours=NULL, sanitize=TRUE) {
  
  if(!is.null(plan)) {
    if(class(plan)=='list') {
      values <- get.values(plan[[1]])
      ct <- get.ct(plan[[1]])
      contours <- get.contours(plan[[1]])
      isocenter <- get.isocenter(plan[[1]])
      plan.name <- rep('', length(plan))
      for(i in 1:length(plan)) {
        plan.name[i] <- plan[[i]]$name
      }
    } else {
      values <- get.values(plan)
      ct <- get.ct(plan)
      contours <- get.contours(plan)
      isocenter <- get.isocenter(plan)
      plan.name <- plan$name
    }
  } else {
    isocenter <- data.frame(x_iso=mean(values$x), y_iso=mean(values$y), z_iso=mean(values$z))
    plan.name <- 'values'
  }
  
  i.plan <- 1 # indice del piano attuale da visualizzare
  
  if(sanitize) {
    values <- sanitize.values(values)
    ct <- sanitize.ct(ct)
  }
  
  # APP
  runApp(
    list(
      
      server=function(input, output)
      { 
        # plot della distribuzione di dose
        output$plot.dose <- renderPlot({
          
          # check per vedere se i.plan è cambiato
          i.plan.new <- which(input$plan==plan.name)
          if(i.plan.new!=i.plan) {
            i.plan <- i.plan.new
            values <- get.values(plan[[i.plan]])
            ct <- get.ct(plan[[i.plan]])
            contours <- get.contours(plan[[i.plan]])
            isocenter <- get.isocenter(plan[[i.plan]])
          }
          
          #my.width <- 640
          #my.height <- 640
          
          if(input$Exit) {stopApp()}
          if(input$plane=='axial (z)') {z <- input$z; y <- NA; x <- NA}
          if(input$plane=='coronal (y)') {z <- NA; y <- input$y; x <- NA}
          if(input$plane=='sagittal (x)') {z <- NA; y <- NA; x <- input$x}
          if(input$show.contours) {my.contours <- contours} else {my.contours <- NULL}
          if(input$show.ct) {
            my.alpha.lower <- input$alpha[1]; my.alpha.upper <- input$alpha[2]
          } else {
            my.alpha.lower <- 1; my.alpha.upper <- 1
          }
          display.slice.all(values=values, ct=ct, contours=my.contours,
                            z=z, x=x, y=y, variable=input$variable,
                            cont=input$show.isolevels, invert.y.axis=input$invert.y.axis,
                            alpha.lower=my.alpha.lower, alpha.upper=my.alpha.upper, HU.window=input$HU.window)       
        })
      },
      
      ui=pageWithSidebar(
        # Application title
        headerPanel('Values (R-Planit interactive GUI)'),
        
        # Sidebar
        sidebarPanel(
          
          # help message
          selectInput("plan", 
                      label = "Plan:",
                      choices = plan.name,
                      selected = plan.name[1]),
          
          hr(),
          
          # variable selection
          h5('Variable selection'),
          selectInput("variable", 
                      label = "Variable to display:",
                      choices = values$variables,
                      selected = values$variables[1]),
          
          checkboxInput("show.isolevels", "Show isolevels",
                        value = FALSE),
          sliderInput("alpha", 
                      label = "alpha value:",
                      min = 0, max = 1, value = c(0.2,1)),
          
          hr(),
          
          # CT display
          h5('CT display setting'),
          checkboxInput("show.ct", "Show CT",
                        value = TRUE),
          checkboxInput("show.contours", "Show contours",
                        value = TRUE),
          sliderInput("HU.window", 
                      label = "HU window:",
                      # min = min(ct$values), max = max(ct$values), value = c(-1000,3000)),
                      min = -1000, max = 3000, value = c(-1000,3000)),
          
          hr(),
          
          # plane selection
          h5('Slice selection'),
          selectInput("plane", 
                      label = "Plane to display:",
                      choices = c("axial (z)", "coronal (y)",
                                  "sagittal (x)"),
                      selected = "axial (z)"),
          
          # slice selection
          sliderInput("z", 
                      label = "z coordinate to display:",
                      min = min(values$z), max = max(values$z), value = isocenter$z_iso),
          sliderInput("y", 
                      label = "y coordinate to display:",
                      min = min(values$y), max = max(values$y), value = isocenter$y_iso),
          sliderInput("x", 
                      label = "x coordinate to display:",
                      min = min(values$x), max = max(values$x), value = isocenter$x_iso),
          
          checkboxInput("invert.y.axis", "Invert y axis",
                        value = TRUE),
          
          hr(),
          
          actionButton('Exit', 'Exit')
          
          
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("plot.dose", width='640px', height='640px')
        )
      )
      
    ), launch.browser=TRUE)
  
}