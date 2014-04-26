# APPS -------------------------------------------------------------------------
# (GUI ed altro)

#' Values visualization WEB-GUI (dataframe ops)
#' 
#' @param plan The plan object
#' @family Apps
#' @export
#' @import shiny
values.app <- function(plan=NULL, values=NULL, ct=NULL, contours=NULL, sanitize=TRUE) {
  
  if(!is.null(plan)) {
    values <- get.values(plan)
    ct <- get.ct(plan)
    contours <- get.contours(plan)
    isocenter <- get.isocenter(plan)
  } else {
    isocenter <- data.frame(x_iso=mean(values$x), y_iso=mean(values$y), z_iso=mean(values$z))
  }
  
  if(sanitize) {
    values <- sanitize.values(values)
    ct <- sanitize.ct(ct)
  }
  
  runApp(
    list(
      
      server=function(input, output)
      {
        # plot della distribuzione di dose
        output$plot.dose <- renderPlot({
          if(input$plane=='axial (z)') {z <- input$z; y <- NA; x <- NA}
          if(input$plane=='coronal (y)') {z <- NA; y <- input$y; x <- NA}
          if(input$plane=='sagittal (x)') {z <- NA; y <- NA; x <- input$x}
          display.slice.all(values=values, ct=ct, contours=contours,
                            z=z, x=x, y=y, variable=input$variable,
                            cont=input$show.isolevels, invert.y.axis=input$invert.y.axis)       
        })
      },
      
      ui=pageWithSidebar(
        # Application title
        headerPanel('R-Planit interactive GUI'),
        
        # Sidebar
        sidebarPanel(
          
          # help message
          h5(plan$name),
          
          hr(),
          
          # variable selection
          h5('Variable selection'),
          selectInput("variable", 
                      label = "Variable to display:",
                      choices = values$variables,
                      selected = values$variables[1]),
          
          checkboxInput("show.isolevels", "Show isolevels",
                        value = FALSE),
          
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
                        value = TRUE)
          
          
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("plot.dose", width = "544px", height='544px')
        )
      )
      
    ), launch.browser=TRUE)
  
}