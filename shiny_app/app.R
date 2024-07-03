library("shiny")
library("ggplot2")
library("sf")
library("fields")
library("sp")
library("patchwork")

load("results-pred-5.RData")
load("useful_data.RData")

theme <- theme(
  panel.background = element_rect(fill = "aliceblue"),
  axis.title.x=element_text(size = 10, margin = margin(t=5)),
  axis.title.y=element_text(size = 10, angle = 90, margin = margin(r=5)),
  plot.title = element_text(hjust = 0.5,size = 10, face = "bold"),
  legend.text=element_text(size=10)
)

geom.prev <-
  (ggplot(data = background) + 
    geom_sf(fill = "antiquewhite") + 
    geom_sf(data = aragon, fill = "antiquewhite") +
    xlab("") + ylab("") + 
    theme
  )

plot_points <- function(z,m_lim,M_lim,limit.x,limit.y,title="",obs.points=NULL)
{
  limits <- st_transform(
    as(
      SpatialPointsDataFrame(
        coords = data.frame(X = limit.x, Y = limit.y), 
        data = data.frame(X = limit.x, Y = limit.y),
        proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
      'sf'
    ),
    2062
  )
  mid <- (m_lim+M_lim)/2
  
  if (!is.null(obs.points)){
    extra.points <- geom_point(data=obs.points,aes(x=X,y=Y,
                                                 fill=temp),
                               color="black",pch=24,size=5,stroke=1)
  } else {
    extra.points <- NULL
  }
  
  # Plot elevation
  geom.prev +
    ggtitle(title) +
    geom_raster(data = grid_aragon, aes(x = points[,1], y = points[,2], fill = z)) +
    scale_fill_gradient2(midpoint = mid, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
                         space = "Lab", limits = c(m_lim,M_lim), name = "ºC") +
    scale_color_gradient2(midpoint = mid, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
                          space = "Lab", limits = c(m_lim,M_lim), name = "ºC") +
    coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2],expand = T) + 
    extra.points
  
}

server <- function(input, output) {
  
  # Reactive expression to compose a data frame containing all of the values
  
  plotGg <- reactive({
    
    input.date <- input$date_slider
    limit.x <- c(-3,2)
    limit.y <- c(39, 43)
    
    day <- input$day_slider
    modes <- as.array(input$sel)
  
      
    date.string <- paste(input$day_slider,input$month_slider,
                         input$year_slider,sep="/")   
    index <- colnames(merged_list$mean)==date.string
    new_plots <- list()
    for(i in 1:length(modes))
    {
      f <- function(i)
      {
        zmin <- ifelse(modes[i]=="sd",input$rangeSd[1],input$rangeTemp[1])
        zmax <- ifelse(modes[i]=="sd",input$rangeSd[2],input$rangeTemp[2])
        z <- merged_list[[modes[i]]][,index]
        if(input$points & modes[i]!="sd"){
          ind <- seq(0,17)*9180 + which(index)
          obs.points <- df.tem[,c("X","Y","temp")]
        } else{
          obs.points <- NULL
        }
        plot_points(z,zmin,zmax,limit.x,limit.y,modes[i],
                    obs.points = obs.points)
      }  
      new_plots[[i]] <- f(i)
    }
    ncol <- ifelse(length(new_plots)<=1,1,2)
    nrow <-ifelse(length(new_plots)<=2,1,2)
    wrap_plots(new_plots) + plot_layout(ncol=ncol,nrow = nrow)
  }) 
  
  
  output$plot <- renderPlot({
    plotGg()
  },width = "auto",height = 700)


}


# Define UI for slider demo application
ui <- fluidPage(
  
  #  Application title
  titlePanel("Maximum Temperature Results"),
  
  # Sidebar with sliders that demonstrate various available options
  sidebarLayout(
    sidebarPanel(
      sliderInput("day_slider", "Select a Day:",
                  min = 1,
                  max = 30,
                  value = 1,
                  step = 1),
      sliderInput("month_slider", "Select a Month:",
                  min = 5,
                  max = 9,
                  value = 5,
                  step = 1),  # You can adjust the step size (e.g., 1 day)
      sliderInput("year_slider", "Select a Year:",
                  min = 1956,
                  max = 2015,
                  value = 1956,
                  step = 1), 
      sliderInput("rangeTemp", "Range Temperatura:",
                  min = -18, max = 55, value = c(-18,55),step=1),
      sliderInput("rangeSd", "Range sd:",
                  min = 0, max = 5, value = c(0,5),step=.1),
      checkboxInput("points",label = "Observed stations",value=F),
      checkboxGroupInput("sel", 
                   "Select data:", 
                   c("mean","sd","q0.025","q0.975"),
                   selected = c("mean")),
      
    ),
    
    
    # Show a table summarizing the values entered
    mainPanel(
      plotOutput("plot")
    )
  )
)

shinyApp(ui=ui,server=server)
