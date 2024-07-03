# ------------------------------------------------------------------------
# Title: plots.R
# Description: R code to implement the different plots of the thesis.
# Date:
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Set environment. Load required data
load("data/saved_data/data_sf.RData")
load("data/saved_data/plots.RData")
load("data/saved_data/data_df.RData")
load("data/saved_data/mesh.RData")
mesh0$loc<-mesh0$loc*10000
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# PACKAGES
library("ggplot2")
library("sf")
library("sp")
library("reshape2")
library("pals")
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
#                                  0. EXAMPLE PLOTS
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# MATÈRN PLOT


mat <- function(kappa,nu)
{
  return(function(d) geoR::matern(d,1/kappa,nu))
}


theme <- theme(
  panel.background = element_rect(fill = "aliceblue"),
  axis.title.x=element_text(size = 15, margin = margin(t=20)),
  axis.title.y=element_text(size = 15, angle = 90, margin = margin(r=20)),
  plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),
  legend.text=element_text(size=15),
  legend.position="inside",
  legend.justification.inside=c(0.9,0.9),
  legend.background=element_rect(fill = "white")
)

f <- function(i)
{
  geom_function(fun=mat(0.5,nu[i]),n=400,linewidth=1.2,aes(col=paste0("$\\nu=",nu[i],"$")))
}

g <- function(i)
{
  geom_function(fun=mat(kappa[i],0.5),n=400,linewidth=1.2,aes(col=paste0("$\\kappa=",nu[i],"$")))
}

nu <- c(0.1,0.5,0.9)

p_nu <- ggplot() +
  xlim(0,10) +
  ylab("$\\rho(d,\\mathbf{\\theta})$") +
  xlab("Distancia") +
  ggtitle("Función de Matèrn con $\\kappa=0.5$") +
  theme + 
  f(1) + f(2) + f(3) +
  scale_color_manual(NULL,values=brewer.pal(n=3,name="Dark2"))


kappa <- c(0.1,0.5,0.9)
p_kappa <- ggplot() +
  xlim(0,10) +
  ylab("$\\rho(d,\\mathbf{\\theta})$") +
  xlab("Distancia") +
  ggtitle("Función de Matèrn con $\\nu=0.5$") +
  theme + 
  g(1) + g(2) + g(3) +
  scale_color_manual(NULL,values=brewer.pal(n=3,name="Dark2"))

# ------------------------------------------------------------------------
# EXAMPLE GRID PLOT

example.grid.plot <- ggplot() + 
  inlabru::gg(inla.mesh.2d(offset=3,max.edge=3,boundary=aragon/10000,cutoff=.5))+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())



# ------------------------------------------------------------------------
#                             1. SPATIAL PLOTTING
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# ELEVATION AND DISTANCE TO COAST PLOT

# Bound box of Aragon
limits <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = data.frame(X = c(-3, 2), Y = c(39, 43)), 
      data = data.frame(X = c(-3, 2), Y = c(39, 43)),
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

# Theme of plot
theme <- theme(
  panel.background = element_rect(fill = "aliceblue"),
  axis.title.x=element_text(size = 15, margin = margin(t=20)),
  axis.title.y=element_text(size = 15, angle = 90, margin = margin(r=20)),
  plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
  legend.text=element_text(size=10),
  legend.position="right",
  legend.background=element_rect(fill = "white")
)

# Elevation plot
plot_el <- ggplot(data = background) + 
  geom_sf(fill = "antiquewhite") + 
  geom_sf(data = aragon, fill = "antiquewhite") +
  xlab("Longitud") + ylab("Latitude") + ggtitle("Elevación") +
  theme + 
  geom_tile(data = grid_aragon, ggplot2::aes(x = st_coordinates(grid_aragon)[,1], y = st_coordinates(grid_aragon)[,2], fill = elevation)) +
  scale_fill_gradient2(midpoint = 1000, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
                       space = "Lab", limits = c(0, 2800), name = "Elevacion (m)") +
  geom_sf(data = df) +
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))

# Distance to coast plot
plot_dist <- ggplot(data = background) + 
  geom_sf(fill = "antiquewhite") + 
  geom_sf(data = aragon, fill = "antiquewhite") +
  xlab("Longitud") + ylab("Latitud") + ggtitle("Distancia a la costa") +
  theme + 
  geom_tile(data = grid_aragon, ggplot2::aes(x = st_coordinates(grid_aragon)[,1], y = st_coordinates(grid_aragon)[,2], fill = dist)) +
  scale_fill_gradient2(midpoint = 100, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
                       space = "Lab", limits = c(0, 230), name = "Distancia (km)") +
  geom_sf(data = df) +
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))


# Distance and elevation arranged plot
ggpubr::ggarrange(plot_el,plot_dist,labels = c("A","B"),ncol=2,nrow=1)


# ------------------------------------------------------------------------
# COVARIABLES EXTRACTION POINTS PLOT

# Bound limits
limits2 <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = data.frame(X = c(-10.5, 6), Y = c(34.5, 45.5)), 
      data = data.frame(X = c(-10.5, 6), Y = c(34.5, 45.5)),
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

# Exact points
cov.points <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = data.frame(X = c(-10,-10,5,5,-10,-2.5,5,-10,-2.5,5), Y = c(35,45,45,35,35,40,45,45,40,35)), 
      data = data.frame(X = c(-10,-10,5,5,-10,-2.5,5,-10,-2.5,5), Y = c(35,45,45,35,35,40,45,45,40,35)),
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

# Plot
ggplot(data = background) + 
  geom_sf(fill = "antiquewhite") + 
  geom_sf(data = aragon, fill = "antiquewhite") +
  xlab("Longitud") + ylab("Latitud") + ggtitle("Puntos covariables") +
  theme +
  geom_point(aes(x=X,y=Y),color="black",size=2,data=df.stations) +
  geom_sf(data=cov.points,color="red",size=4) +
  geom_path(aes(x=X,y=Y),data=as.data.frame(st_coordinates(cov.points)),color="black") +
  coord_sf(xlim = st_coordinates(limits2)[, 1], ylim = st_coordinates(limits2)[, 2]) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))


# ------------------------------------------------------------------------
# ARAGON MESH PLOT

# Bound limits
limits3 <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = data.frame(X = c(-4, 2.5), Y = c(39, 44)), 
      data = data.frame(X = c(-4, 2.5), Y = c(39, 44)),
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

# Plot
ggplot(data = background) + 
  geom_sf(fill = "antiquewhite") + 
  geom_sf(data = aragon, fill = "antiquewhite") +
  xlab("Longitud") + ylab("Latitud") + ggtitle("Triangulación del dominio ") +
  theme +
  coord_sf(xlim = st_coordinates(limits3)[, 1], ylim = st_coordinates(limits3)[, 2]) + 
  inlabru::gg(mesh0) +
  geom_point(aes(x = X, y = Y), color = "lightblue", size = 0.5,data=df.observ) +
  geom_point(aes(x=X,y=Y),color="black",size=2,data=df.stations) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))


# ------------------------------------------------------------------------
#                             2. EXPLORATORY PLOTTING
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# COVARIABLES VS TEMPERATURE 

theme <- theme(
  panel.background = element_rect(fill = "aliceblue"),
  axis.title.x=element_text(size = 15, margin = margin(t=20)),
  axis.title.y=element_text(size = 15, angle = 90, margin = margin(r=20)),
  plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
  legend.text=element_text(size=15),
  legend.background=element_rect(fill = "white")
)

# ---------------------------------
# Data preparation
station_aux <- df.stations

mean.elev <- as.data.frame(with(df.temp,tapply(temp,station,mean)))
colnames(mean.elev) <- c("mean")
mean.elev$station <- rownames(mean.elev)
rownames(mean.elev) <- NULL

mean.elev <- merge(mean.elev,station_aux, by="station",all.x=T)

# ---------------------------------
# Plots

# Elevation plot
elev_plot <- ggplot(mean.elev,aes(x=altitude,y=mean,color=station)) + geom_point(color="black") +
  labs(x = "Elevación (m)", y = "") +
  ggtitle("Elevación") +
  guides(color="none") + 
  geom_smooth(method='lm',se=T,color="black") +
  theme

# Distance plot
dist_plot <- ggplot(mean.elev,aes(x=dist,y=mean,color=station)) + geom_point(color="black") +
  labs(x = "Distancia al mar (km)", y = "") +
  ggtitle("Distancia al mar") +
  guides(color="none") + 
  geom_smooth(method='lm',se=T,color="black") +
  theme

# Latitude plot
lat_plot <- ggplot(mean.elev,aes(x=lat,y=mean,color=station)) + geom_point(color="black") +
  labs(x = "Latitud (\u00B0C)", y = "") +
  ggtitle("Latitud") +
  guides(color="none") + 
  geom_smooth(method='lm',se=T,color="black") +
  theme

# Longitude plot
lon_plot <-ggplot(mean.elev,aes(x=lon,y=mean,color=station)) + geom_point(color="black") +
  labs(x = "Longitud (\u00B0C)", y = "") +
  ggtitle("Longitud") +
  guides(color="none") + 
  geom_smooth(method='lm',se=T,color="black") +
  theme

# Arranged plot
mean_elev_dist_plot <- ggpubr::ggarrange(elev_plot,dist_plot,lon_plot,lat_plot,ncol=2,nrow=2)
mean_elev_dist_plot <- annotate_figure(mean_elev_dist_plot, top = text_grob("Temperatura media(\u00B0C)", 
                                      color = "black", face = "bold", size = 20))


# ------------------------------------------------------------------------
# MEAN AND SD TEMPERATURE FOR EACH MONTH AND STATION


# ---------------------------------
# Data preparation mean

df.temp$index <- ifelse(df.temp$year<=1985,0,1)
month.mean <- as.data.frame(t(as.data.frame(with(df.temp,tapply(temp,list(month,station,index),mean)))))
colnames(month.mean) <- c("May","Jun","Jul","Ago","Sep")
month.mean$station <- rownames(month.mean)
rownames(month.mean) <- NULL

month.mean.1 <- melt(month.mean[1:18,],id=c("station"))
month.mean.2 <- melt(month.mean[19:36,],id=c("station"))

month.mean.diff <- cbind(month.mean.1[,1:2],month.mean.2[,3]-month.mean.1[,3])
colnames(month.mean.diff) <- c("station","variable","value")

# ---------------------------------
# Plots mean

temp_mean_1 <- ggplot(month.mean.1,aes(x=variable,y=value,color=station, group=station)) + geom_line() +
  theme+
  labs(x="", y = "") +
  scale_y_continuous(limits = c(10, 35)) +
  ggtitle("Media (1956-1985)") +
  guides(color="none")

temp_mean_2 <- ggplot(month.mean.2,aes(x=variable,y=value,color=station, group=station)) + geom_line() +
  theme+
  labs(x="", y = "") +
  scale_y_continuous(limits = c(10, 35)) +
  ggtitle("Media (1986-2015)") +
  guides(color="none")


temp_mean_3 <- ggplot(month.mean.diff,aes(x=variable,y=value,color=station, group=station)) + geom_line() +
  theme+
  labs(x="", y = "") +
  ggtitle("Diferencia") +
  guides(color="none") + 
  geom_hline(yintercept=0, color="black",linewidth=1)

# ---------------------------------
# Data preparation sd

df.temp$index <- ifelse(df.temp$year<=1985,0,1)
month.sd <- as.data.frame(t(as.data.frame(with(df.temp,tapply(temp,list(month,station,index),sd)))))
colnames(month.sd) <- c("May","Jun","Jul","Ago","Sep")
month.sd$station <- rownames(month.sd)
rownames(month.sd) <- NULL

month.sd.1 <- melt(month.sd[1:18,],id=c("station"))
month.sd.2 <- melt(month.sd[19:36,],id=c("station"))

month.sd.diff <- cbind(month.sd.1[,1:2],month.sd.2[,3]/month.sd.1[,3])
colnames(month.sd.diff) <- c("station","variable","value")

# ---------------------------------
# Plots sd


temp_sd_1 <- ggplot(month.sd.1,aes(x=variable,y=value,color=station, group=station)) + geom_line() +
  theme +
  labs(x="", y = "") +
  #scale_y_continuous(limits = c(10, 35)) +
  ggtitle("Desviación (1956-1985)") +
  guides(color="none")

temp_sd_2 <- ggplot(month.sd.2,aes(x=variable,y=value,color=station, group=station)) + geom_line() +
  theme +
  labs(x="", y = "") +
  #scale_y_continuous(limits = c(10, 35)) +
  ggtitle("Desviación  (1986-2015)") +
  guides(color="none") 


temp_sd_3 <- ggplot(month.sd.diff,aes(x=variable,y=value,color=station, group=station)) + geom_line() +
  theme+
  labs(x="", y = "") +
  ggtitle("Cociente") +
  guides(color="none") + 
  geom_hline(yintercept=1, color="black",linewidth=1)

# ---------------------------------
# Merged plot

temp_plot_sd <- ggarrange(temp_mean_1,temp_mean_2,temp_mean_3,temp_sd_1,temp_sd_2,temp_sd_3,labels = c("A","B","C","D","E","F"),ncol=3,nrow=2)
temp_plot_sd_ann <- annotate_figure(temp_plot_sd, top = text_grob("Temperatura (\u00B0C) por meses.", 
                                                                            color = "black", face = "bold", size = 20))


# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
#                                   2. RESULTS PLOTS
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Auxiliary Functions

# Load data
load("data/results/results-pred-5.RData")

# Create a ggplot given a field z.
plot_points <- function(z,m_lim,M_lim,theme,title="",points_p=NULL)
{
  library(fields)
  points <- df.observ[,c("X","Y")]
  points <- 10000*points
  limits <- st_transform(
    as(
      SpatialPointsDataFrame(
        coords = data.frame(X = c(-3, 2), Y = c(39, 43)), 
        data = data.frame(X = c(-3, 2), Y = c(39, 43)),
        proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
      'sf'
    ),
    2062
  )
  
  mid <- (m_lim+M_lim)/2
  
  
  if (!is.null(points_p)){
    extra.points <- geom_point(data=points_p,aes(x=X*10000,y=Y*10000,
                                                 fill=last-first),
                               color="black",pch=24,size=2.2,stroke=1.4)
  } else {
    extra.points <- NULL
  }
  

  # Plot elevation
  (ggplot(data = background) + 
      geom_sf(fill = "antiquewhite") + 
      geom_sf(data = aragon, fill = "antiquewhite") +
      xlab("Longitud") + ylab("Latitud") + ggtitle(title) +
      theme +
      geom_raster(data = grid_aragon, ggplot2::aes(x = points[,1], y = points[,2], fill = z)) +
      scale_fill_gradient2(midpoint = mid, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
                           space = "Lab", limits = c(m_lim,M_lim), name = "ºC") +
      scale_color_gradient2(midpoint = mid, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
                            space = "Lab", limits = c(m_lim,M_lim), name = "ºC") +
      coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) + 
      extra.points +
      scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
      scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
    )
  
}

# ------------------------------------------------------------------------
#  MEAN PLOT BY MONTHS

mean_plot <- function()
{
  
  # Station points
  first_d <- df.temp[df.temp$year %in% 1956:1965,c("temp","station","month")]
  mean_fd <- as.data.frame(tapply(first_d$temp,FUN=mean,
                                  INDEX=list(as.factor(first_d$station),
                                             as.factor(first_d$month))))
  colnames(mean_fd)<-c(paste0("first",seq(5,9)))
  
  last_d <- df.temp[df.temp$year %in% 2006:2015,c("temp","station","month")]
  mean_ld <- as.data.frame(tapply(last_d$temp,FUN=mean,
                                  INDEX=list(as.factor(last_d$station),
                                             as.factor(last_d$month))))
  colnames(mean_ld)<-c(paste0("last",seq(5,9)))
  
  l_points <- list()
  for(i in seq(5,9)){
    l_points[[paste0("month",i)]] <- as.data.frame(cbind(first=mean_fd[,i-4],last=mean_ld[,i-4],
                                           X=df.stations$X,Y=df.stations$Y))
  }
  
  
  aux_mean <- merged_list$mean
  
  days_month_summer<-c(31,30,31,31,30)
  
  mean_diff<-data.frame(mayo=1:n_grid_points,junio=1:n_grid_points,
                        julio=1:n_grid_points,agosto=1:n_grid_points,
                        septiembre=1:n_grid_points)
  
  for(i in 5:9){
    
    detect_fd <- paste0("/",i,"/(195[6-9]|196[0-6])$")
    
    mean_fd_v <- aux_mean[,str_detect(colnames(aux_mean),detect_fd)]
    mean_fd <- apply(mean_fd_v,MARGIN=1,FUN=mean)
    
    
    detect_ld <- paste0("/",i,"/(200[5-9]|201[0-5])$")
    
    mean_ld_v <- aux_mean[,str_detect(colnames(aux_mean),detect_ld)]
    mean_ld <- apply(mean_ld_v,MARGIN=1,FUN=mean) 
    
    mean_diff[,i-4] <- mean_ld-mean_fd
  }
  
  library(ggpubr)
  
  theme <- theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.title.x=element_text(size = 10, margin = margin(t=5)),
    axis.title.y=element_text(size = 10, angle = 90, margin = margin(r=5)),
    plot.title = element_text(hjust = 0.5,size = 10, face = "bold"),
    legend.text=element_text(size=10),
    plot.margin = unit(c(0.1,0.1,0.1,0.1),"lines")
  )
  
  pMy <- plot_points(mean_diff[,1],m_lim=-8,M_lim=8,theme=theme,
                     title = "Mayo",points_p = l_points[[1]])
  pJn <- plot_points(mean_diff[,2],m_lim=-8,M_lim=8,theme=theme,
                     title = "Junio",points_p = l_points[[2]])
  pJl <- plot_points(mean_diff[,3],m_lim=-8,M_lim=8,theme=theme,
                     title = "Julio",points_p = l_points[[3]])
  pAg <- plot_points(mean_diff[,4],m_lim=-8,M_lim=8,theme=theme,
                     title = "Agosto",points_p = l_points[[4]])
  pSp <- plot_points(mean_diff[,5],m_lim=-8,M_lim=8,theme=theme,
                     title = "Septiembre",points_p = l_points[[5]])
  
  return(ggarrange(pMy,pJn,pJl,pAg,pSp,nrow=2,ncol=3,common.legend = T,
                   legend="bottom"))
  
}

mean_plot()

# ------------------------------------------------------------------------
#  TEMPORAL SERIES PLOT


mean_temp <- function()
{
  
  # Station points
  first_d <- df.temp[,c("temp","year","month")]
  mean_years <- as.data.frame(tapply(first_d$temp,FUN=mean,
                                  INDEX=list(as.factor(first_d$year),
                                             as.factor(first_d$month))))
  colnames(mean_years)<-c(paste0("month",seq(5,9)))
  
  l_points <- list()
  for(i in seq(5,9)){
    l_points[[paste0("month",i)]] <- as.data.frame(cbind(t=mean_years[,i-4],
                                                         years=1956:2015))
  }
  
  aux_mean <- merged_list$mean
  
  aux_mean <- apply(aux_mean,FUN=mean, MARGIN = 2)
  
  months_days <- c(31,30,31,31,30)
  ind.month <- rep(rep(1:5,months_days),n_years)
  ind.year <- rep(1:n_years,each=n_days_year)
  df_mean <- cbind(t=aux_mean,month=ind.month,year=ind.year)
  
  mean <- t(tapply(aux_mean,FUN=mean,INDEX = list(ind.month,ind.year)))
  
  
  theme <- theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.title.x=element_text(size = 10, margin = margin(t=5)),
    axis.title.y=element_text(size = 10, angle = 90, margin = margin(r=5)),
    plot.title = element_text(hjust = 0.5,size = 10, face = "bold"),
    legend.text=element_text(size=10),
    legend.position="inside",
    legend.justification.inside=c(0.9,0.9),
    legend.background=element_rect(fill = "white")
  )
  
  plots_l <- list()

  names_month <- c("Mayo","Junio","Julio","Agosto","Septiembre")
  for(i in 1:5){
    df_l <- data.frame(year=1956:2015,mean=mean[,i])
    plots_l[[i]] <- ggplot() + 
      geom_line(data=df_l,aes(x=year,y=mean,colour="Predicciones")) + 
      geom_line(data=l_points[[i]],aes(x=years,y=t,colour="Localizaciones")) +
      scale_color_manual(NULL, values = c("Localizaciones" = "#7570B3", 
                                          "Predicciones" = "#D95F02")) +
      theme + ggtitle(names_month[i]) + 
      geom_smooth(data=df_l,aes(x=year,y=mean), method='lm',se=T,color="black")
    
  }
  
  ggarrange(plotlist = plots_l,common.legend = T,legend="bottom")

}

mean_temp()

# ------------------------------------------------------------------------
#  DIFFERENCES QUANTILES PLOT

mean_diff_quantiles <- function()
{
  diff <- merged_list$'q0.975'-merged_list$'q0.025'
  mean_diff <- apply(diff,FUN=mean,MARGIN=1)
  
  theme <- theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.title.x=element_text(size = 10, margin = margin(t=5)),
    axis.title.y=element_text(size = 10, angle = 90, margin = margin(r=5)),
    plot.title = element_text(hjust = 0.5,size = 15),
    legend.text=element_text(size=10),
    legend.background=element_rect(fill = "white")
  )
  
  points <- df.stations[,c("X","Y")]*10000
  
  plot_points(mean_diff,min(mean_diff),max(mean_diff),theme) + 
    geom_point(data=points,aes(x=X,y=Y),color="black",size=3,shape=17) +
    ggtitle("$q_{0.975}(\\mathbf{x})-q_{0.025}(\\mathbf{x})$")
  
}

mean_diff_quantiles()

# ------------------------------------------------------------------------
#  FIELD PLOT - MODEL 3

# Load required data
load("C:/Users/david/Desktop/TFG/data/Results-2/no-covariates.RData")
library("fields")


plot_field <- function()
{

  # Get prediction grid
  loc <- as.matrix(cbind(df.observ$X,df.observ$Y))
  colnames(loc) <- c("X","Y")
  
  # Get square of projections
  xlim <- c(60,110)
  ylim <- c(50,100)
  proj.grid <- inla.mesh.projector(mesh0,xlim=xlim,ylim=ylim,
                                   dim=c(50,50))
  
  field.proj.mean <- inla.mesh.project(projector=proj.grid,
                                       output.3$summary.random$field$mean)
  field.proj.sd <- inla.mesh.project(projector=proj.grid,
                                     output.3$summary.random$field$sd)
  
  # Get boundary data.frame
  bound_aragon <- sf::st_boundary(aragon)/10000
  bound_df <- as.data.frame(bound_aragon[[1]][[1]])
  colnames(bound_df) <- c("x","y")
  
  # Fix it to ARAGON
  xy.in <- splancs::inout(proj.grid$lattice$loc,bound_df)
  field.proj.mean[!xy.in] <- NA
  field.proj.sd[!xy.in] <- NA
  

  df <- data.frame(sd=as.vector(field.proj.sd),
                   mean=as.vector(field.proj.mean),
                   x=proj.grid$lattice$loc[,1]*10000,
                   y=proj.grid$lattice$loc[,2]*10000)
  df <- df[!is.na(df$mean),]
  
  theme <- theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.title.x=element_text(size = 10, margin = margin(t=5)),
    axis.title.y=element_text(size = 10, angle = 90, margin = margin(r=5)),
    plot.title = element_text(hjust = 0.5,size = 15),
    legend.text=element_text(size=10),
    legend.background=element_rect(fill = "white")
  )
  
  limits <- st_transform(
    as(
      SpatialPointsDataFrame(
        coords = data.frame(X = c(-3, 2), Y = c(39, 43)), 
        data = data.frame(X = c(-3, 2), Y = c(39, 43)),
        proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
      'sf'
    ),
    2062
  )
  
  
  
  mean <- (ggplot(data = background) + 
      geom_sf(fill = "antiquewhite") + 
      geom_sf(data = aragon, fill = "antiquewhite") +
      xlab("Longitud") + ylab("Latitud") + 
      theme + ggtitle("$aqui$") +
      geom_raster(data = df, aes(x = x, y = y, fill = mean)) +
      scale_fill_gradient2(midpoint = 0, low = scales::muted("blue"), 
                           mid = "white", high = scales::muted("red"), 
                           space = "Lab", limits = c(-3,3), name = "ºC") +
      coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) + 
      scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
      scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) + 
      geom_point(data=data.frame(x=df.stations$X*10000,y=df.stations$Y*10000),
                 aes(x=x,y=y),color="black",size=2)
  )
  
  sd <- (ggplot(data = background) + 
             geom_sf(fill = "antiquewhite") + 
             geom_sf(data = aragon, fill = "antiquewhite") +
             xlab("Longitud") + ylab("Latitud") +
             theme + ggtitle("$aqui$") +
             geom_raster(data = df, aes(x = x, y = y, fill = sd)) +
             scale_fill_gradient2(midpoint = 2, low = scales::muted("blue"), 
                                  mid = "white", high = scales::muted("red"), 
                                  space = "Lab", limits = c(1,3), name = "ºC") +
             coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) + 
             scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
             scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) + 
             geom_point(data=data.frame(x=df.stations$X*10000,y=df.stations$Y*10000),
                        aes(x=x,y=y),color="black",size=2)
  )
  
  return(ggarrange(mean,sd,common.legend = F))
}

plot_field()

# ------------------------------------------------------------------------
