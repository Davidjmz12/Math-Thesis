# ------------------------------------------------------------------------
# Title: df_creation.R
# Description: R code to create the different data.frames needed for 
#              modeling.
# Date:
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Set environment. Load required data
load("data/saved_data/data_sf.RData")
temp <-read.table("data/text_files/temperature_aragon.txt",header=T)
stations <- read.table("data/text_files/stations_info.txt")
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# PACKAGES
library("sf")
library("sp")
library("NcLibrary")
library("stringr")
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# MAIN

# Function to get a subset of our data.
subset_df <- function(temp,months,year_min,year_max)
{
  return(temp[temp$year>=year_min & temp$year<=year_max & temp$month %in% months,])
}


# Function that returns a data.frame with the temperatures with date and place
get_table_from_df_and_stations <- function(df,stations)
{
  aux_df <- data.frame()
  for(i in 5:dim(df)[2])
  {
    aux_rows <- df[,c(1:4,i)]
    aux_rows$station <- colnames(aux_rows)[5]
    aux_rows <- cbind(aux_rows,stations[which(stations$name==colnames(aux_rows)[5]),][2:4])
    colnames(aux_rows)[5] <- "temp"
    
    aux_df <- rbind(aux_df,aux_rows)
  }
  
  aux_df$station <- as.factor(aux_df$station)
  aux_df$temp <- aux_df$temp/10
  return (aux_df)
}

# Subset study dates
df.temp <- subset_df(temp,c(5,6,7,8,9),1956,2015)
df.temp <- get_table_from_df_and_stations(df.temp,stations)

# Add sea distance
df.temp <- merge(df.temp,as.data.frame(df)[,-2],by.x=c("lon","lat","station","altitude"),by.y=c("x","y","name","elevation"),all.x = T)
df.temp$day.ind.global <- as.integer(rep(1:9180))

# Add day in year index
df.temp$day.year.ind <- 1:(length(unique(df.temp$day.year))-1)

# Create stations data.frame
df.stations <- as.data.frame(st_coordinates(df))
df.stations$station <- df$name

# Add to data.frame X and Y
df.temp <- merge(df.temp, df.stations, by.x= "station", by.y= "station", all.x=T)
df.stations <-merge(df.stations, stations, by.x= "station", by.y= "name")
df.stations$dist <- as.data.frame(df)$dist

# Create data.frame observation
df.observ <- as.data.frame(cbind(grid_aragon,st_coordinates(grid_aragon)))[,-7]
colnames(df.observ)[1:2] <- c("lon","lat")

df.temp$ind.year <- as.integer(df.temp$year-1955)

# Rearrange temp
df.temp <- df.temp[,c(9,1,2,3,13,14,4,5,15,6,7,8,12,10,11)]


# Extract covariables
configuration <- list(dim=list(
                            latitude=list(type="lat", coord_360=T, name="latitude", format=0),
                            longitude=list(type="lon", coord_360=T, name="longitude", format=0)
                      ),
                      time=list(extended=F,time_div=T,name="time"),
                      sep="_",
                      nc_file="data/text_files/covariables.nc")

df.cov <- NcLibrary::read_nc_file(configuration)


# Subset for desired coordinate points and variable
acc <- str_detect(colnames(df.cov),regex("((t|q|v|u)_((-10|5)E_(35|45)N_))|((t|q|v|u)_-3E_40N_)|DAY|MONTH|YEAR|DATE"))
df.cov <- df.cov[,acc]
df.cov <- df.cov[df.cov$MONTH %in% 5:9,]

#Replicate for each station
n_stations <- 18
df.cov <- do.call(rbind, replicate(n_stations, df.cov, simplify = FALSE))


#Renaming the data.frame
M <- matrix(data=c("_83$","_300hPa","_95$","_500hPa","_105$","_700hPa","-10E_45N","TL",
                   "-10E_35N","BL","5E_45N","TR","5E_35N","BR","-3E_40N","C"),
            ncol=2,nrow=8,byrow=T)

for(i in seq_len(nrow(M))){
  colnames(df.cov) <- str_replace(colnames(df.cov),regex(M[i,1]),M[i,2])
}

# Re-scaling the variables

# Temperature
# Kelvin to ºC
temp.ind <- str_detect(colnames(df.cov),regex("t"))
df.cov[,temp.ind] <- df.cov[,temp.ind]-273.15

# Relative Humidity
# Kg * Kg^-1 to g * Kg^-1
hum.ind <- str_detect(colnames(df.cov),regex("q"))
df.cov[,hum.ind] <- df.cov[,hum.ind]*1000

# Wind speed 
# m * s^-1 equal
wind.ind <- str_detect(colnames(df.cov),regex("u|v"))
df.cov[,wind.ind] <- df.cov[,wind.ind]


df.temp <- cbind(df.temp,df.cov[,5:length(df.cov)])
df.temp$day.year <- 121:273


# save data
save(df.temp,df.stations,df.observ,file="data/saved_data/data_df.RData") 

# ------------------------------------------------------------------------