# ------------------------------------------------------------------------
# Title: models.R
# Description: R code to implement the different models.
# Date:
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Set environment. Load required data
setwd("C:/Users/david/Desktop/TFG")
load("data/saved_data/plots.RData")
load("data/saved_data/data_df.RData")
load("data/saved_data/bound.RData")
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# PACKAGES
library("INLA")
library("stringr")
library("ggplot2")
library("sf")
# ------------------------------------------------------------------------
# RESCALE DATA

df.observ[,c("X","Y")] <- df.observ[,c("X","Y")]/10000
df.temp[,c("X","Y")] <-  df.temp[,c("X","Y")]/10000
df.stations[,c("X","Y")] <-  df.stations[,c("X","Y")]/10000

df.observ$altitude <- df.observ$elevation
df.observ$elevation <- NULL

# ------------------------------------------------------------------------
# AUXILIARY FUNCTIONS

# Function to create the .csv of the results
to.csv <- function(output,file,round.n)
{
  df <- round(rbind(output$summary.fixed[-7],output$summary.hyperpar),round.n)
  write.table(df,file=file,sep=",",col.names = T,row.names = F)
}

# Function to get dics table
get.dics <- function(file,round.n){
  outputs <- list(output.1,output.2,output.3,output.4,output.5.1)
  f <- function(x)
  {
    return(c(x$dic$dic,x$dic$p.eff,x$dic$mean.deviance))
  }
  df <- round(t(sapply(outputs, FUN=f)),round.n)
  write.table(df,file=file,sep=",",col.names = T,row.names = F)
}

# Function to send a email to notify process end.
mail <- function(o){

  library(blastula)
  

  email <- compose_email(
    body = md("Finalizado")
  )
  

  email %>%
    smtp_send(
      from = "825068@unizar.es",
      to = "825068@unizar.es",
      subject = ifelse(exists(o), "FINALIZADO CORRECTAMENTE", "FINALIZADO INCORRECTAMENTE"),
      credentials = creds_key("gmail")
    )
}


# Function to merge two simulation results in one
merge_results <- function(res1_file,res2_file,result_file)
{
  
  create_data_file <- function(stack,output)
  {
    index.pred <- inla.stack.index(stack,"pred")$data
    data <- output$summary.fitted.values[index.pred,c(1,2,3,5)]
    return(get_fields_points(days_predict,years_predict,n_grid_points,data))
  }
  
  get_fields_points <- function(days,years,n_points,field)
  {
    names_l <- c("mean","sd","q0.025","q0.975")
    l <- list()
    for (i in 1:4){
      field_aux <- field[,i]
      day_year_to_day_month <- function(x)
      {
        return(paste(df.temp[df.temp$day.year==x,][1,c("day.month","month")],collapse="/"))
      }
      days.formated <- sapply(days,FUN=day_year_to_day_month)
      names <- apply(expand.grid(days.formated,years),MARGIN=1,FUN=function(x) paste(x,collapse="/"))
      n <- length(names)
      aux <- sapply(seq(n),FUN=function(x) field_aux[seq(x,length(field_aux),by=n)])
      colnames(aux) <- names
      l[[names_l[i]]] <- as.data.frame(aux)   
    }
    return(l)
  }
  
  
  days_predict <- seq(1,n_days_year,by=1) + 120
  
  load(res1_file)
  gc()
  
  years_predict <- seq(1,n_years,by=2) + 1955
  list.5.1 <- create_data_file(stack.5,output.5.1,"result-list.5.1.RData")
  
  rm(output.5.1,stack.5)
  gc()
  
  load(res2_file)
  gc()
  
  years_predict <- seq(2,n_years,by=2) + 1955
  list.5.2 <- create_data_file(stack.5,output.5.2,"result-list.5.2.RData")
  
  rm(output.5.2,stack.5)
  gc()
  
  merge_results <- function(aux1_1,aux2_2)
  {
    aux <- cbind(aux1_1[,1:n_days_year],aux2_2[,1:n_days_year])
    for(i in 1:29){
      ind <- i*n_days_year + 1:n_days_year
      aux <- cbind(aux,aux1_1[,ind],aux2_2[,ind])
    }
    return(aux)
  }
  
  names_list <- c("mean","sd","q0.025","q0.975")
  merged_list <- list()
  for(i in 1:4){
    merged_list[[names_list[i]]] <- merge_results(list.5.1[[i]],list.5.2[[i]])
  }
  
  save(merged_list,file=result_file) 
}

# Animate a results df
create_animation <- function(df,dir)
{
  names <- colnames(df)
  for(i in seq(dim(df)[2])){
    plot <- plot_points(df[,i],paste0("Temperature ",names[i]))
    ggsave(filename=paste0(dir,"/",i,".png"),plot=plot,width = 584*4 ,height = 383*4,units="px") 
  }
  
}


# ------------------------------------------------------------------------
#                           DATA PREPARATION
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# FIXED EFFECTS

# Selection of atmospheric covariables. 

regex.sel <- "(xt)_.{2}_(700)hPa$"

sel.ind.atm <- str_detect(colnames(df.temp),regex(regex.sel))
cov.atm <- paste(names(df.temp[,sel.ind.atm]),collapse = "+")

sel.ind <- sel.ind.atm
sel.ind[c(7,8,12,14)] <- TRUE
data.cov <- as.list(df.temp[,sel.ind])

cov <- "1 + year + sin(day.year*2*pi/360) + cos(day.year*2*pi/360) + I(altitude/100) + I(dist)"
fixed.effects <- as.formula(paste0("temp~",
                                   paste(cov,cov.atm,
                                         sep = ifelse(cov.atm=="",""," + "))))


# ------------------------------------------------------------------------
# DATA CONSTANTS

n_stations <- length(df.stations$station)
n_data <- length(df.temp$station)
n_days <- n_data/n_stations
n_days_year <- length(unique(df.temp$day.year))
n_years <- length(unique(df.temp$year))
n_grid_points <- length(df.observ$X)


# Matrix with the coordinates of all the data
coordinates.all.year <- as.matrix(df.stations[df.temp$station,c("X","Y")])
rownames(coordinates.all.year) <- NULL


# ------------------------------------------------------------------------
# MESH CREATION AND MATERN SPECIFICATION

mesh0 <- inla.mesh.2d(offset=c(5,7),max.edge =c(4,7),boundary=aragon/10000,cutoff = 10)

prior.range <- c(sqrt(diff(range(df.observ$X))^2+diff(range(df.observ$Y))^2),0.95)
prior.sigma <- c(3,.01)
spde <- inla.spde2.pcmatern(mesh0,alpha=2,
                            prior.range = prior.range ,
                            prior.sigma = prior.sigma)


# ------------------------------------------------------------------------
# PREDICTION PREPARATION

days_predict_ind <- seq(1,n_days_year,by=1)
years_predict_ind <- seq(1,n_years,by=2)


days_predict <- days_predict_ind + 120
years_predict <- years_predict_ind + 1955

n_days_est <- length(years_predict_ind)*length(days_predict_ind)

st <- expand.grid(days_predict,years_predict)
colnames(st) <- c("day.year","year")
st <- st[rep(1:n_days_est,n_grid_points),]

pred <- df.observ[rep(1:n_grid_points,each=n_days_est),]

atm.pred <- (df.temp[df.temp$year %in% years_predict & 
                       df.temp$day.year %in% days_predict ,sel.ind.atm][rep(1:n_days_est,n_grid_points),])


pred <- cbind(pred,st,atm.pred)

pred.loc <- as.matrix(pred[,c("X","Y")])
pred.group <- as.matrix(pred[,"day.year"])
colnames(pred.group) <- "day.year"
pred.rep <- as.matrix(pred[,"year"])
colnames(pred.rep) <- "year"

pred.cov <- as.list(pred[,!(names(pred) %in% c("X","Y","lat","lon"))])

rm(pred)

# ------------------------------------------------------------------------
#                                 MODELLING
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# MODEL 1
# Just fixed effects
formula.1 <- fixed.effects

output.1 <- inla(formula.1,
                 data= as.list(df.temp),
                 family="gaussian",
                 control.compute=list(dic=F)
)
summary(output.1)


# ------------------------------------------------------------------------
# MODEL 2
# Fixed effects + ar1 in day/year

formula.2 <- update(fixed.effects, . ~ . + 
                      f(day.year,
                        replicate=year,
                        model="ar1",
                        hyper = list(prec = list(initial = log(0.06)))
                      )
)


output.2 <- inla(formula.2,
                 data= as.list(df.temp),
                 family="gaussian",
                 control.compute = list(dic=T)
)

summary(output.2)


# ------------------------------------------------------------------------
# MODEL 3
# Fixed effects + spatial coefficient

Proj.Mat.3 <- inla.spde.make.A(mesh=mesh0,
                               loc=coordinates.all.year)

mesh.index.3 <- inla.spde.make.index(name="field",
                                     n.spde=spde$n.spde)

stack.est.3 <- inla.stack(data=list(temp=df.temp$temp),
                          A=list(Proj.Mat.3,1),
                          effects=list(c(mesh.index.3,list(Intercept=1)),
                                       data.cov),
                          tag="est")

formula.3 <- update(fixed.effects,.~.-1+Intercept + 
                      f(field,model=spde,diagonal = 1e-5))


output.3 <- inla(formula.3,
                 data=inla.stack.data(stack.3,spde=spde),
                 family = "gaussian",
                 control.predictor=list(A=inla.stack.A(stack.3)),
                 control.compute = list(dic=TRUE),
                 verbose=T,
                 inla.mode="experimental",
                 
)

summary(output.3)

# ------------------------------------------------------------------------
# MODEL 4
# Fixed effects + ar1 in day/year + spatial effect

formula.4 <- update(formula.3,.~ . + f(day.year,
                                       replicate=year,
                                       model="ar1",
                                       hyper = list(prec = 
                                                      list(initial = 
                                                             log(0.251)))
)
)

output.4 <- inla(formula.4,
                 data=inla.stack.data(stack.est.3,spde=spde),
                 control.predictor=list(A=inla.stack.A(stack.est.3),compute=T),
                 control.compute = list(dic=TRUE), 
                 verbose=T
)

summary(output.4)



# ------------------------------------------------------------------------
# MODEL 5
# Fixed effects + spatial ar1 in day/year 


# Group is the time index. In our case is the day.total
# n.group is the number of time groups
Proj.Mat.5 <- inla.spde.make.A(mesh=mesh0,
                               loc=coordinates.all.year,
                               group=df.temp$day.year.ind,
                               n.group=n_days_year,
                               repl = df.temp$ind.year,
                               n.repl = n_years)

#Make the temporal index
temporal.index.5 <- inla.spde.make.index(name="spatial.field",
                                         n.spde=spde$n.spde,
                                         n.group=n_days_year,
                                         n.repl=n_years)

stack.est.5 <- inla.stack(data=list(temp=df.temp$temp),
                          A=list(Proj.Mat.5,1),
                          effects=list(c(temporal.index.5,list(Intercept=1)),
                                       data.cov),
                          tag="est")




Proj.Mat.Pred.5 <- inla.spde.make.A(mesh=mesh0,
                                    loc=pred.loc,
                                    group=pred.group-120,
                                    n.group=n_days_year,
                                    repl=pred.rep-1955,
                                    n.repl=n_years)

stack.pred.5 <- inla.stack(data=list(temp=NA),
                         A=list(Proj.Mat.Pred.5,1),
                         effects=list(c(temporal.index.5,list(Intercept=1)),
                                      pred.cov),
                         tag="pred")

stack.5 <- inla.stack(stack.est.5,stack.pred.5)

# Remove to free RAM memory
rm(stack.pred.5,Proj.Mat.5,Proj.Mat.Pred.5,covariate_matrix_std,data.cov)
gc()

# Formula specification
formula.5 <- update(fixed.effects, .~. -1 + Intercept +
                      f(spatial.field,model=spde,group=spatial.field.group,
                        control.group = list(model="ar1"),
                        replicate=spatial.field.repl)
)

# Dummy training to get initial values
initial2 <- inla(formula.5,
                 data = inla.stack.data(stack.5,spde=spde),
                 family="gaussian",
                 control.predictor = list(A=inla.stack.A(stack.5)),
                 control.compute = list(openmp.strategy="huge"),
                 control.inla =  list(diagonal = 10, strategy = "gaussian", int.strategy = "eb"),
                 verbose=T,
                 num.threads=16,
                 inla.mode="experimental"
)

# Real computation. It could take days.
# Model needed 100GB RAM 
output.5 <- inla(formula.5,
                 data = inla.stack.data(stack.5,spde=spde),
                 family="gaussian",
                 control.predictor = list(A=inla.stack.A(stack.5)),
                 control.compute = list(dic=T,openmp.strategy="huge"),
                 verbose=T,
                 control.mode=list(result=initial2,restart=T),
                 num.threads = 16,
                 inla.mode="experimental"
)

mail("output.5")
save(output.5,file="results2-2.RData")
load("data/Results/results2-2.RData")

# ------------------------------------------------------------------------