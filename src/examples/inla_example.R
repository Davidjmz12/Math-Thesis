# ---------------------------------------------------------------------
# Title: example_inla.R
# Description: R code to implement the INLA algorithm example.
# Date:
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# PACKAGES

library("INLA")
library("numDeriv")
library("ggplot")

# ---------------------------------------------------------------------
# Simulation of random data
actual.mean <- 0.5
actual.sd <- 0.7
n <- 200
y <- rnorm(n,actual.mean,actual.sd)


# Likelihood density
likelihood <- function(theta,x)
{
  return(prod(dnorm(y,x,1/sqrt(theta))))
}

# Prior density
prior <- function(theta)
{
  return(dgamma(theta,2,1))
}

# Latent field density
latent.field <- function(x)
{
  return(dnorm(x))
}


# ---------------------------------------------------------------------
# Full posterior on x characterization
mean.full.post.x <- function(theta) return((theta*sum(y))/(n*theta+1))
var.full.post.x <- function(theta) return(1/(n*theta+1))
dens.full.post.x <- function(x,theta) return(dnorm(x,mean.full.post.x(theta),
                                                   sqrt(var.full.post.x(theta))))

# ---------------------------------------------------------------------

# Posterior of theta before normalization
post.theta.no.norm <- function(theta)
{
  mode <- mean.full.post.x(theta)
  
  operation <- prior(theta)*latent.field(mode)*likelihood(theta,mode)*
    sqrt(2*pi*var.full.post.x(theta))
  
  return(as.numeric(operation))
}

post.theta.no.norm.v <- function(theta) sapply(theta,post.theta.no.norm)


post.theta.norm <- integrate(post.theta.no.norm.v,0,10)$value

# Posterior of theta after normalization
dens.post.theta <- function(theta)
{
  return(post.theta.no.norm.v(theta)/post.theta.norm)
}

# Mode of posterior of theta
pos.theta.mode <- optimise(dens.post.theta,interval=c(0,5),maximum=TRUE)$maximum

# ---------------------------------------------------------------------
# Compute grid of thetas

#Negative hessian
H <- -hessian(dens.post.theta,c(pos.theta.mode))
H.inv <- as.numeric(1/H)

# The eigenvalue descomposition of H.inv is exactly 1 * H.inv * 1
z.param <- function(z)
{
  as.numeric(pos.theta.mode+sqrt(H.inv)*z)
}

# Longitude step and limit definitions
longitude.step <- 0.5
limit <- 3

# Grid of thetas
thetas <- c(pos.theta.mode)

log.dens <- function(theta) log(dens.post.theta(theta))

for(dir in c(-1,1)){
  new.dir <- 0+dir*longitude.step
  new.t <- z.param(new.dir)
  diff <-  abs(log.dens(z.param(0))-
                 log.dens(new.t))
  while(diff<limit & new.t>0)
  {
    thetas <- c(thetas,new.t)
    new.dir <- new.dir+dir*longitude.step
    new.t <- z.param(new.dir)
    diff <-  abs(log.dens(z.param(0))-
                   log.dens(new.t))
  }
}
thetas <- sort(thetas)

# ---------------------------------------------------------------------
# Evaluation at grid of thetas. Computing posterior of latent field

dis.integrand <- function(x,theta)
{
  dens.full.post.x(x,theta)*dens.post.theta(theta)
}

aux.fun <- function(x)
{
  mult <- sapply(thetas, function(theta) dis.integrand(x,theta))
  sum(mult)
}
aux.fun.v <- function(x) sapply(x,aux.fun)
pos.x.norm <- integrate(aux.fun.v,-1,1)$value

# Posterior density function of latent field after normalization
dens.post.x <- function(x)
{
  return(aux.fun.v(x)/pos.x.norm)
}

# Posterior latent field mode
pos.x.mode <- optimize(dens.post.x,c(0,1),maximum=T)$maximum

# ---------------------------------------------------------------------
# INLA CODE
inla.output <- inla(y~1,data=data.frame(y=y),
                    control.family=
                      list(hyper=list(prec=list(prior="loggamma",param=c(2,1)))),
                    control.fixed=list(mean.intercept=0,
                                       prec.intercept=1))
# INLA latent density grid
inla.m <- inla.output$marginals.fixed$`(Intercept)`
# INLA hyperparameters density grid
inla.sd <- inla.output$marginals.hyperpar$`Precision for the Gaussian observations`


# ---------------------------------------------------------------------
# PLOTTING

x.latent <- seq(0,1,length.out=100)

df.x <- data.frame()
i <- 1
for(theta_ in thetas)
{
  fun <- function(theta) function(x) dis.integrand(x,theta)/pos.x.norm
  y.mean.x <- sapply(x.latent,fun(theta_))
  df.x <- rbind(df.x,data.frame(x=x.latent,id=i,y=y.mean.x))
  i <- i+1
}


theme <- theme(
  panel.background = element_rect(fill = "aliceblue"),
  axis.title.x=element_text(size = 15, margin = margin(t=20)),
  axis.title.y=element_text(size = 15, angle = 90, margin = margin(r=20)),
  plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
  legend.text=element_text(size=10),
  legend.position="inside",
  legend.justification.inside=c(0.95,0.95),
  legend.background=element_rect(fill = "white")
)

pos.x.mean <- integrate(function(x) x*dens.post.x(x),0,2)$value

# Latent field plot
gx <- ggplot() + theme + xlab("$\\mathbf{x}$") + ylab("") + ggtitle("Posteriori campo latente") +
  geom_line(data=inla.m,aes(x=x,y=y,colour="Salida INLA")) + 
  geom_line(data=data.frame(x=x.latent,y=dens.post.x(x.latent)),
            aes(x=x,y=y,colour="Resultados"),linetype=1) +
  geom_line(data=df.x,aes(x=x,y=y,group=id,colour="Sumandos"))+
  geom_vline(xintercept = actual.mean,color="black",linewidth=0.8) +
  geom_vline(xintercept = inla.output$summary.fixed$mean,color="#D95F02",
             linewidth=0.7,linetype=3) +
  geom_vline(xintercept = pos.x.mean,color= "#7570B3",
             linewidth=0.7,linetype=2) +
  scale_color_manual(NULL, values = c("Salida INLA" = "#D95F02",
                                      "Resultados" = "#7570B3",
                                      "Sumandos"="#66A61E"))


pos.theta.mean <- integrate(function(theta) theta*dens.post.theta(theta),0,5)$value

# Hyperparameters plots
gtheta <-ggplot() + theme +  xlab("$\\mathbf{\\theta}$") + ylab("") + ggtitle("Posteriori hiperparámetros") +
  geom_line(data=inla.sd,aes(x=x,y=y,colour="Salida INLA")) +
  geom_line(data=data.frame(x=thetas,y=dens.post.theta(thetas)),
            aes(x=x,y=y,colour="Resultados")) +
  geom_vline(xintercept = 1/(actual.sd)^2,color="black",linewidth=0.8) +
  geom_vline(xintercept = inla.output$summary.hyperpar$mean,color="#D95F02",
             linewidth=0.7,linetype=3) +
  geom_vline(xintercept = pos.theta.mean,color= "#7570B3",
             linewidth=0.7,linetype=2) +
  scale_color_manual(NULL, values = c("Salida INLA" = "#D95F02",
                                      "Resultados" = "#7570B3"))

# Arranged plot
ggpubr::ggarrange(gx,gtheta,ncol=2,nrow=1)

# ---------------------------------------------------------------------
