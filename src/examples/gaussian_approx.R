# ---------------------------------------------------------------------
# Title: gaussian_approx.R
# Description: R code to implement the Gaussian approximation example.
# Date:
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# PACKAGES
library("ggplot2")

# ---------------------------------------------------------------------

y<-3 # Observed data

# Density function
dens <- function(x)
{
  exp(-0.5*x^2+x*y-exp(x))
}

# Grid of x 
x <- seq(-3,3,length.out=1000)
x.y <- dens(x)

# Get mode of density
mode <- x[which.max(x.y)]


# ---------------------------------------------------------------------
# Gaussian approximation

# Initial mode guess
initial.mode <- -2

# Array of modes
mode.it <- c(initial.mode)
# Array of variances
Q.it <- c(1)
n_it <- 100 # Num iterations

# Function to get taylor coefficients
f.taylor <- function(x) {x*y-exp(x)}

for(i in 1:n_it)
{
  # Taylor coefficients of 1,x,x^2
  tay_cov <- rev(taylor(f.taylor,mode.it[i],n=2))
  # Addapt to our needs
  tay_cov[3] <- -2*tay_cov[3]
  
  # Compute new values
  Q.it <- c(Q.it,Q.it[1]+tay_cov[3])
  new.mode <- solve(Q.it[1]+tay_cov[3],tay_cov[2])
  mode.it <- c(mode.it,new.mode)
}

# ---------------------------------------------------------------------
# PLOTTING

gaussian.df <- function(i)
{
  mu <- mode.it[i]
  sigma <- 1/sqrt(Q.it[i])
  return(data.frame(x,y=dnorm(x,mu,sigma)))
}

theme <- theme(
  panel.background = element_rect(fill = "aliceblue"),
  axis.title.x=element_text(size = 15, margin = margin(t=20)),
  axis.title.y=element_text(size = 15, angle = 90, margin = margin(r=20)),
  plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
  legend.text=element_text(size=10),
  legend.position="inside",
  legend.justification.inside=c(0.1,0.9),
  legend.background=element_rect(fill = "white")
)

# Resulting plot
ggplot() + 
  theme +
  xlab("x") + ylab("") + ggtitle("Ejemplo aproximación Gaussiana") +
  geom_line(data=data.frame(x,x.y),aes(x=x,y=x.y,colour="Densidad exacta")) +
  geom_point(data=data.frame(x=mode.it,y=rep(0,n_it+1)),
             aes(x=x,y=y,colour="Modas aproximadas"),shape=4) +
  geom_point(data=data.frame(x=mode,y=0),aes(x=x,y=y,colour="Moda"),size=2,shape=3) +
  geom_line(data=gaussian.df(1),aes(x=x,y=y,colour="Iteración 1"),
            linetype="dashed") +
  geom_line(data=gaussian.df(n_it),aes(x=x,y=y,colour="Iteración 100"),
            linetype="dashed") +
  scale_color_manual(NULL, values = c("Densidad exacta" = "#7570B3",
                                      "Modas aproximadas" = "red",
                                      "Moda" = "black",
                                      "Iteración 1"="#D95F02",
                                      "Iteración 100"="#66A61E"))

# ---------------------------------------------------------------------

