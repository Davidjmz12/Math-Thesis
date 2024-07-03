# ---------------------------------------------------------------------
# Title: gp_prediction.R
# Description: R code to implement the Gaussian Prediction example.
# Date:
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# PACKAGES
library("ggplot2")
library("RColorBrewer")
library("ggpubr")

# ---------------------------------------------------------------------

loc <- seq(-5,5,length.out=100) # Set locations
d <- function(x) exp(-1/2*x^2) # Set covariance function
M <- d(as.matrix(dist(loc))) # Set covariance matrix
diag(M) <- 1 + 1e-9 # Add to avoid singular matrices
L <- chol(M) # Compute Cholesky descomposition

X <- matrix(rnorm(6*100),ncol=6,nrow=100) # Get 6 samples in matrix


#Get prior sample with covariance matrix M
Z <- t(L) %*% X 

# ---------------------------------------------------------------------
# PLOT PRIOR

Z.melt <- melt(Z)[,-1]
colnames(Z.melt) <- c("id","y")
Z.melt$x <- loc

theme <- theme(
  panel.background = element_rect(fill = "aliceblue"),
  axis.title.x=element_text(size = 15, margin = margin(t=20)),
  axis.title.y=element_text(size = 15, angle = 90, margin = margin(r=20)),
  plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
  legend.text=element_text(size=10),
  legend.position="inside",
  legend.justification.inside=c(0.9,0.9),
  legend.background=element_rect(fill = "white")
)

# Prior plot
prior <- ggplot(data=Z.melt,aes(x=x,y=y,group=as.factor(id))) + 
  geom_line(aes(colour="Muestras priori")) +
  geom_hline(aes(colour="Media GP",yintercept=0), linewidth=1) +
  scale_color_manual(NULL, values = c("Muestras priori" = "#D95F02", "Media GP" = "black")) +
  theme + 
  ylab("$Y(\\mathbf{s})$") +
  xlab("$\\mathbf{s}$") +
  ggtitle("Priori $y^*$")


# ---------------------------------------------------------------------
# POSTERIOR

loc.train <- c(-4.55,-2.55,2.55,4.55)
y.train <- c(4,-2,-3,1)

# Function that returns the posterior of 6 random samples
posterior <- function(loc.train,y.train,loc.r)
{
  n.train <- length(loc.train)
  M.total <- d(as.matrix(dist(c(loc.train,loc.r))))
  M.train <- M.total[1:n.train,1:n.train]
  M.train.r <- M.total[-(1:n.train),1:n.train]
  M.r <- M.total[-(1:n.train),-(1:n.train)]
  
  mean <- M.train.r %*% solve(M.train) %*% y.train
  sd <- M.r - (M.train.r %*% solve(M.train)) %*% t(M.train.r)
  
  return(cbind(mean,t(MASS::mvrnorm(n=6,mu=mean,Sigma=sd))))
}


# Get posterior sampling
M <- posterior(loc.train,y.train,loc)

# ---------------------------------------------------------------------
# POSTERIOR PLOT

mean.df <- data.frame(x=loc,y=M[,1],id=1)
M <- M[,-1]
M.melt <- melt(M)[,-1]
colnames(M.melt) <- c("id","y")
M.melt$x <- loc


point.df <- data.frame(x=loc.train,y=y.train,id=1)

# Posterior plot
post <- ggplot(data=M.melt,aes(x=x,y=y,group=as.factor(id))) + 
  geom_line(aes(colour="Muestras posteriori")) +
  geom_line(data=mean.df,aes(x=x,y=y,colour="Media GP"), linewidth=1) +
  scale_color_manual(NULL, values = c("Muestras posteriori" = "#7570B3", 
                                      "Media GP" = "black",
                                      "Realizaciones GP"="red")) +
  theme + 
  ylab("$Y(\\mathbf{s})$") +
  xlab("$\\mathbf{s}$") +
  geom_point(data=point.df,aes(x=x,y=y,colour="Realizaciones GP"),size=3,shape=4,stroke=1.6) +
  ggtitle("Posteriori $y^*\\mid y$")

# Posterior and priori arranged plot
ggpubr::ggarrange(prior,post,ncol=2,nrow=1,labels = c("A","B"))

# ---------------------------------------------------------------------
