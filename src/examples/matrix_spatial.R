# ---------------------------------------------------------------------
# Title: matrix_spatial.R
# Description: R code to implement matrix spatial plots.
# Date:
# ---------------------------------------------------------------------


library("INLA")
library("ggplot2")


# ---------------------------------------------------------------------
# CODE EXTRACTED FROM :
# https://becarioprecario.bitbucket.io/spde-gitbook/index.html
s <- 3 
n <- 7
points <- cbind(c(0.1, 0.9, 1.5, 2, 2.3, 2, 1), 
             c(1, 1, 0, 0, 1.2, 1.9, 2)) * s 
mesh <- inla.mesh.2d(points[-c(3, 5, 7), ], max.edge = s * 1.7,
                     offset = s / 4, cutoff = s / 2, n = 6) 
m <- mesh$n

# END 
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# PLOTTING


# Plot mesh function
plot_mesh <- function(mesh,points){
  points.df <- data.frame(x=points[,1],y=points[,2],id=1:length(points[,2]))
  theme <- theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.title.x=element_text(size = 15, margin = margin(t=20)),
    axis.title.y=element_text(size = 15, angle = 90, margin = margin(r=20)),
    plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
    legend.text=element_text(size=10),
  )
  p1<-ggplot(data=points.df,aes(x=x,y=y)) + ggtitle("Puntos") +
    inlabru::gg(mesh) + xlab("") + ylab("") +
    theme +
    geom_point(size=6,color="lightblue",fill="white") +
    geom_text(aes(label=id))
  
  axis.df <- data.frame(x=mesh$loc[,1],y=mesh$loc[,2],id=1:length(mesh$loc[,2]))
  p2<-ggplot(data=axis.df,aes(x=x,y=y)) + ggtitle("Vértices")+
    inlabru::gg(mesh) + xlab("") + ylab("") +
    theme +
    geom_point(size=6,color="grey",fill="white") +
    geom_text(aes(label=id))
  
  ggpubr::ggarrange(p1,p2,ncol=2,nrow=1)
  
}

# Function that plots a sparse matrix X
plot_sparse_matrix <- function(X,title="",nums=T,by=1,limits.leyend=NULL,breaks.leyend=NULL)
{
  library(ggplot2)
  # Convert matrix to dataframe
  df <- reshape2::melt(as.matrix(X))
  colnames(df) <- c("row","col","value")
  dim.row <- rev(range(df$row)) + c(0.5,-0.5)
  seq.row <- seq(range(df$row)[1],range(df$row)[2],by=by)
  dim.col <- range(df$col) + c(-0.5,0.5)
  seq.col <- seq(range(df$col)[1],range(df$col)[2],by=by)
  df <- df[round(df$value,2)!=0,]
  
  if(is.null(limits.leyend) | is.null(breaks.leyend)){
    scale <- scale_fill_gradient2(low = "red", mid = "white", 
                                  high = "lightblue", midpoint = 0,name="")
  } else {
    scale <- scale_fill_gradient2(low = "darksalmon", mid = "white", 
                                  high = "lightblue", midpoint = 0,name="",
                                  breaks=breaks.leyend,limits=limits.leyend)
  }
  
  if(nums)
  {
    nums_plot <- geom_text(aes(label = sprintf("%.2f", value)), color = "black")
  } else{
    nums_plot <- NULL
  }
  
  # Create heatmap
  ggplot(df, aes(x=col, y=row)) +
    geom_tile(aes(fill = value), color = "black") +
    scale +
    nums_plot +
    scale_y_reverse(limits=dim.row,breaks=seq.row,expand=c(0,0)) + 
    scale_x_continuous(limits=dim.col,breaks=seq.col,expand=c(0,0)) +
    theme_bw() + xlab("") + ylab("") + ggtitle(title) +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "lightgrey"),
          plot.title = element_text(hjust = 0.5,size = 15, face = "bold")) + 
    coord_fixed(ratio = 1) +
    guides(
      fill = guide_colorbar(
        frame.colour = "black",  # Add a black box around the legend
        frame.linewidth = 0.6,
        ticks.colour = "black",
        ticks.linewidth = 0.6,
        barwidth = 1,  # Adjust the width of the color bar
        barheight = 10  # Adjust the height of the color bar
      )
    )
  
}

# Mesh plot
plot_mesh(mesh,points)

# Projection Matrix plot
A <- inla.spde.make.A(mesh,points)
plot_sparse_matrix(A,title="Matriz de proyección")

# Finite element representation plot
C_G <- inla.mesh.fem(mesh,order=1)
C_plot <- plot_sparse_matrix(C_G$c1,title="C")
G_plot <- plot_sparse_matrix(C_G$g1,title="G")
ggpubr::ggarrange(C_plot,G_plot)


# Function that computes the precision matrix
compute_prec <- function(mesh,tau,kappa)
{
  C_G <- inla.mesh.fem(mesh,order=1)
  C <- C_G$c0
  G <- C_G$g1
  C_inv <- solve(C_G$c0)
  return(tau^2*(kappa^4*C+2*kappa^2*G+G%*%C_inv%*%G))
}

# Precision matrix plot for different hyperparameters
p1.q <-plot_sparse_matrix(compute_prec(mesh,0.5,0.2),
                   title="$\\mathbf{Q}(\\tau=0.5,\\kappa=0.2)$")
p2.q <- plot_sparse_matrix(compute_prec(mesh,0.2,0.5),
                   title="$\\mathbf{Q}(\\tau=0.2,\\kappa=0.5)$")
ggpubr::ggarrange(p1.q,p2.q,ncol=2,nrow=1)


# Precision matrix of a denser mesh
other.mesh <- inla.mesh.2d(loc=data.frame(x=runif(1000,1000,2000),y=runif(1000,0,1000)),
                   cutoff = 150)
plot_sparse_matrix(compute_prec(other.mesh,0.5,0.2),title="Matriz de precisión",F,by=10)

# ---------------------------------------------------------------------