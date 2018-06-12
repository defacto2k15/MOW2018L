suppressPackageStartupMessages({
  library(lattice)
  library(gridExtra)
  library(latticeExtra)
  library(locfit)
  library(cec2017)
  library(lmvar)
  library(earth)
  library(kriging)
  library(sp)
  library(gstat)
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
  library(nnet)
  library(plot3D)
})

config = rev(within(list(), {
  goalFunctionId = 2
  approximation_type = "regression_splines"
  polyLevel= 4
  
  earth = rev(within(list(), {
    nprune   = 20
    nfold    = 10
  }))
  
  kriging = rev(within(list(), {
    variogram.model = "Sph"
    variogram.range = 900
    variogram.kappa = 1
  }))
  
}))
# spline regression





goalFunction <- function(x1, x2)
{
  cec2017(config$goalFunctionId, matrix(
    data = c(x1, x2),
    ncol =2,
    nrow = length(x1)
  ))
}




grid <- replicate(2, runif(1000, min=-100, max=100));
dataPoints.x <- grid[,1];
dataPoints.y <- grid[,2];
ye <- goalFunction(dataPoints.x, dataPoints.y);
panel.real = levelplot( ye~ dataPoints.x* dataPoints.y, contour=TRUE, panel = panel.2dsmoother)

if(config$approximation_type != "kriging"){
  
  polyLevel <- config$polyLevel;
  regressionRelation <- ye ~ poly(dataPoints.x,dataPoints.y,degree=polyLevel)
  
  if( config$approximation_type == "regression_splines"){
    fit <- earth(regressionRelation, pmethod="backward",nprune= config$earth$nprune, nfold=config$earth$nprune);
    grid <- replicate(2, runif(1000, min=-100, max=100));
    testPoints.x <- grid[,1];
    testPoints.y <- grid[,2];
    tinput <- poly(testPoints.x,testPoints.y,degree=polyLevel); 
    predictResults <- predict(fit, as.matrix(tinput));
  } else
    if( config$approximation_type == "polynomial"){
      fit <- lm(regressionRelation, y=TRUE, x=TRUE)
      predictResults <- predict(fit)
    }
  
  
  panel.predicted = levelplot( predictResults~ dataPoints.x* dataPoints.y, contour=TRUE, panel = panel.2dsmoother)
  grid.arrange(panel.predicted, panel.real,  ncol=2)
} else{
  # Function
  
  # dataPoints = data.frame(dataPoints.x, dataPoints.y);
  # dataPoints <- dataPoints[1:200,]
  # 
  # dataPoints$y = myf(dataPoints);
  # ylimits = c( min(dataPoints$y), max(dataPoints$y));
  # dataPoints$y = (dataPoints$y-ylimits[1])/(ylimits[2]-ylimits[1])
  
  # 
  # halves <- split(dataPoints, sample(rep(1:2, nrow(dataPoints)/2)));
  # trainset = halves[[1]];
  # testset = halves[[2]];
  
  set.seed(42)
  x <- seq(-1, 1, by = 0.01)
  y <- seq(-1, 1, by = 0.01)
  grid <- mesh(x, y)
  z    <- with(grid, 1-abs(x+y)-abs(y-x))
  #persp3D(z=z, x=x, y=y, xlab = "X", ylab = "Y", facets=TRUE, theta=10, phi=30)
  
  funct<-function(x,y) {
    goalFunction((x*100),y*100)/1000
    
  }
  
  #wartosci do trenowania sieci
  sample_x <- dataPoints.x/100
  sample_y <- dataPoints.y/100
  
  sample_z <- funct(sample_x,sample_y)
  limits = c(min(sample_z), max(sample_z));
  sample_z =  (sample_z-limits[1])/(limits[2] - limits[1])
  nn <- nnet(data.frame(dataPoints.x/100, dataPoints.y/100), sample_z, size=10, maxit = 30, linout = TRUE)
  
  test_vec = expand.grid(x, y) # gives a 40401x2 test vector
  res=predict(nn, test_vec)
  #res=funct(test_vec$Var1, test_vec$Var2)
  dim(res)<-c(201, 201)
  persp3D(z = res, x = x, y = y, facets=TRUE)
}