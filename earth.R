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
})

config = rev(within(list(), {
  goalFunctionId = 9
  approximation_type = "kriging"
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
  regressionRelation <- ye ~ polym(dataPoints.x,dataPoints.y,degree=polyLevel)
  
  if( config$approximation_type == "regression_splines"){
    fit <- earth(regressionRelation,pmethod="backward",nprune= config$earth$nprune, nfold=config$earth$nprune);
    predictResults <- predict(fit)
  } else
  if( config$approximation_type == "polynomial"){
    fit <- lm(regressionRelation, y=TRUE, x=TRUE)
    predictResults <- predict(fit)
  }
  
  panel.predicted = levelplot( predictResults~ dataPoints.x* dataPoints.y, contour=TRUE, panel = panel.2dsmoother)
  grid.arrange(panel.predicted, panel.real,  ncol=2)
} else{
  r <- data.frame( x= dataPoints.x, y = dataPoints.y, value=ye)
  coordinates(r) <- ~ x + y
  
  grid <- expand.grid(x=seq(-100, 100, by=10), y=seq(-100, 100, by=10))
  coordinates(grid) <- ~ x + y
  
  lzn.vgm <- variogram((value)~1, r) # calculates sample variogram values 
  lzn.fit <- fit.variogram(lzn.vgm, model=vgm(
    model = config$kriging$variogram.model, 
    range = config$kriging$variogram.range,
    kappa = config$kriging$variogram.kappa),
    debug.level=0) # fit model
  
  lzn.kriged <- krige(value ~ 1, r, grid, model=lzn.fit)
  
  p <- lzn.kriged %>% as.data.frame %>%
    ggplot(aes(x=x, y=y)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
    scale_fill_gradient(low = "yellow", high="red") +
    scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
    theme_bw()
  
  
  pp <- levelplot( lzn.kriged$var1.pred~ lzn.kriged$x* lzn.kriged$y, contour=TRUE, panel = panel.2dsmoother)
  grid.arrange(pp, panel.real, ncol=2)  
}