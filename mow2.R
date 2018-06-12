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

source("C:/studiaMagisterskie/projekty/mow/git/my_gaisl.R")


config = rev(within(list(), {
  drawApproximation = T
  goalFunctionId = 6
  approximation_type = "polynomial"
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
dataPoints = data.frame(dataPoints.x, dataPoints.y);

ye <- goalFunction(dataPoints.x, dataPoints.y);
panel.real = levelplot( ye~ dataPoints.x* dataPoints.y, contour=TRUE, panel = panel.2dsmoother)


approximated_target_function_df <- NULL;
if(config$approximation_type != "kriging"){
  
  polyLevel <- config$polyLevel;
  regressionRelation <- ye ~ 1 + dataPoints.x + I(dataPoints.x^2) + dataPoints.y + I(dataPoints.y^2) + I(dataPoints.x+dataPoints.y)
  regressionRelation <- ye ~ poly(dataPoints.x,dataPoints.y,degree=polyLevel)
  
  if( config$approximation_type == "regression_splines"){
    fit <- earth(regressionRelation,pmethod="backward",nprune= config$earth$nprune, nfold=config$earth$nprune);
    predictResults <- predict(fit)
    
  } else 
    if( config$approximation_type == "polynomial"){
      fit <- lm(regressionRelation, y=TRUE, x=TRUE,data=dataPoints)
      predictResults <- predict(fit)
      approximated_target_function_df <- function(df){
         predict(fit, data.frame(dataPoints.x=df[,1], dataPoints.y=df[,2]))
      }
    }
  
  if(config$drawApproximation){
    panel.predicted = levelplot( predictResults~ dataPoints.x* dataPoints.y, contour=TRUE, panel = panel.2dsmoother)
    grid.arrange(panel.predicted, panel.real,  ncol=2)
  }
}


searchRange.x <- c(-100, 100)
searchRange.y <- c(-100, 100)

# rysowanie diagramu pokzaujaca przekrój wartosci
x1 <- seq(searchRange.x[1], searchRange.x[2], by = 1)
x2 <- seq(searchRange.y[1], searchRange.y[2], by = 1)

f <- outer(x1,
           x2,
           function(x,y){
             -goalFunction(x,y)
             #approximated_target_function_df(data.frame(x=x,y=y))
           })


# x - macierz dwuwierszowy, i n-kolumnowy gdzie w każdym wierszu jest n-współrzędnych rodziców
# y - dwuelementowy wektor z wartościami funkcji celu rodziców
# limits - macierz 2xn, ogranicza dziedzinę parametrów. W rzędzie 1 są wartości minimalne, w rzędzie 2 maksymalne
# wymaganiee, aby y[1] > y[2]
# zwraca - macierz 2xn z współrzędnymi potomków
my_crossover_internal <-function(x,y, limits, alternative){
  x1 = x[1,];
  f1 = y[1];
  x2 = x[2,]
  f2 = y[2];
  
  z = x1*x2*(x1-x2);
  a = (x2*f1 -x1*f2)/z;
  b = (x1^2*f2 - x2^2*f1)/z;
  
  xn = -b/(2*a);  
  xn = ifelse(is.na(xn), x1,xn);
  xn = pmax(pmin(xn, limits[2,]), limits[1,])
  
  alpha = runif(1);
  
  xm = x2 + alpha*(x1-x2)*((f1-f2)/f1)
  xm = pmax(pmin(xm, limits[2,]), limits[1,])
  
  children <- matrix(as.double(NA), nrow = 2, ncol = length(x1))
  
  children[1,] <- xn
  children[2,] <- xm
  
  children;
}

ga_spCrossover_R <- function(object, parents)
{
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  a <- runif(1)
  children[1,] <- a*parents[1,] + (1-a)*parents[2,]
  children[2,] <- a*parents[2,] + (1-a)*parents[1,]
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

ga_spCrossover_based <- function(x,y)
{
  n <- ncol(x)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  a <- runif(1)
  return(a*x[1,] + (1-a)*x[2,])
}




# points zawierac bedzie historie znajdywanych punktów
points <- list()

ga_islandCount <- 6;

for(i in 1:ga_islandCount ){ 
  points[[i]] <- list();
}

ga_maxiter <- 100
ga_migrationInterval <- 10;
ga_epochCount <- ga_maxiter/ga_migrationInterval;

current_iteration <- 0;
GA <- my_gaisl(
  type = "real-valued",
  function(x){
      -approximated_target_function_df(x)
      @-goalFunction(x[,1], x[,2])
  },
  min = c(searchRange.x[1], searchRange.y[1]),
  max = c(searchRange.x[2], searchRange.y[2]),
  popSize = 50,
  maxiter = ga_maxiter,
  migrationInterval = ga_migrationInterval,
  run = 100,
  numIslands = ga_islandCount,
  monitor = function(x) {
    ru <<- x;
    current_iteration <<- current_iteration+1;
    
    xs <- vector();
    ys <- vector();
    
    for(i in 1:ga_islandCount ){ 
      xs <-  ru@islands[[i]]@population[,1]
      ys <-  ru@islands[[i]]@population[,2] 
      points[[i]][[current_iteration]] <<- rbind (xs, ys)
      
    }
    gaislMonitor(x)
    
  },
  #crossover = my_crossover,
  parallel = F
)

final_epoch_count <-  nrow(GA@summary[[1]])/GA@migrationInterval;

intensity_colors <- list();
for(island in 1:ga_islandCount){
  max_col <- (ga_islandCount*2)+2;
  intensity_colors[[island]] <- rainbow(n=final_epoch_count, start = (2*island+1)/max_col, end=(2*island+2)/max_col);
}
terrain.colors(final_epoch_count)

filled.contour(x1, x2, f, color.palette = bl2gr.colors,
               plot.axes = {
                 axis(1)
                 axis(2)

                 for(island in 1:ga_islandCount){

                   for (iter in 1:final_epoch_count) {
                     x <- points[[island]][[iter]][1, ] + island;
                     y <- points[[island]][[iter]][2, ] + island;
                     points(x, y, col = intensity_colors[[island]][[iter]])
                   }

                   points(
                     GA@islands[[island]]@solution[, 1],
                     GA@islands[[island]]@solution[, 2],
                     pch = 3,
                     cex = 2,
                     col = intensity_colors[[island]][[1]],
                     lwd = 2
                   )
                 }
               })