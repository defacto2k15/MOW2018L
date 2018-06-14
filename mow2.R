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
  goalFunctionId = 4
  approximation_type = "nn"
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
  
  approximation_integration = rev(within(list(), {
    use_approximation_fit_percent = T;
    appriximation_fit_percent = 1;
    
    use_approximation_population_percent = F;
    appriximation_population_percent = 0.8;
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

goalFunction_dataFrame <- function(x)
{
  cec2017(config$goalFunctionId, matrix(
    data = as.vector(x),
    ncol =ncol(x),
    nrow = nrow(x)
  ))
}

searchRange <- c(-100, 100)


set.seed(42)
grid <- replicate(2, runif(1000, min=-100, max=100));
dataPoints.x <- grid[,1];
dataPoints.y <- grid[,2];
dataPoints = data.frame(dataPoints.x, dataPoints.y);
dataPoints.mat = matrix(c(dataPoints.x, dataPoints.y), ncol=2, nrow=length(dataPoints.x));

ye <- goalFunction_dataFrame(dataPoints.mat);


approximated_target_function_df <- NULL;
targetFunction.real = function(x){
  goalFunction_dataFrame(x);
};
targetFunction.approximated = NULL;

  
polyLevel <- config$polyLevel;
polied <- poly(dataPoints.mat,degree=polyLevel, raw=T);
dfx <- as.data.frame(polied);
regressionRelation <- ye ~ polied;

if( config$approximation_type == "regression_splines"){
  fit <- earth(ye~., data =dfx, pmethod="backward",nprune= config$earth$nprune, nfold=config$earth$nprune);
  predictResults <- predict(fit)
  targetFunction.approximated <- function(df){
    testInput <- poly(df, degree=polyLevel, raw=T);
    as.vector(predict(fit, as.vector(testInput)))
  }
  
} else 
  if( config$approximation_type == "polynomial"){
    fit <- lm(ye~., y=TRUE, x=TRUE,data=dfx)
    predictResults <- predict(fit)
    targetFunction.approximated <- function(df){
      testInput <- poly(df, degree=polyLevel, raw=T);
      predict(fit, testInput)
    }
  }else
    if(config$approximation_type == "nn"){
      
      #wartosci do trenowania sieci
      sample_x <- dataPoints.x/100
      sample_y <- dataPoints.y/100
      #sample_arguments <- 
      
      sample_z <- ye;
      limits = c(min(sample_z), max(sample_z));
      sample_z =  (sample_z-limits[1])/(limits[2] - limits[1])
      nn <- nnet(dataPoints.mat/100, sample_z, size=10, maxit = 30, linout = TRUE)
      
      test_vec = grid/100 # gives a 40401x2 test vector
      predictResults=predict(nn, test_vec)
      targetFunction.approximated <- function(df){
        predict(nn, (df/100))*1000
      }
    }

if(config$drawApproximation){
  grid3 <- expand.grid(seq(searchRange[1], searchRange[2], length.out = 30),
                   seq(searchRange[1], searchRange[2], length.out = 30))
  
  dataPoints3.x <- grid3[,1];
  dataPoints3.y <- grid3[,2];
  pr <- targetFunction.approximated(matrix(c(dataPoints3.x, dataPoints3.y), ncol=2));
  panel.predicted = levelplot( pr~ dataPoints3.x* dataPoints3.y, contour=TRUE, panel = panel.2dsmoother)

  real <- targetFunction.real(matrix(c(dataPoints3.x, dataPoints3.y), ncol=2))
  panel.real = levelplot( real~ dataPoints3.x* dataPoints3.y, contour=TRUE, panel = panel.2dsmoother)
  
  grid.arrange(panel.predicted, panel.real,  ncol=2)
}





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
GA <- my_ga(
  type = "real-valued",
  fitness =  function(x,a){
    if(config$approximation_integration$use_approximation_population_percent){
      if( runif(1) < config$approximation_integration$appriximation_population_percent ){
        return(as.vector(-targetFunction.approximated(x)));
      }else{
        return(-targetFunction.real(x));
      }
    }else if(config$approximation_integration$use_approximation_fit_percent){
      perc <- config$approximation_integration$appriximation_fit_percent;
      splitter <- sample(c(0,1), replace=TRUE, size=nrow(x), prob=c(perc, 1-perc))

      if(length(which(splitter==0)) == 0){
        return -(targetFunction.real(x));
      }
      if(length(which(splitter==1)) == 0){
        return (as.vector(-(targetFunction.approximated(x))));
      }

      splitedPoints <- split(x,splitter);

      results <- list();
      results$`0` <- -targetFunction.approximated(matrix(splitedPoints$`0`, nrow=sum(splitter==0), ncol=2));
      results$`1` <- -targetFunction.real(matrix(splitedPoints$`1`, nrow=sum(splitter==1), ncol=2));


      return(as.vector(unsplit(results, splitter)));
    }else{
      -targetFunction.real(x);
    }
  },
  min = c(searchRange[1], searchRange[1]),
  max = c(searchRange[2], searchRange[2]),
  popSize = 50,
  maxiter = ga_maxiter,
  run = 100,
  monitor = function(x) {
    current_iteration <<- x@iter;
    xs <- x@population[, 1]

    ys <- x@population[, 2]

    points[[x@iter]] <<- rbind (xs, ys)
    gaMonitor(x)

  },
  #crossover = my_crossover_internal
)


intensity_colors <- terrain.colors(ga_maxiter)

# rysowanie diagramu pokzaujaca przekrój wartosci
x1 <- seq(searchRange[1], searchRange[2], by = 1)
x2 <- seq(searchRange[1], searchRange[2], by = 1)

f <- matrix(-targetFunction.approximated(as.matrix(expand.grid(x1,x2))), nrow=length(x1));


filled.contour(x1, x2, f, color.palette = bl2gr.colors,
               plot.axes = {
                 axis(1)
                 axis(2)

                 for (iter in 1:GA@iter) {
                   points(points[[iter]][1, ], points[[iter]][2, ], col = intensity_colors[iter])
                 }

                 points(
                   GA@solution[, 1],
                   GA@solution[, 2],
                   pch = 3,
                   cex = 2,
                   col = "white",
                   lwd = 2
                 )
               })