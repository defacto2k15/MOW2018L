library(GA)
library(cec2017)

goalFunctionId <- 6;

searchRange.x <- c(-100, 100)
searchRange.y <- c(-100, 100)

# funkcja celu
goalFunction <- function(x1, x2)
{
  cec2017(goalFunctionId, matrix(
    data = c(x1, x2),
    ncol = 2,
    nrow = length(x1)
  ))
}

# rysowanie diagramu pokzaujaca przekrój wartosci
x1 <- seq(searchRange.x[1], searchRange.x[2], by = 1)
x2 <- seq(searchRange.y[1], searchRange.y[2], by = 1)

f <- outer(x1,
           x2,
           goalFunction)

# filled.contour(x1, x2, f, color.palette = bl2gr.colors)

# x - macierz dwuwierszowy, i n-kolumnowy gdzie w kazdym wierszu jest n-wspólrzednych rodziców
# y - dwuelementowy wektor z wartosciami funkcji celu rodziców
# limits - macierz 2xn, ogranicza dziedzine parametrów. W rzedzie 1 sa wartosci minimalne, w rzedzie 2 maksymalne
# wymaganiee, aby y[1] > y[2]
# zwraca - macierz 2xn z wspólrzednymi potomków
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

my_crossover <- function(object, parents)
{
  c1 = ga_spCrossover_R(object, parents);
  
  fitness <- object@fitness[parents];
  parents <- object@population[parents,,drop = FALSE]
  if(fitness[1] < fitness[2]){
    t <- fitness[2];
    fitness[2] <- fitness[1];
    fitness[1] <- t;
    
    t <- parents[2,];
    parents[2,] <- parents[1,];
    parents[1,] <- t;
  }

  limits = matrix(c(-100,100,-100,100), nrow=2, ncol=2)
  
  children <- my_crossover_internal(parents, fitness, limits, ga_spCrossover_based);
  
  out <- list(children = (children), fitness = rep(NA,2))
  return(out)
}



B = matrix(c(1,1, 10, 10), nrow=2, ncol=2)
Y = c(10000000,1)
#q <- my_crossover_internal(B,Y)


# points zawierac bedzie historie znajdywanych punktów
points <- list()

ga_maxiter <- 200

current_iteration <- 0;
GA <- ga(
  type = "real-valued",
  fitness =  function(x, a, b){
      if( current_iteration %% 10 ){
        -1000;
      }else{
     -goalFunction(x[1], x[2])
      }
    },
  min = c(searchRange.x[1], searchRange.y[1]),
  max = c(searchRange.x[2], searchRange.y[2]),
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
  crossover = my_crossover
)


intensity_colors <- terrain.colors(ga_maxiter)

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