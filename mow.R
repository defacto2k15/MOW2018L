library(GA)
library(cec2017)

goalFunctionId <- 4;

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

filled.contour(x1, x2, f, color.palette = bl2gr.colors)

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

  }
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