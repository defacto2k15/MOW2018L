
my_ga <- function (type = c("binary", "real-valued", "permutation"), fitness, 
                   ..., min, max, nBits, population = gaControl(type)$population, 
                   selection = gaControl(type)$selection, crossover = gaControl(type)$crossover, 
                   mutation = gaControl(type)$mutation, popSize = 50, pcrossover = 0.8, 
                   pmutation = 0.1, elitism = base::max(1, round(popSize * 0.05)), 
                   updatePop = FALSE, postFitness = NULL, maxiter = 100, run = maxiter, 
                   maxFitness = Inf, names = NULL, suggestions = NULL, optim = FALSE, 
                   optimArgs = list(method = "L-BFGS-B", poptim = 0.05, pressel = 0.5, 
                                    control = list(fnscale = -1, maxit = 100)), keepBest = FALSE, 
                   parallel = FALSE, monitor = if (interactive()) {
                     if (is.RStudio()) 
                       gaMonitor
                     else gaMonitor2
                   } else FALSE, seed = NULL) 
{
  call <- match.call()
  type <- match.arg(type, choices = eval(formals(ga)$type))
  if (!is.function(population)) 
    population <- get(population)
  if (!is.function(selection)) 
    selection <- get(selection)
  if (!is.function(crossover)) 
    crossover <- get(crossover)
  if (!is.function(mutation)) 
    mutation <- get(mutation)
  if (missing(fitness)) {
    stop("A fitness function must be provided")
  }
  if (!is.function(fitness)) {
    stop("A fitness function must be provided")
  }
  if (popSize < 10) {
    warning("The population size is less than 10.")
  }
  if (maxiter < 1) {
    stop("The maximum number of iterations must be at least 1.")
  }
  if (elitism > popSize) {
    stop("The elitism cannot be larger that population size.")
  }
  if (pcrossover < 0 | pcrossover > 1) {
    stop("Probability of crossover must be between 0 and 1.")
  }
  if (is.numeric(pmutation)) {
    if (pmutation < 0 | pmutation > 1) {
      stop("If numeric probability of mutation must be between 0 and 1.")
    }
    else if (!is.function(population)) {
      stop("pmutation must be a numeric value in (0,1) or a function.")
    }
  }
  if (missing(min) & missing(max) & missing(nBits)) {
    stop("A min and max range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!")
  }
  switch(type, binary = {
    nBits <- as.vector(nBits)[1]
    min <- max <- NA
    nvars <- nBits
  }, `real-valued` = {
    min <- as.vector(min)
    max <- as.vector(max)
    nBits <- NA
    if (length(min) != length(max)) {
      stop("min and max must be vector of the same length!")
    }
    nvars <- length(max)
  }, permutation = {
    min <- as.vector(min)[1]
    max <- as.vector(max)[1]
    nBits <- NA
    nvars <- length(seq(min, max))
  })
  if (is.null(suggestions)) {
    suggestions <- matrix(nrow = 0, ncol = nvars)
  }
  else {
    if (is.vector(suggestions)) {
      if (nvars > 1) 
        suggestions <- matrix(suggestions, nrow = 1)
      else suggestions <- matrix(suggestions, ncol = 1)
    }
    else {
      suggestions <- as.matrix(suggestions)
    }
    if (nvars != ncol(suggestions)) 
      stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
  }
  if (is.logical(monitor)) {
    if (monitor) 
      monitor <- gaMonitor
  }
  if (is.null(monitor)) 
    monitor <- FALSE
  if (optim) {
    optimArgs.default <- eval(formals(ga)$optimArgs)
    optimArgs.default$control[names(optimArgs$control)] <- optimArgs$control
    optimArgs$control <- NULL
    optimArgs.default[names(optimArgs)] <- optimArgs
    optimArgs <- optimArgs.default
    rm(optimArgs.default)
    if (any(optimArgs$method == c("L-BFGS-B", "Brent"))) {
      optimArgs$lower <- min
      optimArgs$upper <- max
    }
    else {
      optimArgs$lower <- -Inf
      optimArgs$upper <- Inf
    }
    optimArgs$poptim <- min(max(0, optimArgs$poptim), 1)
    optimArgs$pressel <- min(max(0, optimArgs$pressel), 1)
    optimArgs$control$maxit <- as.integer(optimArgs$control$maxit)
    if (is.null(optimArgs$control$fnscale)) 
      optimArgs$control$fnscale <- -1
    if (optimArgs$control$fnscale > 0) 
      optimArgs$control$fnscale <- -1 * optimArgs$control$fnscale
  }
  parallel <- if (is.logical(parallel)) {
    if (parallel) 
      startParallel(parallel)
    else FALSE
  }
  else {
    startParallel(parallel)
  }
  on.exit(if (parallel) parallel::stopCluster(attr(parallel, 
                                                   "cluster")))
  `%DO%` <- if (parallel && requireNamespace("doRNG", quietly = TRUE)) 
    doRNG::`%dorng%`
  else if (parallel) 
    `%dopar%`
  else `%do%`
  if (!is.null(seed)) 
    set.seed(seed)
  i. <- NULL
  fitnessSummary <- matrix(as.double(NA), nrow = maxiter, ncol = 6)
  colnames(fitnessSummary) <- names(gaSummary(rnorm(10)))
  bestSol <- if (keepBest) 
    vector(mode = "list", length = maxiter)
  else list()
  Fitness <- rep(NA, popSize)
  object <- new("ga", call = call, type = type, min = min, 
                max = max, nBits = nBits, names = if (is.null(names)) 
                  character()
                else names, popSize = popSize, iter = 0, run = 1, maxiter = maxiter, 
                suggestions = suggestions, population = matrix(), elitism = elitism, 
                pcrossover = pcrossover, pmutation = if (is.numeric(pmutation)) 
                  pmutation
                else NA, fitness = Fitness, summary = fitnessSummary, 
                bestSol = bestSol)
  Pop <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  ng <- min(nrow(suggestions), popSize)
  if (ng > 0) {
    Pop[1:ng, ] <- suggestions
  }
  if (popSize > ng) {
    Pop[(ng + 1):popSize, ] <- population(object)[1:(popSize - 
                                                       ng), ]
  }
  object@population <- Pop
  for (iter in seq_len(maxiter)) {
    object@iter <- iter
    if (!parallel) {
      # for (i in seq_len(popSize)) if (is.na(Fitness[i])) {
      #   fit <- fitness(Pop[i, ], ...)
      #   if (updatePop) 
      #     Pop[i, ] <- attributes(fit)[[1]]
      #   Fitness[i] <- fit
      # }
      Fitness <- fitness(Pop, ...);
    }
    else {
      print("TODO131 NOT IMPLEMENTED")
      Fitness <- foreach(i. = seq_len(popSize), .combine = "c") %DO% 
      {
        if (is.na(Fitness[i.])) 
          fitness(Pop[i., ], ...)
        else Fitness[i.]
      }
    }
    object@population <- Pop
    object@fitness <- Fitness
    fitnessSummary[iter, ] <- gaSummary(object@fitness)
    object@summary <- fitnessSummary
    if (is.function(monitor)) {
      monitor(object)
    }
    if (optim & (type == "real-valued")) 
      if (optimArgs$poptim > runif(1)) {
        i <- sample(1:popSize, size = 1, prob = optimProbsel(Fitness, 
                                                             q = optimArgs$pressel))
        opt <- try(suppressWarnings(optim(fn = fitness, 
                                          ..., par = Pop[i, ], method = optimArgs$method, 
                                          lower = optimArgs$lower, upper = optimArgs$upper, 
                                          control = optimArgs$control)), silent = TRUE)
        if (is.function(monitor)) {
          if (!inherits(opt, "try-error")) 
            cat("\b | Local search =", format(opt$value, 
                                              digits = getOption("digits")))
          else cat(" |", opt[1])
          cat("\n")
        }
        if (!inherits(opt, "try-error")) {
          Pop[i, ] <- opt$par
          Fitness[i] <- opt$value
        }
        object@population <- Pop
        object@fitness <- Fitness
        fitnessSummary[iter, ] <- gaSummary(object@fitness)
        object@summary <- fitnessSummary
      }
    if (keepBest) {
      object@bestSol[[iter]] <- unique(Pop[Fitness == max(Fitness, 
                                                          na.rm = TRUE), , drop = FALSE])
    }
    if (is.function(postFitness)) {
      object <- postFitness(object, ...)
      Fitness <- object@fitness
      Pop <- object@population
    }
    if (iter > 1) 
      object@run <- garun(fitnessSummary[seq(iter), 1])
    if (object@run >= run) 
      break
    if (max(Fitness, na.rm = TRUE) >= maxFitness) 
      break
    if (object@iter == maxiter) 
      break
    ord <- order(Fitness, decreasing = TRUE)
    PopSorted <- Pop[ord, , drop = FALSE]
    FitnessSorted <- Fitness[ord]
    if (is.function(selection)) {
      sel <- selection(object)
      Pop <- sel$population
      Fitness <- sel$fitness
    }
    else {
      sel <- sample(1:popSize, size = popSize, replace = TRUE)
      Pop <- object@population[sel, ]
      Fitness <- object@fitness[sel]
    }
    object@population <- Pop
    object@fitness <- Fitness
    if (is.function(crossover) & pcrossover > 0) {
      nmating <- floor(popSize/2)
      mating <- matrix(sample(1:(2 * nmating), size = (2 * 
                                                         nmating)), ncol = 2)
      for (i in seq_len(nmating)) {
        if (pcrossover > runif(1)) {
          parents <- mating[i, ]
          Crossover <- crossover(object, parents)
          Pop[parents, ] <- Crossover$children
          Fitness[parents] <- Crossover$fitness
        }
      }
      object@population <- Pop
      object@fitness <- Fitness
    }
    pm <- if (is.function(pmutation)) 
      pmutation(object)
    else pmutation
    if (is.function(mutation) & pm > 0) {
      for (i in seq_len(popSize)) {
        if (pm > runif(1)) {
          Mutation <- mutation(object, i)
          Pop[i, ] <- Mutation
          Fitness[i] <- NA
        }
      }
      object@population <- Pop
      object@fitness <- Fitness
    }
    if (elitism > 0) {
      ord <- order(object@fitness, na.last = TRUE)
      u <- which(!duplicated(PopSorted, margin = 1))
      Pop[ord[1:elitism], ] <- PopSorted[u[1:elitism], 
                                         ]
      Fitness[ord[1:elitism]] <- FitnessSorted[u[1:elitism]]
      object@population <- Pop
      object@fitness <- Fitness
    }
  }
  if (optim & (type == "real-valued")) {
    optimArgs$control$maxit <- rev(optimArgs$control$maxit)[1]
    i <- which.max(object@fitness)
    opt <- try(suppressWarnings(optim(fn = fitness, ..., 
                                      par = object@population[i, ], method = optimArgs$method, 
                                      lower = optimArgs$lower, upper = optimArgs$upper, 
                                      control = optimArgs$control)), silent = TRUE)
    if (is.function(monitor)) {
      if (!inherits(opt, "try-error")) 
        cat("\b | Final local search =", format(opt$value, 
                                                digits = getOption("digits")))
      else cat(" |", opt[1])
    }
    if (!inherits(opt, "try-error")) {
      object@population[i, ] <- opt$par
      object@fitness[i] <- opt$value
    }
  }
  object@summary <- na.exclude(object@summary)
  attr(object@summary, "na.action") <- NULL
  object@fitnessValue <- max(object@fitness, na.rm = TRUE)
  valueAt <- which(object@fitness == object@fitnessValue)
  solution <- object@population[valueAt, , drop = FALSE]
  if (nrow(solution) > 1) {
    eps <- gaControl("eps")
    solution <- unique(round(solution/eps) * eps, margin = 1)
  }
  colnames(solution) <- parNames(object)
  object@solution <- solution
  if (keepBest) 
    object@bestSol <- object@bestSol[!sapply(object@bestSol, 
                                             is.null)]
  return(object)
}



my_gaisl <-function (type = c("binary", "real-valued", "permutation"), fitness, 
                     ..., min, max, nBits, population = gaControl(type)$population, 
                     selection = gaControl(type)$selection, crossover = gaControl(type)$crossover, 
                     mutation = gaControl(type)$mutation, popSize = 100, numIslands = 4, 
                     migrationRate = 0.1, migrationInterval = 10, pcrossover = 0.8, 
                     pmutation = 0.1, elitism = base::max(1, round(popSize/numIslands * 
                                                                     0.05)), updatePop = FALSE, postFitness = NULL, maxiter = 1000, 
                     run = maxiter, maxFitness = Inf, names = NULL, suggestions = NULL, 
                     optim = FALSE, optimArgs = list(method = "L-BFGS-B", poptim = 0.05, 
                                                     pressel = 0.5, control = list(fnscale = -1, maxit = 100)), 
                     parallel = TRUE, monitor = if (interactive()) {
                       if (is.RStudio()) 
                         gaislMonitor
                       else gaislMonitor2
                     } else FALSE, seed = NULL) 
{
  call <- match.call()
  type <- match.arg(type, choices = eval(formals(gaisl)$type))
  if (!is.function(population)) 
    population <- get(population)
  if (!is.function(selection)) 
    selection <- get(selection)
  if (!is.function(crossover)) 
    crossover <- get(crossover)
  if (!is.function(mutation)) 
    mutation <- get(mutation)
  if (missing(fitness)) {
    stop("A fitness function must be provided")
  }
  if (!is.function(fitness)) {
    stop("A fitness function must be provided")
  }
  if (popSize < 10) {
    warning("The population size is less than 10.")
  }
  if (maxiter < 1) {
    stop("The maximum number of iterations must be at least 1.")
  }
  if (elitism > popSize) {
    stop("The elitism cannot be larger that population size.")
  }
  if (pcrossover < 0 | pcrossover > 1) {
    stop("Probability of crossover must be between 0 and 1.")
  }
  if (is.numeric(pmutation)) {
    if (pmutation < 0 | pmutation > 1) {
      stop("If numeric probability of mutation must be between 0 and 1.")
    }
    else if (!is.function(population)) {
      stop("pmutation must be a numeric value in (0,1) or a function.")
    }
  }
  if (missing(min) & missing(max) & missing(nBits)) {
    stop("A min and max range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!")
  }
  switch(type, binary = {
    nBits <- as.vector(nBits)[1]
    min <- max <- NA
    nvars <- nBits
  }, `real-valued` = {
    min <- as.vector(min)
    max <- as.vector(max)
    nBits <- NA
    if (length(min) != length(max)) {
      stop("min and max must be vector of the same length!")
    }
    nvars <- length(max)
  }, permutation = {
    min <- as.vector(min)[1]
    max <- as.vector(max)[1]
    nBits <- NA
    nvars <- length(seq(min, max))
  })
  if (is.null(suggestions)) {
    suggestions <- matrix(nrow = 0, ncol = nvars)
  }
  else {
    if (is.vector(suggestions)) {
      if (nvars > 1) 
        suggestions <- matrix(suggestions, nrow = 1)
      else suggestions <- matrix(suggestions, ncol = 1)
    }
    else {
      suggestions <- as.matrix(suggestions)
    }
    if (nvars != ncol(suggestions)) 
      stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
  }
  if (is.logical(monitor)) {
    if (monitor) 
      monitor <- gaislMonitor
  }
  if (is.null(monitor)) 
    monitor <- FALSE
  islSize <- max(10, floor(popSize/numIslands))
  migPop <- max(1, floor(migrationRate * islSize))
  numiter <- max(1, floor(maxiter/migrationInterval))
  parallel <- if (is.logical(parallel)) {
    if (parallel) 
      startParallel(numIslands)
    else FALSE
  }
  else {
    startParallel(parallel)
  }
  on.exit(if (parallel) parallel::stopCluster(attr(parallel, 
                                                   "cluster")))
  `%DO%` <- if (parallel && requireNamespace("doRNG", quietly = TRUE)) 
    doRNG::`%dorng%`
  else if (parallel) 
    `%dopar%`
  else `%do%`
  if (!is.null(seed)) 
    set.seed(seed)
  i. <- NULL
  object <- new("gaisl", call = call, type = type, min = min, 
                max = max, nBits = nBits, names = if (is.null(names)) 
                  character()
                else names, popSize = popSize, numIslands = numIslands, 
                migrationRate = migrationRate, migrationInterval = migrationInterval, 
                maxiter = maxiter, run = run, suggestions = suggestions, 
                elitism = elitism, pcrossover = pcrossover, pmutation = pmutation, 
                islands = list(), summary = list(), fitnessValues = list(), 
                solutions = list())
  GAs <- vector(mode = "list", length = numIslands)
  POPs <- rep(list(suggestions), times = numIslands)
  sumryStat <- rep(list(matrix(as.double(NA), nrow = numiter * 
                                 migrationInterval, ncol = 6, dimnames = list(NULL, names(gaSummary(rnorm(10)))))), 
                   numIslands)
  for (iter in seq_len(numiter)) {
    GAs <- foreach(i. = seq_len(numIslands)) %DO% {
      my_ga(type = type, fitness = fitness, ..., min = min, 
            max = max, nBits = nBits, suggestions = POPs[[i.]], 
            population = population, selection = selection, 
            crossover = crossover, mutation = mutation, popSize = islSize, 
            pcrossover = pcrossover, pmutation = pmutation, 
            elitism = elitism, updatePop = updatePop, postFitness = postFitness, 
            maxiter = migrationInterval, run = migrationInterval, 
            maxFitness = maxFitness, names = names, optim = optim, 
            optimArgs = optimArgs, keepBest = FALSE, parallel = FALSE, 
            monitor = FALSE, seed = NULL)
    }
    for (i in seq_len(numIslands)) {
      j <- seq((iter - 1) * migrationInterval + 1, iter * 
                 migrationInterval)
      sumryStat[[i]][j, ] <- GAs[[i]]@summary
      from <- i
      to <- (i%%numIslands) + 1
      j <- order(GAs[[from]]@fitness, decreasing = TRUE)[seq(migPop)]
      migrationPop <- GAs[[from]]@population[j, , drop = FALSE]
      j <- sample(setdiff(seq(GAs[[to]]@popSize), order(GAs[[to]]@fitness, 
                                                        decreasing = TRUE)[seq(elitism)]), size = migPop)
      newpop <- rbind(GAs[[to]]@population[-j, , drop = FALSE], 
                      migrationPop)
      POPs[[to]] <- newpop
    }
    object@islands <- GAs
    object@summary <- sumryStat
    if (is.function(monitor)) {
      monitor(object)
    }
    maxFitnessIslands <- sapply(sumryStat, function(x) apply(x[seq(iter * 
                                                                     migrationInterval), , drop = FALSE], 1, max))
    if (all(apply(maxFitnessIslands, 2, garun) > run)) 
      break
    if (all(maxFitnessIslands >= maxFitness)) 
      break
  }
  fitnessValues <- lapply(object@islands, function(x) x@fitnessValue)
  solutions <- lapply(object@islands, function(x) x@solution)
  object@summary <- lapply(object@summary, function(x) {
    x <- na.exclude(x)
    attr(x, "na.action") <- NULL
    x
  })
  object@fitnessValues <- fitnessValues
  object@solutions <- solutions
  object@fitnessValue <- max(unlist(fitnessValues), na.rm = TRUE)
  object@solution <- solutions[[which(object@fitnessValue == 
                                        unlist(fitnessValues))[1]]]
  return(object)
}