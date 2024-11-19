library(symengine)  # for symbolic derivatives
library(data.table)
library(dplyr)

#--------------
# R Functions for the Paper:
# "Optimal Allocation of Resources Between Control and Surveillance
# for Complex Eradication Scenarios"
#--------------


calc_moments<- function(pgf, parms) {
  # calculate first two moments (i.e mean, variance)
  require(symengine)
  e<- S(pgf)
  parms$s<- 1
  nam<- names(parms)
  args<- list()
  for(i in 1: length(parms)) {
    assign(nam[i],parms[[i]])
    args[[i]]<- S(eval(nam[i]))
  }
  args<- Vector(args)
  d1<- symengine::D(e, "s")
  d2<- symengine::D(d1, "s")
  fd1<- as.function(d1, args=args)
  fd2<- as.function(d2, args=args)
  ex<- do.call(fd1, parms)
  vx<- do.call(fd2, parms) + ex - ex^2
  return(list(Ex=ex,Vx=vx))
}

calc_pmf<- function(pgf, parms, support) {
  # Given PGF, calculate PMF for the support
  require(symengine)
  probs<- rep(NA, length(support))
  parms$s<- 0
  nam<- names(parms)
  args<- list()
  for(i in 1: length(parms)) {
    assign(nam[i],parms[[i]])
    args[[i]]<- S(eval(nam[i]))
  }
  args<- Vector(args)
  dn<- list()
  e<- symengine::S(pgf)
  dn[[1]]<- e
  fe<- as.function(e, args=args)
  probs[1]<- do.call(fe, parms)
  for(i in 2:length(support)) {
    dn[[i]]<- symengine::D(dn[[i-1]], "s")
    fe<- as.function(dn[[i]], args=args)
    ev<- do.call(fe, parms)
    lprob<- log(ev) - lfactorial(support[i])
    probs[i]<- exp(lprob)
  }
  data.frame(n=support,p=probs)
}

# Construct the various PGF's as expressions to facilitate symbolic derivatives

# kappa
kappa <- function(rho, lambda, delta, tr) {
  # calculate kappa for a given probability of absence (tr)
  t <- 1
  pgf<- make_psi_Xk(t)
  f<- as.function(S(pgf))
  P <- f(s=0, rho=rho, lambda=lambda, delta=delta)
  while (P < tr) {
    t <- t + 1
    pgf<- make_psi_Xk(t)
    f<- as.function(S(pgf))
    P <- f(s=0, rho=rho, lambda=lambda, delta=delta)
  }
  return(t)
}


# PGF helper functions

psi_Xk_num <- function(k) {
  if (k == 0) {
    num <- "exp(rho * (s - 1))"
  } else {
    num <- paste0("exp(rho * (", psi_X_num(k - 1), " - 1))")
  }
  num
}


psi_X_num<- function(k) {
  if(k==0){
    num<- "exp(lambda*((1-delta)*s-1))"
  }
  else {
    num<- paste0("exp(lambda*((1-delta)*",psi_X_num(k-1),"-1))")
  }
  num
}

psi_Xk_den <- function(k) {
  if (k == 0) {
    num <- NULL
  } else {
    num <- paste0("exp(rho * (", psi_X_den(k - 1), "))")
  }
  num
}


psi_X_den<- function(k) {
  if(k==0){
    num<- "exp(-lambda*delta)-1"
  }
  else {
    num<- paste0("exp(lambda*((1-delta)*",psi_X_den(k-1),"))-1")
  }
  num
}

make_psi_Xk <- function(k) {
  if (k == 0) {
    pgf <- psi_Xk_num(k = k)
  } else {
    num <- psi_Xk_num(k = k)
    denom <- psi_Xk_den(k = k)
    pgf <- paste0(num, "/", denom)
  }
  pgf
}


make_PND <- function(k) {
  if (k == 1) {
    num <- "exp(rho * (exp(-lambda * delta) - 1))"
  } else {
    num <- paste0("exp(rho * (", psi_X_den(k - 1), "))")
  }
  num
}



psi_D_num<- function(k) {
  if(k==0){
    num<- "exp(lambda*((1-delta)*exp(-lambda*delta)-1))"
  }
  else {
    num<- paste0("exp(lambda*((1-delta)*",psi_D_num(k-1),"-1))")
  }
  num
}

psi_Dk_num <- function(k) {
  if (k == 0) {
    num <- "exp(rho * (exp(-lambda * delta) - 1))"
  } else {
    num <- paste0("exp(rho * (", psi_D_num(k - 1), " - 1))")
  }
  num
}


psi_Dk<- function(k) {
  if(k==0) {
    val<- psi_Dk_num(k=k)
  }
  else{
    num<- psi_Dk_num(k=k)
    denom<- psi_Xk_den(k=k)
    val<- paste0(num,"/",denom)
  }
  val
}


# PK
make_Pk <- function(k) {
  if (k == 1) {
    val <- paste0("1 - ",make_PND(1))
  } else {
    val <- paste0(make_PND(k-1)," * ", "(1 - ",psi_Dk(k-1),")")
  }
  return(val)
}

# Build the various PGF's

make_pgfK <- function(k) {
  tmp <- list()
  for (i in 1:k) {
    tmp[[i]] <- paste0("(", make_Pk(i), ") * s^", i)
  }
  den <- paste0("1 - ", make_PND(k))
  num <- paste0(tmp, collapse = "+")

  return(paste0("(", num, ")/(", den, ")"))
}


make_pgfM <- function(k) {
  t1 <- make_PND(k)
  val <- paste0(t1, " / ", "(1 - (1 - ", t1, ") * s)")

  return(val)
}


make_pgfMk <- function(k) {
  t1 <- make_PND(k)
  t2<- make_pgfK(k)
  val <- paste0(t1," / ","(1 - (1 - ",t1,") * ",t2,")")
  return(val)
}

make_pgfN <- function(k) {
  t1 <- paste0("s^", k)
  t2 <- make_pgfMk(k)
  val <- paste0("(", t1, ") * (", t2, ")")

  return(val)
}



entropyfunc <- function(prob_vector) {
  # Calculate the Shannon entropy for a given probability mass function (PMF)
  prob_vector <- prob_vector[prob_vector > 0]

  entropy_value <- -sum(prob_vector * log2(prob_vector))

  return(entropy_value)
}


calculate_entropy <- function(p_vals, delta_vals, lambda_vals, tr, support) {
  # Calculate the Shannon entropy incorporating prior for model parameters using the survey block model

  grid <- expand.grid(p = p_vals, delta = delta_vals, lambda = lambda_vals)
  grid$rho <- -log(1 - grid$p)  # Calculate rho

  cat("Number of steps (grid size) required to complete the process:", nrow(grid), "\n")

  grid$kappa <- apply(grid, 1, function(row) {
    rho <- as.numeric(row["rho"])
    delta <- as.numeric(row["delta"])
    lambda <- as.numeric(row["lambda"])
    kappa_value <- kappa(rho, lambda, delta, tr)
    return(kappa_value)
  })

  num_rows <- nrow(grid)
  max_support <- max(support)

  column_names <- paste0("pmfN", 0:max_support)
  pmfN_results <- as.data.table(matrix(NA_real_, nrow = num_rows, ncol = length(column_names)))
  setnames(pmfN_results, column_names)

  for (i in 1:num_rows) {
    rho <- grid[i, "rho"]
    lambda <- grid[i, "lambda"]
    delta <- grid[i, "delta"]
    k <- grid[i, "kappa"]

    parms <- list(rho = rho, lambda = lambda, delta = delta)

    pgfN <- make_pgfN(k)

    pmfN <- calc_pmf(pgfN, parms, support = support)

    if (!is.numeric(pmfN$p) || length(pmfN$p) != length(support)) {
      stop("pmfN$p must be a numeric vector of the same length as support.")
    }

    check <- sum(pmfN$p)
    if (check < 1) {
      cat(sprintf("Warning: The exact sum of the PMF is %s.  Calculation continues with normalizing the PMF. You may consider using a larger support range.\n", check))
    }




    for (j in 1:length(pmfN$p)) {
      pmfN_results[i, j] <- pmfN$p[j] / check
    }

    cat("Step No:", i, "\n")
  }

  DistN <- colSums(pmfN_results, na.rm = TRUE) / nrow(pmfN_results)

  entVal <- entropyfunc(DistN)
  ex = sum(support *DistN)
  return(list(Ex=ex,Entropy=entVal))
}


calculate_evpi_p <- function(p_vals, delta_vals, lambda_vals, tr, cs, cm) {
  # Calculate the Expected Value of Perfect Information (EVPI)
  # when p is uncertain and delta is treated as a control variable.

  grid <- expand.grid(p = p_vals, delta = delta_vals, lambda = lambda_vals)
  grid$rho <- -log(1 - grid$p)

  grid$kappa <- apply(grid, 1, function(row) {
    rho <- as.numeric(row["rho"])
    delta <- as.numeric(row["delta"])
    lambda <- as.numeric(row["lambda"])
    kappa_value <- kappa(rho, lambda, delta, tr)
    return(kappa_value)
  })

  maxk = max(grid$kappa)
  grid <- grid %>%
    rowwise() %>%
    mutate(
      EM = {
        parms <- list(rho = rho, lambda = lambda, delta = delta)
        pgfM <- make_pgfM(maxk)
        M_moments <- calc_moments(pgfM, parms)
        M_moments$Ex
      },
      EN = {
        parms <- list(rho = rho, lambda = lambda, delta = delta)
        pgfN <- make_pgfN(maxk)
        N_moments <- calc_moments(pgfN, parms)
        N_moments$Ex
      },
      Utility = cs * EN + cm * EM
    ) %>%
    ungroup()

  dataTable <- as.data.frame(grid)

  utilityMinByP <- dataTable %>%
    group_by(p) %>%
    summarise(MinUtility = min(Utility, na.rm = TRUE), .groups = 'drop')

  EUtilityCurrent <- sum(utilityMinByP$MinUtility, na.rm = TRUE) / nrow(utilityMinByP)

  utilitySumByDelta <- dataTable %>%
    group_by(delta) %>%
    summarise(TotalUtility = sum(Utility, na.rm = TRUE), CountDelta = n(), .groups = 'drop') %>%
    mutate(AverageUtility = TotalUtility / CountDelta)

  EUtilityPerfect <- min(utilitySumByDelta$AverageUtility, na.rm = TRUE)

  EVPI <- EUtilityPerfect - EUtilityCurrent

  return(list(
    EVPI = EVPI
  ))
}

