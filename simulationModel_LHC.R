library(deSolve)
library(tidyverse)

# set working directory
path <- "C:\\Users\\Carolin\\Documents\\Rscripts\\Covid-19\\"
setwd(path)

### the ODE model
covidmodel <- function(t, vars, parms) {
  with(as.list(c(vars, parms)), {
    foi <- 
      (beta1 * I1 +
         beta2 * I2m + beta3 * I3m +
         beta2 * b_par * I2s + b_par * beta3 * I3s +
         beta2 * m_par * Y2 + beta3 * m_par * Y3) / (popsize)
    return(
      list(
        c(
          dS = -foi * S,
          dE = foi * S - sigma * E,
          dI1 = sigma * E - gamma1 * I1,
          dI2m = gamma1 * I1 * (1 - p_par) - gamma2 * I2m,
          dI3m = gamma2 * I2m - gamma3 * I3m,
          dI2s = gamma1 * I1 * p_par - (alpha + gamma2) * I2s,
          dI3s = gamma2 * I2s - gamma3 * I3s,
          dY2 = alpha * I2s - gamma2 * Y2,
          dY3 = gamma2 * Y2 - gamma3 * Y3,
          dInf = foi * S,
          dRep = alpha * I2s
        )
      )
    )
  })
}

### the function doing the simulations
# There are default parameter values for everything, 
# so simply calling the function produces a matrix with simulated variables.
# Alternative parameters can be chosen as argument to the function, 
# e.g. calling "covidsim(alpha = 0.3)"
covidsim <- function(
  parameters = c(
    beta1 = 0.25, beta2 = 0.16, beta3 = 0.016,
    sigma = 1, gamma1 = 0.2, gamma2 = 0.14, gamma3 = 0.14,
    alpha = 0.5, p_par = 0.5, b_par = 0.8, m_par = 0.01,
    start_dist = 60, end_dist = 108, effect_dist = 1
 ),
  popsize = 60e+06, seedsize = 10, burnintime = 30,
  timewindow = 365,
  ...
) {
  altparms <- list(...)
  for(i in names(altparms)) {
   	parameters[i] <- altparms[[i]]
  }
  
  initialstate <- c(
    S = popsize, E = seedsize, I1 = 0,
    I2m = 0, I3m = 0,
    I2s = 0, I3s = 0,
    Y2 = 0, Y3 = 0,
    Infected = 0, Reported = 0
  )
  
  startparms <- c(popsize = popsize, parameters)
  startparms["m_par"] <- startparms["b_par"]
  simres <- ode(
    y = initialstate,
    times = seq(0, burnintime, 1),
    func = covidmodel,
    parms = startparms
  )
  
  nodistparms <- c(popsize = popsize, parameters)
  simres <- ode(
    y = tail(simres, 1)[, -1],
    times = seq(0, min(nodistparms[["start_dist"]], timewindow), 1),
    func = covidmodel,
    parms = nodistparms
  )
  
  if(timewindow <= nodistparms[["start_dist"]]) return(simres)
  distparms <- nodistparms
  distparms[["beta1"]] <- distparms[["beta1"]] * distparms[["effect_dist"]]
  distparms[["beta2"]] <- distparms[["beta2"]] * distparms[["effect_dist"]]
  distparms[["beta3"]] <- distparms[["beta3"]] * distparms[["effect_dist"]]
  simres <- rbind(
    head(simres, -1),
    ode(
      y = tail(simres, 1)[, -1],
      times = seq(distparms[["start_dist"]], min(distparms[["end_dist"]], timewindow), 1),
      func = covidmodel,
      parms = distparms
    )
  )
  
  if(timewindow <= nodistparms[["end_dist"]]) return(simres)
  simres <- rbind(
    head(simres, -1),
    ode(
      y = tail(simres, 1)[, -1],
      times = seq(nodistparms[["end_dist"]], timewindow, 1),
      func = covidmodel,
      parms = nodistparms
    )
  )
  
  return(simres)
}


# READ IN PARAMS FROM FILE ##########################################################

library(XLConnect)

# first drawn an incubation period and then derive
# gamma_1 and sigma under the assumption that
# the exposed time and the time in the infected compartment I1
# are each half the incubation period
# needs to be called before calculating beta
get_gamma1_sigma <- function(params, n)
{
	IP <- runif(n, params$incuPer[2], params$incuPer[3])
	gamma_1 <- 1/(0.5*IP)
	sigma <- 1/(0.5*IP)

	return(list(sigma=sigma, gamma_1=gamma_1))
}

# calculate parameters not defined in file and
# draw n samples from each parameter distribution
getParams <- function(params, n)
{
	R0 <- runif(n, params$R_0[2], params$R_0[3])
	gamma2 <- runif(n, params$gamma_2[2], params$gamma_2[3])
	gamma3 <- runif(n, params$gamma_3[2], params$gamma_3[3])
	p <- runif(n, params$p[2], params$p[3])
	alpha <- runif(n, params$alpha[2], params$alpha[3])
	b <- runif(n, params$b[2], params$b[3])
	m <- runif(n, params$m[2], params$m[3])

	q1 <- rep(params$q_1[1], n)
	q2 <- rep(params$q_2[1], n)
	q3 <- rep(params$q_3[1], n)

	temp <- get_gamma1_sigma(params, n)
	sigma <- temp[[1]]
	gamma1 <- temp[[2]]

	beta <- R0 / (1/gamma1 + (1 - p) * (1/gamma2 + 1/gamma3) + p/gamma2 * (q2 + gamma2*q3/gamma3))
	beta1 <- beta * q1
	beta2 <- beta * q2
	beta3 <- beta * q3

	parameters <- data.frame(R0=R0, beta1=beta1, beta2=beta2, beta3=beta3, 
					sigma=sigma, gamma1=gamma1, gamma2=gamma2, gamma3=gamma3,
						p_par=p, alpha=alpha, b_par=b, m_par=m)

	return(parameters)
}

getParamsCentral <- function(params, n)
{
	R0 <- params$R_0[1]
	gamma2 <- params$gamma_2[1]
	gamma3 <- params$gamma_3[1]
	p <- params$p[1]
	alpha <- params$alpha[1]
	b <- params$b[1]
	m <- params$m[1]

	q1 <- params$q_1[1]
	q2 <- params$q_2[1]
	q3 <- params$q_3[1]

	IP <- params$incuPer[1]
	gamma1 <- 1/(0.5*IP)
	sigma <- 1/(0.5*IP)

	beta <- R0 / (1/gamma1 + (1 - p) * (1/gamma2 + 1/gamma3) + p/gamma2 * (q2 + gamma2*q3/gamma3))
	beta1 <- beta * q1
	beta2 <- beta * q2
	beta3 <- beta * q3

	parameters <- data.frame(R0=R0, beta1=beta1, beta2=beta2, beta3=beta3, 
					sigma=sigma, gamma1=gamma1, gamma2=gamma2, gamma3=gamma3,
						p_par=p, alpha=alpha, b_par=b, m_par=m)

	return(parameters)
}




filename <- "params2.xlsx"

wb <- loadWorkbook(paste0(path, filename))
params.lit <- readWorksheet(wb, sheet="Sheet1")
paramnames <- params.lit$Param
params <- as.data.frame(t(params.lit[,2:4]))
names(params) <- paramnames

n <- 100
parameters <- getParams(params, n)
parameters$start_dist <- rep(60, n)
parameters$end_dist <- rep(365, n)
parameters$effect_dist <- rep(1, n)


# RUN SIMULATIONS ########################################################################


sims <- list()
for(i in 1:nrow(parameters))
{
	sims[[i]] <- as.data.frame(covidsim(parameters[i,], popsize = 60e+04, seedsize = 10, burnintime = 30, timewindow = 365))
	sims[[i]]$rep <- rep(i, nrow(sims[[i]]))
}

df <- do.call(rbind, sims)
df$rep <- as.factor(df$rep)

df$Incidence <- df$Infected - lag(df$Infected)
df$Inc.reported <- df$Reported - lag(df$Reported)
df$Incidence[is.na(df$Incidence)] <- 0
df$Incidence[df$Incidence<0] <- 0
df$Inc.reported[is.na(df$Inc.reported)] <- 0
df$Inc.reported[df$Inc.reported<0] <- 0


p <- ggplot(df, aes(x=time, y=Incidence, colour=rep)) + theme_classic() + geom_line() + 
	theme(legend.position="none") #+ ylim(0, 10000) + xlim(0, 150) 
print(p)


q <- ggplot(df, aes(x=time, y=Inc.reported, colour=rep)) + theme_classic() + geom_line() + 
	theme(legend.position="none") #+ ylim(0, 10000) + xlim(0, 150) 
print(q)


write.table(file="lhc.output.csv", df[, c("time", "Infected", "Reported", "rep")], row.names=FALSE, sep=",")






