# Respirometry models for Allen-Waller & Barott 2023 Symbiosis
# Luella Allen-Waller
# Adapted from Teegan Innis (Innis et al. 2021 GCB)

setwd("~/Box Sync/Barott lab/Data/Lulu Data/Monty pHithon and the Holey Pale/MpHi respirometry/PI model outputs")
library(ggplot2)
library(rjags)
library(snowfall)
library(dplyr)
library(coda)
library(reshape)
library(reshape2)
library(scales)
library(zoo)
library(pipeR)
library(lubridate)
library(readr)
library(tidyverse)
library(ggpubr)
library(patchwork)

# Clear environment
rm(list = ls())
set.seed(3593)

# Import data
Series <- do.call(rbind, lapply(list.files(pattern = "*Raw.csv"), read.csv))

## Subset by timepoint
June <- subset(Series, timepoint == "june")
July <- subset(Series, timepoint == "july")

# Set up the model
## BUGS code
bugs.model <- "
model {
  ## Likelihood
  for(k in 1:N) {
    P[k] ~ dnorm(mu[k],tau)
    mu[k] <- Pmax*(1-exp(-a*I[k]/Pmax))*exp(-b*I[k]/Pmax)+R
  }

  ## Priors
  tau ~ dgamma(100,1)
  Pmax ~ dnorm(0,1)T(0,)
  a ~ dnorm(0,10)T(0,)
  b ~ dnorm(0,10)T(0,)
  R ~ dnorm(0,1)

  ## Prediction
  for(k in 1:N.pr) {
    mu.pr[k] <- Pmax*(1-exp(-a*I.pr[k]/Pmax))*exp(-b*I.pr[k]/Pmax)+R
  }
  
  ## Transformations
  sigma <- 1/sqrt(tau)
  alpha <- a/100
  beta <- b/100
}
"


## JAGS
n.cores <- parallel::detectCores()-1
jags.snowfall <- function(file,data,inits,vars,n.samp,n.thin=1,
                          n.cores=4,n.chains=2,
                          n.adapt=1000,n.update=1000,
                          modules=c("glm")) {
  stopifnot(require(snowfall),require(rjags))
  seeds <- sample(1:1e6, n.cores)
  if (inherits(file, "connection")) {
    model.code <- readLines(file, warn = FALSE)
    file <- tempfile()
    writeLines(model.code, file)
  }
  sfInit(parallel = TRUE, cpus = n.cores)
  sfLibrary(rjags)
  sfExport("file","data","inits","vars","n.thin","n.samp","n.adapt","n.update","modules","seeds")
  on.exit(sfStop())
  s <- sfLapply(1:n.cores, function(i) {
    set.seed(seeds[i])
    inits <- lapply(sample(1:1e6, n.chains),
                    function(s) c(inits,list(.RNG.name="base::Wichmann-Hill",.RNG.seed=s)))
    for(m in modules) load.module(m)
    model <- jags.model(file=file,data=data,inits=inits,
                        n.chains=n.chains,n.adapt=n.adapt)
    update(model,n.update)
    jags.samples(model,var=vars,thin=n.thin,n.iter=n.samp*n.thin)
  })
  lapply(setNames(vars,vars),
         function(v) {
           a <- lapply(s,`[[`,v)
           dm <- dim(a[[1]])
           dm[length(dm)] <- length(a)*dm[length(dm)]
           structure(do.call(c,a),
                     dim=dm,
                     varname=attr(a[[1]],"varname"),
                     class="mcarray")
         })
}

## ---------------------------------------------------------------------------------
# June

# Plot data
ggplot(June, aes(x = rad, y = coefficients.x1, group = ID)) +
  geom_path() + ylab("O2 evolution") + xlab("irradiance") +
  facet_wrap(~ ID)

## Split into separate datasets
dsJune <- split(June, June$ID)

# Set parameters
I.pr <- seq(0, 800, 10)

# Run the model
### lapply
model.june <- lapply(dsJune, function(x) {
  jags.snowfall(textConnection(bugs.model),
                data = list(P = x$coefficients.x1,
                            I = x$rad/100,
                            N = nrow(x),
                            I.pr = I.pr/100,
                            N.pr = length(I.pr)),
                n.chains = 4,
                inits = list(a = 0.1, b = 0.1, Pmax = 1, R = -0.4),
                vars = c("alpha", "beta", "Pmax_µmol_L_min", "R_µmol_L_min", "mu", "sigma", "mu.pr"),
                n.samp = 2000, n.thin = 20)
}) 

# Convert to coda format
coda.june <- lapply(model.june, function(x) {
  mcmc.list(lapply(1:dim(x$alpha)[3], function(ch) mcmc(sapply(x[c("alpha","beta","Pmax_µmol_L_min","R_µmol_L_min")], function(x) x[,,ch]))))
})

# Convert to dataframe
june.summary <- lapply(coda.june, function(x) {
  t(as.data.frame(rbind(colMeans(as.data.frame(as.matrix(x)))))) 
})
summary.first <- as.data.frame(june.summary)
summary.first$variable <- c("alpha", "beta", "Pmax_µmol_L_min", "R_µmol_L_min")
summary.first <- summary.first %>%
  select(variable, everything())
summary.first <- as.data.frame(t(summary.first))
summary.first <- summary.first[-1,]
summary.first$frag <- gsub("X", "", rownames(summary.first))
summary.first$timepoint <- "June"

write.csv(summary.first, "June Resp Mod Summ.csv")

## Posterior predicted mean and credible intervals plots
lapply(dsJune, function(x) {
  plot(coefficients.x1 ~ rad, data = x, pch = 16)
  matlines(I.pr, t(lapply(model.one$mu.pr, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.975))))), lty = 1)
})

########### July

# Plot data
ggplot(July, aes(x = rad, y = coefficients.x1, group = ID)) +
  geom_line() + ylab("O2 evolution") + xlab("irradiance") +
  facet_wrap(~ ID)

## Split into separate datasets
dsJuly <- split(July, July$ID)

# Set parameters
I.pr <- seq(0, 800, 10)

# Run the model
### lapply
model.july <- lapply(dsJuly, function(x) {
  jags.snowfall(textConnection(bugs.model),
                data = list(P = x$coefficients.x1,
                            I = x$rad/100,
                            N = nrow(x),
                            I.pr = I.pr/100,
                            N.pr = length(I.pr)),
                n.chains = 4,
                inits = list(a = 0.1, b = 0.1, Pmax = 1, R = -0.4),
                vars = c("alpha", "beta", "Pmax_µmol_L_min", "R_µmol_L_min", "mu", "sigma", "mu.pr"),
                n.samp = 2000, n.thin = 20)
}) 

# Convert to coda format
coda.july <- lapply(model.july, function(x) {
  mcmc.list(lapply(1:dim(x$alpha)[3], function(ch) mcmc(sapply(x[c("alpha","beta","Pmax_µmol_L_min","R_µmol_L_min")], function(x) x[,,ch]))))
})

# Convert to dataframe
july.summary <- lapply(coda.july, function(x) {
  t(as.data.frame(rbind(colMeans(as.data.frame(as.matrix(x)))))) 
})
summary.second <- as.data.frame(july.summary)
summary.second$variable <- c("alpha", "beta", "Pmax_µmol_L_min", "R_µmol_L_min")
summary.second <- summary.second %>%
  select(variable, everything())
summary.second <- as.data.frame(t(summary.second))
summary.second <- summary.second[-1,]
summary.second$frag <- gsub("X", "", rownames(summary.second))
summary.second$timepoint <- "July"

write.csv(summary.second, "July Resp Mod Summ.csv")

## Posterior predicted mean and credible intervals plots
lapply(dsJuly, function(x) {
  plot(coefficients.x1 ~ rad, data = x, pch = 16)
  matlines(I.pr, t(lapply(model.july$mu.pr, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.975))))), lty = 1)
})

#####---------------------------------------------------------------------
##### Summarize Model #####
coda.summary <- as.data.frame(rbind(summary.first, summary.second))
write.csv(coda.summary, "McapRespModel_summary_output.csv")