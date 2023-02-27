library(mizer)
library(tidyverse)
library(RColorBrewer)

# Resource Response -----------------------------------------------------------

# rebound response function
# m indicates new stable state
# u indicates rate of recovery
# b constrains system maximum and minimum
rebound_resource <- function(m,x,u,b) m-exp(-1*x)*(sin(u*x+b)/sin(b))
#test
m <- 10
u <- -1.5
b <- -.1
x <- seq(0.01,10, length.out = 1000)
plot(x,rebound_resource(m,x,u,b))

# Visualize different resource dynamics -----------------------------------

x <- seq(0.01,10, length.out = 1000)
b <- seq(-.15, -0.05, length.out = 10)

res <- mapply(rebound_resource, x = x, MoreArgs=list( b = b, u = u,m = m))

cols <- brewer.pal(n = 10, name = "Spectral")
par(mar = c(5, 5, 3, 8), xpd=TRUE)

matplot(x, t(res), col=cols, type="l", lty=1, lwd=2, xlab="Time", ylab="Resource Rate")
legend("right",  inset = c(- 0.25, 0), legend=round(b,2), title="Value of b", lwd=2, col=cols, bty="n")

# create time indexed resource dynamics function for Mizer--------------------

#Read in params
params <- readRDS("mizer_output/final_params_exp.RDS")

# Create parameter objects
params_test <- params

other_params(params_test)$t_idx = 1

rebound_forcing <- function (params, n, n_pp, n_other, rates, t, dt, resource_rate, 
          resource_capacity, steps, ...) 
{
  resource_rate <- params@other_params$other$r_pp_array[t + params@other_params$other$t_idx]* params@w_full^(resource_params(params)[["n"]] - 1)
  mur <- resource_rate + rates$resource_mort
  n_steady <- resource_rate * resource_capacity/mur
  n_pp_new <- n_steady + (n_pp - n_steady) * exp(-mur * dt)
  sel <- mur == 0
  n_pp_new[sel] <- n_pp[sel]
  n_pp_new
}

# Attach the new function
params_test@resource_dynamics <- "rebound_forcing"



# Create r_pp array using rebound function
b <- -0.05
x <- seq(0.01,4, length.out = 50)
rebound_rpp <- rebound_resource(m,x,u,b)
rebound_rpp
# Attach resource spectra
other_params(params_test)$r_pp_array <- rebound_rpp

### project with function
test <- project(params_test, t_max = 50)

plotBiomass(test)
# mizer::animateSpectra(test,ylim = c(10e-5,10e10))
mizer::plot(test)

# plot community slope
source("_6_plotting_functions.R")

out <- project(params)

basic_plot(test,out)
