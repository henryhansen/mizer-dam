library(mizer)
library(tidyverse)
library(RColorBrewer)

# Mortality Responses -----------------------------------------------------

# Step for a period of time
## Use a constant value

# decreasing linear for a period of time
## Use seq() from a high value to a low value

# exponential decay
# a indicates y intercept
# k indicates steepness of 
# exp_decay <- function(a,k,x) a*exp(-x/(1/k)) 
# #test
# a <- 1
# k <- 0.5
# x <- seq(1,10,1)
# exp_decay(a,x,k)


# decreasing sigmoid 
#b must be negative
#a must be sufficiently postive
#k is intercept assuming a is large enough)
# dec_sig <- function(a,b,k,x) (-k/(1+exp(a+b*x)))+k
# 
# #test
# a <- 6
# b <- -2.5
# k <- 1
# x <- x
# dec_sig(a,b,k,x)


# Create Master Mortality Function ----------------------------------------

dam_mortality <- function(x, c, a, b, k, type){
  
  #exponential decay function
  exp_decay <- function(a,k,x) a*exp(-x*k) 
  
  #decreasing sigmoid function
  dec_sig <- function(a,b,k,x) (-k/(1+exp(a+b*x)))+k
  
  #if no sequence provided return constant value for mizer
  if(type == "constant"){ 
    return(c)
  }
  
  # if sequence is provided for linear relationship return linear sequence for mizer
  if(type == "linear"){
    return(seq(from = a, to = b, length = length(x)))
  }
  
  # if sequence, A and K are provided, return logistic sequence for mizer
  if(type == "logistic"){
    return(exp_decay(a,k,x))
  }
  
  if(type == "sigmoid"){
    return(dec_sig(a,b,k,x))
  }
}

##### Test if the function work for different types
# dam_mortality(c=0.5, type ="constant") #c works
# plot(dam_mortality(x = c(1:10), a = 1, b = 0, type = "linear")) # linear works
# plot(dam_mortality(x = c(1:10), a = 1, k = 0.5, type = "logistic")) # logistic works
# plot(dam_mortality(x = c(1:10), a = 6, b = -2.5, k = 1, type = "sigmoid")) # sigmoid works

# Plot example responses for logistic and sigmoid  -----------------------------
##### logisitic
# x <- seq(from = 1, to = 50, by = 1)
# k <- seq(0.1, 0.9, by=.1)
# a <- 1
# type <-  "logistic"
# 
# res <- mapply(dam_mortality, x = x, MoreArgs=list( k = k, a = a,type = type))
# 
# cols <- brewer.pal(n = 9, name = "Spectral")
# par(mar = c(5, 5, 3, 8), xpd=TRUE)
# 
# matplot(x, t(res), col=cols, type="l", lty=1, lwd=2, xlab="Time", ylab="Dam Mortality")
# legend("right",  inset = c(- 0.25, 0), legend=k, title="Value of K", lwd=2, col=cols, bty="n")


##### Sigmoid
# x <- seq(from = 1, to = 50, by = 1)
# b <- seq(-0.2, -1, by = -.1)
# k <- 1
# a <- 5
# type <- "sigmoid"
# 
# res <- mapply(dam_mortality, x = x, MoreArgs=list( b = b, k = k, a = a,  type = type))
# 
# cols <- brewer.pal(n = 9, name = "Spectral")
# par(mar = c(5, 5, 3, 8), xpd=TRUE)
# 
# matplot(x, t(res), col=cols, type="l", lty=1, lwd=2, xlab="Time", ylab="Fishing Mortality")
# legend("right",  inset = c(- 0.25, 0), legend=b, title="Value of B", lwd=2, col=cols, bty="n")

# Test dam mortality function (using the logistic function) with Mizer ----------

# # load tuned species_params
# morm <- readRDS("tuned_params.rds") 
# 
# gear_name <- "dam_mortality"
# gear_params(morm)$gear <- rep(gear_name, 17)
# times = seq(from = 1, to = 50, by = 1)
# effort_array <- array(NA, dim = c(length(times), length(gear_name)),
#                       dimnames = list(time = times, gear = gear_name))
# 
# effort_array[,"dam_mortality"] <- dam_mortality(x = times, a = 1, k = 0.1, type = "logistic")
# sim <- project(morm, effort = effort_array, dt = 0.1)
# 
# plot(sim)
# plotYield(sim)
# plotBiomass(sim)
