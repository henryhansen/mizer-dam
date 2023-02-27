### Libraries
library(mizer)
library(mizerExperimental)
source("_6_plotting_functions.R")

# test projection model ---------------------------------------------------

# load params
params <- readRDS("tuned_params_no_migs.rds")

# overview of params
summary(params)

# simple projection with basic params
out <- projectToSteady(params, return_sim = T)
plot(out)

# simple projection with double max recruitment
test <- params
test@species_params$R_max <- test@species_params$R_max*2
test.sim <- project(test)
plot(test.sim)

basic_plot(out,test.sim)



# # Do multiple tests with a variety of r_max values ------------------------
# 
# ### make function to test different r_pp scenarios as a step
# r_max_test <- function(base, rm1, r1t, rm2, r2t) { #given a base community, change r_max 2x
#   values <- base@params@species_params$R_max
#   base@params@species_params$R_max <- values*rm1
#   sim <- project(base, t_max = r1t, dt = 0.1)
#   sim@params@species_params$R_max <- values*rm2
#   simfinal <- project(sim, t_max = r2t,dt = 0.1)
#   return(simfinal)
# }
# 
# values <- out@species_params$R_max
# out@species_params$R_max <- values*rm1
# sim <- project(out, t_max = r1t, dt = 0.1)
# sim@species_params$R_max <- values*rm2
# simfinal <- project(sim, t_max = r2t,dt = 0.1)
# 
# 
# # project different r_max scenarios
# halfR <- r_max_test(out,.5,5,10,10) #half R_max for 5 and recover
# quarR <- r_max_test(out,0.25,5,10,10) #quar R_max for 5 and recover
# tentR <- r_max_test(out,0.2,5,10,10) #10th R_max for 5 and recover
# hadoR <- r_max_test(out,0.5,5,2,10) #half R_max for 5 and double after
# 
# 
# ### Plot Large fish indicator plots for simple tests
# basic_plot(out,halfR)
# basic_plot(out,quarR)
# basic_plot(out,tentR)
# basic_plot(out,hadoR)
# plotBiomass(hadoR)
# 
# 
# # Simulate under 100 different r_pp values from 0.001 to 10 --------------
# r_max_values <- matrix(seq(0.1,5, length.out = 100)) #create range of r_pp values
# 
# # apply test function with 5 year disturbance and 10 year recovery
# many_r_max <- apply(r_max_values, 1, r_max_test, base = out, r1t = 5, rm2 = 1, r2t = 20)
# 
# comm_plot(many_r_max,out,r_max_values)



# build function for dynamic R_max ----------------------------------------

logis <- function(l,k,x,x0) ((l)/ (1+exp(-k*(x-x0))))

l = 2
k = 0.5
x = seq(1,50,1)
x0 = 10

Ls <- list(seq(1,5,length.out = 10)) #new asymptotes that range from 1 to 5

asym <- mapply(Ls, FUN = logis, k=k,x=x,x0=x0) #matrix of asymptotes over time

cols <- brewer.pal(n = 10, name = "Spectral")
par(mar = c(5, 5, 3, 8), xpd=TRUE)
matplot(x= t(matrix(rep(x,10),nrow=10,byrow = T)),y = t(asym),col=cols,type="l", lty=1, lwd=2,)

rmax <- params@species_params$R_max

# given a range of new asymptotes, apply logis function to each rmax matrix
# here is the special function that does the calculation

changeRmax <- function(asym, rmax) {
  nrm <- list()
  for (i in 1:nrow(asym)) {
    nrm[[i]] <- matrix(rep(asym[i,],times = length(rmax)), ncol = ncol(asym), byrow = T)
    nrm[[i]] <- matrix(rep(rmax,50),ncol = 50,byrow = F) + sweep(matrix(unlist(nrm[[i]]),nrow = length(rmax),byrow = F), 1, rmax, `*`)
  }
  return(nrm)
}

dyn_rec <- changeRmax(asym, rmax)

# Dynamic recruitment test function
out <- projectToSteady(params, return_sim = T)


TimeBeverton <- function(rdi, species_params, ...)
{  
  params <- get('params',envir = parent.frame())
  t <- get('t',envir = parent.frame())
  return(rdi/(1 + rdi/unname(params@other_params$other$dyn_rec[, t+ params@other_params$other$t_idx])))
}

dyn_test_recruit <- function(base, dyn_rec) {
  
  
  dyn_rec <- array(dyn_rec,dim = c(16,50))
  
  other_params(base)$t_idx = 1
  
  other_params(base)$dyn_rec = dyn_rec
  
  base <- setRateFunction(base, "RDD", "TimeBeverton")
  ### project with function
  sim <- project(base, t_max = 50)
  return(sim)
}


###### apply test function with dynamic resources
many_rec <- sapply(dyn_rec, dyn_test_recruit, base = params)


comm_plot(many_rec,out,Ls[[1]]) #final plot is broken because lfi values are too small
