library(mizer)

# sta <- list(type ="static",
#      dam_mortality =
#        list(
#          e1 = 0,
#          t1 = 1,
#          e2 = 0,
#          t2 = 49
#        ),
#      resource_rate =
#        list(
#          r1 = 10,
#          t1 = 5,
#          r2 = 10,
#          t2 = 45
#        )
#        ,
#      recruitment_max = 
#        list(
#          rm1 = 0,
#          t1 = 20,
#          rm2 = 0,
#          t2 = 30
#        )
#       )
# 
# dyn <- list(type ="dynamic",
#      dam_mortality =
#        list(
#        a = 1,
#        k = 0.9),
#      resource_rate =
#        list(
#        m = 10,
#        u = -1.5,
#        b = -0.05,
#        start = 0.01,
#        finish = 5),
#      recruitment_max = list(
#        l = 0.001,
#        k = 0.001,
#        x0 = 20
#      ))
# 
# scenario_plot(dyn)
# 
# params <- readRDS("final_params_exp.RDS")

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

TimeBeverton <- function(rdi, species_params, ...) {
  t = get('t',envir = parent.frame())
  params = get('params',envir = parent.frame())
  new_rmax <- params@other_params$other$log_rec[,as.integer(t)]
  # print(list(c(as.integer(t)+1)))
  return(rdi/(1 + rdi/new_rmax))
}



project_dam_removal <- function(params, mechanism_params = NULL, t_max = 50, ...) {
  require(mizer)
  out <- project(params, return_sim = T, t_max = 1)
  
  if (is.null(mechanism_params)) {
    warning("Because you did not input any parameters for dam removal, a mizer::project will be used")
    simulation <- project(params, ...)
  } 
  
  if(mechanism_params$type == "static") {

      #project with dam removal functions
      
      if(is.list(mechanism_params$dam_mortality)) {
        spn <- nrow(species_params(params))
        
        gear_params(params) <- data.frame(gear=c(rep("Dam_Mortality", spn), rep("Migration", spn)),
                                          species = rep(species_params(params)$species,2),
                                          catchability = c(rep(1,spn), rep(0, spn)),
                                          knife_edge_size = c(rep(0.1, 18), gear_params(params)$knife_edge_size))
        
        gear_params(params)$catchability[c(22,36)] <- 1 #replace catchability for migrants back to 1
        
        times = seq(from = 1, to = t_max, by = 1)
        effort_array <- array(NA, dim = c(length(times), 2),
                              dimnames = list(time = times, gear = unique(gear_params(params)$gear)))
        effort_array[,"Dam_Mortality"] <- c(rep(as.double(mechanism_params$dam_mortality["e1"]),
                                                as.double(mechanism_params$dam_mortality["t1"])),
                                            rep(as.double(mechanism_params$dam_mortality["e2"]),
                                                as.double(mechanism_params$dam_mortality["t2"])))
        effort_array[,"Migration"] <- rep(3,t_max)
        # print(effort_array)

      }
    
      if(is.list(mechanism_params$resource_rate)) {
        # #calculate resulting function
        # #add resource array to params and declare extra variables
        ### r_pp changes rebounds using rebound function
        # Create parameter objects
        other_params(params)$t_idx = 1

        # Attach the new function
        params@resource_dynamics <- "rebound_forcing"

        # Create r_pp array
        step_rpp <- unlist(c(rep(mechanism_params$resource_rate["r1"],
                          mechanism_params$resource_rate["t1"]),
                      rep(mechanism_params$resource_rate["r2"],
                          mechanism_params$resource_rate["t2"])))
        print(step_rpp)

        # Attach resource spectra
        other_params(params)$r_pp_array <- step_rpp

    }

      if(is.list(mechanism_params$recruitment_max)) {
        #calculate resulting function
        #add recruitment array to params and declare extra variables
        x <- seq(1,t_max,1)
        time <- x
        values <- params@species_params$R_max
        log_rec <- matrix(ncol = length(time), nrow = length(values), dimnames = list(params@species_params$species,time))
        log_rec[,1] <- values

        rec_chg <- unlist(c(rep(mechanism_params$recruitment_max["rm1"],
                        mechanism_params$recruitment_max["t1"]),
                     rep(mechanism_params$recruitment_max["rm2"],
                         mechanism_params$recruitment_max["t2"])))
        # print(rec_chg)
        for (i in 1:length(time)) {
          log_rec[,i] <- values+values*rec_chg[i]
          log_rec[is.na(log_rec)] <- Inf
        }
          other_params(params)$t_idx = 1

          other_params(params)$log_rec = log_rec
          # print(log_rec)
          params <- setRateFunction(params, "RDD", "TimeBeverton")

      }
    simulation <- project(params, effort = effort_array, t_max = t_max)
  } 
  
  if (mechanism_params$type == "dynamic") {
      if(is.list(mechanism_params$dam_mortality)){
          
        spn <- nrow(species_params(params))
        
        gear_params(params) <- data.frame(gear=c(rep("Dam_Mortality", spn), rep("Migration", spn)),
                                          species = rep(species_params(params)$species,2),
                                          catchability = c(rep(1,spn), rep(0, spn)),
                                          knife_edge_size = c(rep(0.1, 18), gear_params(params)$knife_edge_size))
        
        gear_params(params)$catchability[c(22,36)] <- 1 #replace catchability for migrants back to 1
        
        times = seq(from = 1, to = t_max, by = 1)
        effort_array <- array(NA, dim = c(length(times), 2),
                              dimnames = list(time = times, gear = unique(gear_params(params)$gear)))
        effort_array[,"Migration"] <- rep(3,t_max)
          
        a = as.double(mechanism_params$dam_mortality["a"])
        k = as.double(mechanism_params$dam_mortality["k"])
        b = as.double(mechanism_params$dam_mortality["b"])
        x = seq(1,t_max)

        effort_array[,"Dam_Mortality"] <- c((-k/(1+exp(a+b*x)))+k)
      } else {
      effort_array = 3
    }
    
      if(is.list(mechanism_params$resource_rate)) {
        #calculate resulting function
        #add resource array to params and declare extra variables
        ### r_pp changes rebounds using rebound function
        # Create parameter objects
        other_params(params)$t_idx = 1

        # Attach the new function
        params@resource_dynamics <- "rebound_forcing"

        # Create r_pp array using rebound function
        rebound_resource <- function(m,x,u,b) m-exp(-1*x)*(sin(u*x+b)/sin(b))

        x <- seq(unlist(mechanism_params$resource_rate["start"]),
                 unlist(mechanism_params$resource_rate["finish"]),
                 length.out = t_max)
        print(x)

        rebound_rpp <- rebound_resource(m = as.double(mechanism_params$resource_rate["m"]),
                                        x,
                                        u= as.double(mechanism_params$resource_rate["u"]),
                                        b= as.double(mechanism_params$resource_rate["b"]))
        print(rebound_rpp)
        length(rebound_rpp)

        # Attach resource spectra
        other_params(params)$r_pp_array <- rebound_rpp
    }

      if(is.list(mechanism_params$recruitment_max)) {
        #calculate resulting function
        #add recruitment array to params and declare extra variables
        params <- params

        x <- seq(1,t_max,1)
        time <- x
        values <- params@species_params$R_max
        log_rec <- matrix(ncol = length(time), nrow = length(values), dimnames = list(params@species_params$species,time))
        log_rec[,1] <- values
        l <- as.double(mechanism_params$recruitment_max["l"])
        k <- as.double(mechanism_params$recruitment_max["k"])
        x0 <- as.double(mechanism_params$recruitment_max["x0"])

        rec_chg <- ((l)/ (1+exp(-k*(x-x0))))
        for (i in 1:length(time)) {
          log_rec[,i] <- values+values*rec_chg[i]

        other_params(params)$t_idx = 1

        other_params(params)$log_rec = log_rec

        params <- setRateFunction(params, "RDD", "TimeBeverton")
        }
        
        
      }
    # simulation <- project(params, t_max = t_max)
    simulation <- project(params, effort = effort_array, t_max = t_max)
    # out <- list(x,l,k,x0,log_rec,rec_chg)
  }
  
  return(simulation)
  # return(out)
}


# Tests for function ------------------------------------------------------



# test <- project_dam_removal(params = params, mechanism_params = dyn)
# 
# test <- project_dam_removal(params = params, mechanism_params = sta,t_max = 500)
# plotBiomass(test)

# out <- project(params, t_max = 50)
# 
# many <- list()
# 
# dms <- seq(0,2,length.out =30)
# # dmd <- seq(10,0.1,length.out =30)
# 
# rrs <- seq(10,1,length.out = 30)
# 
# rms <- seq(2, 0.001, length.out = 30)
# 
# for (i in 1:length(dms)) {
#   sta$dam_mortality$e1 <- dms[i]
#   # sta$dam_mortality$t1 <- as.integer(dmd[i])
#   # sta$dam_mortality$t2 <- 50-as.integer(dmd[i])
#   
#   sta$resource_rate$r1 <- rrs[i]
#   
#   sta$recruitment_max$rm2 <- rms[i]
#   
#   many[[i]] <- project_dam_removal(params = params, sta)
# }
# 
# comm_plot(many, out, dms, sz = 15)
# 
# plotlyBiomass(many[[30]])
# # 
# library(viridis)
# 
# mato <- order(many[[15]]@params@species_params$L_inf)
# colors <- viridis(18)
# 
# for (i in 1:length(mato)){
# many[[15]]@params@linecolour[i] <- colors[mato[i]]
# }
# plotlyBiomass(many[[15]])


# basic_plot(out,many[[5]])
# species <- out@params@species_params$species
# wts <- sapply(species, getMeanWeight, sim = many [[15]], min_w <- 1)
# wtdf <- gather(as.data.frame(wts), key = "Species", value = "MeanWeight")
# wtdf$Year <- rep(seq(1,50,1), length(species))
# ggplot(wtdf, aes(x = Year, y = log10(MeanWeight), color = Species)) +
#   geom_line() +
#   theme_minimal()

