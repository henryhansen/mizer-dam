# Scenario Plots ----------------------------------------------------------
scenario_plot <- function(scenario){
  require(functional)
  require(tidyverse)
  if (scenario$type == "static") {
  # Convert values to relative values
    rel_e1 <- scenario$dam_mortality$e1
    rel_t1 <- scenario$dam_mortality$t1
    rel_e2 <- scenario$dam_mortality$e2
    rel_t2 <- scenario$dam_mortality$t2
    
    rel_rr1 <- (scenario$resource_rate$r1-10)/10
    rel_rrt1 <- scenario$resource_rate$t1
    rel_rr2 <- (scenario$resource_rate$r2-10)/10
    rel_rrt2 <- scenario$resource_rate$t2

    rel_rm1 <- scenario$recruitment_max$rm1/10
    rel_rmt1 <- scenario$recruitment_max$t1
    rel_rm2 <- scenario$recruitment_max$rm2/10
    rel_rmt2 <- scenario$recruitment_max$t2
    
  # Plot relative change of values
  years <- -10:50
  baseline <- rep(0,length(years))
  plot(years,
       baseline,
       type="l", 
       lwd = 10, 
       col = alpha("gray",0.3),
       xlab ="Years since dam removal",
       ylab = "Relative Change")
  lines(years, 
        c(rep(0,11),rep(rel_e1,rel_t1),rep(rel_e2,rel_t2)), 
        type="S", 
        col = "#440154FF", 
        lwd=5, 
        lty ="1342")
  lines(years, 
        c(rep(0,11),rep(rel_rr1,rel_rrt1),rep(rel_rr2,rel_rrt2)), 
        type="S", 
        col = "#2A788EFF", 
        lwd=4, 
        lty ="dashed")
  lines(years, 
        c(rep(0,11),rep(rel_rm1,rel_rmt1),rep(rel_rm2,rel_rmt2)), 
        type="S", 
        col = "#7AD151FF", 
        lwd=3, 
        lty ="dotted")
  abline(v = 0, col = alpha("red", alpha = 0.2), lwd=5)
  legend("bottomright",
         legend = c("Dam Mortality","Resource Rate","Recruitment Max"),
         lwd = 3,
         lty = c("1342","dashed","dotted"),
         col = c("#440154FF","#2A788EFF","#7AD151FF"),
         cex = 0.8)
  mtext("Dam Removal", side = 3, adj = 0.15, col = "Red")
  } else  {
    # Convert values to relative values
    rel_a <- scenario$dam_mortality$a
    rel_k <- scenario$dam_mortality$k
    
    rel_m <- scenario$resource_rate$m
    rel_u <- scenario$resource_rate$u
    rel_b <- scenario$resource_rate$b
    rel_s <- scenario$resource_rate$start
    rel_f <- scenario$resource_rate$finish
    
    rel_l <- scenario$recruitment_max$l
    rel_kr <- scenario$recruitment_max$k
    rel_x0 <- scenario$recruitment_max$x0
    # Plot relative change of values
    
    #exponential decay function
    exp_decay <- function(a,k,x) a*exp(-x*k) 
    rebound_resource <- function(m,x,u,b) scale(m-exp(-1*x)*(sin(u*x+b)/sin(b)),
                                                center = T,
                                                scale = 10)
    logis <- function(l,k,x,x0) scale(((l)/(1+exp(-1*k*(x-x0)))),
                                      center = F,
                                      scale = 1)
    
    years <- -10:50
    baseline <- rep(0,length(years))
    plot(years,
         baseline,
         type="l",
         lwd = 10,
         col = alpha("gray",0.3),
         xlab ="Years since dam removal",
         ylab = "Relative Change")
    
    dm_given <- Curry(exp_decay, a = rel_a, k = rel_k)

    curve(dm_given,
          from = 0,
          to = tail(years,1),
          type="S",
          col = "#440154",
          lwd=3,
          lty ="1342",
          add = T)
    
    rr_given <- Curry(rebound_resource, m = rel_m, u = rel_u, b = rel_b)
    # rel_s
    # rel_f
    curve(rr_given,
          from = 0,
          to = tail(years,1),
          type="S",
          col = "#2A788EFF",
          lwd=4,
          lty ="dashed",
          add = T)
    
    rm_given <- Curry(logis, l = rel_l, k = rel_kr, x0 = rel_x0)
    curve(rm_given,
          from = 0,
          to = tail(years,1),
          type="S",
          col = "#7AD151FF",
          lwd=3,
          lty ="dotted",
          add=T)

    abline(v = 0, col = alpha("red", alpha = 0.2), lwd=5)
    legend("bottomright",
           legend = c("Dam Mortality","Resource Rate","Recruitment Max"),
           lwd = 3,
           lty = c("1342","dashed","dotted"),
           col = c("#440154FF","#2A788EFF","#7AD151FF"),
           cex = 0.8)
    mtext("Dam Removal", side = 3, adj = 0.15, col = "Red")

  }
}


# Rank plot ---------------------------------------------------------------
rank_plot <- function(sim, timepts = c(1,20,50,100,200), params, sz,ssz) {
gdslopes <- as.data.frame(getBiomass(sim)) %>% 
  rownames_to_column("Year") %>% 
  slice(timepts) %>% 
  gather("Species","Biomass",-Year) %>% 
  group_by(Year) %>% 
  mutate(Rank = dense_rank(desc(Biomass)))


rp <- ggplot(data = gdslopes, aes(x = as.integer(Year), y = Rank, group = Species)) +
  geom_line(aes(color = Species), alpha = 1, size = 1.5) +
  geom_point(aes(color = Species), alpha = 1, size = 3) +
  scale_colour_manual(values = getColours(params)[gdslopes$Species]) +
  scale_fill_manual(values = getColours(params)[gdslopes$Species])+
  scale_y_reverse(breaks = c(1, 5, 10, 15, 18)) +
  scale_x_continuous(limits = c(-50,200),labels = c("",seq(0,200,50))) +
  geom_text_repel(data = gdslopes %>% filter(Year == 1), 
                  aes(label = paste0(Species)), 
                  hjust = "left", 
                  fontface = "bold", 
                  size = ssz, 
                  nudge_x = .5, 
                  direction = "x") +
  # geom_text_repel(data = gdslopes %>% filter(Year == 200), 
  #                 aes(label = paste0(Species)) , 
  #                 hjust = "right", 
  #                 fontface = "bold", 
  #                 size = ssz, 
  #                 nudge_x = -.75, 
  #                 direction = "x") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = sz)) +
  theme(axis.title = element_text(size = sz+2)) +
  labs(x = "Year", y = "Biomass Rank")

return(rp)
}




### Basic community indicator comparison function for large fish
basic_plot <- function(base_sim,test_sim) {
  # extract data from base community simulation
  final <- idxFinalT(base_sim)
  lfi0 <- getProportionOfLargeFish(base_sim, min_w = 10, max_w = 100e3,threshold_l = 40)[[final]]
  mw0 <- getMeanWeight(base_sim,min_w = 10, max_w = 100e3)[[final]]
  mmw0 <- getMeanMaxWeight(base_sim, min_w = 10, max_w = 100e3)[final, "mmw_biomass"]
  slope0 <- getCommunitySlope(base_sim, min_w = 10, max_w = 100e3)[final, "slope"]
  # extract data from test community simulaiton
  lfi <- getProportionOfLargeFish(test_sim, min_w = 10, max_w = 100e3,threshold_l = 40)
  mw <- getMeanWeight(test_sim, min_w = 10, max_w = 100e3)
  mmw <- getMeanMaxWeight(test_sim, min_w = 10, max_w = 100e3)[, "mmw_biomass"]
  slope <- getCommunitySlope(test_sim, min_w = 10, max_w = 100e3)[, "slope"]
  
  library(ggplot2)
  years <- 1:length(lfi)
  # Simulated data
  community_plot_data <- rbind(
    data.frame(year = years, measure = "LFI", data = lfi),
    data.frame(year = years, measure = "Mean Weight", data = mw),
    data.frame(year = years, measure = "Mean Max Weight", data = mmw),
    data.frame(year = years, measure = "Slope", data = slope))
  # Unexploited data
  community_unfished_data <- rbind(
    data.frame(year = years, measure = "LFI", data = lfi0),
    data.frame(year = years, measure = "Mean Weight", data = mw0),
    data.frame(year = years, measure = "Mean Max Weight", data = mmw0),
    data.frame(year = years, measure = "Slope", data = slope0))
  # Reference level
  community_reference_level <-
    data.frame(year = years, measure = "LFI", data = lfi0 * 0.8)
  # Build up the plot
  ggplot(community_plot_data) + 
    geom_line(aes(x = year, y = data)) +
    facet_wrap(~measure, scales = "free") + 
    geom_line(aes(x = year, y = data), linetype = "dashed",
              data = community_unfished_data) +
    geom_line(aes(x=year,y=data), linetype = "dotted",
              data = community_reference_level)
}




### Community structure comparisons at end
# compare projections with community biomass ------------------------------

comm_plot <- function(simulations, out, values, sz = 12, ptype = "combined", LFI_label = "Dam Mortality") {
  
library(tidyverse)
extract <- function(x) apply(getBiomass(x),1,sum) #create function to extract comm. biomass

many_ex <- lapply(simulations, extract) #extract values

biomasses <- data.frame(matrix(unlist(many_ex), nrow=length(many_ex), byrow=TRUE))#conver to df

gbio <- gather(biomasses,"timestep","biomass") #gather and label columns

gbio$timestep <- parse_number(gbio$timestep)#fix timestep label

gbio$test_value <- rep(values,idxFinalT(simulations[[1]])) #add test_value column

cbio <- ggplot(gbio, aes(x = timestep, y = biomass, color = factor(test_value))) +
  geom_line() +
  scale_colour_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Year")+
  ylab("Biomass (g)")+ 
  theme(text = element_text(size = sz)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), 
                     limits = c(0,1.25e6), breaks = seq(0,1.25e6,100000))


#  compare projections with community slope (getCommunitySlope)---------------
# project base scenario so it's same length as simulation
# out_ext <- project(out, t_max = idxFinalT(simulations[[1]])-idxFinalT(out))

extract_slope <- function(x) getCommunitySlope(x) #create function to extract comm. biomass

many_ex_slope <- lapply(simulations, extract_slope) #extract values

just_slopes <- lapply(many_ex_slope, `[`, 1)#extract only slopes

slopes <- data.frame(matrix(unlist(just_slopes), nrow=length(just_slopes), byrow=TRUE))#conver to df

gslopes <- gather(slopes,"timestep","slopes") #gather and label columns

gslopes$timestep <- parse_number(gslopes$timestep)#fix timestep label

gslopes$test_value <- rep(values,idxFinalT(simulations[[1]])) #add test_value column

cslope <- ggplot(gslopes, aes(x = timestep, y = slopes, color = factor(test_value))) +
  geom_line() +
  scale_colour_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Year")+
  ylab("Size Spectra Slope")+ 
  theme(text = element_text(size = sz)) 


# compare projections with community mean weights------------------------
# project base scenario so it's same length as simulation

extract_meanw <- function(x) getMeanWeight(x) #create function to extract comm. biomass

many_ex_meanw <- lapply(simulations, extract_meanw) #extract values

meanws <- data.frame(matrix(unlist(many_ex_meanw), nrow=length(many_ex_meanw), byrow=TRUE))#conver to df

gmeanws <- gather(meanws,"timestep","mean_weight") #gather and label columns

gmeanws$timestep <- parse_number(gmeanws$timestep)#fix timestep label

gmeanws$test_value <- rep(values,idxFinalT(simulations[[1]])) #add test_value column

cmw <- ggplot(gmeanws, aes(x = timestep, y = mean_weight, color = factor(test_value))) +
  geom_line() +
  scale_colour_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Year")+
  ylab("Mean Weight (g)") + 
  theme(text = element_text(size = sz))

# compare changes within projections for proportion of large fish  ------------------

extract_bf <- function(x) getProportionOfLargeFish(x) #create function to extract comm. biomass

many_ex_bf <- lapply(simulations, extract_bf) #extract values

bigfish <- data.frame(matrix(unlist(many_ex_bf), nrow=length(many_ex_bf), byrow=TRUE))#conver to df

gbig <- gather(bigfish,"timestep","proportion_big_fish") #gather and label columns

gbig$timestep <- parse_number(gbig$timestep)#fix timestep label

gbig$test_value <- rep(values,idxFinalT(simulations[[1]])) #add test_value column


# old plot
# cbig <- ggplot(gbig, aes(x = timestep, y = test_value, z = proportion_big_fish)) +
#   geom_contour_filled() +
#   xlab("Year")+
#   ylab(LFI_label) + 
#   theme(text = element_text(size = sz)) +
#   scale_fill_viridis_d(option = "C",direction = 1, name = "Proportion of Large Fish") 

cbig <- ggplot(gbig, aes(x = timestep, y = proportion_big_fish, color = factor(test_value))) +
  geom_line() +
  scale_colour_viridis_d() +
  xlab("Year")+
  ylab(LFI_label) + 
  theme_minimal() +
  theme(legend.position = "none") +
  theme(text = element_text(size = sz)) 

# Show all community level plots together ---------------------------------
if (ptype == "combined") {
  library(gridExtra)
  return(gridExtra::grid.arrange(cbio, cslope, cmw, cbig))
} else if (ptype == "bio") {
  return(cbio)
} else if (ptype == "slope") {
  return(cslope)
} else if (ptype == "cmw") {
  return(ccmw)
} else if (ptype == "big") {
  return(cbig)
}

}

