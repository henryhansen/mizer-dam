library(mizer)
library(tidyverse)
library(gridExtra)
library(ggrepel)

# Load projection function
source("_5_project_scenario_functions.R")
source("_6_plotting_functions.R")

# Load params
params <- readRDS("mizer_output/final_params_exp.RDS")


# Scenario Run Individual Effects -----------------------------------------

# Set base static scenario
sta <- list(type ="static",
            dam_mortality =
              list(
                e1 = 0,
                t1 = 5,
                e2 = 0,
                t2 = 45
              ),
            resource_rate =
              list(
                r1 = 10,
                t1 = 20,
                r2 = 10,
                t2 = 30
              )
            ,
            recruitment_max = 
              list(
                rm1 = 0,
                t1 = 10,
                rm2 = 0,
                t2 = 40
              )
)
orsta <- sta #keep original static for new runs

out <- project(params, t_max = 50) # create base community

many_dms <- list() #create container for projections


dms <- seq(0,2,length.out =30) #create range of mortality values

for (i in 1:length(dms)) { #Project different mortality values
  sta$dam_mortality$e1 <- dms[i]
  many_dms[[i]] <- project_dam_removal(params = params, sta)
}
# Export model projections
write_rds(many_dms, "mizer_output/dms_scenarios.rds")
# many_dms <- readRDS("mizer_output/dms_scenarios.rds")
# Plot scenarios
dmbio <- comm_plot(many_dms, out, dms, sz = 15, ptype = "bio") #plot projections
dmbio

sta <- orsta #refresh static scenario 

# do same thing for resource rate
many_rrs <- list()

rrs <- seq(10,0.1,length.out = 30)

for (i in 1:length(dms)) {

  sta$resource_rate$r1 <- rrs[i]
  many_rrs[[i]] <- project_dam_removal(params = params, sta)
}
# Export model projections
write_rds(many_rrs, "mizer_output/rrs_scenarios.rds")
# many_rrs <- read_rds("mizer_output/rrs_scenarios.rds")
# Plot scenarios
rrbio <- comm_plot(many_rrs, out, rrs, sz = 15, ptype = "bio") + scale_colour_viridis_d(direction = -1)
rrbio

sta <- orsta

# do same thing for recruitment
many_rms <- list()

rms <- seq(2, 0, length.out = 30)

for (i in 1:length(dms)) {
  sta$recruitment_max$rm2 <- rms[i]
  # print(sta)
  many_rms[[i]] <- project_dam_removal(params = params, sta)
}
# Export model projections
write_rds(many_rms, "mizer_output/rms_scenarios.rds")
# many_rms <- read_rds("mizer_output/rms_scenarios.rds")
# Plot scenarios
rmsbio <- comm_plot(many_rms, out, rms, sz = 15, ptype = "bio") + scale_colour_viridis_d(direction = -1)
rmsbio

sta <- orsta
# now combine all 3 and plot
many_syn <- list()
for (i in 1:length(dms)) {
  sta$dam_mortality$e1 <- dms[i]

  sta$resource_rate$r1 <- rrs[i]
  
  sta$recruitment_max$rm2 <- rms[i]

  many_syn[[i]] <- project_dam_removal(params = params, sta)
}

# Export model projections
write_rds(many_syn, "mizer_output/synergistic_scenarios.rds")
# many_syn <- read_rds("mizer_output/synergistic_scenarios.rds")

synbio <- comm_plot(many_syn, out, dms, sz = 15, ptype = "bio")
synbio

# Plot scenarios
# combined <- grid.arrange)
library(cowplot)
combined <- plot_grid(dmbio, rrbio + theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.title.y = element_blank() ),
                      rmsbio + theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank() ),
                      synbio + theme(axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.title.y = element_blank() ), nrow =1,labels = c("A", "B", "C", "D"))

combined

ggsave2(file = "mizer_figures/individual_effects.jpeg", 
       combined, 
       units = "px", 
       height = 1500, 
       width = 4000)

# Run simulation for Dynamic Long Term Outcomes ---------------------------

# Set base dynamic run
dyn <- list(type ="dynamic",
            dam_mortality =
              list(
                a = 5,
                k = 0.5,
                b = 0),
            resource_rate =
              list(
                m = 10,
                u = -1.5,
                b = -1,
                start = 0.01,
                finish = 10),
            recruitment_max = list(
              l = 0,
              k = .5,
              x0 = 10
            ))

# dyn <- list(type ="dynamic",
#             dam_mortality =
#               list(
#                 a = 1,
#                 k = 0.5,
#                 b = -0),
#             resource_rate = 0,
#             recruitment_max = 0)
# 
# plotBiomass(project_dam_removal(params, dyn))


many_dyn <- list() #Create container for dynamic runs

d_dms <- seq(-3,-0.05,length.out =30) #create a range of values for each function and parameter change

d_rrs <- seq(-1, -0.05, length.out = 30)

d_rms <- seq(3,0.01,length.out = 30)

for (i in 1:length(d_dms)) {
  dyn$dam_mortality$b <- d_dms[i]
  
  dyn$resource_rate$b <- d_rrs[i]
  
  dyn$recruitment_max$l <- d_rms[i]
  
  many_dyn[[i]] <- project_dam_removal(params = params, dyn, t_max = 200)
}

# Export model projections
write_rds(many_dyn, "mizer_output/dynamic_scenarios.rds")

many_dyn <- read_rds("mizer_output/dynamic_scenarios.rds")

# Plot scenarios
dyn_combo <- comm_plot(many_dyn, values = d_dms, LFI_label = "Proportion of big fish",
                       sz = 20)

dyn_combo

ggsave(file = "mizer_figures/dynamic_effects.png",
       dyn_combo,
       units = "px",
       height = 4000,
       width = 4000)



# Compare output from 1st scenario to last scenario -----------------------
rp1 <- rank_plot(many_dyn[[1]], params = params, sz = 20, ssz = 5)
rp30 <- rank_plot(many_dyn[[30]], params = params, sz = 20, ssz = 5)


all_combo <- grid.arrange(dyn_combo,rp1, rp30, 
                          layout_matrix = rbind(c(1, 1), c(2, 3)))

all_combo
# combine dynamic plots together

ggsave(file = "mizer_figures/dynamic_effects.png", 
       all_combo, 
       units = "px", 
       height = 5000, 
       width = 5000)






