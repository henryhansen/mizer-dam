library(mizer)
library(mizerHowTo)
library(mizerExperimental)
library(tidyverse)

# load species_params -----------------------------------------------------
# morm <- read_csv("morrums?n_no_migs.csv", #import data, allow for swedish characters
#                  locale(encoding = "ISO-8859-1"), 
#                  col_names = T,col_types = NULL)  

# morm <- read.csv("morrumsån_no_migs.csv") 
morm <- read.csv("data/morrumsån_w_migs.csv")
# Encoding(morm$Swedish.Name) <- "UTF-8" 
Encoding(morm$Swedish.Name) <- "UTF-8" 


# Plot Species Densities --------------------------------------------------
# 
# dens <- read_csv("../Fish Data/Species-Occurrence.csv", #import data, allow for swedish characters
#                  locale(encoding = "ISO-8859-1"), 
#                  col_names = T,col_types = NULL)  


dens <- read.csv("../Fish Data/Species-Occurrence.csv") 

dens$Art1 <- factor(dens$Art1)#factor fish species

dens$Fiskedatum1 <- as.Date(dens$Fiskedatum1, "%m/%d/%Y") #make sure Y is capitalized!

dens$Density100m2 <- parse_double(dens$Density100m2, locale = locale(decimal_mark = ","))#convert to "." from ","

ggplot(data = dens, aes(x=Fiskedatum1, y = Density100m2)) +
  # geom_line() +
  geom_point() +
  scale_y_log10()+
  facet_wrap(~factor(Art1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# remove salmon for now
# dens_nm <- dens[!grepl("Lax", dens$Art1),]
dens_nm <- dens

# Plot and Calculate average biomass -------------------------------------------

# create geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} 

ggplot(data = dens_nm, aes(x=Art1,y=Density100m2)) +
  geom_boxplot() +
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Extract geometric mean values from before dam removal (July 1, 2020), compared to after

dens_gm <- dens_nm %>% 
  mutate(dam = if_else(Fiskedatum1 < "2020-07-01", "B","A")) %>% #assign B for before and A for after
  group_by(Art1,dam) %>%
  summarise(gm_d = gm_mean(Density100m2), #extract medians
            n = n()) %>%  #extract number of samples
  ungroup() %>% 
  mutate(freq = n / sum(n)) #calculate proportion of samples where fish is seen


### Based on GIS and sample site, the length of the river from Fridafors to
### the mouth is 34139.70 meters and a median width of 30 meters this corresponds
### to a total river area of
area <- 34139.70 * 30

#Read in length data, find median length for each species, convert to mass

# lens <- read_csv("../Fish Data/Individual-Lengths.csv", #import data, allow for swedish characters
#                  locale(encoding = "ISO-8859-1"), 
#                  col_names = T,col_types = NULL)  

lens <- read.csv("../Fish Data/Individual-Lengths.csv") 

lens$Art1 <- factor(lens$Art1)#factor fish species

lens$Fiskedatum <- as.Date(lens$Fiskedatum, "%m/%d/%Y") #make sure Y is capitalized!

# lens_nm <- lens[!grepl("Lax", lens$Art1),]
lens_nm <- lens

lens_gm <- lens_nm %>% group_by(Art1) %>%
  summarise(gm_l = gm_mean(Längd1)/10 , #extract means in cm
            n_samp = n())  #extract number of samples

# prepare join: drop fish with only one sample from median densities and use only before values
# dens_gm_B <- dens_gm[which(dens_gm$dam == "B" & dens_gm$freq < 1),]
dens_gm_B <- dens_gm[which(dens_gm$dam == "B"),]

community <- left_join(dens_gm_B,lens_gm) # join densities and lengths
community$"Swedish.Name" <- community$Art1 #create species column for next join

# join community information with species params
params <- left_join(morm,community,c("Swedish.Name" ="Swedish.Name"))

# params <- params[complete.cases(params), ] # drop empty rows (fish with one value)

params$gm_w <- (params$a*params$gm_l^params$b)

params$biomass_observed_unscaled <- params$gm_w*params$gm_d*area/100 #multiply converted geometric weight by area of river by 100/m2 density
params$biomass_observed <- params$biomass_observed_unscaled*params$freq #multiply biomass by frequency to scale biomass
params$number_observed <- params$gm_d*area/100*params$freq #calculated abundance

# compare original data densities median proportions to biomass proportions-----

dens_orig <- dens_nm %>% 
  mutate(dam = if_else(Fiskedatum1 < "2020-07-01", "B","A")) %>% #assign B for before and A for after
  group_by(Art1,dam) %>%
  summarise(med = median(Density100m2)) %>%  #extract medians
  ungroup() %>% 
  filter(dam == "B") %>% 
  mutate(perc = med / sum(med) *100)

check <- left_join(dens_orig,params[,c(9, 18,19)])
check <- check[complete.cases(check), ]
check$b_perc <- check$biomass_observed/sum(check$biomass_observed)*100

check_tidy <- gather(data = check[,c(1,4,7)], key = "metric", value = "percent",-1)

ggplot(data = check_tidy, aes(x=Art1,y=percent,fill=metric)) +
  geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### clean up columns for mizer
morm_params <- params[,c(1:8,18:19)]
names(morm_params)[1] <- "species"


rm(list=setdiff(ls(), "morm_params")) #drop other objects


# # load for trout and salmon to average size at departure from 2021 tagging data
salmon_smolts <- read.csv("data/Sam-salmon-smolt-tagging-data-2021.csv", row.names = NULL, sep = ";")
trout_smolts <- read.csv("data/Florian-Trout-Smolt-Tagging-Data-2021.csv", row.names = NULL, sep = ";")

#calculate departure weight for trout and salmon
depart_salmon_w <- mean(salmon_smolts$weight..g.,na.rm = T)
depart_trout_w <- mean(as.integer(subset(trout_smolts, species=="Salmo trutta trutta" & life.stage=="smolt")$weight.in.g), na.rm = T)


# Setup species parameters ------------------------------------------------

#create interactions matrix
morm_inter <- matrix(rep(1, length(morm_params$species)^2), 
                     ncol = nrow(morm_params), 
                     dimnames = list(morm_params$species, morm_params$species))


#create species params object

param_untuned <- newMultispeciesParams(morm_params,morm_inter,kappa = 1e5)
species_params(param_untuned)

# Change size at return for eels

param_untuned@species_params[param_untuned@species_params$species == "Eel", "w_min"] <- 0.226
# cite this paper for value @appelbaumGrowthObservationsEuropean1998


# WGBAST report 2020 
#page 176 - Estimates of wild smolt production (*1000) in Baltic salmon rivers updated with 2019 data. 2001-2020
avg_sal_smolts <- mean(38,34,36,38,39,44,33,28,26,25,35,19,30,49,36,44,28,27,19,35)*1000

#updated salmon smolts numbers observed
param_untuned@species_params[param_untuned@species_params$species == "Atlantic Salmon", "number_observed"] <- avg_sal_smolts
param_untuned@species_params[param_untuned@species_params$species == "Atlantic Salmon", "biomass_observed"] <- avg_sal_smolts*depart_salmon_w

#page 239 -  sea trout smolt estimates 2002-2019
avg_trout_smolts<- mean(6995,3526,5086,5517,10220,6867,3612,5298,3461,3173,2126)

#update trout smolts numbers observed
param_untuned@species_params[param_untuned@species_params$species == "Brown trout", "number_observed"] <- avg_trout_smolts
param_untuned@species_params[param_untuned@species_params$species == "Brown trout", "number_observed"] <- avg_trout_smolts*depart_trout_w
# add size based mortality to trout and salmon over departure length and all others 0
gear_params(param_untuned)[,"catchability"] <- 0
gear_params(param_untuned)[gear_params(param_untuned) == "Brown trout","catchability"] <- 1
gear_params(param_untuned)[gear_params(param_untuned) == "Atlantic Salmon","catchability"] <- 1

gear_params(param_untuned)[gear_params(param_untuned) == "Brown trout","knife_edge_size"] <- depart_trout_w
gear_params(param_untuned)[gear_params(param_untuned) == "Atlantic Salmon","knife_edge_size"] <- depart_salmon_w

# param_untuned@initial_effort[["knife_edge_gear"]] <- 1
param_untuned@initial_effort[["knife_edge_gear"]] <- 3

# # Add migrant reproduction size -----------------------------------------
# 
params <- param_untuned
# # params <- steady(param_untuned)#first look at the size spectra
# plotlySpectra(params, power = 2)
# plotBiomassObservedVsModel(params)
# 
#adjust egg density and size at maturity for trout and salmon so numbers match
species_params(params)[c("Brown trout","Atlantic Salmon"),"min_w_mat"] <- 13
species_params(params)[c("Brown trout","Atlantic Salmon"),"w_mat25"] <- 15



# interactions' -----------------------------------------------------------
coded_params <- params
### Remove predators from diet of salmon, trout, and minnow
btr <- which(row.names(getInteraction(coded_params)) == "Brown trout")
cmw <- which(row.names(getInteraction(coded_params)) == "Common Minnow")
ats <- which(row.names(getInteraction(coded_params)) == "Atlantic Salmon")
lpreds <- c("Carps","Pike","Pikeperch","Burbot")
newinter <- getInteraction(coded_params)
newinter[c(btr,cmw,ats),lpreds] <- 0
coded_params <- setInteraction(coded_params,newinter)

saveRDS(coded_params,"data/noMort.rds")

# # # set external mortality rate to 0.15 for all species
morts <- getExtMort(coded_params)
morts[morts<1] <- 0.15
coded_params@mu_b <- morts

# matchbiomass ------------------------------------------------------------

exp_params <- matchBiomasses(coded_params)

plot(exp_params)
plotBiomassObservedVsModel(exp_params)
getReproductionLevel(exp_params)


steady_params <- steady(exp_params)
plot(steady_params)
plotBiomassObservedVsModel(steady_params)
getReproductionLevel(steady_params)



# cycle through calibration -----------------------------------------------


cycle_params <- steady_params |> calibrateBiomass() |> matchBiomasses() |> steady() |>
  calibrateBiomass() |> matchBiomasses() |> steady() |>
  calibrateBiomass() |> matchBiomasses() |> steady() |>
  calibrateBiomass() |> matchBiomasses() |> steady()


plotBiomassObservedVsModel(cycle_params)
cycle_params@species_params$erepro
plotGrowthCurves(cycle_params, species_panel = TRUE)
getReproductionLevel(cycle_params)
plot(x = log10(species_params(cycle_params)$w_inf),
     y = log10(getRDI(cycle_params)/getRDD(cycle_params)))
text(x = log10(species_params(cycle_params)$w_inf),
     y = log10(getRDI(cycle_params)/getRDD(cycle_params)),
     labels=cycle_params@species_params$species)

plot(cycle_params)
plotBiomass(project(cycle_params, t_max = 500))


# # setbeverton holt --------------------------------------------------------

#consider fixing erepro and preserving it and tuning only the ones that are going extinct (Eel and bullhead)

test <- setBevertonHolt(cycle_params, erepro = c(10, #
                                                 1,
                                                 0.1,
                                                 0.001,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 0.1,
                                                 1,
                                                 0.001))  # original value 0.1



test <- steady(test, preserve = c("erepro"))

plotBiomassObservedVsModel(test)
test@species_params$erepro
plot(x = log10(species_params(test)$w_inf),
     y = log10(getRDI(test)/getRDD(test)))
text(x = log10(species_params(test)$w_inf),
     y = log10(getRDI(test)/getRDD(test)),
     labels=test@species_params$species)
round(getReproductionLevel(test),3)

plot(project(test, effort = 3))



# final param -------------------------------------------------------------

final_params <- projectToSteady(test, t_max = 250, progress_bar = F, tol = 0.001) 

round(getReproductionLevel(final_params), 5)

plotBiomassObservedVsModel(final_params)
plotBiomass(project(final_params))
plot(final_params)
plot(x = log10(species_params(final_params)$w_inf),
     y = log10(getRDI(final_params)/getRDD(final_params)),
     xlab ="W_inf",
     ylab = "Log10(RDI/RDD)",
     main = "RDI/RDD plot",
     xlim =c(0,5)
)
text(x = log10(species_params(final_params)$w_inf),
     y = log10(getRDI(final_params)/getRDD(final_params)),
     labels=final_params@species_params$species)

# add color scheme
# set colors based on Linf but smolts are changed to reflect departure size
linf <- species_params(final_params)$L_inf
linf[4] <- mean(trout_smolts$length.in.mm, na.rm = T)/10
linf[18] <- mean(salmon_smolts$length..mm., na.rm = T)/10
sizeorder <- order(linf)

colors <- viridis::viridis(18, begin = 0.02, end = 0.98)
for (i in 1:length(sizeorder)) {
  final_params@linecolour[i] <- colors[which(sizeorder == i)]
}

final_params@linetype[19] <- "dashed"

saveRDS(final_params, "mizer_output/final_params_exp.RDS")

final_params <- readRDS("mizer_output/final_params_exp.RDS")



# Export Calibration ------------------------------------------------------

library(gridExtra)
library(ggrepel)


bmp <- plotBiomassObservedVsModel(final_params) + 
  theme_minimal() +
  theme(legend.position="none")
bmp

# create data for rdi/rdi plot
ddd <- data.frame(species = species_params(final_params)$species,
                  RDI = getRDI(final_params),
                  RDD = getRDD(final_params),
                  L_inf = linf)

dddplot <- ggplot(ddd, aes(x = log10(linf), y = log10(RDI/RDD))) +
  geom_point(aes(color = species)) + 
  theme_minimal() + 
  labs(x = expression(paste("log10(","L",infinity," [cm])")), 
       y = "log10(RDI/RDD)")  +
  theme(legend.position="none") +
  geom_label_repel(aes(label = species,
                   color = species),
                   fill = 'white',
                   size = 3.5,
                   min.segment.length = 0.1) +
  scale_colour_manual(values = getColours(final_params)[ddd$species]) +
  scale_fill_manual(values = getColours(final_params)[ddd$species])+
  scale_y_continuous(limits = c(-1,7))

combined <- gridExtra::grid.arrange(bmp, dddplot, ncol = 2)

ggsave(file = "mizer_figures/final_calibration.png", 
       combined, 
       units = "px", 
       height = 3000, 
       width = 4000)


# Plot coexistance and spectra ---------------------------------------------------
library(RColorBrewer)
final_params@linecolour[1:18] <- c(brewer.pal(11, name = "Set3"),brewer.pal(7, name = "Set1"))

top <- plotBiomass(project(final_params)) +
  theme_minimal() +
  theme(legend.position="none") 

bottom <- plotSpectra(final_params) + 
  scale_y_log10(limits = c(1e-06, 1e10)) +
  theme_minimal() 

tb <- gridExtra::grid.arrange(top,bottom, ncol = 1)

ggsave(file = "mizer_figures/coexist_spectra.png", 
       tb, 
       units = "px", 
       height = 3000, 
       width = 4000)

# old parameterization with custom erepro for all species -----------------

# test <- setBevertonHolt(cycle_params, erepro = c(10, #
#                                                  1,
#                                                  0.00001,
#                                                  0.001,
#                                                  0.1,
#                                                  0.001,
#                                                  0.1,
#                                                  0.1,
#                                                  0.00001,
#                                                  0.001,
#                                                  0.001,
#                                                  0.00001,
#                                                  0.1,
#                                                  0.1,
#                                                  0.1,
#                                                  0.1,
#                                                  1,
#                                                  0.001))  # original value 0.1
# test <- steady(test, preserve = c("R_max"))#original value R_max


