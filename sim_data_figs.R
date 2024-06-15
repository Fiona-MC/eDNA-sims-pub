# load sim data for 25k samples for one run

##################################
# for the actual data runs
###################################

library(ggplot2)
library(dplyr)
library(tidyr)

data_dir <- "multiSim_10sp_random_moreSamples"
run <- 10
numSamples <- 10000

sim_figs_dir <- paste0("/space/s1/fiona_callahan/", data_dir, "/simDataFigs/")
dir.create(sim_figs_dir)

sitetab <- read.csv(paste0("/space/s1/fiona_callahan/", data_dir, "/randomRun", run, "/sim_sitetab_sampled", numSamples, ".csv"))
sitetab_abd <- read.csv(paste0("/space/s1/fiona_callahan/", data_dir, "/randomRun", run, "/sim_sitetab_readAbd_sampled", numSamples, ".csv"))

names(sitetab_abd)

# Find the most common Lat Long pair
most_common_latlong <- sitetab_abd %>%
  group_by(Lat, Long) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  slice(1)

print(most_common_latlong)

thisLat <- most_common_latlong[1, 1]
thisLong <- most_common_latlong[1, 2]

p1 <- ggplot(sitetab_abd, aes(x = Long, y = Lat, color = Sp2)) +
  geom_point(size = 30, shape = 15) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Spatial Distribution of Species 1 Abundance",
       x = "Longitude",
       y = "Latitude",
       color = "Abundance") +
  theme_minimal()
p1

ggsave(plot = p1, filename = paste0(sim_figs_dir, "abd_byLocation_sampled", numSamples, "_run", run, ".pdf"))

sitetab_site1 <- sitetab_abd[sitetab_abd$Lat == as.numeric(thisLat) & sitetab_abd$Long == as.numeric(thisLong), ]

p2 <- ggplot(sitetab_site1, aes(x = Age, y = Sp3, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "Temporal Trend of Species 1 Abundance",
       x = "Time",
       y = "Abundance") +
  theme_minimal()
p2

#ggsave(plot = p2, filename = paste0(sim_figs_dir, "abd_byTime_", numSamples, "_run", run, ".pdf"))


# Reshape the data into long format
sitetab_long_site1 <- sitetab_site1 %>%
  pivot_longer(cols = starts_with("Sp"), names_to = "Species", values_to = "Abundance") %>%
  filter(Species %in% c("Sp1", "Sp2", "Sp3", "Sp4", "Sp5"))

# Create the plot
p3 <- ggplot(sitetab_long_site1, aes(x = Age, y = Abundance, color = Species, group = Species)) +
  geom_line() +
  geom_point() +
  labs(title = "Temporal Trend of Species Abundance",
       x = "Time",
       y = "Abundance") +
  theme_minimal()
p3
ggsave(plot = p3, filename = paste0(sim_figs_dir, "abd_byTime_", numSamples, "_run", run, ".pdf"), width = 20, height = 10)



# Reshape the data into long format
sitetab_long <- sitetab_abd %>%
  pivot_longer(cols = starts_with("Sp"), names_to = "Species", values_to = "Abundance") %>%
  filter(Species %in% c("Sp1", "Sp2", "Sp3", "Sp4", "Sp5"))

# Create the plot
p4 <- ggplot(sitetab_long, aes(x = Age, y = Abundance, color = Species, group = Species)) +
  geom_line() +
  geom_point() +
  labs(title = "Temporal Trend of Species Abundance",
       x = "Time",
       y = "Abundance") +
  theme_minimal()
p4
ggsave(plot = p4, filename = paste0(sim_figs_dir, "abd_byTime_allLoc_", numSamples, "_run", run, ".pdf"), width = 20, height = 10)










##################################
# for the example ones that I simulated just for this
###################################
library(ggplot2)
library(gridExtra)

data_dir <- "/space/s1/fiona_callahan/multiSim_examineSims/10sp_random/randomRun1/"

sim_data <- readRDS(paste0(data_dir, "sim_data.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))
locList <- readRDS(paste0(data_dir, "locList.Rdata"))

locPlots <- c(19, 8, 6, 37, 64, 46)
outDir <- data_dir
sitetab_sampled <- read.csv(paste0(data_dir, "sim_sitetab_sampled.csv"))
spList <- 1:10
readAbdMode <- TRUE

#wrangle data
# this creates a list of length of number of locations and in each index there is a DF with colums time, species1, species2, species3
species_abundance <- list() 
for (s_index in seq_along(locList)) {
    # make data frame for that location
    species_abundance[[s_index]] <- data.frame(time = 1:params$num_gens)
    for (k in spList) {
        spAbund <- lapply(sim_data, FUN = function(x) {return(x$N[[s_index]][k])}) # nolint: brace_linter.
        species_abundance[[s_index]][[paste0("species", k)]] <- unlist(spAbund)

        spReadAbund <- lapply(sim_data, FUN = function(x) {return(x$readAbd[[s_index]][k])}) # nolint: brace_linter.
        species_abundance[[s_index]][[paste0("species", k, "reads")]] <- unlist(spReadAbund)

        spPresAbs <- lapply(sim_data, FUN = function(x) {return(x$y[[s_index]][k])}) # nolint: brace_linter.
        species_abundance[[s_index]][[paste0("species", k, "pres")]] <- unlist(spPresAbs)
    }
}

#wrangle data
carrying_capacities <- list()
for (s_index in seq_along(locList)){
    # make data frame for that location
    carrying_capacities[[s_index]] <- data.frame(time = 1:params$num_gens)
    for (k in spList){
    Kcapacity <- lapply(sim_data, FUN = function(x) {return(x$K[[s_index]][k])})
    carrying_capacities[[s_index]][[paste0("species", k)]] <- unlist(Kcapacity)
    }
}


readsVsAbd <- ggplot(species_abundance[[46]], aes(x = species1, y = species1reads)) +
  geom_point(size = 3) +
  labs(x = "Species 1 Abundance", 
        y = "Species 1 read count") 
readsVsAbd
ggsave(readsVsAbd, filename = paste0(outDir, "readsVsAbd.pdf"))


species_abundance[[46]]$species1pres <- as.factor(species_abundance$species1pres)
# Create the box plot 
readsVsPres <- ggplot(species_abundance[[46]], aes(x = species1pres, y = species1, group = species1pres)) +
  geom_boxplot() +
  theme_minimal()

readsVsPres

ggsave(readsVsPres, filename = paste0(outDir, "readsVsPres.pdf"))


#plot species against each other
s_indexL <- locPlots
abd_plotList <- list()
cc_plotList <- list()
for (i in seq_along(s_indexL)){
    s_index <- s_indexL[i]
    # species_abundance[[s_index]]
    df <- reshape::melt(species_abundance[[s_index]][species_abundance[[s_index]]$time %in% 100:1000, ], id.vars = "time")
    names(df) <- c("time", "species", "abundance")
    abd_plot <- ggplot(df, aes(x = time, y = abundance, color = species)) + # nolint: object_usage_linter.
    geom_line() +
    ggtitle(paste("Species abundance for location", s_index))
    
    carrying_capacities[[s_index]]
    df <- reshape::melt(carrying_capacities[[s_index]][carrying_capacities[[s_index]]$time %in% 100:1000, ], id.vars = "time")
    names(df) <- c("time", "species", "carrying_capacity")
    cc_plot <- ggplot(df, aes(x = time, y = carrying_capacity, color = species)) + # nolint: object_usage_linter.
    geom_line() +
    ggtitle(paste("Carrying capacities for location", s_index))
    abd_plotList[[i]] <- abd_plot
    cc_plotList[[i]] <- cc_plot
}

abd_plots <- arrangeGrob(abd_plotList[[1]], abd_plotList[[2]], abd_plotList[[3]], abd_plotList[[4]],
                        abd_plotList[[5]], abd_plotList[[6]], nrow = 2)
ccplots <- arrangeGrob(cc_plotList[[1]], cc_plotList[[2]], cc_plotList[[3]], cc_plotList[[4]], 
                        cc_plotList[[5]], cc_plotList[[6]], nrow = 2)
ggsave(abd_plots, filename = paste0(outDir, "abdPlots.pdf"), width = 16, height = 8)
ggsave(ccplots, filename = paste0(outDir, "ccPlots.pdf"), width = 16, height = 8)

#Plot average popn over space
#data wrangling
x_coords <- lapply(locList, FUN = function(x) {return(x[1])})
y_coords <- lapply(locList, FUN = function(x) {return(x[2])}) 

sp_average <- lapply(species_abundance, FUN = function(sp_abd) {
    return(apply(X = sp_abd, MARGIN = 2, FUN = mean))
})

abd_av_df <- data.frame(s_index = seq_along(locList), x = unlist(x_coords), y = unlist(y_coords))
for (k in spList){
    #just cause the first column isnt real
    sp_i <- k + 1
    sp_av <- lapply(sp_average, FUN = function(spav) {return(spav[sp_i])})
    abd_av_df[[paste0("species", k, "AverageAbundance")]] <- unlist(sp_av)
}
#abd_av_df

aesNamesString <- lapply(spList, FUN = function(sp) {paste0("species", sp, "AverageAbundance")})
#plot average abundance of each species across space
avgAbdPlotsL <- lapply(aesNamesString, FUN = function(varName) {
    p <- ggplot(abd_av_df, aes(x = x, y = y, label = s_index, color = !!sym(varName))) + # nolint: object_usage_linter.
    geom_text() +
    scale_color_gradient(low = "blue", high = "green")
    #return(p)
})

avgAbdPlot <- arrangeGrob(grobs = avgAbdPlotsL)
ggsave(avgAbdPlot, file = paste0(outDir, "avgAbdPlot.pdf"), height = 10, width = 20)

# plot covs over time for spatial locations 1 through 8 (first row in current setup)
covPlotsL <- list()
for(cov in 1:params$numCovs){
    covariate <- paste0("Cov", cov)
    covPlotsL[[cov]] <- ggplot(data = sitetab_sampled[sitetab_sampled$labID <= 8, ], aes_string(x = "Age", y = covariate, group = "labID", color = "labID")) + # nolint: object_usage_linter, line_length_linter.
        geom_point() +
        geom_line()
}
covplots <- arrangeGrob(grobs = covPlotsL)
ggsave(covplots, file = paste0(outDir, "covPlots.pdf"), height = 8, width = 10)
