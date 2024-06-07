# load sim data for 25k samples for one run

##################################
# for the actual data runs
###################################

library(ggplot2)
library(dplyr)

data_dir <- "multiSim_10sp"
run <- 10

sitetab <- read.csv(paste0("/space/s1/fiona_callahan/", data_dir, "/randomRun", run, "/sim_sitetab_sampled10000.csv"))
sitetab_abd <- read.csv(paste0("/space/s1/fiona_callahan/", data_dir, "/randomRun", run, "/sim_sitetab_readAbd_sampled10000.csv"))

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


sitetab_site1 <- sitetab_abd[sitetab_abd$Lat == as.numeric(thisLat) & sitetab_abd$Long == as.numeric(thisLong), ]

p2 <- ggplot(sitetab_site1, aes(x = Age, y = Sp2, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "Temporal Trend of Species 1 Abundance",
       x = "Time",
       y = "Abundance") +
  theme_minimal()
p2




##################################
# for the example ones that I simulated just for this
###################################
library(ggplot2)
library(gridExtra)

data_dir <- "/space/s1/fiona_callahan/multiSim_examineSims/100sp_random/randomRun1/"

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
ggsave(abd_plots, filename = paste0(outDir, "abdPlots.png"), width = 16, height = 8)
ggsave(ccplots, filename = paste0(outDir, "ccPlots.png"), width = 16, height = 8)

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
ggsave(avgAbdPlot, file = paste0(outDir, "avgAbdPlot.png"), height = 10, width = 20)

# plot covs over time for spatial locations 1 through 8 (first row in current setup)
covPlotsL <- list()
for(cov in 1:10){
    covariate <- paste0("Cov", cov)
    covPlotsL[[cov]] <- ggplot(data = sitetab_sampled[sitetab_sampled$labID <= 8, ], aes_string(x = "Age", y = covariate, group = "labID", color = "labID")) + # nolint: object_usage_linter, line_length_linter.
        geom_point() +
        geom_line()
}
covplots <- arrangeGrob(grobs = covPlotsL)
ggsave(covplots, file = paste0(outDir, "covPlots.png"), height = 8, width = 10)
