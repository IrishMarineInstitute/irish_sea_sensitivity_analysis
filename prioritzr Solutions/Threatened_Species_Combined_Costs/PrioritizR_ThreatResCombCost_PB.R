# ---------------------------------------------------
# Threatened species reserves. All features, Combined costs 
# Author: Patricia Breen 
# ---------------------------------------------------

# Libraries
library(raster)
library(sp)
library(prioritizr)
library(sf)
library(gurobi)
library(dplyr)
library(tibble)
library(scales)
library(ggplot2)
library(topsis)
library(withr)
library(purrr)
library(ggspatial)
library(rgdal)

# ---------------------------------
# 1.0 Get data
# ---------------------------------
# ---------------------------------
# 1.1 Make a raster stack of available feature layers
# ---------------------------------

setwd("./mpa_ecological_sensitivity/Test_Scenarios")

feature_tifs <- list.files(pattern = "\\.tif$")
features <- lapply(feature_tifs, raster)
# Set the names of the list elements to the file names
names(features) <- feature_tifs
# Create a raster stack from the list of rasters
feat_stack_1km <- stack(features)
names(feat_stack_1km)

# ---------------------------------
# 1.2 Make a raster stack of available cost layers
# ---------------------------------

setwd("../data/PreppedCostLayers")

cost_tifs <- list.files(pattern = "\\.tif$")
costs <- lapply(cost_tifs, raster)
# Set the names of the list elements to the file names
names(costs) <- cost_tifs
# Create a raster stack from the list of rasters
cost_stack_1km <- stack(costs)
names(cost_stack_1km)

# ---------------------------------
# 1.3 Get planning units and count number of cells
# ---------------------------------

setwd("../planning_units")
pu_3km <- raster("planning_units_3km.tif")

# ---------------------------------
# 2.0 Aggregate to 3km scale
# ---------------------------------

# Aggregation method averages values across aggregated cells
feat_stack_3km <-aggregate(feat_stack_1km, fact=3)
# Aggregation method sums values across aggregated cells to maintain total costs of layer
cost_stack_3km <-aggregate(cost_stack_1km, fact=3, fun = sum)

cellno <- sum(!is.na(getValues(cost_stack_3km[[1]])))
# ---------------------------------
# 3.0 Set up scenario and solve
# ---------------------------------

# ---------------------------------
# 3.1 Set up cost layer, features and targets 
# ---------------------------------

setwd("../")
setwd("../Solutions")

scenario_name <- "ThreatRes_CombCost_PB" #Each scenario needs a unique name
scenario_group <- "Threatened Species Reserves"

date <- as.character(Sys.time())

# Choose cost layer to apply and adjust slice accordingly.  If no costs, use "pu_3km"
cost <- cost_stack_3km[[1]]
cost_name <- names(cost)
cost_value <- cellStats(cost, sum, na.rm = TRUE)

dir.create(scenario_name)
setwd(paste0("", scenario_name))

# List features and store for later tracking
feature_names <- names(feat_stack_1km)
feature_list <- as.list(feature_names)
included_features <- paste(feature_names, collapse = "; ")

# Set targets and store for later tracking file, also 
#targets <- as.list(rep(0.2, length(feature_names)))
targets <- (rep(0.28, length(feature_names))) #targets 0.28  
targets[c(11, 12, 15, 16, 32, 33)] <- 0.6 #target of 0.6 set for the fronts and for the 4 threatened species
targets_used <- paste(targets, collapse = "; ")

#write target for each feature to text file
featurextarget <- cbind(feature_list, targets)
write.table(featurextarget, file = paste0(scenario_name,"_feature_targets.txt"), sep = ",", row.names = FALSE)

# ---------------------------------
# 3.2 Generate boundary length data for the planning units
# ---------------------------------
bld_3km <- boundary_matrix(pu_3km)
bld_3km@x <-rescale(bld_3km@x, to = c(0.01,100))

# ---------------------------------
# 3.3 Define a range of penalty values
# ---------------------------------
prelim_lower <- -1   # change this for your own data
prelim_upper <- 3  # change this for your own data
prelim_penalty <- round(10^seq(prelim_lower, prelim_upper, length.out = 9), 3)

# print penalty values
print(prelim_penalty)

# ---------------------------------
# 3.4 Pose a 3km-scale Priortizr problem with boundary penalties
# ---------------------------------

p1 <- problem(cost, feat_stack_3km) %>%
  add_min_set_objective() %>%  # minimize the cost of the solution whilst ensuring that all targets are met
  add_relative_targets(targets) %>%
  add_binary_decisions()

save(p1, file = paste0(scenario_name, "_problem.RData"))

# ---------------------------------
# 3.5 Generate Prelims 
# ---------------------------------
# generate preliminary prioritizations based on each penalty
# note that we specify a relaxed gap and time limit for the solver
# adjust solver as required and according to licence constraints

prelim <- lapply(prelim_penalty, function(x) {
  s <-
    p1 %>%
    add_boundary_penalties(penalty = x, data = bld_3km) %>%
    add_neighbor_constraints(k = 2) %>%
    add_gurobi_solver(threads = parallel::detectCores(TRUE) - 1, gap = 0.2, time_limit = 120, verbose = TRUE) %>%
    solve()
  s <- stack(s)
  names(s) <- paste0("bp_", round(x,2))
  return(s)
})

save(prelim, file = paste0(scenario_name, "_prelim.RData")) 

png(paste0(scenario_name, "_prelims.png"), width = 1200, height = 1200, res = 300)
plot(stack(prelim), axes = FALSE, box = FALSE, legend = FALSE)
dev.off()

# ---------------------------------
# 3.6 Generate Solutions
# ---------------------------------

# define a new set of penalty values that ranges below the two preferred prelim values
# This time a linear scale is applied
penalty <- round(seq(prelim_penalty[8], prelim_penalty[9], length.out = 9), 3)

solution <- lapply(penalty, function(x) {
  s <-
    p1 %>%
    add_boundary_penalties(penalty = x, data = bld_3km) %>%
    add_neighbor_constraints(k = 2) %>%
    add_gurobi_solver(threads = parallel::detectCores(TRUE) - 1, gap = 0.1, time_limit = 1200, verbose = TRUE) %>%
    solve()
  s <- stack(s)
  names(s) <- paste0("bp_", round(x,2))
  return(s)
})

save(solution, file = paste0(scenario_name, "_solution.RData")) 
solutions <-stack(solution)

#double check working directory
png(paste0(scenario_name, "_solutions.png"), width = 1200, height = 1200, res = 300)
spplot(solutions, axes = FALSE, box = FALSE)
dev.off()

for (i in 1:nlayers(solutions)) {  
  # Write files  
  soln <- solutions[[i]]  
  writeRaster(soln, paste0(scenario_name, "_bp_", penalty[[i]], ".tif"), format = "GTiff", overwrite = TRUE)}

# ---------------------------------
# 4.0 Evaluate solutions
# ---------------------------------

# ---------------------------------
# 4.1.1 Evaluate feature representations and export to CSV
# ---------------------------------

# Evaluate feature representation in the list of raster solutions
featurereps <- lapply(solution, function(x) {
  eval_feature_representation_summary(p1, x)
})

# Name the tibbles in the list according to the penalty applied
names(featurereps) <- paste0(scenario_name, "_bp_", penalty)

# Examine them in the console
print(featurereps)

# stack the fifth column of all tibbles and bind the data
feat_df <- featurereps %>% 
  map_dfc(~select(.x, relative_held)) %>% 
  set_names(nm = names(featurereps)) %>% 
  mutate(feature = c(featurereps[[1]][[2]]))

# convert to dataframe and assign feature as row names
feat_df <- as.data.frame(feat_df)
rownames(feat_df) <- feat_df$feature
feat_df <- feat_df[, -10]
# reduce to 3 decimal places
feat_df <- round(feat_df,3)
# print the resulting table
feat_df
# Export the summary to CSV
write.csv(feat_df, file=paste0(scenario_name, "_solutions_feature_reps.csv"))

# ---------------------------------
# 4.1.2 Evaluate whether any features have missed their target
# ---------------------------------

# Initialize an empty list to hold the results
OK_solns <- list()

# Loop over each column in soln_df
for (i in 1:ncol(feat_df)) {
  # Extract the values in the current column
  column_values <- feat_df[[i]]
  # Check if all values in the current column are greater than or equal to the corresponding values in targets
  is_OK <- all(column_values >= targets)
  # Append the result to the OK_solns list
  OK_solns[[i]] <- is_OK
}

# Print the resulting list of TRUE/FALSE values
OK_solns

# Iterate through the list to identify those items which are false:

# Initialize an empty list to hold the indices of any FALSE values
false_indices <- list()
# Loop over each element in OK_solns
for (i in seq_along(OK_solns)) {
  # Check if the current element is FALSE
  if (!OK_solns[[i]]) {
    # If it is, append its index to the false_indices list
    false_indices[[length(false_indices) + 1]] <- i
  }
}

# Print the resulting list of false indices
false_indices


# ---------------------------------
# 4.2 Pick best solution
# ---------------------------------

# ---------------------------------
# 4.2.1 Convert 'solutions' (a list of raster stacks) to a list of rasters
# ---------------------------------

# Create an empty list to store the raster layers
solution_layers <- list()
# Loop over each raster stack in the list and extract the first layer
for (i in seq_along(solution)) {
  solution_layers[[i]] <- solution[[i]][[1]]
}
names(solution_layers) <- names(solutions)

# ---------------------------------
# 4.2.2 Remove layers which did not hit targets
# ---------------------------------

# Get the indices of the elements to keep (i.e., the complement of false_indices)
keep_indices <- setdiff(seq_along(solution_layers), unlist(false_indices))
# Subset solution_layers using the keep_indices
OK_layers <- solution_layers[keep_indices]
# Print the resulting filtered list
OK_layers

# ---------------------------------
# 4.2.3 Generate dataframe of relevant parameters with lapply
# ---------------------------------

soln_summary <- lapply(solution_layers, function(x) {   # To apply to ALL layers then use <lapply(soln_layers, function(x)> instead
  clumps <- clump(x)
  clump_sizes <- freq(clumps)
  total_clumps <- length(unique(clumps))
  orphans <- which(clump_sizes[, 2] ==4)
  soln_cost = eval_cost_summary(p1, x)$cost
  boundary <- eval_boundary_summary(p1, x)$boundary
  pu <- eval_n_summary(p1, x)$n
  data.frame(clumps = total_clumps, orphans = length(orphans),
    cost = soln_cost, pu = pu, boundary = boundary, area = pu/cellno)
})

# Subset the penalty list using keep_indices (so as to correctly apply names)
penalty_filtered <- penalty[keep_indices]
penalty_filtered
names(soln_summary) <- paste0("bp_", penalty_filtered)

# Combine the tables into a single dataframe
sum_df<- bind_rows(soln_summary)
sum_df <- as.data.frame(sum_df)
rownames(sum_df) <- names(soln_summary)
sum_df <- round(sum_df,3)
sum_df

# Export the summary to CSV
write.csv(sum_df, file=paste0(scenario_name, "solution_comparisons_3km.csv"))


# ---------------------------------
# 4.2.4 Pick best solution (based on minimizing orphans and cost)
# ---------------------------------

# Get the minimum value of 'orphans' in the 'sum_df'
min_orphans <- min(sum_df$orphans)
# Create a new dataframe 'low_orph_df' that only includes rows where 'orphans' is equal to the minimum or the minimum plus 1
low_orph_df <- subset(sum_df, orphans == min_orphans | orphans == min_orphans + 1)
# Print the resulting dataframe
low_orph_df

# Find the lowest cost in 'low_orph_df' (regardless of the value of 'orphans')
min_cost_all <- min(low_orph_df$cost)
# Calculate the ratio of the lowest cost to 'cost_value'
ratio_all <- min_cost_all / cost_value
ratio_all
# Find the lowest cost in 'low_orph_df' where 'orphans' equals 'min_orphans'
min_cost_min_orphans <- min(low_orph_df$cost[low_orph_df$orphans == min_orphans])
# Calculate the ratio of the lowest cost for 'min_orphans' to 'cost_value' minus 0.02
ratio_min_orphans <- (min_cost_min_orphans / cost_value) - 0.02
ratio_min_orphans

# Determine which ratio is lower and select the corresponding row into 'best_soln_df'
if (ratio_all < ratio_min_orphans) {
  best_soln_df <- low_orph_df[low_orph_df$cost == min_cost_all,]
} else {
  best_soln_df <- low_orph_df[low_orph_df$cost == min_cost_min_orphans & low_orph_df$orphans == min_orphans,]
}

# Print the resulting dataframe and write to a text file
best_soln_df
write.csv(best_soln_df, file = paste0(scenario_name, "_bestsolution.csv"))

# -------------------------------
# 4.2.5 Plot best solution
# -------------------------------

best_solution <- solution_layers[[rownames(best_soln_df)]]
#best_solution <- solution_layers[[1]] #This bit should be automated. See code above. But I had to do this bit manually because 
#somewhere along the line the layer names lost a decimal place and I couldn't find the bug and given time constraints I just had to 
#move on. 

sol_df <- rasterToPoints(best_solution)
sol_df <- data.frame(sol_df)
sol_df[ ,3] <- as.factor(sol_df[ ,3]) 

# ------------------
# 3.8 Plot for final solutions
# ------------------
land <- readOGR("../../data/GIS_layers/IREUK_poly_ITM.shp") #in github folder GIS_layers
sa <- readOGR("../../data/GIS_layers/irishsea_aoi.shp") #in github folder GIS_layers

ggplot() +
  geom_tile(sol_df, mapping = aes(x = x, y = y, fill = sol_df[ ,3]))+
  layer_spatial(sa, fill = NA, colour = "grey80")+
  layer_spatial(land, fill = "grey85", colour = NA)+
  scale_fill_manual(name = "Status", values = c("white","green4"), na.value = NA, labels = c("Not selected", "Selected"))+
  coord_sf(xlim = c(700000,780000), ylim = c(598000,810000))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        rect = element_blank(), 
        legend.position = c(0.2,0.5),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=6),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  )

ggsave(paste0(scenario_name, "_bestsolution.png"), width = 8, height = 6, dpi = 300)

writeRaster(best_solution, paste0(scenario_name, "_bestsolution.tif"), format = "GTiff", overwrite = TRUE)

# ---------------------------------
# 5.0 Evaluate solution costs on a sectoral basis
# ---------------------------------

# Initialize a matrix to store the results
result_matrix <- matrix(0, nrow = nlayers(solutions) * nlayers(cost_stack_3km), ncol = 3, 
                        dimnames = list(NULL, c("Solution Layer", "Cost Layer", "Cost Sum")))

# Loop through each layer in the "solutions" stack
count <- 0
for (i in 1:nlayers(solutions)) {
  
  # Extract the current "solution" layer and its name
  current_solution <- solutions[[i]]
  current_solution_name <- names(solutions)[i]
  
  # Loop through each layer in the "costs" stack
  for (j in 1:nlayers(cost_stack_3km)) {
    # Extract the current "cost" layer and its name
    current_cost <- cost_stack_3km[[j]]
    current_cost_name <- names(cost_stack_3km)[j]
    # multiply rasters to get costs of proposed protected area (areas outside protection are multiplied by zero)
    soln_sect_cost <- current_cost * current_solution
    # Sum the values of the selected cells
    tot_soln_sect_cost <- sum(soln_sect_cost[!is.na(soln_sect_cost)])
    # Update the result matrix with the current results
    count <- count + 1
    result_matrix[count, ] <- c(current_solution_name, current_cost_name, tot_soln_sect_cost)
  }
}

write.csv(result_matrix, file=paste0(scenario_name, "_sector_cost_eval.csv"))

# Convert data matrix to a data frame
result_df <- data.frame(result_matrix)
result_df$Cost.Sum <- as.numeric(result_df$Cost.Sum)

# Remove "X1km_bp_" and ".tif" from Solution.Layer
result_df$Solution.Layer <- gsub("bp_", "", result_df$Solution.Layer) 

# Convert Solution.Layer to numeric
result_df$Solution.Layer <- as.numeric(result_df$Solution.Layer)

# Plot using ggplot
ggplot(result_df, aes(x = Solution.Layer, y = Cost.Sum, color = Cost.Layer, group = Cost.Layer)) +
  geom_line() + 
  geom_point() +
  labs(x = "Boundary penalty", y = "Cost to Sectors", fill = "Sector") +
  theme_classic()

# Save the plot as a PNG file
ggsave(paste0(scenario_name, "_sector_cost_eval.png"), width = 8, height = 6, dpi = 300)

# ---------------------------------
# 6.0 Write to Scenarios.csv
# ---------------------------------

setwd("../")
final_penalty <- rownames(best_soln_df)
final_penalty <- gsub("bp_", "", final_penalty)

# Write details to CSV file for tracking purposes
tracking <- read.csv(file = "Scenarios.csv")   

tracking[nrow(tracking)+1,] <- list (scenario_group, scenario_name, included_features, cost_name,
                                     targets_used, date, final_penalty, best_soln_df$boundary,
                                     best_soln_df$cost, round(best_soln_df$cost/cost_value, 3),
                                     best_soln_df$clumps, best_soln_df$orphans, best_soln_df$pu, best_soln_df$area) 
tracking

write.csv(tracking, file = "Scenarios.csv", row.names = FALSE)

#END