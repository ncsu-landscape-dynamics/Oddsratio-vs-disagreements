## Data creation for comparing calibration and validation methods 4 datasets
## create 200x200 grid of 250 m resolution pixels (dataset 1)
devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "feature/oddsratio", force = TRUE)
library(raster)
library(PoPS)

ncol <- 200
nrow <- 200
number <- nrow *ncol
xmax <- nrow * 250
ymax <- ncol * 250
# create random host map with all_plants being a raster with all values set to the max that was generated in host map
host <- matrix(round(runif(number, 0, 100)), nrow = nrow, ncol = ncol)
cols <- round(runif(number*.4, 1, ncol))
rows <- round(runif(number*.4, 1, nrow))
for (zeros in 1:number*.4) {
  host[cols[zeros],rows[zeros]] <- 0
}
all_plants <- matrix(max(host), nrow = nrow, ncol = ncol)
## create initial infections for single and multiple infection scenarios
initial_infection <- matrix(0, nrow = nrow, ncol = ncol)
column <- round(runif(1, 1, ncol))
row <- round(runif(1, 1, nrow))
initial_infection[column, row] <- initial_infection[column, row] + round(runif(1, 1, 8))
initial_infection_multi <- initial_infection
# Add five new location of initial infection for multiple initial infection scenario.
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))

## convert to raster for saving
host <- raster(host, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
all_plants <- raster(all_plants, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection <- raster(initial_infection, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection_multi <- raster(initial_infection_multi, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)

writeRaster(host, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/host.tif", overwrite = TRUE)
writeRaster(all_plants, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/all_plants.tif", overwrite = TRUE)
writeRaster(initial_infection, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_single.tif", overwrite = TRUE)
writeRaster(initial_infection_multi, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_multi.tif", overwrite = TRUE)


host_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/host.tif"
total_plants_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/all_plants.tif"
temperature_file <- ""
temperature_coefficient_file <- ""
precipitation_coefficient_file <-""
use_lethal_temperature <- FALSE
temp <- FALSE
precip <- FALSE
season_month_start <- 1
season_month_end <- 12
time_step <- "month"
start_time <- 2012
end_time <- 2015
dispersal_kern <- "cauchy"
percent_short_distance_dispersal <- 1.0
long_distance_scale <- 0.0
lethal_temperature <- 0
lethal_temperature_month <- 1
wind_dir <- "NONE"
kappa <- 0
random_seed <- 42
treatments_file <- ""
treatment_years <- c(0)
management <- FALSE
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0

## reproductive rate = 0.7 distance scale = 21 at extent of 200 by 200 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_single.tif"
short_distance_scale <- 21
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr0.7_sd21_extent200_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 21 at extent of 200 by 200 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_multi.tif"
short_distance_scale <- 21
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr0.7_sd21_extent200_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 83 at extent of 200 by 200 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_single.tif"
short_distance_scale <- 83
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr0.7_sd83_extent200_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 83 at extent of 200 by 200 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_multi.tif"
short_distance_scale <- 83
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr0.7_sd83_extent200_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 21 at extent of 200 by 200 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_single.tif"
short_distance_scale <- 21
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr1.5_sd21_extent200_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 21 at extent of 200 by 200 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_multi.tif"
short_distance_scale <- 21
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr1.5_sd21_extent200_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 83 at extent of 200 by 200 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_single.tif"
short_distance_scale <- 83
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr1.5_sd83_extent200_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 83 at extent of 200 by 200 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_multi.tif"
short_distance_scale <- 83
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr1.5_sd83_extent200_multi.tif", overwrite = TRUE)
plot(infected_stack)

## create 1000 x 1000 grid of 250 m resolution pixels (dataset 2)
ncol <- 1000
nrow <- 1000
number <- nrow *ncol
xmax <- nrow * 250
ymax <- ncol * 250
# create random host map with all_plants being a raster with all values set to the max that was generated in host map
host <- matrix(round(runif(number, 0, 100)), nrow = nrow, ncol = ncol)
cols <- round(runif(number*.4, 1, ncol))
rows <- round(runif(number*.4, 1, nrow))
for (zeros in 1:number*.4) {
  host[cols[zeros],rows[zeros]] <- 0
}
all_plants <- matrix(max(host), nrow = nrow, ncol = ncol)
## create initial infections for single and multiple infection scenarios
initial_infection <- matrix(0, nrow = nrow, ncol = ncol)
column <- round(runif(1, 1, ncol))
row <- round(runif(1, 1, nrow))
initial_infection[column, row] <- initial_infection[column, row] + round(runif(1, 1, 8))
initial_infection_multi <- initial_infection
# Add five new location of initial infection for multiple initial infection scenario.
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))

## convert to raster for saving
host <- raster(host, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
all_plants <- raster(all_plants, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection <- raster(initial_infection, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection_multi <- raster(initial_infection_multi, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)

writeRaster(host, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/host.tif", overwrite = TRUE)
writeRaster(all_plants, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/all_plants.tif", overwrite = TRUE)
writeRaster(initial_infection, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_single.tif", overwrite = TRUE)
writeRaster(initial_infection_multi, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_multi.tif", overwrite = TRUE)

host_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/host.tif"
total_plants_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/all_plants.tif"
temperature_file <- ""
temperature_coefficient_file <- ""
precipitation_coefficient_file <-""
use_lethal_temperature <- FALSE
temp <- FALSE
precip <- FALSE
season_month_start <- 1
season_month_end <- 12
time_step <- "month"
start_time <- 2012
end_time <- 2015
dispersal_kern <- "cauchy"
percent_short_distance_dispersal <- 1.0
long_distance_scale <- 0.0
lethal_temperature <- 0
lethal_temperature_month <- 1
wind_dir <- "NONE"
kappa <- 0
random_seed <- 42
treatments_file <- ""
treatment_years <- c(0)
management <- FALSE
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0

## reproductive rate = 0.7 distance scale = 21 at extent of 1000 by 1000 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_single.tif"
short_distance_scale <- 21
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/infected_years_rr0.7_sd21_extent1000_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 21 at extent of 1000 by 1000 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_multi.tif"
short_distance_scale <- 21
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/infected_years_rr0.7_sd21_extent1000_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 83 at extent of 1000 by 1000 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_single.tif"
short_distance_scale <- 83
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/infected_years_rr0.7_sd83_extent1000_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 83 at extent of 1000 by 1000 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_multi.tif"
short_distance_scale <- 83
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/infected_years_rr0.7_sd83_extent1000_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 21 at extent of 1000 by 1000 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_single.tif"
short_distance_scale <- 21
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/infected_years_rr1.5_sd21_extent1000_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 21 at extent of 1000 by 1000 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_multi.tif"
short_distance_scale <- 21
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/infected_years_rr1.5_sd21_extent1000_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 83 at extent of 1000 by 1000 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_single.tif"
short_distance_scale <- 83
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/infected_years_rr1.5_sd83_extent1000_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 83 at extent of 1000 by 1000 with resoultion of 250 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/initial_infections_multi.tif"
short_distance_scale <- 83
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/1000/infected_years_rr1.5_sd83_extent1000_multi.tif", overwrite = TRUE)
plot(infected_stack)

## create 200x200 grid of 1000 m resolution pixels (dataset 3)
ncol <- 200
nrow <- 200
number <- nrow *ncol
xmax <- nrow * 1000
ymax <- ncol * 1000
# create random host map with all_plants being a raster with all values set to the max that was generated in host map
host <- matrix(round(runif(number, 0, 100)), nrow = nrow, ncol = ncol)
cols <- round(runif(number*.4, 1, ncol))
rows <- round(runif(number*.4, 1, nrow))
for (zeros in 1:number*.4) {
  host[cols[zeros],rows[zeros]] <- 0
}
all_plants <- matrix(max(host), nrow = nrow, ncol = ncol)
## create initial infections for single and multiple infection scenarios
initial_infection <- matrix(0, nrow = nrow, ncol = ncol)
column <- round(runif(1, 1, ncol))
row <- round(runif(1, 1, nrow))
initial_infection[column, row] <- initial_infection[column, row] + round(runif(1, 1, 8))
initial_infection_multi <- initial_infection
# Add five new location of initial infection for multiple initial infection scenario.
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))

## convert to raster for saving
host <- raster(host, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
all_plants <- raster(all_plants, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection <- raster(initial_infection, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection_multi <- raster(initial_infection_multi, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)

writeRaster(host, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/host.tif", overwrite = TRUE)
writeRaster(all_plants, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/all_plants.tif", overwrite = TRUE)
writeRaster(initial_infection, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_single.tif", overwrite = TRUE)
writeRaster(initial_infection_multi, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_multi.tif", overwrite = TRUE)


host_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/host.tif"
total_plants_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/all_plants.tif"
temperature_file <- ""
temperature_coefficient_file <- ""
precipitation_coefficient_file <-""
use_lethal_temperature <- FALSE
temp <- FALSE
precip <- FALSE
season_month_start <- 1
season_month_end <- 12
time_step <- "month"
start_time <- 2012
end_time <- 2015
dispersal_kern <- "cauchy"
percent_short_distance_dispersal <- 1.0
long_distance_scale <- 0.0
lethal_temperature <- 0
lethal_temperature_month <- 1
wind_dir <- "NONE"
kappa <- 0
random_seed <- 42
treatments_file <- ""
treatment_years <- c(0)
management <- FALSE
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0

## reproductive rate = 0.7 distance scale = 21 at extent of 200 by 200 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_single.tif"
short_distance_scale <- 21
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/infected_years_rr0.7_sd21_extent200_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 21 at extent of 200 by 200 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_multi.tif"
short_distance_scale <- 21
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/infected_years_rr0.7_sd21_extent200_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 83 at extent of 200 by 200 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_single.tif"
short_distance_scale <- 83
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/infected_years_rr0.7_sd83_extent200_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 83 at extent of 200 by 200 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_multi.tif"
short_distance_scale <- 83
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/infected_years_rr0.7_sd83_extent200_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 21 at extent of 200 by 200 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_single.tif"
short_distance_scale <- 21
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/infected_years_rr1.5_sd21_extent200_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 21 at extent of 200 by 200 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_multi.tif"
short_distance_scale <- 21
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/infected_years_rr1.5_sd21_extent200_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 83 at extent of 200 by 200 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_single.tif"
short_distance_scale <- 83
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/infected_years_rr1.5_sd83_extent200_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 83 at extent of 200 by 200 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/initial_infections_multi.tif"
short_distance_scale <- 83
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/200/infected_years_rr1.5_sd83_extent200_multi.tif", overwrite = TRUE)
plot(infected_stack)

## create 1000 x 1000 grid of 1000 m resolution pixels(dataset 4)
ncol <- 1000
nrow <- 1000
number <- nrow *ncol
xmax <- nrow * 1000
ymax <- ncol * 1000
# create random host map with all_plants being a raster with all values set to the max that was generated in host map
host <- matrix(round(runif(number, 0, 100)), nrow = nrow, ncol = ncol)
cols <- round(runif(number*.4, 1, ncol))
rows <- round(runif(number*.4, 1, nrow))
for (zeros in 1:number*.4) {
  host[cols[zeros],rows[zeros]] <- 0
}
all_plants <- matrix(max(host), nrow = nrow, ncol = ncol)
## create initial infections for single and multiple infection scenarios
initial_infection <- matrix(0, nrow = nrow, ncol = ncol)
column <- round(runif(1, 1, ncol))
row <- round(runif(1, 1, nrow))
initial_infection[column, row] <- initial_infection[column, row] + round(runif(1, 1, 8))
initial_infection_multi <- initial_infection
# Add five new location of initial infection for multiple initial infection scenario.
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))
col1 <- column + round(runif(1, -10, 10))
row1 <- row + round(runif(1, -10, 10))
initial_infection_multi[col1, row1] <- initial_infection_multi[col1, row1] + round(runif(1, 0, 8))

## convert to raster for saving
host <- raster(host, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
all_plants <- raster(all_plants, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection <- raster(initial_infection, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection_multi <- raster(initial_infection_multi, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)

writeRaster(host, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/host.tif", overwrite = TRUE)
writeRaster(all_plants, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/all_plants.tif", overwrite = TRUE)
writeRaster(initial_infection, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_single.tif", overwrite = TRUE)
writeRaster(initial_infection_multi, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_multi.tif", overwrite = TRUE)

host_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/host.tif"
total_plants_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/all_plants.tif"
temperature_file <- ""
temperature_coefficient_file <- ""
precipitation_coefficient_file <-""
use_lethal_temperature <- FALSE
temp <- FALSE
precip <- FALSE
season_month_start <- 1
season_month_end <- 12
time_step <- "month"
start_time <- 2012
end_time <- 2015
dispersal_kern <- "cauchy"
percent_short_distance_dispersal <- 1.0
long_distance_scale <- 0.0
lethal_temperature <- 0
lethal_temperature_month <- 1
wind_dir <- "NONE"
kappa <- 0
random_seed <- 42
treatments_file <- ""
treatment_years <- c(0)
management <- FALSE
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0

## reproductive rate = 0.7 distance scale = 21 at extent of 1000 by 1000 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_single.tif"
short_distance_scale <- 21
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/infected_years_rr0.7_sd21_extent1000_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 21 at extent of 1000 by 1000 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_multi.tif"
short_distance_scale <- 21
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/infected_years_rr0.7_sd21_extent1000_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 83 at extent of 1000 by 1000 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_single.tif"
short_distance_scale <- 83
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/infected_years_rr0.7_sd83_extent1000_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 0.7 distance scale = 83 at extent of 1000 by 1000 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_multi.tif"
short_distance_scale <- 83
reproductive_rate <- 0.7

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/infected_years_rr0.7_sd83_extent1000_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 21 at extent of 1000 by 1000 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_single.tif"
short_distance_scale <- 21
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/infected_years_rr1.5_sd21_extent1000_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 21 at extent of 1000 by 1000 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_multi.tif"
short_distance_scale <- 21
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/infected_years_rr1.5_sd21_extent1000_multi.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 83 at extent of 1000 by 1000 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_single.tif"
short_distance_scale <- 83
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/infected_years_rr1.5_sd83_extent1000_single.tif", overwrite = TRUE)
plot(infected_stack)

## reproductive rate = 1.5 distance scale = 83 at extent of 1000 by 1000 with resoultion of 1000 m with single initial infection location
infected_file <- "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/initial_infections_multi.tif"
short_distance_scale <- 83
reproductive_rate <- 1.5

data <- PoPS::pops(infected_file, host_file, total_plants_file, reproductive_rate,
                   use_lethal_temperature, temp, precip, management, mortality_on,
                   temperature_file, temperature_coefficient_file,
                   precipitation_coefficient_file, treatments_file,
                   season_month_start, season_month_end, time_step,
                   start_time, end_time, treatment_years,
                   dispersal_kern, percent_short_distance_dispersal,
                   short_distance_scale, long_distance_scale,
                   lethal_temperature, lethal_temperature_month,
                   mortality_rate, mortality_time_lag,
                   wind_dir, kappa, random_seed)

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (j in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[j]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}
writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/1000/1000/infected_years_rr1.5_sd83_extent1000_multi.tif", overwrite = TRUE)
plot(infected_stack)