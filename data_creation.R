## Data creation for comparing calibration and validation methods 2 datasets
## 250 x 250 data creation (dataset 1)
devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "feature/oddsratio")
library(raster)
library(PoPS)
ncol = 250
nrow = 250
number = nrow *ncol
xmax = nrow * 100
ymax = ncol * 100
host <- matrix(round(runif(number, 0, 200)), nrow = nrow, ncol = ncol)
all_plants <- matrix(max(host), nrow = nrow, ncol = ncol)
initial_infection <- matrix(0, nrow = nrow, ncol = ncol)
column <- round(runif(1, 1, ncol))
row <- round(runif(1, 1, nrow))
initial_infection[column, row] <- initial_infection[column, row] + round(runif(1, 0, 20))
initial_infection[column, row] <- initial_infection[column + round(runif(1, -10, 10)), row + round(runif(1, -10, 10))] + round(runif(1, 0, 20))
initial_infection[column, row] <- initial_infection[column + round(runif(1, -10, 10)), row + round(runif(1, -10, 10))] + round(runif(1, 0, 20))
initial_infection[column, row] <- initial_infection[column + round(runif(1, -10, 10)), row + round(runif(1, -10, 10))] + round(runif(1, 0, 20))
initial_infection[column, row] <- initial_infection[column + round(runif(1, -10, 10)), row + round(runif(1, -10, 10))] + round(runif(1, 0, 20))

## convert to raster for saving
host <- raster(host, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
all_plants <- raster(all_plants, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)
initial_infection <- raster(initial_infection, xmn = 0, ymn = 0, xmx = xmax, ymx = ymax)

writeRaster(host, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/host.tif", overwrite = TRUE)
writeRaster(all_plants, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/all_plants.tif", overwrite = TRUE)
writeRaster(initial_infection, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/initial_infections.tif", overwrite = TRUE)


## reproductive rate = 0.3 distance scale = 45 with management
infected_file = "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/initial_infections.tif"
host_file = "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/host.tif"
total_plants_file = "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/all_plants.tif"
temperature_file = ""
temperature_coefficient_file = ""
precipitation_coefficient_file =""
use_lethal_temperature = FALSE
temp = FALSE
precip = FALSE
season_month_start = 1
season_month_end = 12
time_step = "month"
start_time = 2012
end_time = 2017
dispersal_kern = "cauchy"
percent_short_distance_dispersal = 1.0
short_distance_scale = 38
long_distance_scale = 0.0
lethal_temperature = 0
lethal_temperature_month = 1
wind_dir = "NONE"
kappa = 0
random_seed = 42
reproductive_rate = 0.8
treatments_file = ""
treatment_years = c(0)
management = FALSE
mortality_on = FALSE
mortality_rate = 0
mortality_time_lag = 0

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


ncol = 250
nrow = 250
number = nrow *ncol
xmax = nrow * 100
ymax = ncol * 100

infected_stack <- stack(raster(data$infected_before_treatment[[1]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
for (i in 2: length(data$infected_before_treatment)) {
  infected_stack <- stack(infected_stack, raster(data$infected_before_treatment[[i]], xmn = 0, ymn = 0, xmx = xmax, ymx = ymax))
}

writeRaster(infected_stack, "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/reproductive_rate_0_3_distance_scale_45/infected_years_with_management.tif", overwrite = TRUE)

