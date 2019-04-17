## Calibrate the 2 synthetic datasets using quantity, allocation, and configuration disagreement metrics
devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "feature/parallel_calibration", force = TRUE)
library(PoPS)
library(raster)

infected_years_file = "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/infected_years_rr0.7_sd21_extent200_single.tif"
num_iterations = 1000
start_reproductive_rate = 1.0
start_short_distance_scale = 40
sd_reproductive_rate = 0.2
sd_short_distance_scale = 1

infected_file = "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/initial_infections_single.tif"
host_file = "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/host.tif"
total_plants_file = "H:/Team Drives/APHIS  Projects/Papers/oddsratio vs disagreements/data/250/200/all_plants.tif"
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
end_time = 2014
dispersal_kern = "cauchy"
percent_short_distance_dispersal = 1.0
short_distance_scale = 21
long_distance_scale = 0.0
lethal_temperature = 0
lethal_temperature_month = 1
wind_dir = "NONE"
kappa = 0
random_seed = 42
reproductive_rate = 0.5
treatments_file = ""
treatment_years = c(0)
management = FALSE
mortality_on = FALSE
mortality_rate = 0
mortality_time_lag = 0

params_3year <- calibrate(infected_years_file, num_iterations, start_reproductive_rate, 
                          start_short_distance_scale, sd_reproductive_rate, sd_short_distance_scale,
                          infected_file, host_file, total_plants_file, reproductive_rate,
                          use_lethal_temperature, temp, precip, management, mortality_on,
                          temperature_file, temperature_coefficient_file, 
                          precipitation_coefficient_file, treatments_file,
                          season_month_start, season_month_end, time_step,
                          start_time, end_time, treatment_years,
                          dispersal_kern, percent_short_distance_dispersal,
                          short_distance_scale, long_distance_scale,
                          lethal_temperature, lethal_temperature_month,
                          mortality_rate, mortality_time_lag,
                          wind_dir, kappa)
