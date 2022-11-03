## AWESOME PHYSICS MODEL (APM)
require(remotes)
remotes::install_github("robertladwig/LakeModelR")
remotes::install_github("aemon-j/gotmtools")

library(LakeModelR)

aquatic_sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |>
  dplyr::filter(aquatics == 1)

lake_sites <- aquatic_sites %>%
  filter(field_site_subtype == 'Lake')

require(neon4cast)

vars = c("air_temperature", "air_pressure", "relative_humidity", "surface_downwelling_longwave_flux_in_air", "surface_downwelling_shortwave_flux_in_air",
         "eastward_wind", "northward_wind")

#read in the targets data
targets <- read_csv('https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz')

targets <- targets %>%
  filter(site_id %in% lake_sites$field_site_id)


# New forecast only available at 5am UTC the next day 

forecast_date <- Sys.Date() 
noaa_date <- forecast_date - days(2)

df_future <- neon4cast::noaa_stage2()

noaa_future_allvars <- df_future |> 
  dplyr::filter(reference_datetime == noaa_date,
                datetime >= forecast_date,
                site_id %in% lake_sites$field_site_id,
                variable %in% vars) |> 
  dplyr::collect()

# past stacked weather
df_past <- neon4cast::noaa_stage3()

#Other variable names can be found at https://projects.ecoforecast.org/neon4cast-docs/Shared-Forecast-Drivers.html#stage-2

noaa_past_allvars <- df_past |> 
  dplyr::filter(site_id %in% lake_sites$field_site_id,
                datetime >= ymd('2017-01-01'),
                variable %in% vars) |> 
  dplyr::collect()

require(gotmtools)


future_weather = noaa_future_allvars %>% 
  group_by(datetime, site_id, variable) %>% 
  summarise(prediction = mean(prediction, na.rm = T),
            longitude = first(longitude),
            latitude = first(latitude)) %>% 
  ungroup %>% 
  pivot_wider(values_from = prediction, names_from = variable)

names(future_weather)

future_weather_onelake = future_weather %>% 
  filter(site_id == "BARC")
# select(datetime,
#        air_temperature,
#        relative_humidity,
#        surface_downwelling_shortwave_flux_in_air,
#        latitude,
#        longitude,
#        surface_downwelling_longwave_flux_in_air,
#        site_id)


future_weather_onelake$surface_downwelling_shortwave_flux_in_air[793] = 0

cloud_cover = calc_cc(date = future_weather_onelake$datetime,
                      airt = future_weather_onelake$air_temperature - 273.15,
                      relh = future_weather_onelake$relative_humidity,
                      swr = zoo::na.approx(future_weather_onelake$surface_downwelling_shortwave_flux_in_air),
                      lat = mean(future_weather_onelake$latitude),
                      lon = mean(future_weather_onelake$longitude),
                      elev = 10)

future_weather_onelake$ea = future_weather_onelake$relative_humidity * 10^(9.28603523 - 2322.37885/(future_weather_onelake$air_temperature))

future_weather_onelake = future_weather_onelake %>% 
  mutate(wind = sqrt(eastward_wind^2 + northward_wind^2))

future_weather_onelake$cc = cloud_cover
daily_meteo = future_weather_onelake
sigma = 5.67*10^(-8)
emissivity = 0.97
eps = 0.97

targets %>%
  filter(variable == 'temperature', site_id == "BARC") %>% 
  slice_max(datetime)

df_xiao = matrix(NA, nrow =  nrow(future_weather_onelake), ncol = 10)


for (j in 1: ncol(df)){
  temp = 23.3
  for (n in 1:nrow(future_weather_onelake)) {
    Q <- (
      longwave(cc = daily_meteo[n, "cc"], sigma = sigma, Tair = (daily_meteo[n, "air_temperature"] - 273.15) %>% pull(air_temperature), ea = daily_meteo[n, "ea"], emissivity = emissivity, Jlw = daily_meteo[n, "surface_downwelling_longwave_flux_in_air"]) +
        backscattering(emissivity = emissivity, sigma = sigma, Twater = temp[n], eps = eps) +
        latent(Tair = (daily_meteo[n, "air_temperature"] - 273.15) %>% pull(air_temperature), Twater = temp[n], Uw = daily_meteo[n, "wind"] %>% pull(wind), p2 = 1, pa = daily_meteo[n, "air_pressure"] %>% pull(air_pressure), ea=daily_meteo[n, "ea"] %>% pull(ea), RH = daily_meteo[n, "relative_humidity"] %>% pull(relative_humidity), A = 0.13 * 1000000, Cd = 0.0013) +
        sensible(Tair = (daily_meteo[n, "air_temperature"] - 273.15) %>% pull(air_temperature), Twater = temp[n], Uw = daily_meteo[n, "wind"] %>% pull(wind), p2 = 1, pa = daily_meteo[n, "air_pressure"] %>% pull(air_pressure), ea=daily_meteo[n, "ea"] %>% pull(ea), RH = daily_meteo[n, "relative_humidity"] %>% pull(relative_humidity), A = 0.13 * 1000000, Cd = 0.0013))
    
    H =  (1- 0.9) * (daily_meteo[n, "surface_downwelling_shortwave_flux_in_air"])
    
    temp = append(temp,as.numeric(temp[n] + (Q + H)/(4184 * calc_dens(temp[n])) * 3600 + rnorm(mean = 0, n = 1, sd = 0.005)))
    
    print(n)
  }
  df_xiao[, j] = temp[-1]
}

apply(df_xiao, 1, sd) %>% plot()
apply(df_xiao, 1, mean) %>% plot()

pred = daily_meteo %>% select(datetime, site_id) %>% 
  mutate(family = "normal",
         variable = "temperature",
         sigma = apply(df_xiao, 1, sd),
         mu = apply(df_xiao, 1, mean),
         model_id = "GLEON_physics"
  ) %>% 
  pivot_longer(names_to = "parameter",
               values_to = "prediction",
               cols = c(sigma, mu)) %>% 
  mutate(reference_datetime = Sys.Date(),
         datetime = as.Date(datetime)) %>% 
  group_by(datetime, variable, parameter, family, model_id, reference_datetime, site_id) %>% 
  summarise(
    prediction = mean(prediction)
  ) %>% 
  ungroup

temp_lm_forecast_EFI = pred

# Start by writing the forecast to file
theme <- 'aquatics'
date <- temp_lm_forecast_EFI$reference_datetime[1]
forecast_name_1 <- paste0(temp_lm_forecast_EFI$model_id[1], ".csv")
forecast_file_1 <- paste(theme, date, forecast_name_1, sep = '-')
forecast_file_1

write_csv(temp_lm_forecast_EFI, forecast_file_1)

neon4cast::forecast_output_validator(forecast_file_1)

neon4cast::submit(forecast_file = forecast_file_1,
                  ask = F) # if ask = T (default), it will produce a pop-up box asking if you

message('AHOI!')