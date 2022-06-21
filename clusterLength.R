#### 1. Install libraries ####

# install.packages("extRemes")
# install.packages("lubridate")
# install.packages("ncdf4")
# install.packages("PCICt")
# install.packages("tidyverse")

#### 2. Load libraries ####

library(extRemes)  # for decluster
library(lubridate) # for as.Date
library(ncdf4)     # for nc_open
library(PCICt)     # for as.PCICt
library(tidyverse) # for case_when

#### 3. Remove all objects and set the working directory ####

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### 4. Functions ####

# define function to calculate season dates
season.dates <- function(start.date, end.date) {
  season.dates <- as.Date(start.date)
  while (tail(season.dates, n=1) < end.date){
    tmp <- ymd(as.Date(tail(season.dates, n=1))) %m+% months(3)
    season.dates <- c(season.dates, tmp)
  }
  season.dates[length(season.dates)] <- end.date
  season.dates <- as.POSIXct(season.dates, tz ="UTC", format = "%Y-%m-%d")
  attr(season.dates, "tzone") <- "UTC"
  return(season.dates)
}

# define function to calculate cluster lengths
cluster.lengths <- function(data, threshold, step.size) {
  e <- extremalindex(data, threshold = threshold)
  n.clusters <- e[[2]]
  run.length <- e[[3]]
  d <- decluster(data, threshold = threshold, method = "runs", r = run.length)
  # count the time between the first and last exceedance in each cluster
  cluster.id <- attributes(d)$clusters
  exceedance.index <- which(data > threshold)
  cluster.lengths <- integer(n.clusters)
  for (i in 1:n.clusters){
    index.min <- min(which(cluster.id==i))
    index.max <- max(which(cluster.id==i))
    cluster.lengths[i] <- (exceedance.index[index.max] - 
                             exceedance.index[index.min] + 1)*step.size
  }
  return(cluster.lengths)
}

# define function to calculate seasonal cluster lengths
seasonal.cluster.lengths <- function(grid, datetime, dates, data, event, 
                                     threshold, step.size, offset) {
  seasonal.cluster.lengths <- data.frame(length = double(), season = integer())
  n.seasons <- length(dates) - 1
  for (i in 1:n.seasons){
    # subset data
    if (grid == TRUE){
      start.time <- which(as.character(datetime) == dates[i])
      end.time <- which(as.character(datetime) == dates[i+1])
    } else {
      start.time <- which(datetime == dates[i])
      end.time <- which(datetime == dates[i+1])
    }
    subset.data <- data[start.time:(end.time-1),]
    # calculate cluster lengths
    if (length(which(subset.data[[event]] > threshold)) > 1){
      # error occurs when only 0 or 1 event exceeds threshold
      # warning occurs when only 2 events exceed threshold
      cluster.lengths <- cluster.lengths(subset.data[[event]], threshold, 
                                         step.size)
      tmp.data <- data.frame(length = cluster.lengths, season = (i+offset)%%4)
    } else {
      tmp.data <- data.frame(length = NA, season = (i+offset)%%4)
    }
    seasonal.cluster.lengths <- rbind(seasonal.cluster.lengths, tmp.data)
  }
  return(seasonal.cluster.lengths)
}

#### 5. Definitions ####

# ramp event = maximum change in capacity factor, within a given time window
# drought event = percentage of time capacity factor below capacity threshold, 
# within a given time window

#### 6. Parameters ####

window.size <- 12 # hours
ramp.threshold <- 0.9
drought.threshold <- 0.9
capacity.threshold <- 0.05

#### 7. Read power curve and its corrections ####

file.name <- 'data/power_curve/IEA_curve.csv'
IEA.curve <- read.csv(file.name, header = T)
file.name <- 'data/power_curve/IEA_TIC.csv'
IEA.TIC <- read.csv(file.name, header = T)

#### 8. Set modelled and observed data ####

# modelled data:
modelled <- 0
while ((modelled != 1) && (modelled != 2)){
  modelled <- readline(prompt="Which modelled data would you like to use?
                       Enter 1 for EURO-CORDEX or 2 for UKCP18: ")
}
if (modelled == 1){
  modelled <- "CORDEX"
} else {
  modelled <- "UKCP18"
}

# observed data:
observed <- 0
while ((observed != 1) && (observed != 2)){
  observed <- readline(prompt="Which observed data would you like to use?
                       Enter 1 for Europlatform or 2 for K13: ")
}
if (observed == 1){
  observed <- "Europlatform"
} else {
  observed <- "K13"
}

if (modelled == "CORDEX"){
  modelled.long <- "EURO-CORDEX"
} else {
  modelled.long <- "UKCP18"
}

comparison <- case_when(
  modelled == 'CORDEX' && observed == 'Europlatform' ~ "cordex_epf",
  modelled == 'CORDEX' && observed == 'K13' ~ "cordex_k13",
  modelled == 'UKCP18' && observed == 'Europlatform' ~ "ukcp18_epf",
  modelled == 'UKCP18' && observed == 'K13' ~ "ukcp18_k13"
)

dir.create(paste0('plots/', comparison))

#### 9. Create and save data frame of modelled data ####

# create data frame of modelled data time and wind speed 
# at nearest neighbour to observed data site
file.name <- case_when(
  comparison == 'cordex_epf' ~ 
    paste0('sfcWind_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_',
           'MOHC-HadREM3-GA7-05_v1_3hr_Europlatform.nc'),
  comparison == 'cordex_k13' ~ 
    paste0('sfcWind_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_',
           'MOHC-HadREM3-GA7-05_v1_3hr_K13.nc'),
  comparison == 'ukcp18_epf' ~ 
    'sfcWind_rcp85_land-cpm_uk_2.2km_01_1hr_Europlatform.nc',
  comparison == 'ukcp18_k13' ~ 
    'sfcWind_rcp85_land-cpm_uk_2.2km_01_1hr_K13.nc'
)

if (modelled == "CORDEX"){
  tmp <- nc_open(paste0('data/euro-cordex/', file.name))
  grid.time <- ncvar_get(tmp, varid = 'time') # days since 1949-12-01
} else {
  tmp <- nc_open(paste0('data/ukcp18/', file.name))
  grid.time <- ncvar_get(tmp, varid = 'time') # hours since 1970-01-01
}

grid.speed <- ncvar_get(tmp, varid = 'sfcWind')
nc_close(tmp)
grid.data <- data.frame(time = grid.time, speed = grid.speed)
rm(tmp, grid.time, grid.speed)
head(grid.data)

# calculate capacity factor for each wind speed
grid.data$capacity <- 
  approx(IEA.curve$WS, IEA.curve$CF, xout = grid.data$speed, rule = 2)$y*
  approx(IEA.TIC$WS, IEA.TIC$TIC, xout = grid.data$speed, rule = 2)$y
# set max capacity factor = 1
grid.data$capacity <- pmin(grid.data$capacity, 1)
head(grid.data)

# calculate ramp event
grid.step <- as.numeric(grid.data$time[2] - grid.data$time[1])*24 # hours
n.elements <- window.size/grid.step
n.dates <- nrow(grid.data)
grid.data$ramp <- integer(n.dates)
for (i in n.elements:n.dates){
  grid.data$ramp[i] <- 
    max(grid.data$capacity[(i - (n.elements - 1)):i]) - 
    min(grid.data$capacity[(i - (n.elements - 1)):i])
}
head(grid.data)

# calculate drought event
grid.data$drought <- integer(n.dates)
for (i in n.elements:n.dates){
  capacity.factors <- grid.data$capacity[(i - (n.elements - 1)):i]
  grid.data$drought[i] <-sum(capacity.factors < capacity.threshold)/n.elements
}
head(grid.data)

# save modelled data
file.name <- paste0('data/rds/', comparison, '.rds')
saveRDS(grid.data, file.name)

#### 10. Create and save data frame of observed data ####

# create data frame of observed data time and wind speed
if (observed == "Europlatform"){
  file.name <- 'data/obs_knmi/europlt_processed.csv'
} else {
  file.name <- 'data/obs_knmi/k13_processed.csv'
}
site.data <- read.csv(file.name, skip = 9, header = T)
# convert Datetime from type character to type double
site.data$Datetime <- as.POSIXct(site.data$Datetime, format="%d/%m/%Y %H:%M", 
                                 tz="UTC")
range(site.data$Datetime) # "2003-04-01 00:10:00 UTC" "2022-02-01 00:00:00 UTC"
# when Flag==0, replace FF_10M_10 with NA
n.dates <- nrow(site.data)
for (i in 1:n.dates){
  if (site.data$Flag[i] == 0){
    site.data$FF_10M_10[i] <- NA
  }
}
# keep only datetime and speed columns
site.data <- site.data[c("Datetime", "FF_10M_10")]
# rename columns
colnames(site.data) <- c("datetime", "speed")
head(site.data)

# calculate capacity factor for each wind speed
site.data$capacity <- 
  approx(IEA.curve$WS, IEA.curve$CF, xout = site.data$speed, rule = 2)$y*
  approx(IEA.TIC$WS, IEA.TIC$TIC, xout = site.data$speed, rule = 2)$y
# set max capacity factor = 1
site.data$capacity <- pmin(site.data$capacity, 1)
head(site.data)

# calculate ramp event
site.step <- as.numeric(site.data$datetime[2] - site.data$datetime[1])/60 # hours
n.elements <- window.size/site.step
site.data$ramp <- integer(n.dates)
for (i in n.elements:n.dates){
  site.data$ramp[i] <- 
    max(site.data$capacity[(i - (n.elements - 1)):i], na.rm = TRUE) - 
    min(site.data$capacity[(i - (n.elements - 1)):i], na.rm = TRUE)
}
sum(is.na(site.data$ramp)) # 0
# replace -Inf with 0 (set min ramp = 0)
site.data$ramp <- pmax(site.data$ramp, 0)
sum(is.infinite(site.data$ramp)) # 0
site.data[70:75,]

# calculate drought event
site.data$drought <- integer(n.dates)
for (i in n.elements:n.dates){
  capacity.factors <- site.data$capacity[(i - (n.elements - 1)):i]
  site.data$drought[i] <-sum(capacity.factors < capacity.threshold, 
                             na.rm = TRUE)/n.elements
}
sum(is.na(site.data$drought)) # 0
site.data[70:75,]

# save observed data
file.name <- paste0('data/rds/', observed, '.rds')
saveRDS(site.data, file.name)

#### 11. Read modelled data ####

file.name <- paste0('data/rds/', comparison, '.rds')
grid.data <- readRDS(file.name)

#### 12. Read observed data ####

file.name <- paste0('data/rds/', observed, '.rds')
site.data <- readRDS(file.name)

# # modify k13 data so that end date = '2019-01-01'
# end.date <- '2019-01-01'
# end.row <- which(site.data$datetime == end.date)
# site.data <- site.data[1:end.row,]
# tail(site.data)

#### 13. Calculate modelled data datetime ####

if (modelled == "CORDEX"){
  grid.datetime <- as.PCICt(grid.data$time*3600*24, origin = '1949-12-01', 
                            cal="360_day")
  range(grid.datetime) # "1980-01-01" "2005-12-01"
} else {
  grid.datetime <- as.PCICt(grid.data$time*3600, origin = '1970-01-01', 
                            cal="360_day")
  range(grid.datetime) # "1980-12-01 00:00:00" "2000-11-30 23:00:00"
}

#### 14. Calculate modelled data seasons ####

if (modelled == "CORDEX"){
  start.date <- '1979-12-01'
  end.date <- '2005-12-01'
  grid.dates <- season.dates(start.date, end.date)
  grid.dates[1] <- as.Date('1980-01-01')
  grid.dates <- format(as.Date(grid.dates), "%Y-%m-%d %H:%M:%S")
} else {
  start.date <- '1980-12-01'
  end.date <- '2000-11-30'
  grid.dates <- season.dates(start.date, end.date)
  grid.dates <- format(as.Date(grid.dates), "%Y-%m-%d %H:%M:%S")
  grid.dates[length(grid.dates)] <- format(as_datetime("2000-11-30 23:00:00"), 
                                           "%Y-%m-%d %H:%M:%S")
}

# first season = winter => add 2 to season to make winter = third season
grid.offset <- 2

#### 15. Calculate modelled data cluster lengths ####

grid.step <- as.numeric(grid.data$time[2] - grid.data$time[1])*24 # hours

# ramp events
grid.ramps <- cluster.lengths(grid.data$ramp, ramp.threshold, grid.step)

# drought events
grid.droughts <- cluster.lengths(grid.data$drought, drought.threshold, 
                                 grid.step)

#### 16. Calculate modelled data seasonal cluster lengths ####

# ramp events
grid.ramp <- seasonal.cluster.lengths(TRUE, grid.datetime, grid.dates, 
                                      grid.data, event=4, ramp.threshold, 
                                      grid.step, grid.offset)
head(grid.ramp)

# drought events
grid.drought <- seasonal.cluster.lengths(TRUE, grid.datetime, grid.dates, 
                                         grid.data, event=5, drought.threshold, 
                                         grid.step, grid.offset)
head(grid.drought)

#### 17. Calculate observed data seasons ####

if (observed == "Europlatform"){
  end.date <- '2022-02-01'
} else {
  end.date <- '2022-02-01'
  # end.date <- '2019-01-01'
}

start.date <- '2003-03-01'
site.dates <- season.dates(start.date, end.date)
site.dates[1] <- as.POSIXct("2003-04-01 00:10:00", tz ="UTC", 
                            format = "%Y-%m-%d %H:%M:%S")

# first season = spring => add 3 to season to make spring = fourth season
site.offset <- 3

#### 18. Calculate observed data cluster lengths ####

site.step <- as.numeric(site.data$datetime[2] - site.data$datetime[1])/60 # hours

# ramp events
site.ramps <- cluster.lengths(site.data$ramp, ramp.threshold, site.step)

# drought events
site.droughts <- cluster.lengths(site.data$drought, drought.threshold, 
                                 site.step)

#### 19. Calculate observed data seasonal cluster lengths ####

# ramp events
site.ramp <- seasonal.cluster.lengths(FALSE, site.data$datetime, site.dates, 
                                      site.data, event=4, ramp.threshold, 
                                      site.step, site.offset)
head(site.ramp)

# drought events
site.drought <- seasonal.cluster.lengths(FALSE, site.data$datetime, site.dates, 
                                         site.data, event=5, drought.threshold, 
                                         site.step, site.offset)
head(site.drought)

#### 20. Boxplots comparing modelled and observed cluster lengths ####

# ramp events with outliers
file.name <- paste0('plots/', comparison, '/boxplot_ramp_outlier.png')
png(file.name)
boxplot(grid.ramps, site.ramps, notch = TRUE, horizontal = FALSE, 
        outline = TRUE, ylab = "Ramp length (hours)", 
        names = c(modelled.long, observed))
dev.off()

# ramp events without outliers
file.name <- paste0('plots/', comparison, '/boxplot_ramp.png')
png(file.name)
boxplot(grid.ramps, site.ramps, notch = TRUE, horizontal = FALSE, 
        outline = FALSE, ylab = "Ramp length (hours)", 
        names = c(modelled.long, observed))
dev.off()

# drought events with outliers
file.name <- paste0('plots/', comparison, '/boxplot_drought_outlier.png')
png(file.name)
boxplot(grid.droughts, site.droughts, notch = TRUE, horizontal = FALSE, 
        outline = TRUE, ylab = "Drought length (hours)", 
        names = c(modelled.long, observed))
dev.off()

# drought events without outliers
file.name <- paste0('plots/', comparison, '/boxplot_drought.png')
png(file.name)
boxplot(grid.droughts, site.droughts, notch = TRUE, horizontal = FALSE, 
        outline = FALSE, ylab = "Drought length (hours)", 
        names = c(modelled.long, observed))
dev.off()

#### 21. Boxplots comparing modelled and observed seasonal cluster lengths ####

# create data frame with modelled and observed seasonal ramp cluster lengths
site.ramp$season <- site.ramp$season + 0.5
both.ramp <- rbind(grid.ramp, site.ramp)

# ramp events with outliers
file.name <- paste0('plots/', comparison, '/boxplot_ramp_seasonal_outlier.png')
png(file.name, width = 811, height = 505)
boxplot(length~season, data=both.ramp, notch = TRUE, horizontal = FALSE, 
        xlab = "", ylab = "Ramp length (hours)", outline = TRUE, 
        axes = FALSE)
axis(1, at = 1:8, cex.axis = 1, mgp=c(3,2,0), 
     labels = c(paste0("Spring\n", modelled), paste0("Spring\n", observed), 
                paste0("Summer\n", modelled), paste0("Summer\n", observed), 
                paste0("Autumn\n", modelled), paste0("Autumn\n", observed), 
                paste0("Winter\n", modelled), paste0("Winter\n", observed)))
axis(2)
box()
dev.off()

# ramp events without outliers
file.name <- paste0('plots/', comparison, '/boxplot_ramp_seasonal.png')
png(file.name, width = 811, height = 505)
boxplot(length~season, data=both.ramp, notch = TRUE, horizontal = FALSE, 
        xlab = "", ylab = "Ramp length (hours)", outline = FALSE, 
        axes = FALSE)
axis(1, at = 1:8, cex.axis = 1, mgp=c(3,2,0), 
     labels = c(paste0("Spring\n", modelled), paste0("Spring\n", observed), 
                paste0("Summer\n", modelled), paste0("Summer\n", observed), 
                paste0("Autumn\n", modelled), paste0("Autumn\n", observed), 
                paste0("Winter\n", modelled), paste0("Winter\n", observed)))
axis(2)
box()
dev.off()

# create data frame with modelled and observed seasonal drought cluster lengths
site.drought$season <- site.drought$season + 0.5
both.drought <- rbind(grid.drought, site.drought)

# drought events with outliers
file.name <- paste0('plots/', comparison, 
                   '/boxplot_drought_seasonal_outlier.png')
png(file.name, width = 811, height = 505)
boxplot(length~season, data=both.drought, notch = TRUE, horizontal = FALSE, 
        xlab = "", ylab = "Drought length (hours)", outline = TRUE, 
        axes = FALSE)
axis(1, at = 1:8, cex.axis = 1, mgp=c(3,2,0), 
     labels = c(paste0("Spring\n", modelled), paste0("Spring\n", observed), 
                paste0("Summer\n", modelled), paste0("Summer\n", observed), 
                paste0("Autumn\n", modelled), paste0("Autumn\n", observed), 
                paste0("Winter\n", modelled), paste0("Winter\n", observed)))
axis(2)
box()
dev.off()

# drought events without outliers
file.name <- paste0('plots/', comparison, '/boxplot_drought_seasonal.png')
png(file.name, width = 811, height = 505)
boxplot(length~season, data=both.drought, notch = TRUE, horizontal = FALSE, 
        xlab = "", ylab = "Drought length (hours)", outline = FALSE, 
        axes = FALSE)
axis(1, at = 1:8, cex.axis = 1, mgp=c(3,2,0), 
     labels = c(paste0("Spring\n", modelled), paste0("Spring\n", observed), 
                paste0("Summer\n", modelled), paste0("Summer\n", observed), 
                paste0("Autumn\n", modelled), paste0("Autumn\n", observed), 
                paste0("Winter\n", modelled), paste0("Winter\n", observed)))
axis(2)
box()
dev.off()

#### 22. Find outlier in cordex_epf data ramp events cluster lengths ####

which(grid.ramp$length > 800) # 601
nrow(grid.ramp) # 635
tail(grid.ramp, 40)
range(grid.datetime) # "1980-01-01" "2005-12-01"
# create subset of modelled data in spring 2004
start.time <- which(as.character(grid.datetime) == "2004-03-01 00:00:00")
end.time <- which(as.character(grid.datetime) == "2004-06-01 00:00:00")
subset.data <- grid.data[start.time:end.time,]
head(subset.data)
# plot time series of modelled data ramp events in spring 2004
file.name <- paste0('plots/', comparison, '/timeseries_ramp_spring2004.png')
png(file.name, width = 811, height = 505)
plot(subset.data$ramp, type="o")
abline(h=ramp.threshold, col='red')
dev.off()
# calculate run length of modelled data ramp clusters in spring 2004
e <- extremalindex(subset.data$ramp, threshold = ramp.threshold)
run.length <- e[[3]] # 318
