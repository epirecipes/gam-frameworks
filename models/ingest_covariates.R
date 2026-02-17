library(data.table)

# precipitation in mm

load = function(path)
{
    files = list.files(path, full.names = TRUE)
    data = rbindlist(lapply(files, fread))
    
    data = data[, .(adminCode, date, value)]
    names(data)[3] = "precip"
    return (data[])
}

precip_path = "~/Documents/lassa_sentinel/output/covariates/precip_chirps/daily"
precip = load(precip_path)
precip = precip[date >= "2010-01-01"]
# saveRDS(precip, "./data/precip.rds")


# SPI (note: only available once per month)

spi_path = "~/Documents/lassa_sentinel/output/covariates/spi_chirps/daily/spi(NN)/spi(NN)_chirps_admin1_daily.csv"

result = list()
for (n in c(1, 3, 6, 12)) {
    path = stringr::str_replace_all(spi_path, "\\(NN\\)", as.character(n))
    data = fread(path)
    data = data[date >= "2010-01-01", .(adminCode, date, value)]
    names(data)[3] = paste0("spi", n)
    result[[length(result) + 1]] = data
}
spi = merge(result[[1]], result[[2]], by = c("adminCode", "date"), all = TRUE)
spi = merge(spi, result[[3]], by = c("adminCode", "date"), all = TRUE)
spi = merge(spi, result[[4]], by = c("adminCode", "date"), all = TRUE)
# saveRDS(spi, "./data/spi.rds")


# EVI

evi_path = "~/Documents/lassa_sentinel/output/covariates/evi_modis/daily"
files = list.files(evi_path, full.names = TRUE)
evi = lapply(files, fread)
evi = rbindlist(evi)

evi = evi[date >= "2010-01-01", .(adminCode, date, evi = value)]
# saveRDS(evi, "./data/evi.rds")

# temperature

temp_path = "~/Documents/lassa_sentinel/output/covariates/temp_modis/daily"
files = list.files(temp_path, full.names = TRUE)
temp1 = lapply(files, fread)
temp1 = rbindlist(temp1)

temp1 = temp1[date >= "2010-01-01", .(adminCode, date, temp1 = value)]
# saveRDS(temp, "./data/temp.rds")

# alt. temperature
temp2_path = "~/Documents/lassa_sentinel/output/covariates/temp_era5land/monthly"
files = list.files(temp2_path, full.names = TRUE)
temp2 = lapply(files, fread)
temp2 = rbindlist(temp2)

temp2 = temp2[date >= "2010-01-01", .(adminCode, date, temp = value)]
temp2[, adminCode := stringr::str_sub(adminCode, 1, 5)]
# Aggregate to state by taking mean over each LGA. Ideally this would be a 
# statewide temperature instead.
temp2 = temp2[, .(temp2 = mean(temp)), by = .(adminCode, date)]





# merge all climate

data = precip |>
  merge(spi, by = c("adminCode", "date"), all = TRUE) |>
  merge(evi, by = c("adminCode", "date"), all = TRUE) |>
  merge(temp1, by = c("adminCode", "date"), all = TRUE) |>
  merge(temp2, by = c("adminCode", "date"), all = TRUE)

saveRDS(data, "./data/climate.rds")

data[, sum(is.na(precip))]
data[, sum(is.na(spi1))]
data[, sum(is.na(spi3))]
data[, sum(is.na(spi6))]
data[, sum(is.na(spi12))]
data[, sum(is.na(evi))]
data[, sum(is.na(temp1))]
data[, sum(is.na(temp2))]


# take a look
library(ggplot2)
ggplot(data[adminCode == "NG020"]) +
    geom_line(aes(x = date, y = temp1, colour = "temp1"))
