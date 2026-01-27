## Load data and process
library(data.table)

# Case data frame
df = fread("./data/lassa_states_endemic.tsv")

# Restrict variables
df = df[, .(adminName, adminCode, date, epiWeek, year, logpop, polyid, region, y)]

# Set date column to date format
df[, date := as.Date(date)]

# Get set of admin codes for use
admins = df[, unique(adminCode)]

# Climate data
climate = readRDS("./data/climate.rds")

# Merge datasets
df = merge(df, climate, by = c("adminCode", "date"), all = TRUE)

# Exclude non-focal states
df = df[adminCode %in% admins]

# Up to latest date for data
latest_date = df[!is.na(y), max(date)]
df = df[date <= latest_date]

# TODO TEMP: only use earlier case line when multiples are there
df[, dummy := 1:.N, by = .(adminCode, date)]
df = df[dummy == 1]
df$dummy = NULL

# Fill in polyid
df[, polyid := as.numeric(stringr::str_sub(adminCode, 3))]

# Fill in region
regions = df[!is.na(region), .(region = unique(region)), by = adminCode]
df[regions, region := i.region]

# Fill in adminName
adminNames = df[!is.na(adminName), .(adminName = unique(adminName)), by = adminCode]
df[adminNames, adminName := i.adminName]

# # Exclude missing (forecast) datapoints
# df = df[!is.na(df$y), ]


# Neighbour graph file
nb_loc = "./data/states_nbmatrix"


# Result
return (list(df = df, nb_loc = nb_loc))
