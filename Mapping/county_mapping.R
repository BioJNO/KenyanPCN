# Step 1: load the packages you will need ------------------------------------

# To install a package
install.packages("ggplot2")
install.packages("maptools")

# load required libraries
# TO DO: remove unused
library(ggplot2)
library(rgdal)
library(maptools)
library(ggmap)
library(dplyr)
library(reshape2)
library(gridExtra)
library(grid)

# Step 2: import shapefile data ----------------------------------------------
# Download a shapefile for kenya and save it in the directory (folder) you
# are working in. 
# The shapefile we are using is available here: 
# https://africaopendata.org/dataset/kenya-counties-shapefile

# Read in shapefile data.
kenya = readOGR(dsn="KenyanCountyShapefile",
                layer="County")
kenya@data$id = rownames(kenya@data)
kenya.points = fortify(kenya, region="id")

# store the shapefile data as a data frame (table)
kenya.df = left_join(kenya.points, kenya@data, by="id")

# subset shapefile data for particular counties
nyandarua.df <- kenya.df[kenya.df$COUNTY=="Nyandarua",]
keiyo.marakwet.df <- kenya.df[kenya.df$COUNTY=="Keiyo-Marakwet",]

# Step 3: import and merge your data -----------------------------------------
# Your data folders also need to be in the directory you're working from
# read in data tables for plotting
cyst_data <- read.csv("cyst_collection_data.csv")
metadata <- read.csv("sampling_metadata.csv")

# data tables can be merged so long as they share a column with the same
# heading which contains common values e.g. sample IDs.
# Note: this is case sensitive! 
all.data <- merge(cyst_data, metadata, by = "Sample.code")

# write out the merged file. 
write.csv(all.data, "collection_merged_with_metadata.csv")

# select rows in your dataset based on the value in a specified column.
nyandarua_cyst <- all.data[all.data$County.x=="Nyandarua",]
Taita.Taveta.cyst <- cyst_data[cyst_data$County=="Taita Taveta",]
Keiyo.Marakwet.cyst <- cyst_data[cyst_data$County=="Elgeyo Marakwet",]

# Add a column to your dataset which highlights counties by potato production
# give this column a default value of FALSE. This kind of value (TRUE/FALSE)
# is called a "Boolean" value. This means it is possible for it to hold one
# of two values. 
kenya.df$potato.production <- FALSE

# change the value of the "potato.production" column to "TRUE"
# based on the value of the "COUNTY" column. 
kenya.df$potato.production[kenya.df$COUNTY == "Nyandarua"] <- TRUE
kenya.df$potato.production[kenya.df$COUNTY == "Taita Taveta"] <- TRUE
kenya.df$potato.production[kenya.df$COUNTY == "Keiyo-Marakwet"] <- TRUE

# Step 4: plot a base map ----------------------------------------------------
# Plot a map of Kenya with counties coloured by the value of the 
# potato production we created above.
map <- ggplot() + geom_polygon(data = kenya.df,
                               aes(x = long,
                                   y = lat,
                                   group = group,
                                   fill = potato.production),
                               colour = "black")
map

# Save as a scalable vector graphic.
# svg is a text format you can view in inkscape.
# You can alter and save in inkscape (https://inkscape.org/)
# To include in a word document or powerpoint save as an .emf 
# Use vector graphics for presentations, especially posters, 
# to avoid losing resolution when the image is enlarged.
ggsave("Kenya_county_potato_production_map.svg", dpi=600, width=7, height =6)

# Plot the same map without axes. 
map <- ggplot() + geom_polygon(data = kenya.df,
                               aes(x = long,
                                   y = lat,
                                   group = group,
                                   fill = potato.production),
                               colour = "black") +
                               theme_void()

map

# Step 4: add points to your base map ----------------------------------------
# As long as you have latitude/longitude values you can add these points to
# your base map.
p1 = map + geom_point(mapping = aes(Longitude_2, Latitude_2),
                          color = "white",
                          data = cyst_data)
p1

# Save as a scalable vector graphic.
ggsave("Kenya_county_sampled_map.svg", dpi=600, width=7, height =6)

# Plot only one county -------------------------------------------------------
# As above except use a subset of the data frame of the shape file. 
nyandarua.map <- ggplot() + geom_polygon(data = nyandarua.df,
                               aes(x = long,
                                   y = lat,
                                   group = group),
                               colour = "black",
                               fill = "gray75")
nyandarau.map

# Plot all sample positions in a specific county
p1 = nyandarua.map + geom_point(mapping = aes(Longitude_2, Latitude_2),
                                color = "white",
                                data = nyandarua_cyst)
p1


# layer cyst count on top of sample locations
# colour by potato variety 
# size by number of cysts recovered
# TO DO: sort cbfriendly discrete colour scheme
p2 = p1 + geom_point(mapping = aes(Longitude_2, Latitude_2,
                                   colour = Potato.Variety.1,
                                   fill = Potato.Variety.1,
                                   group = 1,
                                   size = No_of_Cysts),
                     alpha = 0.8,
                     shape = 21,
                     data=nyandarua_cyst)
p2

