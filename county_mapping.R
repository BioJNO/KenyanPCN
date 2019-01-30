# Load required packages
# some of these packages we don't need
library(ggplot2)
library(rgdal)
library(maptools)
library(ggmap)
library(dplyr)
library(reshape2)
library(gridExtra)
library(grid)

# Plot a base R map of Kenya ----------------------------------------------
# Read shapefile data.
kenya = readOGR(dsn="KenyanCountyShapefile",
                layer="County")
kenya@data$id = rownames(kenya@data)
kenya.points = fortify(kenya, region="id")
kenya.df = left_join(kenya.points, kenya@data, by="id")

# Plot the base map
basemap = ggplot(kenya.df, aes(long, lat, group=group))
basemap = basemap + geom_polygon(fill="gray75") + theme(panel.background = element_rect(fill = "gray95"))
basemap = basemap + coord_equal()
basemap
