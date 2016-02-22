#Learning Spatial Distribution Modelling
#install.packages(c("raster", "rgdal", "dismo", "rJava"))

require(dismo)
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file, header=TRUE, sep=",")
head(bradypus)

#Load the workspace downlaoded from GBIF
load("PumaRecords.rda")

#Puma <- gbif("puma", "concolor*", geo=FALSE)
#Pumageo <- subset(Puma, !is.na(lon) & !is.na(lat))

#save(Pumageo, file="PumaRecords.rda")

dim(Pumageo)

#Plots but takes a while - some are clearly wrong - they're in the ocean
require(leaflet)
m <- leaflet(N) %>% addTiles() 
m %>% setView(-72.519848, 42.373157, zoom = 1)
m %>% addMarkers(~ lon, ~ lat, clusterOptions = markerClusterOptions(), popup=as.character(N$lon))

#Would be nice to have leaflet display the lat and lon data for each marker.


New<-subset(Pumageo, subset=lon<(-40)&lat>(-65))

dups2 <- duplicated(New[, c("lon", "lat")])
Ne <- New[!dups2, ]

N<-Ne[!(Ne$lon=="-48.86401" | Ne$lon=="-62.54101" | Ne$lon=="-61.98063" | Ne$lon=="-61.9617" | Ne$lon=="-79.95644" | Ne$lon=="-79.95622"),]
#-48.86401
#-62.54101
#-61.98063
#-61.9617
#-79.95644
#-79.95622




require(ggmap)

#Find America
America<-get_map(location = c(lon = -80, lat = 10), source="google", zoom=2)
#Save the map
AmericaMap<-ggmap(America, extent = "device", legend = "topright")


#Contour plot on a static map
AmericaMap + 
  stat_density2d(aes(x = lon, y = lat,
                     fill = ..level..,alpha=..level..), geom = "polygon", data = N) +
  scale_fill_gradient(low = "black", high = "red") + ggtitle("Puma Density")




dim(Ne)

library(sp)
library(maptools)
coordinates(Pumageo) <- ~lon+lat
data(wrld_simpl)
crs(Pumageo) <- crs(wrld_simpl)
ovr <- over(Pumageo, wrld_simpl)
cntr <- ovr$NAME
i <- which(is.na(cntr))

j <- which(cntr != Pumageo$country)
# for the mismatches, bind the country names of the polygons and points
cbind(cntr, Pumageo$country)[j,]

plot(Pumageo)
plot(wrld_simpl, add=T, border="blue", lwd=2)
points(Pumageo[j, ], col="red", pch=20, cex=2)
