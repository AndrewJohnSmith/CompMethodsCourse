# Spatial analyses

#----------------------------------------------------
# Mantel test
#----------------------------------------------------
# Mantel's test of significance is evaluated via permutation procedures,
# in which the rows and/or columns of the distance matrices are randomly rearranged.

# Test the relationship between fruit and leaf consumption in primates using a Mantel test
library(vegan)

# Import the fruit data
data_fruit <- read.csv("data_fruit.csv")

# Remove the species label column
data_fruit2 = data_fruit[,-1]

# Convert the data to a Euclidian distance matrix
data_fruit2_dist = dist(data_fruit2)

# Import the leaves data
data_leaves <- read.csv("data_leaves.csv")

# Remove the species label column
data_leaves2 = data_leaves[,-1]

# Convert the data to a Euclidian distance matrix
data_leaves2_dist = dist(data_leaves2)

# Run the test

mantel01 = mantel(data_leaves2_dist, data_fruit2_dist, method="pearson", permutations=9999)

mantel01

#-----------------------------------------------------------
# Partial Mantel test
#-----------------------------------------------------------

# Partial Mantel statistic uses partial correlation conditioned on the third matrix.

# Now account for body mass

# Import the body mass data
data_mass <- read.csv("data_mass.csv")

# Remove the species label column
data_mass2 = data_mass[,-1]

# Convert the data to a Euclidian distance matrix
data_mass_dist = dist(data_mass2)

# Run the partial Mantel test

partial_mantel01 = mantel.partial(data_leaves2_dist, data_fruit2_dist, data_mass_dist, method="pearson", permutations=9999)

partial_mantel01

##############################################
# Run a Mantel test examining the relationship between community species composition and environmental factors
# Use the  community_species_mantel.csv and community_environ_mantel.csv datasets

# Then run a partial Mantel test with the same datasets, while controlling for geographic distance
# Use the community_geography_mantel.csv file as the third matrix

# Email the results to me
##############################################

#------------------------------------
# Moran's I
#------------------------------------

# Examine whether there is spatial autocorrelation in the annual rainfall data
# associated with the communities you worked with above

library(ape)

# Import the location data
community_geography_mantel <- read.csv("community_geography_mantel.csv")

# Remove the site label column
community_geography_mantel2 = community_geography_mantel[,-1]

# Convert the data to a Euclidian distance matrix;
# note that ideally you would compute spherical distances but planar distance are good enough for this exercise
# Also use the as.matrix function because you need it to be a full matrix with diagonals
# *Using the dist function alone only produces a lower half matrix with no diagonal
community_geography_mantel_dist = as.matrix(dist(community_geography_mantel2))

# Take inverse of the matrix values and replace the diagonal entries with zero

community_geography_mantel_dist_inv <- 1/community_geography_mantel_dist

diag(community_geography_mantel_dist_inv) <- 0

# Import the site data

community_rain <- read.csv("community_rain.csv")

rain_MoranI = Moran.I(community_rain$Ann_rain, community_geography_mantel_dist_inv)

rain_MoranI

# observed	the computed Moran's I.
# expected	the expected value of I under the null hypothesis.
# sd	the standard deviation of I under the null hypothesis.
# p.value	the P-value of the test of the null hypothesis against the alternative hypothesis specified in alternative.


#################################################################################
# YOU CAN USE R AS A GIS PLATFORM
#################################################################################
# Here are a couple of exercises (from Natalie Cooper) to demontrate some very basic functions in R
#################################################################################

#Natalie Cooper July 2013

#This code comes from Cooper & Nunn 2013. Identifying future zoonotic disease threats: Where are the gaps in our understanding of primate infectious diseases?
#Evolution, Medicine and Public Health. 2013. 26-37. DOI: 10.1093/emph/eot001.

#The purpose was to use a database containing records of primates with a disease
#and then to map the amount of disease sampling in different countries across the world.

#-------------------------------------------------------------------------
#Load required packages
#-------------------------------------------------------------------------

library(maptools)
gpclibPermit()
library(rgdal)
library(sp)
library(raster)

library(PBSmapping)
library(rworldmap)
library(RColorBrewer)
library(rgeos)

#-------------------------------------------------------------------------
#Read in parasite/primate host data
#-------------------------------------------------------------------------

ds<-read.delim("HostParasiteTestData.txt", header = TRUE)

#Note that the data have been pre-cleaned to simplify the code
#I have corrected the taxonomies
#removed non-matching data
#and removed non accurate localities
#The full dataset is available online if you want it

#-------------------------------------------------------------------------
#convert to R Data Structure 'SpatialPointsDataFrame'
#-------------------------------------------------------------------------

dsSPDF<-SpatialPointsDataFrame(ds[,4:5],data.frame(ds[,1:10]))
#creates dataframe - first thing needs to be long then lat, then the second is all the data
#order is CRUCIAL: LONGITUDE then LATITUDE

proj4string(dsSPDF)<-CRS("+proj=longlat")
#make sure the projection is lat long
#to match the world map

#-------------------------------------------------------------------------
#Create ID row 
#-------------------------------------------------------------------------

dsSPDF@data$ID<-1:(length(dsSPDF@data$ParasiteCorrectedName))

#This isn't needed for all analyses
#Here I need an ID for each host-parasite combination per locality per paper
#This is my measure of "sampling effort"
#Could do unique number of primate species or similar

#-------------------------------------------------------------------------
#Add world map data
#-------------------------------------------------------------------------

data(wrld_simpl)

#-------------------------------------------------------------------------
#Plot points on world map
#-------------------------------------------------------------------------

plot(wrld_simpl)
points(dsSPDF, col = "red", cex = 0.5, pch = 16)

#You can modify the plot window using xlim and ylim
#x axis runs from -179 to 179
# y axis runs from -80 to 80

#-------------------------------------------------------------------------
#Transform wrld_simpl and dsSPDF so they have the same projection
#-------------------------------------------------------------------------

dsSPDF<-spTransform(dsSPDF, CRS("+proj=longlat +ellps=WGS84"))
wrld_simpl<-spTransform(wrld_simpl, CRS("+proj=longlat +ellps=WGS84"))

#-------------------------------------------------------------------------
#Work out per country sampling effort
#-------------------------------------------------------------------------

countrysamples<-over(dsSPDF, wrld_simpl)
#This looks in each country's polygon
#and extracts how many sampling points there are in the country

dsSPDF@data$ISO3<-countrysamples$ISO3
#adds a column to dsSPDF containing country polygons
#So shows which country the samples are present in
#ISO3 is just a code in wrld_simpl for each country

#-------------------------------------------------------------------------
#Calculate number of host-parasite sampling occasions 
#per country polygon
#-------------------------------------------------------------------------

samples<-with(dsSPDF@data, aggregate(dsSPDF@data$ID, by = list(dsSPDF@data$ISO3), FUN = length))
#This aggregates the number of sampling occassions in each country polygon

#aggregate creates a dataframe with names x and Group.1. 
#Group.1 is the factor the data were aggregated by (here = country)
#x is the number of samples in the country
#Rename these to make it easier to follow your code:

samples$ISO3<-samples$Group.1
samples$sampocc<-samples$x

#remove x and Group.1 so you don't have any confusing variables hanging around!
samples<-samples[,-c(1:2)] 

#-------------------------------------------------------------------------
#Link data to map
#-------------------------------------------------------------------------

spdata <- joinCountryData2Map(samples
              , joinCode = "ISO3"
              , nameJoinColumn = "ISO3"
              , projection = "none")

#This will give the following warning:
#57 codes from your data successfully matched countries in the map
#0 codes from your data failed to match with a country code in the map
#189 codes from the map weren't represented in your data

#This just means that there are 57 countries with primate disease sampling
#No primate disease samples that aren't within a country (which is good!)
#189 countries without primate disease sampling
#Check this is what you expect for your data!

#-------------------------------------------------------------------------
#Plot map of sampling effort for each country
#-------------------------------------------------------------------------

mapCountryData(spdata, colourPalette = c("red", "orange", "yellow"), nameColumnToPlot="sampocc", borderCol = "black", mapTitle = "Sampling effort per country")

#This has the warning message:
#In rwmGetColours(colourPalette, numColours) :
  #3 colours specified and 7 required, using interpolation to calculate colours

#This just means that to evenly plot the colours on the map it needs 7 different colours
#I only gave it 3 so it interpolates the colours in between
#i.e. various shades of orange

#You can vary the colours, use RColorBrewer etc

mapCountryData(spdata, colourPalette = rainbow(7), nameColumnToPlot="sampocc", borderCol = "black", mapTitle = "Sampling effort per country")
mapCountryData(spdata, colourPalette = heat.colors(7), nameColumnToPlot="sampocc", borderCol = "black", mapTitle = "Sampling effort per country")
mapCountryData(spdata, colourPalette = terrain.colors(7), nameColumnToPlot="sampocc", borderCol = "black", mapTitle = "Sampling effort per country")

###################################################################################################

#Natalie Cooper July 2013

#This code shows you quickly how to manipulate IUCN Mammal GR maps and also WORLDCLIM env data.
#Make sure you use the proper citations for the maps and WORLDCLIM if you use this

#-------------------------------------------------------------------------
#Load required packages
#-------------------------------------------------------------------------

library(maptools)
gpclibPermit()
library(rgdal)
library(sp)
library(raster)

library(PBSmapping)
library(rworldmap)
library(RColorBrewer)
library(rgeos)

#-------------------------------------------------------------------------
#Read in IUCN range maps
#These are spatial polygons
#Available in a zipped file from the IUCN website
#Note you need the polygons so you must download terrestrial mammals only
#the all mammals file is in a different format which R can't read
#Note that shapefiles like this come in multiple files
#so you will have files called hc_grid1d.dbf, hc_grid1d.prj, hc_grid1d.shp and hc_grid1d.shx
#You need ALL of them in your working directory to make this work
#-------------------------------------------------------------------------

maps<-readShapeSpatial("MAMMTERR")#This is SLOW; it took about 3 mins to run on my fast computer

proj4string(maps) <-CRS("+proj=longlat +datum=WGS84")
#ensure the projection is correct***This is really important! check your map projection before use!

#If you need to correct the taxonomy so it matches the data:
maps@data$BINOMIAL<-gsub("Mico argentatus","Callithrix argentata", maps@data$BINOMIAL)

#You can plot range maps either individually or on the world map:

plot(maps[which(maps@data$BINOMIAL == "Pan troglodytes"),])

#or:

data(wrld_simpl)
plot(wrld_simpl)
plot(maps[which(maps@data$BINOMIAL == "Pan troglodytes"),], add = TRUE, col = "green")

#-------------------------------------------------------------------------
#To look at overlap/intersections between two species:
#-------------------------------------------------------------------------
#****NOTE: Creating the "intersects" object didn't work for me but go ahead with the other
#****lines of code and see what you get***********************************

intersects<-gIntersection(maps[which(maps@data$BINOMIAL == "Pan troglodytes"),], 
		maps[which(maps@data$BINOMIAL == "Colobus guereza"),])

plot(wrld_simpl)
plot(maps[which(maps@data$BINOMIAL == "Pan troglodytes"),], add = TRUE, col = "green")
plot(maps[which(maps@data$BINOMIAL == "Colobus guereza"),], add = TRUE, col = "blue")
plot(intersects, add= TRUE, col = "yellow")

#-------------------------------------------------------------------------
#Manipulating IUCN range maps
#This is easy, you can subset as you would a normal dataframe
#You just need to use maps@data$VariableName, rather than the standard maps$Variable Name
#because now your dataframe is just part of the whole maps object
#-------------------------------------------------------------------------
#For example, to replace spaces in species names with _ 

maps@data$BINOMIAL<-gsub(" ", "_", maps@data$BINOMIAL)

#-----------------------------------------------
#Plotting colourful maps of environmental data
#----------------------------------------------

#Read in worldclim data (available directly through R)
#extract WORLDCLIM data. Object is a RasterStack with a layer for each variable
#Note you can vary the resolution of the data you extract

bioclim<-getData("worldclim", var = "bio", res = 10)#takes a while the first time

bioclim<-unstack(bioclim)#unstacks raster so you can deal with each variable individually

plot(bioclim[[6]], zlim=c(-100,250), axes = F, col=rev(heat.colors(25)), ext = matrix(c(150,50,-120,-40), nrow = 2))
#zlim defines the range of the variable (temperature in this case = bioclim 6)
#extent defines the lat and long you want to display
#you can use any of R's colour palettes eg rainbow, terrain.colors or a user defined set

plot(bioclim[[6]], zlim=c(-100,250), axes = F, col=rev(rainbow(25)), ext = matrix(c(150,50,-120,-40), nrow = 2))
plot(bioclim[[6]], zlim=c(-100,250), axes = F, col=rev(heat.colors(25)), ext = matrix(c(150,50,-120,-40), nrow = 2))

#Note that "rev" reverses the colour scale, because I prefer red = hot. 

#You can add the parasite location points into this etc.

###################################################################################################

## You can also do species distribution modeling in R using the dismo package; the associated vignette is helpful
## http://cran.r-project.org/web/packages/dismo/

# A list of R packages that work with spatial data is here: http://cran.r-project.org/web/views/Spatial.html

# A couple of other good website for working with spatial data in R are:

# http://www.people.fas.harvard.edu/~zhukov/spatial.html
# http://spatial.ly/r/

###################################################################################################







