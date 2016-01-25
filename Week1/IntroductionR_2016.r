#############################################
### Modified from several sources including:
### VertNet Biodiversity Informatics Training Workshop
### https://github.com/mkoo/BITW
### and
### Rogery Mundry's Quick Intro to R lecture at Anthrotree 2014
#
### Introduction to R
#
### for R.3.0.1
#
#################################################

##R is a calculator (try your own):
555/66
2+3
#Respects laws of precedence
2+3*4
(2+3)*4
###
#Creating objects: operator is "<-" OR "="
#R is case sensitive so q and  Q are different variables.
x <- 2

y = 3

z <- 5

x+y*z

result <- x+y*z

result

#Objects are stored in the workspace.
ls()         #see them all
rm(x)      #remove 'x'
ls()
#also this
objects()

# In RStudio a list of the objects can be found in the upper right frame


#############
###       Vectors
xvec <- 1:5
xvec

yvec <-10:15
#see yvec

zvec <- c(10,17,12,8,5,1)
zvec

#    [] square brackets allows returning values by indexing
zvec[3]
yvec[2:4]
#     What happens when the index is out of range?
xvec[8]

#To convert one atomic value to another these are handy:
#as.character(x)
#as.integer(x)
#as.logical(x)
#as.numeric(x)
#as.factor(x)

#########################
###      Working directories and Datasets - Import/Export
###      Installing packages, getting down to business!
#
#setwd(), dir(), getwd(), list.files(), install.packages()
#Majority of the time you will likely be importing data already in a table
# In Windows, directory paths / need to be escaped with double//. Can copy/paste from
# Explorer address bar, then switch the slashes

getwd() # Where is the working directory?
setwd("~/Desktop/BITW-Workshop")  #set for a Mac OS

# Easier to set working directory if you use RStudio
# Point & click in RStudio: Tools --> Global Options
# Copy the LH_data1.csv file into your working directory

dir() #See what is in directory
list.files() #does same as above but may take longer

# Install a package using RStudio
# Click the Packages button in the lower right frame
# Click Install Pacakages
# Search for your package of interest
# Click install; be sure that install dependencies is checked
# Installing a package put the package on your computer but you cannot use it in R until you load it
# Load a package using RStudio by clicking next to the package name in your package list window (lower right frame)
# Now you can use the package


# Import a csv file
Data <- read.csv("LH_data1.csv", header=T)

# You can import  other file types. For instance, if you had a tab delimted file you could use
# data=read.table(file="file_1.txt",header=T,sep="\t")

# Import data using RStudio
# Click Import dataset button in upper right frame
# Notice that this will name the object LH_data1

# See a sample of the data you just read into the working space:
colnames(LH_data1)  #see column names in the header row
head(LH_data1)       # First 6 rows of data frame


# Other handy ways to quickly peek at your data:
tail(LH_data1)      # Last 6 rows-- I like to check that no errant comma(s) or tab(s) has messed up the fields

# You can see all your data in RStudio in the upper left frame!

# To check the structure of your data, use str()
str(LH_data1)

# Some more info about the data can be obtained using summary()
summary(LH_data1)

# To create a correlation matrix of all your numerical data use cor()
cor(LH_data1[,2:9]) # The square brackets are needed because the first variable is not numerical (species) and would produce an error. The numbers after the comma refer to columns 2 through 9.

# To create a plot of the correlations use pairs()
pairs(cor(LH_data1[,2:9]))

# You can find some basic functions (calculate mean, median, max, min, etc.) in R here http://www.statmethods.net/management/functions.html

# Delete a row from the dataset
data_new<-LH_data1[-3,]

# Save the new dataset as a csv file
write.csv(data_new, "data_new.csv", row.names=FALSE)



#########################
# R is a great graphing environment!

# Exercise for making a histogram and modifying a plot; from http://www.r-bloggers.com/basics-of-histograms/

# Create some fake normally distributed data
BMI<-rnorm(n=1000, m=24.2, sd=2.2) 

# Make a histogram of the data; RStudio creates plots in the lower right frame
hist(BMI)

# To understand how R created the histogram:
#Save the histogram as an object
histinfo<-hist(BMI)

# Produce the output of info used by R
histinfo

# Change the number of bins
hist(BMI, breaks=20, main="Breaks=20")
hist(BMI, breaks=5, main="Breaks=5")

# Set the exact  start and end point of each bin
hist(BMI, breaks=seq(17,32,by=3), main="Breaks is vector of breakpoints")

# Instead of counting the number of datapoints per bin, R can give the probability densities using the freq=FALSE option:
hist(BMI, freq=FALSE, main="Density plot")
lines(density(BMI), col="red", lwd=2)

#Plot just the densities
plot(density(BMI), main="Density plot")

# This code does several things: adjusts the x axis label; creates a title for the graph; creates histogram bars that are light green; modifies the scales of the x and y axes
# These functions can be used with other types of plots
hist(BMI, freq=FALSE, xlab="Body Mass Index", main="Distribution of Body Mass Index", col="lightgreen", xlim=c(15,35),  ylim=c(0, .20))

# If you want to get a list of R color options, then type
colors()

# Create a scatterplot using your imported dataset; plot body mass vs. longveity
# pch sets the symbol for the samples; see this site for more symbol options (as well as other graphical parameters) http://www.statmethods.net/advgraphs/parameters.html
# To tell R which variables you want to plot, you have to specify the dataset and then the variable in the dataset. This is done with the dollar sign.
plot(LH_data1$AdultBodyMass_g, LH_data1$MaxLongevity_m, main="My First Scatterplot", xlab="Body Mass ", ylab="Max. Longevity", pch=19)
#Add colours
plot(LH_data1$AdultBodyMass_g, LH_data1$MaxLongevity_m, main="My First Scatterplot", xlab="Body Mass ", ylab="Max. Longevity", pch=19, col=unclass(Data$Clade))

# Log10 transforming the data would make the plot nicer. If you use log (as opposed to 1og10), then R will use the natural log.
plot(log10(LH_data1$AdultBodyMass_g), log10(LH_data1$MaxLongevity_m), main="My First Scatterplot", xlab="Body Mass ", ylab="Max. Longevity", pch=19)
#Adding colours
plot(log10(LH_data1$AdultBodyMass_g), log10(LH_data1$MaxLongevity_m), main="My First Scatterplot", xlab="Body Mass ", ylab="Max. Longevity", pch=19, col=unclass(Data$Clade))

# In RStudio, click the Zoom button in the lower right frame to get a better view of your plot
# You can save your plot by clicking the Export button
# The file is saved in your R working directory

# You can also partition the graphics window so more than one plot can be display simultaneously
par(mfrow=c(2,1))

# Recreate the unlogged and logged plots
plot(LH_data1$AdultBodyMass_g, LH_data1$MaxLongevity_m, main="My First Scatterplot", xlab="Body Mass ", ylab="Max. Longevity", pch=19)
plot(log10(LH_data1$AdultBodyMass_g), log10(LH_data1$MaxLongevity_m), main="My First Scatterplot", xlab="Body Mass ", ylab="Max. Longevity", pch=19)

# reset the graphics frame to one plot
par(mfrow=c(1,1))

# Create a boxplot of social group size for each clade in the dataset
boxplot(SocialGrpSize~Clade,data=LH_data1, main="Differences in Social Group Size", xlab="Clade", ylab="Social Group Size")

# The boxplot shows the following aspects of the data for each group

# the lower quartile, defining the lower limit of the box in each figure;
# the sample median, represented by the heavy line inside each box;
# the upper quartile, defining the upper limit of each box;
# the upper end of the nominal data range is defined as the upper quartile plus 1.5 times the IQD  (IQD = interquartile distance, i.e. the distance from the upper to lower quartile)
# the lower end of the nominal data range is defined as the lower quartile minus 1.5 times the IQD
# values outside the nominal data range are defined as "outliers" and noted with an open circle

# Log10 transform social group size to get the plot in log space
boxplot(log10(SocialGrpSize)~Clade,data=LH_data1, main="Differences in Social Group Size", xlab="Clade", ylab="Social Group Size")

#Do a basic statistical test
t.test(log10(Data$SocialGrpSize)~Data$Clade)

# Save your R workspace (i.e. all the objects and data from this session)
# In RStudio, click the disk icon in the upper right frame
# You can open your saved workspace later by clicking the open folder icon next to the disk icon

#####DONE FOR TODAY