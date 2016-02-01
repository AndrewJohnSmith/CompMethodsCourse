# Bivariate regression and allometry
# open brain_bodymass.csv dataset
# this dataset contains three variables, clade, female endocranial volume, and female body mass

#Set working directory
setwd("/Users/Home/Dropbox/Amherst/Courses/CompMethods/CompMethodsCourse/Week2")

brain_bodymass<-read.table("brain_bodymass.txt", header=T, row.names=1)

#NOT run
#DataCSV<-read.csv("brain_bodymass.txt", sep="\t")
  
# Explore the data
str(brain_bodymass)

summary(brain_bodymass)

cor(brain_bodymass[,3:4])

plot(brain_bodymass[,3:4], col=unclass(brain_bodymass$clade), pch=20)

# Examine the relationship between ECV and body mass using OLS model

model_01<-lm(ecv_female~mass_female, data=brain_bodymass)
model_02<-lm(ecv_female~mass_female+I(mass_female^2), data=brain_bodymass)


summary(model_01) # shows an output table of your model

# Plot the regression

plot(ecv_female~mass_female, data=brain_bodymass, main="Body Mass vs. ECV - Raw Values", pch=20, col=unclass(brain_bodymass$clade))

# Draw a best fit line

abline(model_01, col="red", lwd=2)
abline(model_02, col="blue", lwd=2)

# Examine some model diagnostics
# First, split the graphics windown into four sections

par(mfrow=c(2,2))

plot(model_01)   # Data that may be outliers are labeled by their case number

# Can also create a Cook's distance plot more directly
# First, reset the graphics window

par(mfrow=c(1,1))

plot(cooks.distance(model_01))
  #A cooks of >1 is normally not good.

################################################

# Run the same variables but log10 transform them this time

model_02<-lm(log10(ecv_female)~log10(mass_female), data=brain_bodymass)

summary(model_02) # shows an output table of your model


# Take a look at the model diagnostics
par(mfrow=c(2,2))
plot(model_02)

# How do these diagnostics compare to the model using untransformed data?
  #Looks much better Q-Q is looking like a positive relationship.
  # Other plot lines are flat (which is good).

# Calculate a 95% confidence interval
confint(model_02, 'log10(mass_female)',level=0.95)

# You can also use the predict () function to get the fitted values and confidence interval for all samples
predict1<-predict(model_02, interval="confidence")

# You can view the fitted values by typing predict1. Though it will be a long list.
# You can also save the values to a file and open it in Excel
write.csv(predict1, "fitted_values.csv")

par(mfrow=c(1,1))
# Plot the regression
plot(log10(ecv_female)~log10(mass_female), data=brain_bodymass, main="Body Mass vs. ECV - Log10 Values", pch=20, col=unclass(brain_bodymass$clade))

# Draw a best fit line

abline(model_02, col="blue", lwd=2)

# Save your plot by clicking the Export button in RStudio

# You can generate the model residuals. Save them as an object first.
my_resid<-resid(model_02)

# You can save them to a csv file
write.csv(my_resid, "raw_resid.csv", row.names=F)

# Open the csv file and view the residuals

# Though, standardized residuals are more useful
my_st_resid<-rstandard(model_02)

#View what standardized vs unstandized residials are.
  boxplot(my_resid)
  boxplot(my_st_resid)

# You can save them to a csv file
write.csv(my_st_resid, "st_resid.csv", row.names=F)

# Open the csv file and view the standardized residuals

###########################################

# Rerun the last regression but switch the independent and dependent variables
# How do the results differ?
# Plot X vs. Y. How does the plot differ?
# Email me both of the plots (one with body mass as the predictor and the other with brain size as predictor)

par(mfrow=c(1,2))
plot(log10(ecv_female)~log10(mass_female), data=brain_bodymass, main="Original Body Mass vs. ECV - Log10 Values", pch=20, col=unclass(brain_bodymass$clade))
plot(log10(mass_female)~log10(ecv_female), data=brain_bodymass, main="Switched Body Mass vs. ECV - Log10 Values", pch=20, col=unclass(brain_bodymass$clade))
par(mfrow=c(1,1))
###########################################

# Standard (reduced) major axis regression
# Install and load lmodel2 package

install.packages("lmodel2")
require(lmodel2)
?lmodel2
rma_model01<-lmodel2(log10(ecv_female)~log10(mass_female), data=brain_bodymass, nperm=999)

rma_model01

plot(rma_model01, "SMA")      # Plots best fit line and 95% conf limits

###########################################
# Run an ANOVA using the anova() function

anova_model01<-aov(log10(mass_female)~clade, data=brain_bodymass)

anova_model01

summary(anova_model01)
#View multiple comparisons
TukeyHSD(anova_model01)

par(mfrow=c(2,2))

plot(anova_model01)      # view diagnostic plots

par(mfrow=c(1,1))

boxplot(log10(mass_female)~clade, data= brain_bodymass, ylab="log female mass")


##########################################
# You can also use the linear model function to conduct an ANOVA, i.e. conduct a linear model with a categorical predictor

anova_model02<-lm(log10(mass_female)~clade, data=brain_bodymass)

summary(anova_model02) # Notice that R lists each group within clade separately as a predictor (in alphabetical order) and sets the first group as the intercept.
# The mean value of the first group (Catarrhines) is the estimate of the "intercept". The mean values of the other groups can be calculated by adding or substracting
# their estimate values (depending on the sign) to the estimate for the "intercept"
# How do the results compare running a model using the lm() function vs. the anova() function?

par(mfrow=c(2,2))

plot(anova_model02)     # view diagnostic plots

# the package ggplot2 has lots of nice features for graphics
# http://docs.ggplot2.org/current/


# Save your workspace by clicking the disk icon in the upper right window of RStudio

# DONE FOR TODAY

# Interacting with a scatterplot
attach(brain_bodymass)
plot(log10(ecv_female)~log10(mass_female), data=brain_bodymass, main="Original Body Mass vs. ECV - Log10 Values", pch=20, col=unclass(brain_bodymass$clade)) # scatterplot
identify(log10(ecv_female)~log10(mass_female), labels=row.names(brain_bodymass)) # identify points
coords <- locator(type="l") # add lines
coords # display list 
detach(brain_bodymass)
