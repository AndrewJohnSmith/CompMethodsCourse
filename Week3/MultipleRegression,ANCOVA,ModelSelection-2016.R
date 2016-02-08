# Multiple regression, ANCOVA, model selection

# Multiple Regression##############################################

# Import the mult_regr_data.csv data file
setwd("/Users/Home/Dropbox/Amherst/Courses/CompMethods/CompMethodsCourse/Week3")
require(car)
mult_regr_data<-read.csv("mult_regr_data.csv", row.names=1)

# Explore the data

str(mult_regr_data)

summary(mult_regr_data)

# Create a correlation matrix, similar to the one produced in Week 1
cor(mult_regr_data)
pairs(mult_regr_data)

# Create a multiple regression using two predictor variables
mult_reg_model01<-lm(MaxLongevity_m~Precip_Mean_mm+AdultBodyMass_g, data=mult_regr_data)

summary(mult_reg_model01)

# Check the model assumptions ; first change the graphics window to display 4 graphs (2 X 2)

par(mfrow=c(2,2))

plot(mult_reg_model01)

# The model would probably fit the assumptions better if the variables were log transformed

# Re-run the model using log transformed data

mult_reg_model02<-lm(log10(MaxLongevity_m)~log10(Precip_Mean_mm)+log10(AdultBodyMass_g), data=mult_regr_data)

summary(mult_reg_model02)

plot(mult_reg_model02)

# Check the variance inflation factor values
# first install the car package
require(car)

vif(mult_reg_model02)
#Close to one so variables are not colinear

# ANCOVA ####################################################
# Use the brain_bodymass dataset for this section

# ANCOVA is a special type of multiple regression (i.e. linear model) with at least one categorical predictor and a second predictor (continuous or categorical)
# Tests for differences in the intercept (mean value) of groups. This is often used to test for grade shifts among clades.

# Import the brain_bodymass.csv dataset
brain_bodymass<-read.csv("brain_bodymass.csv", row.names=1)

# Test for differences in the Y value among clades in the relationship between body mass and brain size.
# In other words, for a given body mass, do different clades have different brain sizes?

# Use an interaction term to test for homogeneity of slopes. An interaction term is designated by an asterisk between variables.

ancova_model01<-lm(log10(ecv_female)~log10(mass_female)*clade, data = brain_bodymass)

# Interaction plots can be made in R using the command:
# interaction.plot(factorA, factorB, Response)
# first reset your graphics window

par(mfrow=c(1,1))

interaction.plot(log10(brain_bodymass$mass_female), brain_bodymass$clade, log10(brain_bodymass$ecv_female))

# The plot is a bit funky. There are probably some options to make it nice here
# https://stat.ethz.ch/R-manual/R-patched/library/stats/html/interaction.plot.html

# Switch the order of factor A and B to alter which variable is plotted on the x-axis, 
# though in this case switching the variables in the plot would make it worse

# If the interaction term is significant, then the slopes differ among clades.
# This would violates the homogeneity of slopes assumption and differences in the Y values cannot be examined.
# Since we have three groups in the clade variable, we have two interaction terms

summary(ancova_model01)

# How would you interpret the results?
  #Significant difference in slopes b/w mass and Platyrrhine to Caterrine
# If the interaction term is not signficant, then re-run the analysis without the interaction term.
# In this example, assume no interaction term is significant

ancova_model02<-lm(log10(ecv_female)~log10(mass_female)+clade, data = brain_bodymass)

summary(ancova_model02)

# How would you interpret the results?
  #Sig dif in mass between catherines and platyrrhines
# Check the model assumptions with plot (); though first modify your graphics window
par(mfrow=c(2,2))

plot(ancova_model02)

#########################################################################

# Model selection #

# Use the mult_regr_data dataset for this section

# If you have several independent variables, then you may not want to include all of them in the model,
# especially if some of the variables are poor predictors - this can lead to overfitting the model.


# stepwise model selection by AIC  

# Install the MASS package
require(MASS)
fit<-lm(MaxLongevity_m~SocialGrpSize+Precip_Mean_mm+AdultBodyMass_g+SexualMaturityAge_d, data=mult_regr_data)

step <- stepAIC(fit, direction="both")

step$anova # display results

# Lower AIC models indicate a better model
  #Precipitation predictor is not a decent addition to the model

#########################################################################

# This approach creates all possible models using your predictor variables and
# uses AICc to judge the best models and best predictor variables 

# Install MuMIn package
require(MuMIn)
model01x<-lm(log10(MaxLongevity_m)~log10(SocialGrpSize)+log10(Precip_Mean_mm)+log10(AdultBodyMass_g)+log10(SexualMaturityAge_d), data=mult_regr_data)

options(na.action=na.fail) # Only performs analysis if no missing values are in dataset

dd<-dredge(model01x) # Produces a set of models based on all possible combinations of predictor variables

dd    # View a summary of the models

dd_models<-get.models(dd, subset="TRUE")  # Produces a list of which predictors are in each model

dd_models  # This will list the specific variables in each model  

# Models within 2 AICc values of the best model are considered equally good
# The "delta" column lists the difference in AICc value between that model and the best model

importance(dd)  # Sum of Akaike weights over all models that include the explanatory variable
# Variables with an Akaike weight closer to 1 are more important than those with values close to 0

# Which model(s) is the best? Which variable(s) is most important? Email me your answers.
  #Model 5 was the best, adult body mass was the 'most important' variable
