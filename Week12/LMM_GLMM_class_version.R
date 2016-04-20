# Linear mixed models (LMM) and generalized linear mixed models (GLMM)

library(lme4)

library(lmerTest)

library(MuMIn)

#-------------------------------------------------------------
# LMM PREDICTING GROUP SPREAD USING FOCAL INDIVIDUAL AND DATE AS RANDOM EFFECTS 

# Mixed models are similar to regular linear models except that they include at least one random effect. 
# These random effects essentially give structure to the error term of the mode.
# In the case of our models below, we add a random effect for each focal individual and each date. 
# These random effect quantifies idiosyncratic variation that is due to different individuals and different dates.
# For more info on mixed models a brief and good tutorial can be found here
# http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf
#-------------------------------------------------------------

new_data_climate_gpspread_nomiss <- read.csv("new_data_climate_gpspread_nomiss.csv")

View(new_data_climate_gpspread_nomiss)

# Create the LMM using the lmer function
# Fixed effects are included in the model in the same way as a typical linear model
# Random effects are included in the model using the (1 | variable) term
# The "1" indicates that we are modeling individual differences by assuming different random intercepts for each subject.
# That is, each subject (in this case the Focal Individual and the Date) is assigned a different intercept value, and the mixed model estimates these intercepts for you.

gp_spread_model001=lmer(GP_SPRD ~ Comp.1 + Comp.2 + REPRO + ACTIVITY + GP_TYPE + (1 | FOCAL) + (1 | DATE), data=new_data_climate_gpspread_nomiss)

# Produce the model output
summary(gp_spread_model001)

# In the Random Effects table of the output
# Examine the Std Dev column. 
# This shows the amount of variation in Y explained by Date, Focal, and Residual variation accounting for the random effects

# You will notice that p values are generated for each predictor but not for the overall model.
# You can examine the statistical significance of the overall model by using a log likelihood ratio test.
# Create a null model only using the random effects.
# Compare each model by using the anova function, for instance
# anova(gp_spread_model001_null, gp_spread_model001)
# If you get a signigncant p value, then the model with predictors is significantly better than the null model

# Now use AICc to determine the best combination of predictors (fixed effects)

# Change the na.action setting in preparation for the multi-model selection analysis
options(na.action=na.fail)

# Run the AICc analysis
gp_spread_model001_Aicc=dredge(gp_spread_model001)

# Produce the results
gp_spread_model001_Aicc

# Which predictor is most important?
importance(gp_spread_model001_Aicc)

##### Re-run the model above but only include FOCAL as a random effect
##### Re-run the model above but only include DATE as a random effect
##### Re-run the model above but include FOCAL and DATE as interacting random effects

##### Based on your analyses, how should you model the random effects? 
##### Base your answer on the AICc values of the full model from each analysis

# Is colinearity among the predictors a problem?
### Checking VIF using mer-utils.R code from the weblog of the Human Language Processing (HLP) lab at the University of Rochester.  https://hlplab.wordpress.com/2011/02/24/diagnosing-collinearity-in-lme4/

source('C:/Users/jkamilar/Desktop/R working folder/mer-utils.R')

vif.mer(MODEL_NAME)


#----------------------------------------------------------------------
#**********************************************************************
#----------------------------------------------------------------------


#-------------------------------------------------------------
# GLMM PREDICTING GROUP SIZE USING FOCAL INDIVIDUAL AND DATE AS RANDOM EFFECTS
#-------------------------------------------------------------

new_data_climate_gpsize_nomiss <- read.csv("new_data_climate_gpsize_nomiss.csv")

View(new_data_climate_gpsize_nomiss)

# Note that we are now using the glmer function (instead of the lmer function)
# Now you need to specify a distribution family. In this case we are modeling a poisson distribution.

gp_sz_model001=glmer(GP_SZ ~ Comp.1 + Comp.2 + REPRO + ACTIVITY + GP_TYPE + (1 | FOCAL) + (1 | DATE), family="poisson", data=new_data_climate_gpsize_nomiss)

summary(gp_sz_model001)

options(na.action=na.fail)

gp_sz_model001_AICc=dredge(gp_sz_model001)

gp_sz_model001_AICc

importance(gp_sz_model001_AICc)

# Create a null model only using the random effects.
# Compare each model by using the anova function

##### Re-run the model above but only include FOCAL as a random effect
##### Re-run the model above but only include DATE as a random effect
##### Re-run the model above but include FOCAL and DATE as interacting random effects

##### Based on your analyses, how should you model the random effects? 
##### Base your answer on the AICc values of the full model from each analysis

gp_sz_model002=glmer(GP_SZ ~ Comp.1 + Comp.2 + REPRO + ACTIVITY + GP_TYPE + (1 | FOCAL), family="poisson", data=new_data_climate_gpsize_nomiss)

summary(gp_sz_model002)

options(na.action=na.fail)

gp_sz_model002_AICc=dredge(gp_sz_model002)

gp_sz_model002_AICc

importance(gp_sz_model002_AICc)

#Importance of ACTIVITY switches with Comp.1
#> importance(gp_sz_model001_AICc)
#GP_TYPE REPRO Comp.2 Comp.1 ACTIVITY
#Importance:          1.00    0.97  0.82   0.66   0.04    
#N containing models:   16      16    16     16     16    
#> importance(gp_sz_model002_AICc)
#GP_TYPE REPRO Comp.2 ACTIVITY Comp.1
#Importance:          1.00    1.00  1.00   1.00     0.82  
#N containing models:   16      16    16     16       16  
