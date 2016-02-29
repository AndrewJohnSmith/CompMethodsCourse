setwd("~/Dropbox/Amherst/Courses/CompMethods/CompMethodsCourse/Week6")

# Modeling Trait Evolution

library(ape)

library (phytools)

primate_traits <- read.csv("primate_traits.csv")
primate_tree05 <- read.tree("primate_tree05.txt")

female_mass <- primate_traits$female_mass
names(female_mass) <- primate_traits$species


#----------------------------------------------------
# Model trait evolution via Brownian motion
#----------------------------------------------------

library(geiger)

modelBM<-fitContinuous(primate_tree05, female_mass, model = "BM")

modelBM

#GEIGER-fitted comparative model of continuous data
# fitted ?BM? model parameters:
#	sigsq = 8062695.905460
#	z0 = 8922.262375

# model summary:
#	log-likelihood = -336.438483
#	AIC = 676.876966
#	AICc = 677.305538
#	free parameters = 2

#Convergence diagnostics:
#	optimization iterations = 100
#	failed iterations = 0
#	frequency of best fit = 1.00

# object summary:
#	'lik' -- likelihood function
#	'bnd' -- bounds for likelihood search
#	'res' -- optimization iteration summary
#	'opt' -- maximum likelihood parameter estimates

modelBM$opt

#$sigsq
#[1] 8062696
#$z0
#[1] 8922.262
#$lnL
#[1] -336.4385
#$method
#[1] "Brent"
#$k
#[1] 2
#$aic
#[1] 676.877
#$aicc
#[1] 677.3055

#---------------------------------------------------------------
# Model trait evolution via Ornstein-Uhlenbeck model
#---------------------------------------------------------------

modelOU<-fitContinuous(primate_tree05, female_mass, model = "OU")

modelOU

#GEIGER-fitted comparative model of continuous data
# fitted ?OU? model parameters:
#	alpha = 0.017382
#	sigsq = 11200378.611272
#	z0 = 9224.008374

# model summary:
	log-likelihood = -335.942007
	AIC = 677.884013
	AICc = 678.772902
	free parameters = 3

Convergence diagnostics:
	optimization iterations = 100
	failed iterations = 0
	frequency of best fit = 0.04

 object summary:
	'lik' -- likelihood function
	'bnd' -- bounds for likelihood search
	'res' -- optimization iteration summary
	'opt' -- maximum likelihood parameter estimates
> summary(modelOU)
    Length Class      Mode
lik   1    bm         function
bnd   2    data.frame list
res 400    -none-     numeric
opt   8    -none-     list

modelOU$opt

$alpha
[1] 0.01738177

$sigsq
[1] 11200379

$z0
[1] 9224.008

$lnL
[1] -335.942

$method
[1] "subplex"

$k
[1] 3

$aic
[1] 677.884

$aicc
[1] 678.7729


#------------------------------------
# Model trait evolution via an early burst model
#------------------------------------

modelEB<-fitContinuous(primate_tree05, female_mass, model = "EB")

modelEB

GEIGER-fitted comparative model of continuous data
 fitted ?EB? model parameters:
	a = -0.000001
	sigsq = 8063208.908645
	z0 = 8922.254795

 model summary:
	log-likelihood = -336.438512
	AIC = 678.877024
	AICc = 679.765913
	free parameters = 3

Convergence diagnostics:
	optimization iterations = 100
	failed iterations = 0
	frequency of best fit = 0.04

 object summary:
	'lik' -- likelihood function
	'bnd' -- bounds for likelihood search
	'res' -- optimization iteration summary
	'opt' -- maximum likelihood parameter estimates


###### You can compare the models by examining their AICc values  #########

#**********************************************************
# Re-run the three analyses above using log female mass

# What is the most likely model of trait evolution?
#**********************************************************

library(qpcR)

modelLBM<-fitContinuous(primate_tree05, log(female_mass), model = "BM")
modelLOU<-fitContinuous(primate_tree05, log(female_mass), model = "OU")
modelLEB<-fitContinuous(primate_tree05, log(female_mass), model = "EB")

model.primate<-matrix(,3,4,dimnames = list(c("Brownian Motion", "Early Burst", "Ornstein-Uhlenbeck"),c("log likelihood", "AICc", "Delta AICc", "AICc Weights")))

model.primate[,1]<-c(modelLBM$opt$lnL, modelLEB$opt$lnL, modelLOU$opt$lnL)
model.primate[,2]<-c(modelLBM$opt$aicc, modelLEB$opt$aicc, modelLOU$opt$aicc)

aic.all.primate<-as.matrix(model.primate[,2])

scor.wts.primate<-akaike.weights(aic.all.primate)

model.primate[,3]<-scor.wts.primate$deltaAIC
model.primate[,4]<-scor.wts.primate$weights
model.primate

#                   log likelihood    AICc  Delta AICc   AICc Weights
#Brownian Motion         -32.04379 68.51615   0.000000    0.6310529
#Early Burst             -32.04324 70.97536   2.459217    0.1845243
#Ornstein-Uhlenbeck      -32.04379 70.97646   2.460317    0.1844228

#---------------------------------------------------------
# Rate of Discrete Trait Evolution using the diversitree package
# note: code based on the tutorial accompanying diversitree
#---------------------------------------------------------

install.packages("diversitree")

library(diversitree)

# Create a "birth-death" tree; birth rate = 0.1; death rate = 0.03
# In a "real" situation you would import your phylogeny instead

phy <- tree.bd(c(.1, .03), max.taxa=50)

plot(ladderize(phy))

# Simulate a discrete binary trait on the tree
# The trait starts in state 0; x0=0
# Simulates the evolution of a trait with an average rate of transition from state 0 to 1 is 0.1; average rate of transition from 1 to 0 is 0.2
# In a "real" situation you would import your data (binary trait) instead

states <- sim.character(phy, c(.1, .2), x0=0, model="mk2")

#----------------------------------------------------------
# We will use a maximum likelihood approach
#----------------------------------------------------------

# Create a likelihood function of trait evolution given the data and tree
# This is based on Pagel 1994 method
# Should produce similar results to the ace function in the ape package

lik.mk2 <- make.mk2(phy, states)

# Run a maximum likelihood analysis that yields the rate of trait evolution (i.e. rate of transition from 0 to 1 and 1 to 0) 

fit.mk2 <- find.mle(lik.mk2, c(.1, .1), method="subplex")  # the c(.1, .1) sets the inital parameters; setting different inital parameters should yield the same result

coef(fit.mk2)

      q01       q10 
0.1210342 0.2575020     # Our results yield values that are close to our simulation setting, which was q01=0.1 and q10=0.2

# why is there a difference between our simulation setting values and the maximum likelihood results?

# Creates a likelihood function where the transition rates (q values) are constrained to be equal

lik.mk1 <- constrain(lik.mk2, q10 ~ q01)

fit.mk1 <- find.mle(lik.mk1, .1, method="subplex")  # Only one transition rate (.1) is needed here because the q values are equal

# Can run a ANOVA function to test if two different transition rates is a better fit than one transition rate (i.e. rate from 0 to 1 and 1 to 0 is equal)
# Conducts a log likelihood ratio test

anova(fit.mk2, mk1=fit.mk1)

     Df   lnLik    AIC  ChiSq Pr(>|Chi|)  
full  2 -29.912 63.825                    
mk1   1 -32.848 67.695 5.8703     0.0154 *
---
Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1

# These results suggest that there is a significant difference between an equal rates of trait transition vs. unequal rates

# Why are the degrees of freedom different for the two models (full vs. mk1)?
  #Two transition rates vs. 1?

#**** What does it mean if the transition rate is higher in one direction than the other?? *******************
  #More likely to observed transtions in one direction (i.e., gain of function rather than loos of function).

#----------------------------------------------------
# Sample the likelihood distribution using a Bayesian approach (MCMC)
#----------------------------------------------------

prior.exp <- make.prior.exponential(10) # sets the prior distribution

# The prior distrbution defines the uncertainty before the data are taken into account

# Definition 1: prior represents a population of possible parameter values, from which the current value has been drawn
# Definition 2: prior expresses our knowledge of the parameter as if the value of the parameter is a random realization from the prior distribution

# The exponential distribution is defined from 0 to infinity, which is good because the rate of trait transition cannot be negative

# Diversitree also allows for a uniform prior; this is an non-informative prior (i.e. you don't have much prior information)

samples <- mcmc(lik.mk2, c(.1, .1), nsteps=10000, prior=prior.exp, w=.1, print.every=1000)
# lik.mk2 is the likelihood function from above
# c (.1, .1) is the starting point for the Markov chain
# nsteps is the number of steps or the number of times the parameter space is searched by the chain
# in a real analysis, a non-informative prior is often used. in this case, the number of steps is much higher, often in the millions (1-5 million)
# prior is the prior distrbution we set above
# w is the tuning parameter;  related to the "acceptance rate" for each step, which you usually want to be 20-40%
# print.every = records the parameters every X number of steps

samples <- subset(samples, i > 1000)
# From the 5000 steps, only consider steps >500
# The first 500 steps are considered the "burnin" phase
# During the burnin phase, the chain will often search areas of parameter space that are not useful
# A conservative approach is to not use the first 10% of the total steps 

mean(samples$q01 > samples$q10)

# The frequency of samples where the q01 parameter is greater than the q10 parameter
  #0.1146667 - just over 10%

par(mfrow=c(1,1))

col <- c("#004165", "#eaab00")
# sets the colors for parameters q01 and q10

profiles.plot(samples[c("q01", "q10")], col.line=col, las=1, legend.pos="topright")

# Posterior probability distributions for the parameters of the Mk2 model. 

abline(v=c(.1, .2), col=col)

# True values are indicated by the solid vertical lines.

# Save your plot


#-------------------------------------
# NOW USE "REAL" DATA
#-------------------------------------

primate_binary <- read.csv("primate_binary.csv") # this dataset quantifies the mandible of primates as robust (1) or not robust (0)

robustmandible <- primate_binary$robustmandible   # extract the variable to a new object

names(robustmandible) <- primate_binary$species  # assign species names to each row

tree_10 <- read.tree("primate_tree05.txt") # read the phylogeny

lik.mk20 <- make.mk2(tree_10, robustmandible)

fit.mk20 <- find.mle(lik.mk20, c(.1, .1), method="subplex")

#coef(fit.mk20)
#       q01        q10 
#0.03168182 0.03240891 

lik.mk10 <- constrain(lik.mk20, q10 ~ q01)

fit.mk10 <- find.mle(lik.mk10, .1, method="subplex") 

anova(fit.mk20, mk10=fit.mk10)

     Df   lnLik    AIC     ChiSq Pr(>|Chi|)
full  2 -18.583 41.165                     
mk10  1 -18.584 39.167 0.0019786     0.9645

# The rate of evolution from non-robust to robust is similar to robust to non-robust

# Sample using MCMC

prior.exp <- make.prior.exponential(10)

samples20 <- mcmc(lik.mk20, c(.1, .1), nsteps=10000, prior=prior.exp, w=.1, print.every=1000)

samples20 <- subset(samples20, i > 1000)

mean(samples20$q01 > samples20$q10)

profiles.plot(samples20[c("q01", "q10")], col.line=terrain.colors(3), las=1, legend.pos="topright")


#The rate of evolution from non-robust to robust is similar to robust to non-robust
  #Therefore we would expect the plots to be on top of each other.

#A bit of ancestral state recon. here
fitER <- rerootingMethod(tree_10, robustmandible, model = "ER")
print(fitER)

plotTree(tree_10, setEnv = TRUE, offset = 0.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, piecol = c("blue", "red", "yellow"), cex = 0.6)
tiplabels(pie = to.matrix(robustmandible, sort(unique(robustmandible))), piecol = c("blue", "red", "yellow"), cex = 0.3)

#Doesn't seem to get that the bottom-most clade is more likely to be all red (including the nodes).