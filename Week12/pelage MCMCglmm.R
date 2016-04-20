# Project using a phylogenetic Markov Chain Monte Carlo GLMM to test whether primate coloration
# blends into background color and does it differ across different parts of the body
# and across species with different types of color vision

# Location on the body and visual system are fixed effects and species is designated as a random effect
# Species is a random effect because each species was measured multiple times (at different locations on their body)


library("ape")
library("MCMCglmm")

#################################
#####IMPORT DATA & PHYLOGENY#####
#################################
data <- read.csv("data for comparative analyses.csv",header=TRUE)
str(data)

tree <- read.nexus("pelage project consensus tree.nex")
str(tree)
plot.phylo(tree)
is.binary.tree(tree) #true
is.ultrametric(tree) #true

##############################
#####STATISTICAL ANALYSES#####
##############################
#PRIOR
prior <- list(R=list(V=1,nu=0,fix=1),G=list(G1=list(V=1,nu=0)))

#MODEL INCLUDING LOCATION & VISION
model.LocVis <- MCMCglmm(Overlap~Location+Vision,random=~Species,family="categorical",pedigree=tree,data=data,
                         prior=prior,nitt=1050000,thin=1000,burnin=50000)
plot(model.LocVis$Sol)
plot(model.LocVis$VCV)
autocorr(model.LocVis$VCV)

summary(model.LocVis)


#Location effects: Overlap ~ Location + Vision 
#
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
#(Intercept)     14.6147   7.6159  22.5870    50.65 <0.001 ***
#LocationBelly    2.6612  -0.8392   5.8616   821.48  0.078 .  
#LocationCap     -1.1843  -3.2348   1.1627  1000.00  0.288    
#LocationTail    -1.5398  -3.6721   0.5147  1000.00  0.140    
#VisionTetra     -4.4978  -7.7914  -1.4535   655.35  0.002 ** 
#VisionTri       -3.8291  -7.1597  -0.6499   651.26  0.002 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
	#When using discrete variables the intercept is actually the 'back' (study looks at hair on the back).
	#The negative or positive values are relative to the 'back' value of 14.6.
		#Therefore the posterior mean for VisionTri is ~11, or -3.8 below the 'back' value of 14.6.