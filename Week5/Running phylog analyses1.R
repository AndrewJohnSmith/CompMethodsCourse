#-----------------------------------------------
# Running phylogenetic analyses 1
#-----------------------------------------------
# Lots of code borrowed from Anthrotree wiki
#-----------------------------------------------

# Load the packages: phytools and caper

library(phytools)

library(caper)

setwd("/Users/Home/Dropbox/Amherst/Courses/CompMethods/CompMethodsCourse/Week5")

# Import a phylogenetic tree; flex_species_tree.txt

flex_data_tree<-read.tree("flex_species_tree.txt")

# Import your dataset; flex_data_final.csv
# Use the species names in your dataset as your row names
flex_data_final <- read.csv("flex_data_final.csv", row.names=1)

# Do your taxon names match in the phylogeny and dataset?
# Use the name.check function in the geiger package
# Load geiger

library(geiger)

name.check(flex_data_tree, flex_data_final)

# This returns output for species that are in your tree but not your dataset
# And species in your dataset that are not in your tree

# Drop a taxon from the phylogeny
flex_data_tree_no_cc <- drop.tip(flex_data_tree, "Cebus_capucinus")

# View your tree
plot(ladderize(flex_data_tree_no_cc), cex=0.5)
axisPhylo()

# Computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.
cophenetic(flex_data_tree)["Papio_hamadryas", "Propithecus_diadema"]

# Computes the pairwise distances among all taxa
cophenetic_all <- cophenetic(flex_data_tree)

# This output is large so write it to a csv file to view

write.csv(cophenetic_all, "cophenetic_all.csv")

# Open the csv file in Excel


#-----------------------------------------------------------------------------
# Calculating Phylogenetically Independent Contrasts
#-----------------------------------------------------------------------------

# use caper
library(caper)

# Import phylogey; this file contains more than one tree
tree_list <- read.nexus("ProbSet5.tree.nex")

# Import associated dataset
ds <- read.table("ProbSet5.data.txt", header=T)

# We will use the second tree only
tree <- tree_list[[2]]

plot(ladderize(tree))
axisPhylo()

# Create a comparative data object using the caper package
# This object is used by caper and contains both a phylogeny
# and an associated dataset

primate <- comparative.data(phy = tree, data = ds, names.col = Species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

# Some notes about the comparative data object: it matches the taxa in your dataset and the taxa in your phylogeny automatically
# vcv=TRUE indicates that a variance-covariance matrix is created
# na.omit = FALSE stops the function from removing species without data for certain variables
# warn.dropped = TRUE will tell you if there are any species which are not in the tree and the dataset and are therefore dropped from the data object
# You can see the list of dropped species by typing primate$dropped
# Make sure you check this list is what you expected, it may reveal typos in your species names

primate$dropped
#Everything was there

# The crunch function calculates phylogenetically independent contrasts
model.ic <- crunch(log(HomeRange) ~ log(FemMass), data = primate)

summary(model.ic)

# You can examine the contrasts using the function 'caic.table' with the model as output
caic.table(model.ic)

# The table provides the standardized contrasts for each of the variables, along with the 
# variance (the branch length connecting the species or nodes), the number of descendent lineages 
# from each node, the node depth, and the studentized residuals (used to identify outliers).

# Examine the model diagnostics

par(mfrow=c(2,2)) #so you can see all 4 plots at once
plot(model.ic)
par(mfrow=c(1,1)) #to reset the graphic window to just one graph

# caper has an inbuilt option that automatically removes outliers with studentized residuals >?3 (Jones and Purvis 1997), using the argument "robust=TRUE"

model.ic2<-crunch(log(HomeRange) ~ log(FemMass), data = primate, robust = TRUE)
summary(model.ic2)

# Note that one outlier has automatically been removed


#---------------------------------------------------------------------------
# Running a phylogenetic generalized least squares model (PGLS) using the caper package
#---------------------------------------------------------------------------

library (caper)

# load the tab delimited file

primatedata <-read.table("Primatedata.txt", sep = "\t", header = TRUE)
head(primatedata)

# load the nexus file

primatetree <-read.nexus("consensusTree_10kTrees_Version2.nex")

# Note that phylogenies in R cannot have spaces in the tip labels, instead R uses underscores to separate Genus and species (Genus_species). The names of the species in the tree must match those in the data, therefore in this dataset the spaces in species names have been replaced with underscores.

primatedata$Binomial<-gsub(" ", "_", primatedata$Binomial)

# Create the comparative data object
primate <- comparative.data(phy = primatetree, data = primatedata, names.col = Binomial, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
primate$dropped

# NOTE regarding the arguments to this function: vcv = TRUE stores a variance covariance matrix of your tree, which you will need for pgls. na.omit = FALSE stops the function from removing species without data for certain variables. warn.dropped = TRUE will tell you if there are any species which are not in the tree and the dataset and are therefore dropped from the data object.

# If you drop species, you can see the list of dropped species by typing primate$dropped. Make sure you check this list is what you expected, it may reveal typos in your species names. If you want to turn this off use warn.dropped = FALSE. Note that the comparative.data function does all the matching of data and tree for you. Note that we expect to drop some species here because not all the species in primatetree are in primatedata and vice versa.

# The function for PGLS analyses in caper is pgls. To fit a model of log gestation length against log body mass which uses the maximum likelihood estimate of lambda we use the following code:

model.pgls<-pgls(log(GestationLen_d)~ log(AdultBodyMass_g), data = primate, lambda ='ML')
  # note that lambda can also equal a fixed number, e.g. you can use lambda = 1 if you want a model that is equal to using independent contrasts

summary(model.pgls)

# The relatoinship between the variables can be plotted

plot(log(GestationLen_d)~ log(AdultBodyMass_g), data = primatedata, pch=20)
abline(model.pgls, col="red", lwd=2)

par(mfrow=c(2,2))  # Change your graphics window
plot(model.pgls) # Check the model diagnostics
par(mfrow=c(1,1))#to reset the graphic window to just one graph

# Check if you have any sample that have a phylogenetic standardized residual > 3 or < -3

res1 <- residuals(model.pgls, phylo = TRUE)

res1 <- res1/sqrt(var(res1))[1]

rownames(res1)<-rownames(model.pgls$residuals)

# This will list any taxa that have large residuals, if any
rownames(res1)[(abs(res1)>3)]

# You can then remove these species from the data and tree and redo the analysis.

primate_nooutliers<-primate[-which(abs(res1)>3),]

model.pgls_nooutliers<-pgls(log(GestationLen_d)~ log(AdultBodyMass_g), data = primate_nooutliers, lambda='ML')

summary(model.pgls_nooutliers)
#Kappa and delta are branch length transformers.
#The kappa parameter differentially stretches or compresses individual phylogenetic branch lengths and can be used to test for a punctuational versus gradual mode of trait evolution. Kappa > 1.0 stretches long branches more than shorter ones, indicating that longer branches contribute more to trait evolution (as if the rate of evolution accelerates within a long branch). Kappa < 1.0 compresses longer branches more than shorter ones. In the extreme of Kappa = 0.0, trait evolution is independent of the length of the branch. Kappa = 0.0 is consistent with a punctuational mode of evolution.
#The parameter delta scales overall path lengths in the phylogeny – the distance from the root to the species, as well as the shared path lengths. It can detect whether the rate of trait evolution has accelerated or slowed over time as one moves from the root to the tips, and can find evidence for adaptive radiations.  If the estimate of Delta < 1.0, this says that shorter paths (earlier evolution in the phylogeny) contribute disproportionately to trait evolution – this is the signature of an adaptive radiation: rapid early evolution followed by slower rates of change among closely related species. Delta > 1.0 indicates that longer paths contribute more to trait evolution. This is the signature of accelerating evolution as time progresses. Seen this way, delta is a parameter that detects differential rates of evolution over time and re-scales the phylogeny to a basis in which the rate of evolution is constant.
#The parameter lambda reveals whether the phylogeny correctly predicts the patterns of covariance among species on a given trait. This important parameter in effect indicates whether one of the key assumptions underlying the use of comparative methods – that species are not independent – is true for a given phylogeny and trait. If a trait is in fact evolving among species as if they were independent, this parameter will take the value 0.0 and indicate that phylogenetic correction can be dispensed with. A lambda value of 0.0 corresponds to the tree being represented as a star or big-bang phylogeny.  If traits are evolving as expected given the tree topology and the random walk model, lambda takes the value of 1.0.  Values of lambda = 1.0 are consistent with the constant-variance model (sometimes called Brownian motion) being a correct representation of the data. Intermediate values of lambda arise when the tree topology over-estimates the covariance among species.


# We can also plot the likelihood surface of lambda. If the surface is relatively flat or bimodal then
# lambda may be relatively unstable

profile_lambda <- pgls.profile(model.pgls, which="lambda") # vary lambda
plot(profile_lambda)

#You want a single peak - multiple peaks or a flat line is bad and is typically caused by low sample size
# There are cases where you will get an error related to the optimization of lambda
# If so, then you may need to adjust the bounds on the search of the maximum likelihood space
# The default values for the bounds are lambda: 1e-6 to 1. 
# To change the bounds, use the "bounds" argument within your pgls model, for example:

bnds <- list(lambda = c(0.01, 1), kappa = c(0, 3), delta = c(0, 3))
#that changes the boundaries of lambda to a minimum of 0.01. Don't worry about the other parameters now.

model.pgls2<-pgls(log(GestationLen_d)~ log(AdultBodyMass_g), data = primate, lambda ="ML", bounds = bnds)

#Don't really need to change the bounds as optim worked, lambda doesn't change as I haven't changed the data
profile_lambda2 <- pgls.profile(model.pgls2, which="lambda") # vary lambda
plot(profile_lambda2)

#-----------------------------------------------------------------
# Phylogenetic RMA regression using Phytools
#-----------------------------------------------------------------

# Load data and tree

mass_brain2 <- read.csv("mass_brain2.csv", row.names=1)
head(mass_brain2)

tree_new <- read.tree("tree_new.txt")

# be sure that the taxa are the same the dataset and tree (may need to prune tree)

name.check(tree_new, mass_brain2)

# separate each variable into its own object
brain <- mass_brain2$log_ecv_female
mass <- mass_brain2$log_mass_female

# Add names to each vector
names(brain) <- rownames(mass_brain2)
names(mass) <- rownames(mass_brain2)

# Make Phylogenetic RMA model
rma_model1 <- phyl.RMA(mass,brain,tree_new,method="lambda")

# Produce all the output (see below for more details); you may want to write this to a csv file

rma_model1

# RMA.beta a vector of RMA regression coefficients.
# V a VCV matrix for the traits.
# lambda fitted value of lambda (method="lambda" only).
# logL log-likelihood (method="lambda" only).
# test a vector containing results for hypothesis tests on beta.
# resid a vector of residuals for y given x.
              

#-----------------------------------------------------------------
# Phylogenetic ANOVA using pgls function in caper
#-----------------------------------------------------------------

brain_tree <- read.tree("flex_species_tree.txt")

brain_data <- read.csv("flex_data_final.csv")
head(brain_data)

brain_comp <- comparative.data(phy = brain_tree, data = brain_data, names.col = species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

phyl_anova01 <- pgls(log10(abs_female_ECV)~catarrhine, data = brain_comp, lambda ="ML")

summary(phyl_anova01)

#Can also use the phyANOVA function in phytools

#-----------------------------------------------------------------
# Calculating phylogenetic signal using the phytools package
#-----------------------------------------------------------------

# Load dataset and tree

flex_species_tree<-read.tree("flex_species_tree.txt")

socorgflex <- read.csv("soc_org_flex_phyl_signal.csv")
socorgflex2 <- socorgflex$SocOrgFlex
socorgflex2 <- matrix(socorgflex2)
rownames(socorgflex2) <- socorgflex$Species

# Using Blomberg's K; output yields K value and p value
phylosig(flex_species_tree, socorgflex2, method="K", test=TRUE, nsim=1000)
#$K
#[1] 0.363947
#$P
#[1] 0.034
  #Significantly weak signal

# Using Pagel's lambda; Output yields lambda value, the optimized likelihood (logL),
# the likelihood if lambda=0 (logL0), and the p value (tests if optimized lambda is greater than 0)

phylosig(flex_species_tree, socorgflex2, method="lambda", test=TRUE, nsim=1000)

#-------------------------------------------------------------------------------------
# Plot the trait on the tree using contMap function in phytools
#-------------------------------------------------------------------------------------

# Extract the variable as a vector

socorgflex2=socorgflex$SocOrgFlex

# Assign a name to each value in the vector; based on the species names in the dataset

names(socorgflex2)=socorgflex$Species

# Create the plot

obj <- contMap(flex_species_tree, socorgflex2, fsize=0.5)

# Can also create different types of plots

## plot leftward
plot(obj, direction = "leftwards")

## plot fan style tree
plot(obj, type = "fan", lwd = 5)  

# More plotting options can be found in Liam Revell's file: Plotting phylogenies.docx

#----------------------------------------------------------------

# Some useful info here:   https://www.zoology.ubc.ca/~schluter/R/phylogenetic/

# and here:   http://bodegaphylo.wikispot.org/Phylogenetics_and_Comparative_Methods_in_R



