# Ordination and related methods in R

# Some text and code from Dolph Schluter - https://www.zoology.ubc.ca/~schluter/R/multivariate/

#-----------------------------------------------
# Principal Components Analysis
#-----------------------------------------------

# read dataset into R
setwd("/Users/Home/Dropbox/Amherst/Courses/CompMethods/CompMethodsCourse/Week7")

pca_data1 <- read.csv("pca_data1.csv")

# If your variables are measured on different scales then you should use a
# correlation matrix to extract components

# If species or sample names are in the first column then you must specify the columns to use
# in this case, columns 2 through 6 (unless you set rownames)
# scale = TRUE sets the analysis to use a covariation matrix
pca_model01 <- prcomp(pca_data1 [,2:6], scale=TRUE)

# A biplot displays the data points along two principal components (the first and second components, by default)
# and adds arrows to indicate the contributions of each trait to these principal components.
# The graphs can get messy if there are too many variables.
# Use the cex option as well as the xlim and ylim options to help fit the labels onto the graph.
#Rotation is the same as loadings

biplot(pca_model01, cex = 0.7)

# Change the scale of the X and Y axes
biplot(pca_model01, cex = 0.7, xlim=c(-0.4, 0.6), ylim=c(-0.4, 0.6))

# Save your  biplot as a jpg

# Plot components 3 and 4

biplot(pca_model01, cex = 0.7, choices=c(3,4))

# Results are extracted from the prcomp object, e.g.
summary(pca_model01)                 # variance explained by the components

plot(pca_model01,type="lines")       # "scree" plot of eigenvalues

screeplot(pca_model01, type="lines") # same
#Greater than 1 indicates the new axis is contrubuting to variation more so than the original variable.
#Less than 1 indicates variable is less than the original variation.

pca_model01$rotation                 # eigenvectors for each component(with trait loadings)

pca_model01$sdev^2                   # eigenvalues for each component(variances)

predict(pca_model01)                 # PCA scores
  pca_model01$x #can also do this
# export PCA scores to csv file

write.csv(predict(pca_model01), "pca_model01_scores.csv")

#-----------------------------------------------------
# Non-metric multidimentional scaling
#-----------------------------------------------------

# This method preserves only the rank order of distances between points, as much as possible given that the data have been reduced to a small number of dimensions. 
# The command is in the MASS library, so load it to begin

library(MASS)

# Read in your data
mada_occur_data02 <- read.csv("mada_occur_data02.csv")

# Remove the site label column
mada_occur_data10 <- mada_occur_data02[,-1]

# Convert the presence - absence matrix to a Euclidian distance matrix
mada_dist10 <- dist(mada_occur_data10)

# Run the analysis
NMDS01 <- isoMDS(mada_dist10, k = 2)
# K indicates the number of axes you wish to return

# Output the results
NMDS01

# Plot the results 
plot(NMDS01$points, xlab="Coordinate 1", ylab="Coordinate 2", main="Nonmetric MDS")
text(NMDS01$points, labels = mada_occur_data02$site, cex=.7, pos=4)

#-------------------------------------------------------
# Correspondence analysis
#-------------------------------------------------------

library(MASS)

#Use the mada_occur_data10 data object
#To analyze, put the species as separate columns of a data frame and the sites as separate rows.
#The correspondence analysis is carried out as follows.
  #nf refers to the number of axes to be extracted, usually just 1 or 2. Results are stored in the correspondence analysis object

correspond_model01 <- corresp(mada_occur_data10, nf = 2)

# If nf = 2, the following command creates a scatter plot of both sites and species, with their scales arranged to correspond.

# Species that are close to one another in the plot tend to occur at the same sites. Likewise, nearby sites in the plot will tend to have similar species compositions (or species abundances, if abundance data are analyzed).
#Species located next to sites in the graph occur predominantly there, whereas species falling between sites tend to occur in two or more sites.

plot(correspond_model01)

# The column and row scores indicate the contributions of individual species and sites to the correspondence axes.

correspond_model01$rscore
  
correspond_model01$cscore


#---------------------------------------------
# Phylogenetic principal components analysis
#---------------------------------------------

library(phytools)

primate_tree05 <- read.tree("primate_tree05.txt") 
primate_traits <- read.csv("primate_traits.csv")
rownames(primate_traits) <- primate_traits$species
primate_traits_m <- data.matrix(primate_traits)
primate_traits_m2 <- primate_traits_m[,-1]

# There is a zero in the dataset so we have to add a constant while log transforming
primate_traits_m2 <- log((primate_traits_m2)+0.01)

#Uses lambda
primate_traits_phyl_pca <- phyl.pca(primate_tree05, primate_traits_m2, method="lambda", mode="corr")

#Can also use BM
primate_traits_phyl_pca <- phyl.pca(primate_tree05, primate_traits_m2, method="BM", mode="corr")

phylomorphospace(primate_tree05, primate_traits_phyl_pca$S, label="horizontal")

# Output options:

# A list with the following components:
# Eval 	diagonal matrix of eigenvalues.
# Evec 	matrix with eigenvectors in columns.
# S 	matrix with scores.
# L 	matrix with loadings.
# lambda 		fitted value of lambda (method="lambda" only).
# logL 	log-likelihood for lambda model (method="logL" only).

primate_traits_phyl_pca$Eval

primate_traits_phyl_pca$Evec

primate_traits_phyl_pca$S

primate_traits_phyl_pca$lambda

#***************************************
# Use two of the above methods on your data (or a subet of your data)
# Email the results to me
#***************************************