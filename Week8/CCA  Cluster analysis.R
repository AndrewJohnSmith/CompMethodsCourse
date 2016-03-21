
#---------------------------------------------
# Canonical Correspondence Analysis
#---------------------------------------------

library (vegan)

# Import the environmental variables. These will be the predictors.
africa_environ <- read.csv("africa_environ.csv")

View(africa_environ)

# Import the species presence/absence matrix. This will serve as the "dependent matrix"
africa_primates <- read.csv("africa_primates.csv")

View(africa_primates)

# The model structure is similar to a linear model, with the dependent data matrix (in this case species presence/absence) to the left of ~
# and the predictors to the right. You may notice that not all the environmental variables are being used.
# Also, the "data" argument only refers to the predictor data set.
cca.model1 <- cca(africa_primates ~ TotalRain + MinTempMean + MaxTempMean, data=africa_environ)

# get some results
cca.model1

#Call: cca(formula = africa_primates ~ TotalRain + MinTempMean + MaxTempMean, data = africa_environ)
#
#              Inertia Proportion Rank
#Total          4.5024     1.0000
#Constrained    0.6757     0.1501    3
#Unconstrained  3.8267     0.8499   28
#Inertia is mean squared contingency coefficient
#
#Eigenvalues for constrained axes:
#  CCA1   CCA2   CCA3
#0.4273 0.1786 0.0698
#
#Eigenvalues for unconstrained axes:
#   CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8
#0.6467 0.5335 0.4732 0.3293 0.2849 0.2657 0.2273 0.1880
#(Showed only 8 of all 28 unconstrained eigenvalues)

# You can create a biplot of the CCA model
row.names(cca.model1$CCA$wa)<-as.character(africa_environ$Site) #Add in site names
plot(cca.model1)


#############
# You can test for the significance of the model
#############

anova.cca(cca.model1)

#Permutation test for cca under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: cca(formula = africa_primates ~ TotalRain + MinTempMean + MaxTempMean, data = africa_environ)
#         Df ChiSquare     F Pr(>F)
#Model     3    0.6757 2.472  0.001 ***
#Residual 42    3.8267
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1

#############
# You can also examine the effects of individual predictors in the model. This test is sequential: the terms are analysed in the order they happen to be in the model.
#############

anova(cca.model1, by="term", permutations=999)

#Permutation test for cca under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999
#
#Model: cca(formula = africa_primates ~ TotalRain + MinTempMean + MaxTempMean, data = africa_environ)
#            Df ChiSquare      F Pr(>F)
#TotalRain    1    0.3939 4.3227  0.001 ***
#MinTempMean  1    0.0943 1.0351  0.576
#MaxTempMean  1    0.1875 2.0583  0.052 .
#Residual    42    3.8267
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1

###############
# Here is another way to examine the effects of individual predictors in the model based on marginal effects (?Type III effects?). This is a better approach thatn the previous way
# because the order of the predictor in the model doesn't matter.
###############

anova(cca.model1, by="mar", permutations=999)

#Permutation test for cca under reduced model
#Marginal effects of terms
#Permutation: free
#Number of permutations: 999
#
#Model: cca(formula = africa_primates ~ TotalRain + MinTempMean + MaxTempMean, data = africa_environ)
#            Df ChiSquare      F Pr(>F)
#TotalRain    1    0.2301 2.5257  0.014 *
#MinTempMean  1    0.1722 1.8905  0.091 .
#MaxTempMean  1    0.1875 2.0583  0.044 *
#Residual    42    3.8267
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1

#Basically, African primate position (or abundance) is influenced by total rain and maximum temperature

########
# Also, it is possible to analyze significance of each axis
########

anova(cca.model1, by="axis", permutations=999)

#Permutation test for cca under reduced model
#Marginal tests for axes
#Permutation: free
#Number of permutations: 999
#
#Model: cca(formula = africa_primates ~ TotalRain + MinTempMean + MaxTempMean, data = africa_environ)
#         Df ChiSquare      F Pr(>F)
#CCA1      1    0.4273 4.6897  0.001 ***
#CCA2      1    0.1786 1.9600  0.067 .
#CCA3      1    0.0698 0.7665  0.833
#Residual 42    3.8267
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1

#First axis is really driving the trend (and to some extent the second axis)

##################################################################################

#---------------------------------------------
# Cluster Analysis
#---------------------------------------------

# Use the Euclidian distance matrix

# Create the data matrix using the primate traits data from DFA

primate_traits_dfa <- read.csv("primate_traits_dfa.csv")

View(primate_traits_dfa)

primate_traits_cluster <- primate_traits_dfa[,3:6]

# Create a Euclidian distance matrix

primate_traits_dist <- dist(primate_traits_cluster)

# View the distance matrix
primate_traits_dist

# Run the cluster function, which creates a hierarchical cluster; it uses the complete linkage method by default; also known as farthest neighbour clustering
cluster01 <- hclust(primate_traits_dist)

# View the dendrogram
plot(cluster01)

#----------------------------------------------------
# Run a UPGMA cluster analysis; this type is usually the best to use
#----------------------------------------------------

# Install and load the cluster package

library(cluster)

cluster_upgma01 <- agnes(primate_traits_dist, method="average")
# This package can also use a raw data matrix and calculate various distance metrics within the function; it can also use different methods to contruct the cluster analysis; refer to the package manual for details

plot(cluster_upgma01)

# For a brief definition of different hierarchical clustering methods see:
# http://nlp.stanford.edu/IR-book/completelink.html

#----------------------------------------------------------------------
# A cluster analysis with bootstrap support
#----------------------------------------------------------------------

# Use the recluster package

install.packages("recluster")

library(recluster)

# Read in data

africa_primates2 <- read.csv("africa_primates2.csv")

rownames(africa_primates2) <- africa_primates2$Site

africa_primates3 <- africa_primates2[,-1]

# This function creates a series of trees by resampling the order of sites in the dissimilarity matrix.
# Then, it computes a consensus among them. The resulting tree is unaffected by original row order.
# tr = The number of trees to be used for the consensus;  p = proportion for a clade to be represented in the consensus tree;
# dist = type of distance metric to use; method = the one specified is UPGMA; phylo = only needs to be specified if conducing a phylogenetic betadiversity analysis

africa_primates2_cluster <- recluster.cons(africa_primates3, phylo = NULL, tr = 100, p = 0.5, dist = "sorensen", method = "average", blenghts=TRUE, select=FALSE)

# Bootstraps the consensus tree; boot = number of replicates; level = use 1 to keep the same number of sites in the replicates

africa_primates2_boot01 <- recluster.boot(africa_primates2_cluster$cons, africa_primates3, phylo = NULL, tr = 100, p = 0.5, dist = "sorensen", method = "average", boot = 100, level = 1)

# Plot the bootstrapped tree

recluster.plot(africa_primates2_cluster$cons,africa_primates2_boot01,direction="downwards")

##PCA on the environmental data
AfEnvironPCA<-prcomp(africa_environ[,2:6])
summary(AfEnvironPCA)
AfEnvScore<-AfEnvironPCA$x

cca.model2 <- cca(africa_primates ~ PC1 + PC2, data=as.data.frame(AfEnvScore))

row.names(cca.model2$CCA$wa)<-as.character(africa_environ$Site) #Add in site names
plot(cca.model2)
