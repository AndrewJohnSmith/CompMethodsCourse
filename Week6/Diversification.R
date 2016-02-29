# Diversification

# Based on code from Roi Holzman and Samantha Price

library(laser)
library(geiger)
library(ape)
library(apTreeshape)

#--------------------------------
# Lineage through time plot
#--------------------------------

plot(ladderize(primate_tree05), cex=0.6); axisPhylo()

#lets plot a species through time plot. we will use a plotting command from ape
ltt.plot(primate_tree05)

#under constant birth/death rates, we expect an exponential increase in species numbers with time. this should translate to a straight line on a log scale. two plotting commands enable a qualitative examination of the pattern
ltt.plot(primate_tree05,log='y')

# in ltt.plot (from the ape package), the input is the phylogeny, and the handle log='y' denotes a logarithmic scale for the y axis. this command enables control over titles, axes scale and graphic handles

# we can also look at the distribution of branching times in our tree
hist(branching.times(primate_tree05))

#--------------------------------------------------
### Simulating trees ##
#--------------------------------------------------
# simulating trees is a powerful feature of R which facilitates simulations and allow test of hypothesis and formulation of null models.
# Here we will use Geiger's sim.bdtree function
# Starting from a root node the function simulates the growth of a phylogenetic tree under a uniform, time-homogeneous birth-death process.
# If birth is greater than death, then the number of lineages is expected to grow exponentially.
# The arguments for the function include birth rate (b = Per-lineage birth (speciation) rate), death rate (d = Per-lineage death (extinction) rate), 
# the run time (t = after how many years the simulation has to stop), n = the maximum number of taxa, and a seed. 
# the latter is useful in case you want to run your simulations again and get the exact same result. 
# lastly, the function can be set up to return trees with 'all extinct' taxa if extinct=T

# We can start with a pure-birth tree (where d = 0)

simulated.tree1<-sim.bdtree(b=0.2, d=0, t=10)
plot(ladderize(simulated.tree1))

# its interesting to note the variation in tree shape and # species:
# we will plot 9 trees on one panel using the par command to prepare the graphic object

par(mfrow=c(3,3))

#then simulate 9 trees under the same parameters; notice that this code includes a loop function

for (i in 1:9){
	simulated.tree2<-sim.bdtree(b=0.2, d=0, t=10)
	plot(simulated.tree2)
	}

# We can add a death rate and check out the tree. The extinct=F flag prevents trees with no survivors
for (i in 1:9){
simulated.tree3<-sim.bdtree(b=0.2, d=0.05, t=15, extinct=F)
plot(simulated.tree3)
}

# note that not all taxa get to the current time- these are extinct taxa, that can be removed with
par(mfrow=c(1,1))

simulated.tree4=drop.extinct(simulated.tree3)

plot(simulated.tree4)






