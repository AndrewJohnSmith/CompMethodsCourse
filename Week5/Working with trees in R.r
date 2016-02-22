#----------------------------------------------
# Working with trees in R
#----------------------------------------------
# Some code borrowed from: https://www.zoology.ubc.ca/~schluter/R/phylogenetic/ and
# Phylogenetic comparative methods wiki at NEScent http://www.r-phylo.org/wiki/Main_Page
#----------------------------------------------------------------------------------------

# install the packages: ape, geiger, phytools

# Input a phylogeny into R; download a consensus tree of your choice from the 10K Trees website
# 10K Trees will save it as a nexus file. Try to convert it to Newick format using Figtree.

# For a Newick formatted tree use; this file format can also use a .txt extension
mytree <- read.tree("MyNewickTreefile.tre")
# For a Nexus formatted tree use
mytree <- read.nexus("MyNexusTreefile.nex")

# some of the information stored in the phylo object "mytree"
mytree               # prints basic tree information
mytree$tip.label     # species names
mytree$tip           # shows the numbering of the tips
mytree$edge.length   # lengths of all tree branches
mytree$edge          # identity of branches (from node, to node)


# Use the plot.phylo command to view a phylogeny stored in the object "mytree"
# Include the "cex" option to help reduce overlap of species labels
plot.phylo(mytree, cex=0.6)

# the plot.phylo function has several option for displaying trees
# use help (plot.phylo) to see options and try to modify your tree

# Resolving polytomies

# break polytomies in random order.
dichotomousphylogeny <- multi2di(phylogeny, random = TRUE)

# collapse these zero length branches into 'true multichotomies'
dichotomousphylogeny <- multi2di(phylogeny, random = FALSE)

# Check if your phylogeny is fully resolved (i.e. no polytomies).
is.binary.tree(dichotomousphylogeny)


# Write Newick formatted tree
write.tree(mytree_edit1, file="MyNewickTreefile_edit1.tre")
#  OR
write.tree(mytree_edit1, file="MyNewickTreefile_edit1.txt")

# Write Nexus formatted tree
write.nexus(mytree_edit2, file="MyNexusTreefile_edit2.nex")


