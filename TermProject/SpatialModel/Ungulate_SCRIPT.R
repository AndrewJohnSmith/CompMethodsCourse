#Run this when you first download R.
install.packages(c("geiger", "qpcR", "phytools", "ouch", "ape", "OUwie", "picante", "parallel", "surface", "auteur"))
#Once all are installed (R will do this for you) you no longer have to worry about them.

#For Centrarchids as listed as factors
Data[,4]<-as.numeric(as.character(Data[,4]))
Data[,5]<-as.numeric(as.character(Data[,5]))
Data[,6]<-as.numeric(as.character(Data[,6]))


setwd("/Users/Home/Dropbox/Amherst/Courses/CompMethods/CompMethodsCourse/TermProject")
source("Ungulate_SOURCE.R")

Files<-MakeFiles()


#Your tree and data is stored as 'Tree' and 'Data' respectively.
Tree<-Files[[1]]
Data<-Files[[2]]



Trait<-as.numeric(Data[,1])
names(Trait)<-rownames(Data)

Data[,1] <- ordered(Data[,1], levels = c("EX", "EW", "CR", "EN", "VU", "NT", "LC"))

Trees<-di2multi(Tree, tol=1e-01)
contMap(Trees,Trait,res=300, fsize=0.2)



#Standardize variables against standard length
Standardize(Data, "Measurement_1", "StLength")


#Geiger Evolutionary Models#
#Single variable
MyModels<-GeigerModels(Tree, Data, "preorb.siz")
#View results
MyModels

#OUCH Evolutionary Models#
#Single Variable
MyOUCHModels<-OuchModels(Tree, Data, "preorb.siz", "Habitat")
#View results
MyOUCHModels

#OUwie Evolutionary Models#
#Generate trees with stocastically mapped characters.
MySimTrees<-SimmapTrees(Tree, Data, "Diet", 5)

#View the first four of those character change trees.
PlotSimTreeDiet(MySimTrees, 0.6)
PlotSimTreeHabitat(MySimTrees, 0.6)

#Perform the OU analysis with your recently generated trees.
#Run the analysis over multiple cores to speed it up.
MyOUwieModels<-OUwieModelsMC(MySimTrees, Data, "preorb.siz", "Diet", "OUM")
#View results
MyOUwieModels

##Diet##
FinalOutput<-CollateOUwieNormDiet(MyOUwieModels)

##Habitat##
FinalOutput<-CollateOUwieNormHabitat(MyOUwieModels)


#Cannot use "Trait" here, need to insert the column number(s) that your data is in.
#For just a single trait use Data[,columnnumber] e.g...
DispData<-Disparity(Tree, Data, Data[,1])

#Or
#use Data[,column-number-start:column-number-end]
DispDatas<-Disparity(Tree, Data, Data[,1])


#Correlations in evolutionary rate along the tree
Rate.Shifts(Tree, Data, "ma", Ngens=1000)
Rate.Correlation(Tree, Data, "ma", "preorb.siz")
#Sometimes Trait.Correlation wont run - try a couple of different trait combinations.
Trait.Correlation(Tree, Data, "ma", "preorb.siz")

#Plotting Extras
AncestralMap(Trees, Data, "CentDataPrune")
PMorphospace(Tree, Data, Habitat, preorb.siz, ma)


#Difficult to use
ConvergentTest<-SURFACE(Tree, Data, Data[,3:4])



## :(  ##
#Problem scripts#

#Advanced OUwie models. Many of these will saddle.
MyOUwieModels<-OUwieModelsMC(MySimTrees, Data, "preorb.siz", "Diet", "OUMV")
#View results
MyOUwieModels

#Saving results - use diet if you used the diet regime.
#Use if objective function may be at a saddle point
FinalOutput<-CollateOUwieEigHabitat(MyOUwieModels)
FinalOutput<-CollateOUwieEigDiet(MyOUwieModels)

#You may be able to use this if all converged.
FinalOutput<-CollateOUwieNormDiet(MyOUwieModels)
FinalOutput<-CollateOUwieNormHabitat(MyOUwieModels)

#Descripts of the models you can run in place of "OUM".#
#model=OU1  a single peak Ornstein-Uhlenbeck model across the entire tree.
#model=OUM  a multi-peak Ornstein-Uhlenbeck model with different optima (theta) for each regime.
#model=OUMV  a multi-peak Ornstein-Uhlenbeck model with different optima (theta) and different Brownian rate parameter (sigma2) for each regime.
#model=OUMA  a multi-peak Ornstein-Uhlenbeck model with different optima (theta) and different strength of selection parameter (alpha) for each regime.
#OUMA didn't seem to run to completion.
#model=OUMVA  a multi-peak Ornstein-Uhlenbeck model with different optima (theta), different Brownian rate parameter (sigma2), and different strength of selection parameter (alpha) for each regime.
#OUMVA rarely has enough information to work.






#Tree set-up
#Data engineering
	#Read in the tree

Taxa<-Tr$tip.label

Taxa<-Taxa[c(-28, -41, -100, -117, -161, -224, -295)]


status<-NULL
for(i in 1:250){
input = gsub("_", "-", Taxa[i])
h <- htmlParse(paste("http://api.iucnredlist.org/go/", input, sep = ""))

status[i] <- xpathSApply(h, '//div[@id="red_list_category_code"]', xmlValue)
 }

status



#Remove data deficient taxa 
Poor<-Taxa[status=="DD"]
#Need to run the tree in again
Taxa<-Tr$tip.label
Taxon<-Taxa[c(28, 41, 100, 117, 161, 224, 295)]
New<-c(Poor, Taxon)

#Drop DD taxa from tree
NewTaxa<-drop.tip(Tr, New)

Tx<-NewTaxa$tip.label
#Rerun once errors are removed to check everthing removed

#Save the new files
DF<-cbind(Tx, status)
write.csv(DF, "UngulateFile.csv", quote=F)
write.tree(NewTaxa, "UngulateTree.tre")