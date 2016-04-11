# DFA and MANOVA

setwd("~/Dropbox/Amherst/Courses/CompMethods/CompMethodsCourse/Week9")

#-----------------------------------------------
# Discriminant function analysis
#-----------------------------------------------

library(MASS)

# read the data into R   (these are log10 transformed data)

primate_traits_dfa <- read.csv("primate_traits_dfa.csv")

# Set the rownames as the species names

rownames(primate_traits_dfa)=primate_traits_dfa$species

# Examine the mean and standard deviation of each variable

sapply(primate_traits_dfa[3:6],mean)  # the sapply function can be used to apply some other function to each column in a data frame, in this case the mean function

sapply(primate_traits_dfa[3:6],sd)

#This appears almost identical to the apply(DATA[,1:n], 2, FUN) argument

# Run the DFA model

dfa_model01 = lda(clade ~ abs_female_ECV + female_mass + groupsize + per_leaves, data = primate_traits_dfa)

# You can plot the results
plot(dfa_model01)

# Produce some of the model results

print(dfa_model01)

# Output:
Call:
lda(clade ~ abs_female_ECV + female_mass + groupsize + per_leaves,
    data = primate_traits_dfa)

Prior probabilities of groups:    # This is based on the number of samples per group
     cata    platyr     strep
0.5161290 0.3225806 0.1612903

Group means:
       abs_female_ECV female_mass groupsize per_leaves
cata         2.019215    3.946186 1.2866047   1.352307
platyr       1.709137    3.571151 1.1011565   1.204701
strep        1.466894    3.575161 0.8952527   1.417832

Coefficients of linear discriminants:    # Describes how the original variables are related to the discriminant function axes
                      LD1        LD2
abs_female_ECV -9.6796097  4.4403014
female_mass     6.1609366 -5.7604476
groupsize      -0.2728523 -0.6223835
per_leaves     -0.5676214  0.5161332

Proportion of trace:       # Variance explained by each DFA axis
   LD1    LD2
0.8884 0.1116

# You can also produce the discriminant function scores for each sample

predict(dfa_model01)$x
#Any reason why it gives only 2 axes from four variables?

# Output
                                 LD1         LD2
Alouatta_caraya           0.33524595  0.43622481
Alouatta_palliata         0.96307938 -0.22814749
Alouatta_pigra            1.53272684 -0.49650156
Alouatta_seniculus        0.49500242  0.22749848
Ateles_belzebuth         -1.01266468 -0.31584744
Ateles_geoffroyi         -1.42458318  0.12539256
Cebus_apella             -1.91013844  1.80900751
Cercocebus_agilis        -1.82909688  0.48387690
Cercopithecus_ascanius   -0.95755024  1.11966695
Cercopithecus_mitis      -0.49811090  0.41383453
Chlorocebus_aethiops     -0.79627753  0.91882368
Colobus_guereza           0.28961305 -0.16219735
Erythrocebus_patas       -0.34258719 -0.40980227
Eulemur_fulvus            1.68210981  0.60104243
Gorilla_gorilla          -1.31104279 -2.54078074
Hylobates_lar            -1.76458334  1.36191636
Indri_indri               3.23246946 -1.09528865
Lagothrix_lagotricha     -0.99998591 -0.08035621
Lemur_catta               2.43292825  0.33101348
Lophocebus_albigena      -0.79594579  0.09277325
Macaca_sylvanus          -0.15365486 -0.91262286
Nasalis_larvatus          0.24156315 -0.60812618
Pan_troglodytes          -1.94903454 -1.97126513
Papio_hamadryas          -1.45565795 -0.78231642
Piliocolobus_badius       0.35754396 -0.56808856
Pithecia_pithecia        -0.07302894  1.98675097
Propithecus_diadema       2.30687244 -0.70064701
Saguinus_fuscicollis      2.63641642  2.17053557
Semnopithecus_entellus   -0.47369856 -0.70459235
Symphalangus_syndactylus -0.89822925  0.11995498
Varecia_variegata         2.14029987 -0.62173223

# You can assess the accuracy of the prediction by examining the percent correct for each category. This time, re-run the model using the CV=TRUE argument, which generates jacknifed (i.e., leave one out) predictions.

dfa_model02 = lda(clade ~ abs_female_ECV + female_mass + groupsize + per_leaves, data = primate_traits_dfa, CV=TRUE)

jacknifed_pred <- table(primate_traits_dfa$clade, dfa_model02$class)

diag(prop.table(jacknifed_pred, 1))

# Percent correctly predicted for all samples
#cata platyr  strep 
#0.75   0.30   0.80 
  #Some of these are pretty shitty.

sum(diag(prop.table(jacknifed_pred)))

# Now examine how well each species was classifed
  #Many don't have that awesome percentages (i.e., <50% certainty)
dfa_model02$posterior

# You can write these results to a csv file

write.csv(dfa_model02$posterior, "posterior_prob.csv")

#******************************
# Email me an example of a poorly classifed catarrhine, platyrrhine, and strepsirrhine species
#******************************

# Classify new or unknown samples according to the DFA model

# Read in a data file of unclassifed taxa (e.g., a fossil taxon)
primate_traits_dfa_unk <- read.csv("primate_traits_dfa_unk.csv")

#Change to CV = FALSE tp give an 'lda' object
dfa_model03 = lda(clade ~ abs_female_ECV + female_mass + groupsize + per_leaves, data = primate_traits_dfa, CV=F)

dfa_model03_class <- predict(dfa_model03, primate_traits_dfa_unk)

# Obtain the classification results
  #Rows indicate the probability that the 3 'fossil' taxa fit into one of the primate groupings

dfa_model03_class$posterior

# Obtain the DFA scores for the unknown samples

dfa_model03_class$x

#Plot them
plot(dfa_model01)
points(dfa_model03_class$x, pch=20, col="red")

#------------------------------------------------------
# Run a MANOVA to obtain the p value of the DFA analysis (which should be identical to a DFA except the idependent and dependent variables are switched)
#------------------------------------------------------

# Extract the dependent variables from the dataset
dep_var_matrix01=primate_traits_dfa[,3:6]

# Transform the dataset into a matrix

dep_var_matrix01_m=as.matrix(dep_var_matrix01)

# Run the MANOVA model

manova_model01 <- manova(dep_var_matrix01_m ~ primate_traits_dfa$clade)

# Obtain the model's p value based on Wilks' lambda

summary(manova_model01, test="Wilks")
#Crazy significant - the four variables are different depending on the clade.

##############################################################

