# Ancestral State Reconstruction

# Reconstructing Continuous Traits

library(ape)

#-------------------------------------------------
# Using the ape package
#-------------------------------------------------

MLreconstruction2 <- ace(log(female_mass), primate_tree05, type="continuous", method="ML")

MLreconstruction2

    Ancestral Character Estimation

Call: ace(x = log(female_mass), phy = primate_tree05, type = "continuous", 
    method = "ML")

    Log-likelihood: 101.879 

$ace
       32        33        34        35        36        37        38 
 8.474102  8.557320  9.102226  8.964412  8.794472  8.513742  8.355311 
       39        40        41        42        43        44        45 
 8.502138  8.854828  8.862567  9.000189  9.023924  8.980581  9.089969 
       46        47        48        49        50        51        52 
 9.605884  9.070550 10.587464  7.853360  7.842406  8.339132  8.529147 
       53        54        55        56        57        58        59 
 8.639080  8.498902  8.665969  8.932169  7.739842  8.347854  8.583144 
       60        61 
 8.065168  7.922234 

$sigma2
[1] 0.012327438 0.001617496

$CI95
        [,1]      [,2]
32  7.718239  9.229965
33  8.014990  9.099650
34  8.693271  9.511182
35  8.638564  9.290260
36  8.526330  9.062613
37  8.271690  8.755793
38  8.117924  8.592697
39  8.257079  8.747198
40  8.593592  9.116063
41  8.599386  9.125749
42  8.761587  9.238792
43  8.730196  9.317653
44  8.678429  9.282734
45  8.787843  9.392095
46  9.223223  9.988546
47  8.799165  9.341935
48 10.285030 10.889897
49  7.499768  8.206951
50  7.510437  8.174375
51  8.023584  8.654680
52  8.349020  8.709274
53  8.507704  8.770455
54  8.336092  8.661712
55  8.380669  8.951269
56  8.776927  9.087411
57  7.399784  8.079899
58  7.852468  8.843239
59  8.139146  9.027143
60  7.662124  8.468211
61  7.553788  8.290680


plot(ladderize(primate_tree05), cex=0.6); axisPhylo()

nodelabels(pch = 21, cex=(MLreconstruction2$ace*0.33))  # 0.33 is a multiplier (i.e. scales) of the reconstructed value to better display on the tree. You can modify this value.
 
write.csv(MLreconstruction2$ace, "MLreconstruction.csv")
write.csv(MLreconstruction2$CI95, "MLreconstruction_CI95.csv")

# What is the reconstructed value for the last common ancestor of gorillas and chimpanzees?
mrca(primate_tree05)["Gorilla_gorilla", "Pan_troglodytes"]
  #48

# What is the reconstructed value for the last common ancestor of ring tailed lemurs and ruffed lemurs?
mrca(primate_tree05)["Lemur_catta", "Varecia_variegata"]
  #60

#----------------------------------------------------------
# Using phytools
#----------------------------------------------------------

# Option 1

aa <- fastAnc(primate_tree05, female_mass, CI = TRUE)

# Note that female mass is not log transformed

aa

$ace
       32        33        34        35        36        37        38 
 8922.262 11195.869 17132.090 11242.297  8251.880  5858.007  4779.014 
       39        40        41        42        43        44        45 
 5661.930  8289.524  8161.404  8915.342  9879.350  9045.248 10037.715 
       46        47        48        49        50        51        52 
27929.275 12068.341 52664.350  4783.861  4606.927  5625.679  5282.131 
       53        54        55        56        57        58        59 
 5734.735  5053.866  6674.561  7735.889  4262.100  5472.658  6042.684 
       60        61 
 3781.372  3141.573 

$CI95
          [,1]     [,2]
32 -17747.7774 35592.30
33  -7939.7909 30331.53
34   2702.4698 31561.71
35   -254.9263 22739.52
36  -1209.1828 17712.94
37  -2682.5264 14398.54
38  -3596.9722 13155.00
39  -2984.7619 14308.62
40   -927.8818 17506.93
41  -1124.6811 17447.49
42    496.4480 17334.24
43   -484.5851 20243.28
44  -1615.9400 19706.44
45   -622.5296 20697.96
46  14427.3745 41431.18
47   2492.7266 21643.95
48  41993.2093 63335.49
49  -7692.2076 17259.93
50  -7106.1890 16320.04
51  -5508.1332 16759.49
52  -1073.5029 11637.77
53   1099.2532 10370.22
54   -690.7370 10798.47
55  -3392.0020 16741.12
56   2258.2926 13213.49
57  -7736.4332 16260.63
58 -12006.6468 22951.96
59  -9623.4664 21708.83
60 -10439.7242 18002.47
61  -9858.7753 16141.92

# Option 2

anc_BM <- anc.ML(primate_tree05, female_mass, model = "BM")

anc_BM

$sig2
        
4165726 

$ace
       32        33        34        35        36        37        38 
 8922.262 11195.869 17132.090 11242.297  8251.880  5858.007  4779.014 
       39        40        41        42        43        44        45 
 5661.930  8289.524  8161.404  8915.342  9879.350  9045.248 10037.715 
       46        47        48        49        50        51        52 
27929.275 12068.341 52664.350  4783.861  4606.927  5625.679  5282.131 
       53        54        55        56        57        58        59 
 5734.735  5053.866  6674.561  7735.889  4262.100  5472.658  6042.684 
       60        61 
 3781.372  3141.573 

$logLik
[1] -599.9908

$counts
function gradient 
       1        1 

$convergence
[1] 0

$message
[1] "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"

$model
[1] "BM"

attr(,"class")
[1] "anc.ML"


anc_BM <- anc.ML(primate_tree05, log(female_mass), model = "BM")


#------------------------------------------------------------------
# Reconstructing states using maximum likelihood and a OU model
#------------------------------------------------------------------

# Revell warns that the OU modeol "has not be thoroughly tested & some bugs were reported for an earlier version"

anc_OU <- anc.ML(primate_tree05, female_mass, CI = TRUE, model = "OU")
anc_OU

$sig2
[1] 4236332

$alpha
[1] 1e-08

$ace
       32        33        34        35        36        37        38 
 8922.262 11195.869 17132.090 11242.297  8251.880  5858.007  4779.014 
       39        40        41        42        43        44        45 
 5661.930  8289.524  8161.404  8915.342  9879.350  9045.248 10037.715 
       46        47        48        49        50        51        52 
27929.275 12068.341 52664.350  4783.861  4606.927  5625.679  5282.131 
       53        54        55        56        57        58        59 
 5734.735  5053.866  6674.561  7735.889  4262.100  5472.658  6042.684 
       60        61 
 3781.372  3141.573 

$logLik
[1] -599.995

$counts
function gradient 
       2        2 

$convergence
[1] 0

$message
[1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

$model
[1] "OU"

attr(,"class")
[1] "anc.ML"
 
#A few things from phytools
phenogram(primate_tree05, log(female_mass), spread.labels = TRUE)
fancyTree(primate_tree05, type = "phenogram95", x = log(female_mass))
